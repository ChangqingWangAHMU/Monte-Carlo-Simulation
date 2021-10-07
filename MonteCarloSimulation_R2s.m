%Implementation of the Monte carlo simulation for relaxivity-iron calibration in hepatic iron overload
%Software: Matlab2016b
%Author: ChangqingWang
%Date: 09-17-2021
%References: [1] Wang C, Reeder SB, Hernando D. Relaxivity-iron calibration in hepatic iron overload: Reproducibility and extension of a Monte Carlo model. NMR in Biomedicine. 2021;e4604. doi:10.1002/nbm.4604 
%            [2] Ghugre NR, Doyle EK, Storey P, Wood JC. Relaxivity-iron calibration in hepatic iron overload: predictions of a Monte Carlo model. Magn Reson Med. 2015;74(3):879-883.
%            [3] Ghugre NR, Wood JC. Relaxivity-iron calibration in hepatic iron overload: probing underlying biophysical mechanisms using a Monte Carlo model. Magn Reson Med. 2011;65(3):837-847.
%Note that unit of distance is um, and unit of time is msec

clear all;clc;
addpath(genpath(pwd));
if isempty(gcp('nocreate'))%start a parallel pool of 10 workers
    parpool(10);
else
    disp('Already initialized');
end
tic;
name_distribution=input('Please input name of distribuion:r1, r2, r3 or r4:\n','s');%r1:random distribution, r2:cellular anisotropy, r3:cellular anisotropy with nearest neighbor, r4:cellular anisotropy with nearest neighbor along with sinusoids
HIC=input('Please input HIC (mg/g):\n');%hepatic iron concentration
D=input('Please input D:\n');%diffusion coefficient, in um^2/msec
Loadinputs(HIC);
load(['Inputs_HIC' num2str(HIC) '.mat']);

%% %%%%%%%%%%%%%%%%%%%%%%%%Virtual Liver Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (name_distribution)
    case 'r1'
        % %(a) random spheres distribtuion
        r=GaussianRandomDistribution(N_spheres,size_sphere,HIC);
    case 'r2'
        % %(b) cellular anisotropy
        r=CellularAnisotropy(amount_iron,size_sphere,HIC); 
    case 'r3'
        % % (c) cellular anisotropy with nearest neighbor 
        [r,distance_nn]=CellularAnisotropy_NN(N_spheres,amount_iron,size_sphere,distance_nn,HIC);
    case 'r4'
        % %(d) cellular anisotropy with nearest neighbor along with sinusoids
        [r,distance_nn]=CellularAnisotropy_NNS(amount_iron_he,amount_iron_sin,size_sphere,distance_nn,HIC);
    otherwise
        disp('This is error!');
end
load([name_distribution '_Sphere_HIC' num2str(HIC) '.mat']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic field generation%%%%%%%%%%%%%%%%%%%%
chi_L=(HIC/WDR)*chi_F/VF_target;%particle susceptibility
[X,Y,Z]=ndgrid(-40:0.5:40,-40:0.5:40,-40:0.5:40);
X=X(:);Y=Y(:);Z=Z(:);
deltaB=zeros(1,161^3);
parfor i=1:161^3
    tem=(r-repmat([X(i) Y(i) Z(i)],N_spheres,1)).^2;
    deltaB(i)=sum(B0/3*chi_L*size_sphere.^3'.*(2*tem(:,3)-tem(:,1)-tem(:,2))./(sum(tem,2)).^(5/2));
end
X=reshape(X,161,161,161);Y=reshape(Y,161,161,161);Z=reshape(Z,161,161,161);deltaB=reshape(deltaB,161,161,161);
F=griddedInterpolant(X,Y,Z,deltaB,'cubic');%interpolating function of magnetic dipole field
save([name_distribution '_F_HIC' num2str(HIC) '.mat'],'F')
toc
clear X Y Z
% load([name_distribution '_F_HIC' num2str(HIC) '.mat']);%option for loading the generated data
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%Proton movement%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_protons=10*500;
time_interval=30;
time_step=0.5*1e-3;
N_timesteps=time_interval/time_step+1;
%proton mobolity:mobility
spmd
for ip=1:N_protons/10
    ip
    while 1%the first proton position (not in sphere)
        r_protons=-40+80*rand(1,3);
        insphere=find(sum((repmat(r_protons,N_spheres,1)-r).^2,2)<size_sphere.^2');
        if isempty(insphere)
            break;
        end
    end
    %Brownian Motion
    x_tem=zeros(N_timesteps,3);
    dx=zeros(3,N_timesteps-1);
    [x_tem,dx]=brownian_motion_simulation(N_timesteps,D,time_step);
    r_temp=repmat(r_protons,N_timesteps,1)+x_tem;
    r_protons_m(:,:,ip)=r_temp;
    
    %nearby spheres table
    table_nearbysphere=cell(time_interval/0.2+1,1);
    for it=1:size(table_nearbysphere,1)
        table_nearbysphere{it,1}=find(sum((repmat(r_protons_m(400*(it-1)+1,:,ip),N_spheres,1)-r).^2,2)-size_sphere.^2'<4)';
    end
    
    %proton-sphere collision
    for it=1:N_timesteps-1
        index_len40=r_protons_m(it+1,:,ip)<-40;
        index_mt40=r_protons_m(it+1,:,ip)>40;
        if sum(index_len40>0)||sum(index_mt40>0)
            r_protons_m(it+1,index_len40,ip)=mod(r_protons_m(it+1,index_len40,ip),40);%periodic boundary
            r_protons_m(it+2:end,index_len40,ip)=repmat(r_protons_m(it+1,index_len40,ip),N_timesteps-it-1,1)+cumsum(dx(index_len40,it+1:end)');
            r_protons_m(it+1,index_mt40,ip)=mod(r_protons_m(it+1,index_mt40,ip),40)-40;%periodic boundary
            r_protons_m(it+2:end,index_mt40,ip)=repmat(r_protons_m(it+1,index_mt40,ip),N_timesteps-it-1,1)+cumsum(dx(index_mt40,it+1:end)');
            table_nearbysphere{floor((it-1)/400)+1}=find(sum((repmat(r_protons_m(it+1,:,ip),N_spheres,1)-r).^2,2)-size_sphere.^2'<4)';%update sphere table
            continue;
        end

        for itable=table_nearbysphere{floor((it-1)/400)+1}
            A=norm(r_protons_m(it+1,:,ip)-r_protons_m(it,:,ip))^2;%intersection between line segment and sphere
            B=2*(r_protons_m(it+1,:,ip)-r_protons_m(it,:,ip))*(r_protons_m(it,:,ip)-r(itable,:))';
            C=norm(r_protons_m(it,:,ip)-r(itable,:))^2-size_sphere(itable)^2;
            DD=B^2-4*A*C;
            if DD<0
                continue;
            elseif (-B-sqrt(DD))/(2*A)>0&&(-B-sqrt(DD))/(2*A)<1
                r_tem=r_protons_m(it,:,ip)+(-B-sqrt(DD))/(2*A)*(r_protons_m(it+1,:,ip)-r_protons_m(it,:,ip));
                r_protons_m(it+1:end,:,ip)=r_protons_m(it+1:end,:,ip)-repmat((r_protons_m(it+1,:,ip)-r_tem),N_timesteps-it,1);
            end
        end
    end
end
end
toc

for i=1:10
    r_protons_m_total(:,:,:,i)=r_protons_m{i};
end
r_protons_m=reshape(r_protons_m_total,N_timesteps,3,N_protons);
clear r_protons_m_total;

%% %%%%%%%%%%%%%%%%%%%%%%%%%Phase Accrual%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ratio_gyromagnetic=2.675*1e+5;%ms^(-1).T^(-1)
Phase=zeros(size(r_protons_m,3),N_timesteps);
parfor it=1:N_timesteps
    Phase(:,it)=B0+F(squeeze(r_protons_m(it,:,:))');
end
Phase=cumsum(Phase,2);
Phase=Ratio_gyromagnetic*time_step*Phase;%Equation [6] in reference [3]
clear r_protons_m;
delete(gcp);

%% %%%%%%%%%%%%%%%%%%%%%%%%%MR Signal Synthesis%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T20=50;%msec
t=1:time_step*1e+3:time_interval;
S0=1;

for i=1:size(Phase,1)
    s(i,:)=S0*exp(-t/T20+1i.*Phase(i,t/time_step+1));
end
fun=@(x,t) x(1)*exp(-t*x(2));
x = lsqcurvefit(fun,[1 1],t,abs(sum(s,1)))
M0=x(1);R2s=x(2);
Signal_R2s=sum(s,1);
save([name_distribution '_Phase_D' num2str(D/0.19) '_HIC' num2str(HIC) '.mat'],'Phase','Signal_R2s','-v7.3');
save([name_distribution '_R2s_Result_D' num2str(D/0.19) '_HIC' num2str(HIC) '.mat'],'x');
clear s;
run MonteCarloSimulation_R2.m
toc