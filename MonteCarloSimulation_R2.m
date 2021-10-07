%Monte Carlo Simulatvion For R2 based on Phase Generated by R2* Simulations
% clear all;clc;
% addpath(genpath(pwd));
% tic;
% name_distribution=input('Please input name of distribuion:r1, r2, r3 or r4:\n','s');
% HIC=input('Please input HIC (mg/g):\n');
% D=input('Please input D:\n');
% load([name_distribution '_Phase_Res_D' num2str(D/0.19) '_HIC' num2str(HIC) '.mat']);
Phase=[Phase(:,1),Phase(:,2:end)-Phase(:,1:end-1)];%phase at time step K for each proton
TE_T2=exp(linspace(log(0.1),log(30),6));%TE

time_interval=30;
time_step=0.5*1e-3;
N_timesteps=time_interval/time_step+1;

for i=1:length(TE_T2)
    Phase_T2(:,i)=-sum(Phase(:,1:floor(TE_T2(i)/2/time_step)),2)'+sum(Phase(:,floor(TE_T2(i)/2/time_step)+1:floor(TE_T2(i)/time_step)),2)';
end

T20=50;
S0=1;

for i=1:size(Phase,1)
    s(i,:)=S0*exp(-TE_T2/T20+1i.*Phase_T2(i,:));
end
fun=@(x,t) x(1)*exp(-t*x(2))+x(3);
x1 = lsqcurvefit(fun,[abs(sum(s(:,1),1)), 1, 0],TE_T2,abs(sum(s,1)))
M0=x1(1);R2=x1(2);
save([name_distribution '_R2_Result_D' num2str(D/0.19) '_HIC' num2str(HIC) '.mat'],'x1');