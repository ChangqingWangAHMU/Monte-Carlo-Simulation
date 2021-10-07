%Load inputs for MonteCarloSimulation. 
%Created by Changqing Wang (10/23/2015)
function Loadinputs(HIC)
B0=1.5;%magentic field, in Tesla
% HIC=20;%hepatic iron concentration, in mg/g
WDR=4.1;%tissue wet-to-dry weight ratio
chi_F=1.2e-6;%intrinsic suscepibility or iron, in m^3/Kg_Fe
N_spheres=0;%number of iron particles

%Iron Morphology Description in Table 1
VF=0.00127*HIC+0.00368;%Volume fraction
beta_sd=0.0002965*HIC+0.09353;%Size distribution(MRM paper)
upsilon_sd=2.188*beta_sd+2.695;
% beta_sd=0.000733*HIC+0.0917;%Size distribution (Dissertation)
% upsilon_sd=1.933*beta_sd+2.724;

upsilon_nn=exp(0.4764*log(HIC)-0.1309);%Nearest neighbor
beta_nn=exp(-1.4029*log(upsilon_nn)+0.7404);
upsilon_ca=0.01787*HIC+1.0245;%Cellular anisotropy
beta_ca=exp(0.3751*log(HIC)-3.3663);
u_sd=0.05;u=0;

%spere size distribution
VF_target=0;
while VF_target<VF
    N_spheres=N_spheres+200;
    rng(1,'v5uniform');x_cdf_sd=rand(1,N_spheres);
    size_sphere=zeros(1,N_spheres);
    for i=1:N_spheres
        size_sphere(1,i)=fzero(@(x) gamcdf(x-u_sd,upsilon_sd,beta_sd)-x_cdf_sd(1,i),0);
    end
    VF_target=4/3*pi*sum(size_sphere.^3)/80^3;%volume fraction
end

A_he=5.2;B_he=0.57;C_he=36;%Hepatocyte (MRM paper)
A_tot=5.1;B_tot=0.57;C_tot=60;%Total (MRM paper)
% A_he=6.79;B_he=0.69;C_he=36;%Hepatocyte (Dissertation)
% A_tot=5.98;B_tot=0.63;C_tot=60;%Total (Dissertation)

ironscore_hep=C_he/(1+A_he*exp(-B_he*log(HIC)));
ironscore_tot=C_tot/(1+A_tot*exp(-B_tot*log(HIC)));
VF_he=round(ironscore_hep/ironscore_tot*100)/100;

%iron loading distribution and nearst neighbor distance distribution
amount_iron=CA_GDF(N_spheres,u,upsilon_ca,beta_ca);
distance_nn=NN_GDF(N_spheres,u,upsilon_nn,beta_nn);
amount_iron_sin=floor(N_spheres*(1-VF_he)/18*ones(1,18));
amount_iron_he=CA_GDF(N_spheres-sum(amount_iron_sin),u,upsilon_ca,beta_ca);
save(['Inputs_HIC' num2str(HIC) '.mat']);


