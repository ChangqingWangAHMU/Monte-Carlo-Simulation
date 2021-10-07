function amount_iron=CA_GDF(N_spheres,u,upsilon_ca,beta_ca)
%The amount of iron per hepatocyte
%input:N_spheres--number of spheres;
%      u--location parameter in GDF;
%      upsilon_ca--shape parameter in GDF;
%      beta_ca--scale parameter in GDF
%output:amount_iron--amount of iron per hepatocyte

% x_ca=0:0.1:1;%GDF_CA
% GDF_CA=gampdf(x_ca-u,upsilon_ca,beta_ca);
% GDF_CA=GDF_CA/sum(GDF_CA);
% figure;plot(x_ca,GDF_CA,'*');xlabel('Normailized iron loading');ylabel('Relative frequency');title('GDF\_CA');

rng(3,'v5uniform');x_cdf_ca=rand(1,64);
amount_iron=zeros(1,64);
for i=1:64
    amount_iron(1,i)=fzero(@(x) gamcdf(x-u,upsilon_ca,beta_ca)-x_cdf_ca(1,i),0);
end
% figure;hist(amount_iron,10);xlabel('Normailized iron loading');ylabel('Frequency');title('Cellular Anisotropy');

amount_iron=round(amount_iron/sum(amount_iron)*N_spheres);
amount_iron(end)=amount_iron(end)+(N_spheres-sum(amount_iron(:)));