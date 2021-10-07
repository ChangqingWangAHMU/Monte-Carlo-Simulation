function distance_nn=NN_GDF(N_spheres,u,upsilon_nn,beta_nn)
%Nearest neighbor distance distribution
%input:N_spheres--number of spheres;
%      u--location parameter in GDF;
%      upsilon_nn--shape parameter in GDF;
%      beta_nn--scale parameter in GDF
%output:distance_nn--nearest neighbor distance

% x_nn=0:0.5:8;%GDF_NN
% GDF_NN=gampdf(x_nn-u,upsilon_nn,beta_nn);
% GDF_NN=GDF_NN/sum(GDF_NN);
% figure;plot(x_nn,GDF_NN,'*');xlabel('Nearest Neighbor Distance (um)');ylabel('Relative frequency');title('GDF\_NN');

rng(4,'v5uniform');x_cdf_nn=rand(1,N_spheres);%small HIC: multiply 1-1e-3
distance_nn=zeros(1,N_spheres);
for i=1:N_spheres
    distance_nn(1,i)=fzero(@(x) gamcdf(x-u,upsilon_nn,beta_nn)-x_cdf_nn(1,i),0);
end
% figure;hist(distance_nn,100);xlabel('Nearest Neighbor Distance (um)');ylabel('Frequency');title('Cellular Anisotropy with NN');