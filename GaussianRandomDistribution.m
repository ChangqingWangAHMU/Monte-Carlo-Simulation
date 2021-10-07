function r=GaussianRandomDistribution(N_spheres,size_sphere,HIC)
%Gaussian random sphere distribution
%input:N_spheres--number of spheres;
%      size_sphere--radius of spheres;
%output:r--gaussian random distributed spheres' position without collision

% r=mvnrnd(zeros(1,3),15^2*eye(3),length(size_sphere));%gaussian random distribution (mean and standard deviation unknown here)
rng(2,'v5uniform');r=-40+80.*rand(N_spheres,3);%gaussian random distributed spheres position (by uniform distribution here)

% %sphere collision
% temp=[-40:20:40];
% for i=1:4
%     for j=1:4
%         for k=1:4
%             index_he=find(r(:,1)>temp(i)&r(:,1)<temp(i+1)&r(:,2)>temp(j)&r(:,2)<temp(j+1)&r(:,3)>temp(k)&r(:,3)<temp(k+1));
%             r_temp=r(index_he,:);
%             r_temp=CheckCollision(r_temp,size_sphere(index_he),i,j,k);
%             r(index_he,:)=r_temp;
%         end
%     end
% end

%sphere collision is checked directly 
r=CheckCollision(r,size_sphere);
save(['r1_Sphere_HIC' num2str(HIC) '.mat']);
scatter3sphere(r,size_sphere);title('Random distribution');