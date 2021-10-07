function r=CellularAnisotropy(amount_iron,size_sphere,HIC)
%Cellular anisotropy sphere distribution
%input:amount_iron--amount of iron per hepatocyte
%       size_sphere--radius of spheres;
%output:cellular anisotropy distributed spheres' position without collision

index_ini=cumsum(amount_iron)+1;
index_ini=[1 index_ini(1:end-1)];

tem=[-40:20:20];    
p=0;%index of hepatocyte
r=[];
hwait=waitbar(0,'Processing');
for i=1:4
    for j=1:4
        for k=1:4
            p=p+1;
            r_temp=zeros(amount_iron(1,p),3);
            r_temp=[tem(i)+20.*rand(amount_iron(1,p),1) tem(j)+20.*rand(amount_iron(1,p),1) tem(k)+20.*rand(amount_iron(1,p),1)];
            r_temp=CheckCollision(r_temp,size_sphere(index_ini(p):index_ini(p)+size(r_temp,1)-1),i,j,k);
            r=[r;r_temp];
            waitbar(p/64,hwait,[num2str(fix(p/64*100)) '% completed']);
        end
    end
end
close(hwait);
save(['r2_Sphere_HIC' num2str(HIC) '.mat']);
scatter3sphere(r,size_sphere);title('Cellular Anisotropy');