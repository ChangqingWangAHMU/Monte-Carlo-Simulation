function [r,distance_nn]=CellularAnisotropy_NN(N_spheres,amount_iron,size_sphere,distance_nn,HIC)
%Cellular anisotropy sphere distribution with nearest neighbor
%input:amount_iron--number of spheresfor each hepatocyte;
%      size_sphere--radius of spheres;
%      distance_nn--nearest neighbor distance
%output:cellular anisotropy with nearest neighbor distributed spheres' position without collision

p=0;%index of hepatocyte
m=0;%%index of sphere
tem=[-40:20:20];    
r=[];
index_ini=cumsum(amount_iron)+1;
index_ini=[1 index_ini(1:end-1)];
hwait=waitbar(0,'Processing');
for i=1:4
    for j=1:4
        for k=1:4
            p=p+1;
            r_temp=zeros(amount_iron(p),3);
            r_temp(1,:)=[tem(i)+20.*rand(1) tem(j)+20.*rand(1) tem(k)+20.*rand(1)];%generate the 1st reference sphere in each hepatocyte
            if amount_iron(p)==0
                continue;
            end
            for l=2:amount_iron(p)%generate subsequent spheres in each hepatocyte and check collision
                while 1
                    theta=pi*rand(1);phi=2*pi*rand(1);
                    r_temp(l,:)=r_temp(l-1,:)+[distance_nn(1,index_ini(p)+l-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+l-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+l-2)*cos(theta)];
                    if (r_temp(l,1)>tem(i))&&(r_temp(l,1)<tem(i)+20)&&(r_temp(l,2)>tem(j))&&(r_temp(l,2)<tem(j)+20)&&(r_temp(l,3)>tem(k))&&(r_temp(l,3)<tem(k)+20)
                        break;
                    end
                end
                
                collision_index=find(size_sphere(index_ini(p)+l-1).*ones(1,l-1)+size_sphere(index_ini(p):(index_ini(p)+l-2))...
                                     -sqrt(sum((repmat(r_temp(l,:),l-1,1)-r_temp(1:l-1,:)).^2,2))'>1e-10);
                m=0;
                while ~isempty(collision_index)
                    kc=collision_index(1);
                    while 1%step 1:placing the concerned sphere at the surface of the colliding one
                        radius1=size_sphere(index_ini(p)+kc-1);radius2=size_sphere(index_ini(p)+l-1);
                        theta=pi*rand(1);phi=2*pi*rand(1);
                        r_temp(l,:)=r_temp(kc,:)+[(radius1+radius2)*sin(theta)*cos(phi),(radius1+radius2)*sin(theta)*sin(phi),(radius1+radius2)*cos(theta)];
                        if (r_temp(l,1)>tem(i))&&(r_temp(l,1)<tem(i)+20)&&(r_temp(l,2)>tem(j))&&(r_temp(l,2)<tem(j)+20)&&(r_temp(l,3)>tem(k))&&(r_temp(l,3)<tem(k)+20)
                            break;
                        end
                    end
       
                    collision_index1=find(size_sphere(index_ini(p)+l-1).*ones(1,l-1)+size_sphere(index_ini(p):(index_ini(p)+l-2))...
                                     -sqrt(sum((repmat(r_temp(l,:),l-1,1)-r_temp(1:l-1,:)).^2,2))'>1e-10);%step 2:recheck collision
                    if ~(length(collision_index1)<length(collision_index))%step 3:collision unsolved, then new position
                        m=m+1;
                        while 1%subsequent sphere position with fixed diatance in hepatocyte
                            theta=pi*rand(1);phi=2*pi*rand(1);
                            r_temp(l,:)=r_temp(l-1,:)+[distance_nn(1,index_ini(p)+l-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+l-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+l-2)*cos(theta)];
                            if (r_temp(l,1)>tem(i))&&(r_temp(l,1)<tem(i)+20)&&(r_temp(l,2)>tem(j))&&(r_temp(l,2)<tem(j)+20)&&(r_temp(l,3)>tem(k))&&(r_temp(l,3)<tem(k)+20)
                                break;
                            end
                        end
                         
                        if m==50%new distance
                            distance_nn(1,index_ini(p)+l-2)=8*rand(1);
                            while 1
                                theta=pi*rand(1);phi=2*pi*rand(1);
                                r_temp(l,:)=r_temp(l-1,:)+[distance_nn(1,index_ini(p)+l-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+l-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+l-2)*cos(theta)];
                                if (r_temp(l,1)>tem(i))&&(r_temp(l,1)<tem(i)+20)&&(r_temp(l,2)>tem(j))&&(r_temp(l,2)<tem(j)+20)&&(r_temp(l,3)>tem(k))&&(r_temp(l,3)<tem(k)+20)
                                    break;
                                end
                            end
                        elseif m==1000%new reference
                            while 1
                                r_temp(l,:)=[tem(i)+20.*rand(1) tem(j)+20.*rand(1) tem(k)+20.*rand(1)];
                                if (r_temp(l,1)>tem(i))&&(r_temp(l,1)<tem(i)+20)&&(r_temp(l,2)>tem(j))&&(r_temp(l,2)<tem(j)+20)&&(r_temp(l,3)>tem(k))&&(r_temp(l,3)<tem(k)+20)
                                    break;
                                end
                            end
                            break;
                        end
                        
                        collision_index=find(size_sphere(index_ini(p)+l-1).*ones(1,l-1)+size_sphere(index_ini(p):(index_ini(p)+l-2))...
                                     -sqrt(sum((repmat(r_temp(l,:),l-1,1)-r_temp(1:l-1,:)).^2,2))'>1e-10);
                        continue;
                    end
                    collision_index=collision_index1; 
                end
            end
            r=[r;r_temp]; 
            waitbar(p/64,hwait,[num2str(fix(p/64*100)) '% completed']);
        end
    end
end
close(hwait);
save(['r3_Sphere_HIC' num2str(HIC) '.mat']);
% scatter3sphere(r,size_sphere);title('Cellular Anisotropy with NN');