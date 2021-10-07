function [r,distance_nn]=CellularAnisotropy_NNS(amount_iron_he,amount_iron_sin,size_sphere,distance_nn,HIC)
%Cellular anisotropy sphere distribution with nearest neighbor along
%sinusoids
%input:amount_iron_he--number of spheresfor each hepatocyte;
%      amount_iron_sin--number of spheresfor each sinusoid;
%      size_sphere--radius of spheres;
%      distance_nn--nearest neighbor distance
%output:cellular anisotropy with nearest neighbor along sinusoids distributed spheres' position without collision

tem2=[-20 0 -40 -20;20 0 -40 -20;0 -20 -40 -20;0 20 -40 -20; ...
        -20 -20 -20 0;-20 20 -20 0;20 -20 -20 0;20 20 -20 0;0 0 -20 0; ...
        -20 0 0 20;20 0 0 20;0 -20 0 20;0 20 0 20; ...
        -20 -20 20 40;-20 20 20 40;20 -20 20 40;20 20 20 40;0 0 20 40];%cylinders centers and height
%% %%%%%%%%%%%%%%%%%%%%%Spheres in Hepatocyte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_he=[];
p=0;%index of hepatocyte
tem1=[-40:20:20];%hepatocytes' position
index_ini=cumsum(amount_iron_he)+1;
index_ini=[1 index_ini(1:end-1)];
hwait=waitbar(0,'Processing Hepatocyte');
for i=1:4
    if i==1%cylinder position in ith slice
        temp3=tem2(1:4,:);
    elseif i==2
        temp3=tem2(5:9,:);
    elseif i==3
        temp3=tem2(10:13,:);
    elseif i==4 
        temp3=tem2(14:18,:);
    end
    for j=1:4
        for k=1:4
            p=p+1;
            q=1;
            r_temp=zeros(amount_iron_he(p),3);
            
            while 1%generate 1st reference sphere in each hepatocyte (not within cylinder)
                r_temp(1,:)=[tem1(k)+20.*rand(1) tem1(j)+20.*rand(1) tem1(i)+20.*rand(1)];
                incylinder=find(sum((repmat(r_temp(1,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                if isempty(incylinder)
                    break;
                end
            end
            
            while q<amount_iron_he(p)%generate subsequent sphere position and check collision
                q=q+1;
                
                while 1%subsequent sphere position with fixed diatance in hepatocyte while not in cylinder
                    theta=pi*rand(1);phi=2*pi*rand(1);
                    r_temp(q,:)=r_temp(q-1,:)+[distance_nn(1,index_ini(p)+q-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+q-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+q-2)*cos(theta)];
                    incylinder=find(sum((repmat(r_temp(q,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                    if isempty(incylinder)&&r_temp(q,1)>tem1(k)&&r_temp(q,1)<tem1(k)+20&&r_temp(q,2)>tem1(j)&&r_temp(q,2)<tem1(j)+20&&r_temp(q,3)>tem1(i)&&r_temp(q,3)<tem1(i)+20
                        break;
                    end
                end
                
                collision_index=find(size_sphere(index_ini(p)+q-1).*ones(1,q-1)+size_sphere(index_ini(p):(index_ini(p)+q-2))...
                                     -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);%check collision
                m=0;
                while ~isempty(collision_index)
                    kc=collision_index(1);
                    while 1%step 1:placing the concerned sphere at the surface of the colliding one 
                        radius1=size_sphere(index_ini(p)+kc-1);radius2=size_sphere(index_ini(p)+q-1);
                        theta=pi*rand(1);phi=2*pi*rand(1);
                        r_temp(q,:)=r_temp(kc,:)+[(radius1+radius2)*sin(theta)*cos(phi),(radius1+radius2)*sin(theta)*sin(phi),(radius1+radius2)*cos(theta)];
                        incylinder=find(sum((repmat(r_temp(q,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                        if isempty(incylinder)&&r_temp(q,1)>tem1(k)&&r_temp(q,1)<tem1(k)+20&&r_temp(q,2)>tem1(j)&&r_temp(q,2)<tem1(j)+20&&r_temp(q,3)>tem1(i)&&r_temp(q,3)<tem1(i)+20
                            break;
                        end
                    end
                            
                    collision_index1=find(size_sphere(index_ini(p)+q-1).*ones(1,q-1)+size_sphere(index_ini(p):(index_ini(p)+q-2))...
                                     -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);%step 2: recheck collision
                    if ~(length(collision_index1)<length(collision_index))%step 3:collision unsolved
                        m=m+1;
                        while 1%subsequent sphere position with fixed diatance in hepatocyte while not in cylinder
                            theta=pi*rand(1);phi=2*pi*rand(1);
                            r_temp(q,:)=r_temp(q-1,:)+[distance_nn(1,index_ini(p)+q-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+q-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+q-2)*cos(theta)];
                            incylinder=find(sum((repmat(r_temp(q,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                            if isempty(incylinder)&&r_temp(q,1)>tem1(k)&&r_temp(q,1)<tem1(k)+20&&r_temp(q,2)>tem1(j)&&r_temp(q,2)<tem1(j)+20&&r_temp(q,3)>tem1(i)&&r_temp(q,3)<tem1(i)+20
                                break;
                            end
                        end
                            
                        if m==50%new distance
                            distance_nn(1,index_ini(p)+q-2)=8*rand(1);
                            while 1
                                theta=pi*rand(1);phi=2*pi*rand(1);
                                r_temp(q,:)=r_temp(q-1,:)+[distance_nn(1,index_ini(p)+q-2)*sin(theta)*cos(phi),distance_nn(1,index_ini(p)+q-2)*sin(theta)*sin(phi),distance_nn(1,index_ini(p)+q-2)*cos(theta)];
                                incylinder=find(sum((repmat(r_temp(q,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                                if isempty(incylinder)&&r_temp(q,1)>tem1(k)&&r_temp(q,1)<tem1(k)+20&&r_temp(q,2)>tem1(j)&&r_temp(q,2)<tem1(j)+20&&r_temp(q,3)>tem1(i)&&r_temp(q,3)<tem1(i)+20
                                    break;
                                end
                            end
                        elseif m==1000%new reference sphere
                            while 1
                                r_temp(q,:)=[tem1(k)+20.*rand(1) tem1(j)+20.*rand(1) tem1(i)+20.*rand(1)];
                                incylinder=find(sum((repmat(r_temp(q,1:2),size(temp3,1),1)-temp3(:,1:2)).^2,2)<25);
                                if isempty(incylinder)&&r_temp(q,1)>tem1(k)&&r_temp(q,1)<tem1(k)+20&&r_temp(q,2)>tem1(j)&&r_temp(q,2)<tem1(j)+20&&r_temp(q,3)>tem1(i)&&r_temp(q,3)<tem1(i)+20
                                    break;
                                end
                            end
                            break;
                        end
                            
                        collision_index=find(size_sphere(index_ini(p)+q-1).*ones(1,q-1)+size_sphere(index_ini(p):(index_ini(p)+q-2))...
                                                -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);
                        continue;
                    end
                    collision_index=collision_index1;         
                end
            end
            r_he=[r_he;r_temp];
            waitbar(p/64,hwait,[num2str(fix(p/64*100)) '% completed']);
        end
    end
end
close(hwait);

%% %%%%%%%%%%%%%%%%%%%%%Spheres in Sinusoid%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_sin=[];
size_sphere_sin=size_sphere(sum(amount_iron_he)+1:end);
distance_nn_sin=distance_nn(sum(amount_iron_he)-64+1:end);
hwait=waitbar(0,'Processing Sinusoid');
for i=1:size(tem2,1)
    r_temp=zeros(amount_iron_sin(i),3);
    q=1;
    
	phi_cy=2*pi*rand(1);%generate first sphere in each sinusoid
	rho_cy=5*rand(1);
	z_cy=20*rand(1);
	r_temp(1,:)=[tem2(i,1)+rho_cy.*cos(phi_cy) tem2(i,2)+rho_cy.*sin(phi_cy) tem2(i,3)+z_cy];
    while q<amount_iron_sin(i)
        q=q+1;
        
        while 1%subsequent sphere position with fixed distance in sinsoid
            theta=pi*rand(1);phi=2*pi*rand(1);
            r_temp(q,:)=r_temp(q-1,:)+[distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*cos(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*sin(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*cos(theta)];
            incylinder=find(sum((r_temp(q,1:2)-tem2(i,1:2)).^2,2)<25);
            if ~isempty(incylinder)&&r_temp(q,3)>tem2(i,3)&&r_temp(q,3)<tem2(i,4)
                break;
            end
        end
        
        collision_index=find(size_sphere_sin((i-1)*amount_iron_sin(i)+q).*ones(1,q-1)+size_sphere_sin((i-1)*amount_iron_sin(i)+1:(i-1)*amount_iron_sin(i)+q-1)...
                                     -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);%check collision
        m=0;
        while ~isempty(collision_index)
            kc=collision_index(1);
            while 1%step 1:placing the concerned sphere at the surface of the colliding one 
                radius1=size_sphere_sin((i-1)*amount_iron_sin(i)+kc);radius2=size_sphere_sin((i-1)*amount_iron_sin(i)+q);
                theta=pi*rand(1);phi=2*pi*rand(1);
                r_temp(q,:)=r_temp(kc,:)+[(radius1+radius2)*sin(theta)*cos(phi),(radius1+radius2)*sin(theta)*sin(phi),(radius1+radius2)*cos(theta)];
                incylinder=find(sum((r_temp(q,1:2)-tem2(i,1:2)).^2,2)<25);
                if ~isempty(incylinder)&&r_temp(q,3)>tem2(i,3)&&r_temp(q,3)<tem2(i,4)
                    break;
                end
            end
                
            collision_index1=find(size_sphere_sin((i-1)*amount_iron_sin(i)+q).*ones(1,q-1)+size_sphere_sin((i-1)*amount_iron_sin(i)+1:(i-1)*amount_iron_sin(i)+q-1)...
                                     -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);%step 2: recheck collision after moving the reference spheres
            if ~(length(collision_index1)<length(collision_index))%step 3:collision unsolved
                m=m+1;
                while 1%subsequent sphere position with fixed distance in sinsoid
                    theta=pi*rand(1);phi=2*pi*rand(1);
                    r_temp(q,:)=r_temp(q-1,:)+[distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*cos(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*sin(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*cos(theta)];
                    incylinder=find(sum((r_temp(q,1:2)-tem2(i,1:2)).^2,2)<25);
                    if ~isempty(incylinder)&&r_temp(q,3)>tem2(i,3)&&r_temp(q,3)<tem2(i,4)
                        break;
                    end
                end
                    
                if m==50%new distance
                    distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)=8*rand(1);
                    while 1
                        theta=pi*rand(1);phi=2*pi*rand(1);
                        r_temp(q,:)=r_temp(q-1,:)+[distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*cos(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*sin(theta)*sin(phi),distance_nn_sin(1,(i-1)*amount_iron_sin(i)+q-i)*cos(theta)];
                        incylinder=find(sum((r_temp(q,1:2)-tem2(i,1:2)).^2,2)<25);
                        if ~isempty(incylinder)&&r_temp(q,3)>tem2(i,3)&&r_temp(q,3)<tem2(i,4)
                            break;
                        end
                    end
                elseif m==1000%new reference sphere
                    phi_cy=2*pi*rand(1);
                    rho_cy=5*rand(1);
                    z_cy=20*rand(1);
                    r_temp(q,:)=[tem2(i,1)+rho_cy.*cos(phi_cy) tem2(i,2)+rho_cy.*sin(phi_cy) tem2(i,3)+z_cy];
                    break;
                end
                    
                collision_index=find(size_sphere_sin((i-1)*amount_iron_sin(i)+q).*ones(1,q-1)+size_sphere_sin((i-1)*amount_iron_sin(i)+1:(i-1)*amount_iron_sin(i)+q-1)...
                                     -sqrt(sum((repmat(r_temp(q,:),q-1,1)-r_temp(1:q-1,:)).^2,2))'>1e-10);
                continue;
            end
            collision_index=collision_index1;
        end  
    end
    r_sin=[r_sin;r_temp];
	waitbar(i/18,hwait,[num2str(fix(i/18*100)) '% completed']);
end
close(hwait);
distance_nn(sum(amount_iron_he)-64+1:end)=distance_nn_sin;
%% %%%%%%%%%%%%%%%%%%%%%%%%%Result%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=[r_he;r_sin];
save(['r4_Sphere_HIC' num2str(HIC) '.mat']);
scatter3sphere(r,size_sphere);title('Cellular Anisotropy with NNS');
