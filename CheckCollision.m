function r=CheckCollision(r,size_sphere,ii,jj,kk)
%Check collision between spheres
%input:r-sphere position;
%      size_sphere-sphere radius
%output:r--sphere position after checking collisions
temp=[-40:20:40];
if nargin==2%random gaussian distribution
    ii=1;jj=1;kk=1;
    spacing=80;
else 
    spacing=20;
end

for i=2:size(r,1)
    collision_index=find(size_sphere(i).*ones(1,i-1)+size_sphere(1:i-1)-sqrt(sum((repmat(r(i,:),i-1,1)-r(1:i-1,:)).^2,2))'>1e-10);
    while ~isempty(collision_index)
        k=collision_index(1);%step1:placing the concerned sphere at the surface of the colliding one
        radius1=size_sphere(k);radius2=size_sphere(i);
        theta=pi*rand(1);phi=2*pi*rand(1);
        r(i,:)=r(k,:)+[(radius1+radius2)*sin(theta)*cos(phi),(radius1+radius2)*sin(theta)*sin(phi),(radius1+radius2)*cos(theta)];
        collision_index1=find(size_sphere(i).*ones(1,i-1)+size_sphere(1:i-1)-sqrt(sum((repmat(r(i,:),i-1,1)-r(1:i-1,:)).^2,2))'>1e-10);%step 2:recheck collision
            
        if ~(length(collision_index1)<length(collision_index))||min(r(i,:))<-40||max(r(i,:))>40%step 3:collision unsolved, then new position
            r(i,:)=[temp(ii)+spacing.*rand(1) temp(jj)+spacing.*rand(1) temp(kk)+spacing.*rand(1)];
            collision_index=find(size_sphere(i).*ones(1,i-1)+size_sphere(1:i-1)-sqrt(sum((repmat(r(i,:),i-1,1)-r(1:i-1,:)).^2,2))'>1e-10);
            continue;
        end
        collision_index=collision_index1;
    end    
end