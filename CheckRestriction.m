function  [r_protons_m,dx]=CheckRestriction(r_protons,dx)
%check restriction for proton mobility (not pass through boundary)
syms x4 y4
tem1=[-40:20:40];
tem2=[-20 0 -40 -20;20 0 -40 -20;0 -20 -40 -20;0 20 -40 -20; ...
        -20 -20 -20 0;-20 20 -20 0;20 -20 -20 0;20 20 -20 0;0 0 -20 0; ...
        -20 0 0 20;20 0 0 20;0 -20 0 20;0 20 0 20; ...
        -20 -20 20 40;-20 20 20 40;20 -20 20 40;20 20 20 40;0 0 20 40];% tem2 description of sinusoid
tem3(:,:,1)=[1 3 4 1;3 7 6 4;2 6 7 5;1 2 5 1];%tem3 description of inetersection sinusoid and hepatocyte
tem3(:,:,2)=[3 4 3 4;2 6 7 5;3 7 6 4;2 5 2 5];
tem3(:,:,3)=[1 3 4 1;3 7 6 4;2 6 7 5;1 2 5 1];
tem3(:,:,4)=[3 4 3 4;2 6 7 5;3 7 6 4;2 5 2 5];
[x,y,z]=meshgrid(1:4,1:4,1:4);
xyz=16*(x-1)+4*(y-1)+z;
    
r_protons_1=r_protons(1,:);
for index_he=1:64
    i=x(find(xyz==index_he));j=y(find(xyz==index_he));k=z(find(xyz==index_he));
    if r_protons_1(1)>tem1(i)&&r_protons_1(1)<tem1(i+1)&&r_protons_1(2)>tem1(j)...
            &&r_protons_1(2)<tem1(j+1)&&r_protons_1(3)>tem1(k)&&r_protons_1(3)<tem1(k+1)%decide whether in hepatocyte
        incylinder=find(sum((repmat(r_protons_1(1,1:2),size(tem2,1),1)-tem2(:,1:2)).^2,2)<25 & ...
                        r_protons_1(1,3)>tem2(:,3) & r_protons_1(1,3)<tem2(:,4));
        if  ~isempty(incylinder) %if in sinusiod
%             disp(['In sinusoid (' num2str(incylinder) ')']);
            index=find(sum((repmat(tem2(incylinder(1),1:2),size(r_protons,1),1)-r_protons(:,1:2)).^2,2)>25 | ...
                       r_protons(:,3)<tem2(incylinder(1),3) | r_protons(:,3)>tem2(incylinder(1),4),1);
            while ~isempty(index)
                if sum((tem2(incylinder(1),1:2)-r_protons(index,1:2)).^2,2)>25%reflection in cylinder, x-y case
                    A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                    B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-tem2(incylinder(1),1:2))';
                    C=norm(tem2(incylinder(1),1:2)-r_protons(index-1,1:2))^2-5^2;
                    u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                    if u1>0 && u1<1
                        u=u1;
                    else
                        u=u2;
                    end
                    x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                    y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                    eq1=(r_protons(index,2)+y4)/2*(y2-tem2(incylinder(1),2))==... 
                        -(x2-tem2(incylinder(1),1))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem2(incylinder(1),2));
                    eq2=(y4-r_protons(index,2))*(x2-tem2(incylinder(1),1))==(x4-r_protons(index,1))*(y2-tem2(incylinder(1),2));
                    [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                    r_protons(index,1)=double(solx4);
                    r_protons(index,2)=double(soly4);
                end
                if r_protons(index,3)<tem2(incylinder(1),3)%reflection in cylinder, z case
                    r_protons(index,3)=2*tem2(incylinder(1),3)-r_protons(index,3);
                end
                if r_protons(index,3)>tem2(incylinder(1),4)
                    r_protons(index,3)=2*tem2(incylinder(1),4)-r_protons(index,3);
                end
                r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                index=find(sum((repmat(tem2(incylinder(1),1:2),size(r_protons,1),1)-r_protons(:,1:2)).^2,2)>25 | ...
                            r_protons(:,3)<tem2(incylinder(1),3) | r_protons(:,3)>tem2(incylinder(1),4),1);
            end
            break;
        else%if in hepatocyte
%             disp(['In hepatocyte (' num2str(i) ',' num2str(j) ',' num2str(k) ')']);
            switch tem3(i,j,k)  
                case 1 %for case without intersection with cylidner
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                       r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                       r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1),1);
                    while ~isempty(index)
                    if r_protons(index,1)<tem1(i)%reflection in cube, x case
                        r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                    end
                    if r_protons(index,1)>tem1(i+1)
                        r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                    end
                    if r_protons(index,2)<tem1(j)%reflection in cube, y case
                        r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                    end
                    if r_protons(index,2)>tem1(j+1)
                        r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                    end
                    if r_protons(index,3)<tem1(k) %reflection in cube, z case
                        r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                    end
                    if r_protons(index,3)>tem1(k+1)
                        r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                    end
                    r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                                r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                                r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1),1);
                    end
                    break;
                case 2 %for case upper-left intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                        r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                        r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                        sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i) tem1(j+1)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i) tem1(j+1)])';
                            C=norm([tem1(i) tem1(j+1)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j+1))==... 
                                -(x2-tem1(i))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j+1));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i))==(x4-r_protons(index,1))*(y2-tem1(j+1));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                            r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                            r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                            sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    end
                    break;
                case 3 %for case upper-right intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                        r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                        r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                        sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i+1) tem1(j+1)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i+1) tem1(j+1)])';
                            C=norm([tem1(i+1) tem1(j+1)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j+1))==... 
                                -(x2-tem1(i+1))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j+1));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i+1))==(x4-r_protons(index,1))*(y2-tem1(j+1));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                            r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                            r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                            sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    end
                    break;
                case 4 %for case lower-right intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                       r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                       r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                       sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i+1) tem1(j)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i+1) tem1(j)])';
                            C=norm([tem1(i+1) tem1(j)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j))==... 
                                -(x2-tem1(i+1))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i+1))==(x4-r_protons(index,1))*(y2-tem1(j));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                            r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                            r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                            sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                    end
                    break;
                case 5 %for case lower-left intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                       r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                       r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                       sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i) tem1(j)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i) tem1(j)])';
                            C=norm([tem1(i) tem1(j)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j))==... 
                                -(x2-tem1(i))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i))==(x4-r_protons(index,1))*(y2-tem1(j));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                           r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                           r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                           sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                    end
                    break;
                case 6 %for case lower-left and upper-right intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                        r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                        r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                        sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j)],size(r_protons,1),1)).^2,2)<25 |...
                        sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i) tem1(j)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i) tem1(j)])';
                            C=norm([tem1(i) tem1(j)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j))==... 
                                -(x2-tem1(i))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i))==(x4-r_protons(index,1))*(y2-tem1(j));
                            [solx4,soly4]=solve([eq1,eq2],x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i+1) tem1(j+1)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i+1) tem1(j+1)])';
                            C=norm([tem1(i+1) tem1(j+1)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j+1))==... 
                                -(x2-tem1(i+1))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j+1));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i+1))==(x4-r_protons(index,1))*(y2-tem1(j+1));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                            r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                            r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                            sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j)],size(r_protons,1),1)).^2,2)<25 |...
                            sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j+1)],size(r_protons,1),1)).^2,2)<25,1);
                    end
                    break;
                case 7 %for case upper-left and lower-right intersection
                    index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                        r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                        r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                        sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j+1)],size(r_protons,1),1)).^2,2)<25 |...
                        sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                    while ~isempty(index)
                        if r_protons(index,1)<tem1(i)%reflection in cube, x case
                            r_protons(index,1)=2*tem1(i)-r_protons(index,1);
                        end
                        if r_protons(index,1)>tem1(i+1)
                            r_protons(index,1)=2*tem1(i+1)-r_protons(index,1);
                        end
                        if r_protons(index,2)<tem1(j)%reflection in cube, y case
                            r_protons(index,2)=2*tem1(j)-r_protons(index,2);
                        end
                        if r_protons(index,2)>tem1(j+1)
                            r_protons(index,2)=2*tem1(j+1)-r_protons(index,2);
                        end
                        if r_protons(index,3)<tem1(k) %reflection in cube, z case
                            r_protons(index,3)=2*tem1(k)-r_protons(index,3);
                        end
                        if r_protons(index,3)>tem1(k+1)
                            r_protons(index,3)=2*tem1(k+1)-r_protons(index,3);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i) tem1(j+1)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i) tem1(j+1)])';
                            C=norm([tem1(i) tem1(j+1)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j+1))==... 
                                -(x2-tem1(i))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j+1));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i))==(x4-r_protons(index,1))*(y2-tem1(j+1));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        if norm(r_protons(index,1:2)-[tem1(i+1) tem1(j)])<5% intersection of sinusoid and hepatocyte
                            A=norm(r_protons(index-1,1:2)-r_protons(index,1:2))^2;
                            B=2*(r_protons(index,1:2)-r_protons(index-1,1:2))*(r_protons(index-1,1:2)-[tem1(i+1) tem1(j)])';
                            C=norm([tem1(i+1) tem1(j)]-r_protons(index-1,1:2))^2-5^2;
                            u1=(-B+sqrt(B^2-4*A*C))/(2*A);u2=(-B-sqrt(B^2-4*A*C))/(2*A);
                            if u1>0 && u1<1
                                u=u1;
                            else
                                u=u2;
                            end
                            x2=r_protons(index-1,1)+u*(r_protons(index,1)-r_protons(index-1,1));
                            y2=r_protons(index-1,2)+u*(r_protons(index,2)-r_protons(index-1,2));%(x2,y2)the intersection point
                            eq1=(r_protons(index,2)+y4)/2*(y2-tem1(j))==... 
                                -(x2-tem1(i+1))*((r_protons(index,1)+x4)/2-x2)+y2*(y2-tem1(j));
                            eq2=(y4-r_protons(index,2))*(x2-tem1(i+1))==(x4-r_protons(index,1))*(y2-tem1(j));
                            [solx4,soly4]=solve(eq1,eq2,x4,y4);%point after reflection
                            r_protons(index,1)=double(solx4);
                            r_protons(index,2)=double(soly4);
                        end
                        r_protons(index+1:end,:)=repmat(r_protons(index,:),size(r_protons,1)-index,1)+cumsum(dx(:,index:end),2)';
                        index=find(r_protons(:,1)<tem1(i) | r_protons(:,1)>tem1(i+1) | ...
                            r_protons(:,2)<tem1(j) | r_protons(:,2)>tem1(j+1) | ...
                            r_protons(:,3)<tem1(k) | r_protons(:,3)>tem1(k+1) | ...
                            sum((r_protons(:,1:2)-repmat([tem1(i) tem1(j+1)],size(r_protons,1),1)).^2,2)<25 |...
                            sum((r_protons(:,1:2)-repmat([tem1(i+1) tem1(j)],size(r_protons,1),1)).^2,2)<25,1);
                   end
                   break;
            end
        end
    end
end    
r_protons_m=r_protons;
dx=(r_protons_m(2:end,:)-r_protons_m(1:end-1,:))';