function [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator( vx,vy,ele_size,order,end_cpt_factor,trim_factor,fig )
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\C_files');
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
addpath('C:\Users\Soham\Box Sync\Fibrin Modeling\German Formulation\C_files');
addpath('C:\Users\Soham\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
%% Description
% Network_Mesh_Generator : Generates mesh for random network.
%% OUTPUT

%% INPUT
% npt: No of ploygons in Voronoi diagram.
% V : Coordinates of the junction point and the ends of the straight lines meeting at the junction. (#rows = 4 and %columns = 3)

%% Code
MarkerSize=25;LineWidth=2;
%% generate mesh.
% Identify interior points
int_pt=[];key=1;
for j=1:size(vx,2)
    if (vx(1,j)~=0&&vx(1,j)~=1)&&(vy(1,j)~=0&&vy(1,j)~=1)
        int_pt(key,1)=vx(1,j); int_pt(key,2)=vy(1,j);key=key+1;
    end
    if (vx(2,j)~=0&&vx(2,j)~=1)&&(vy(2,j)~=0&&vy(2,j)~=1)
        int_pt(key,1)=vx(2,j); int_pt(key,2)=vy(2,j);key=key+1;
    end
end
int_pt=unique(int_pt,'rows','stable');
% Trim vertices at the interior points.
L=zeros(size(vx,2),1);
for i=1:size(vx,2)
    L(i)=norm([vx(1,i)-vx(2,i),vy(1,i)-vy(2,i)]);
end
% Len=min(L);% shortest vetrex.
Vx=[];Vy=[];key=1;
for i=1:size(vx,2)
    x1=vx(1,i); x2=vx(2,i); y1=vy(1,i); y2=vy(2,i);
    Len=norm([x1-x2,y1-y2]);
    for j=1:size(int_pt,1)
        if x1==int_pt(j,1) && y1==int_pt(j,2)
           key1=1;
           break
        else
            key1=0;
        end
    end
    for j=1:size(int_pt,1)
        if x2==int_pt(j,1) && y2==int_pt(j,2)
           key2=1;
           break
        else
           key2=0;
        end
    end
    if key1==1
        Vx(1,key)= x1+trim_factor*Len*(x2-x1)/L(i);
        Vy(1,key)= y1+trim_factor*Len*(y2-y1)/L(i);
    else
        Vx(1,key)= x1;
        Vy(1,key)= y1;
    end
    if key2==1
        Vx(2,key)= x2+trim_factor*Len*(x1-x2)/L(i);
        Vy(2,key)= y2+trim_factor*Len*(y1-y2)/L(i);
    else
        Vx(2,key)= x2;
        Vy(2,key)= y2;
    end
    key=key+1;
end
if fig=='on'
    figure
    plot(Vx,Vy,'b-','LineWidth',LineWidth);hold on
end


% Discretize segments with NURBS and store information of each patch.
Patch=struct;
Cpt=[];key_cpt=1;
for i=1:size(Vx,2)
    p1=[Vx(1,i);Vy(1,i)]; p2=[Vx(2,i);Vy(2,i)];
    l=norm(p2-p1);
    if ele_size<=l
        nele=ceil(norm(p1-p2)/ele_size);
    else
        nele=ceil(norm(p1-p2)/ele_size)+1;% if length of the fiber smaller than requested element size then 2 elements generated on the fiber.
    end
    nos=order+nele;
    % Form know vector.
    kv=zeros(1,nos+order+1);
    for j=size(kv,2)-order:size(kv,2)
        kv(1,j)=1;
    end
%     kv(1,order+2)=ele_size/10;kv(1,size(kv,2)-order-1)=1-ele_size/10;
    for j=2:nele
        kv(1,order+j)=(j-1)/nele;
    end
    % Form weights
    wts=ones(1,nos);
    Wts(1,key_cpt:key_cpt+nos-1)=wts;
    % Form Control points.
    cpt=zeros(4,nos);
    % Cpt next to the ends is added based on the length of the fiber.
    len=min(ele_size,l);
    cpt(1:2,1)=p1;cpt(1:2,nos)=p2;
    % Cpt next to the ends is added based on the length of the fiber.
%     cpt(1:2,2)=p1+(len*end_cpt_factor)*(p2-p1)/l;cpt(1:2,nos-1)=p2+(len*end_cpt_factor)*(p1-p2)/l;
%     for j=1:nos-4
%         cpt(1:2,2+j)=p1+(j/((nos+1)-4))*(p2-p1);
%     end
    for j=1:nos-2
        cpt(1:2,1+j)=p1+(j/((nos+1)-2))*(p2-p1);
    end
    num_cpt=[key_cpt:key_cpt+nos-1];
    Cpt(:,key_cpt:key_cpt+nos-1)=cpt;key_cpt=key_cpt+nos;
    if fig=='on'
        plot(cpt(1,:),cpt(2,:),'r.','MarkerSize',MarkerSize);hold on;
    end
    %Form patch structure
    Patch(i).kv=kv;Patch(i).wts=wts;Patch(i).cpt=cpt;Patch(i).num_cpt=num_cpt;Patch(i).order=order;
end
key_patch=i+1;


% Form junctions
A=[];B=[];% end points of the segments junction wise.
jun2patch=[]; % patches at a junction.
for i=1:size(int_pt,1)
    temp_B=[];temp_A=[];key=1;
    for j=1:size(vx,2)
        if (int_pt(i,1)==vx(1,j))&&(int_pt(i,2)==vy(1,j)) 
            temp_B(key,:)=[Vx(2,j),Vy(2,j)];
            temp_A(key,:)=[Vx(1,j),Vy(1,j)];
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,1,2];
            key=key+1;
        end
        if (int_pt(i,1)==vx(2,j))&&(int_pt(i,2)==vy(2,j)) 
            temp_A(key,:)=[Vx(2,j),Vy(2,j)];
            temp_B(key,:)=[Vx(1,j),Vy(1,j)];
            n1=Patch(j).num_cpt;
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,size(n1,2),size(n1,2)-1];
            key=key+1;
        end
    end
%     A(1,2*(i-1)+1:2*(i-1)+2)=int_pt(i,:);B(1,2*(i-1)+1:2*(i-1)+2)=int_pt(i,:);
    A(1:3,2*(i-1)+1:2*(i-1)+2)=temp_A;B(1:3,2*(i-1)+1:2*(i-1)+2)=temp_B;
end

% Form arcs at junction, discretize them and form associated bending strip.
Patch_Bs=struct;key_patchbs=1;cyclic=[1,2;2,3;3,1];
for i=1:size(int_pt,1)
    %% Find additional 3 cpts that define arcs
    c=zeros(4,3);
    c(1:2,1)=A(1,2*(i-1)+1:2*(i-1)+2);c(1:2,2)=int_pt(i,:);c(1:2,3)=A(2,2*(i-1)+1:2*(i-1)+2);
    V1=B(1,2*(i-1)+1:2*(i-1)+2)-A(1,2*(i-1)+1:2*(i-1)+2);V2=B(2,2*(i-1)+1:2*(i-1)+2)-A(2,2*(i-1)+1:2*(i-1)+2);
    V1=V1/norm(V1); V2=V2/norm(V2);
    w=[1,((dot(V1,V2)+1)/2)^0.5,1];
    c(4,:)=w(1,:);kv=[0 0 0 1 1 1];
    arc=nrbmak(c,kv); 
    arc=nrbkntins(arc,0.5);
    coefs=arc.coefs;
    cpt=coefs(1:3,:);wts=coefs(4,:);
    %%%% 1st two cpt at junction %%%%
    cpt_arc(1:3,1:2)=cpt(1:3,2:3);
    wts_arc(:,1:2)=wts(:,2:3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c=zeros(4,3);
    c(1:2,1)=A(2,2*(i-1)+1:2*(i-1)+2);c(1:2,2)=int_pt(i,:);c(1:2,3)=A(3,2*(i-1)+1:2*(i-1)+2);
    V1=B(2,2*(i-1)+1:2*(i-1)+2)-A(2,2*(i-1)+1:2*(i-1)+2);V2=B(3,2*(i-1)+1:2*(i-1)+2)-A(3,2*(i-1)+1:2*(i-1)+2);
    V1=V1/norm(V1); V2=V2/norm(V2);
    w=[1,((dot(V1,V2)+1)/2)^0.5,1];
    c(4,:)=w(1,:);kv=[0 0 0 1 1 1];
    arc=nrbmak(c,kv); 
    arc=nrbkntins(arc,0.5);
    coefs=arc.coefs;
    cpt=coefs(1:3,:);wts=coefs(4,:);
    %%%% 3rd cpt at junction %%%%
    cpt_arc(1:3,3)=cpt(1:3,3);
    wts_arc(:,3)=wts(:,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Wts(1,key_cpt:key_cpt+2)=wts_arc(1,:);Cpt(1:3,key_cpt:key_cpt+2)=cpt_arc(:,:);new_num_cpt=[key_cpt:key_cpt+2];key_cpt=key_cpt+3;
    if fig=='on'
        plot(cpt_arc(1,:),cpt_arc(2,:),'k.','MarkerSize',MarkerSize);hold on
    end
    %% Define arcs
    for j=1:3
        patch_num=[jun2patch(cyclic(j,1),3*(i-1)+1),jun2patch(cyclic(j,2),3*(i-1)+1)];
        patch_end=[jun2patch(cyclic(j,1),3*(i-1)+2),jun2patch(cyclic(j,2),3*(i-1)+2)];
        n1=Patch(patch_num(1)).num_cpt;n2=Patch(patch_num(2)).num_cpt;
        num_cpt=[n1(patch_end(1)),new_num_cpt(cyclic(j,1)),new_num_cpt(cyclic(j,2)),n2(patch_end(2))];
        cpt=Cpt(:,num_cpt);wts=Wts(:,num_cpt);kv=[0 0 0 0.5 1 1 1];
        Patch(key_patch).kv=kv;Patch(key_patch).wts=wts;Patch(key_patch).cpt=cpt;Patch(key_patch).num_cpt=num_cpt;Patch(key_patch).order=2;
        key_patch=key_patch+1;
        % Define bending strips
        bs_end=jun2patch(cyclic(j,1),3*(i-1)+3);
        num_cpt_bs=[n1(bs_end),num_cpt(1),num_cpt(2)];
        kv_bs=[0 0 0 1 1 1];wts_bs=Wts(1,num_cpt_bs);cpt_bs=Cpt(:,num_cpt_bs);
        Patch_Bs(key_patchbs).kv=kv_bs;Patch_Bs(key_patchbs).wts=wts_bs;Patch_Bs(key_patchbs).cpt=cpt_bs;Patch_Bs(key_patchbs).num_cpt=num_cpt_bs;
        Patch_Bs(key_patchbs).order=2;key_patchbs=key_patchbs+1;
    end   
end
%% Form Reference triads for each patch and bending strip
A0=zeros(3,3*size(Patch,2));
for i=1:size(Patch,2)
    cpt=Patch(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    A0(2,3*(i-1)+1:3*(i-1)+3)=[0,0,1];
    A0(3,3*(i-1)+1:3*(i-1)+3)=cross(A0(1,3*(i-1)+1:3*(i-1)+3),A0(2,3*(i-1)+1:3*(i-1)+3));
end
A0_Bs=zeros(3,3*size(Patch_Bs,2));
for i=1:size(Patch_Bs,2)
    cpt=Patch_Bs(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0_Bs(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    A0_Bs(2,3*(i-1)+1:3*(i-1)+3)=[0,0,1];
    A0_Bs(3,3*(i-1)+1:3*(i-1)+3)=cross(A0_Bs(1,3*(i-1)+1:3*(i-1)+3),A0_Bs(2,3*(i-1)+1:3*(i-1)+3));
end
%% Plot arcs
if fig=='on'
    pts=20;
    % Find NURBS shape fun.
    order_arc=2;kv_arc=[0 0 0 0.5 1 1 1];
    for i=size(vx,2)+1:size(Patch,2)
        w_arc=Patch(i).wts;
        for k=1:4
            for j=1:pts+1
                xi=(1/pts)*(j-1);
                [N(j,k),dN(j,k)]=NURBSbasis(k,order_arc,xi,kv_arc,w_arc);
            end
        end
        P=Patch(i).cpt;
        P=P(1:3,:);
        L0=zeros(pts+1,3);
        for j=1:4
            L0(:,1)=N(:,j)*P(1,j)+L0(:,1);L0(:,2)=N(:,j)*P(2,j)+L0(:,2);L0(:,3)=N(:,j)*P(3,j)+L0(:,3);
        end
        plot3(L0(:,1),L0(:,2),L0(:,3),'g','LineWidth',LineWidth);hold on
    end
end
end










