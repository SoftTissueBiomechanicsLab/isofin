function [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_3D( segment,junction,ele_size,order,end_cpt_factor,trim_factor,fig )
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')
%% Description
% Generates mesh for 3D random network (Generate using single fractal tree).

%% Trim segments at the interior points (Junctions).
L=zeros(size(segment,1),1);
for i=1:size(segment,1)
    x=segment(i,1)-segment(i,4);
    y=segment(i,2)-segment(i,5);
    z=segment(i,3)-segment(i,6);
    L(i)=norm([x,y,z]);
end
Len=min(L);% shortest segment.
SEGMENT=[];
% Loop index
key=1;
for i=1:size(segment,1) % for all segments
    x1=segment(i,1); x2=segment(i,4);
    y1=segment(i,2); y2=segment(i,5);
    z1=segment(i,3); z2=segment(i,6);
    for j=1:size(junction,1) % Find if start is part of a junction.
        if x1==junction(j,1) && y1==junction(j,2) && z1==junction(j,3)
           key1=1;
           break
        else
           key1=0;
        end
    end
    for j=1:size(junction,1)% Find if end is part of a junction.
        if x2==junction(j,1) && y2==junction(j,2) && z2==junction(j,3)
           key2=1;
           break
        else
           key2=0;
        end
    end
    if key1==1 % trim the segment from the start
        % Trimming based on the length of the shortest segment in the network.
        SEGMENT(key,1)= x1+trim_factor*Len*(x2-x1)/L(i);
        SEGMENT(key,2)= y1+trim_factor*Len*(y2-y1)/L(i);
        SEGMENT(key,3)= z1+trim_factor*Len*(z2-z1)/L(i);
    else       % don't trim the segment from the start
        SEGMENT(key,1)= x1;
        SEGMENT(key,2)= y1;
        SEGMENT(key,3)= z1;
    end
    if key2==1 % trim the segment from the end
        % Trimming based on the length of the shortest segment in the network.
        SEGMENT(key,4)= x2+trim_factor*Len*(x1-x2)/L(i);
        SEGMENT(key,5)= y2+trim_factor*Len*(y1-y2)/L(i);
        SEGMENT(key,6)= z2+trim_factor*Len*(z1-z2)/L(i);
    else       % don't trim the segment from the end
        SEGMENT(key,4)= x2;
        SEGMENT(key,5)= y2;
        SEGMENT(key,6)= z2;
    end
    key=key+1;
end
if fig=='on'
    figure
    for i=1:size(SEGMENT,1)
        x1=SEGMENT(i,1);  x2=SEGMENT(i,4);y1=SEGMENT(i,2);  y2=SEGMENT(i,5);z1=SEGMENT(i,3);  z2=SEGMENT(i,6);
        plot3([x1,x2],[y1,y2],[z1,z2],'b-','LineWidth',1);hold on
    end
end


%% Discretize segments with NURBS and store information of each patch.
Patch=struct;
Cpt=[];key_cpt=1;
for i=1:size(SEGMENT,1)
    p1=[SEGMENT(i,1);SEGMENT(i,2);SEGMENT(i,3)]; 
    p2=[SEGMENT(i,4);SEGMENT(i,5);SEGMENT(i,6)];
    l=norm(p2-p1);
    if ele_size<=l
        nele=ceil(norm(p1-p2)/ele_size);
    else
        % if length of the fiber smaller than requested element size then 
        % 2 elements generated on the fiber.
        nele=ceil(norm(p1-p2)/ele_size)+1;
    end
    nos=order+nele;
    % Form knot vector.
    kv=zeros(1,nos+order+1);
    for j=size(kv,2)-order:size(kv,2)
        kv(1,j)=1;
    end
    for j=2:nele
        kv(1,order+j)=(j-1)/nele;
    end
    % Form weights
    wts=ones(1,nos);
    Wts(1,key_cpt:key_cpt+nos-1)=wts;
    % Form control points.
    cpt=zeros(4,nos);
    cpt(1:3,1)=p1;
    cpt(1:3,nos)=p2;
    % Control point next to the ends is added based on the length of the fiber.
    len=min(ele_size,l);
    cpt(1:3,2)=p1+(len*end_cpt_factor)*(p2-p1)/l;
    cpt(1:3,nos-1)=p2+(len*end_cpt_factor)*(p1-p2)/l;
    % With placing second last cpts (from both ends) close to the respective end cpts.
    for j=1:nos-4 
        cpt(1:3,2+j)=p1+(j/((nos+1)-4))*(p2-p1);
    end

    num_cpt=key_cpt:key_cpt+nos-1;
    Cpt(:,key_cpt:key_cpt+nos-1)=cpt;key_cpt=key_cpt+nos;
    % if fig=='on'
    %     plot3(cpt(1,:),cpt(2,:),cpt(3,:),'r.','MarkerSize',15);hold on;
    % end
    %Form patch structure
    Patch(i).kv=kv;Patch(i).wts=wts;Patch(i).cpt=cpt;Patch(i).num_cpt=num_cpt;Patch(i).order=order;
end
key_patch=i+1;


%% Form new junctions
A=[]; B=[];% end points of the segments junction wise.
jun2patch=[]; % patches at a junction.
for i=1:size(junction,1)
    temp_B=[]; temp_A=[]; key=1;
    for j=1:size(segment,1)
        x1=segment(j,1);  x2=segment(j,4);
        y1=segment(j,2);  y2=segment(j,5);
        z1=segment(j,3);  z2=segment(j,6);
        X1=SEGMENT(j,1);  X2=SEGMENT(j,4);
        Y1=SEGMENT(j,2);  Y2=SEGMENT(j,5);
        Z1=SEGMENT(j,3);  Z2=SEGMENT(j,6);
        if x1==junction(i,1) && y1==junction(i,2) && z1==junction(i,3)
            temp_B(key,:)=[X2,Y2,Z2];
            temp_A(key,:)=[X1,Y1,Z1];
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,1,2];
            key=key+1;
        end
        if x2==junction(i,1) && y2==junction(i,2) && z2==junction(i,3)
            temp_A(key,:)=[X2,Y2,Z2];
            temp_B(key,:)=[X1,Y1,Z1];
            n1=Patch(j).num_cpt;
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,size(n1,2),size(n1,2)-1];
            key=key+1;
        end
    end
    A(1:3,3*(i-1)+1:3*(i-1)+3)=temp_A;B(1:3,3*(i-1)+1:3*(i-1)+3)=temp_B;
end

%% Form arcs at junction, discretize them and form associated bending strip.
Patch_Bs=struct;key_patchbs=1;cyclic=[1,2;2,3;3,1];
for i=1:size(junction,1)
    %% Find additional 3 cpts that define arcs
    c=zeros(4,3);
    c(1:3,1)=A(1,3*(i-1)+1:3*(i-1)+3);
    c(1:3,2)=junction(i,1:3);
    c(1:3,3)=A(2,3*(i-1)+1:3*(i-1)+3);
    V1=B(1,3*(i-1)+1:3*(i-1)+3)-A(1,3*(i-1)+1:3*(i-1)+3);
    V2=B(2,3*(i-1)+1:3*(i-1)+3)-A(2,3*(i-1)+1:3*(i-1)+3);
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
    c(1:3,1)=A(2,3*(i-1)+1:3*(i-1)+3);
    c(1:3,2)=junction(i,1:3);
    c(1:3,3)=A(3,3*(i-1)+1:3*(i-1)+3);
    V1=B(2,3*(i-1)+1:3*(i-1)+3)-A(2,3*(i-1)+1:3*(i-1)+3);
    V2=B(3,3*(i-1)+1:3*(i-1)+3)-A(3,3*(i-1)+1:3*(i-1)+3);
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
    % if fig=='on'
    %     plot3(cpt_arc(1,:),cpt_arc(2,:),cpt_arc(3,:),'k.','MarkerSize',15);hold on
    % end
    %% Define arcs
    for j=1:3
        patch_num=[jun2patch(cyclic(j,1),3*(i-1)+1),jun2patch(cyclic(j,2),3*(i-1)+1)];
        patch_end=[jun2patch(cyclic(j,1),3*(i-1)+2),jun2patch(cyclic(j,2),3*(i-1)+2)];
        n1=Patch(patch_num(1)).num_cpt;n2=Patch(patch_num(2)).num_cpt;
        num_cpt=[n1(patch_end(1)),new_num_cpt(cyclic(j,1)),new_num_cpt(cyclic(j,2)),n2(patch_end(2))];
        cpt=Cpt(:,num_cpt);wts=Wts(:,num_cpt);kv=[0 0 0 0.5 1 1 1];
        Patch(key_patch).kv=kv;
        Patch(key_patch).wts=wts;
        Patch(key_patch).cpt=cpt;
        Patch(key_patch).num_cpt=num_cpt;
        Patch(key_patch).order=2;
        key_patch=key_patch+1;
        % Define bending strips
        bs_end=jun2patch(cyclic(j,1),3*(i-1)+3);
        num_cpt_bs=[n1(bs_end),num_cpt(1),num_cpt(2)];
        kv_bs=[0 0 0 1 1 1];cpt_bs=Cpt(:,num_cpt_bs);
%         wts_bs=Wts(1,num_cpt_bs);
        wts_bs=ones(1,size(num_cpt_bs,2));
        Patch_Bs(key_patchbs).kv=kv_bs;
        Patch_Bs(key_patchbs).wts=wts_bs;
        Patch_Bs(key_patchbs).cpt=cpt_bs;
        Patch_Bs(key_patchbs).num_cpt=num_cpt_bs;
        Patch_Bs(key_patchbs).order=2;
        key_patchbs=key_patchbs+1;
    end   
end
%% Form Reference triads for each patch and bending strip
% Unit vectors of the triad are stored in rows, thus we need to transpose them in Network_KR!
A0=zeros(3,3*size(Patch,2)); 
for i=1:size(Patch,2)
    cpt=Patch(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    ran_vec=[0,1,0];
    while ran_vec==A0(1,3*(i-1)+1:3*(i-1)+3)
        ran_vec=randn(1,3);ran_vec=ran_vec/norm(ran_vec);
    end
    A0(2,3*(i-1)+1:3*(i-1)+3)=cross(A0(1,3*(i-1)+1:3*(i-1)+3),ran_vec);
    A0(2,3*(i-1)+1:3*(i-1)+3)=A0(2,3*(i-1)+1:3*(i-1)+3)/norm(A0(2,3*(i-1)+1:3*(i-1)+3));
    A0(3,3*(i-1)+1:3*(i-1)+3)=cross(A0(1,3*(i-1)+1:3*(i-1)+3),A0(2,3*(i-1)+1:3*(i-1)+3));
end
A0_Bs=zeros(3,3*size(Patch_Bs,2));
for i=1:size(Patch_Bs,2)
    cpt=Patch_Bs(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0_Bs(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    ran_vec=[0,1,0];
    while ran_vec==A0_Bs(1,3*(i-1)+1:3*(i-1)+3)
        ran_vec=randn(1,3);ran_vec=ran_vec/norm(ran_vec);
    end
    A0_Bs(2,3*(i-1)+1:3*(i-1)+3)=cross(A0_Bs(1,3*(i-1)+1:3*(i-1)+3),ran_vec);
    A0_Bs(2,3*(i-1)+1:3*(i-1)+3)=A0_Bs(2,3*(i-1)+1:3*(i-1)+3)/norm(A0_Bs(2,3*(i-1)+1:3*(i-1)+3));
    A0_Bs(3,3*(i-1)+1:3*(i-1)+3)=cross(A0_Bs(1,3*(i-1)+1:3*(i-1)+3),A0_Bs(2,3*(i-1)+1:3*(i-1)+3));
end
%% Plot arcs
if fig=='on'
    title ('IGA Mesh')
    xlabel('X Axis');ylabel('Y Axis');zlabel('Z Axis');
    % Bottom horizontal edges
    plot3([1,1],[0,1],[0,0],'k','LineWidth',1);hold on;
    plot3([0,0],[0,1],[0,0],'k','LineWidth',1);hold on;
    plot3([0,1],[1,1],[0,0],'k','LineWidth',1);hold on;
    plot3([0,1],[0,0],[0,0],'k','LineWidth',1);hold on;
    % Top horizontal edges
    plot3([1,0],[1,1],[1,1],'k','LineWidth',1);hold on;
    plot3([0,1],[0,0],[1,1],'k','LineWidth',1);hold on;
    plot3([1,1],[0,1],[1,1],'k','LineWidth',1);hold on;
    plot3([0,0],[1,0],[1,1],'k','LineWidth',1);hold on;
    % Vetrical edges
    plot3([0,0],[0,0],[0,1],'k','LineWidth',1);hold on;
    plot3([1,1],[0,0],[0,1],'k','LineWidth',1);hold on;
    plot3([1,1],[1,1],[0,1],'k','LineWidth',1);hold on;
    plot3([0,0],[1,1],[0,1],'k','LineWidth',1);hold on;
    pts=20;
    % Find NURBS shape fun.
    order_arc=2;kv_arc=[0 0 0 0.5 1 1 1];
    for i=size(segment,1)+1:size(Patch,2)
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
        plot3(L0(:,1),L0(:,2),L0(:,3),'g-','LineWidth',1);hold on
    end
end
end










