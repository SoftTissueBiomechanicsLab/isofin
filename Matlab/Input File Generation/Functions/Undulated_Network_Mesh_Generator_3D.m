function [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Undulated_Network_Mesh_Generator_3D( segment,Segment_Period,junction,Rotation,ele_size,order,period,Amp,trim_factor,fig )
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')
%% Description
% Network_Mesh_Generator : Generates mesh for 3D random network (Generate using single fractal tree).

%% OUTPUT
% Cpt:
% Wts:
% Patch:
% Patch_Bs:
% A0:
% A0_Bs:
%% INPUT
% segment
% Segment_Period
% junction
% Rotation
% ele_size
% order
% period : Row vector containing max, min and
% smooth_factor
% trim_factor
% fig

%% Code
%% generate mesh.
%% Trim segments at the interior points (Junctions).
MarkerSize=20;LineWidth=2;
L=zeros(size(segment,1),1);
for i=1:size(segment,1)
    x=segment(i,1)-segment(i,4);y=segment(i,2)-segment(i,5);z=segment(i,3)-segment(i,6);
    L(i)=norm([x,y,z]);
end
Len=min(L);% shortest segment.
SEGMENT=zeros(size(segment));key=1;
for i=1:size(segment,1) % for all segments
    x1=segment(i,1);  x2=segment(i,4);y1=segment(i,2);  y2=segment(i,5);z1=segment(i,3);  z2=segment(i,6);
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
        SEGMENT(key,1)= x1+trim_factor*Len*(x2-x1)/L(i);
        SEGMENT(key,2)= y1+trim_factor*Len*(y2-y1)/L(i);
        SEGMENT(key,3)= z1+trim_factor*Len*(z2-z1)/L(i);
    else       % don't trim the segment from the start
        SEGMENT(key,1)= x1;
        SEGMENT(key,2)= y1;
        SEGMENT(key,3)= z1;
    end
    if key2==1 % trim the segment from the end
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
% if fig=='on'
%     for i=1:size(SEGMENT,1)
%         x1=SEGMENT(i,1);  x2=SEGMENT(i,4);y1=SEGMENT(i,2);  y2=SEGMENT(i,5);z1=SEGMENT(i,3);  z2=SEGMENT(i,6);
%         plot3([x1,x2],[y1,y2],[z1,z2],'k-','LineWidth',1);hold on
%     end
% end


%% Discretize segments with NURBS and store information of each patch (Introduce sin undulation in local y and z direction).
Patch=struct;Patch_Bs=struct;key_patchbs=1;
Cpt=[];Wts=[];key_cpt=1;
for i=1:size(SEGMENT,1)
%% Find cpts for main sin curve part of the fiber in local coordinate system.
    % Find the no of elements on the fiber based on element size.
    p1=[SEGMENT(i,1);SEGMENT(i,2);SEGMENT(i,3)]; p2=[SEGMENT(i,4);SEGMENT(i,5);SEGMENT(i,6)];
    l=norm(p2-p1);
    if ele_size<=l
        nele=ceil(norm(p1-p2)/ele_size);
    else
        nele=ceil(norm(p1-p2)/ele_size)+1;% if length of the fiber smaller than requested element size then 2 elements generated on the fiber.
    end
    % Find points to be interpolated in local coordinate system.
    ny=period(Segment_Period(i,1));nz=period(Segment_Period(i,2));
    npt=order+nele;
    X=zeros(3,npt);X(1,1)=0;X(1,npt)=1;
    X(1,2:npt-1)=(1/(npt-1))*(1:npt-2);
    X(2,:)=Amp*sin((X(1,:))*(2*ny*pi));X(3,:)=Amp*sin((X(1,:))*(2*nz*pi));
    % Find cpts, wts, kv in local coordinate system. 
    der(:,1:2)=[1,0,0;1,0,0]';
    [kv,wts,cpt] = NURBS_Interpolate_Der(X,order,der);
    nos=size(cpt,2);
%% Find cpts in global coordinate system.
    cpt=cpt*l; % shrink the curve according to the length of the fiber.
    X=X*l;
    R_L2G=Rotation{i};
    for j=1:nos
        cpt(:,j)=R_L2G*cpt(:,j); % rotate the curve to align with global orientation of the curve.
        cpt(:,j)=cpt(:,j)+p1;% shift it to the start of the fiber in global triad.
    end
    for j=1:size(X,2)
        X(:,j)=R_L2G*X(:,j); % rotate the curve to align with global orientation of the curve.
        X(:,j)=X(:,j)+p1;% shift it to the start of the fiber in global triad.
    end
    cpt(4,:)=0; % torsion dof
    % Add cpts and weights of main curve of fiber to the respective global matrix.
    Wts(1,key_cpt:key_cpt+nos-1)=wts; % Add to global weights of the network.
    Cpt(:,key_cpt:key_cpt+nos-1)=cpt; % Add to global cpts of the network.
    num_cpt=[key_cpt:key_cpt+nos-1];key_cpt=key_cpt+nos;
%     plot nurbs shape of the complete fiber.
    if fig == 'on'
        N=[];dN=[];nos=size(cpt,2);pts=200;
        for k=1:nos
            for j=1:pts+1
                xi=(1/pts)*(j-1);
                [N(j,k),dN(j,k)]=NURBSbasis(k,order,xi,kv,wts);
            end
        end
        q=cpt(1:3,:);
        L0=zeros(pts+1,3);
        for j=1:nos
            L0(:,1)=N(:,j)*q(1,j)+L0(:,1);L0(:,2)=N(:,j)*q(2,j)+L0(:,2);L0(:,3)=N(:,j)*q(3,j)+L0(:,3);
        end
        plot3(X(1,:),X(2,:),X(3,:),'.k','MarkerSize',MarkerSize);hold on;
        plot3(cpt(1,:),cpt(2,:),cpt(3,:),'--.r','MarkerSize',MarkerSize,'LineWidth',LineWidth);hold on;
%         for j=1:nos
%              plot3(cpt(1,:),cpt(2,:),cpt(3,:),'.r','MarkerSize',MarkerSize);hold on;
%         end
        plot3(L0(:,1),L0(:,2),L0(:,3),'b','LineWidth',LineWidth);hold on;
    end
    %% Form patch structure
% Main patch and smoothening curves
    % Middle(Main curve) of the fiber.
    Patch(i).kv=kv;Patch(i).wts=wts;Patch(i).cpt=cpt;Patch(i).num_cpt=num_cpt;Patch(i).order=order;  
end
key_patch=i+1;

%% Form new junctions
A=[];B=[];% end points of the segments junction wise.
jun2patch=[]; % patches at a junction.
for i=1:size(junction,1)
    temp_B=[];temp_A=[];key=1;
    for j=1:size(segment,1)
        x1=segment(j,1);  x2=segment(j,4);y1=segment(j,2);  y2=segment(j,5);z1=segment(j,3);  z2=segment(j,6);
        X1=SEGMENT(j,1);  X2=SEGMENT(j,4);Y1=SEGMENT(j,2);  Y2=SEGMENT(j,5);Z1=SEGMENT(j,3);  Z2=SEGMENT(j,6);
        if x1==junction(i,1) && y1==junction(i,2) && z1==junction(i,3)
            temp_B(key,:)=[X2,Y2,Z2];
            temp_A(key,:)=[X1,Y1,Z1];
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,1,2];
%             n1=Patch(3*(j-1)+2).num_cpt;
%             jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[3*(j-1)+2,size(n1,2),size(n1,2)-1];
            key=key+1;
        end
        if x2==junction(i,1) && y2==junction(i,2) && z2==junction(i,3)
            temp_A(key,:)=[X2,Y2,Z2];
            temp_B(key,:)=[X1,Y1,Z1];
            n1=Patch(j).num_cpt;
            jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[j,size(n1,2),size(n1,2)-1];
%             n2=Patch(3*(j-1)+3).num_cpt;
%             jun2patch(key,3*(i-1)+1:3*(i-1)+3)=[3*(j-1)+3,size(n2,2),size(n2,2)-1];
            key=key+1;
        end
    end
%     A(1,2*(i-1)+1:2*(i-1)+2)=int_pt(i,:);B(1,2*(i-1)+1:2*(i-1)+2)=int_pt(i,:);
    A(1:3,3*(i-1)+1:3*(i-1)+3)=temp_A;B(1:3,3*(i-1)+1:3*(i-1)+3)=temp_B;
end

%% Form arcs at junction, discretize them and form associated bending strip.
Patch_Bs=struct;key_patchbs=1;cyclic=[1,2;2,3;3,1];
for i=1:size(junction,1)
    cpt_arc=[];wts_arc=[];
    %% Find additional 3 cpts that define arcs
    c=zeros(4,3);
    c(1:3,1)=A(1,3*(i-1)+1:3*(i-1)+3);c(1:3,2)=junction(i,1:3);c(1:3,3)=A(2,3*(i-1)+1:3*(i-1)+3);
    V1=B(1,3*(i-1)+1:3*(i-1)+3)-A(1,3*(i-1)+1:3*(i-1)+3);V2=B(2,3*(i-1)+1:3*(i-1)+3)-A(2,3*(i-1)+1:3*(i-1)+3);
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
    c(1:3,1)=A(2,3*(i-1)+1:3*(i-1)+3);c(1:3,2)=junction(i,1:3);c(1:3,3)=A(3,3*(i-1)+1:3*(i-1)+3);
    V1=B(2,3*(i-1)+1:3*(i-1)+3)-A(2,3*(i-1)+1:3*(i-1)+3);V2=B(3,3*(i-1)+1:3*(i-1)+3)-A(3,3*(i-1)+1:3*(i-1)+3);
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
%         plot3(cpt_arc(1,:),cpt_arc(2,:),cpt_arc(3,:),'k.','MarkerSize',MarkerSize);hold on
          for t=1:size(cpt_arc,2)
%               plot3(cpt_arc(1,t),cpt_arc(2,t),cpt_arc(3,t),'.k','MarkerSize',MarkerSize);hold on;
%               text(cpt_arc(1,t),cpt_arc(2,t),cpt_arc(3,t),num2str(new_num_cpt(1,t)));
          end
    end
    %% Define arcs
    for j=1:3
        patch_num=[jun2patch(cyclic(j,1),3*(i-1)+1),jun2patch(cyclic(j,2),3*(i-1)+1)];
        patch_end=[jun2patch(cyclic(j,1),3*(i-1)+2),jun2patch(cyclic(j,2),3*(i-1)+2)];
        n1=Patch(patch_num(1)).num_cpt;n2=Patch(patch_num(2)).num_cpt;
        num_cpt=[n1(patch_end(1)),new_num_cpt(cyclic(j,1)),new_num_cpt(cyclic(j,2)),n2(patch_end(2))];
        cpt=Cpt(:,num_cpt);wts=Wts(:,num_cpt);kv=[0 0 0 0.5 1 1 1];
%         plot3(cpt(1,:),cpt(2,:),cpt(3,:),'ok','MarkerSize',15);
        Patch(key_patch).kv=kv;Patch(key_patch).wts=wts;Patch(key_patch).cpt=cpt;Patch(key_patch).num_cpt=num_cpt;Patch(key_patch).order=2;
        key_patch=key_patch+1;
        % Define bending strips
        bs_end=jun2patch(cyclic(j,1),3*(i-1)+3);
        num_cpt_bs=[n1(bs_end),num_cpt(1),num_cpt(2)];
        kv_bs=[0 0 0 1 1 1];cpt_bs=Cpt(:,num_cpt_bs);
%         wts_bs=Wts(1,num_cpt_bs);
        wts_bs=ones(1,size(num_cpt_bs,2));
        Patch_Bs(key_patchbs).kv=kv_bs;Patch_Bs(key_patchbs).wts=wts_bs;Patch_Bs(key_patchbs).cpt=cpt_bs;Patch_Bs(key_patchbs).num_cpt=num_cpt_bs;
        Patch_Bs(key_patchbs).order=2;key_patchbs=key_patchbs+1;
%         plot3(cpt_bs(1,:),cpt_bs(2,:),cpt_bs(3,:),'.k','MarkerSize',15);
    end   
end
%% Form Reference triads for each patch and bending strip
A0=zeros(3,3*size(Patch,2));
for i=1:size(Patch,2)
    cpt=Patch(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    ran_vec=[0,1,0];
    cross_vec=cross(ran_vec,A0(1,3*(i-1)+1:3*(i-1)+3));
    while norm(cross_vec)==0
        ran_vec=randn(1,3);ran_vec=ran_vec/norm(ran_vec);
        cross_vec=cross(ran_vec,A0(1,3*(i-1)+1:3*(i-1)+3));
    end
    A0(2,3*(i-1)+1:3*(i-1)+3)=cross(A0(1,3*(i-1)+1:3*(i-1)+3),ran_vec);A0(2,3*(i-1)+1:3*(i-1)+3)=A0(2,3*(i-1)+1:3*(i-1)+3)/norm(A0(2,3*(i-1)+1:3*(i-1)+3));
    A0(3,3*(i-1)+1:3*(i-1)+3)=cross(A0(1,3*(i-1)+1:3*(i-1)+3),A0(2,3*(i-1)+1:3*(i-1)+3));
end
A0_Bs=zeros(3,3*size(Patch_Bs,2));
for i=1:size(Patch_Bs,2)
    cpt=Patch_Bs(i).cpt; 
    p1=cpt(1:3,1);p2=cpt(1:3,2);
    A0_Bs(1,3*(i-1)+1:3*(i-1)+3)= (p2-p1)/norm(p2-p1);
    ran_vec=[0,1,0];
    cross_vec=cross(ran_vec,A0_Bs(1,3*(i-1)+1:3*(i-1)+3));
    while norm(cross_vec)==0
        ran_vec=randn(1,3);ran_vec=ran_vec/norm(ran_vec);
        cross_vec=cross(ran_vec,A0_Bs(1,3*(i-1)+1:3*(i-1)+3));
    end
    A0_Bs(2,3*(i-1)+1:3*(i-1)+3)=cross(A0_Bs(1,3*(i-1)+1:3*(i-1)+3),ran_vec);A0_Bs(2,3*(i-1)+1:3*(i-1)+3)=A0_Bs(2,3*(i-1)+1:3*(i-1)+3)/norm(A0_Bs(2,3*(i-1)+1:3*(i-1)+3));
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
        N=[];dN=[];
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
        plot3(L0(:,1),L0(:,2),L0(:,3),'g-','LineWidth',LineWidth);hold on
    end
end
end










