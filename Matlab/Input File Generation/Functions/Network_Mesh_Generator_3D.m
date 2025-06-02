function [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_3D( segment,junction,ele_size,order,end_cpt_factor,trim_factor,fig )
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')
%% Description
% Generates mesh for 3D random network (Generate using single fractal tree).

%% Trim segments at the interior points (Junctions).
fiber_length=zeros(size(segment,1),1);
for ii=1:size(segment,1)
    x=segment(ii,1)-segment(ii,4);
    y=segment(ii,2)-segment(ii,5);
    z=segment(ii,3)-segment(ii,6);
    fiber_length(ii)=norm([x,y,z]);
end
min_fiber_length=min(fiber_length);% shortest segment.

% Initialize new segment matrix as zero matrix size of original segment
SEGMENT=zeros(size(segment));

% Loop through each fiber in segment
for ii=1:size(segment,1) % for all segments
    % Get coordinates of each fiber node
    x1=segment(ii,1); x2=segment(ii,4);
    y1=segment(ii,2); y2=segment(ii,5);
    z1=segment(ii,3); z2=segment(ii,6);
    % Loop through each row of junction and determine if fiber node is in
    % junction (i.e., if multiple fibers intercept at the node in segment)
    % Do this loop for first node (x1,...) and second node (x2,...)
    for j=1:size(junction,1) % Find if start is part of a junction.
        if x1==junction(j,1) && y1==junction(j,2) && z1==junction(j,3)
           % Key for if first fiber node is a junction
           key1=1;
           break
        else
           key1=0;
        end
    end
    for j=1:size(junction,1)% Find if end is part of a junction.
        if x2==junction(j,1) && y2==junction(j,2) && z2==junction(j,3)
           % Key for if second fiber node is a junction
           key2=1;
           break
        else
           key2=0;
        end
    end
    if key1==1 % trim the segment from the start
        % Trimming based on the length of the shortest segment in the network.
        SEGMENT(ii,1)= x1+trim_factor*min_fiber_length*(x2-x1)/fiber_length(ii);
        SEGMENT(ii,2)= y1+trim_factor*min_fiber_length*(y2-y1)/fiber_length(ii);
        SEGMENT(ii,3)= z1+trim_factor*min_fiber_length*(z2-z1)/fiber_length(ii);
    else       % don't trim the segment from the start
        SEGMENT(ii,1)= x1;
        SEGMENT(ii,2)= y1;
        SEGMENT(ii,3)= z1;
    end
    if key2==1 % trim the segment from the end
        % Trimming based on the length of the shortest segment in the network.
        SEGMENT(ii,4)= x2+trim_factor*min_fiber_length*(x1-x2)/fiber_length(ii);
        SEGMENT(ii,5)= y2+trim_factor*min_fiber_length*(y1-y2)/fiber_length(ii);
        SEGMENT(ii,6)= z2+trim_factor*min_fiber_length*(z1-z2)/fiber_length(ii);
    else       % don't trim the segment from the end
        SEGMENT(ii,4)= x2;
        SEGMENT(ii,5)= y2;
        SEGMENT(ii,6)= z2;
    end

end


%% Discretize segments with NURBS and store information of each patch.
Patch=struct;
Cpt=[];key_cpt=1;
for ii=1:size(SEGMENT,1)
    % Get first and second node of fiber i
    p1=[SEGMENT(ii,1);SEGMENT(ii,2);SEGMENT(ii,3)]; 
    p2=[SEGMENT(ii,4);SEGMENT(ii,5);SEGMENT(ii,6)];
    temp_fiber_length=norm(p2-p1);
    if ele_size<=temp_fiber_length
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
    len=min(ele_size,temp_fiber_length);
    cpt(1:3,2)=p1+(len*end_cpt_factor)*(p2-p1)/temp_fiber_length;
    cpt(1:3,nos-1)=p2+(len*end_cpt_factor)*(p1-p2)/temp_fiber_length;
    % With placing second last cpts (from both ends) close to the respective end cpts.
    for j=1:nos-4 
        cpt(1:3,2+j)=p1+(j/((nos+1)-4))*(p2-p1);
    end

    num_cpt=key_cpt:key_cpt+nos-1;
    Cpt(:,key_cpt:key_cpt+nos-1)=cpt;key_cpt=key_cpt+nos;
    % if fig=='on'
    % plot3(cpt(1,:),cpt(2,:),cpt(3,:),'.','MarkerSize',15);hold on;
    % end
    %Form patch structure
    Patch(ii).kv=kv;
    Patch(ii).wts=wts;
    Patch(ii).cpt=cpt;
    Patch(ii).num_cpt=num_cpt;
    Patch(ii).order=order;
end
key_patch=ii+1;


%% Form new junctions
A=[]; B=[];% end points of the segments junction wise.
jun2patch=[]; % patches at a junction.
for ii=1:size(junction,1)
    temp_B=[]; temp_A=[]; key=1;
    for j=1:size(segment,1)
        % Get coordinates for original and new fiber segments
        x1=segment(j,1);  x2=segment(j,4);
        y1=segment(j,2);  y2=segment(j,5);
        z1=segment(j,3);  z2=segment(j,6);
        X1=SEGMENT(j,1);  X2=SEGMENT(j,4);
        Y1=SEGMENT(j,2);  Y2=SEGMENT(j,5);
        Z1=SEGMENT(j,3);  Z2=SEGMENT(j,6);
        % Check if original fiber nodes are junctions. If so, store the new
        % nodal coordinates in temporary arrays. temp_A stores the new
        % nodes that originally were at junction  i. temp_B stores the
        % other end of the fiber
        if x1==junction(ii,1) && y1==junction(ii,2) && z1==junction(ii,3)            
            temp_A(key,:)=[X1,Y1,Z1];
            temp_B(key,:)=[X2,Y2,Z2];
            % If first node of fiber is a junction, save fiber number along
            % with 1 and 2 to index the control points 
            jun2patch(key,3*(ii-1)+(1:3))=[j,1,2];
            key=key+1;
        end
        if x2==junction(ii,1) && y2==junction(ii,2) && z2==junction(ii,3)
            temp_A(key,:)=[X2,Y2,Z2];
            temp_B(key,:)=[X1,Y1,Z1];
            n1=Patch(j).num_cpt;
            % If second node of fiber is a junction, save fiber number along
            % with number of control points (and # cpts-1) to index the 
            % control points 
            jun2patch(key,3*(ii-1)+(1:3))=[j,size(n1,2),size(n1,2)-1];
            key=key+1;
        end
    end
   
    % A is an arrays that store the new coordinates of fiber ends which
    % originally met at a given junction. B stores the coordinates of the
    % other end of the fiber corresponding to A. If a network has a
    % different
    num_fibers = size(temp_A,1);
    A(1:num_fibers,3*(ii-1)+(1:3))=temp_A;
    B(1:num_fibers,3*(ii-1)+(1:3))=temp_B;
end

%% Form arcs at junction, discretize them and form associated bending strip.

% Generalized for arbitrary number of fibers meeting at each junction
Patch_Bs=struct; key_patchbs=1;
for ii = 1:size(junction,1)
    % Temporary variables to hold values of A and B to be manipulated in
    % correctly creating arcs
    temp_A = A(:,3*(ii-1)+(1:3));
    temp_B = B(:,3*(ii-1)+(1:3));    
    temp_jun2patch_1 = jun2patch(:,3*(ii-1)+1);
    temp_jun2patch_2 = jun2patch(:,3*(ii-1)+2);
    temp_jun2patch_3 = jun2patch(:,3*(ii-1)+3);


    % Determine how many rows are nonzero to determine number of fibers at
    % a junction. Dividing the sum by 3 and ceiling that value gives binary
    % values - 1 if that row corresponds to a fiber, 0 if it is empty. This
    % is my current workaround for dealing with junctions that have
    % different numbers of fibers. 
    nonzero_rows = ceil(sum(temp_A ~= [0 0 0],2)/3);
    num_fibers = sum(nonzero_rows);

    % Sort temp_A and temp_B so that nodes are ordered clockwise
    if num_fibers > 3
        % Be sure to only include points corresponding to fibers. If there
        % are varied numbers of fibers at each junction, there will be zero
        % rows in the variables, which we do not want to use
        [temp_A,sort_id] = Sort_Points(temp_A(1:num_fibers,:));
        temp_B = temp_B(sort_id,:);
        temp_jun2patch_1 = temp_jun2patch_1(sort_id);
        temp_jun2patch_2 = temp_jun2patch_2(sort_id);
        temp_jun2patch_3 =  temp_jun2patch_3(sort_id);

    end

    % Using num_fibers, create the cyclic_variable
    cyclic = [1:num_fibers; 2:num_fibers 1]';
    cpt_arc = zeros(3,num_fibers);
    wts_arc = zeros(1,num_fibers);
    for jj = 1:num_fibers-1
        % temp_cpt is the temporary control point variable to store the 
        % control points used in generating arcs between fibers
        temp_cpt = zeros(4,3);
        % First column corresponds to junction-end of one fiber
        temp_cpt(1:3,1) = temp_A(jj,:);%A(jj,3*(ii-1)+(1:3));
        % Second column corresponds to junction coordinates
        temp_cpt(1:3,2) = junction(ii,1:3);
        % Third column corresponds to junction-end of other fiber
        temp_cpt(1:3,3) = temp_A(jj+1,:);%A(jj+1,3*(ii-1)+(1:3));
        % V1 and V2 are used to determine the tangent vector of two fibers
        % meeting at a junction
        V1 = temp_B(jj,:) - temp_A(jj,:);%B(jj,3*(ii-1)+(1:3)) - A(jj,3*(ii-1)+(1:3));
        V2 = temp_B(jj+1,:) - temp_A(jj+1,:);%B(jj+1,3*(ii-1)+(1:3)) - A(jj+1,3*(ii-1)+(1:3));
        % Normalize vectors
        V1 = V1/norm(V1); V2 = V2/norm(V2);
        % Calculate weight for the NURBS and assign to 4th row of temp_cpt
        temp_cpt(4,:) = [1,((dot(V1,V2)+1)/2)^0.5,1];
        % Knot vector for the arc
        kv = [0 0 0 1 1 1];
        % Calculate coefficients for the arc to connect the two fibers
        arc=nrbmak(temp_cpt,kv); 
        % Insert a knot into the spline to create two control points along
        % the arc
        arc=nrbkntins(arc,0.5);
        coefs=arc.coefs; 
        % First three rows of coefs are control point coordinates, while 
        % the fourth row holds the control point weights
        cpt = coefs(1:3,:); wts = coefs(4,:);
        % If on first loop to create the arc, store columns 2 and 3 of cpt 
        % and w, otherwise, only store column 3 of cpt and w
        if jj==1
            cpt_arc(1:3,1:2) = cpt(:,2:3);
            wts_arc(1:2) = wts(2:3);
        else
            cpt_arc(1:3,jj+1) = cpt(1:3,3);
            wts_arc(jj+1) = wts(3);
        end
    end

    Wts(1,key_cpt:(key_cpt+num_fibers-1)) = wts_arc;
    Cpt(1:3,(key_cpt:key_cpt+num_fibers-1)) = cpt_arc;
    new_num_cpt = key_cpt:(key_cpt+num_fibers-1);
    key_cpt=key_cpt+num_fibers;

    for kk=1:num_fibers
        % Get patch number (corresponds to fiber number)
        % patch_num=[jun2patch(cyclic(kk,1),3*(ii-1)+1),jun2patch(cyclic(kk,2),3*(ii-1)+1)];
        patch_num=[temp_jun2patch_1(cyclic(kk,1)),temp_jun2patch_1(cyclic(kk,2))];

        % Get control points of end of patch near a junction
        % patch_end=[jun2patch(cyclic(kk,1),3*(ii-1)+2),jun2patch(cyclic(kk,2),3*(ii-1)+2)];
        patch_end=[temp_jun2patch_2(cyclic(kk,1)),temp_jun2patch_2(cyclic(kk,2))];

        % Get control points for the patches
        n1=Patch(patch_num(1)).num_cpt;
        n2=Patch(patch_num(2)).num_cpt;
        % Store the control points for the arc (to be created) in num_cpt.
        % The arc will start at the first/last control point of patch_num(1),
        % go through two of the new control points, and end at the
        % first/last control point of patch_num(2) (Depends on the
        % direction the patches/fibers run in).
        num_cpt=[n1(patch_end(1)),new_num_cpt(cyclic(kk,1)),new_num_cpt(cyclic(kk,2)),n2(patch_end(2))];
        % Get the control points and weights for the arcs
        cpt=Cpt(:,num_cpt);
        wts=Wts(:,num_cpt);
        % Create a new knot vector for this arc
        kv=[0 0 0 0.5 1 1 1];
        % Add new patches for the arcs
        Patch(key_patch).kv=kv;
        Patch(key_patch).wts=wts;
        Patch(key_patch).cpt=cpt;
        Patch(key_patch).num_cpt=num_cpt;
        Patch(key_patch).order=2;
        key_patch=key_patch+1;
        % Define bending strips - start by defining the control points for
        % the bending strip. bs_end gives one end of the bending strip
        % which connects the end of a fiber to the arc. This connection 
        % allows for bending moments to transfer across fibers
        % bs_end=jun2patch(cyclic(kk,1),3*(ii-1)+3);
        bs_end=temp_jun2patch_3(cyclic(kk,1));

        % Get control point numbers for the bending strip
        num_cpt_bs=[n1(bs_end),num_cpt(1),num_cpt(2)];
        % Define the knot vector for the bending strip
        kv_bs=[0 0 0 1 1 1];
        % Get the control points for the bending strip and add them to a
        % structure containing the bending strips
        cpt_bs=Cpt(:,num_cpt_bs);
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
for ii=1:size(Patch,2)
    % Get control points for Patch i
    cpt=Patch(ii).cpt; 
    % Get first two control points of Patch i, and use them to get the
    % vector tangent to the patch
    p1=cpt(1:3,1);
    p2=cpt(1:3,2);
    % First row of A0 holds the tangent vectors to the patches
    A0(1,3*(ii-1)+(1:3))= (p2-p1)/norm(p2-p1);
    % Initalize a random vector. Modify it if has the same direction as the
    % tangent to the patch
    ran_vec=[0,1,0];
    while ran_vec==A0(1,3*(ii-1)+(1:3))
        ran_vec=randn(1,3);
        ran_vec=ran_vec/norm(ran_vec);
    end
    % Find a vector normal to the tangent vecctor using the random vector
    A0(2,3*(ii-1)+(1:3))=cross(A0(1,3*(ii-1)+(1:3)),ran_vec);
    % Normalize new normal 
    A0(2,3*(ii-1)+(1:3))=A0(2,3*(ii-1)+(1:3))/norm(A0(2,3*(ii-1)+(1:3)));
    % Using 2 normalized vectors of the triad, find a third that is normal
    % to the plane of the 2 vectors.
    A0(3,3*(ii-1)+(1:3))=cross(A0(1,3*(ii-1)+(1:3)),A0(2,3*(ii-1)+(1:3)));
end

% Repeat above process for bending strips
A0_Bs=zeros(3,3*size(Patch_Bs,2));
for ii=1:size(Patch_Bs,2)
    cpt=Patch_Bs(ii).cpt; 
    p1=cpt(1:3,1);
    p2=cpt(1:3,2);
    A0_Bs(1,3*(ii-1)+(1:3))= (p2-p1)/norm(p2-p1);
    ran_vec=[0,1,0];
    while ran_vec==A0_Bs(1,3*(ii-1)+(1:3))
        ran_vec=randn(1,3);
        ran_vec=ran_vec/norm(ran_vec);
    end
    A0_Bs(2,3*(ii-1)+(1:3))=cross(A0_Bs(1,3*(ii-1)+(1:3)),ran_vec);
    A0_Bs(2,3*(ii-1)+(1:3))=A0_Bs(2,3*(ii-1)+(1:3))/norm(A0_Bs(2,3*(ii-1)+(1:3)));
    A0_Bs(3,3*(ii-1)+(1:3))=cross(A0_Bs(1,3*(ii-1)+(1:3)),A0_Bs(2,3*(ii-1)+(1:3)));
end
end










