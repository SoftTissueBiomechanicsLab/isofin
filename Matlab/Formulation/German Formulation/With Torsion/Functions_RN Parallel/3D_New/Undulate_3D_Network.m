function [] = Undulate_3D_Network(network_path,new_network_path,period,select_period)
%% Undulates existing 3D straight fiber networks.
%% Input
%% Output
%% Code
load(network_path) % load mat file containing straight fiber 3D network.
num_segments=size(Segment,1);
Segment_Period=zeros(num_segments,2);
Rotation={};Length=zeros(num_segments,1);
for i=1:num_segments
%% Assign periods of sin function in local y and z.
    % NOTE: The longer the fiber the more it has choice to select the
    % periods.(This is done as we control the mesh size globally and not decide the number of elements per fiber)
    % Find choices the fiber has based on it's length.
    n1=(Segment(i,4:6)-Segment(i,1:3))'; 
    Length(i) = norm(n1);% length of the fiber
    key=1;
    for j=1:size(select_period,2)-1
        if Length(i) > select_period(j) && Length(i) < select_period(j+1)
            break;
        else
            key=key+1;
        end      
    end
    Segment_Period(i,:)=ceil(key*rand(1,2));
%% Assign local orientation of the fiber wrt global coordinate system.
   n1=n1/Length(i); % unit vector in the direction of the fiber.(% 1st local direction. )
   ranvec=randn(3,1);ranvec=ranvec/norm(ranvec);
   while ranvec == n1
       ranvec=randn(3,1);ranvec=ranvec/norm(ranvec);
   end
   n2=cross(n1,ranvec);n2=n2/norm(n2);n3=cross(n1,n2);% 2nd and 3rd in local direction.
   R_L2G(:,1)=n1;R_L2G(:,2)=n2;R_L2G(:,3)=n3; %Rotation Matrix Local to Global.
   Rotation{i}=R_L2G;
end
save(new_network_path,'Junction','Segment','Segment_Period','Rotation','period')
end















