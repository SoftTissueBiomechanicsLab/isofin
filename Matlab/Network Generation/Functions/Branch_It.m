function [] = Branch_It(n,SegEnds,range_theta,range_phi,range_R,NumIter,key,filepath)
%% RECURSIVE FUNCTION TO FORM A TREE
%% Input
% n : number of children fibers
% SegEnds: Trunk segment ends.
% range_theta: range of angle that branch can make wrt the trunk.
% range_phi: range of angle that segment connecting end of the new branch
% and centre of the bottom face of the cone can make with the n2 direction.
% range_R: range of length of the branch.

%% Code
stop_growth=zeros(1,n);
%% Establish Local coordinate system.
% Unit vector in local 1 direction pointing towards progressive branching.
e1=(SegEnds(:,2)-SegEnds(:,1))/norm(SegEnds(:,2)-SegEnds(:,1)); 
ranVec= randn(3,1);% vector in any random direction.
while ranVec==e1
ranVec=randn(3,1); % vector in any random direction.
end
% Unit vector in local 2 and local 1 direction respectively.
e2=cross(e1,ranVec);e2=e2/norm(e2);e3=cross(e1,e2);e3=e3/norm(e3); 
%% Find end of n new branches/segments
% Find local coordinates
R=range_R(1)+rand(n,1)*(range_R(2)-range_R(1));
theta=range_theta(1)+rand(n,1)*(range_theta(2)-range_theta(1));
phi=range_phi(1)+rand(n,1)*(range_phi(2)-range_phi(1));
x(1,:)=R.*cos(theta); 
x(2,:)=R.*sin(theta).*sin(phi);
x(3,:)=R.*sin(theta).*cos(phi);

% Find global coordinates
% Global to Local transformation matrix.
T_G2L(1,:)=e1;T_G2L(2,:)=e2;T_G2L(3,:)=e3; 
% Local to Global transformation matrix.
X=zeros(3,n);
for i=1:n % Rotation
    X(:,i)=T_G2L\x(:,i);
end

for i=1:n % Translation
    X(1,i)=SegEnds(1,2)+X(1,i);X(2,i)=SegEnds(2,2)+X(2,i);X(3,i)=SegEnds(3,2)+X(3,i);
end
%% Check if segments within the box
for i=1:n
[value,FaceOut] = PtInBox(X(:,i));
if value == 0
    stop_growth(1,i)=1;
    X(:,i)=terminate(SegEnds(:,2),X(:,i),FaceOut); 
end
end
% Check if new segments smaller than required.
key_small=0;
for i=1:n % Check is new segments smaller than required.
    if norm(SegEnds(:,2)-X(:,i)) < range_R(1)
        key_small=1;
        break;
    else
        key_small=0; 
    end
end

%% Store coordinates of the junctions and segments.
if key==1
    Segment(1,1:3)=SegEnds(:,1);Segment(1,4:6)=SegEnds(:,2);
    Segment(2,1:3)=SegEnds(:,2);Segment(2,4:6)=X(:,1);
    Segment(3,1:3)=SegEnds(:,2);Segment(3,4:6)=X(:,2);
    Junction(1,1:3)=SegEnds(:,2);Junction(1,4:6)=SegEnds(:,1);
    Junction(1,7:9)=X(:,1);Junction(1,10:12)=X(:,2);
    save(filepath,'Junction','Segment')
else
    if key_small==0
        load(filepath,'Junction','Segment');
        s1=size(Junction,1);s2=size(Segment,1);
        Segment(s2+1,1:3)=SegEnds(:,2);Segment(s2+1,4:6)=X(:,1);
        Segment(s2+2,1:3)=SegEnds(:,2);Segment(s2+2,4:6)=X(:,2);
        Junction(s1+1,1:3)=SegEnds(:,2);Junction(s1+1,4:6)=SegEnds(:,1);
        Junction(s1+1,7:9)=X(:,1);Junction(s1+1,10:12)=X(:,2);
        save(filepath,'Junction','Segment')
    end
end
%% Call the same function at these ends.
if key==NumIter
    return 
else
    key=key+1;
    if key_small==0
        if stop_growth(1,1)==0
            SegEnds1(:,1)=SegEnds(:,2);SegEnds1(:,2)=X(:,1);
            Branch_It(n,SegEnds1,range_theta,range_phi,range_R,NumIter,key,filepath)
        end
        if stop_growth(1,2)==0
            SegEnds2(:,1)=SegEnds(:,2);SegEnds2(:,2)=X(:,2);
            Branch_It(n,SegEnds2,range_theta,range_phi,range_R,NumIter,key,filepath)
        end
    end
end
end

