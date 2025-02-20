function [value,face] = Isotropy_Criteria(Segment,NumFiber)
%% Checks Isotropy Criteria for the network.
%% INPUT
% Segment: Contains coordinates of the segments.
% tol: tolerence for acceptance of variance in no of fiber ends on the
% faces.
%% OUTPUT
% Criterion 1: No of fiber ends on all the faces of the cube don't vary much.
value = 0;
face=zeros(1,6); % stores no of fiber ends on the faces.
for i=1:size(Segment,1)
    for j=1:2
        X=Segment(i,3*(j-1)+1:3*(j-1)+3);
        if X(1)==1 % fiber on face 1
            face(1,1)=face(1,1)+1;
        end
        if X(1)==0 % fiber on face 2
            face(1,2)=face(1,2)+1;
        end
        if X(2)==1 % fiber on face 3
            face(1,3)=face(1,3)+1;
        end
        if X(2)==0 % fiber on face 4
            face(1,4)=face(1,4)+1;
        end
        if X(3)==1 % fiber on face 5
            face(1,5)=face(1,5)+1;
        end
        if X(3)==0 % fiber on face 6
            face(1,6)=face(1,6)+1;
        end
    end
end
if max(face)<NumFiber(2) && min(face)>NumFiber(1)
   value=1;
end
end

