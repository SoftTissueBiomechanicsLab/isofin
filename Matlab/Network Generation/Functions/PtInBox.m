function [value,FaceOut] = PtInBox(X)
%% Function to check if the point lies in the box. (Box in positive Octant)
%% Code
% Out of face 1
FaceOut=[0];value=1;key=1;
if  X(1)>=1
   FaceOut(key)=1;value = 0;key=key+1;
end
% Out of face 2
if  X(1)<=0
   FaceOut(key)=2;value = 0;key=key+1;
end
% Out of face 3
if  X(2)>=1
   FaceOut(key)=3;value = 0;key=key+1;
end
% Out of face 4
if  X(2)<=0
   FaceOut(key)=4;value = 0;key=key+1;
end
% Out of face 5
if  X(3)>=1
   FaceOut(key)=5;value = 0;key=key+1;
end
% Out of face 6
if  X(3)<=0
   FaceOut(key)=6;value = 0;key=key+1;
end
end

