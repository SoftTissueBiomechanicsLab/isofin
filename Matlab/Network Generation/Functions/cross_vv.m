function [V] = cross_vv(v,u)
%% Calculates cross product of two 3D vectors, in the following order (v x u).
%% Code
V(1,1)=v(2)*u(3)-v(3)*u(2);V(2,1)=v(3)*u(1)-v(1)*u(3);V(3,1)=v(1)*u(2)-v(2)*u(1);
end

