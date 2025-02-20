function [ vector ] = tangent_d1r( a1,a11,a1r,a11r,norm_a1 )
%Derivative of tangent 
%% INPUT
% r = rth dof in a particular element.
%% Code
n=norm_a1;vector=a11r/n-dot_p(a1,a1r)*a11/n^3+3*dot_p(a1,a11)*a1*dot_p(a1,a1r)/n^5-(dot_p(a1,a11r)+dot_p(a1r,a11))*a1/n^3-dot_p(a1,a11)*a1r/n^3;
end

