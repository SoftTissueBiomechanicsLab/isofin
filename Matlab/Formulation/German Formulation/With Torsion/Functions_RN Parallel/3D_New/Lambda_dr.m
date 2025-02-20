function [ matrix ] = Lambda_dr( T,t,tr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);
% dTt=dot_p(T,t);% Scalar
% dTtr=dot_p(T,tr);
% CTt=cross(T,t);% Column vector.
% CTtr=cross(T,tr);
% matrix=dTtr*I+cross_vM(CTtr,I)-(dTtr/(1+dTt)^2)*CTt*transpose(CTt)+(1/(1+dTt))*(CTt*transpose(CTtr)+CTtr*transpose(CTt));

I=[1,0,0;0,1,0;0,0,1];dTt=T(1)*t(1)+T(2)*t(2)+T(3)*t(3);dTtr=T(1)*tr(1)+T(2)*tr(2)+T(3)*tr(3);CTt=cross(T,t);CTtr=cross(T,tr);matrix=dTtr*I+cross_vM(CTtr,I)-(dTtr/(1+dTt)^2)*CTt*transpose(CTt)+(1/(1+dTt))*(CTt*transpose(CTtr)+CTtr*transpose(CTt));

end

