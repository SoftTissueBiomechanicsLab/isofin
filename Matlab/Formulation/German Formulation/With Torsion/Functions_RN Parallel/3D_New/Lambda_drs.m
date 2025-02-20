function [ matrix ] = Lambda_drs( T,t,tr,ts,trs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);
I=[1,0,0;0,1,0;0,0,1];dTt=T(1)*t(1)+T(2)*t(2)+T(3)*t(3);CTt=cross(T,t);%Identity matrix  Scalar and a Column vector.

% m1=dot_p(T,trs)*I+cross_vM(cross(T,trs),I)+(2*dot_p(T,tr)*dot_p(T,ts)/(1+dTt)^3-dot_p(T,trs)/(1+dTt)^2)*CTt*transpose(CTt);
% m2=-(dot_p(T,tr)/(1+dTt)^2)*(CTt*transpose(cross(T,ts))+cross(T,ts)*transpose(CTt));
% m3=-(dot_p(T,ts)/(1+dTt)^2)*(CTt*transpose(cross(T,tr))+cross(T,tr)*transpose(CTt));
% m4=(1/(1+dTt))*(CTt*transpose(cross(T,trs))+cross(T,ts)*transpose(cross(T,tr))+cross(T,tr)*transpose(cross(T,ts))+cross(T,trs)*transpose(CTt));
% matrix=m1+m2+m3+m4;
matrix=dot_p(T,trs)*I+cross_vM(cross(T,trs),I)+(2*dot_p(T,tr)*dot_p(T,ts)/(1+dTt)^3-dot_p(T,trs)/(1+dTt)^2)*CTt*transpose(CTt)-(dot_p(T,tr)/(1+dTt)^2)*(CTt*transpose(cross(T,ts))+cross(T,ts)*transpose(CTt))-(dot_p(T,ts)/(1+dTt)^2)*(CTt*transpose(cross(T,tr))+cross(T,tr)*transpose(CTt))+(1/(1+dTt))*(CTt*transpose(cross(T,trs))+cross(T,ts)*transpose(cross(T,tr))+cross(T,tr)*transpose(cross(T,ts))+cross(T,trs)*transpose(CTt));
end
