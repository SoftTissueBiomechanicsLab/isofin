function [ matrix ] = Lambda_d1r( T,T1,t,t1,tr,t1r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);
I=[1,0,0;0,1,0;0,0,1];dTt=T(1)*t(1)+T(2)*t(2)+T(3)*t(3);CTt=cross(T,t);%Identity matrix Scalar and a Column vector.

% m1=(dot_p(T1,tr)+dot_p(T,t1r))*I+cross_vM((cross(T1,tr)+cross(T,t1r)),I);
% m2=2*(((dot_p(T1,t)+dot_p(T,t1))*dot_p(T,tr))/(1+dTt)^3)*CTt*transpose(CTt);
% m3=-((dot_p(T1,tr)+dot_p(T,t1r))/(1+dTt)^2)*CTt*transpose(CTt);
% m4=-((dot_p(T1,t)+dot_p(T,t1))/(1+dTt)^2)*(CTt*transpose(cross(T,tr))+cross(T,tr)*transpose(CTt));
% m5=-(dot_p(T,tr)/(1+dTt)^2)*(CTt*transpose(cross(T1,t)+cross(T,t1))+(cross(T1,t)+cross(T,t1))*transpose(CTt));
% m61=CTt*transpose(cross(T1,tr)+cross(T,t1r))+cross(T,tr)*transpose(cross(T1,t)+cross(T,t1));
% m62=(cross(T1,t)+cross(T,t1))*transpose(cross(T,tr))+(cross(T1,tr)+cross(T,t1r))*transpose(CTt);
% m6=(1/(1+dTt))*(m61+m62);
% matrix=m1+m2+m3+m4+m5+m6;
matrix=(dot_p(T1,tr)+dot_p(T,t1r))*I+cross_vM((cross(T1,tr)+cross(T,t1r)),I)+2*(((dot_p(T1,t)+dot_p(T,t1))*dot_p(T,tr))/(1+dTt)^3)*CTt*transpose(CTt)-((dot_p(T1,tr)+dot_p(T,t1r))/(1+dTt)^2)*CTt*transpose(CTt)-((dot_p(T1,t)+dot_p(T,t1))/(1+dTt)^2)*(CTt*transpose(cross(T,tr))+cross(T,tr)*transpose(CTt))-(dot_p(T,tr)/(1+dTt)^2)*(CTt*transpose(cross(T1,t)+cross(T,t1))+(cross(T1,t)+cross(T,t1))*transpose(CTt))+(1/(1+dTt))*(CTt*transpose(cross(T1,tr)+cross(T,t1r))+cross(T,tr)*transpose(cross(T1,t)+cross(T,t1))+(cross(T1,t)+cross(T,t1))*transpose(cross(T,tr))+(cross(T1,tr)+cross(T,t1r))*transpose(CTt));
end
