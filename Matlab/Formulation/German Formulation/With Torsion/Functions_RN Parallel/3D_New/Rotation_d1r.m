function [ matrix ] = Rotation_d1r( t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);sn=sin(psi);cs=cos(psi);CtI=cross_vM(t,I);
% matrix=psi_1r*(-sn*I+cs*CtI)-psi_1*psi_r*(cs*I+sn*CtI)+psi_1*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(t1,I)+sn*cross_vM(t1r,I);
I=[1,0,0;0,1,0;0,0,1];sn=sin(psi);cs=cos(psi);CtI=cross_vM(t,I);matrix=psi_1r*(-sn*I+cs*CtI)-psi_1*psi_r*(cs*I+sn*CtI)+psi_1*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(t1,I)+sn*cross_vM(t1r,I);
end
