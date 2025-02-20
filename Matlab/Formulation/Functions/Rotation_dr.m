function [ matrix ] = Rotation_dr( t,tr,psi,psi_r )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);sn=sin(psi);
% matrix=psi_r*(-sn*I+cos(psi)*cross_vM(t,I))+sn*cross_vM(tr,I);

I=[1,0,0;0,1,0;0,0,1];sn=sin(psi);matrix=psi_r*(-sn*I+cos(psi)*cross_vM(t,I))+sn*cross_vM(tr,I);
end

