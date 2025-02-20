function [ matrix ] = Rotation_d1(t,t1,psi,psi_1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Code
% sn=sin(psi);I=eye(3);
% matrix=psi_1*(-sn*I+cos(psi)*cross_vM(t,I))+sn*cross_vM(t1,I);  
sn=sin(psi);I=[1,0,0;0,1,0;0,0,1];matrix=psi_1*(-sn*I+ cos(psi)*cross_vM(t,I)) + sn*cross_vM(t1,I); 
end
