function [ matrix ] = Rotation_drs( t,tr,ts,trs,psi,psi_r,psi_s )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% psi_rs=0;
% I=eye(3);sn=sin(psi);cs=cos(psi);CtI=cross_vM(t,I);
% matrix=psi_rs*(-sn*I+cs*CtI)-psi_s*psi_r*(cs*I+sn*CtI)+psi_s*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(ts,I)+sn*cross_vM(trs,I);

psi_rs=0;I=[1,0,0;0,1,0;0,0,1];sn=sin(psi);cs=cos(psi);CtI=cross_vM(t,I);matrix=psi_rs*(-sn*I+cs*CtI)-psi_s*psi_r*(cs*I+sn*CtI)+psi_s*cs*cross_vM(tr,I)+psi_r*cs*cross_vM(ts,I)+sn*cross_vM(trs,I);
end

