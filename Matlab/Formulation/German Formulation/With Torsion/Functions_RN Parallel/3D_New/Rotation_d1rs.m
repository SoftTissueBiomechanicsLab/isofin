function [ matrix ] = Rotation_d1rs( t,t1,tr,ts,t1r,t1s,trs,t1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);sn=sin(psi);cs=cos(psi);
% psi_rs=0;psi_1rs=0;
I=[1,0,0;0,1,0;0,0,1];sn=sin(psi);cs=cos(psi);psi_rs=0;psi_1rs=0;
% m1=(-psi_1rs*sn-psi_1r*psi_s*cs-psi_1s*psi_r*cs-psi_1*psi_rs*cs+psi_1*psi_r*psi_s*sn)*I;
% m2=(psi_1rs*cs-psi_1r*psi_s*sn-psi_1s*psi_r*sn-psi_1*psi_rs*sn-psi_1*psi_r*psi_s*cs)*cross_vM(t,I);
% m3=(psi_1s*cs-psi_1*psi_s*sn)*cross_vM(tr,I)+(psi_1r*cs-psi_1*psi_r*sn)*cross_vM(ts,I)+(psi_rs*cs-psi_r*psi_s*sn)*cross_vM(t1,I);
% m4=(psi_1*cs)*cross_vM(trs,I)+(psi_r*cs)*cross_vM(t1s,I)+(psi_s*cs)*cross_vM(t1r,I)+(sn)*cross_vM(t1rs,I);
% matrix=m1+m2+m3+m4;
matrix=(-psi_1rs*sn-psi_1r*psi_s*cs-psi_1s*psi_r*cs-psi_1*psi_rs*cs+psi_1*psi_r*psi_s*sn)*I+(psi_1rs*cs-psi_1r*psi_s*sn-psi_1s*psi_r*sn-psi_1*psi_rs*sn-psi_1*psi_r*psi_s*cs)*cross_vM(t,I)+(psi_1s*cs-psi_1*psi_s*sn)*cross_vM(tr,I)+(psi_1r*cs-psi_1*psi_r*sn)*cross_vM(ts,I)+(psi_rs*cs-psi_r*psi_s*sn)*cross_vM(t1,I)+(psi_1*cs)*cross_vM(trs,I)+(psi_r*cs)*cross_vM(t1s,I)+(psi_s*cs)*cross_vM(t1r,I)+(sn)*cross_vM(t1rs,I);
end
