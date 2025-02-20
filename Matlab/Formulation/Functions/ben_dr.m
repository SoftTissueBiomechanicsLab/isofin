function [ scalar ] = ben_dr( LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,a1,a1r,A0,alpha )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% LT0T=Lambda(T0,T);LTt=Lambda(T,t);RT=Rotation(T,PSI);Rt=Rotation(t,psi);
% LT0T1=Lambda_d1(T0,T01,T,T1);LTt1=Lambda_d1(T,T1,t,t1);RT1=Rotation_d1(T,T1,PSI,PSI_1);Rt1=Rotation_d1(t,t1,psi,psi_1);
% LTtr=Lambda_dr(T,t,tr);Rtr=Rotation_dr(t,tr,psi,psi_r);
% LTt1r=Lambda_d1r(T,T1,t,t1,tr,t1r);Rt1r=Rotation_d1r(t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r);
v1=(Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T+Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T+Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T+Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1)*A0(:,alpha);v2=(Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,alpha);

% scalar=dot_p(v1,a1)+dot_p(v2,a1r);
scalar=v1(1)*a1(1)+v1(2)*a1(2)+v1(3)*a1(3)+v2(1)*a1r(1)+v2(2)*a1r(2)+v2(3)*a1r(3);
end
