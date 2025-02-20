function [ scalar ] = tor_dr( LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,A0,beta,alpha )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% LT0T=Lambda(T0,T);LTt=Lambda(T,t);RT=Rotation(T,PSI);Rt=Rotation(t,psi);
% LT0T1=Lambda_d1(T0,T01,T,T1);LTt1=Lambda_d1(T,T1,t,t1);RT1=Rotation_d1(T,T1,PSI,PSI_1);Rt1=Rotation_d1(t,t1,psi,psi_1);
% LTtr=Lambda_dr(T,t,tr);Rtr=Rotation_dr(t,tr,psi,psi_r);
% LTt1r=Lambda_d1r(T,T1,t,t1,tr,t1r);Rt1r=Rotation_d1r(t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r);

% m1=Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T;
% m2=Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T;
% m3=Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T;
% m4=Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1;
% v1=(m1+m2+m3+m4)*A0(:,beta);
u1=(Rt*LTt*RT*LT0T)*A0(:,alpha);v1=(Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T+Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T+Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T+Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1)*A0(:,beta);


% m1=Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T;
% m2=Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1;
% v2=(m1+m2)*A0(:,beta);
u2=(Rtr*LTt*RT*LT0T+Rt*LTtr*RT*LT0T)*A0(:,alpha);v2=(Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,beta);

% scalar=dot_p(v1,u1)+dot_p(v2,u2);
scalar=v1(1)*u1(1)+v1(2)*u1(2)+v1(3)*u1(3)+v2(1)*u2(1)+v2(2)*u2(2)+v2(3)*u2(3);
end















