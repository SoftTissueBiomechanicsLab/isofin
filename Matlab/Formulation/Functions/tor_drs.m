function [ scalar ] = tor_drs( LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,A0,beta,alpha )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% LT0T=Lambda(T0,T);LTt=Lambda(T,t);RT=Rotation(T,PSI);Rt=Rotation(t,psi);
% LT0T1=Lambda_d1(T0,T01,T,T1);LTt1=Lambda_d1(T,T1,t,t1);RT1=Rotation_d1(T,T1,PSI,PSI_1);Rt1=Rotation_d1(t,t1,psi,psi_1);
% LTtr=Lambda_dr(T,t,tr);Rtr=Rotation_dr(t,tr,psi,psi_r);LTts=Lambda_dr(T,t,ts);Rts=Rotation_dr(t,ts,psi,psi_s);
% LTt1r=Lambda_d1r(T,T1,t,t1,tr,t1r);Rt1r=Rotation_d1r(t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r);
% LTt1s=Lambda_d1r(T,T1,t,t1,ts,t1s);Rt1s=Rotation_d1r(t,t1,ts,t1s,psi,psi_1,psi_s,psi_1s);
% LTtrs=Lambda_drs(T,t,tr,ts,trs);Rtrs=Rotation_drs(t,tr,ts,trs,psi,psi_r,psi_s);
% LTt1rs=Lambda_d1rs(T,T1,t,t1,tr,ts,trs,t1r,t1s,t1rs);Rt1rs=Rotation_d1rs(t,t1,tr,ts,t1r,t1s,trs,t1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s);

% m1=Rt1rs*LTt*RT*LT0T+Rt1s*LTtr*RT*LT0T+Rt1r*LTts*RT*LT0T+Rt1*LTtrs*RT*LT0T;
% m2=Rtrs*LTt1*RT*LT0T+Rts*LTt1r*RT*LT0T+Rtr*LTt1s*RT*LT0T+Rt*LTt1rs*RT*LT0T;
% m3=Rtrs*LTt*RT1*LT0T+Rts*LTtr*RT1*LT0T+Rtr*LTts*RT1*LT0T+Rt*LTtrs*RT1*LT0T;
% m4=Rtrs*LTt*RT*LT0T1+Rts*LTtr*RT*LT0T1+Rtr*LTts*RT*LT0T1+Rt*LTtrs*RT*LT0T1;
% v1=(m1+m2+m3+m4)*A0(:,beta);
u1=(Rt*LTt*RT*LT0T)*A0(:,alpha);v1=(Rt1rs*LTt*RT*LT0T+Rt1s*LTtr*RT*LT0T+Rt1r*LTts*RT*LT0T+Rt1*LTtrs*RT*LT0T+Rtrs*LTt1*RT*LT0T+Rts*LTt1r*RT*LT0T+Rtr*LTt1s*RT*LT0T+Rt*LTt1rs*RT*LT0T+Rtrs*LTt*RT1*LT0T+Rts*LTtr*RT1*LT0T+Rtr*LTts*RT1*LT0T+Rt*LTtrs*RT1*LT0T+Rtrs*LTt*RT*LT0T1+Rts*LTtr*RT*LT0T1+Rtr*LTts*RT*LT0T1+Rt*LTtrs*RT*LT0T1)*A0(:,beta);

% m1=Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T;
% m2=Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T;
% m3=Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T;
% m4=Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1;
% v2=(m1+m2+m3+m4)*A0(:,beta);
u2=(Rts*LTt*RT*LT0T+Rt*LTts*RT*LT0T)*A0(:,alpha);v2=(Rt1r*LTt*RT*LT0T+Rt1*LTtr*RT*LT0T+Rtr*LTt1*RT*LT0T+Rt*LTt1r*RT*LT0T+Rtr*LTt*RT1*LT0T+Rt*LTtr*RT1*LT0T+Rtr*LTt*RT*LT0T1+Rt*LTtr*RT*LT0T1)*A0(:,beta);

% m1=Rt1s*LTt*RT*LT0T+Rts*LTt1*RT*LT0T+Rt1*LTts*RT*LT0T+Rt*LTt1s*RT*LT0T;
% m2=Rts*LTt*RT1*LT0T+Rts*LTt*RT*LT0T1+Rt*LTts*RT1*LT0T+Rt*LTts*RT*LT0T1;
% v3=(m1+m2)*A0(:,beta);
u3=(Rtr*LTt*RT*LT0T+Rt*LTtr*RT*LT0T)*A0(:,alpha);v3=(Rt1s*LTt*RT*LT0T+Rts*LTt1*RT*LT0T+Rt1*LTts*RT*LT0T+Rt*LTt1s*RT*LT0T+Rts*LTt*RT1*LT0T+Rts*LTt*RT*LT0T1+Rt*LTts*RT1*LT0T+Rt*LTts*RT*LT0T1)*A0(:,beta);

% m1=Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T;
% m2=Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1;
% v4=(m1+m2)*A0(:,beta);
u4=(Rtrs*LTt*RT*LT0T+Rts*LTtr*RT*LT0T+Rtr*LTts*RT*LT0T+Rt*LTtrs*RT*LT0T)*A0(:,alpha);v4=(Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,beta);

% scalar=dot_p(v1,u1)+dot_p(v2,u2)+dot_p(v3,u3)+dot_p(v4,u4);

scalar=v1(1)*u1(1)+v1(2)*u1(2)+v1(3)*u1(3)+v2(1)*u2(1)+v2(2)*u2(2)+v2(3)*u2(3)+v3(1)*u3(1)+v3(2)*u3(2)+v3(3)*u3(3)+v4(1)*u4(1)+v4(2)*u4(2)+v4(3)*u4(3);
end













