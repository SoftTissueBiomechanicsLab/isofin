function [ matrix ] = Lambda_d1rs( T,T1,t,t1,tr,ts,t1r,t1s,trs,t1rs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);
I=[1,0,0;0,1,0;0,0,1];

%% Scalars
% dTt=dot_p(T,t);% Scalar
% dT1t=dot_p(T1,t);dTt1=dot_p(T,t1);
% dTts=dot_p(T,ts);dTtr=dot_p(T,tr);dTtrs=dot_p(T,trs);
% dTt1s=dot_p(T,t1s);dTt1r=dot_p(T,t1r);dTt1rs=dot_p(T,t1rs);
% dT1ts=dot_p(T1,ts);dT1tr=dot_p(T1,tr);dT1trs=dot_p(T1,trs);
dTt=T(1)*t(1)+T(2)*t(2)+T(3)*t(3);dT1t=T1(1)*t(1)+T1(2)*t(2)+T1(3)*t(3);dTt1=T(1)*t1(1)+T(2)*t1(2)+T(3)*t1(3);dTts=T(1)*ts(1)+T(2)*ts(2)+T(3)*ts(3);dTtr=T(1)*tr(1)+T(2)*tr(2)+T(3)*tr(3);dTtrs=T(1)*trs(1)+T(2)*trs(2)+T(3)*trs(3);dTt1s=T(1)*t1s(1)+T(2)*t1s(2)+T(3)*t1s(3);dTt1r=T(1)*t1r(1)+T(2)*t1r(2)+T(3)*t1r(3);dTt1rs=T(1)*t1rs(1)+T(2)*t1rs(2)+T(3)*t1rs(3);dT1ts=T1(1)*ts(1)+T1(2)*ts(2)+T1(3)*ts(3);dT1tr=T1(1)*tr(1)+T1(2)*tr(2)+T1(3)*tr(3);dT1trs=T1(1)*trs(1)+T1(2)*trs(2)+T1(3)*trs(3);

%% Column Vetors
% CTt=cross(T,t);% Column vector.
% CT1t=cross(T1,t);CTt1=cross(T,t1);
% CTts=cross(T,ts);CTtr=cross(T,tr);CTtrs=cross(T,trs);
% CTt1s=cross(T,t1s);CTt1r=cross(T,t1r);CTt1rs=cross(T,t1rs);
% CT1ts=cross(T1,ts);CT1tr=cross(T1,tr);CT1trs=cross(T1,trs);
CTt=cross(T,t);CT1t=cross(T1,t);CTt1=cross(T,t1);CTts=cross(T,ts);CTtr=cross(T,tr);CTtrs=cross(T,trs);CTt1s=cross(T,t1s);CTt1r=cross(T,t1r);CTt1rs=cross(T,t1rs);CT1ts=cross(T1,ts);CT1tr=cross(T1,tr);CT1trs=cross(T1,trs);

%% Matrices
% m1=(dT1trs+dTt1rs)*I+cross_vM((CT1trs+CTt1rs),I);
% m2=-6*(((dT1t+dTt1)*dTtr*dTts)/(1+dTt)^4)*CTt*transpose(CTt);
% m3=((2*(dT1ts+dTt1s)*dTtr+2*(dT1t+dTt1)*dTtrs)/(1+dTt)^3)*CTt*transpose(CTt);
% m4=(2*(dT1t+dTt1)*dTtr/(1+dTt)^3)*(CTt*transpose(CTts)+CTts*transpose(CTt));
% m5=(2*(dT1tr+dTt1r)*dTts/(1+dTt)^3-(dT1trs+dTt1rs)/(1+dTt)^2)*CTt*transpose(CTt);
% m6=-((dT1tr+dTt1r)/(1+dTt)^2)*(CTt*transpose(CTts)+CTts*transpose(CTt));
% m7=(2*(dT1t+dTt1)*dTts/(1+dTt)^3-(dT1ts+dTt1s)/(1+dTt)^2)*(CTt*transpose(CTtr)+CTtr*transpose(CTt));
% m8=-((dT1t+dTt1)/(1+dTt)^2)*(CTt*transpose(CTtrs)+CTts*transpose(CTtr)+CTtr*transpose(CTts)+CTtrs*transpose(CTt));
% m9=(2*dTtr*dTts/(1+dTt)^3-dTtrs/(1+dTt)^2)*(CTt*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTt));
% m10=-(dTtr/(1+dTt)^2)*(CTt*transpose(CT1ts+CTt1s)+CTts*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTts)+(CT1ts+CTt1s)*transpose(CTt));
% m11=-(dTts/(1+dTt)^2)*(CTt*transpose(CT1tr+CTt1r)+CTtr*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTtr)+(CT1tr+CTt1r)*transpose(CTt));
% m12=(1/(1+dTt))*(CTt*transpose(CT1trs+CTt1rs)+CTts*transpose(CT1tr+CTt1r)+CTtr*transpose(CT1ts+CTt1s)+CTtrs*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTtrs)+(CT1ts+CTt1s)*transpose(CTtr)+(CT1tr+CTt1r)*transpose(CTts)+(CT1trs+CTt1rs)*transpose(CTt));
% matrix=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12;

matrix=(dT1trs+dTt1rs)*I+cross_vM((CT1trs+CTt1rs),I)-6*(((dT1t+dTt1)*dTtr*dTts)/(1+dTt)^4)*CTt*transpose(CTt)+((2*(dT1ts+dTt1s)*dTtr+2*(dT1t+dTt1)*dTtrs)/(1+dTt)^3)*CTt*transpose(CTt)+(2*(dT1t+dTt1)*dTtr/(1+dTt)^3)*(CTt*transpose(CTts)+CTts*transpose(CTt))+(2*(dT1tr+dTt1r)*dTts/(1+dTt)^3-(dT1trs+dTt1rs)/(1+dTt)^2)*CTt*transpose(CTt)-((dT1tr+dTt1r)/(1+dTt)^2)*(CTt*transpose(CTts)+CTts*transpose(CTt))+(2*(dT1t+dTt1)*dTts/(1+dTt)^3-(dT1ts+dTt1s)/(1+dTt)^2)*(CTt*transpose(CTtr)+CTtr*transpose(CTt))-((dT1t+dTt1)/(1+dTt)^2)*(CTt*transpose(CTtrs)+CTts*transpose(CTtr)+CTtr*transpose(CTts)+CTtrs*transpose(CTt))+(2*dTtr*dTts/(1+dTt)^3-dTtrs/(1+dTt)^2)*(CTt*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTt))-(dTtr/(1+dTt)^2)*(CTt*transpose(CT1ts+CTt1s)+CTts*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTts)+(CT1ts+CTt1s)*transpose(CTt))-(dTts/(1+dTt)^2)*(CTt*transpose(CT1tr+CTt1r)+CTtr*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTtr)+(CT1tr+CTt1r)*transpose(CTt))+(1/(1+dTt))*(CTt*transpose(CT1trs+CTt1rs)+CTts*transpose(CT1tr+CTt1r)+CTtr*transpose(CT1ts+CTt1s)+CTtrs*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTtrs)+(CT1ts+CTt1s)*transpose(CTtr)+(CT1tr+CTt1r)*transpose(CTts)+(CT1trs+CTt1rs)*transpose(CTt));
end
