function [ matrix ] = Lambda_d1( T,T1,t,t1 )

% if T=T0 and t=T then derivative of T0 equal to zeros i.e give T1= zero
% vector as input.

I=eye(3);
dTt=T(1)*t(1)+T(2)*t(2)+T(3)*t(3);
dTt1=T(1)*t1(1)+T(2)*t1(2)+T(3)*t1(3);
dT1t=T1(1)*t(1)+T1(2)*t(2)+T1(3)*t(3);
CTt=cross(T,t); CTt1=cross(T,t1); CT1t=cross(T1,t);
matrix=(dTt1+dT1t)*I+cross_vM((CTt1+CT1t),I)-((dTt1+dT1t)/(1+dTt)^2)*CTt*...
    transpose(CTt)+(1/(1+dTt))*...
    (CTt*transpose(CT1t+CTt1)+(CT1t+CTt1)*transpose(CTt));
end
