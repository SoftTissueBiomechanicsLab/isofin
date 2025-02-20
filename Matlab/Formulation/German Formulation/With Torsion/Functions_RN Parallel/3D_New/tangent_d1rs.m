function [ vector ] = tangent_d1rs( a1,a11,a1r,a1s,a11r,a11s,norm_a1 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% n=norm_a1;
% v1=-a11r*dot_p(a1,a1s)/n^3-(dot_p(a1,a1r)*a11s+dot_p(a1s,a1r)*a11)/n^3+3*dot_p(a1,a1r)*a11*dot_p(a1,a1s)/n^5;
% v2=3*(dot_p(a1s,a1r)*dot_p(a1,a11)*a1+dot_p(a1,a1r)*(dot_p(a1s,a11)+dot_p(a1,a11s))*a1)/n^5;
% v3=3*dot_p(a1,a1r)*dot_p(a1,a11)*a1s/n^5-15*dot_p(a1,a1r)*dot_p(a1,a11)*a1*dot_p(a1,a1s)/n^7;
% v4=-(dot_p(a1s,a11r)+dot_p(a11s,a1r))*a1/n^3;
% v5=-(dot_p(a1,a11r)+dot_p(a11,a1r))*a1s/n^3+3*(dot_p(a1,a11r)+dot_p(a11,a1r))*a1*dot_p(a1,a1s)/n^5;
% v6=-(dot_p(a1s,a11)+dot_p(a1,a11s))*a1r/n^3+3*dot_p(a1,a11)*a1r*dot_p(a1,a1s)/n^5;
% vector=v1+v2+v3+v4+v5+v6;

n=norm_a1;vector=-a11r*dot_p(a1,a1s)/n^3-(dot_p(a1,a1r)*a11s+dot_p(a1s,a1r)*a11)/n^3+3*dot_p(a1,a1r)*a11*dot_p(a1,a1s)/n^5+3*(dot_p(a1s,a1r)*dot_p(a1,a11)*a1+dot_p(a1,a1r)*(dot_p(a1s,a11)+dot_p(a1,a11s))*a1)/n^5+3*dot_p(a1,a1r)*dot_p(a1,a11)*a1s/n^5-15*dot_p(a1,a1r)*dot_p(a1,a11)*a1*dot_p(a1,a1s)/n^7-(dot_p(a1s,a11r)+dot_p(a11s,a1r))*a1/n^3-(dot_p(a1,a11r)+dot_p(a11,a1r))*a1s/n^3+3*(dot_p(a1,a11r)+dot_p(a11,a1r))*a1*dot_p(a1,a1s)/n^5-(dot_p(a1s,a11)+dot_p(a1,a11s))*a1r/n^3+3*dot_p(a1,a11)*a1r*dot_p(a1,a1s)/n^5;
end

