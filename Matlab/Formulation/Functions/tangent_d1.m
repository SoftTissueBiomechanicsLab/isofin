function [ vector ] = tangent_d1( a1,a11,norm_a1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% n=norm_a1;vector=a11/n-((dot_p(a1,a11))*a1)/n^3;
n=norm_a1;vector=a11/n-((a1(1)*a11(1)+a1(2)*a11(2)+a1(3)*a11(3))*a1)/n^3;
vec1=a11/n;vec2=-((a1(1)*a11(1)+a1(2)*a11(2)+a1(3)*a11(3))*a1)/n^3;
end