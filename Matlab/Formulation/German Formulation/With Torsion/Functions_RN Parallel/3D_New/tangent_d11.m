function [ vector ] = tangent_d11( a1,a11,a111,norm_a1 )
%UNTITLED5 Summary of this function goes here
%% INPUT

%% Code
n=norm_a1;
vector=a111/n-dot_p(a1,a11)*a11/n^3+3*(dot_p(a1,a11)^2)*a1/n^5-(dot_p(a1,a111)+dot_p(a11,a11))*a1/n^3-dot_p(a1,a11)*a11/n^3;
end

