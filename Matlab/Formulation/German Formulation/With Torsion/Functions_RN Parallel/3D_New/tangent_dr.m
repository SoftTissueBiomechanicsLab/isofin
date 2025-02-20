function [ vector ] = tangent_dr( a1,a1r,norm_a1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Code
n=norm_a1;vector=a1r/n-((dot_p(a1,a1r))*a1)/n^3;
end