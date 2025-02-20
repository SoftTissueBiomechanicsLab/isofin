function [ vector ] = tangent_drs( a1,a1r,a1s,norm_a1 )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%% Code
n=norm_a1;vector=-a1r*dot_p(a1,a1s)/n^3-(a1s*dot_p(a1,a1r)+a1*dot_p(a1s,a1r))/n^3+3*dot_p(a1,a1r)*a1*dot_p(a1,a1s)/n^5;
end

