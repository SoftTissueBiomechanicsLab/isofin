function [ scalar ] = dot_p( v1,v2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
scalar=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3);
% scalar=transpose(v1)*v2; % Slower
end

