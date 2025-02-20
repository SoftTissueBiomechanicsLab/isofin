function [ matrix ] = cross_vM( v,M )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Code
% matrix(1,1)=v(2)*M(3,1)-v(3)*M(2,1);
% matrix(1,2)=v(2)*M(3,2)-v(3)*M(2,2);
% matrix(1,3)=v(2)*M(3,3)-v(3)*M(2,3);
% matrix(2,1)=v(3)*M(1,1)-v(1)*M(3,1);
% matrix(2,2)=v(3)*M(1,2)-v(1)*M(3,2);
% matrix(2,3)=v(3)*M(1,3)-v(1)*M(3,3);
% matrix(3,1)=v(1)*M(2,1)-v(2)*M(1,1);
% matrix(3,2)=v(1)*M(2,2)-v(2)*M(1,2);
% matrix(3,3)=v(1)*M(2,3)-v(2)*M(1,3);

matrix(1,1)=v(2)*M(3,1)-v(3)*M(2,1);matrix(1,2)=v(2)*M(3,2)-v(3)*M(2,2);matrix(1,3)=v(2)*M(3,3)-v(3)*M(2,3);matrix(2,1)=v(3)*M(1,1)-v(1)*M(3,1);matrix(2,2)=v(3)*M(1,2)-v(1)*M(3,2);matrix(2,3)=v(3)*M(1,3)-v(1)*M(3,3);matrix(3,1)=v(1)*M(2,1)-v(2)*M(1,1);matrix(3,2)=v(1)*M(2,2)-v(2)*M(1,2);matrix(3,3)=v(1)*M(2,3)-v(2)*M(1,3);

end


