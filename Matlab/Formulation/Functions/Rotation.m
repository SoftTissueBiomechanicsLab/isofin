function [ matrix ] = Rotation( t,psi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Code
% I=eye(3);
% matrix=I*cos(psi)+sin(psi)*cross_vM(t,I);% Rotation matrix.
I=[1,0,0;0,1,0;0,0,1];matrix=I*cos(psi)+sin(psi)*cross_vM(t,I);% Rotation matrix.
end

