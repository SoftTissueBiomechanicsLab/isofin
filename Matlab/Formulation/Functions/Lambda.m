function [ matrix ] = Lambda( N0,N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Code
% d=dot_p(N0,N);% Scalar
% C=cross(N0,N);% Column vector.
% matrix=d*eye(3)+cross_vM(C,eye(3))+(1/(1+d))*(C*transpose(C));% Transformation matrix.
I=[1,0,0;0,1,0;0,0,1];d=N0(1)*N(1)+N0(2)*N(2)+N0(3)*N(3);C=cross(N0,N);matrix=d*I+cross_vM(C,I)+(1/(1+d))*(C*transpose(C));% Transformation matrix.
end

