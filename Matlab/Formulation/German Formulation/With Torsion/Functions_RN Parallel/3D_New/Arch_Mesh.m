function [ kv,cpt,wts,nos arch ] = Arch_Mesh( radius,center,angle,nele,order )
% Generates mesh for an arch.
%% OUTPUT
% kv : knot vector
% cpt : control Points
% wts : weights
% nos : number of shape functions
%% INPUT
% center : row vector containing coordinates of the centre of the arch. [x,y,z]
% radius : radius of the arch.
% angle : row vector containing start and end angle of the arch. [sang eang] 
% nele : number of elements required.
% order : order of NURBS shape functions required. 
%% Code
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
% Generating basic arch structure.
sang=angle(1);eang=angle(2);
arch = nrbcirc(radius,center,sang,eang); % single element order 2 NURBS
% arch = nrbcirc(); % single element order 2 NURBS
% degree elevation
ntimes=order-2; % # of times we need to elevate.
arch = nrbdegelev(arch, ntimes);
% knot insertion
kv_in=zeros(1,nele-1);
for i=1:nele-1
    kv_in(1,i)=i/nele;
end
arch = nrbkntins(arch,kv_in);
% extract all required parameters from the structure.
kv=arch.knots;cpt=arch.coefs(1:3,:);wts=arch.coefs(4,:);nos=arch.number;
end

