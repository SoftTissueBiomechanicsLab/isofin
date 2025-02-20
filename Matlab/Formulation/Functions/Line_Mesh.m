function [ kv,cpt,wts,nos line ] = Line_Mesh( p1,p2,nele,order )
% Generates mesh for an arch.
%% OUTPUT
% kv : knot vector
% cpt : control Points
% wts : weights
% nos : number of shape functions
%% INPUT
% p1 : Start point of the line.
% p2 : End point of the line.
% nele : number of elements required.
% order : order of NURBS shape functions required. 
%% Code
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
% Generating basic line structure.
line = nrbline(p1,p2); % single element order 1 NURBS
% degree elevation
ntimes=order-1; % # of times we need to elevate.
line = nrbdegelev(line, ntimes);
% knot insertion
kv_in=zeros(1,nele-1);
for i=1:nele-1
    kv_in(1,i)=i/nele;
end
line = nrbkntins(line,kv_in);
% extract all required parameters from the structure.
kv=line.knots;cpt=line.coefs(1:3,:);wts=line.coefs(4,:);nos=line.number;
end

