function [ kv,Cpt,w,nos ] = Mesh_Refinement( kv,Cpt,w,old_deg,new_deg,nele )
% Mesh_Refinement : Manipulates degree of shape functions and number of elements.
%% OUTPUT
%% INPUT
%% Code
addpath('C:\Users\smm5969\Desktop\NURBS_MESH');
% Create structure for unrefined mesh.
coefs(1:3,:)=Cpt;coefs(4,:)=w;
nurbs=nrbmak(coefs,kv);
% NOTE: Degree elevation and knot insertion do not commute. Here we first elevate the degree then insert knots.
% Degree elevation.
ntimes=new_deg-old_deg;
nurbs=nrbdegelev(nurbs, ntimes);
% Knot insertion.
kv_in=zeros(1,nele-1);
for i=1:nele-1
    kv_in(1,i)=i/nele;
end
nurbs = nrbkntins(nurbs,kv_in);
% extract all required parameters from the structure.
kv=nurbs.knots;Cpt=nurbs.coefs(1:3,:);w=nurbs.coefs(4,:);nos=nurbs.number;
end

