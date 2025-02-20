function [ Cpt,kv,w] = Network_Mesh( SL_Ends,Nele,p )
% Network_Mesh: Provides mesh for required network geometry.
% Details: Provides control points, know vectors and weights for each cable
% in the network. 
%% OUTPUT
% Cpt : Control points of all cables involved. 
% kv: knotvector of all cables
% w: weights associated with control points for all cables.
%% INPUT
% SL_Ends : Provides of ends of all the straight lines in the network.
% Nele : No of elements in each cable.
% p : order of NURBS.
%% CODE
ncable=size(Sl_Ends,1);
Cpt=zeros(4,(p+Nele)*ncable);
for i=1:ncable
    
end
end

