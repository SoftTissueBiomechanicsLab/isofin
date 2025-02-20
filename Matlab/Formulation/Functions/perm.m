function [ value] = perm( i,j,k )
%% perm: permutation operator
%% output: 
% value: returns 1,0, or -1.
%% input:
% inputs are indices.
value= -1;
if i==j||j==k||k==i
    value=0;
end
if i==1&&j==2&&k==3||i==2&&j==3&&k==1||i==3&&j==1&&k==2
   value=1;
end
end

