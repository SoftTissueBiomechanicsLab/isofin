function [ Value ] = Is_PosDef( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Value=1;
a=eig(A);
if isreal(a) ==0
    disp('Eigen values of the matrix are not real!');Value=0;
else
    for i=1:length(a) 
       if a(i)<0
           Value=0;
       end
    end
end
end


