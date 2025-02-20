function [kv,wts,cpt] = NURBS_Interpolate(X,order)
%% NOTE: Chord method for parametrization and averaging for knot vector are the best possible methods for interpolation.(NURBS Book)
%% Code.
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\C_files');
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
%% Find the parametrization using chord method.
num_cpt=size(X,2);dim=size(X,1);
s=zeros(1,num_cpt);
for i=2:num_cpt % each of point in s is the cumulative chord distance upto that point.
    s(i)=s(i-1)+norm([X(1,i)-X(1,i-1),X(2,i)-X(2,i-1)]); 
end
s=s/max(s);
%% Find knotvector using averaging method.
kv=zeros(1,2*(order+1)+num_cpt-order-1);
kv(1:order+1)=0;kv(size(kv,2)-order:size(kv,2))=1;
for i=2:num_cpt-(order)
    kv(i+order)=(1/order)*sum(s(i:i+order-1));
end
wts=ones(1,num_cpt); % We don't play with weigths! (Yes I know it sucks!).
%% Evaluate control points for the interpolation.
% Form the matrix to solve linear equations to get cpts.
A=zeros(num_cpt);
for i=1:num_cpt
    N=zeros(1,num_cpt);dN=N;
    for j=1:num_cpt
        [N(j),dN(j)]=NURBSbasis(j,order,s(i),kv,wts);% Finding jth NURBS basis fun at s(i) parameter.
    end
    A(i,:)=N;
end
% Find cpts
cpt=zeros(dim,num_cpt);
for i=1:dim
    cpt(i,:)=(A\X(i,:)')';
end
end













