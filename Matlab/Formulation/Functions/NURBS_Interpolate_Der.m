function [kv,wts,cpt] = NURBS_Interpolate_Der(X,order,der)
%% NOTE: Chord method for parametrization and averaging for knot vector are the best possible methods for interpolation.(NURBS Book)
%% This function interpolates curves at given points with perscibed derivates on both ends of the curve.
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
m=num_cpt+order+2; %with two end derivatives
kv=zeros(1,m+1);kv(1:order+1)=0;kv(m+1-order:m+1)=1;
for i=1:num_cpt-(order)+1
    kv(i+order+1)=(1/order)*sum(s(i:i+order-1));
end
wts=ones(1,num_cpt+2); % We don't play with weigths! (Yes I know it sucks!).
%% Evaluate control points for the interpolation.
A=zeros(num_cpt+2);
for i=1:num_cpt
    N=zeros(1,num_cpt+2);dN=N;
    for j=1:num_cpt+2
        [N(j),dN(j)]=NURBSbasis(j,order,s(i),kv,wts);% Finding jth NURBS basis fun at s(i) parameter.
    end
    A(i,:)=N;
end
% Modify A to incorporate derivative requirements at ends
der=der*1;
% Insert equation in second row
A(3:size(A,1)-1,:)=A(2:size(A,1)-2,:);A(2,:)=0;A(2,1)=-1;A(2,2)=1;
X_new=zeros(dim,num_cpt+2);X_new(1:dim,1)=X(1:dim,1);
X_new(1:dim,3:num_cpt+1)=X(1:dim,2:num_cpt);
X_new(1:dim,2)=(kv(order+2)/order)*der(:,1);
% Insert equation in last row
A(num_cpt+2,num_cpt+1)=-1;A(num_cpt+2,num_cpt+2)=1;
X_new(1:dim,num_cpt+2)=((1-kv(m-order-1))/order)*der(:,1);
% Find cpts
cpt=zeros(dim,num_cpt+2);
for i=1:dim
    cpt(i,:)=(A\X_new(i,:)')';
end
end













