function [ K_nc,R_nc,t,f ] = Build_KR_nc(Q,Ele2Cp,num_ele,num_cp,kv,order,w,s,mag,angle,DOF)
%% Build_KR_nc: Gives Stiffness(Global) matrix component due to point follower load at the end of the cantilevered beam.
% Beam under goes fluttering response due to the follower tip load.


%% Output
% K_nc: Tangent stiffness matrix contribution of non-conservative forces.
% R_nc: Residue vector contribution of non-conservative forces.
%% Input
% NOS: Unused variable!
% Nele: No of elements in a patch.
% Dim: Dimension of the problem.
% P0: cpt matrix containing x-y-z of all cpts of undeformed configuration. Each column of P0 has one cpt 
% P: cpt matrix containing x-y-z of all cpts at some newton iteration. Each column of P has one cpt
% End_pt: Matrix conatining end pts of all the elements in the patch. Each row corresponds to an element. 
% Ele2Cp: Matrix containing indices of control pts. Each row corresponds to an element.
% DOF: Matrix containing indices of dofs of cpts. Each row corresponds to a cpt.
%% Code
%% Form complete stiffness matrix and residue vector for a patch. 
ncp=size(Q,2);npt=size(s,2); % number of points where point load is applied.
K_nc=zeros(ncp*4);R_nc=zeros(ncp*4,1);% Only unconstrained dof are considered.
%% CODE
%%  Calculate tangent vector at the parametric point where load is applied.
[N,dN1,dN2]=NURBS1DBasis2ndDers(s,order,kv,w);
nos=size(N,2);
for j=1:size(Ele2Cp,2)
Q_ele(:,j)=Q(:,Ele2Cp(num_ele,j));
end
a1=zeros(3,1);
for m=1:nos
  a1=a1+dN1(m)*Q_ele(1:3,m);
end
t=a1/norm(a1);
% Calculate the direction of the force.
Rot=[cos(angle),-sin(angle);sin(angle),cos(angle)];% Rotation matrix.
f=Rot*t(1:2,1)*mag; % direction of load.

%% Calculate residue due to the non-conservative load.
R_nc(DOF(num_cp,1):DOF(num_cp,2),1)=f;
force_mag=zeros(order+1,1);
force_mag(order+1,1)=mag;

%% Calculate tangent stiffness matrix due to the non-conservative load.
% Form K_temp and find contribution to K_nc.
K_temp=zeros(nos*4);
nor_a1=norm(a1);
for m=1:nos
    for m1=1:2
        for n=1:nos
            for n1=1:2
                a_temp=Rot*a1(1:2,:);
                K_temp((m-1)*4+m1,(n-1)*4+n1)=force_mag(m)*((dN1(n)*Rot(m1,n1)/nor_a1)-a_temp(m1,1)*dN1(n)*a1(n1,1)/nor_a1^3);
                K_nc(DOF(Ele2Cp(num_ele,m),m1),DOF(Ele2Cp(num_ele,n),n1))=K_nc(DOF(Ele2Cp(num_ele,m),m1),DOF(Ele2Cp(num_ele,n),n1))+K_temp((m-1)*4+m1,(n-1)*4+n1);
            end
        end
    end
end
end






