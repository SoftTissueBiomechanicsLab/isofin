function [ K,R,Q,Reaction ] = Build_KR(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,A0,ngp,D_DOF,F_DOF)
%% K_patch: Gives Stiffness(Global) matrix of a Patch.

%% Output
% K_Patch: Tangent stiffness matrix for a patch of elements.
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
%% Modify Q as per Dirichlet's BC
NOS_Total=size(P,2);
for i=1:size(D_DOF)
    for j=1:4
        if D_DOF(i,j+1)~=Inf
        % Modify Q
           Q(j,D_DOF(i,1))=P(j,D_DOF(i,1))+D_DOF(i,j+1);
        end
    end
end
%% Form complete stiffness matrix and residue vector for a patch. 
% nos=size(P,2);
K=zeros(NOS_Total*4,NOS_Total*4);R=zeros(NOS_Total*4,1);% Only unconstrained dof are considered.
%% Apply residual forces before displacement/force increment.
parfor i=1:Nele % Parallel for loop here.
  K_temp=zeros(NOS_Total*4,NOS_Total*4);R_temp=zeros(NOS_Total*4,1);% Only unconstrained dof are considered.
  [ K_ele,R_ele ]=Element_KR(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
  for j=1:size(Ele2Cp(i,:),2)
      for j1=1:4
          r=4*(j-1)+j1;% rth dof (local).
          R_temp(DOF(Ele2Cp(i,j),j1),1)=R_temp(DOF(Ele2Cp(i,j),j1),1)+R_ele(r,1);
          for k=1:size(Ele2Cp(i,:),2)
             for k1=1:4
                 s=4*(k-1)+k1;% sth dof (local).
                 if r<=s
                    a=DOF(Ele2Cp(i,j),j1);b=DOF(Ele2Cp(i,k),k1);
                     K_temp(a,b)=K_temp(a,b)+K_ele(r,s);
                     if a~=b
                        K_temp(b,a)=K_temp(b,a)+K_ele(r,s);
                     end
                 end
             end
         end
      end
  end
  K=K+K_temp;R=R+R_temp;
end
Reaction=R;
%% Apply Neumann BC (Force BC)
for i=1:size(F_DOF)
    for j=1:4
        % Modify R
        key=DOF(F_DOF(i,1),j);
        R(key,1)= R(key,1)-F_DOF(i,j+1);
    end
end
%% Apply Dirichlet's BC (Displacement BC)
for i=1:size(D_DOF)
    for j=1:4
        if D_DOF(i,j+1)~=Inf
            key=DOF(D_DOF(i,1),j);
        % Modify K
           K(key,:)=0;K(key,key)=1;
        % Modify R
           R(key,:)=0;
        end
    end
end
end






