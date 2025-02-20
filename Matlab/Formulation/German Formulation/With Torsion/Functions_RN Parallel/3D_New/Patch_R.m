function [ Reaction ] = Patch_R(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,A0,ngp)
%% Description: Script use to calculate reaction of a patch.
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
% Note: Constrained dofs are numbered highest.
% NC: Smallest index of constrained dof.
%% Code

NOS_Total=size(P,2);
Reaction=zeros(NOS_Total*4,1);% Only unconstrained dof are considered.
parfor i=1:Nele % parfor
  [ R_ele ]=Element_R(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
  R_temp=zeros(NOS_Total*4,1);
  for j=1:size(Ele2Cp(i,:),2)
      for j1=1:4
          r=4*(j-1)+j1;% rth dof (local).
          R_temp(DOF(Ele2Cp(i,j),j1),1)=R_temp(DOF(Ele2Cp(i,j),j1),1)+R_ele(r,1);
      end
  end
  Reaction=Reaction+R_temp;
end
end
