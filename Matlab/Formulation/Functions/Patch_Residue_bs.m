function [ F_Patch ] = Patch_Residue_bs(F_Patch,Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,NC,A0,ngp)
%% R_Patch: Gives Residual(Global) vector of a patch.
%% Output
% R_Patch: Residue Vector for a patch of elements.
%% Input
% NOS:
% Nele: No of elements in a patch.
% Dim: Dimension of the problem.
% P0: cpt matrix containing x-y-z of all cpts of undeformed configuration. Each column of P0 has one cpt 
% P: cpt matrix containing x-y-z of all cpts at some newton iteration. Each column of P has one cpt
% End_pt: Matrix conatining end pts of all the elements in the patch. Each row corresponds to an element. 
% Ele2Cp: Matrix containing indices of control pts. Each row corresponds to an element.
% DOF: Matrix containing indices of dofs of cpts. Each row corresponds to a cpt.
% Note: Constrained dofs are numbered highest.
% NC: Smallest index of constrained dof.
for i=1:Nele
      F_ele=Element_Residue_bs(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
      for j=1:size(Ele2Cp(i,:),2)
          for j1=1:4
              r=4*(j-1)+j1;% rth dof (local).
              if DOF(Ele2Cp(i,j),j1)<NC
                F_Patch(DOF(Ele2Cp(i,j),j1),1)=F_Patch(DOF(Ele2Cp(i,j),j1),1)+F_ele(r,1);
              end
          end                
      end
end
end

