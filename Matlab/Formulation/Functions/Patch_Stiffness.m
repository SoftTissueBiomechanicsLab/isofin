function [ K_Patch ] = Patch_Stiffness(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,NC,A0,ngp)
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
K_Patch=zeros(NC-1,NC-1);% Only unconstrained dof are considered.
for i=1:Nele 
      K_ele=Element_Stiffness(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
      for j=1:size(Ele2Cp(i,:),2)
          for k=1:size(Ele2Cp(i,:),2)
                 for t1=1:4
                     r=4*(j-1)+t1;% rth dof (local).
                     for t2=1:4
                         s=4*(k-1)+t2;% sth dof (local).
                         if r<=s
                             a=DOF(Ele2Cp(i,j),t1);b=DOF(Ele2Cp(i,k),t2);
                             if max(a,b)<NC
                               K_Patch(a,b)=K_Patch(a,b)+K_ele(r,s);
                               if a~=b
                                   K_Patch(b,a)=K_Patch(b,a)+K_ele(r,s);
                               end
                             end
                         end
                     end
                 end
          end
      end
end
% K_Patch_L=transpose(K_Patch);
% K_Patch(logical(eye(size(K_Patch,2))))=0;
% K_Patch_U=K_Patch;
% K_Patch=K_Patch_L+K_Patch_U;
end
