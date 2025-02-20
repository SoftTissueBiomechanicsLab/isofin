function [ K_Patch,R_Patch ] = Patch_KR(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,NC,A0,ngp)
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

K_Patch=zeros(NC-1,NC-1);R_Patch=zeros(NC-1,1);% Only unconstrained dof are considered.
parfor i=1:Nele % parfor
  [ K_ele,R_ele ]=Element_KR(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
   K_temp=zeros(NC-1,NC-1);R_temp=zeros(NC-1,1);
  for j=1:size(Ele2Cp(i,:),2)
      for j1=1:4
          r=4*(j-1)+j1;% rth dof (local).
          if DOF(Ele2Cp(i,j),j1)<NC
               R_temp(DOF(Ele2Cp(i,j),j1),1)=R_temp(DOF(Ele2Cp(i,j),j1),1)+R_ele(r,1);
%              R_Patch(DOF(Ele2Cp(i,j),j1),1)=R_Patch(DOF(Ele2Cp(i,j),j1),1)+R_ele(r,1);
          end
          for k=1:size(Ele2Cp(i,:),2)
             for k1=1:4
                 s=4*(k-1)+k1;% sth dof (local).
                 if r<=s
                     a=DOF(Ele2Cp(i,j),j1);b=DOF(Ele2Cp(i,k),k1);
                     if max(a,b)<NC
                         K_temp(a,b)=K_temp(a,b)+K_ele(r,s);
%                        K_Patch(a,b)=K_Patch(a,b)+K_ele(r,s);
                       if a~=b
                        K_temp(b,a)=K_temp(b,a)+K_ele(r,s);
%                       K_Patch(b,a)=K_Patch(b,a)+K_ele(r,s);
                       end
                     end
                 end
             end
         end
      end
  end
  R_Patch=R_Patch+R_temp;K_Patch=K_Patch+K_temp;
end
end
