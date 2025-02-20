function [ R_patch ] = PatchN_R_bs_New(Nele,P,Qf,End_pt,Ele2Cp,knotVector,order,weights,Mat,DOF,A0,ngp)
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
%% Form complete stiffness matrix and residue vector for a patch. 
nos=size(P,2);
K_patch=sparse(nos*4,nos*4);R_patch=sparse(nos*4,1);% Only unconstrained dof are considered.
for i=1:Nele %% MAKE PARALLEL HERE
    % get K and R for ith element.
    [ R_ele ]=Element_R_bs_New(P,Qf,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat,A0,ngp);
    % forming indices and values for sparse matrix.
    n=size(R_ele,1);
    row_index_K=zeros(n*(n+1)/2,1);column_index_K=row_index_K;key_K=1;
    row_index_R=zeros(n,1);column_index_R=ones(n,1);key_R=1;
    for j=1:size(Ele2Cp(i,:),2)
      for j1=1:4
          row_index_R(key_R,1)=DOF(Ele2Cp(i,j),j1);key_R=key_R+1;
% % %           for k=1:size(Ele2Cp(i,:),2)
% % %              for k1=1:4
% % %                  column_index_K(key_K,1)=DOF(Ele2Cp(i,j),j1);row_index_K(key_K,1)=DOF(Ele2Cp(i,k),k1);key_K=key_K+1;
% % %              end
% % %          end
      end
    end
    % Form sparse K_patch and R_patch.
% % %     K_patch=K_patch+sparse(row_index_K,column_index_K,K_ele,4*nos,4*nos);
    R_patch=R_patch+sparse(row_index_R,column_index_R,R_ele,4*nos,1);
end
end






