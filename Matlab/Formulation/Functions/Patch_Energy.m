function [ Eps_Patch,Kappa_Patch,ME_Patch,BE_Patch,TE_Patch,SE_Patch] = Patch_Energy(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,A0,ngp)
%% Element_Energy: Calculates strain, curvature, membrane,bending and torsion energy of the patch.

%% Output
% Eps_Patch: Membrane strain (1x(no of element) matix)
% Kappa_Patch: Bending curvatures (2x(no of element) matrix)
% ME_Patch: Membrane strain energy (1x(no of element) matix)
% BE_Patch: Bending strain energy (1x(no of element) matix)
% TE_Patch: Torsioin strain energy (1x(no of element) matix)

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

% % K_Patch=zeros(NC-1,NC-1);R_Patch=zeros(NC-1,1);% Only unconstrained dof are considered.
Eps_Patch= zeros(Nele,1);Kappa_Patch= zeros(Nele,2);ME_Patch = zeros(Nele,1);BE_Patch = zeros(Nele,1);TE_Patch = zeros(Nele,1);SE_Patch = zeros(Nele,1);
for i=1:Nele % parfor
% %   [ K_ele,R_ele ]=Element_KR(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat(i,:),A0,ngp);
  [ Eps_Element,Kappa_Element,ME_Element,BE_Element,TE_Element,SE_Element ] = Element_Energy(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat,A0,ngp);
  Eps_Patch(i,1)=Eps_Element;Kappa_Patch(i,:)=Kappa_Element';
  ME_Patch(i,1)=ME_Element;BE_Patch(i,1)=BE_Element;TE_Patch(i,1)=TE_Element;SE_Patch(i,1)=SE_Element;
end
end
