function [ Eps_Patch,Kappa_Patch,ME_Patch,BE_Patch,TE_Patch,SE_Patch] = Patch_Energy(Nele,P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,A0,ngp)
%% Element_Energy: Calculates strain, curvature, membrane,bending and torsion energy of the patch.

%% Output
% Eps_Patch: Membrane strain (1x(no of element) matix)
% Kappa_Patch: Bending curvatures (2x(no of element) matrix)
% ME_Patch: Membrane strain energy (1x(no of element) matix)
% BE_Patch: Bending strain energy (1x(no of element) matix)
% TE_Patch: Torsioin strain energy (1x(no of element) matix)

%% Input
% Nele: No of elements in a patch.
% Dim: Dimension of the problem.
% P0: control point matrix containing x-y-z of all control points of undeformed configuration. 
% Each column of P0 has one control point 
% P: control point matrix containing x-y-z of all control point at some newton iteration. 
% Each column of P has one control point
% End_pt: Matrix conatining end points of all the elements in the patch. 
% Each row corresponds to an element. 
% Ele2Cp: Matrix containing indices of control points. Each row corresponds to an element.
% DOF: Matrix containing indices of dofs of control points. 
% Each row corresponds to a control point.
% Note: Constrained dofs are numbered highest.
%% Code

Eps_Patch= zeros(Nele,1);
Kappa_Patch= zeros(Nele,2);
ME_Patch = zeros(Nele,1);
BE_Patch = zeros(Nele,1);
TE_Patch = zeros(Nele,1);
SE_Patch = zeros(Nele,1);
for i=1:Nele 
  [ Eps_Element,Kappa_Element,ME_Element,BE_Element,TE_Element,SE_Element ] = ...
      Element_Energy(P,Q,End_pt(i,:),Ele2Cp(i,:),knotVector,order,weights,Mat,A0,ngp);
  % Store energies for each element
  Eps_Patch(i,1)=Eps_Element;
  Kappa_Patch(i,:)=Kappa_Element';
  ME_Patch(i,1)=ME_Element;
  BE_Patch(i,1)=BE_Element;
  TE_Patch(i,1)=TE_Element;
  SE_Patch(i,1)=SE_Element;
end
end
