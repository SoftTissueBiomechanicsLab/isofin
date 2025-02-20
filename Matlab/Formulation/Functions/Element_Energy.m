function [ Eps_Element,Kappa_Element,ME_Element,BE_Element,TE_Element,SE_Element ] = Element_Energy(P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,A0,ngp)
%% Element_Energy: Calculates strain, curvature, membrane,bending and torsion energy of the element.

%% Output
% Eps_Element: Membrane strain  (1x1 matix) (mean of all int point values)
% Kappa_Element: Bending curvatures (2x1 matix) (mean of all int point values)
% ME_Element: Membrane strain energy (1x1 matix) 
% BE_Element: Bending strain energy (1x1 matix)
% TE_Element: Torsioin strain energy (1x1 matix)

%% Inputs
% P0: cpt matrix containing x-y-z of all cpts of undeformed configuration. Each column of P has one cpt 
% P: cpt matrix containing x-y-z of all cpts at some newton iteration. Each column of P has one cpt
% End_pt: row vector containing end points of an element.#
% Parent Space : e; Parametric Space : E; Physical Space : x and X.
% Ele2Cp: a row vector containing index of control point (cpt indices are same as shape fun indices).
% dof: matrix with degrees of freedom corresponding to the control pts.
% row: x-y-z components of cpt(no of rows = 3 for 3D) and column = index of cpt (no of columns = total number of cpts).
% knotVector: Knot vector for a patch to which the element belongs.
% order: order of NURBS shape functions. 
% weights: wieghts associated with the knot vector.
% Mat: Material parameters.
% A0: Refrence triad vectors. (3x3 matrix of those vectors).
% T0: Initial
%%
T0=A0(:,1)/dot_p(A0(:,1),A0(:,1))^0.5; % tangent vector to the centreline in refrence triad.
Area=Mat(1);E=Mat(2);G=Mat(3);I2=Mat(4);I3=Mat(5);Ip=Mat(6);
%% Form element cpt matrix.
P_ele=zeros(4,size(Ele2Cp,2));Q_ele=P_ele;
for i=1:size(Ele2Cp,2)
    P_ele(:,i)=P(:,Ele2Cp(1,i));  % P_ele : cpt matrix containing x-y-z of all cpts within element.
    Q_ele(:,i)=Q(:,Ele2Cp(1,i));
end
nos=order+1;
% K_Element=zeros(size(Ele2Cp,2)*4);F_Element=zeros(size(Ele2Cp,2)*4,1);
%% Compute residue components and tangent stiffness components for each cp and assemble residue vector.
% Gauss point declaration and shape function computation.
[ Gauss_e, Gauss_W ] = Gauss_quad( ngp );
% Compute shape fun derivatives wrt parameter at s for all control points.
N=zeros(ngp,nos);dN1=N;dN2=N;
Ept2=End_pt(1,2);Ept1=End_pt(1,1);
J1=0.5*(End_pt(1,2)-End_pt(1,1));% Jacobian for parent-parameter pair. Note: For 1D no need of determinant.
for i=1:ngp
    s = 0.5*((Ept2-Ept1)*Gauss_e(i)+(Ept2+Ept1));% s is Gauss pt in parametric space.
    [N(i,:),dN1(i,:),dN2(i,:)]=NURBS1DBasis2ndDers(s,order,knotVector,weights);
end
ME_Element=0;BE_Element=0;TE_Element=0;eps=zeros(1,ngp);kappa=zeros(2,ngp);
for i=1:ngp % Gauss integration
    %% Energy calculation
    % Forming required vectors.
    % A1,a1 = tangent vector of the centre line for undeformed and deformed respectively.
    xc=zeros(4,1);Xc=zeros(4,1);a1=zeros(4,1);A1=zeros(4,1);a11=a1;A11=A1;
    Q_temp=Q_ele;P_temp=P_ele;N_temp=N;dN1_temp=dN1;dN2_temp=dN2;
    for m1=1:nos
      xc=xc+N(i,m1)*Q_temp(:,m1);Xc=Xc+N_temp(i,m1)*P_temp(:,m1);
      a1=a1+dN1_temp(i,m1)*Q_temp(:,m1);A1=A1+dN1_temp(i,m1)*P_temp(:,m1);
      a11=a11+dN2_temp(i,m1)*Q_temp(:,m1);A11=A11+dN2_temp(i,m1)*P_temp(:,m1);
    end
    norm_A1=dot_p(A1(1:3,1),A1(1:3,1))^0.5;norm_a1=dot_p(a1(1:3,1),a1(1:3,1))^0.5;
    PSI=Xc(4,1);PSI_1=A1(4,1);psi=xc(4,1);psi_1=a1(4,1);
    % Normalized tangent vectors.
    T=tangent(A1(1:3,1));t=tangent(a1(1:3,1));T01=zeros(3,1);
    T1=tangent_d1(A1(1:3,1),A11(1:3,1),norm_A1);t1=tangent_d1(a1(1:3,1),a11(1:3,1),norm_a1);             
    % Required matrices.      
    Ben=zeros(2,1); ben=zeros(2,1);Tor=zeros(2,1); tor=zeros(2,1);
    LT0T=Lambda(T0,T);LTt=Lambda(T,t);RT=Rotation(T,PSI);Rt=Rotation(t,psi);
    LT0T1=Lambda_d1(T0,T01,T,T1);LTt1=Lambda_d1(T,T1,t,t1);
    RT1=Rotation_d1(T,T1,PSI,PSI_1);Rt1=Rotation_d1(t,t1,psi,psi_1);
    for m1=2:3
      % For bending.
      Ben(m1-1,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,m1),A1(1:3,1));
      ben(m1-1,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,m1),a1(1:3,1));
    end
    % For torsion.
    Tor(1,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,3),RT*LT0T*A0(:,2));
    tor(1,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,3),Rt*LTt*RT*LT0T*A0(:,2));
    Tor(2,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,2),RT*LT0T*A0(:,3));
    tor(2,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,2),Rt*LTt*RT*LT0T*A0(:,3));
    % Strain and curvatures
    eps(i)=0.5*(dot_p(a1(1:3,1),a1(1:3,1))-dot_p(A1(1:3,1),A1(1:3,1))); % Stretching(Membrane) strain;
    kappa(:,i)=[(ben(1)-Ben(1)),(ben(2)-Ben(2))];
    % Membrane energy
    me = (0.5*E*Area*(eps(i)^2/(norm_A1^3)))*J1*Gauss_W(i);
    ME_Element = me + ME_Element;
    % Bending energy
    be = (0.5*(E/(norm_A1^3))*(I3*kappa(1,i)^2+I2*kappa(2,i)^2))*J1*Gauss_W(i);
    BE_Element = be + BE_Element;
    % TOrsion energy
    te = (0.25*(G/(norm_A1^2))*Ip*((tor(1)-Tor(1))^2+(tor(2)-Tor(2))^2))*J1*Gauss_W(i);
    TE_Element = te + TE_Element;
end
SE_Element=ME_Element+BE_Element+TE_Element;
Eps_Element=mean(eps);Kappa_Element=mean(kappa,2);
end




