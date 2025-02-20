function [ F_Element ] = Element_R(P,Q,End_pt,Ele2Cp,knotVector,order,weights,Mat,A0,ngp)
%% Element_Residual: Gives residual vector for an element.

%% Output
% R_Element: Matrix containing residual vector. Each column corresponds to
% cpt in the element. No of rows = dimension of the problem.
% Provides residual vector of corresponding dof of all control points in an element.
% Membrane and Bending contributions considered.

%% Inputs
% P0: cpt matrix containing x-y-z of all cpts of undeformed configuration. Each column of P has one cpt 
% P: cpt matrix containing x-y-z of all cpts at some newton iteration. Each column of P has one cpt
% End_pt: row vector containing end points of an element.
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
% K_Element=zeros(size(Ele2Cp,2)*4);
F_Element=zeros(size(Ele2Cp,2)*4,1);
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
for i=1:ngp %% MAKE PARALLEL HERE
%     s = 0.5*((Ept2-Ept1)*Gauss_e(i)+(Ept2+Ept1));% s is Gauss pt in parametric space.
%     [N(i,:),dN1(i,:),dN2(i,:)]=NURBS1DBasis2ndDers(s,order,knotVector,weights);
%     K_temp=zeros(size(Ele2Cp,2)*4);
    F_temp=zeros(size(Ele2Cp,2)*4,1);
    for j=1:size(Ele2Cp,2)
        for j1=1:4
            r=4*(j-1)+j1;% rth dof (local).
            % Forming required vectors.
            % A1,a1 = tangent vector of the centre line for undeformed and deformed respectively.
              xc=zeros(4,1);Xc=zeros(4,1);a1=zeros(4,1);A1=zeros(4,1);a11=a1;A11=A1;
              ar=a1;a1r=a1;a11r=a1;
              Q_temp=Q_ele;P_temp=P_ele;N_temp=N;dN1_temp=dN1;dN2_temp=dN2;
              for m1=1:nos
                  xc=xc+N(i,m1)*Q_temp(:,m1);Xc=Xc+N_temp(i,m1)*P_temp(:,m1);
                  a1=a1+dN1_temp(i,m1)*Q_temp(:,m1);A1=A1+dN1_temp(i,m1)*P_temp(:,m1);
                  a11=a11+dN2_temp(i,m1)*Q_temp(:,m1);A11=A11+dN2_temp(i,m1)*P_temp(:,m1);
              end
              norm_A1=dot_p(A1(1:3,1),A1(1:3,1))^0.5;norm_a1=dot_p(a1(1:3,1),a1(1:3,1))^0.5;
              ar(j1,1)=N(1,j);a1r(j1,1)=dN1(i,j);a11r(j1,1)=dN2(i,j);
              PSI=Xc(4,1);PSI_1=A1(4,1);psi=xc(4,1);psi_1=a1(4,1);
              psi_r=ar(4,1);psi_1r=a1r(4,1);
            % Normalized tangent vectors.
              T=tangent(A1(1:3,1));t=tangent(a1(1:3,1));T01=zeros(3,1);
              T1=tangent_d1(A1(1:3,1),A11(1:3,1),norm_A1);t1=tangent_d1(a1(1:3,1),a11(1:3,1),norm_a1);
              tr=tangent_dr(a1(1:3,1),a1r(1:3,1),norm_a1);
              t1r=tangent_d1r(a1(1:3,1),a11(1:3,1),a1r(1:3,1),a11r(1:3,1),norm_a1);                
            % Required matrices.      
              Ben=zeros(2,1); ben=zeros(2,1);ben_r=zeros(2,1);
              Tor=zeros(2,1); tor=zeros(2,1);tor_r=zeros(2,1);
              LT0T=Lambda(T0,T);LTt=Lambda(T,t);RT=Rotation(T,PSI);Rt=Rotation(t,psi);
              LT0T1=Lambda_d1(T0,T01,T,T1);LTt1=Lambda_d1(T,T1,t,t1);
              RT1=Rotation_d1(T,T1,PSI,PSI_1);Rt1=Rotation_d1(t,t1,psi,psi_1);
              LTtr=Lambda_dr(T,t,tr);Rtr=Rotation_dr(t,tr,psi,psi_r);
              LTt1r=Lambda_d1r(T,T1,t,t1,tr,t1r);Rt1r=Rotation_d1r(t,t1,tr,t1r,psi,psi_1,psi_r,psi_1r);
            for m1=2:3
                % For bending.
                Ben(m1-1,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,m1),A1(1:3,1));
                ben(m1-1,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,m1),a1(1:3,1));
                ben_r(m1-1,1)=ben_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,a1(1:3,1),a1r(1:3,1),A0,m1);
            end
            % For torsion.
              Tor(1,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,3),RT*LT0T*A0(:,2));
              tor(1,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,3),Rt*LTt*RT*LT0T*A0(:,2));
              tor_r(1,1)=tor_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,A0,3,2);
              Tor(2,1)=dot_p((RT1*LT0T+RT*LT0T1)*A0(:,2),RT*LT0T*A0(:,3));
              tor(2,1)=dot_p((Rt1*LTt*RT*LT0T+Rt*LTt1*RT*LT0T+Rt*LTt*RT1*LT0T+Rt*LTt*RT*LT0T1)*A0(:,2),Rt*LTt*RT*LT0T*A0(:,3));
              tor_r(2,1)=tor_dr(LT0T,LTt,LT0T1,LTt1,LTtr,LTt1r,RT,Rt,RT1,Rt1,Rtr,Rt1r,A0,2,3);

            %% Membrane Contribution
              Fm=((Area*E*0.5*(dot_p(a1(1:3,1),a1(1:3,1))-dot_p(A1(1:3,1),A1(1:3,1)))*dot_p(a1(1:3,1),a1r(1:3,1)))/(dot_p(A1(1:3,1),A1(1:3,1)))^1.5);           
            %% Bending Contribution  
              Fb=(E/(dot_p(A1(1:3,1),A1(1:3,1)))^1.5)*(I3*(ben(1)-Ben(1))*ben_r(1)+I2*(ben(2)-Ben(2))*ben_r(2));
            %% Torsion Contribution  
              Ft=(0.5*G*Ip/dot_p(A1(1:3,1),A1(1:3,1))^0.5)*((tor(1)-Tor(1))*tor_r(1)+(tor(2)-Tor(2))*tor_r(2));
            %% Total element residue vector.  
              F_temp(r,1)=(Fm+Fb+Ft)*J1*Gauss_W(i);
%               for k=1:size(Ele2Cp,2)
%                 for k1=1:4
%                     s=4*(k-1)+k1;% sth dof (local). 
%                     if r<=s                      
%                         % Forming required vectors.
%                         % A1,a1 = tangent vector of the centre line for undeformed and deformed respectively.
%                         as=zeros(4,1);a1s=zeros(4,1);a11s=a1s;
%                         as(k1,1)=N(1,k);
%                         a1s(k1,1)=dN1(i,k);
%                         a11s(k1,1)=dN2(i,k);
%                         psi_s=as(4,1);psi_1s=a1s(4,1); 
%                         % Normalized tangent vectors.
%                         ts=tangent_dr(a1(1:3,1),a1s(1:3,1),norm_a1);t1s=tangent_d1r(a1(1:3,1),a11(1:3,1),a1s(1:3,1),a11s(1:3,1),norm_a1);
%                         trs=tangent_drs(a1(1:3,1),a1r(1:3,1),a1s(1:3,1),norm_a1);
%                         t1rs=tangent_d1rs(a1(1:3,1),a11(1:3,1),a1r(1:3,1),a1s(1:3,1),a11r(1:3,1),a11s(1:3,1),norm_a1);
%                         % Required matrices.      
%                         ben_s=zeros(2,1);ben_rs=zeros(2,1);
%                         tor_s=zeros(2,1);tor_rs=zeros(2,1);
%                         LTts=Lambda_dr(T,t,ts);Rts=Rotation_dr(t,ts,psi,psi_s);
%                         LTt1s=Lambda_d1r(T,T1,t,t1,ts,t1s);Rt1s=Rotation_d1r(t,t1,ts,t1s,psi,psi_1,psi_s,psi_1s);
%                         LTtrs=Lambda_drs(T,t,tr,ts,trs);Rtrs=Rotation_drs(t,tr,ts,trs,psi,psi_r,psi_s);
%                         LTt1rs=Lambda_d1rs(T,T1,t,t1,tr,ts,t1r,t1s,trs,t1rs);Rt1rs=Rotation_d1rs(t,t1,tr,ts,t1r,t1s,trs,t1rs,psi,psi_1,psi_r,psi_s,psi_1r,psi_1s);
%                         % Bending
%                         for m1=2:3
%                             ben_s(m1-1,1)=ben_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,a1(1:3,1),a1s(1:3,1),A0,m1);
%                             ben_rs(m1-1,1)=ben_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,a1(1:3,1),a1r(1:3,1),a1s(1:3,1),A0,m1);
%                         end
%                         % Torsion              
%                         tor_s(1,1)=tor_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,A0,3,2);
%                         tor_rs(1,1)=tor_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,A0,3,2);                
% 
%                         tor_s(2,1)=tor_dr(LT0T,LTt,LT0T1,LTt1,LTts,LTt1s,RT,Rt,RT1,Rt1,Rts,Rt1s,A0,2,3);
%                         tor_rs(2,1)=tor_drs(LT0T,LTt,LT0T1,LTt1,LTtr,LTts,LTt1r,LTt1s,LTtrs,LTt1rs,RT,Rt,RT1,Rt1,Rtr,Rts,Rt1r,Rt1s,Rtrs,Rt1rs,A0,2,3);
%                         %% Membrane Contribution
%                           Km=Area*E*(dot_p(a1(1:3,1),a1s(1:3,1))*dot_p(a1(1:3,1),a1r(1:3,1))+0.5*(dot_p(a1(1:3,1),a1(1:3,1))-dot_p(A1(1:3,1),A1(1:3,1)))*dot_p(a1s(1:3,1),a1r(1:3,1)))/dot_p(A1(1:3,1),A1(1:3,1))^1.5;          
%                         %% Bending Contribution  
%                           Kb=(E/dot_p(A1(1:3,1),A1(1:3,1))^1.5)*(I3*ben_s(1)*ben_r(1)+I2*ben_s(2)*ben_r(2)+I3*(ben(1)-Ben(1))*ben_rs(1)+I2*(ben(2)-Ben(2))*ben_rs(2));                 
%                         %% Torsion Contribution  
%                           Kt=(0.5*G*Ip/dot_p(A1(1:3,1),A1(1:3,1))^0.5)*(tor_r(1)*tor_s(1)+tor_r(2)*tor_s(2)+(tor(1)-Tor(1))*tor_rs(1)+(tor(2)-Tor(2))*tor_rs(2));
%                         %% Total Stiffness  
%                         K_temp(r,s)=(Km+Kb+Kt)*J1*Gauss_W(i);
%                     end      
%                 end
%              end
        end
    end
%     K_Element=K_Element+K_temp;
    F_Element=F_Element+F_temp;
end
end




