function [Stop_Key] = Stop_Criterion_Network_3D(P,Q0,Q,d,R,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,DOF_Bs,tol_D,tol_R,choice)
%% Stopping Criterion: Based on residue or displacement or both
%% Input
% P :
% Q1 :
% Q2:
% R :
% tol_D :
% tol_R :
%% Output
%% Code
% Area=Mat_Bs(1);E_bs=Mat_Bs(2);Ip=Mat_Bs(6);r=10^-3;
%% 
ncp=size(DOF,1); % no of control points
switch choice
%% Residue only
    case 1
        Ratio_R=Norm_R(length(Norm_R))/Norm_R(1);
        if Ratio_R<tol_R
            Stop_Key=1;
        else
            Stop_Key=0;
        end

%% Displacement Only
    case 2
        d0_disp=zeros((3/4)*size(R,1),1);d0_rot=zeros((1/4)*size(R,1),1);
        d_disp=zeros((3/4)*size(R,1),1);d_rot=zeros((3/4)*size(R,1),1);
        for j=1:size(Q0,2)
            d0_disp(3*(j-1)+1,1)=Q0(1,j)-Q(1,j);d0_disp(3*(j-1)+2,1)=Q0(2,j)-Q(2,j);d0_disp(3*(j-1)+3,1)=Q0(3,j)-Q(3,j);
            d0_rot(j,1)=Q0(4,j)-Q(4,j); 
            d_disp(3*(j-1)+1,1)=d(DOF(j,1),1);d_disp(3*(j-1)+2,1)=d(DOF(j,2),1);d_disp(3*(j-1)+3,1)=d(DOF(j,3),1);
            d_rot(j,1)=d(DOF(j,4),1);    
        end
        Ratio_D_disp = norm(d_disp)/norm(d0_disp);
%         Ratio_D_disp = max(abs(d_disp))/max(abs(d0_disp))
        if Ratio_D_disp < tol_D
            Stop_Key=1;
        else
            Stop_Key=0;
        end

%% Both displacement and residue
    case 3
        % Comparing residual to the fluxes for residual criterion
        % Checking how far we have come from assumed solution 
            [ Flux ] = Network_R_New(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,DOF_Bs);Flux=full(Flux);
            flux_disp=zeros(3*ncp,1);flux_rot=zeros(ncp,1);r_disp=zeros(ncp*3,1);r_rot=zeros(ncp,1);
            d_disp=zeros(3*ncp,1);d_rot=zeros(ncp,1);
            for k=1:ncp
                flux_disp(3*(k-1)+1:3*(k-1)+3,1)=Flux(DOF(k,1:3),1);flux_rot(k,1)=Flux(DOF(k,4),1);
                r_disp(3*(k-1)+1:3*(k-1)+3,1)=R(DOF(k,1:3),1);r_rot(k,1)=R(DOF(k,4),1);
                d_disp(3*(k-1)+1:3*(k-1)+3,1)=d(DOF(k,1:3),1);d_rot(k,1)=d(DOF(k,4),1);
            end
            avg_flux_disp=sum(abs(flux_disp))/length(flux_disp);avg_flux_rot=sum(abs(flux_rot))/length(flux_rot);
            ratio_R_disp=max(abs(r_disp))/avg_flux_disp
            ratio_R_rot=max(abs(r_rot))/avg_flux_rot
            ratio_U_disp=max(abs(d_disp))/max(max(abs(Q0(:,1:3)-Q(:,1:3))));
            ratio_U_rot=max(abs(d_rot))/max(max(abs(Q0(:,4)-Q(:,4))));
            
            % Ratio_D=max(ratio_disp,ratio_rot); Ratio_R=max(ratio_disp,ratio_rot); %both displacement and torsion
            
            Ratio_R=ratio_R_disp % only for displacement field
            Ratio_D=ratio_U_disp % only for displacement field
        if Ratio_D < tol_D && Ratio_R < tol_R
            Stop_Key=1;
        else
            Stop_Key=0;
        end  
end
end















