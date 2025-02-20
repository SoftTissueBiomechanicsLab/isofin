function [Stop_Key] = Stop_Criterion(P,Q,R,Norm_R0,Norm_R,d,Norm_D0,Norm_D,DOF,tol_D,tol_R,choice)
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
switch choice
%% Residue only
    case 1
% %        R_disp=zeros((3/4)*size(R,1),1);
% %        for j=1:size(Q0,2)
% %            R_disp(3*(j-1)+1,1)=R(DOF(j,1),1);R_disp(3*(j-1)+2,1)=R(DOF(j,2),1);R_disp(3*(j-1)+3,1)=R(DOF(j,3),1);  
% %        end
% %        Norm_R_disp=norm(R_disp)
% %        Ratio_R_disp=Norm_R_disp/Norm_R0_disp
       if Norm_R0(2)==0
          Ratio_R_disp=Norm_R(1)/Norm_R0(1)
          if Ratio_R_disp < tol_R
             Stop_Key=1;
          else
             Stop_Key=0;
          end
       else
          Ratio_R_disp=Norm_R(1)/Norm_R0(1)
          Ratio_R_rot=Norm_R(2)/Norm_R0(2)
          if Ratio_R_disp < tol_R && Ratio_R_rot < tol_R
             Stop_Key=1;
          else
             Stop_Key=0;
          end  
       end
       
%% Displacement Only
    case 2
% %         d1_disp=zeros((3/4)*size(R,1),1);d2_disp=zeros((3/4)*size(R,1),1);d1_rot=zeros((1/4)*size(R,1),1);d2_rot=zeros((1/4)*size(R,1),1);
% %         d_disp=zeros((3/4)*size(R,1),1);d_rot=zeros((3/4)*size(R,1),1);
% %         for j=1:size(Q0,2)
% %             d1_disp(3*(j-1)+1,1)=Q0(1,j)-P(1,j);d1_disp(3*(j-1)+2,1)=Q0(2,j)-P(2,j);d1_disp(3*(j-1)+3,1)=Q0(3,j)-P(3,j);
% %             d1_rot(j,1)=Q0(4,j)-P(4,j); 
% % %             d2_disp(3*(j-1)+1,1)=Q2(1,j)-P(1,j);d2_disp(3*(j-1)+2,1)=Q2(2,j)-P(2,j);d2_disp(3*(j-1)+3,1)=Q2(3,j)-P(3,j);
% % %             d2_rot(j,1)=Q2(4,j)-P(4,j);
% %             d_disp(3*(j-1)+1,1)=d(DOF(j,1),1);d_disp(3*(j-1)+2,1)=d(DOF(j,2),1);d_disp(3*(j-1)+3,1)=d(DOF(j,3),1);
% %             d_rot(j,1)=d(DOF(j,4),1);    
% %         end
% %        Ratio_D_disp = norm(d_disp)/norm(d1_disp)
% %        if Ratio_D_disp < tol_D
% %            Stop_Key=1;
% %        else
% %            Stop_Key=0;
% %        end
        if Norm_D0(2)==0
              Ratio_D_disp=Norm_D(1)/Norm_D0(1)
              if Ratio_D_disp < tol_D
                 Stop_Key=1;
              else
                 Stop_Key=0;
              end
        else
              Ratio_D_disp=Norm_D(1)/Norm_D0(1)
              Ratio_D_rot=Norm_D(2)/Norm_D0(2)
              if Ratio_D_disp < tol_D && Ratio_D_rot < tol_D
                 Stop_Key=1;
              else
                 Stop_Key=0;
              end  
       end
%% Both
    case 3
end
end















