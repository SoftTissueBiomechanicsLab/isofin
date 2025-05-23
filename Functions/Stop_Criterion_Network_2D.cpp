#include <iostream>
#include <cmath>
#include"../Eigen/Eigen/Dense"
#include"../Eigen/Eigen/Sparse"
#include <tuple>
#include <vector>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
int Stop_Criterion_Network_2D(MatrixXd P,MatrixXd Q_pre,MatrixXd Q, MatrixXd R,MatrixXd Fluxes,MatrixXd d,MatrixXi DOF,MatrixXd D_DOF,double u0,double tol_D,double tol_R,int choice )
{   int Stop_Key;
	if (choice==1) // Residue criterion
	{
		// double Ratio_R;
        // Ratio_R = R.norm()/R0.norm();
        // cout << "Norm_R:  " << R.norm() << endl;
        // cout << "Ratio_R:  " << Ratio_R << endl;
        // if (Ratio_R < tol_R)
        // {
            	// Stop_Key=1;
		// }
        // else
        // {
           // Stop_Key=0;	
	    // }

	}
	if (choice==2) // Displacement criterion
	{   
		// int num_disp,num_rot;
		// num_disp = 3*P.cols();num_rot = P.cols();
		// RowVectorXd d0_disp(num_disp),d0_rot(num_rot),d_disp(num_disp),d_rot(num_rot);
		// for (int j=0;j<Q0.cols();j++)
        // {   
        	// d0_disp(3*j)=Q0(0,j)-Q(0,j);d0_disp(3*j+1)=Q0(1,j)-Q(1,j);d0_disp(3*j+2)=Q0(2,j)-Q(2,j);
            // d0_rot(j)=Q0(3,j)-Q(3,j); 
            // d_disp(3*j)=d(DOF(j,0)-1,0);d_disp(3*j+1)=d(DOF(j,1)-1,0);d_disp(3*j+2)=d(DOF(j,2)-1,0);
            // d_rot(j)=d(DOF(j,3)-1,0); 
		// }
		// double Ratio_D_disp,Ratio_D_rot;
        // Ratio_D_disp = d_disp.norm()/d0_disp.norm();
        // Ratio_D_rot = d_rot.norm()/d0_rot.norm();
        // //cout << "d0_disp:  " << d0_disp << endl;
        // //cout << "d0_rot:  " << d0_rot << endl;
        // cout << "Ratio_D_disp:  " << Ratio_D_disp << endl;
        // cout << "Ratio_D_rot:  " << Ratio_D_rot << endl;
        // if (Ratio_D_disp < tol_D)
        // {
            	// Stop_Key=1;
		// }
        // else
        // {
           // Stop_Key=0;	
	    // } 

	} 
	if (choice==3) // Both displacement and residue criterion
	{
        MatrixXd flux_disp,r_disp,d_disp,Q_diff_disp;
        Q_diff_disp.resize(2,2*P.cols()); 
        double avg_flux_disp,ratio_R_disp,Ratio_D,ratio_U_disp,Ratio_R;
        double key_disp,temp_disp;
        
        flux_disp.resize(2*P.cols(),1);r_disp.resize(2*P.cols(),1);d_disp.resize(2*P.cols(),1);
                
        // Comparing residual to the fluxes for residual criterion
        // Checking how far we have come from assumed solution
        //Separate different types of dofs (displacement and rotation)
        key_disp=0;temp_disp=0;
        for (int i=0;i<P.cols();i++)
        {
                for (int j=0;j<2;j++)
            {
                if (D_DOF(i,j)!=u0)
                {
                    temp_disp+=pow(Fluxes(DOF(i,j)-1,0),2);
                }	
            } 
            flux_disp.block(2*i,0,2,1)=Fluxes.block(4*i,0,2,1);
            r_disp.block(2*i,0,2,1)=R.block(4*i,0,2,1);
            d_disp.block(2*i,0,2,1)=d.block(4*i,0,2,1);
            Q_diff_disp.block(0,i,2,1)=Q_pre.block(0,i,2,1)-Q.block(0,i,2,1);
        }
        flux_disp=flux_disp.cwiseAbs();r_disp=r_disp.cwiseAbs();
        d_disp=d_disp.cwiseAbs();avg_flux_disp=flux_disp.norm();
        // cout << "avg_flux_disp" << avg_flux_disp << endl;
        // avg_flux_disp=temp_disp/key_disp;avg_flux_rot=temp_rot/key_rot;
        // avg_flux_disp=pow(temp_disp,0.5);avg_flux_rot=pow(temp_rot,0.5);
        // cout << "avg_flux_disp" << avg_flux_disp << endl;
        // cout << "r_disp.maxCoeff()" << r_disp.maxCoeff() << endl;
        ratio_R_disp=r_disp.maxCoeff()/avg_flux_disp;
        
        Q_diff_disp=Q_diff_disp.cwiseAbs();d_disp=d_disp.cwiseAbs();
        ratio_U_disp=d_disp.maxCoeff()/Q_diff_disp.maxCoeff();


        // only considering displacement dof no rotation
        Ratio_D=ratio_U_disp;Ratio_R=ratio_R_disp;
        
        // cout << " ratio_R_disp " << ratio_R_disp << endl;
        // cout << " ratio_R_rot " << ratio_R_rot << endl;
        // cout << " ratio_U_disp " << ratio_U_disp << endl;
        // cout << " ratio_U_rot " << ratio_U_rot << endl;

        
        if (Ratio_D < tol_D && Ratio_R < tol_R)
        {
            Stop_Key=1;
            cout << "Ratio_D " << Ratio_D << endl;
            cout << "Ratio_R " << Ratio_R << endl;
        }
        else
        {
            Stop_Key=0;
            cout << "Ratio_D " << Ratio_D << endl;
            cout << "Ratio_R " << Ratio_R << endl;
        }
	}
	return Stop_Key;
}











