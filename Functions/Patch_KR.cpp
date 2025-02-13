#include <iostream>
#include <cmath>
#include"../Eigen/Eigen/Dense"
#include <tuple>
using namespace std;
using namespace Eigen;

// Function Description
// to be oontinued ......
tuple<MatrixXd,MatrixXd>Patch_KR(int Nele,MatrixXd &P,MatrixXd &Q,MatrixXd &End_pt,MatrixXd &Ele2Cp,RowVectorXd &knotVector,int order,RowVectorXd &weights,MatrixXd &Mat,MatrixXd &DOF,int NC,Matrix3d &A0,int ngp)
{ 
    int i,j,j1,k,k1,a,b,nos,r,s,dof;
    dof=4;
	MatrixXd K_patch(NC-1,NC-1),K_ele(4*nos,4*nos); // initializing patch K(consistent tangent) and R (resiude).
	VectorXd R_patch(NC-1),R_ele(4*nos); // initializing patch K(consistent tangent) and R (resiude).
	K_patch = MatrixXd::Zero(NC-1,NC-1);R_patch = VectorXd::Zero(NC-1);
	nos = order+1;
	K_ele = MatrixXd::Zero(4*nos,4*nos);R_ele = VectorXd::Zero(4*nos);
	RowVector2d vec1;RowVectorXd vec2(order+1),vec3(6);
	for ( i=0;i<Nele;i++ )
	{   cout << "Element: " << i+1 << endl << endl;  
	    vec1=End_pt.row(i);vec2=Ele2Cp.row(i);vec3=Mat.row(i);
		tie(K_ele,R_ele)=Element_KR(P,Q,vec1,vec2,knotVector,order,weights,vec3,A0,ngp);
		for ( j=0;j<nos;j++ )
		{   
			for ( j1=0;j1<nos;j1++ )
		    {   
		    	r=dof*j+j1; // rth dof (local).
		    	if (DOF(Ele2Cp(i,j)-1,j1)<NC)
		    	{
		    		R_patch(DOF(Ele2Cp(i,j)-1,j1)-1) += R_ele(r);
				}
			    for (k=0;k<nos;k++)
			    {   
			    	for ( k1=0;k1<nos;k1++ )
			    	{   
			    		s=dof*k+k1;// sth dof (local).
			    		if (r<=s)
			    		{
			    			a=DOF(Ele2Cp(i,j)-1,j1);b=DOF(Ele2Cp(i,k)-1,k1);
			    		   	if(max(a,b)<NC)
			    		   	{
			    		   		K_patch(a-1,b-1)+=K_ele(r,s);
			    		   		if(a!=b)
			    		   	    {
			    		   		   K_patch(b-1,a-1)+=K_ele(r,s);
			    		   		
							    }
							}
						}
					}
				}
			
			
		    }
		 } 
	}
	return make_tuple(K_patch,R_patch);
}
