#include <iostream>
#include <cmath>
#include"../Eigen/Eigen/Dense"
#include"../Eigen/Eigen/Sparse"
#include <tuple>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.
typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

// Function Description
// to be oontinued ......
tuple<SpMat,MatrixXd>PatchN_KR_bs(int Nele,MatrixXd P,MatrixXd Q,MatrixXd End_pt,MatrixXi Ele2Cp,RowVectorXd knotVector,int order,RowVectorXd weights,MatrixXd Mat,MatrixXi DOF,Matrix3d A0,int ngp)
{ 
    int i,j,j1,k,k1,a,b,nos,NOS,r,s,dof,num_dof,num_cp;
    dof=4; // dofs per cpt.
    nos = order+1; // no of shape funs per element.
	num_cp=(Nele+order);// total no of cpts in a patch.
	num_dof=dof*(num_cp); // total dofs in a patch.
    
	SpMat K_patch(num_dof,num_dof); // initializing patch K(consistent tangent)(SPARSE MATRIX).
	
	MatrixXd K_ele(4*nos,4*nos),R_patch(num_dof,1),R_ele(4*nos,1);; // initializing element K(consistent tangent),R and patch R (DENSE MATRIX).
	R_patch = MatrixXd::Zero(num_dof,1);
	K_ele = MatrixXd::Zero(4*nos,4*nos);R_ele = MatrixXd::Zero(4*nos,1);
	RowVector2d vec1;RowVectorXd vec3(6);
	RowVectorXi vec2(order+1);
	
	
	for ( i=0;i<Nele;i++ )
	{   
	    vec1=End_pt.row(i);vec2=Ele2Cp.row(i);vec3=Mat;
		tie(K_ele,R_ele)=Element_KR_bs(P,Q,vec1,vec2,knotVector,order,weights,vec3,A0,ngp);
		for ( j=0;j<nos;j++ )
		{   
			for ( j1=0;j1<dof;j1++ )
		    {   
		    	r=dof*j+j1; // rth dof (local).
		    	R_patch(DOF(Ele2Cp(i,j)-1,j1)-1) += R_ele(r);
			    for (k=0;k<nos;k++)
			    {   
			    	for ( k1=0;k1<dof;k1++ )
			    	{   
			    		s=dof*k+k1;// sth dof (local).
			    		if (r<=s)
			    		{
			    			a=DOF(Ele2Cp(i,j)-1,j1);b=DOF(Ele2Cp(i,k)-1,k1);
			    		   	K_patch.coeffRef(a-1,b-1)+= K_ele(r,s);
						}
					}
				}
			
			
		    }
		 } 
	}
	return make_tuple(K_patch,R_patch);
}
