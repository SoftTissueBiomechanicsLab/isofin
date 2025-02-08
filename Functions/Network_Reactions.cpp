#include <iostream>
#include <cmath>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Dense>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Sparse>
#include <tuple>
#include <vector>
#include <omp.h>
#include <time.h>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace chrono;

typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.
typedef SparseMatrix<double,ColMajor> SpMat_CM; // declares a column-major sparse matrix type of double
typedef SparseMatrix<double, RowMajor> SpMat_RM; // declares a row-major sparse matrix type of double
typedef SparseVector<double,RowMajor> SpVec_RM; // declares a row-major sparse vector type of double

// Function Description
// to be oontinued ......
tuple<MatrixXd>Network_Reactions(MatrixXd &P,MatrixXd &Q,vector<Patch_Structure> &Patch,vector<Patch_Structure> &Patch_Bs, MatrixXd &A0, MatrixXd &A0_Bs, RowVectorXd &Mat,RowVectorXd &Mat_Arcs, RowVectorXd &Mat_Bs, MatrixXi &DOF,MatrixXd &D_DOF,MatrixXd &F_DOF,double u0,int num_threads)
{   
	//Scalars
	int i,j;
	int order,num_cpt,Nele,ngp,len_kv,patch_size,patch_bs_size;
	patch_size=Patch.size();
	patch_bs_size=Patch_Bs.size();    
	// Vectors
    RowVectorXd kv,wts;
	
	MatrixXd R_Patch,R_Patch_Bs;
	MatrixXd R,Reactions;
	R.resize(4*P.cols(),1);R=MatrixXd::Zero(4*P.cols(),1);
	Reactions.resize(4,P.cols());Reactions=MatrixXd::Zero(4,P.cols());
    
	// Matrices
    MatrixXd p,q,End_pt;
    Matrix3d a0; // matrix containing triad of a patch.
    MatrixXi n; // matrix containing IDs of cpt in a patch.
    MatrixXi dof,Ele2Cp; //  matrix containing local dof of a patch.
 	
    //// Form K and R without bending strip for the whole network. ////
    // Form vector of K_Patch and R_Patch
	vector <MatrixXd> R_vector; 
	R_vector.resize(patch_size);
	    
    //***************  PARALLEL PART OF THE FUNCTION  ***************//
	omp_set_num_threads(num_threads);  //   START PARALLEL THREADS 
    #pragma omp parallel 
	{   
		int i;
		//**** Form K and R without bending strip for the whole network. ****//
		#pragma omp for schedule(dynamic) //  START PARALLEL FOR LOOP 
	    for(i=0;i<patch_size; i++) // Get Patch stiffness and residue
	    {  
			///// Variables defined inside the threads.
			// Scalars
			int order,Nele,num_cpt,ngp,len_kv;
			
			// Vectors
			MatrixXi n; // matrix containing IDs of cpt in a patch.
			RowVectorXd kv,wts;
			
			// Matrices
			MatrixXd p,q,End_pt;
			MatrixXi dof,Ele2Cp; //  matrix containing local dof of a patch.
			Matrix3d a0; // matrix containing triad of a patch.
			
			MatrixXd R_Patch;
			
			num_cpt=Patch[i].num_cpt;len_kv=Patch[i].len_kv;order=Patch[i].order;Nele=num_cpt-order;
			n.resize(1,num_cpt);
			n=Patch[i].patch2cpt; // Cpt IDs of the patch.
			kv.resize(1,len_kv);kv=Patch[i].kv; // knotvector of the patch.
			wts.resize(1,num_cpt);wts=Patch[i].wts; // weights of the patch.
			ngp=order+1; // no of gauss points
			// Control point coordinates of the patch.
			p.resize(4,num_cpt);p=Patch[i].cpt;
			q.resize(4,num_cpt);q=p;
	        for (int j1=0; j1<num_cpt; j1++) // select q from Q 
			{   
				q.col(j1)=Q.col(n(0,j1)-1);
			}					
			End_pt.resize(Nele,2); // End points of the elements.
			for (int j=1;j<=Nele;j++)
			{
				End_pt(j-1,0)=1.0*(j-1)/Nele;End_pt(j-1,1)=1.0*j/Nele;
			}
			Ele2Cp.resize(Nele,order+1); // ID of control points for elements in the patch.
			for(int j=1;j<=Nele;j++)
			{
				for(int k=1;k<=order+1;k++)
				{
					Ele2Cp(j-1,k-1)= j+k-1;
				}	
			}
			int key_K=1;dof.resize(num_cpt,4);
		    for (int j=0;j<num_cpt;j++) // local ID of dofs
			{
				dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
			}
			// get triad of the patch.
			a0.row(0)=A0.col(3*i);a0.row(1)=A0.col(3*i+1);a0.row(2)=A0.col(3*i+2); // A0 stores transpose of the triad. (Mesh genearator convention!)
			// Get K and R for patch.
			// Sparse Matrices.
			R_Patch.resize(4*num_cpt,1);
			if( i < patch_size - patch_bs_size)
			{   
				tie(R_Patch) = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
			}
			else
			{   
				tie(R_Patch) = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Arcs,dof,a0,ngp);
			}	
	        R_vector[i]=R_Patch;			
	   }   // END PARALLEL FOR LOOP  	   
   }   // END PARALLEL THREADS 
   //**** Assemble contribution of R_Patch to network R. ****//
   int r,r1,s,s1,a,b;
   for(int i=0;i<patch_size; i++) // demonstarte that making vector of K_Patch works.
   {	
		num_cpt=Patch[i].num_cpt; order=Patch[i].order; Nele=num_cpt-order;
		n.resize(1,num_cpt);n=Patch[i].patch2cpt; // Cpt IDs of the patch.
   	    int key_K=1;dof.resize(num_cpt,4);
	    for (int j=0;j<num_cpt;j++) // local ID of dofs
		{
			dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
		}
		// Get K_Patch and R_Patch for ith patch.
		R_Patch.resize(4*num_cpt,1);
		R_Patch = R_vector[i];
		//  Add Patch contribution to R (Since K_Patch is sparse only iterate over non-zero elements)
		// R
		for (r=0;r<num_cpt;r++)
		   for (r1=0;r1<4;r1++)
			{   
			    a = n(0,r);
				R(DOF(a-1,r1)-1,0) += R_Patch(dof(r,r1)-1,0);
		    }	
   }
   //**** Calculate Reactions (fiber ends)
   for (i=0;i<P.cols();i++)
   {
   		for (j=0;j<4;j++)
   		{
   			if (D_DOF(i,j)==0)
			{
			   	Reactions(j,i)=R(DOF(i,j)-1,0);
			}	
		}
   }
   return make_tuple(Reactions);
}


















