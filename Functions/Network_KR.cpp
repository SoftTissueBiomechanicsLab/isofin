#include <iostream>
#include <cmath>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Sparse>
#include <tuple>
#include <vector>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.
typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef SparseMatrix<double, RowMajor> SpMat2; // declares a column-major sparse matrix type of double

// Function Description
// to be oontinued ......


tuple<SpMat,MatrixXd,MatrixXd>Network_KR(MatrixXd &P,MatrixXd &Q,vector<Patch_Structure> &Patch,vector<Patch_Structure> &Patch_Bs, MatrixXd &A0, MatrixXd &A0_Bs, RowVectorXd &Mat,RowVectorXd &Mat_Arcs, RowVectorXd &Mat_Bs, MatrixXi &DOF,MatrixXd &D_DOF,MatrixXd &F_DOF,double u0)
{ 
    //Scalars
	int i,j,order,num_cpt,Nele,ngp;
    // Vectors
    RowVectorXd kv,wts;
	MatrixXd R_Patch,R_Patch_Bs,R;
	R.resize(4*P.cols(),1);R=MatrixXd::Zero(4*P.cols(),1);
    // Matrices
    MatrixXd p,q,End_pt;
    Matrix3d a0; // matrix containing triad of a patch.
    MatrixXi n; // matrix containing IDs of cpt in a patch.
    MatrixXi dof,Ele2Cp; //  matrix containing local dof of a patch.
	
	
	// Sparse Matrices.
	SpMat K_Patch_Bs,K_Patch;
	SpMat2 K;
	K.resize(4*P.cols(),4*P.cols());
	// Modify Q as per boundry conditions.
	
	for(int i=0;i<P.cols();i++)
	{   
		for(int j=0;j<4;j++)
		{
			if(D_DOF(i,j)!=u0 )
			{   	
				// Modify Q
				Q(j,i) = P(j,i) + D_DOF(i,j);
			}
		}
	}
//	double len_kv;
	int len_kv;
	
    //***************** Form K and R without bending strip for the whole network. *****************//
    // Form vector of K_Patch and R_Patch
	vector <SpMat> K_vector;vector <MatrixXd> R_vector; 
	K_vector.resize(Patch.size());R_vector.resize(Patch.size());

    for(int i=0;i<Patch.size(); i++) // Get Patch stiffness and residue
    {  
		num_cpt=Patch[i].num_cpt;len_kv=Patch[i].len_kv;order=Patch[i].order;Nele=num_cpt-order;
		n.resize(1,num_cpt);
		n=Patch[i].patch2cpt; // Cpt IDs of the patch.
		kv.resize(1,len_kv);kv=Patch[i].kv; // knotvector of the patch.
		wts.resize(1,num_cpt);wts=Patch[i].wts; // weights of the patch.
		Nele = num_cpt-order; // no of elements in a patch.
		ngp=order+1; // no of gauss points
		// Control point coordinates of the patch.
		p.resize(4,num_cpt);p=Patch[i].cpt;
		q.resize(4,num_cpt);q=p;
        for (int j1=0; j1<num_cpt; j1++) // select q from Q 
		{   
			q.col(j1)=Q.col(n(0,j1)-1);
		}		
//		cout << "q: " << q << endl;
		
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
		K_Patch.resize(4*num_cpt,4*num_cpt);R_Patch.resize(4*num_cpt,1);
		if( i < Patch.size() - Patch_Bs.size())
		{   
//		    cout << "STIFFNESS OF FIBER" << endl << endl;
			tie(K_Patch,R_Patch) = PatchN_KR(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
		}
		else
		{   
//		    cout << "STIFFNESS OF ARC" << endl << endl;
			tie(K_Patch,R_Patch) = PatchN_KR(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Arcs,dof,a0,ngp);
		}
        K_vector[i]=K_Patch;R_vector[i]=R_Patch;			
   }
  
  
  
  
   // Assemble contribution of K_Patch and R_Patch to network K and R.
   for(i=0;i<Patch.size(); i++) // demonstarte that making vector of K_Patch works.
   {
   	    num_cpt=Patch[i].num_cpt; 
		n.resize(1,num_cpt);n=Patch[i].patch2cpt; // Cpt IDs of the patch.
   	    int key_K=1;dof.resize(num_cpt,4);
	    for (int j=0;j<num_cpt;j++) // local ID of dofs
		{
			dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
		}
		// Get K_Patch and R_Patch for ith patch.
		K_Patch.resize(4*num_cpt,4*num_cpt);R_Patch.resize(4*num_cpt,1);K_Patch = K_vector[i];R_Patch = R_vector[i];
		//  Add to K and R
		for (int r=0;r<num_cpt;r++)
		{
			for (int r1=0;r1<4;r1++)
			{   
				int a; a = n(0,r);
				R(DOF(a-1,r1)-1,0) += R_Patch(dof(r,r1)-1,0);
				for (int s=0;s<num_cpt;s++)
				{
					for (int s1=0;s1<4;s1++)
					{
						int b = n(0,s);
						if(dof(r,r1)<=dof(s,s1))
						{
							K.coeffRef(DOF(a-1,r1)-1,DOF(b-1,s1)-1) += K_Patch.coeffRef(dof(r,r1)-1,dof(s,s1)-1);
							if (DOF(a-1,r1)!=DOF(b-1,s1))
							{
								K.coeffRef(DOF(b-1,s1)-1,DOF(a-1,r1)-1) += K_Patch.coeffRef(dof(r,r1)-1,dof(s,s1)-1);
							}
						}
					}
				}
			}
		}
   }

   //***************** Adding bending strips *****************//
   // Form vector of K_Patch_Bs and R_Patch_Bs 
	K_vector.resize(Patch_Bs.size());R_vector.resize(Patch_Bs.size());

    for(int i=0;i<Patch_Bs.size(); i++) // Get Patch stiffness and residue
    {  
		num_cpt=Patch_Bs[i].num_cpt;len_kv=Patch_Bs[i].len_kv;order=Patch_Bs[i].order;Nele=num_cpt-order;
		n.resize(1,num_cpt);
		n=Patch_Bs[i].patch2cpt; // Cpt IDs of the patch.
		kv.resize(1,len_kv);kv=Patch_Bs[i].kv; // knotvector of the patch.
		wts.resize(1,num_cpt);wts=Patch_Bs[i].wts; // weights of the patch.
		Nele = num_cpt-order; // no of elements in a patch.
		ngp=order+1; // no of gauss points
		// Control point coordinates of the patch.
		p.resize(4,num_cpt);p=Patch_Bs[i].cpt;
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
		a0.row(0)=A0_Bs.col(3*i);a0.row(1)=A0_Bs.col(3*i+1);a0.row(2)=A0_Bs.col(3*i+2); // A0 stores transpose of the triad. (Mesh genearator convention!)
		// Get K and R for patch.
		// Sparse Matrices.
		K_Patch_Bs.resize(4*num_cpt,4*num_cpt);R_Patch_Bs.resize(4*num_cpt,1);
	
		tie(K_Patch_Bs,R_Patch_Bs) = PatchN_KR_bs(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Bs,dof,a0,ngp);
        K_vector[i]=K_Patch_Bs;R_vector[i]=R_Patch_Bs;			
   }
   // Assemble contribution of K_Patch_Bs and R_Patch_Bs to network K and R.
   for(i=0;i<Patch_Bs.size(); i++) // demonstarte that making vector of K_Patch works.
   {
   	    num_cpt=Patch_Bs[i].num_cpt; 
		n.resize(1,num_cpt);n=Patch_Bs[i].patch2cpt; // Cpt IDs of the patch.
   	    int key_K=1;dof.resize(num_cpt,4);
	    for (int j=0;j<num_cpt;j++) // local ID of dofs
		{
			dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
		}
		// Get K_Patch and R_Patch for ith patch.
		K_Patch_Bs.resize(4*num_cpt,4*num_cpt);R_Patch_Bs.resize(4*num_cpt,1);K_Patch_Bs = K_vector[i];R_Patch_Bs = R_vector[i];
		//  Add to K and R
		for (int r=0;r<num_cpt;r++)
		{
			for (int r1=0;r1<4;r1++)
			{   
				int a; a = n(0,r);
				R(DOF(a-1,r1)-1,0) += R_Patch_Bs(dof(r,r1)-1,0);
				for (int s=0;s<num_cpt;s++)
				{
					for (int s1=0;s1<4;s1++)
					{
						int b = n(0,s);
						if(dof(r,r1)<=dof(s,s1))
						{
							K.coeffRef(DOF(a-1,r1)-1,DOF(b-1,s1)-1) += K_Patch_Bs.coeffRef(dof(r,r1)-1,dof(s,s1)-1);
							if (DOF(a-1,r1)!=DOF(b-1,s1))
							{
								K.coeffRef(DOF(b-1,s1)-1,DOF(a-1,r1)-1) += K_Patch_Bs.coeffRef(dof(r,r1)-1,dof(s,s1)-1);
							}
						}
					}
				}
			}
		}
   }

   //// Apply Boundary conditions 
   // Apply Neumann BC
   int key_K;
   for(int i=0;i<P.cols();i++)
   {
   	    for(int j=0; j<4; j++)
   	    {
   	    	key_K = DOF(i,j); R(key_K-1,0) -= F_DOF(i,j);  		
	    }
   }
   
   // Apply Dirichlet's BC
   SparseVector<double,RowMajor> Z1(P.cols()*4);
   for(int i=0;i<P.cols();i++)
   {
   	    for(int j=0; j<4; j++)
   	    {
   	    	if(D_DOF(i,j)!=u0)
   	    	{
   	    		key_K=DOF(i,j);
   	    		// Modify K
   	    		K.row(key_K-1) = Z1; // All the elements in the row set to zero.
				K.coeffRef(key_K-1,key_K-1) = 1;
				// Modify R
   	    		R(key_K-1,0)=0;
			}			   		
	    }
   }
   return make_tuple(K,R,Q);
}


















