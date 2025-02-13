#include <iostream>
#include <cmath>
#include"../Eigen/Eigen/Dense"
#include"../Eigen/Eigen/Sparse"
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
tuple<SpMat_RM,MatrixXd,MatrixXd,MatrixXd>Network_KR_Parallel_Debug(MatrixXd &P,MatrixXd &Q,vector<Patch_Structure> &Patch,vector<Patch_Structure> &Patch_Bs, MatrixXd &A0, MatrixXd &A0_Bs, RowVectorXd &Mat,RowVectorXd &Mat_Arcs, RowVectorXd &Mat_Bs, MatrixXi &DOF,MatrixXd &D_DOF,MatrixXd &F_DOF,double u0,int num_threads)
{   
	//Scalars
	int i,j;
	int order,num_cpt,Nele,ngp,len_kv,patch_size,patch_bs_size;
	patch_size=Patch.size();
	patch_bs_size=Patch_Bs.size();    
	// Vectors
    RowVectorXd kv,wts;
	
	MatrixXd R_Patch,R_Patch_Bs;
	MatrixXd R,Fluxes;
	R.resize(4*P.cols(),1);R=MatrixXd::Zero(4*P.cols(),1);
	Fluxes.resize(4*P.cols(),1);Fluxes=MatrixXd::Zero(4*P.cols(),1);
    
	// Matrices
    MatrixXd p,q,End_pt;
    Matrix3d a0; // matrix containing triad of a patch.
    MatrixXi n; // matrix containing IDs of cpt in a patch.
    MatrixXi dof,Ele2Cp; //  matrix containing local dof of a patch.
	
	
	// Sparse Matrices.
	SpMat_CM K_Patch_Bs,K_Patch;
	SpMat_RM K(4*P.cols(),4*P.cols());
	double sparcity=0.003; // Max no of elements in a column ( choosing this properly saves time in assembly)
    // Note: Sparcity is observed in different meshes of same network case.
    // It might be possible to calculate max no of elements in a row/column, I haven't done it yet.
    K.reserve(VectorXi::Constant(4*P.cols(),pow(sparcity,0.5)*4*P.cols())); // This step saves time while assembling K. 'coeffRef' function takes much time if space is reserves. 
	
	// Modify Q as per boundry conditions.
	cout << "Modify Q as per boundry conditions." << endl << endl;
    auto t1=steady_clock::now();
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
    auto t2=steady_clock::now();
	double del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
	cout << "Time required : " << del_t << " seconds" << endl << endl;
	
    //// Form K and R without bending strip for the whole network. ////
    // Form vector of K_Patch and R_Patch
	vector <SpMat_CM> K_vector,K_vector_bs;
	vector <MatrixXd> R_vector,R_vector_bs; 
	K_vector.resize(patch_size);R_vector.resize(patch_size);
	K_vector_bs.resize(patch_bs_size);R_vector_bs.resize(patch_bs_size);
    
    //***************  PARALLEL PART OF THE FUNCTION  ***************//
   cout << "Calculating K and R for fibers" << endl;
    t1=steady_clock::now();

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
			SpMat_CM K_Patch;
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
			K_Patch.resize(4*num_cpt,4*num_cpt);R_Patch.resize(4*num_cpt,1);
			if( i < patch_size - patch_bs_size)
			{   
				tie(K_Patch,R_Patch) = PatchN_KR(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
			}
			else
			{   
				tie(K_Patch,R_Patch) = PatchN_KR(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Arcs,dof,a0,ngp);
			}	
	        K_vector[i]=K_Patch;R_vector[i]=R_Patch;			
	   }   // END PARALLEL FOR LOOP 
	   #pragma omp single
	   {
	    t2=steady_clock::now();
	    del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
	    cout << "Time required : " << del_t << " seconds" << endl << endl;
	   }
	   
	   //**** Adding bending strips ****//
           #pragma omp single
          cout << "Calculating K and R for bending strips" << endl;
           t1=steady_clock::now();
	    
        #pragma omp for schedule(dynamic) //  START PARALLEL FOR LOOP
	    for(i=0;i<patch_bs_size; i++) // Get Patch stiffness and residue
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
			SpMat_CM K_Patch_Bs;
			MatrixXd R_Patch_Bs;
	        
			num_cpt=Patch_Bs[i].num_cpt;len_kv=Patch_Bs[i].len_kv;order=Patch_Bs[i].order;Nele=num_cpt-order;
			n.resize(1,num_cpt);
			n=Patch_Bs[i].patch2cpt; // Cpt IDs of the patch.
			kv.resize(1,len_kv);kv=Patch_Bs[i].kv; // knotvector of the patch.
			wts.resize(1,num_cpt);wts=Patch_Bs[i].wts; // weights of the patch.
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

	        K_vector_bs[i]=K_Patch_Bs;R_vector_bs[i]=R_Patch_Bs;	
	   }  // END PARALLEL FOR LOOP 	   
   }   // END PARALLEL THREADS 
   
    t2=steady_clock::now();
    del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
	cout << "Time required : " << del_t << " seconds" << endl << endl;
     
  cout << "Assembly of K and R of fibers (NEW)" << endl;
   t1=steady_clock::now();
   //**** Assemble contribution of K_Patch and R_Patch to network K and R. ****//
   int r,r1,s,s1,a,b;
   for(int i=0;i<patch_size; i++) // demonstrate that making vector of K_Patch works.
   {	
		num_cpt=Patch[i].num_cpt; order=Patch[i].order; Nele=num_cpt-order;
		n.resize(1,num_cpt);n=Patch[i].patch2cpt; // Cpt IDs of the patch.
   	    int key_K=1;dof.resize(num_cpt,4);
	    for (int j=0;j<num_cpt;j++) // local ID of dofs
		{
			dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
		}
		// Get K_Patch and R_Patch for ith patch.
		K_Patch.resize(4*num_cpt,4*num_cpt);R_Patch.resize(4*num_cpt,1);
		K_Patch = K_vector[i];R_Patch = R_vector[i];
		
		//  Add Patch contribution to K and R (Since K_Patch is sparse only iterate over non-zero elements)
		// R
		for (r=0;r<num_cpt;r++)
		   for (r1=0;r1<4;r1++)
			{   
			    a = n(0,r);
				R(DOF(a-1,r1)-1,0) += R_Patch(dof(r,r1)-1,0);
		    }

		//K
		double value;
		for (int col_index=0; col_index<K_Patch.outerSize(); ++col_index) // Searching in each column
         { 
		   s=col_index/4;s1=col_index%4;b = n(0,s);
		   for (SpMat_CM::InnerIterator it(K_Patch,col_index); it; ++it) //looking for rows with nnz
            {
			    int row_index=it.row();value=it.value();
				r=row_index/4;r1=row_index%4;a=n(0,r);
				K.coeffRef(DOF(a-1,r1)-1,DOF(b-1,s1)-1) += value;
			    if (DOF(a-1,r1)!=DOF(b-1,s1))
			    {
				  K.coeffRef(DOF(b-1,s1)-1,DOF(a-1,r1)-1) += value;
			    }
		    }
		}	
   }
   t2=steady_clock::now();
   del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
  cout << "Time required : " << del_t << " seconds" << endl << endl;
      
  cout << "Assembly of K and R of bending strips (NEW)" << endl;
   t1=steady_clock::now();
   //**** Assemble contribution of K_Patch_Bs and R_Patch_Bs to network K and R. ****//
   for(int i=0;i<patch_bs_size; i++) // demonstarte that making vector of K_Patch works.
   {
//   	    cout << "Patch no : " << i << endl << endl;
            num_cpt=Patch_Bs[i].num_cpt; 
		n.resize(1,num_cpt);n=Patch_Bs[i].patch2cpt; // Cpt IDs of the patch.
   	    int key_K=1;dof.resize(num_cpt,4);
	    for (int j=0;j<num_cpt;j++) // local ID of dofs
		{
			dof(j,0)=key_K;dof(j,1)=key_K+1;dof(j,2)=key_K+2;dof(j,3)=key_K+3;key_K=key_K+4;
		}
		// Get K_Patch and R_Patch for ith patch.
		K_Patch_Bs.resize(4*num_cpt,4*num_cpt);R_Patch_Bs.resize(4*num_cpt,1);
		K_Patch_Bs = K_vector_bs[i];R_Patch_Bs = R_vector_bs[i];
		
		//  Add Patch contribution to K and R (Since K_Patch is sparse only iterate over non-zero elements)
		// R
		for (r=0;r<num_cpt;r++)
		   for (r1=0;r1<4;r1++)
			{   
			    a = n(0,r);
				R(DOF(a-1,r1)-1,0) += R_Patch_Bs(dof(r,r1)-1,0);
		    }
               // cout << "R assembly done" << endl << endl;

		//K
		double value;
		for (int col_index=0; col_index<K_Patch_Bs.outerSize(); ++col_index) // Searching in each column
         { 
		   s=col_index/4;s1=col_index%4;b = n(0,s);
		   for (SpMat_CM::InnerIterator it(K_Patch_Bs,col_index); it; ++it) //looking for rows with nnz
            {
			    int row_index=it.row();value=it.value();
				r=row_index/4;r1=row_index%4;a=n(0,r);
				K.coeffRef(DOF(a-1,r1)-1,DOF(b-1,s1)-1) += value;
			    if (DOF(a-1,r1)!=DOF(b-1,s1))
			    {
				  K.coeffRef(DOF(b-1,s1)-1,DOF(a-1,r1)-1) += value;
			    }
		    }
		}
        // cout << "K assembly done" << endl << endl;

   }
   t2=steady_clock::now();
   del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
   cout << "Time required : " << del_t << " seconds" << endl << endl;
    
   //**** Store Fluxes Before modifying R (Residual) as per boundary conditions ****//
   Fluxes=R;

	//**** Apply Boundary conditions ****//
   // Apply Neumann BC  
   cout << "Apply Neumann BC" << endl;
    t1=steady_clock::now();
    int key_K;
	for(int i=0;i<P.cols();i++)
    {
		for(int j=0; j<4; j++)
   	    {
			key_K = DOF(i,j);
			R(key_K-1,0) -= F_DOF(i,j);  		
	    }
    }
    t2=steady_clock::now();
    del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
   cout << "Time required : " << del_t << " seconds" << endl << endl;
   
   // Apply Dirichlet's BC
   cout << "Apply Dirichlet's BC" << endl;
   t1=steady_clock::now();
    SpVec_RM Z1(P.cols()*4); // This vector used to set row of tangent stiffness zero.
	for(int i=0;i<P.cols();i++)
    {
		for(int j=0; j<4; j++)
   	    {
   	    	if(D_DOF(i,j)!=u0)
   	    	{   
   	    		key_K=DOF(i,j);
   	    		//Modify K}
                K.row(key_K-1) = Z1; // All the elements in the row set to zero.
				K.coeffRef(key_K-1,key_K-1) = 1; // This element at row and column of corresponsing dof is set to 1.
				// Modify R
   	    		R(key_K-1,0)=0;
			}			   		
	    }
    }
   t2=steady_clock::now(); 
   del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
  cout << "Time required : " << del_t << " seconds" << endl << endl;    
   K.makeCompressed();
   return make_tuple(K,R,Fluxes,Q);
}


















