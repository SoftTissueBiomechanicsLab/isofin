#include </home/smm5969/Desktop/cpp/Functions/AllFunctions_JS.cpp>
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\John Snow\Functions\AllFunctions.cpp>
#include <iostream>
#include <cmath>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
//#include<Eigen/Dense>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Sparse>
//#include<Eigen/Sparse>
#include <tuple>
#include <math.h>
#include <fstream>
#include <time.h>
#include <omp.h>
#include <chrono>

// EIGENVALUE CALCULATION.
//#include <Eigen/Eigenvalues> 

// SOLVER
// Direct 
//#include <Eigen/SparseCholesky> 
//#include <Eigen/SparseLU> 
#include </home/smm5969/Desktop/cpp/Eigen/Eigen/SparseQR> 
//#include <Eigen/OrderingMethods>

// Iterative
#include </home/smm5969/Desktop/cpp/Eigen/Eigen/IterativeLinearSolvers>



using namespace std;
using namespace Eigen; 
using namespace chrono;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.
typedef SparseMatrix<double, RowMajor> SpMat_RM; // declares a row-major sparse matrix type of double
typedef SparseMatrix<double,ColMajor> SpMat_CM; // declares a column-major sparse matrix type of double
typedef SparseVector<double,RowMajor> SpVec_RM; // declares a row-major sparse vector type of double

int main ()
{    
	// String
    string InputFileDir,MeshFileDir,MeshFilePath,MeshFileName,ResultFileDir,ResultFileName;
    InputFileDir="//home//smm5969//Desktop//cpp//Debug//Input//";// directory where mesh data and input file is stored.
    MeshFileDir="//home//smm5969//Desktop//cpp//Debug//Input//MeshData//";// directory where mesh data and input file is stored.
	ResultFileDir="//home//smm5969//Desktop//cpp//Debug//Results//"; // directory where result data and input file is stored.
	// Scalars
	double Max_Disp;
	double u0; // Big number to identify unassigned displacement dofs.
    double r,Area,E,E_Bs,G,G_Arcs,v,I2,I3,Ip;
	// Vectors
	RowVectorXd Wts; // Weights of the network.
	RowVectorXd Mat(6),Mat_Arcs(6),Mat_Bs(6);// Row vectors containing material properties of the fibers.
	RowVectorXi Vec1(4); // Analysis parameters;
	RowVectorXd Vec2(2); // Analysis parameters;
	
	// Matrices
	MatrixXd Cpt,A0,A0_Bs,P,Q,Q_serial,Q_pre; // Cpts and Triads of the network.
	Matrix3d a0; // matrix containing triad of a patch.
	MatrixXd Face_D(6,4),Face_F(6,4); // 6 faces and 4 dofs
	MatrixXd D_DOF,F_DOF;
	SpMat_RM K;
	MatrixXd R,Fluxes,Reactions;
	// vector of Structures
	vector<Patch_Structure> Patch, Patch_Bs; // Patch and Patch_Bs structures of the network.

    //// Read Analysis Inputs.
	ReadInputFile(InputFileDir,MeshFileDir,MeshFileName,MeshFilePath,Mat,Mat_Bs,Mat_Arcs,Max_Disp,u0,Face_D,Face_F,Vec1,Vec2);
	cout << "Reading Input File Data done!" << endl;
	// cout << "Mat: "<< Mat << endl;cout << "Mat_Bs: "<< Mat_Bs << endl;cout << "Mat_Arcs: "<< Mat_Arcs << endl;
	Area=Mat(0);E=Mat(1);G=Mat(2);I2=Mat(3);I3=Mat(4);Ip=Mat(5);E_Bs=Mat_Bs(1);
	G_Arcs=Mat_Arcs(2);
	
	//cout<< "Soham" << endl;
	//// Read Network geometry from the text file. (Generated by matlab) ////
	Read_IGAMeshData(MeshFilePath,Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs);
    cout << "Reading IGA Mesh Data done!" << endl;
	int num_cpt= Cpt.cols();
	
  //// Dof's for fibers and arcs and bending strips ////
	MatrixXi DOF,DOF_Bs;
	DOF.resize(num_cpt,4);DOF_Bs.resize(num_cpt,4);DOF_Bs=MatrixXi::Zero(num_cpt,4);
	for (int i=0;i<num_cpt;i++)
	{
		DOF(i,0)=4*i+1;DOF(i,1)=4*i+2;DOF(i,2)=4*i+3;DOF(i,3)=4*i+4;
	}
	
	for (int i=0;i<Patch_Bs.size();i++)
	{
		for (int j=0;j<Patch_Bs[i].num_cpt;j++)
		{
			DOF_Bs.row(Patch_Bs[i].patch2cpt(j)-1)=DOF.row(Patch_Bs[i].patch2cpt(j)-1);
		}
	}
	
	//***** Boundary Conditions *****//
    //Uniaxial Tension in X
	D_DOF.resize(num_cpt,4),F_DOF.resize(num_cpt,4);
	
	P.resize(Cpt.rows(),Cpt.cols());
	P=Cpt;Q=P;Q_pre=Q;
    K.resize(4*P.cols(),4*P.cols());R.resize(4*P.cols(),1);
	Reactions.resize(4,P.cols());R=MatrixXd::Zero(4,P.cols());
	
		
	//**************** Analysis ****************//
    cout << "***************** Starting Analysis ******************" << endl;
    P=Cpt;Q=P;Q_pre=Q;
	//cout << "P: "<< P << endl << endl;
    // Write unloaded state to the result files.
    // // ResultFileName=MeshFileName+"_Result";
	// // WriteResultFile(ResultFileDir,ResultFileName,P,Reactions);
	//// Newton Iteration ////
    int num_threads,num_inc,max_numiter,max_attempts;
	num_threads=Vec1(0);num_inc=Vec1(1);max_numiter=Vec1(2);max_attempts=Vec1(3);
	double disp_ratio,disp_inc,max_disp_inc,tol_D,tol_R;
	tol_D=Vec2(0);tol_R=Vec2(1);
    disp_ratio=0.0;max_disp_inc=1.0/float(num_inc);disp_inc=max_disp_inc;
    

    int ngp,order,num_iter,attempts,Stop_Key,break_key,num_disp,num_rot;
	num_disp = 3*P.cols();num_rot = P.cols();
    

    MatrixXd d_DOF(D_DOF.rows(),D_DOF.cols()),R0_disp(num_disp,1),R0_rot(num_rot,1),Norm_R0(2,1),d(R.rows(),1);
    MatrixXd R_disp(num_disp,1),R_rot(num_rot,1),Norm_R(2,50);
    Norm_R=MatrixXd::Zero(2,50);
    double conv_order_disp,conv_order_rot;
    MatrixXd Q0(Q.rows(),Q.cols()); 
     
	int inc_num = 1; attempts=0;
	auto t1_total = steady_clock::now();
	
	while(disp_ratio < 1.0)
    {
	  disp_ratio+=disp_inc;
      cout<<"Strain applied in this increment: "<<Max_Disp*100*(disp_ratio)<<"%; Final Strain :  "<<Max_Disp*100<<"% "<< endl << endl;
      ngp=order+1;
      tie(d_DOF,F_DOF) = RN_3D_BC( num_cpt,Cpt,Face_D,Face_F,u0,disp_ratio ); // applied boundary conditions on a unit cube.
	  num_iter=1;Stop_Key=0;
      auto t1_inc = steady_clock::now();
      while(Stop_Key==0)
        {
          auto t1_iter = steady_clock::now();
          cout << "Iteration No: " << num_iter << endl << endl;
          auto t1=steady_clock::now();
		  tie(K,R,Fluxes,Q)=Network_KR_Parallel_Debug(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,d_DOF,F_DOF,u0,num_threads);
//		  tie(K,R,Q)=Network_KR_Parallel(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,d_DOF,F_DOF,u0,num_threads);
//		  tie(K,R,Q)=Network_KR(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,d_DOF,F_DOF,u0);

//		  tie(K_serial,R_serial,Q_serial)=Network_KR(P,Q_serial,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,d_DOF,F_DOF,u0);
//		  if((K-K_serial).norm()<pow(10,-6)) cout << endl <<"Tangent Stiffness from serial and parallel EQUAL" << endl << endl;
//		  else cout << endl <<"Tangent Stiffness from serial and parallel NOT EQUAL ; (K-K_serial).norm() = " << (K-K_serial).norm() << endl << endl; 
//		  if((R-R_serial).norm()<pow(10,-6)) cout << endl <<"Tangent Stiffness from serial and parallel EQUAL" <<  endl << endl;
//		  else cout <<"Residue from serial and parallel NOT EQUAL ; (R-R_serial).norm() = " << (R-R_serial).norm() << endl << endl;
		  
		  auto t2=steady_clock::now();
	      double del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
////	      cout << "Total time required for assembly of K and R : " << del_t << " seconds" << endl << endl;
//		  cout << "n : " << K.cols() << endl;
//		  cout << "nnz : " << K.nonZeros() << endl;
          if(num_iter==1)
        	{
        	  Q0=Q; // store initial Q0 for applying displacement convergence criterion.
			  // R0=R; // Normalized residue vector for first iteration, used for residue criterion.
              // Normalize R0
        	  // for (int j=0;j<Q.cols();j++)
	            // {   
		        	// R0(3*j,0)=R(DOF(j,0)-1,0)/(Area*E_Bs);R0(3*j+1,0)=R(DOF(j,1)-1,0)/(Area*E_Bs);
					// R0(3*j+2,0)=R(DOF(j,2)-1,0)/(Area*E_Bs); R0(3*j+3,0)=R(DOF(j,3)-1,0)/(pow(Ip,0.75)*E_Bs); 
				// }
	       }
	      
		   //Calculate change in position
		  // Direct Solvers
////		  t1=steady_clock::now();
////		  SparseQR <SparseMatrix<double>,NaturalOrdering<int>> solver1;
////          solver1.compute(K);
////		  if(solver1.info()!=Success) 
////			{
////                cout << "decomposition failed" << endl << endl;
////                break;
////            }
////		  d =  -solver1.solve(R);
////		  if(solver1.info()!=Success) 
////			{
////                cout << "solving failed" << endl << endl;
////                break;
////            }
////          t2=steady_clock::now();
////	      del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
////	      cout << "Time required for solving linear system (SparseQR): " << del_t << " seconds" << endl << endl;
////          cout << " d.norm(): " << d.norm() << endl << endl ;
          
		  // Iterative Solvers
		  t1=steady_clock::now();
//		  BiCGSTAB< SparseMatrix<double>,IdentityPreconditioner>  solver2;
//		  BiCGSTAB< SparseMatrix<double>,DiagonalPreconditioner<double> >  solver2;
		  int solver_threads;
		  
		  BiCGSTAB< SparseMatrix<double>,IncompleteLUT<double> >  solver2;
          solver2.compute(K);
		  d =  -solver2.solve(R);
          t2=steady_clock::now();
	      del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
////	      cout << "Time required for solving linear system (BiCGSTAB): " << del_t << " seconds" << endl << endl;
//	      cout << "#iterations:     " << solver2.iterations() << endl;
//          cout << "estimated error: " << solver2.error()      << endl;
////	      t1=steady_clock::now();
//////		  BiCGSTAB< SparseMatrix<double>,IdentityPreconditioner>  solver2;
//////		  BiCGSTAB< SparseMatrix<double>,DiagonalPreconditioner<double> >  solver2;
////		  
////		  solver_threads=1;
//////		  setNbThreads(solver_threads);
////		  omp_set_num_threads(solver_threads);
////		  cout << "No of threads Eigen using for BiCGSTAB: " << nbThreads() <<endl;
////		  #pragma omp parallel
////		  {
////		    BiCGSTAB< SparseMatrix<double>,IncompleteLUT<double> >  solver2;
////            solver2.compute(K);
////			d =  -solver2.solve(R);
////		  }
////          t2=steady_clock::now();
////	      del_t=duration_cast<microseconds>(t2-t1).count()/pow(10,6);
////	      cout << "Time required for solving linear system (BiCGSTAB): " << del_t << " seconds" << endl << endl;
		  
		  // // // Normalize R
    	  // // for (int j=0;j<Q.cols();j++)
            // // {   
	        	// // R(DOF(j,0)-1,0)=R(DOF(j,0)-1,0)/(Area*E_Bs);R(DOF(j,1)-1,0)=R(DOF(j,1)-1,0)/(Area*E_Bs);
				// // R(DOF(j,2)-1,0)=R(DOF(j,2)-1,0)/(Area*E_Bs); R(DOF(j,3)-1,0)=R(DOF(j,3)-1,0)/(pow(Ip,0.75)*E_Bs); 
			// // }
	       // // Norm_R(0,num_iter-1)=R.norm(); // vector containing norm of residue for all iterations.
	      // // // Rate of Convergence
	      // // if (num_iter>2)
	      // // {
	      	  
			// // conv_order=log(Norm_R(0,num_iter-1)/Norm_R(0,num_iter-2))/log(Norm_R(0,num_iter-2)/Norm_R(0,num_iter-3));
			// // cout << "conv_order: " << conv_order << endl<< endl; 
		  // // }
		  //// Stopping Criterion
           if(num_iter>1)
            {
//            	int choice=2; // Displacemnt criterion only
            	// int choice=1; // Residue criterion only
            	int choice=3; // Both residual and displacement criterion 
            	Stop_Key = Stop_Criterion_Network_3D(P,Q_pre,Q,R,Fluxes,d,DOF,d_DOF,u0,tol_D,tol_R,choice);
//            	cout << "Stop_Key" << Stop_Key <<endl;
			}
		  // Find new control points
		    for (int j=0; j<Q.cols();j++)
		    {   
		    	Q(0,j)+=d(DOF(j,0)-1);
				Q(1,j)+=d(DOF(j,1)-1);
				Q(2,j)+=d(DOF(j,2)-1);
				Q(3,j)+=d(DOF(j,3)-1);
            }
     
             auto t2_iter = steady_clock::now();
             double del_t_iter = duration_cast<microseconds>(t2_iter - t1_iter).count()/pow(10,6);
////             cout << "Time required for iteration no "<< num_iter << " : " <<del_t_iter << " seconds" <<  endl << endl;	
//             if (num_iter==1) break; // Check one newton iteration.
			 num_iter+=1;
			 if (num_iter>max_numiter)
			 {  cout<<endl<<endl<<"MAX ITERATIONS OF "<<max_numiter<<" REACHED, REDUCING DISPLACEMENT INCREMENT"<<endl<<endl; 
			    break_key=1;
				break;
			 }
			 else
			 {
			 	break_key=0;
			 }
	    } // END OF ITERATION
      if(break_key==0)
      {
      	// cout << "Norm_R: " << Norm_R << endl<< endl;
		Q_pre=Q;
      	// Reactions (Yet to be calculated)
		tie(Reactions)=Network_Reactions(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,d_DOF,F_DOF,u0,num_threads);
		// Save the solution
        // // ResultFileName=MeshFileName+"_Result";
	    // // WriteResultFile(ResultFileDir,ResultFileName,Q,Reactions);
	    
	    auto t2_inc = steady_clock::now();
        double del_t_inc = duration_cast<microseconds>(t2_inc - t1_inc).count()/pow(10,6);
////        cout<<"Time required for increment no "<<inc_num<<" : "<<del_t_inc<<" seconds"<<endl<<endl;
	    inc_num+=1;
	    attempts=0;
	    if(1.0-disp_ratio>disp_inc)
	  	{
	  		if(disp_inc<max_disp_inc) disp_inc=2*disp_inc;
		}
		else
		{
			disp_inc=1.0-disp_ratio;
		}
	    
	  }
	  else
	  {
	  	Q=Q_pre;
		disp_ratio-=disp_inc; // back to previous disp ratio.
		disp_inc=disp_inc/2;
		attempts+=1;  
		cout << "attempts: "<<attempts<<endl;
		if (attempts==max_attempts)
	    {
		   cout<<endl<<endl<<"STOPPING THE ANALYSIS: TOO MANY ATTEMPTS MADE FOR THE LAST INCREMENT!!!"<<endl<<endl;
		   break;
	    }
	  } 
	  
      auto t2_inc = steady_clock::now();
      double del_t_inc = duration_cast<microseconds>(t2_inc - t1_inc).count()/pow(10,6);
//      cout << "Time required for increment no "<< inc_num << " : " <<del_t_inc << " seconds" <<  endl << endl; 	
      
      // Total time elapsed since start of the analysis
      auto t2_total = steady_clock::now();
      double del_t_total = duration_cast<microseconds>(t2_total - t1_total).count()/pow(10,6);
      if ( del_t_total <=60)
      cout << endl << "Total time elapsed "<< inc_num << " : " << del_t_total << " seconds" <<  endl << endl << endl;
      if ( del_t_total > 60 && del_t_total < 3600 )
      cout << endl << "Total time elapsed "<< inc_num << " : " << del_t_total/60.0 << " minutes" <<  endl << endl << endl;
      if ( del_t_total > 3600 && del_t_total < 3600*24 )
      cout << endl << "Total time elapsed "<< inc_num << " : " << del_t_total/3600.0 << " hours" <<  endl << endl << endl;
      if ( del_t_total > 3600*24)
      cout << endl << "Total time elapsed " << " : " << del_t_inc/(3600.0*24) << " days" <<  endl << endl << endl;
      //cout << Q << endl;
	  if (inc_num==2) break; // Check one increment.
    } //END OF ANALYSIS
    
	return 0;
}


























