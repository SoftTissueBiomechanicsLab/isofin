#include <iostream>
#include <cmath>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Sparse>
#include <tuple>
#include <vector>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

tuple<MatrixXd,MatrixXd> RN_3D_BC_SimpleShear(int &num_cpt,MatrixXd &Cpt,MatrixXd &Face_D,MatrixXd &Face_F,double u0,double disp_ratio)
{ // This function applies boundary conditions on a unit cube;
	vector<int> Face_x1,Face_x0,Face_y1,Face_y0,Face_z1,Face_z0,key;
	key.assign(6,0); // set key size = 6 and value of each element = 0;
	//// Identify cpts on the faces ////
	for (int i=0;i<num_cpt;i++)
	{   
		// face 1 (x=1)
		if(Cpt(0,i)==1 && Cpt(1,i)!=0 && Cpt(1,i)!=1 && Cpt(2,i)!=0 && Cpt(2,i)!=1 )
		{   
			Face_x1.resize(Face_x1.size()+1);Face_x1[key[0]]=i;key[0]=key[0]+1;
		}
		// face 2 (x=0)
		if(Cpt(0,i)==0 && Cpt(1,i)!=0 && Cpt(1,i)!=1 && Cpt(2,i)!=0 && Cpt(2,i)!=1 )
		{
			Face_x0.resize(Face_x0.size()+1);Face_x0[key[1]]=i;key[1]+=1;
		}
		// face 3 (y=1)
		if(Cpt(1,i)==1 && Cpt(0,i)!=0 && Cpt(0,i)!=1 && Cpt(2,i)!=0 && Cpt(2,i)!=1 )
		{
			Face_y1.resize(Face_y1.size()+1);Face_y1[key[2]]=i;key[2]+=1;
		}
		// face 4 (y=0)
		if(Cpt(1,i)==0 && Cpt(0,i)!=0 && Cpt(0,i)!=1 && Cpt(2,i)!=0 && Cpt(2,i)!=1 )
		{
			Face_y0.resize(Face_y0.size()+1);Face_y0[key[3]]=i;key[3]+=1;
		}
		// face 5 (z=1)(Fix rotations in all directions i.e fix second last control point of the fibers)
		if((1-Cpt(2,i))<=0.05 && Cpt(1,i)!=0 && Cpt(1,i)!=1 && Cpt(0,i)!=0 && Cpt(0,i)!=1 )
		{
			Face_z1.resize(Face_z1.size()+1);Face_z1[key[4]]=i;key[4]+=1;
		}
		// face 6 (z=0) (Fix rotations in all directions i.e fix second last control point of the fibers)
		if((Cpt(2,i)-0)<=0.05  && Cpt(1,i)!=0 && Cpt(1,i)!=1 && Cpt(0,i)!=0 && Cpt(0,i)!=1 )
		{
			Face_z0.resize(Face_z0.size()+1);Face_z0[key[5]]=i;key[5]+=1;
		}	
    }        
	
	// Form D_DOF (for applying displacement boundary conditions on the edges) ////
	MatrixXd D_DOF;D_DOF.resize(num_cpt,4);
	for (int i=0;i<num_cpt;i++)
	{
		D_DOF.row(i) << u0,u0,u0,u0;
	}
	// face 1 (x=1)
	for (int i=0;i<Face_x1.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(0,j)!=u0)
			{
				D_DOF(Face_x1[i],j) = disp_ratio*Face_D(0,j);
			}
		}
	}
	// face 2 (x=0)
	for (int i=0;i<Face_x0.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(1,j)!=u0)
			{
				D_DOF(Face_x0[i],j) = disp_ratio*Face_D(1,j);
			}
		}
	}
	// face 3 (y=1)
	for (int i=0;i<Face_y1.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(2,j)!=u0)
			{
				D_DOF(Face_y1[i],j) = disp_ratio*Face_D(2,j);
			}
		}
	}
	// face 4 (y=0)
	for (int i=0;i<Face_y0.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(3,j)!=u0)
			{
				D_DOF(Face_y0[i],j) = disp_ratio*Face_D(3,j);
			}
		}
	}
	// face 5 (z=1)
	for (int i=0;i<Face_z1.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(4,j)!=u0)
			{
				D_DOF(Face_z1[i],j) = disp_ratio*Face_D(4,j);
			}
		}
	}	
	// face 6 (z=0)
	for (int i=0;i<Face_z0.size();i++)
	{
		for(int j=0;j<Face_D.cols();j++)
		{
			if(Face_D(5,j)!=u0)
			{
				D_DOF(Face_z0[i],j) = disp_ratio*Face_D(5,j);
			}
		}
	}
	
	// Form F_DOF (for applying force boundary conditions on the edges) ////
	MatrixXd F_DOF;F_DOF.resize(num_cpt,4);F_DOF=MatrixXd::Zero(num_cpt,4);
	// face 1 (x=1)
	for (int i=0;i<Face_x1.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_x1[i],j) = Face_F(0,j)/Face_x1.size();
		}
	}
	// face 2 (x=0)
	for (int i=0;i<Face_x0.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_x0[i],j) = Face_F(1,j)/Face_x0.size();
		}
	}
	// face 3 (y=1)
	for (int i=0;i<Face_y1.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_y1[i],j) = Face_F(2,j)/Face_y1.size();
		}
	}
	// face 4 (y=0)
	for (int i=0;i<Face_y0.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_y0[i],j) = Face_F(3,j)/Face_y0.size();
		}
	}
	// face 5 (z=1)
	for (int i=0;i<Face_z1.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_z1[i],j) = Face_F(4,j)/Face_z1.size();
		}
	}	
	// face 6 (z=0)
	for (int i=0;i<Face_z0.size();i++)
	{
		for(int j=0;j<Face_F.cols();j++)
		{
				F_DOF(Face_z0[i],j) = Face_F(5,j)/Face_z0.size();
		}
	}
	
    return make_tuple(D_DOF,F_DOF);
}





























