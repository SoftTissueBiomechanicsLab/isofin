// Function to read IGA mesh data stored in text file.
#include <iostream>
#include <fstream>
#include"../Eigen/Eigen/Dense"
#include<string>
#include <cmath>
#include <vector>
#include <tuple>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
// Structure for patch.
struct Patch_Structure
{   
    int order,len_kv,num_cpt;
	RowVectorXd kv,wts;
	MatrixXi patch2cpt; // contains integer
	MatrixXd cpt;
};

void Read_InputFile(string &filepath,MatrixXd &Cpt,RowVectorXd &Wts,vector<Patch_Structure> &Patch,vector<Patch_Structure> &Patch_Bs,MatrixXd &A0,MatrixXd &A0_Bs)
{
	// Scalars
	int Num_Cpt,i,j,k,num_patch,order,len_kv,num_cpt,num_patch_bs;
	double value;
	// Vectors
	// Matrices
//	MatrixXd Cpt,A0,A0_Bs;
	ifstream file;
	file.open(filepath);
	file >> Num_Cpt;
	Cpt.resize(4,Num_Cpt);Wts.resize(1,Num_Cpt);
	
	// Read control point coordinates.
	for (i=0; i<Num_Cpt; i++) 
	{   for(j=0;j<4;j++)
		{
			file >> value; Cpt(j,i)=value;
		}
    }
    
    // Read weights.
    for (i=0; i<Num_Cpt; i++) 
	{  
			file >> value; Wts(0,i)=value;
    }
    
    // Read Patch data.
    file >> num_patch;
	Patch.resize(num_patch);
    for (i=0; i<num_patch; i++) 
	{       file >> len_kv;file >> num_cpt;file >> order;
	        Patch[i].len_kv=len_kv;Patch[i].num_cpt=num_cpt;Patch[i].order=order;
			Patch[i].kv.resize(1,len_kv);Patch[i].wts.resize(1,num_cpt);
			Patch[i].patch2cpt.resize(1,num_cpt);Patch[i].cpt.resize(4,num_cpt);
			for (j=0; j<len_kv; j++) 
			{ 
				file >> value; Patch[i].kv(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{ 
				file >> value; Patch[i].wts(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{ 
				file >> value; Patch[i].patch2cpt(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{   for(k=0;k<4;k++)
				{
					file >> value; Patch[i].cpt(k,j)=value;
				}
		    }
    }
    
  
    // Read Patch_Bs data.
    file >> num_patch_bs;
    Patch_Bs.resize(num_patch_bs); 
    for (i=0; i<num_patch_bs; i++) 
	{       file >> len_kv;file >> num_cpt;file >> order;
	        Patch_Bs[i].len_kv=len_kv;Patch_Bs[i].num_cpt=num_cpt;Patch_Bs[i].order=order;
			Patch_Bs[i].kv.resize(1,len_kv);Patch_Bs[i].wts.resize(1,num_cpt);
			Patch_Bs[i].patch2cpt.resize(1,num_cpt);Patch_Bs[i].cpt.resize(4,num_cpt);
			for (j=0; j<len_kv; j++) 
			{ 
				file >> value; Patch_Bs[i].kv(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{ 
				file >> value; Patch_Bs[i].wts(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{ 
				file >> value; Patch_Bs[i].patch2cpt(0,j)=value;
    		}
    		for (j=0; j<num_cpt; j++) 
			{   for(k=0;k<4;k++)
				{
					file >> value; Patch_Bs[i].cpt(k,j)=value;
				}
		    }
    }
    
    // Read Triads A0
    A0.resize(3,3*num_patch);
    for (i=0;i<num_patch;i++)
    {
    	for (j=0;j<3;j++)
    	{
    		file >> value;A0(j,3*i)=value;file >> value;A0(j,3*i+1)=value;file >> value;A0(j,3*i+2)=value;
		}
	}
    
	// Read Triads A0_Bs
    A0_Bs.resize(3,3*num_patch_bs);
    for (i=0;i<num_patch_bs;i++)
    {
    	for (j=0;j<3;j++)
    	{
    		file >> value;A0_Bs(j,3*i)=value;file >> value;A0_Bs(j,3*i+1)=value;file >> value;A0_Bs(j,3*i+2)=value;
		}
	}
    file.close();
}
