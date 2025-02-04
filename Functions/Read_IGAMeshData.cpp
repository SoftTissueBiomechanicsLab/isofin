// Function to read IGA mesh data stored in text file.
#include <iostream>
#include <fstream>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
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

void Read_IGAMeshData(string &filepath,MatrixXd &Cpt,RowVectorXd &Wts,vector<Patch_Structure> &Patch,vector<Patch_Structure> &Patch_Bs,MatrixXd &A0,MatrixXd &A0_Bs)
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
//	cout << "Num_Cpt: " << Num_Cpt <<endl << endl;
	Cpt.resize(4,Num_Cpt);Wts.resize(1,Num_Cpt);
	
	// Read control point coordinates.
//	cout << "Read control point coordinates" << endl << endl;
	for (i=0; i<Num_Cpt; i++) 
	{   for(j=0;j<4;j++)
		{
			file >> value; Cpt(j,i)=value;
		}
//		cout << "i : " << i << endl;
//		cout << "Cpt(:,i): " << Cpt(0,i) << Cpt(1,i) << Cpt(2,i) << Cpt(3,i) << endl << endl;
    }
//    cout << "Cpt(0,Num_Cpt-1): " << Cpt(0,Num_Cpt-1) <<endl << endl;cout << "Cpt(1,Num_Cpt-1): " << Cpt(1,Num_Cpt-1) <<endl << endl;
//    cout << "Cpt(2,Num_Cpt-1): " << Cpt(2,Num_Cpt-1) <<endl << endl;cout << "Cpt(3,Num_Cpt-1): " << Cpt(3,Num_Cpt-1) <<endl << endl;
    // Read weights.
//    cout << "Read weights " << endl << endl;
    for (i=0; i<Num_Cpt; i++) 
	{  
			file >> value; 
//			cout << "value: " << value <<endl << endl;
			Wts(0,i)=value;
    }
    
//    cout << "value: " << value <<endl << endl;
//    cout << "Wts(0,Num_Cpt-2): " << Wts(0,Num_Cpt-2) <<endl << endl;
//    cout << "Wts(0,Num_Cpt-1): " << Wts(0,Num_Cpt-1) <<endl << endl;
    
	// Read Patch data.
//    cout << "Read Patch data" << endl << endl;
    file >> num_patch;
//	cout << "num_patch: "<< num_patch << endl;
	Patch.resize(num_patch);
//    cout << "Soham" << endl;

	for (i=0; i<num_patch; i++) 
	{       
//	        cout << "Patch No : " << i << endl << endl;
			file >> len_kv;file >> num_cpt;file >> order;
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
//    cout << "Read Patch_Bs data" << endl << endl;
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
//    cout << "Read Triads A0" << endl << endl;
    A0.resize(3,3*num_patch);
    for (i=0;i<num_patch;i++)
    {
    	for (j=0;j<3;j++)
    	{
    		file >> value;A0(j,3*i)=value;file >> value;A0(j,3*i+1)=value;file >> value;A0(j,3*i+2)=value;
		}
	}
    
	// Read Triads A0_Bs
//	cout << "Read Triads A0_Bs" << endl << endl;
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
