// Function to read IGA mesh data stored in text file.
#include <iostream>
#include <fstream>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Dense>
#include<string>
#include <cmath>
#include <vector>
#include <tuple>
#include <cmath>
#include <math.h>
using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.

void ReadInputFile(string &InputFileDir,string &MeshFileDir,string &MeshFileName,string &MeshFilePath,RowVectorXd &Mat,RowVectorXd &Mat_Bs,RowVectorXd &Mat_Arcs,double &Max_Disp,double &u0,MatrixXd &Face_D,MatrixXd &Face_F,RowVectorXi &Vec1,RowVectorXd &Vec2)
{
	
	// Read input and mesh file names.
	ifstream InputFile;
    string InputFileName,InputFilePath;
    cout<<endl<<"Please enter file name of the Input file"<<endl<<endl; // Ask for name of input file.
    cin>>InputFileName;
	InputFilePath = InputFileDir+InputFileName+".txt";
	InputFile.open(InputFilePath); // Open input file. 
	InputFile>>MeshFileName; // Read mesh file name.
	MeshFilePath=MeshFileDir+MeshFileName+".txt";
	
	
	// Read Material Properties;
	double r,Area,E,E_Bs,G,G_Arcs,E_Arcs,v,I2,I3,Ip;
	InputFile>>r; // Radius of the fibers.
	Area=M_PI*pow(r,2);// c/s area.
	InputFile>>E; // Young's Modulus of fibers.
	E_Bs = E*pow(10,2); // Young's Modulus of bending strip.
	E_Arcs = E*pow(10,2);
    G=0.5*(1+v)*E;// Shear Modulus of fibers.
	G_Arcs = G*pow(10,2); // Shear Modulus of arcs at junctions.
	InputFile>>v; // Poisson's Ratio .
	I2=(M_PI/4)*pow(r,4);I3=I2; // Area moment of inertia.
	Ip=I3+I2; // Polar moment of inertia.
	Mat(0)=Area;Mat(1)=E;Mat(2)=G;Mat(3)=I2;Mat(4)=I3;Mat(5)=Ip;
	Mat_Bs(0)=Area;Mat_Bs(1)=E_Bs;Mat_Bs(2)=G;Mat_Bs(3)=I2;Mat_Bs(4)=I3;Mat_Bs(5)=Ip;
	Mat_Arcs(0)=Area;Mat_Arcs(1)=E_Arcs;Mat_Arcs(2)=G_Arcs;Mat_Arcs(3)=I2;Mat_Arcs(4)=I3;Mat_Arcs(5)=Ip;

	//Read Boundary Conditions.
	InputFile>>Max_Disp; // max displacement of the face.
	InputFile>>u0; // Big number to identify unassigned displacement dofs.
	int i,j;
	double value;
	for (i=0; i<6; i++) 
	{   for(j=0;j<4;j++)
		{
			InputFile >> value; Face_D(i,j)=value;//Displacement BC matrix for faces.
		}
    }
    for (i=0; i<6; i++) 
	{   for(j=0;j<4;j++)
		{
			InputFile >> value; Face_F(i,j)=value;// Force BC matrix for faces.
		}
    }
////    cout << Max_Disp << endl;cout << u0 << endl;cout << Face_D << endl;cout << Face_F << endl;
	// Read Analysis parameters
	InputFile>>Vec1(0); // no of openmp threads.
	InputFile>>Vec1(1); // total increments requested.
	InputFile>>Vec1(2); // max on of iterations per increment.
	InputFile>>Vec1(3); // max attempts to reduce applied disp (Disp halved at evert attempt)
	InputFile>>Vec2(0); // displacement stopping critetion tolerence
	InputFile>>Vec2(1); // force stopping critetion tolerence
////	cout << Vec1(0) << endl;cout << Vec1(1) << endl;cout << Vec1(2) << endl;
////	cout << Vec1(3) << endl;cout << Vec2(0) << endl;cout << Vec2(1) << endl;
	InputFile.close();
	
}
