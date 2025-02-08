#include <iostream>
#include <fstream>
#include <cmath>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Dense>
//#include<Eigen/Dense>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Sparse>
//#include<Eigen/Sparse>
#include <tuple>
#include <math.h>
#include <fstream>
#include <time.h>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace Eigen;
typedef Matrix <int, Dynamic, Dynamic> MatrixXi; // define matrix for integers.
typedef Matrix <int, 1, Dynamic> RowVectorXi; // define row vector for integers.

void WriteResultFile(string &ResultFileDir,string &ResultFileName, MatrixXd &Q, MatrixXd &Reactions)
{
	string ResultFilePath,ReactionFilePath;
	ResultFilePath = ResultFileDir+ResultFileName+".txt";
	ReactionFilePath = ResultFileDir+ResultFileName +"_Reactions.txt";
	ofstream ResultFile,ResultFileReactions;
	// Deformed control points
	ResultFile.open(ResultFilePath,std::ios_base::app);
	ResultFile << Q << endl << endl << endl;
	ResultFile.close();	
	// Reactions
	ResultFileReactions.open(ReactionFilePath,std::ios_base::app);
	ResultFileReactions << Reactions << endl << endl << endl;
	ResultFileReactions.close();
	// Log file ( Time for each iteration, residue norm at each iteration, rate of convergence)
	
}
