#include<iostream>
#include<cmath>
#include</Users/sotiriskakaletsis/Documents/GitHub/isofin/Eigen/Eigen/Dense>
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Rotation( Vector3d &t, double psi )
{
	Matrix3d matrix,I; 
	I = Matrix3d::Identity();
	matrix = I*cos(psi)+sin(psi)*cross_vM(t,I);
	return matrix;
}
