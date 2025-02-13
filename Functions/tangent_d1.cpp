#include<iostream>
#include"../Eigen/Eigen/Dense"
using namespace Eigen;
using namespace std;
Vector3d tangent_d1( Vector3d &a1,Vector3d &a11, double norm_a1 )
{
	Vector3d t,t1,t2;
	double n = norm_a1;
	t = a11/n -(a1.dot(a11))*a1/pow(n,3);
	return t;
}
