#include<iostream>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
using namespace Eigen;
using namespace std;
Vector3d tangent( Vector3d &v )
{
	Vector3d t;
	t = v/v.norm();
	return t;
}
