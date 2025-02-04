#include<iostream>
#include<cmath>
#include</home/smm5969/Desktop/cpp/Eigen/Eigen/Dense>
//#include <C:\Users\smm5969\Box Sync\Fibrin CPP\Functions\cross_vM.cpp>
using namespace Eigen;
using namespace std;
Matrix3d Lambda_d1rs( Vector3d &T, Vector3d &T1, Vector3d &t, Vector3d &t1, Vector3d &tr, Vector3d &ts, Vector3d &t1r, Vector3d &t1s, Vector3d &trs, Vector3d &t1rs )
{
	double dTt,dT1t,dTt1,dTts,dTtr,dTtrs,dTt1s,dTt1r,dTt1rs,dT1ts,dT1tr,dT1trs; 
	Vector3d CTt,CT1t,CTt1,CTts,CTtr,CTtrs,CTt1s,CTt1r,CTt1rs,CT1ts,CT1tr,CT1trs;
	Vector3d vec_sum1;
	Matrix3d matrix,I,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12; 
	
	dTt=T.dot(t);dT1t=T1.dot(t);dTt1=T.dot(t1);dTts=T.dot(ts);dTtr=T.dot(tr);dTtrs=T.dot(trs);
	dTt1s=T.dot(t1s);dTt1r=T.dot(t1r);dTt1rs=T.dot(t1rs);dT1ts=T1.dot(ts);dT1tr=T1.dot(tr);dT1trs=T1.dot(trs);
	
	
	CTt=T.cross(t);CT1t=T1.cross(t);CTt1=T.cross(t1);CTts=T.cross(ts);CTtr=T.cross(tr);CTtrs=T.cross(trs);
	CTt1s=T.cross(t1s);CTt1r=T.cross(t1r);CTt1rs=T.cross(t1rs);CT1ts=T1.cross(ts);CT1tr=T1.cross(tr);CT1trs=T1.cross(trs);
	
	
	I = Matrix3d::Identity(); 
	
	vec_sum1 = (CT1trs+CTt1rs);
	
	m1=(dT1trs+dTt1rs)*I+cross_vM(vec_sum1,I);
	m2=-6*(((dT1t+dTt1)*dTtr*dTts)/pow(1+dTt,4))*CTt*CTt.transpose();
	m3=((2*(dT1ts+dTt1s)*dTtr+2*(dT1t+dTt1)*dTtrs)/pow(1+dTt,3))*CTt*CTt.transpose();
	m4=(2*(dT1t+dTt1)*dTtr/pow(1+dTt,3))*(CTt*CTts.transpose()+CTts*CTt.transpose());
	m5=(2*(dT1tr+dTt1r)*dTts/pow(1+dTt,3)-(dT1trs+dTt1rs)/pow(1+dTt,2))*CTt*CTt.transpose();
	m6=-((dT1tr+dTt1r)/pow(1+dTt,2))*(CTt*CTts.transpose()+CTts*CTt.transpose());
	m7=(2*(dT1t+dTt1)*dTts/pow(1+dTt,3)-(dT1ts+dTt1s)/pow(1+dTt,2))*(CTt*CTtr.transpose()+CTtr*CTt.transpose());
	m8=-((dT1t+dTt1)/pow(1+dTt,2))*(CTt*CTtrs.transpose()+CTts*CTtr.transpose()+CTtr*CTts.transpose()+CTtrs*CTt.transpose());
	m9=(2*dTtr*dTts/pow(1+dTt,3)-dTtrs/pow(1+dTt,2))*(CTt*(CT1t+CTt1).transpose()+(CT1t+CTt1)*CTt.transpose());
	m10=-(dTtr/pow(1+dTt,2))*(CTt*(CT1ts+CTt1s).transpose()+CTts*(CT1t+CTt1).transpose()+(CT1t+CTt1)*CTts.transpose()+(CT1ts+CTt1s)*CTt.transpose());
	m11=-(dTts/pow(1+dTt,2))*(CTt*(CT1tr+CTt1r).transpose()+CTtr*(CT1t+CTt1).transpose()+(CT1t+CTt1)*CTtr.transpose()+(CT1tr+CTt1r)*CTt.transpose());
	m12=(1/(1+dTt))*(CTt*(CT1trs+CTt1rs).transpose()+CTts*(CT1tr+CTt1r).transpose()+CTtr*(CT1ts+CTt1s).transpose()+CTtrs*(CT1t+CTt1).transpose()+(CT1t+CTt1)*CTtrs.transpose()+(CT1ts+CTt1s)*CTtr.transpose()+(CT1tr+CTt1r)*CTts.transpose()+(CT1trs+CTt1rs)*CTt.transpose());
	
	matrix=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12;

	return matrix;
}









