// Nurbs basis and related function.
#include <Functions\FindSpan.cpp>
#include <Functions\DersBasisFun.cpp>
#include <Functions\NURBS1DBasis2ndDers.cpp>
// Gaussian quadrature points.
#include <Functions\Gauss_quad.cpp>

// tangents and derivatives
#include <Functions\tangent.cpp>
#include <Functions\tangent_d1.cpp>
#include <Functions\tangent_d11.cpp>
#include <Functions\tangent_dr.cpp>
#include <Functions\tangent_d1r.cpp>
#include <Functions\tangent_drs.cpp>
#include <Functions\tangent_d1rs.cpp>

// cross product of a vector and a matrix
#include <Functions\cross_vM.cpp>

// rotation matrix and derivatives
#include <Functions\Rotation.cpp>
#include <Functions\Rotation_d1.cpp>
#include <Functions\Rotation_dr.cpp>
#include <Functions\Rotation_d1r.cpp>
#include <Functions\Rotation_drs.cpp>
#include <Functions\Rotation_d1rs.cpp>

// lambda matrices
#include <Functions\Lambda.cpp>
#include <Functions\Lambda_d1.cpp>
#include <Functions\Lambda_dr.cpp>
#include <Functions\Lambda_d1r.cpp>
#include <Functions\Lambda_drs.cpp>
#include <Functions\Lambda_d1rs.cpp>

// bending derivatives
#include <Functions\ben_dr.cpp>
#include <Functions\ben_drs.cpp>


// torsion derivatives
#include <Functions\tor_dr.cpp>
#include <Functions\tor_drs.cpp>

// Element_KR & Element_KR_bs
#include <Functions\Element_KR.cpp>
#include <Functions\Element_KR_bs.cpp>

// PatchN_KR & PatchN_KR_bs
#include <Functions\PatchN_KR.cpp>
#include <Functions\PatchN_KR_bs.cpp>

// Read_IGAMeshData
#include <Functions\Read_IGAMeshData.cpp>

// Network_KR 
#include <Functions\Network_KR.cpp>
#include <Functions\Network_KR_Parallel.cpp>

// Modify K and R as per boundary conditions.
#include <Functions\RN_3D_BC.cpp>

// Stopping Criteion
#include <Functions\Stop_Criterion_Network_3D.cpp>



