#ifndef CURVATURE_H
#define CURVATURE_H

#include "mesh_definitions.h"

using Eigen::Vector3f;
using Eigen::Matrix3f;

struct CurvatureInfo {
	OpenMesh::Vec3f directions[2];
    
	// The position of curvature information need to lie on a vertex
	Eigen::Vector3f pos;

	// Matrix of curvature
	Eigen::Matrix3f M;
	Eigen::Matrix3d MD;

	Eigen::Vector3f dir[2];
	double curvatures[2];

	void ComputeSLERP(CurvatureInfo & v1, float t, CurvatureInfo & result)
	{
		float theta = acosf(pos.dot(v1.pos));
		float thetaDiv = 1.f/theta;
		result.pos = ((thetaDiv*sinf((1.f -t)*theta))*pos + (thetaDiv*sinf(t*theta))*v1.pos);
		result.dir[0] = ((thetaDiv*sinf((1.f -t)*theta))*dir[0] + (thetaDiv*sinf(t*theta))*v1.dir[0]).normalized();
		result.dir[1] = ((thetaDiv*sinf((1.f -t)*theta))*dir[1] + (thetaDiv*sinf(t*theta))*v1.dir[1]).normalized();
	}
};

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature);
void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative);

inline void ComputeBarycentricCoords(Vector3f & pt, Vector3f & a, Vector3f & b, Vector3f & c, float & u, float & v, float & w)
{
	// Real Time Collision Detection
	Vector3f v0 = b-a, v1 = c-a, v2 = pt-a;
	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);
	float denom = d00 * d11 - d01 * d01;
	v = (d11*d20 - d01*d21)/denom;
	w = (d00*d21 - d01*d20)/denom;
	u = 1.f - v - w;
}

// Make this nice (put it in curvautreinfo)
/*inline void ComputeBarycentricSLERP(CurvatureInfo & a, CurvatureInfo & b,
	CurvatureInfo & c, Vector3f & pt, CurvatureInfo & result)
{
	float u,v,w;
	ComputeBarycentricCoords(pt, a.pos, b.pos, c.pos, u, v, w);

	result.pos = u*a.pos + v*b.pos + w*c.pos;
	result.dir[0] = (u*a.dir[0] + v*b.dir[0] + w*c.dir[0]).normalized();
	result.dir[1] = (u*a.dir[1] + v*b.dir[1] + w*c.dir[1]).normalized();


	/*CurvatureInfo tmp0, tmp1;
	a.ComputeSLERP(b, u+v, tmp0);
	a.ComputeSLERP(c, u+v, tmp1);
	tmp0.ComputeSLERP(tmp1, v/(u+v), result);
}*/

void InterpolateCurvature(CurvatureInfo & a, CurvatureInfo & b,
	CurvatureInfo & c, Vector3f & pt, CurvatureInfo & result);

#endif

