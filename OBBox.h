
#ifndef OBBOX_H
#define OBBOX_H

#include "mesh_definitions.h"
#include "Triangle.h"

using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Matrix4f;

struct OBBox {
public:
	// Local space
	Vector3f centerL;			//The center of the OBB
	Vector3f halfWidthsL;		//The axis length of the OBB
	Vector3f rotVecL[3];		//The axis vectors of the OBB

	// World space
	Vector3f centerW;			//The center of the OBB
	Vector3f halfWidthsW;		//The axis length of the OBB
	Vector3f rotVecW[3];		//The axis vectors of the OBB
public:
	OBBox();
	//~OBBox() {}

	void ComputeOBB(const Vector3f & primaryDir, const Vector3f & nor, const Vector3f & startingPt, float length, float width);

	bool OBBoxIntersectW(OBBox & oBB);
	// Computes approximate distance between two OBB's
	// They should be closely aligend along one axis
	float ComputeClosestDistanceAprox(OBBox & oBB);
	
	//TransformOBB is done as follows:
	//centerW = center,
	//rotVecW[i] = rot*rotVecW[i],
	//halfWidthsW -= fScl*halfWidthsW
	//The transformation is clearly from world space to world space
	void TransformOBBW(Vector3f & center, Matrix4f & rot, float fScl);
	void TransformOBBW(Matrix4f & rT);
	bool IntersectTriangle(Vector3f & v1, Vector3f & v2, Vector3f & v3);

	float VolumeW();
	float SurfaceAreaW();

	void DrawBox();

	/* Converts the axis vectors into a rotation matrix of column form */
	void RotationW(Matrix3f & rot);
	/* Converts the axis vectors into a rotation matrix of row form */
	void RotationWT(Matrix3f & rot);

	OBBox & operator=(OBBox & ref);
};

inline float OBBox::VolumeW()
{
	return 8.f*
		halfWidthsW.x()*
		halfWidthsW.y()*
		halfWidthsW.z();
}

inline float OBBox::SurfaceAreaW()
{
	return 8.f*(
		halfWidthsW.x()*halfWidthsW.y() + 
		halfWidthsW.y()*halfWidthsW.z() + 
		halfWidthsW.x()*halfWidthsW.z());
}

inline void OBBox::TransformOBBW(Vector3f & center, Matrix4f & rot4x4, float fScl)
{
	Eigen::Matrix3f rot(rot4x4.block(0, 0, 3, 3));
	centerW = rot*centerL + center;
	rotVecW[0] = rot*rotVecL[0];
	rotVecW[1] = rot*rotVecL[1];
	rotVecW[2] = rot*rotVecL[2];

	halfWidthsW = fScl*halfWidthsL;
}

inline void OBBox::TransformOBBW(Matrix4f & rT)
{
	Eigen::Matrix3f rot(rT.block(0, 0, 3, 3));
	centerW = rot*centerL + Vector3f(rT(0, 3), rT(1, 3), rT(2, 3));
	rotVecW[0] = rot*rotVecL[0];
	rotVecW[1] = rot*rotVecL[1]; 
	rotVecW[2] = rot*rotVecL[2];
}

inline void OBBox::RotationW(Matrix3f & rot)
{
	rot(0, 0) = rotVecW[0].x(), rot(1, 0) = rotVecW[0].y(), rot(2, 0) = rotVecW[0].z();
	rot(0, 1) = rotVecW[1].x(), rot(1, 1) = rotVecW[1].y(), rot(2, 1) = rotVecW[1].z();
	rot(0, 2) = rotVecW[2].x(), rot(1, 2) = rotVecW[2].y(), rot(2, 2) = rotVecW[2].z();
}

inline void OBBox::RotationWT(Matrix3f & rot)
{
	rot(0, 0) = rotVecW[0].x(), rot(0, 1) = rotVecW[0].y(), rot(0, 2) = rotVecW[0].z();
	rot(1, 0) = rotVecW[1].x(), rot(1, 1) = rotVecW[1].y(), rot(1, 2) = rotVecW[1].z();
	rot(2, 0) = rotVecW[2].x(), rot(2, 1) = rotVecW[2].y(), rot(2, 2) = rotVecW[2].z();
}


#endif
