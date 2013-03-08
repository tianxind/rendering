
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "mesh_definitions.h"

using Eigen::Vector3f;

float DistanceToTrianglePlane(const Vector3f & p1, const Vector3f & p2, const Vector3f & p3, const Vector3f & point);

bool PointInTriangle(const Vector3f & p1, const Vector3f & p2, const Vector3f & p3, Vector3f & point);

// Assume that ray is starting from within the triangle
bool TriangleEdgesSegmentIntersect(const Vector3f & p1, const Vector3f & p2, const Vector3f p3,
	const Vector3f & e0, const Vector3f & e1, Vector3f & intersectionPos, int & edgeIdx);


bool triAABBOverlap(Vector3f & boxcenter, Vector3f & boxhalfsize, Vector3f * triverts);



#endif