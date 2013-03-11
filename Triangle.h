
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "mesh_definitions.h"
#include "curvature.h"

using Eigen::Vector3f;

struct Triangle {
	CurvatureInfo v[3];
	Vector3f n;
};

float DistanceToTrianglePlane(const Triangle & tri, const Vector3f & point);
bool PointInTriangle(Triangle & tri, Vector3f & point);
bool TriangleEdgesSegmentIntersect(const Triangle & tri, const Vector3f & e0, const Vector3f & e1, 
	Vector3f & intersectionPos, int & edgeIdx, Mesh & mesh, Mesh::FaceHandle & currentFace, Mesh::FaceHandle & previousFace);

inline void ProjectToTrianglePlane(const Triangle & tri, Vector3f & v)
{
	float dist = tri.n.dot(v);
	v -= dist*tri.n;
}
bool triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

#endif