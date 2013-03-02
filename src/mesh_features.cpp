#include "mesh_features.h"
#include <Eigen/Core>
using namespace OpenMesh;
using namespace Eigen;

double computeDirection(Mesh &mesh, const Mesh::FaceHandle &fh, Vec3f
                        cameraPos) {
  Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh);
  Vector3d c(cameraPos[0], cameraPos[1], cameraPos[2]);
  Vec3f p_center(0, 0, 0);
  for (; fv_it; ++fv_it) {
    p_center += mesh.point(fv_it);
  }
  p_center /= 3.0;
  Vector3d p(p_center[0], p_center[1], p_center[2]);
  Vector3d v = c - p;
  Vec3f normal = mesh.normal(fh);
  Vector3d n(normal[0], normal[1], normal[2]);
  return c.dot(n);
}

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
  // CHECK IF e IS A SILHOUETTE HERE
  Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
  Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
  Mesh::FaceHandle fh_0 = mesh.opposite_face_handle(h0);
  Mesh::FaceHandle fh_1 = mesh.opposite_face_handle(h1);
  double d0 = computeDirection(mesh, fh_0, cameraPos); 
  double d1 = computeDirection(mesh, fh_1, cameraPos);
  return d0 * d1 < 0;
  // ------------------------------------------------
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
  // CHECK IF e IS SHARP HERE
  Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
  Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
  Mesh::FaceHandle fh_0 = mesh.opposite_face_handle(h0);
  Mesh::FaceHandle fh_1 = mesh.opposite_face_handle(h1);
  Vec3f n0 = mesh.normal(fh_0);
  Vector3d normal0(n0[0], n0[1], n0[2]);   
  Vec3f n1 = mesh.normal(fh_1);
  Vector3d normal1(n1[0], n1[1], n1[2]); 
  return normal0.dot(normal1) < 1/2;
  // ------------------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
	return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

