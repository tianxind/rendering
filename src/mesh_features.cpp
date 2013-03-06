#include "mesh_features.h"
#include <Eigen/Core>
using namespace OpenMesh;
using namespace Eigen;

double computeDirection(Mesh &mesh, const Mesh::FaceHandle &fh,
                        Vector3d v) {
  Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh);
  Vec3f normal = mesh.normal(fh);
  Vector3d n(normal[0], normal[1], normal[2]);
  return v.dot(n);
}

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
  // CHECK IF e IS A SILHOUETTE HERE
  Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(e, 0);
  Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(e, 1);
  Mesh::FaceHandle fh_0 = mesh.opposite_face_handle(h0);
  Mesh::FaceHandle fh_1 = mesh.opposite_face_handle(h1);
  Vec3f v0 = mesh.point(mesh.to_vertex_handle(h0));
  Vec3f v1 = mesh.point(mesh.to_vertex_handle(h1));
  Vec3f mid = (v0 + v1) * 0.5;
  Vector3d c(cameraPos[0], cameraPos[1], cameraPos[2]);
  Vector3d p(mid[0], mid[1], mid[2]);
  Vector3d v = c - p;
  double d0 = computeDirection(mesh, fh_0, v); 
  double d1 = computeDirection(mesh, fh_1, v);
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
  // -------------------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
  return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

