
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) 
{
	Eigen::Matrix3d I = Matrix3d::Identity();
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		const Vec3f & vertex_iF = mesh.point(v_it.handle());
		const Vec3f & vertexNormal_iF = mesh.normal(v_it.handle());

		Vector3d v_i(vertex_iF[0], vertex_iF[1], vertex_iF[2]);
		Vector3d N_i(vertexNormal_iF[0], vertexNormal_iF[1], vertexNormal_iF[2]);

		double sumWeight_i = 0.0;
		Eigen::Matrix3d M_i = Matrix3d::Zero();
		for (Mesh::VertexOHalfedgeIter voh_v_it = mesh.voh_iter(v_it); voh_v_it; ++voh_v_it) {
			const Vec3f & vertex_jF = mesh.point(mesh.to_vertex_handle(voh_v_it));
			Vector3d v_j(vertex_jF[0], vertex_jF[1], vertex_jF[2]);

			// Calculate tangent vector
			Vector3d T_ij = (I - N_i * N_i.transpose()) * (v_i - v_j);
			T_ij.normalize();

			// Calculate curvature
			Vector3d v_ji = v_j - v_i;
			double kappa_ij = N_i.dot(v_ji) * 2.0 / v_ji.squaredNorm();

			// Calculate weight
			double w_ij = 0.0;
			// Warning: assume mesh is closed
			w_ij += mesh.calc_sector_area(voh_v_it.handle());
			w_ij += mesh.calc_sector_area(mesh.opposite_halfedge_handle(voh_v_it.handle()));

			M_i += (w_ij * kappa_ij) * (T_ij * T_ij.transpose());
			sumWeight_i += w_ij;
		}
		// Normalize M_i
		M_i /= sumWeight_i;

		// Find the two eigenvectors
		EigenSolver<Matrix3d> solver(M_i);
		CurvatureInfo info;
		/*for (int i=0, j=0; i<3; ++i) {
			double eig = real(solver.eigenvalues()(i));
			if (abs(eig) > 1e-6) {
				info.curvatures[j] = eig;
				Vector3d v = solver.pseudoEigenvectors().block(0, i, 3, 1);
				info.directions[j] = Vec3f(v[0], v[1],v[2]);
				info.dir[j++] = Vector3f(v[0], v[1], v[2]);
			}
		}*/

		Vector3d vvv(0,1,0);
		Vector3d curvDir = -N_i.cross(vvv);
		info.directions[0] = Vec3f(curvDir[0], curvDir[1], curvDir[2]);
		info.directions[1] = Vec3f(curvDir[0], curvDir[1], curvDir[2]);
		info.dir[0] = Vector3f(curvDir[0], curvDir[1], curvDir[2]);
		info.dir[1] = Vector3f(curvDir[0], curvDir[1], curvDir[2]);

		info.pos = Vector3f(v_i[0], v_i[1], v_i[2]);
		mesh.property(curvature, v_it) = info;
	}
}

void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative) {
// Compute vector to viewer and project onto tangent plane, then use components in principal directions
  // to find curvature
  computeCurvature(mesh, curvature);
  for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it !=
         mesh.vertices_end(); ++v_it) {
    Vec3f v = camPos - mesh.point(v_it.handle());
    CurvatureInfo c_info = mesh.property(curvature, v_it.handle());
    Vector3d view(v[0], v[1], v[2]);
    Vec3f t1 = c_info.directions[0];
    Vector3d T1(t1[0], t1[1], t1[2]);
    Vec3f t2 = c_info.directions[1];
    Vector3d T2(t2[0], t2[1], t2[2]);
    double sin = T1.dot(view);
    double cos = T2.dot(view);
    double kw = sin * sin * c_info.curvatures[0] + cos * cos * c_info.curvatures[1];
    mesh.property(viewCurvature, v_it) = kw;
  }
  // -------------------------------------------------------------------------------------------------

  // We'll use the finite elements piecewise hat method to find per-face gradients of the view curvature
  // CS 348a doesn't cover how to differentiate functions on a mesh (Take CS 468! Spring 2013!) so we provide code here

  for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
    double c[3];
    Vec3f p[3];
    
    Mesh::ConstFaceVertexIter fvIt = mesh.cfv_iter(it);
    for (int i = 0; i < 3; i++) {
      p[i] = mesh.point(fvIt.handle());
      c[i] = mesh.property(viewCurvature,fvIt.handle());
      ++fvIt;
    }
    
    Vec3f N = mesh.normal(it.handle());
    double area = mesh.calc_sector_area(mesh.halfedge_handle(it.handle()));
    
    mesh.property(viewCurvatureDerivative,it) = (N%(p[0]-p[2]))*(c[1]-c[0])/(2*area) + (N%(p[1]-p[0]))*(c[2]-c[0])/(2*area);
  }
}
