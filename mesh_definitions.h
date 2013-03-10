#ifndef MESH_DEFINITIONS_H
#define MESH_DEFINITIONS_H

#include <cassert>
#define _USE_MATH_DEFINES
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#define PRINTV(v) \
	std::cout << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl << std::endl;

struct MyTraits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> Mesh;

#endif