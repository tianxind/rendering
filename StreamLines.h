
#ifndef STREAMLINES_H
#define STREAMLINES_H

#include "mesh_definitions.h"
#include "curvature.h"
#include "OBBox.h"
#include <map>
#include <queue>

struct OBBoxEx : public OBBox {
	// What streamline does it belong to
	int oBBoxIdx;
	// What index is it in the streamline
	int streamLineIdx;
};

class StreamLines {
private:
	struct Edge {
		Edge() {}
		Edge(const Vector3f & p1, const Vector3f & p2) {pos[0] = p1; pos[1] = p2; }
		Vector3f pos[2];
		double dist[2];
	};
	struct StreamLine {
		std::vector<Edge> edges;
	};
	struct Seed {
		Seed() {}
		Seed(Vector3f & v, Mesh::FaceHandle & _f_h):f_h(_f_h) {curv.pos = v; }

		CurvatureInfo curv;
		Mesh::FaceHandle f_h;
	};
	struct StreamLinePoint { //MAKE SURE BOTH THESE Structs are changed if a change is made to one
		StreamLinePoint() {}
		StreamLinePoint(Vector3f & v, Mesh::FaceHandle & _f_h):f_h(_f_h) {curv.pos = v; }

		CurvatureInfo curv;
		Mesh::FaceHandle f_h;
	};

	// Keep a reference of the mesh data
	Mesh * m_pMesh;
	OpenMesh::VPropHandleT<CurvatureInfo> * m_pCurvature;

	// Integration step size
	float m_StepSize;

	// Essentially the halfwidth of the OBBox representing the stream line
	float m_StreamLineSep;

	int m_NumFacesVisited;

	int m_NumOBBoxes;

	/*struct OBBToAdd {
		OBBToAdd(OBBoxEx & _obb, Mesh::FaceHandle & _face):obb(_obb),face(_face) {}
		OBBoxEx obb;
		Mesh::FaceHandle face;
	};
	std::vector<OBBToAdd> m_OBBsToAdd;*/

	std::vector<StreamLine> m_StreamLines;

	std::vector<Mesh::FaceHandle> m_FacesTouched;

	Mesh::FaceHandle FindLargestFace();
	void FindNextTrinagle(std::queue<Seed> & seeds, Mesh::FaceHandle & f_h, OBBoxEx & obb);
	bool ComputeStreamLineGroup(Mesh::FaceHandle & f_h, int streamLineIdx);
	bool IntegrationStep(StreamLinePoint & slPoint, OBBoxEx & obb, float & stepSize, float dir, Mesh::FaceHandle & previousFace);
	bool IntegrateAlongDirectionField(std::queue<Seed> & seeds, float dir, int streamLineIdx);

	void AddOBBToFaces(OBBoxEx & obb, Mesh::FaceHandle & f_h, std::map<unsigned int, bool> & visitedFaces);
	bool WrapToNextFace(StreamLinePoint & slPoint, float & partialStepSize, float dir, Mesh::FaceHandle & previousFace);

	// Integration 
	void IntegrateEuler(Vector3f & pos, const Vector3f & dir, float stepSize) {
		pos = pos + stepSize*dir;
	}
public:
	StreamLines();
	//~StreamLines() {}

	void ComputeStreamLines(Mesh & _mesh, OpenMesh::VPropHandleT<CurvatureInfo> & _curvature, float streamLineSeparation);
	void DrawStreamLines();
};



#endif