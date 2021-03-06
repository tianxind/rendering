
#include "StreamLines.h"
#include "OBBox.h"
#include <queue>
#include "glut.h"

using namespace OpenMesh;

// Good thing eigen and openmesh have exactly the same size for 3d vectors
#define TRIANGLE(___f_h) (m_pMesh->property(faceDataHandle, ___f_h).tri)
#define CURVATURE(___v_h) (m_pMesh->property(*m_pCurvature, ___v_h))

struct FaceData {
	FaceData():bVisited(0) {}
	//~FaceData() {}

	bool bVisited;
	Triangle tri;
	std::vector<OBBoxEx> oBBoxes;
};

FPropHandleT<FaceData> faceDataHandle;

StreamLines::StreamLines()
{
}

Mesh::FaceHandle StreamLines::FindLargestFace()
{
	/*Mesh::FaceIter face;
	float largestArea = 0;
	for (Mesh::FaceIter f_it = m_pMesh->faces_begin(), end = m_pMesh->faces_end(); f_it != end; ++f_it) {
		float area = m_pMesh->calc_sector_area(m_pMesh->halfedge_handle(f_it));
		if (area > largestArea) {
			largestArea = area;
			face = f_it;
		}
	}
	return face.handle();*/

	Mesh::FaceHandle f_h;
	return f_h;
}

void StreamLines::ComputeStreamLines(Mesh & _mesh, VPropHandleT<CurvatureInfo> & _curvature, float streamLineSeparation)
{
	m_pMesh = &_mesh;
	m_pCurvature = &_curvature;
	m_pMesh->add_property(faceDataHandle);
	// Add the triangles
	for (Mesh::FaceIter f_it = m_pMesh->faces_begin(), end=m_pMesh->faces_end(); f_it != end; ++f_it) {
		FaceData & faceData = m_pMesh->property(faceDataHandle, f_it);
		
		Mesh::FaceVertexIter fv_it = m_pMesh->fv_begin(f_it);
		faceData.tri.v[0] = CURVATURE(fv_it);
		faceData.tri.v[1] = CURVATURE(++fv_it);
		faceData.tri.v[2] = CURVATURE(++fv_it);
		faceData.tri.n = VT_EIGEN(m_pMesh->calc_face_normal(f_it));
	}

	m_StepSize = .01f;
	m_StreamLineSep = .01;//streamLineSeparation;

	m_NumFacesVisited = 0;
	int numStreamLines = 0;
	//while (m_NumFacesVisited != m_pMesh->n_faces()) {
	int counter = 0;
	m_NumOBBoxes = 0;

	m_StreamLines.clear();

	Mesh::FaceIter f_it = m_pMesh->faces_begin();
	//for (int i=0; i<400; i++) ++f_it;

	for (Mesh::FaceIter end = m_pMesh->faces_end(); f_it != end; ++f_it) {
			//if (counter++ == 64) {
			Mesh::FaceHandle f_h = f_it.handle();
			if (m_pMesh->property(faceDataHandle, f_h).bVisited) {
				continue;
			}

			// Consider if this should not be incremented if no streamline is found
			m_StreamLines.push_back(StreamLine());
			if (!ComputeStreamLineGroup(f_h, numStreamLines++))
				m_StreamLines.pop_back();
		//}
	}
	//}

	//Reset stuff
	for (int i=0, end=m_FacesTouched.size(); i<end; i++) {
		m_pMesh->property(faceDataHandle, m_FacesTouched[i]).bVisited = 0;
		m_pMesh->property(faceDataHandle, m_FacesTouched[i]).oBBoxes.clear();
	}
	m_FacesTouched.clear();
}

bool StreamLines::ComputeStreamLineGroup(Mesh::FaceHandle & f_h, int streamLineIdx)
{	
	Vector3f centroid;
	//Vector3f nor;
	//nor = VT_EIGEN(m_pMesh->calc_face_normal(f_h));
	m_pMesh->calc_face_centroid(f_h, VT_OM(centroid));
	//centroid.y() += .01f;

			/*glBegin(GL_LINES);
			glColor3f(1,0,0);

			glVertex3f(centroid.x(), centroid.y(), centroid.z());
			Vector3f newPt = centroid + .1f*nor;
			glVertex3f(newPt.x(), newPt.y(), newPt.z());

			glEnd();*/

	OBBoxEx obb;
	CurvatureInfo newCurv;
	std::queue<Seed> seeds;
	seeds.push(Seed(centroid, f_h));
	//while (!seeds.empty()) {
	//	Seed & seed = seeds.front();
		
		bool r = 0;
		r |= IntegrateAlongDirectionField(seeds, 1.f, streamLineIdx); //consider dir when looking into orthogonal streamlines
		r |= IntegrateAlongDirectionField(seeds, -1.f, streamLineIdx);
		
		// Find the next seed in which there are 2 possible locations and 2 directions
		/*ComputeIntegrationBox(seed, newCurv, obb, 1.f);
		FindNextTrinagle(seeds, seed.f_h, obb);
		ComputeIntegrationBox(seed, newCurv, obb, -1.f);
		FindNextTrinagle(seeds, seed.f_h, obb);*/

		seeds.pop();
	//}

	return r;
}
/*
void StreamLines::FindNextTrinagle(std::queue<Seed> & seeds, Mesh::FaceHandle & f_h, OBBoxEx & obb)
{
	for (Mesh::FaceHalfedgeIter fh_it = m_pMesh->fh_begin(f_h), end = m_pMesh->fh_end(f_h); fh_it != end; ++fh_it) {
		Mesh::FaceHandle nextFaceH = m_pMesh->face_handle(fh_it);
		FaceData & face = m_pMesh->property(faceDataHandle, nextFaceH);

		
		
		if (face.oBBoxes[i].OBBoxIntersectW(obb)) {
		}
	}
}
*/

bool StreamLines::IntegrateAlongDirectionField(std::queue<Seed> & seeds, float dir, int streamLineIdx)
{
	Seed & seed = seeds.front();
	StreamLinePoint slPoint;
	slPoint.curv = seed.curv;
	slPoint.f_h = seed.f_h;

	float totalStreamLineLen = 0;
	float stepSize = m_StepSize;
	int counter = 0;
	m_NumOBBoxes = 0;

	Mesh::FaceHandle previousFace;
	previousFace.invalidate();
	while (counter++ < 40) {
		if (counter >= 20)
			int a=0;
		Mesh::FaceHandle f_h = slPoint.f_h;
		FaceData face = m_pMesh->property(faceDataHandle, f_h);
		if (!face.bVisited) m_NumFacesVisited++;
		face.bVisited = 1;

		OBBoxEx obb;
		StreamLinePoint newSLPoint = slPoint;
		if (!IntegrationStep(newSLPoint, obb, stepSize, dir, previousFace)) return 0;
		obb.streamLineIdx = streamLineIdx;
		obb.oBBoxIdx = m_NumOBBoxes++;

		float closestDist = FLT_MAX;
		for (int i=0, end=face.oBBoxes.size(); i<end; i++) {
			// Do not consider the last obb added as it will certianly be intersecting
			Vector3f newPos = obb.centerW + obb.halfWidthsW.x()*obb.rotVecW[0];
			Vector3f newPos2 = obb.centerW - obb.halfWidthsW.x()*obb.rotVecW[0];
			if (face.oBBoxes[i].oBBoxIdx >= m_NumOBBoxes-3 && face.oBBoxes[i].streamLineIdx == streamLineIdx) continue;
			
			if (face.oBBoxes[i].OBBoxIntersectW(obb)) { // End the streamline
				goto end;
			} else {									// Find the closest box to the box
				//float dist = face.oBBoxes[i].ComputeClosestDistanceAprox(obb);
				//if (dist < closestDist) closestDist = dist;
			}
		}
		
		// Add the obb to the faces
		std::map<unsigned int, bool> visitedFaces;
		//obb.DrawBox();
		AddOBBToFaces(obb, f_h, visitedFaces);
		m_pMesh->property(faceDataHandle, f_h).oBBoxes.push_back(obb);
		m_FacesTouched.push_back(f_h);

		/*glBegin(GL_LINES);
		glColor3f(1,0,0);

		glVertex3f(slPoint.curv.pos[0], slPoint.curv.pos[1], slPoint.curv.pos[2]);
		glVertex3f(newSLPoint.curv.pos[0], newSLPoint.curv.pos[1], newSLPoint.curv.pos[2]);

		glEnd();*/

		m_StreamLines.back().edges.push_back(Edge(slPoint.curv.pos, newSLPoint.curv.pos));

		//seeds.push(Seed());
		slPoint = newSLPoint;
	}
end:
	//return totalStreamLineLen > minLen;
	return 1;
}

void StreamLines::AddOBBToFaces(OBBoxEx & obb, Mesh::FaceHandle & f_h, std::map<unsigned int, bool> & visitedFaces)
{
	if (visitedFaces.find(f_h.idx()) != visitedFaces.end()) return;
	visitedFaces.insert(std::pair<unsigned int, bool>(f_h.idx(), true));

	for (Mesh::FaceFaceIter ff_it = m_pMesh->ff_begin(f_h), end = m_pMesh->ff_end(f_h); ff_it != end; ++ff_it) {
		Mesh::FaceVertexIter fv_it = m_pMesh->fv_begin(ff_it);
		Vec3f v1 = m_pMesh->point(fv_it);
		Vec3f v2 = m_pMesh->point(++fv_it);
		Vec3f v3 = m_pMesh->point(++fv_it);
		// if intersect obb tri
		Vector3f newPos = obb.centerW + obb.halfWidthsW.x()*obb.rotVecW[0];
		Vector3f newPos2 = obb.centerW - obb.halfWidthsW.x()*obb.rotVecW[0];
		if (obb.IntersectTriangle(VT_EIGEN(v1), VT_EIGEN(v2), VT_EIGEN(v3))) {
			m_pMesh->property(faceDataHandle, ff_it).oBBoxes.push_back(obb);
			m_FacesTouched.push_back(ff_it.handle());
			// Find all the adjacent triangles 
			//AddOBBToFaces(obb, ff_it.handle(), visitedFaces);
		}
	}
}

bool StreamLines::WrapToNextFace(StreamLinePoint & slPoint, float & stepSize, float dir, Mesh::FaceHandle & previousFace)
{
	int edgeIdx = 0;
	Vector3f intersectPos;
	// Note that the openmesh edge index is -1 + the standard index system

	if (!TriangleEdgesSegmentIntersect(TRIANGLE(slPoint.f_h), 
		slPoint.curv.pos, slPoint.curv.pos + (dir*stepSize)*slPoint.curv.dir[0], intersectPos, edgeIdx, *m_pMesh, slPoint.f_h, previousFace)) return 0;
	
	Mesh::FaceHandle nextFaceHandle;
	Mesh::FaceFaceIter ff_it = m_pMesh->ff_begin(slPoint.f_h);
	for (int i=0; i<edgeIdx; i++) {
		++ff_it;
	}
	nextFaceHandle = ff_it.handle();

	// Flip the direction vector
	Vec3f n1 = m_pMesh->normal(slPoint.f_h);
	Vec3f n2 = m_pMesh->normal(nextFaceHandle);
	Vector3f m(VT_EIGEN((n1 + n2).normalized()));
	Vector3f a1 = dir * slPoint.curv.dir[0];
	Vector3f a2 = a1 - 2.f*m*(a1.dot(m));

	//std::cout << "Intersection Positions\n\n";
	//PRINTV(slPoint.curv.pos);
	//PRINTV(intersectPos)
	
	// Compute the partial step size, that is, the step size to the intersection position
	float dist = (intersectPos - slPoint.curv.pos).norm();
	//if ((intersectPos - slPoint.curv.pos).norm() > stepSize) 
	//	__debugbreak();
	stepSize = (intersectPos - slPoint.curv.pos).norm();

	// Flip because it will be multipied by dir later
	slPoint.curv.dir[0] = dir*(a2.normalized());
	slPoint.curv.pos = intersectPos;
	slPoint.f_h = nextFaceHandle;

	return 1;
}

bool StreamLines::IntegrationStep(StreamLinePoint & slPoint, OBBoxEx & obb, float & stepSize, float dir, Mesh::FaceHandle & previousFace)
{
	// Triangle of face
	Triangle & tri = TRIANGLE(slPoint.f_h);
	/*Mesh::FaceVertexIter fv_it = m_pMesh->fv_iter(slPoint.f_h);
	CurvatureInfo p1 = CURVATURE(fv_it);
	CurvatureInfo p2 = CURVATURE(++fv_it);
	CurvatureInfo p3 = CURVATURE(++fv_it);
	Vector3f nor(VT_EIGEN(m_pMesh->calc_face_normal(slPoint.f_h)));*/

	// Determine if we are not on an edge
	// Normally, floating point comparison does not work
	if (stepSize == m_StepSize) {
		//Update curv and obb
		// Compute curvature at the slPoint
		/*std::cout << p1.dir[0] << std::endl;
		std::cout << std::endl;

		std::cout << p2.dir[0] << std::endl;
		std::cout << std::endl;

		std::cout << p3.dir[0] << std::endl;
		std::cout << std::endl;*/

		// may be necessary to project this crap to the normal plane
		InterpolateCurvature(tri, slPoint.curv.pos, slPoint.curv);
		ProjectToTrianglePlane(tri, slPoint.curv.dir[0]);
		//PRINTV(slPoint.curv.dir[0]);

		//float d = abs(slPoint.curv.dir[0].dot(slPoint.curv.dir[1]));
		//if (abs(slPoint.curv.dir[0].dot(slPoint.curv.dir[1])) > .05f) {
		//	std::cout << "Bad orthonormal frame: " << d << std::endl;
		//}
	} // else the CurvatureInfo is not updated

	// Compute the new position
	StreamLinePoint newSLPoint = slPoint;
	IntegrateEuler(newSLPoint.curv.pos, dir*newSLPoint.curv.dir[0], stepSize);
	//PRINTV(slPoint.curv.pos);
	//PRINTV(newSLPoint.curv.pos);
	float distance = DistanceToTrianglePlane(tri, newSLPoint.curv.pos);
	// May need to project to the plane
	//std::cout << "Dist: " << distance << std::endl;

	// Check if we are still in the face and compute the allowable step size
	float partialStepSize = stepSize;
	//PRINTV(p1.pos);
	//PRINTV(p2.pos);
	//PRINTV(p3.pos);
	//PRINTV(slPoint.curv.pos);
	bool partialStep = 0;
	if (!PointInTriangle(tri, newSLPoint.curv.pos)) {
		newSLPoint = slPoint;
		if (!WrapToNextFace(newSLPoint, partialStepSize, dir, previousFace)) return 0;
		partialStep = 1;
	}

	// Compute a slighlty larger obb
	const float expansion = 0*partialStepSize;
	obb.ComputeOBB(dir*slPoint.curv.dir[0], tri.n, slPoint.curv.pos, partialStepSize+expansion, m_StreamLineSep);

	//PRINTV(slPoint.curv.pos);
	//std::cout << '\n';

	previousFace = slPoint.f_h;
	slPoint = newSLPoint;
	// partialStepSize may be larger than stepSize in certain situations with WrapToNextFace
	stepSize = partialStep ? abs(stepSize-partialStepSize) : m_StepSize;

	return 1;
}

void StreamLines::DrawStreamLines()
{
	for (int i=0, end=m_StreamLines.size(); i<end; i++) {
		std::vector<Edge> & edges = m_StreamLines[i].edges;
		for (int j=0, end2=m_StreamLines[i].edges.size(); j < end2; j++) {
			glBegin(GL_LINES);
			glColor3f(0,0,0);

			glVertex3f(edges[j].pos[0].x(), edges[j].pos[0].y(), edges[j].pos[0].z());
			glVertex3f(edges[j].pos[1].x(), edges[j].pos[1].y(), edges[j].pos[1].z());

			glEnd();
		}
	}
}

