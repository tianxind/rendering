
#include "OBBox.h"
#include "Triangle.h"
#include "glut.h"

using namespace Eigen;

//For use in the OBBox Intersect Test
#define Abs(x) (((x) < 0 ? -(x) : (x)))
#define _Abs(x, t) ((t=x) < 0 ? -t : t)
#define DotXYZ(a, b) (a(0)*b(0) + a(1)*b(1) + a(2)*b(2))

OBBox::OBBox()
{
}

void OBBox::ComputeOBB(const Vector3f & primaryDir, const Vector3f & nor, const Vector3f & startingPt, float length, float width)
{
	rotVecL[0] = primaryDir;
	rotVecL[2] = nor;
	rotVecL[1] = primaryDir.cross(nor);
	centerL = startingPt + (.5f*length)*rotVecL[0];
	halfWidthsL[0] = .5f*length;
	halfWidthsL[1] = .5f*width;
	halfWidthsL[2] = .5f*width;

	rotVecW[0] = rotVecL[0];
	rotVecW[1] = rotVecL[1];
	rotVecW[2] = rotVecL[2];
	centerW = centerL;
	halfWidthsW = halfWidthsL;
}

bool OBBox::OBBoxIntersectW(OBBox & oBB)
{
	//See Ericson, Christer, Real Time Collision Detection, page 101, and
	//See Gottschalk, Stefan, Collision Queries using Oriented Bounding Boxes

	//a is assumed to be the this OBBox
	//b is assumed to be oBB
	float fRA=0, fRB=0, t=0;
	const float ep = 1e-4f;
	Vector3f * rVA = rotVecW;
	Vector3f * rVB = oBB.rotVecW;
	Vector3f & hWA = halfWidthsW;
	Vector3f & hWB = oBB.halfWidthsW;

	//Rotation matrices in a's coordinate frame
    Matrix4f rot, aRt;
	rot(0, 0) = DotXYZ(rVA[0], rVB[0]); aRt(0, 0) = Abs(rot(0, 0)) + ep;
	rot(0, 1) = DotXYZ(rVA[0], rVB[1]); aRt(0, 1) = Abs(rot(0, 1)) + ep;
	rot(0, 2) = DotXYZ(rVA[0], rVB[2]); aRt(0, 2) = Abs(rot(0, 2)) + ep;

	rot(1, 0) = DotXYZ(rVA[1], rVB[0]); aRt(1, 0) = Abs(rot(1, 0)) + ep;
	rot(1, 1) = DotXYZ(rVA[1], rVB[1]); aRt(1, 1) = Abs(rot(1, 1)) + ep;
	rot(1, 2) = DotXYZ(rVA[1], rVB[2]); aRt(1, 2) = Abs(rot(1, 2)) + ep;

	rot(2, 0) = DotXYZ(rVA[2], rVB[0]); aRt(2, 0) = Abs(rot(2, 0)) + ep;
	rot(2, 1) = DotXYZ(rVA[2], rVB[1]); aRt(2, 1) = Abs(rot(2, 1)) + ep;
	rot(2, 2) = DotXYZ(rVA[2], rVB[2]); aRt(2, 2) = Abs(rot(2, 2)) + ep;

    //Translation vector T in a's coordinate frame
	Vector3f T(oBB.centerW - centerW);
	T = Vector3f(DotXYZ(T, rVA[0]), DotXYZ(T, rVA[1]), DotXYZ(T, rVA[2]));

	//Axis rotVecA[0]
	fRB = hWB(0) * aRt(0, 0) + hWB(1) * aRt(0, 1) + hWB(2) * aRt(0, 2);
	if (Abs(T(0)) > (hWA(0) + fRB)) return 0;

	//Axis rotVecA[1]
	fRB = hWB(0) * aRt(1, 0) + hWB(1) * aRt(1, 1) + hWB(2) * aRt(1, 2);
	if (Abs(T(1)) > (hWA(1) + fRB)) return 0;

	//Axis rotVecA[2]
	fRB = hWB(0) * aRt(2, 0) + hWB(1) * aRt(2, 1) + hWB(2) * aRt(2, 2);
	if (Abs(T(2)) > (hWA(2) + fRB)) return 0;

	//Axis rotVecB[0]
	fRA = hWA(0) * aRt(0, 0) + hWA(1) * aRt(1, 0) + hWA(2) * aRt(2, 0);
	if (_Abs(T(0) * rot(0, 0) + T(1) * rot(1, 0) + T(2) * rot(2, 0), t) > (fRA + hWB(0))) return 0;

	//Axis rotVecB[1]
	fRA = hWA(0) * aRt(0, 1) + hWA(1) * aRt(1, 1) + hWA(2) * aRt(2, 1);
	if (_Abs(T(0) * rot(0, 1) + T(1) * rot(1, 1) + T(2) * rot(2, 1), t) > (fRA + hWB(1))) return 0;

	//Axis rotVecB[2]
	fRA = hWA(0) * aRt(0, 2) + hWA(1) * aRt(1, 2) + hWA(2) * aRt(2, 2);
	if (_Abs(T(0) * rot(0, 2) + T(1) * rot(1, 2) + T(2) * rot(2, 2), t) > (fRA + hWB(2))) return 0;

	//Axes rotVecA[0] X rotVecB[0];
	fRA = hWA(1) * aRt(2, 0) + hWA(2) * aRt(1, 0);
	fRB = hWB(1) * aRt(0, 2) + hWB(2) * aRt(0, 1);
	if (_Abs(T(2) * rot(1, 0) - T(1) * rot(2, 0), t) > (fRA + fRB)) return 0;

    //Axes rotVecA[0] X rotVecB[1];
	fRA = hWA(1) * aRt(2, 1) + hWA(2) * aRt(1, 1);
	fRB = hWB(0) * aRt(0, 2) + hWB(2) * aRt(0, 0);
	if (_Abs(T(2) * rot(1, 1) - T(1) * rot(2, 1), t) > (fRA + fRB)) return 0;

	//Axes rotVecA[0] X rotVecB[2];
	fRA = hWA(1) * aRt(2, 2) + hWA(2) * aRt(1, 2);
	fRB = hWB(0) * aRt(0, 1) + hWB(1) * aRt(0, 0);
	if (_Abs(T(2) * rot(1, 2) - T(1) * rot(2, 2), t) > (fRA + fRB)) return 0;

	//Axes rotVecA[1] X rotVecB[0];
	fRA = hWA(0) * aRt(2, 0) + hWA(2) * aRt(0, 0);
	fRB = hWB(1) * aRt(1, 2) + hWB(2) * aRt(1, 1);
	if (_Abs(T(0) * rot(2, 0) - T(2) * rot(0, 0), t) > (fRA + fRB)) return 0;

	//Axes rotVecA[1] X rotVecB[1];
	fRA = hWA(0) * aRt(2, 1) + hWA(2) * aRt(0, 1);
	fRB = hWB(0) * aRt(1, 2) + hWB(2) * aRt(1, 0);
	if (_Abs(T(0) * rot(2, 1) - T(2) * rot(0, 1), t) > (fRA + fRB)) return 0;

	//Axes rotVecA[1] X rotVecB[2];
	fRA = hWA(0) * aRt(2, 2) + hWA(2) * aRt(0, 2);
	fRB = hWB(0) * aRt(1, 1) + hWB(1) * aRt(1, 0);
	if (_Abs(T(0) * rot(2, 2) - T(2) * rot(0, 2), t) > (fRA + fRB)) return 0;

	//Axes rotVecA[2] X rotVecB[0];
    fRA = hWA(0) * aRt(1, 0) + hWA(1) * aRt(0, 0);
    fRB = hWB(1) * aRt(2, 2) + hWB(2) * aRt(2, 1);
    if (_Abs(T(1) * rot(0, 0) - T(0) * rot(1, 0), t) > fRA + fRB) return 0;

    //Axes rotVecA[2] X rotVecB[1];
    fRA = hWA(0) * aRt(1, 1) + hWA(1) * aRt(0, 1);
    fRB = hWB(0) * aRt(2, 2) + hWB(2) * aRt(2, 0);
    if (_Abs(T(1) * rot(0, 1) - T(0) * rot(1, 1), t) > fRA + fRB) return 0;

    //Axes rotVecA[2] X rotVecB[2];
    fRA = hWA(0) * aRt(1, 2) + hWA(1) * aRt(0, 2);
    fRB = hWB(0) * aRt(2, 1) + hWB(1) * aRt(2, 0);
    if (_Abs(T(1) * rot(0, 2) - T(0) * rot(1, 2), t) > fRA + fRB) return 0;

	return 1;
}

inline void PointInOBB(Vector3f pt, Vector3f hWB, Vector3f * rVB, Vector3f & centerB)
{
	pt -= centerB;
	float cmp1 = DotXYZ(pt, rVB[0]);
	if (cmp1 < -hWB.x() || cmp1 > hWB.x()) return;
	float cmp2 = DotXYZ(pt, rVB[1]);
	if (cmp2 < -hWB.y() || cmp2 > hWB.y()) return;
	float cmp3 = DotXYZ(pt, rVB[2]);
	if (cmp3 < -hWB.z() || cmp3 > hWB.z()) return;
}

float OBBox::ComputeClosestDistanceAprox(OBBox & oBB)
{
	float hX = halfWidthsW.x();
	float hY = halfWidthsW.y();
	float hZ = halfWidthsW.z();
	// Test this in oBB
	PointInOBB(centerW + (+hX*rotVecW[0]+hY*rotVecW[1]+hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (+hX*rotVecW[0]+hY*rotVecW[1]-hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (+hX*rotVecW[0]-hY*rotVecW[1]+hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (+hX*rotVecW[0]-hY*rotVecW[1]-hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	
	PointInOBB(centerW + (-hX*rotVecW[0]+hY*rotVecW[1]+hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (-hX*rotVecW[0]+hY*rotVecW[1]-hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (-hX*rotVecW[0]-hY*rotVecW[1]+hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);
	PointInOBB(centerW + (-hX*rotVecW[0]-hY*rotVecW[1]-hZ*rotVecW[2]), oBB.halfWidthsW, oBB.rotVecW, oBB.centerW);

	// Test oBB in this
	PointInOBB(oBB.centerW + (+hX*rotVecW[0]+hY*rotVecW[1]+hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (+hX*rotVecW[0]+hY*rotVecW[1]-hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (+hX*rotVecW[0]-hY*rotVecW[1]+hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (+hX*rotVecW[0]-hY*rotVecW[1]-hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	
	PointInOBB(oBB.centerW + (-hX*rotVecW[0]+hY*rotVecW[1]+hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (-hX*rotVecW[0]+hY*rotVecW[1]-hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (-hX*rotVecW[0]-hY*rotVecW[1]+hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);
	PointInOBB(oBB.centerW + (-hX*rotVecW[0]-hY*rotVecW[1]-hZ*rotVecW[2]), halfWidthsW, rotVecW, centerW);

	return 1.f;
}

bool OBBox::IntersectTriangle(Vector3f & v1, Vector3f & v2, Vector3f & v3)
{
	Matrix3f rot;
	RotationW(rot);
	Vector3f relV[3] = {rot*(v1-centerW), rot*(v2-centerW), rot*(v3-centerW) };
	float triVert[3][3] = { 
		{relV[0].x(), relV[0].y(), relV[0].z()}, 
		{relV[1].x(), relV[1].y(), relV[1].z()},
		{relV[2].x(), relV[2].y(), relV[2].z()}
	};
	float zero[3] = {0};
	float halfWidths[3] = {halfWidthsW.x(), halfWidthsW.y(), halfWidthsW.z()};
	return triBoxOverlap(zero, halfWidths, triVert);

	//return triAABBOverlap((float*)
	//triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float * triverts[3])
	//return triAABBOverlap(Vector3f(0,0,0), halfWidthsW, relV);
}

void OBBox::DrawBox()
{

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  /* A step backward, then spin the cube */
  glTranslatef(centerW.x(), centerW.y(), centerW.z());
  float mat4x4[] = {rotVecW[0].x(), rotVecW[0].y(), rotVecW[0].z(), 0,
				    rotVecW[1].x(), rotVecW[1].y(), rotVecW[1].z(), 0,
				    rotVecW[2].x(), rotVecW[2].y(), rotVecW[2].z(), 0,
					0,			  0,			  0,			  1 };
  glMultMatrixf((float*)mat4x4);
  glScalef(halfWidthsW.x(),halfWidthsW.y(),halfWidthsW.z());
  //glRotatef(alpha, 0, 1, 0);

  /* We tell we want to draw quads */
  glBegin(GL_QUADS);

  /* Every four calls to glVertex, a quad is drawn */
  glColor3f(0, 0, 0); glVertex3f(-1, -1, -1);
  glColor3f(0, 0, 1); glVertex3f(-1, -1,  1);
  glColor3f(0, 1, 1); glVertex3f(-1,  1,  1);
  glColor3f(0, 1, 0); glVertex3f(-1,  1, -1);

  glColor3f(1, 0, 0); glVertex3f( 1, -1, -1);
  glColor3f(1, 0, 1); glVertex3f( 1, -1,  1);
  glColor3f(1, 1, 1); glVertex3f( 1,  1,  1);
  glColor3f(1, 1, 0); glVertex3f( 1,  1, -1);

  glColor3f(0, 0, 0); glVertex3f(-1, -1, -1);
  glColor3f(0, 0, 1); glVertex3f(-1, -1,  1);
  glColor3f(1, 0, 1); glVertex3f( 1, -1,  1);
  glColor3f(1, 0, 0); glVertex3f( 1, -1, -1);

  glColor3f(0, 1, 0); glVertex3f(-1,  1, -1);
  glColor3f(0, 1, 1); glVertex3f(-1,  1,  1);
  glColor3f(1, 1, 1); glVertex3f( 1,  1,  1);
  glColor3f(1, 1, 0); glVertex3f( 1,  1, -1);

  glColor3f(0, 0, 0); glVertex3f(-1, -1, -1);
  glColor3f(0, 1, 0); glVertex3f(-1,  1, -1);
  glColor3f(1, 1, 0); glVertex3f( 1,  1, -1);
  glColor3f(1, 0, 0); glVertex3f( 1, -1, -1);

  glColor3f(0, 0, 1); glVertex3f(-1, -1,  1);
  glColor3f(0, 1, 1); glVertex3f(-1,  1,  1);
  glColor3f(1, 1, 1); glVertex3f( 1,  1,  1);
  glColor3f(1, 0, 1); glVertex3f( 1, -1,  1);

  glEnd();

  glPopMatrix();
}

OBBox & OBBox::operator=(OBBox & ref)
{
	memcpy(this, &ref, sizeof(OBBox));
	return *this;
}


