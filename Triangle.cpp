
#include "Triangle.h"

#define Abs(x) ((x < 0 ? -x : x))

float DistanceToTrianglePlane(const Vector3f & p1, const Vector3f & p2, const Vector3f & p3, const Vector3f & point)
{
	Vector3f v1(p2-p1);
	Vector3f v2(p3-p1);
	Vector3f n(v1.cross(v2));
	n.normalize();
	float d = -(n.dot(p1));

	return (n.dot(point) + d);
}

bool PointInTriangle(const Vector3f & p1, const Vector3f & p2, const Vector3f & p3, Vector3f & point)
{
	Vector3f a(p1 - point);
	Vector3f b(p2 - point);
	Vector3f c(p3 - point);

	Vector3f u = b.cross(c);
	Vector3f v = c.cross(a);

	if (u.dot(v) < 0) return 0;

	Vector3f w = a.cross(b);

	if (u.dot(w) < 0) return 0;

	return 1;
}

float Clamp(float n, float min, float max) {
	if (n < min) return min;
	if (n > max) return max;
	return n;
}

bool FindClosestPoint(const Vector3f & p1, const Vector3f & q1, const Vector3f & p2, const Vector3f & q2, Vector3f & point)
{
	float s,t;
	const float epsUpper = 1+1e-2;
	const float epsLower = -1e-2;

	Vector3f d1(q1-p1);
	Vector3f d2(q2-p2);
	Vector3f r(p1-p2);
	float a = d1.dot(d1);
	float e = d2.dot(d2);
	float f = d2.dot(r);
	float c = d1.dot(r);
	
	float b = d1.dot(d2);
	float denom = a*e-b*b;

	if (denom != 0)
		s = (b*f - c*e)/denom;
	else s = 0;
	if (s < epsLower || s > epsUpper) return 0;

	t = (b*s + f)/e;
	if (t < epsLower) {
		return 0;
		
		//t = 0;
		//s = Clamp(-c/a, 0, 1.f);
	} else if (t > epsUpper) {
		return 0;

		//t = 1.f;
		//s = Clamp((b-c)/a, 0.f, 1.f);
	}

	point = p2 + d2*t;

	return 1;
}

bool TriangleEdgesSegmentIntersect(const Vector3f & p1, const Vector3f & p2, const Vector3f p3, const Vector3f & e0, const Vector3f & e1,
	Vector3f & intersectionPos, int & edgeIdx, Mesh & mesh, Mesh::FaceHandle & currentFace, Mesh::FaceHandle & previousFace)
{
	int r[3] = {0};
	//float d1 = DistanceToTrianglePlane(p1, p2, p3, e0);
	//float d2 = DistanceToTrianglePlane(p1, p2, p3, e1);
	//Vector3f dif = e1 - e0;
	//std::cout << '\n' << std::endl;
	//PRINTV(e0);
	//PRINTV(e1);
	//PRINTV(p1);
	//PRINTV(p2);
	//PRINTV(p3);
	Vector3f interPos[3];
	r[0] = FindClosestPoint(e0, e1, p3, p1, interPos[0]);
	r[1] = FindClosestPoint(e0, e1, p1, p2, interPos[1]);
	r[2] = FindClosestPoint(e0, e1, p2, p3, interPos[2]);
	if (r[0]+r[1]+r[2] == 3) {
		__debugbreak();	//Triangle too small
	}
	if (r[0]+r[1]+r[2] == 2) {
		edgeIdx = 0;
		/*for (Mesh::FaceFaceIter ff_it=mesh.ff_begin(currentFace), end=mesh.ff_end(currentFace); ff_it != end; ++ff_it) {

			Mesh::FaceHandle ffff = ff_it.handle();
			int nodthing=0;
		}*/
		for (Mesh::FaceFaceIter ff_it=mesh.ff_begin(currentFace), end=mesh.ff_end(currentFace); ff_it != end; ++ff_it, ++edgeIdx) {
			// SERIOUS PROBLEM IF THIS HAPPENS
			// Never should happen
			if (edgeIdx >= 3 || edgeIdx < 0)
				__debugbreak();

			//Mesh::FaceHandle ffff = ff_it.handle();
			if (r[edgeIdx] != 0 && ff_it.handle() != previousFace)  break;
		}
		intersectionPos = interPos[edgeIdx];
		return 1;
	}
	if (r[0]+r[1]+r[2] == 0)
		return 0;//__debugbreak(); //Never should happen

	edgeIdx = 0*r[0] + 1*r[1] + 2*r[2];
	intersectionPos = interPos[edgeIdx];
	
	return 1;
}

float DistanceToTrianglePlane(const Triangle & tri, const Vector3f & point) 
{
	return DistanceToTrianglePlane(tri.v[0].pos, tri.v[1].pos, tri.v[2].pos, point);
}

bool PointInTriangle(Triangle & tri, Vector3f & point)
{
	return PointInTriangle(tri.v[0].pos, tri.v[1].pos, tri.v[2].pos, point);
}

bool TriangleEdgesSegmentIntersect(const Triangle & tri, const Vector3f & e0, const Vector3f & e1, 
	Vector3f & intersectionPos, int & edgeIdx, Mesh & mesh, Mesh::FaceHandle & currentFace, Mesh::FaceHandle & previousFace)
{
	return TriangleEdgesSegmentIntersect(tri.v[0].pos, tri.v[1].pos, tri.v[2].pos, e0, e1, intersectionPos, edgeIdx, mesh, currentFace, previousFace);
}














		/*Vector3f p[3] = {p1,p2,p3};
		Vector3f e[3] = {(p1-p3).normalized(), (p2-p1).normalized(), (p3-p2).normalized() };
		Vector3f nor(e[0].cross(-e[2]));
		nor.normalize();

		int edgesIdx[2] = {r[0] ? 0 : 1, r[2] ? 2 : 1 };
		// Which edge is in the direction of the edge ray
		Vector3f v(e1-e0);
		int signs = 1*(e[edgesIdx[0]].dot(v) > 0) + 2*(e[edgesIdx[1]].dot(v) > 0) + 3*(e[0].dot(v) > 0);

		if (signs == 0 || signs >= 3)
			__debugbreak();
		else edgeIdx = edgesIdx[signs-1];
		return 1;*/
































/********************************************************/

/* AABB-triangle overlap test code                      */

/* by Tomas Akenine-Möller                              */

/* Function: int triBoxOverlap(float boxcenter[3],      */

/*          float boxhalfsize[3],float triverts[3][3]); */

/* History:                                             */

/*   2001-03-05: released the code in its first version */

/*   2001-06-18: changed the order of the tests, faster */

/*                                                      */

/* Acknowledgement: Many thanks to Pierre Terdiman for  */

/* suggestions and discussions on how to optimize code. */

/* Thanks to David Hunt for finding a ">="-bug!         */

/********************************************************/

#include <math.h>

#include <stdio.h>

#include "Triangle.h"

#define X 0

#define Y 1

#define Z 2

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 



#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])



#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 



#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;



int planeBoxOverlap(float normal[3], float vert[3], float maxbox[3])	// -NJMP-

{

  int q;

  float vmin[3],vmax[3],v;

  for(q=X;q<=Z;q++)

  {

    v=vert[q];					// -NJMP-

    if(normal[q]>0.0f)

    {

      vmin[q]=-maxbox[q] - v;	// -NJMP-

      vmax[q]= maxbox[q] - v;	// -NJMP-

    }

    else

    {

      vmin[q]= maxbox[q] - v;	// -NJMP-

      vmax[q]=-maxbox[q] - v;	// -NJMP-

    }

  }

  if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-

  if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-

  

  return 0;

}





/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;



/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;



bool triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3])

{



  /*    use separating axis theorem to test overlap between triangle and box */

  /*    need to test for overlap in these directions: */

  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */

  /*       we do not even need to test these) */

  /*    2) normal of the triangle */

  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */

  /*       this gives 3x3=9 more tests */

   float v0[3],v1[3],v2[3];

//   float axis[3];

   float min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed

   float normal[3],e0[3],e1[3],e2[3];



   /* This is the fastest branch on Sun */

   /* move everything so that the boxcenter is in (0,0,0) */

   SUB(v0,triverts[0],boxcenter);

   SUB(v1,triverts[1],boxcenter);

   SUB(v2,triverts[2],boxcenter);



   /* compute triangle edges */

   SUB(e0,v1,v0);      /* tri edge 0 */

   SUB(e1,v2,v1);      /* tri edge 1 */

   SUB(e2,v0,v2);      /* tri edge 2 */



   /* Bullet 3:  */

   /*  test the 9 tests first (this was faster) */

   fex = fabsf(e0[X]);

   fey = fabsf(e0[Y]);

   fez = fabsf(e0[Z]);

   AXISTEST_X01(e0[Z], e0[Y], fez, fey);

   AXISTEST_Y02(e0[Z], e0[X], fez, fex);

   AXISTEST_Z12(e0[Y], e0[X], fey, fex);



   fex = fabsf(e1[X]);

   fey = fabsf(e1[Y]);

   fez = fabsf(e1[Z]);

   AXISTEST_X01(e1[Z], e1[Y], fez, fey);

   AXISTEST_Y02(e1[Z], e1[X], fez, fex);

   AXISTEST_Z0(e1[Y], e1[X], fey, fex);



   fex = fabsf(e2[X]);

   fey = fabsf(e2[Y]);

   fez = fabsf(e2[Z]);

   AXISTEST_X2(e2[Z], e2[Y], fez, fey);

   AXISTEST_Y1(e2[Z], e2[X], fez, fex);

   AXISTEST_Z12(e2[Y], e2[X], fey, fex);



   /* Bullet 1: */

   /*  first test overlap in the {x,y,z}-directions */

   /*  find min, max of the triangle each direction, and test for overlap in */

   /*  that direction -- this is equivalent to testing a minimal AABB around */

   /*  the triangle against the AABB */



   /* test in X-direction */

   FINDMINMAX(v0[X],v1[X],v2[X],min,max);

   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;



   /* test in Y-direction */

   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);

   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;



   /* test in Z-direction */

   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);

   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;



   /* Bullet 2: */

   /*  test if the box intersects the plane of the triangle */

   /*  compute plane equation of triangle: normal*x+d=0 */

   CROSS(normal,e0,e1);

   // -NJMP- (line removed here)

   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-



   return 1;   /* box and triangle overlaps */

}

/*bool triAABBOverlap(Vector3f & boxcenter, Vector3f & boxhalfsize, Vector3f * triverts)
{
	return triBoxOverlap((float*)&boxcenter, (float*)&boxhalfsize, (float**)&triverts);
}*/


