#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GL/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
#include "mesh_features.h"
#include "image_generation.h"
#include "decimate.h"
using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<double> viewCurvature;
FPropHandleT<Vec3f> viewCurvatureDerivative;
VPropHandleT<CurvatureInfo> curvature;
Mesh mesh;

bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = {0,0,2,1};
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
bool showSurface = true, showAxes = true, showCurvature = false, showNormals = false;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float shininess[] = { 50.0 };

/* Returns 1 if kw is postivie and -1 if kw is negative */
int KwSign(double kw){
  if(kw < 0) return -1;
  else return 1;
}

/* If Zero crossing exists, stores pairs of vertices to 
 * interpolate between in double* interp and returns true.
 * Otherwise, returns false.
 */
bool ZeroCrossingExists(double* kw, size_t* interp){
  size_t index = 0;
  for(size_t i = 1; i < 3; ++i){
    if(KwSign(kw[i]) != KwSign(kw[i-1])){
      interp[index] = i-1;
      interp[index+1] = i;
      index += 2;
    } 
  }
  if(index == 2) {
    interp[index] = 2;
    interp[index+1] = 0;
  }
  return index != 0;
}

/* Returns Dwkw where gradient of kw is multiplied by
 * w interpolated at the centroid of the face
 */
double ComputeDirectionalGradient(Mesh::FaceIter f_It, Vec3f* vertex,
                                  Vec3f camPos){
  Vec3f kw_gradient = mesh.property(viewCurvatureDerivative,f_It.handle());
  Vec3f centroid = vertex[0] + vertex[1] + vertex[2]; 
  centroid /= 3;

  Vec3f v = camPos - centroid;
  Vec3f n = mesh.normal(f_It.handle());
  Vec3f w = v - n*(n|v);
  
  return kw_gradient|w;
}

/* Interpolates the point coordinates where kw = 0 */
Vec3f InterpZeroCrossingPt(Vec3f* vertex, double* kw, int i1, int i2){
  double t = -1 * kw[i1] / (kw[i2]-kw[i1]);
  return vertex[i1]*(1-t) + vertex[i2]*t;
}

/* Finds Zero Crossing and stores their coordinates in Vec3f* zero_x */
bool FindZeroCrossings(Vec3f* vertex, double* kw, Vec3f* zero_x, 
                       Mesh::FaceIter f_It, Vec3f actualCamPos){
  size_t interp[4];                  // stores pairs of vertex indices
  if(ZeroCrossingExists(kw, interp)){
    double dwkw = ComputeDirectionalGradient(f_It, vertex, actualCamPos);
    if(dwkw > 0){
      Vec3f p1 = InterpZeroCrossingPt(vertex, kw, interp[0], interp[1]);
      Vec3f p2 = InterpZeroCrossingPt(vertex, kw, interp[2], interp[3]);
      zero_x[0] = p1;
      zero_x[1] = p2;
      return true;
    }
  } 
  return false;
}

void renderSuggestiveContours(Vec3f actualCamPos) { // use this camera position to account for panning etc.
  glColor3f(.5,.5,.5);

  // --------------ALGORITHM------------------
  // for each face
  //     interpolate kw  between each pair of vertices: v0->v1 v1->v2 v2->v0
  //          if(kw == 0 AND Dwkw > 0 for any two pts along these edges)
  //               connect these points
  // NOTE:  Dwkw = viewCurvatureDerivative * w 
  // QUESTION: what is w? do we interpolate it?
  // -----------------------------------------

  for(Mesh::FaceIter it = mesh.faces_begin(); it !=mesh.faces_end(); ++it){
    Mesh::ConstFaceVertexIter fvIt = mesh.cfv_iter(it);
    Vec3f vertex[3];           // each vertex
    double kw[3];              // kw's of each vertex
    Vec3f zero_x[2];           // pts to connect on 2 diff edges of a triangle
    bool zero_x_found = false; // whether 2 pts have been found

    for(int i = 0; i < 3; i++){
      vertex[i] = mesh.point(fvIt.handle());
      kw[i] = mesh.property(viewCurvature, fvIt.handle());
      ++fvIt;
    }
          
    // Find points on different edges where kw = 0 & dwkw >0
    zero_x_found = FindZeroCrossings(vertex, kw, zero_x, it, actualCamPos);
    Vec3f pt1 = zero_x[0];
    Vec3f pt2 = zero_x[1];

    // Connect these points
    if(zero_x_found){
      glBegin(GL_LINES);
      glVertex3f(pt1[0], pt1[1], pt1[2]);
      glVertex3f(pt2[0], pt2[1], pt2[2]);
      glEnd();
     }
  }
}

void drawTriangles()
{
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
		Mesh::ConstFaceVertexIter cfvIt;
		cfvIt = mesh.cfv_iter(f_it.handle());

		glBegin(GL_TRIANGLES);
		//glColor3f(0.1, 0.2, 0.3);
 
		Vec3f pointA = mesh.point(cfvIt.handle());
		Vec3f nA = mesh.normal(cfvIt.handle());
		glNormal3d(nA[0], nA[1], nA[2]);
		glVertex3f(pointA[0], pointA[1], pointA[2]);

		Vec3f pointB = mesh.point((++cfvIt).handle());
		Vec3f nB = mesh.normal(cfvIt.handle());
		glNormal3d(nB[0], nB[1], nB[2]);
		glVertex3f(pointB[0], pointB[1], pointB[2]);

		Vec3f pointC = mesh.point((++cfvIt).handle());
		Vec3f nC = mesh.normal(cfvIt.handle());
		glNormal3d(nC[0], nC[1], nC[2]);
		glVertex3f(pointC[0], pointC[1], pointC[2]);

		glEnd();
	}	
}

void renderMesh() {
	if (!showSurface) glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); // render regardless to remove hidden lines
	
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);

	glDepthRange(0.001,1);
	glEnable(GL_NORMALIZE);
	
	// draw the filled polygons
    glLineWidth(1.1f);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glEnable( GL_POLYGON_OFFSET_FILL );
    glPolygonOffset(1, 1);
    drawTriangles();
    glDisable( GL_POLYGON_OFFSET_FILL );
	
    // draw the wireframe
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);
	drawTriangles();
   	glEnable(GL_LIGHTING);
	
	if (!showSurface) glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	
	glDisable(GL_LIGHTING);
	glDepthRange(0,0.999);
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
		renderSuggestiveContours(actualCamPos);
	
	// We'll be nice and provide you with code to render feature edges below
	glBegin(GL_LINES);
	glColor3f(0,0,0);
	glLineWidth(2.0f);
	for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
		if (isFeatureEdge(mesh,*it,actualCamPos)) {
			Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
			Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
			Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
			Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));
			glVertex3f(source[0],source[1],source[2]);
			glVertex3f(target[0],target[1],target[2]);
		}
	glEnd();
	
		if (showCurvature) {
		  glBegin(GL_LINES);
		  glColor3f(0,1,0);
		  for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
			CurvatureInfo info = mesh.property(curvature, it.handle());
			Vec3f p = mesh.point(it.handle());
			Vec3f e1 = info.directions[0];
			Vec3f e2 = info.directions[1];
			
			float vecLen = .01f;
			Vec3f d1Forward = p + e1 * vecLen;
			Vec3f d1Back = p - e1 * vecLen;
			glVertex3f(d1Forward[0], d1Forward[1], d1Forward[2]);
			glVertex3f(d1Back[0], d1Back[1], d1Back[2]);
		
			Vec3f d2Forward = p + e2 * vecLen;
			Vec3f d2Back = p - e2 * vecLen;
			glVertex3f(d2Forward[0], d2Forward[1], d2Forward[2]);
			glVertex3f(d2Back[0], d2Back[1], d2Back[2]);
		  }
		  glEnd();
		}
	
	if (showNormals) {
		glBegin(GL_LINES);
		glColor3f(0,1,0);
		for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
			Vec3f n = mesh.normal(it.handle());
			Vec3f p = mesh.point(it.handle());
			Vec3f d = p + n*.01;
			glVertex3f(p[0],p[1],p[2]);
			glVertex3f(d[0],d[1],d[2]);
		}
		glEnd();
	}
	
	glDepthRange(0,1);
}

void display() {
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
	glEnable(GL_LIGHT0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,windowWidth,windowHeight);
	
	float ratio = (float)windowWidth / (float)windowHeight;
	gluPerspective(50, ratio, .01, 10); // 50 degree vertical viewing angle, zNear = 1, zFar = 1000
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPos[0]+pan[0], cameraPos[1]+pan[1], cameraPos[2]+pan[2], pan[0], pan[1], pan[2], up[0], up[1], up[2]);
	
	// Draw mesh
	renderMesh();

	// Draw axes
	if (showAxes) {
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glLineWidth(1);
			glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(1,0,0); // x axis
			glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,1,0); // y axis
			glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,1); // z axis
		glEnd(/*GL_LINES*/);
	}
	
	glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
	else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
	else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);
	
	lastPos[0] = x;
	lastPos[1] = y;
}

void mouseMoved(int x, int y) {
	int dx = x - lastPos[0];
	int dy = y - lastPos[1];
	Vec3f curCamera(cameraPos[0],cameraPos[1],cameraPos[2]);
	Vec3f curCameraNormalized = curCamera.normalized();
	Vec3f right = up % curCameraNormalized;

	if (leftDown) {
		// Assume here that up vector is (0,1,0)
		Vec3f newPos = curCamera - 2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
		newPos = newPos.normalized() * curCamera.length();
		
		up = up - (up | newPos) * newPos / newPos.sqrnorm();
		up.normalize();
		
		for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
	}
	else if (rightDown) for (int i = 0; i < 3; i++) cameraPos[i] /= pow(1.1,dy*.1);
	else if (middleDown) {
		pan += -2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
	}

	
	lastPos[0] = x;
	lastPos[1] = y;
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
	computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);
	
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);

	if (key == 's' || key == 'S') showSurface = !showSurface;
	else if (key == 'a' || key == 'A') showAxes = !showAxes;
	else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
	else if (key == 'n' || key == 'N') showNormals = !showNormals;
	else if (key == 'w' || key == 'W') writeImage(mesh, windowWidth, windowHeight, "renderedImage.svg", actualCamPos);
	else if (key == 'q' || key == 'Q') exit(0);
	glutPostRedisplay();
}

void reshape(int width, int height) {
	windowWidth = width;
	windowHeight = height;
	glutPostRedisplay();
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " mesh_filename\n";
		exit(0);
	}
	
	IO::Options opt;
	opt += IO::Options::VertexNormal;
	opt += IO::Options::FaceNormal;
	
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	
	cout << "Reading from file " << argv[1] << "...\n";
	if ( !IO::read_mesh(mesh, argv[1], opt )) {
		cout << "Read failed.\n";
		exit(0);
	}

	cout << "Mesh stats:\n";
	cout << '\t' << mesh.n_vertices() << " vertices.\n";
	cout << '\t' << mesh.n_edges() << " edges.\n";
	cout << '\t' << mesh.n_faces() << " faces.\n";
	
	//simplify(mesh,.1f);
	
	mesh.update_normals();
	
	mesh.add_property(viewCurvature);
	mesh.add_property(viewCurvatureDerivative);
	mesh.add_property(curvature);
	
	// Move center of mass to origin
	Vec3f center(0,0,0);
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) center += mesh.point(vIt);
	center /= mesh.n_vertices();
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) -= center;

	// Fit in the unit sphere
	float maxLength = 0;
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) maxLength = max(maxLength, mesh.point(vIt).length());
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) /= maxLength;
	
	computeCurvature(mesh,curvature);

	up = Vec3f(0,1,0);
	pan = Vec3f(0,0,0);
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
	computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);

	glutInit(&argc, argv); 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 
	glutInitWindowSize(windowWidth, windowHeight); 
	glutCreateWindow(argv[0]);

	glutDisplayFunc(display);
	glutMotionFunc(mouseMoved);
	glutMouseFunc(mouse);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);

	glutMainLoop();
	
	return 0;
}
