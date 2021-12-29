#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <tuple>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

Polyhedron* poly;
std::vector<PolyLine> lines;
std::vector<icVector3> points;

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = tranlate y, 2 = rotate
int N = 30;

// IBFV related variables (Van Wijk 2002)
//https://www.win.tue.nl/~vanwijk/ibfv/
#define NPN		64
#define SCALE	4.0
#define ALPHA	8
float tmax = win_width / (SCALE * NPN);
float dmax = SCALE / win_width;
unsigned char* pixels;

// Tuples
std::tuple<Quad**, bool> quads_sorted_by_x;
std::tuple<Quad**, bool> quads_sorted_by_y;

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);
void initIBFV();

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void displayIBFV();
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*functions for element picking*/
void display_vertices(GLenum mode, Polyhedron* poly);
void display_quads(GLenum mode, Polyhedron* poly);
void display_selected_vertex(Polyhedron* poly);
void display_selected_quad(Polyhedron* poly);

/*display vis results*/
void display_polyhedron(Polyhedron* poly);

icVector3 grad_at_point(Quad quad, float x, float y);
float bilinear_interpolate_at_point(Quad* quad, float x, float y);
Quad get_quad_by_position(float x, float y, int index_xlist_last_quad = 0, int index_ylist_last_quad = 0); // Gets quad and sorts list
Vertex* find_with_coords(double x, double y, int vert_count, Vertex** vlist);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[])
{
	/*load mesh from ply file*/
	FILE* this_file = fopen("../data/scalar_data/r4.ply", "r");
	poly = new Polyhedron(this_file);
	fclose(this_file);
	
	/*initialize the mesh*/
	poly->initialize(); // initialize the mesh
	poly->write_info();


	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Scientific Visualization");


	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	poly->finalize();	// finalize everything
	free(pixels);
	return 0;
}


/******************************************************************************
Find with coordinates
******************************************************************************/
Vertex* find_with_coords(double x, double y, int v_count, Vertex** vlist) {
	for (int i = 0; i < v_count; i++) {
		if (vlist[i]->x == x && vlist[i]->y == y) {
			return vlist[i];
		}
	}
	return NULL;
}


/******************************************************************************
Gradient at Point
******************************************************************************/
icVector3 grad_at_point(Quad quad, float x, float y)
{


}

/******************************************************************************
Bilinear Interpolation at a Point
******************************************************************************/
float bilinear_interpolate_at_point(Quad* quad, float x0, float y0)
{
		Vertex* v1 = (*quad).verts[0];
		Vertex* v2 = quad->verts[1];
		Vertex* v3 = quad->verts[2];
		Vertex* v4 = quad->verts[3];

		double x1 = std::min(std::min(v3->x, v1->x), std::min(v2->x, v4->x));
		double x2 = std::max(std::max(v3->x, v1->x), std::max(v2->x, v4->x));
		double y1 = std::min(std::min(v3->y, v1->y), std::min(v2->y, v4->y));
		double y2 = std::max(std::max(v3->y, v1->y), std::max(v2->y, v4->y));

		Vertex* x1y1 = find_with_coords(x1, y1, 4, quad->verts);
		Vertex* x2y1 = find_with_coords(x2, y1, 4, quad->verts);
		Vertex* x1y2 = find_with_coords(x1, y2, 4, quad->verts);
		Vertex* x2y2 = find_with_coords(x2, y2, 4, quad->verts);

		double fx1y1 = x1y1->scalar;
		double fx2y1 = x2y1->scalar;
		double fx1y2 = x1y2->scalar;
		double fx2y2 = x2y2->scalar;

		if (x1 <= x0 && x2 >= x0) {
			if (y1 <= y0 && y2 >= y0) {
				double p1 = ((x2 - x0) / (x2 - x1)) * ((y2 - y0)/(y2 - y1));
				double p2 = ((x0-x1)/(x2 - x1))*((y2 - y0)/(y2 - y1));
				double p3 = ((x2 - x0)/(x2 - x1)) * ((y0 - y1)/(y2 - y1));
				double p4 = ((x0- x1)/(x2 - x1)) * ((y0 - y1)/(y2 - y1));
				return (p1 * fx1y1 + p2 * fx2y1 + p3 * fx1y2 + p4 * fx2y2);//TODO: scale at the end. *HEIGHT_AMPLIFY*max_bounds;
			}
		} else {
			return NULL;
		}
}

/******************************************************************************
Get Quad by Position
******************************************************************************/
Quad get_quad_by_position(float x, float y, int index_xlist_last_quad = 0, int index_ylist_last_quad = 0)
{

}
/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode)
{
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	glScalef(0.9 / poly->radius, 0.9 / poly->radius, 0.9 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Initialize IBFV patterns
******************************************************************************/

void initIBFV()
{
	pixels = (unsigned char*)malloc(sizeof(unsigned char) * win_width * win_height * 3);
	memset(pixels, 255, sizeof(unsigned char) * win_width * win_height * 3);

	tmax = win_width / (SCALE * NPN);
	dmax = SCALE / win_width;

	int lut[256];
	int phase[NPN][NPN];
	GLubyte pat[NPN][NPN][4];
	int i, j, k;

	for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
	for (i = 0; i < NPN; i++)
		for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256;

	for (i = 0; i < NPN; i++)
	{
		for (j = 0; j < NPN; j++)
		{
			pat[i][j][0] =
				pat[i][j][1] =
				pat[i][j][2] = lut[(phase[i][j]) % 255];
			pat[i][j][3] = ALPHA;
		}
	}
	
	glNewList(1, GL_COMPILE);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
	glEndList();
}

/******************************************************************************
Pick objects from the scene
******************************************************************************/

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	ptr = (GLuint*)buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	return seed_id;
}

/******************************************************************************
Diaplay all quads for selection
******************************************************************************/

void display_quads(GLenum mode, Polyhedron* this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	for (i = 0; i < this_poly->nquads; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		Quad* temp_q = this_poly->qlist[i];
		
		glBegin(GL_POLYGON);
		for (j = 0; j < 4; j++) {
			Vertex* temp_v = temp_q->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

/******************************************************************************
Diaplay all vertices for selection
******************************************************************************/

void display_vertices(GLenum mode, Polyhedron* this_poly)
{
	for (int i = 0; i < this_poly->nverts; i++) {
		if (mode == GL_SELECT)
			glLoadName(i + 1);

		CHECK_GL_ERROR();

		Vertex* temp_v = this_poly->vlist[i];
		drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15);
	}
	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay selected quad
******************************************************************************/

void display_selected_quad(Polyhedron* this_poly)
{
	if (this_poly->selected_quad == -1)
	{
		return;
	}

	unsigned int i, j;

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	Quad* temp_q = this_poly->qlist[this_poly->selected_quad];

	glBegin(GL_POLYGON);
	for (j = 0; j < 4; j++) {
		Vertex* temp_v = temp_q->verts[j];
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(temp_v->x, temp_v->y, temp_v->z);
	}
	glEnd();
}

/******************************************************************************
Diaplay selected vertex
******************************************************************************/

void display_selected_vertex(Polyhedron* this_poly)
{
	if (this_poly->selected_vertex == -1)
	{
		return;
	}

	Vertex* temp_v = this_poly->vlist[this_poly->selected_vertex];
	drawDot(temp_v->x, temp_v->y, temp_v->z, 0.15, 1.0, 0.0,0.0);

	CHECK_GL_ERROR();
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

/******************************************************************************
Min and max
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	// clear out lines and points
	lines.clear();
	points.clear();

	switch (key) { 
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	// solid color display with lighting
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// wireframe display
		display_mode = 2;
		glutPostRedisplay();
		break;

	case '3':	// checkerboard display
	{
		display_mode = 3;

		double L = (poly->radius * 2) / 30;
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			for (int j = 0; j < 4; j++) {

				Vertex* temp_v = temp_q->verts[j];

				temp_v->R = int(temp_v->x / L) % 2 == 0 ? 1 : 0;
				temp_v->G = int(temp_v->y / L) % 2 == 0 ? 1 : 0;
				temp_v->B = 0.0;
			}
		}
		glutPostRedisplay();
	}
	break;

	case '4':	// Drawing points and lines created by the dots_and_lines_example() function
		display_mode = 4;
		dots_and_lines_example(&points, &lines);
		glutPostRedisplay();
		break;

	case '5':	// IBFV vector field display
		display_mode = 5;
		glutPostRedisplay();
		break;

	case '6':	// Problem 2a -> Displaying contours
		display_mode = 6;
		{
			//Finding min and max scalars 
			double N_scalar = (poly->max_scalar - poly->min_scalar) / N;
			double contour = poly->min_scalar;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < poly->nverts; j++) {
					Vertex* vtemp = poly->vlist[j];
					if (vtemp->scalar >= contour) {
						vtemp->level = 1;
					}
					else {
						vtemp->level = 0;
					}
				}
				for (int i = 0; i < poly->nedges; i++) {
					Edge* etemp = poly->elist[i];
					Vertex* v1 = etemp->verts[0];
					Vertex* v2 = etemp->verts[1];

					// delete any old crossing points
					if (etemp->crossing) {
						delete etemp->crossing;
						etemp->crossing = NULL;
					}

					if (v1->level != v2->level) {
						double scalar1 = v1->scalar;
						double scalar2 = v2->scalar;
						
						double Alpha = (contour - scalar1) / (scalar2 - scalar1);
						double cross_x = Alpha * (v2->x - v1->x) + v1->x;
						double cross_y = Alpha * (v2->y - v1->y) + v1->y;
						double cross_z = Alpha * (v2->z - v1->z) + v1->z;

						icVector3* cross_point = new icVector3(cross_x, cross_y, cross_z);
						etemp->crossing = cross_point;
					}
				}
				// Getting edge crossing points

				PolyLine N_contour;

				for (int i = 0; i < poly->nquads; i++) {
					Quad* qtemp = poly->qlist[i];

					std::vector<icVector3*> crossing_points;
					for (int j = 0; j < 4; j++) {
						Edge* etemp = qtemp->edges[j];
						if (etemp->crossing) {
							crossing_points.push_back(etemp->crossing);
						}
					}
					
					switch (crossing_points.size()) {
					case 2: {
						LineSegment line(*crossing_points[0], *crossing_points[1]);
						N_contour.push_back(line);  // find out how to get things to show up correctly 
					}
					case 4: {
						break;
					}
					}
				}
				lines.push_back(N_contour);

				contour += N_scalar;
			}
			glutPostRedisplay();

			break;
		}
	case '7': //Problem 2b -> displaying contours with solid colors
		display_mode = 7;
		{
			//Finding min and max scalars 
			double N_scalar = (poly->max_scalar - poly->min_scalar) / N;
			double contour = poly->min_scalar;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < poly->nverts; j++) {
					Vertex* vtemp = poly->vlist[j];
					if (vtemp->scalar >= contour) {
						vtemp->level = 1;
					}
					else {
						vtemp->level = 0;
					}
				}
				for (int i = 0; i < poly->nedges; i++) {
					Edge* etemp = poly->elist[i];
					Vertex* v1 = etemp->verts[0];
					Vertex* v2 = etemp->verts[1];

					// delete any old crossing points
					if (etemp->crossing) {
						delete etemp->crossing;
						etemp->crossing = NULL;
					}

					if (v1->level != v2->level) {
						double scalar1 = v1->scalar;
						double scalar2 = v2->scalar;

						double Alpha = (contour - scalar1) / (scalar2 - scalar1);
						double cross_x = Alpha * (v2->x - v1->x) + v1->x;
						double cross_y = Alpha * (v2->y - v1->y) + v1->y;
						double cross_z = contour;

						icVector3* cross_point = new icVector3(cross_x, cross_y, cross_z);
						etemp->crossing = cross_point;
					}
				}
				// Getting edge crossing points

				PolyLine N_contour;

				for (int i = 0; i < poly->nquads; i++) {
					Quad* qtemp = poly->qlist[i];

					std::vector<icVector3*> crossing_points;
					for (int j = 0; j < 4; j++) {
						Edge* etemp = qtemp->edges[j];
						if (etemp->crossing) {
							crossing_points.push_back(etemp->crossing);
						}
					}

					switch (crossing_points.size()) {
					case 2: {
						LineSegment line(*crossing_points[0], *crossing_points[1]);
						N_contour.push_back(line);  // find out how to get things to show up correctly 
					}
					case 4: {
						break;
					}
					}
				}
				lines.push_back(N_contour);

				contour += N_scalar;
			}
			glutPostRedisplay();
			break;
		}
	case '8':
		display_mode = 8;
		{
			//Finding min and max scalars 
			double N_scalar = (poly->max_scalar - poly->min_scalar) / N;
			double contour = poly->min_scalar;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < poly->nverts; j++) {
					Vertex* vtemp = poly->vlist[j];
					if (vtemp->scalar >= contour) {
						vtemp->level = 1;
					}
					else {
						vtemp->level = 0;
					}
				}
				for (int i = 0; i < poly->nedges; i++) {
					Edge* etemp = poly->elist[i];
					Vertex* v1 = etemp->verts[0];
					Vertex* v2 = etemp->verts[1];

					// delete any old crossing points
					if (etemp->crossing) {
						delete etemp->crossing;
						etemp->crossing = NULL;
					}

					if (v1->level != v2->level) {
						double scalar1 = v1->scalar;
						double scalar2 = v2->scalar;

						double Alpha = (contour - scalar1) / (scalar2 - scalar1);
						double cross_x = Alpha * (v2->x - v1->x) + v1->x;
						double cross_y = Alpha * (v2->y - v1->y) + v1->y;
						double cross_z = contour;

						icVector3* cross_point = new icVector3(cross_x, cross_y, cross_z);
						etemp->crossing = cross_point;
					}
				}
				// Getting edge crossing points

				PolyLine N_contour;

				for (int i = 0; i < poly->nquads; i++) {
					Quad* qtemp = poly->qlist[i];

					std::vector<icVector3*> crossing_points;
					for (int j = 0; j < 4; j++) {
						Edge* etemp = qtemp->edges[j];
						if (etemp->crossing) {
							crossing_points.push_back(etemp->crossing);
						}
					}

					switch (crossing_points.size()) {
					case 2: {
						LineSegment line(*crossing_points[0], *crossing_points[1]);
						N_contour.push_back(line);  // find out how to get things to show up correctly 
					}
					case 4: {
						break;
					}
					}
				}
				lines.push_back(N_contour);

				contour += N_scalar;

				//Critical points

				for (int i = 0; i < poly->nquads; i++) {
					Quad* qtemp = poly->qlist[i];
					Vertex* v1 = qtemp->verts[0];
					Vertex* v2 = qtemp->verts[1];
					Vertex* v3 = qtemp->verts[2];
					Vertex* v4 = qtemp->verts[3];

					double x1 = v3->x;
					double x2 = v1->x;
					double y1 = v3->y;
					double y2 = v1->y;

					double fx1y1 = v3->scalar;
					double fx2y1 = v4->scalar;
					double fx1y2 = v2->scalar;
					double fx2y2 = v1->scalar;

					
					double x = ((x2 * fx1y1) - (x1 * fx2y1) - (x2 * fx1y2) + (x1 * fx2y2)) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);
					double y = ((y2 * fx1y1) - (y2 * fx2y1) - (y1 * fx1y2) + (y1 * fx2y2)) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);

					if (x >= x1 && x <= x2 && y >= y1 && y <= y2) {
							
							double fxy = ( (((x2 - x) / (x2 - x1)) * ((y2 - y) / (y2 - y1))) * fx1y1) +
										 ( (((x - x1) / (x2 - x1)) * ((y2 - y) / (y2 - y1))) * fx2y1) +
								         ( (((x2 - x) / (x2 - x1)) * ((y - y1) / (y2 - y1))) * fx1y2) +
									      (((x - x1) / (x2 - x1)) * ((y - y1) / (y2 - y1))) * fx2y2 ;

							

				
							double z0 = 0;

							points.push_back(icVector3(x, y, z0));
					}
					
				}
				
			}


			glutPostRedisplay();
			break;
		}
	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;
		}
	}

	/******************************************************************************
	Callback function for dragging mouse
	******************************************************************************/

	void motion(int x, int y) {
		float r[4];
		float s, t;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		switch (mouse_mode) {
		case 2:

			Quaternion rvec;

			mat_to_quat(rotmat, rvec);
			trackball(r, s_old, t_old, s, t);
			add_quats(r, rvec, rvec);
			quat_to_mat(rvec, rotmat);

			s_old = s;
			t_old = t;

			display();
			break;

		case 1:

			translation[0] += (s - s_old);
			translation[1] += (t - t_old);

			s_old = s;
			t_old = t;

			display();
			break;
		}
	}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_SHIFT) {  // build up the selection feedback mode

				/*select face*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_quads(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_quad = processHits(hits, selectBuf);
				printf("Selected quad id = %d\n", poly->selected_quad);
				glutPostRedisplay();

				CHECK_GL_ERROR();

			}
			else if (button == GLUT_LEFT_BUTTON && key == GLUT_ACTIVE_CTRL)
			{
				/*select vertex*/

				GLuint selectBuf[512];
				GLint hits;
				GLint viewport[4];

				glGetIntegerv(GL_VIEWPORT, viewport);

				glSelectBuffer(win_width, selectBuf);
				(void)glRenderMode(GL_SELECT);

				glInitNames();
				glPushName(0);

				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();

				/*  create 5x5 pixel picking region near cursor location */
				gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y), 1.0, 1.0, viewport);

				set_view(GL_SELECT);
				set_scene(GL_SELECT, poly);
				display_vertices(GL_SELECT, poly);

				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				glFlush();

				glMatrixMode(GL_MODELVIEW);

				hits = glRenderMode(GL_RENDER);
				poly->selected_vertex = processHits(hits, selectBuf);
				printf("Selected vert id = %d\n", poly->selected_vertex);
				glutPostRedisplay();

			}

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);

	// reset IBFV pixels buffer
	free(pixels);
	initIBFV();
}

/******************************************************************************
Display IBFV vector field visualization (used for Project 3)
******************************************************************************/

void displayIBFV()
{
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_BLEND);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.5, 0.5, 0.5, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw the mesh using pixels and use vector field to advect texture coordinates
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	double modelview_matrix[16], projection_matrix[16];
	int viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix);
	glGetIntegerv(GL_VIEWPORT, viewport);

	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];

		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];

			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);

			tx = tx / win_width;
			ty = ty / win_height;
			
			icVector2 dp = icVector2(vtemp->vx, vtemp->vy);
			normalize(dp);
			dp *= dmax;

			double dx = -dp.x;
			double dy = -dp.y;

			float px = tx + dx;
			float py = ty + dy;

			glTexCoord2f(px, py);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}

	glEnable(GL_BLEND);

	// blend in noise pattern
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0);
	glScalef(2.0, 2.0, 1.0);

	glCallList(1);

	glBegin(GL_QUAD_STRIP);

	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	glTexCoord2f(0.0, tmax); glVertex2f(0.0, 1.0);
	glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
	glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	glDisable(GL_BLEND);

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// draw the mesh using pixels without advecting texture coords
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* qtemp = poly->qlist[i];
		glBegin(GL_QUADS);
		for (int j = 0; j < 4; j++)
		{
			Vertex* vtemp = qtemp->verts[j];
			double tx, ty, dummy;
			gluProject((GLdouble)vtemp->x, (GLdouble)vtemp->y, (GLdouble)vtemp->z,
				modelview_matrix, projection_matrix, viewport, &tx, &ty, &dummy);
			tx = tx / win_width;
			ty = ty / win_height;
			glTexCoord2f(tx, ty);
			glVertex3d(vtemp->x, vtemp->y, vtemp->z);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
}

/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, poly);

	/*display the mesh*/
	display_polyhedron(poly);

	/*display selected elements*/
	display_selected_vertex(poly);
	display_selected_quad(poly);


	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_polyhedron(Polyhedron* poly)
{
	unsigned int i, j;

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	case 3:	// checkerboard pattern display
	{
		glDisable(GL_LIGHTING);
		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 4:// points and lines drawing example
	case 6:
	
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		// draw lines
		for (int k = 0; k < lines.size(); k++)
		{
			drawPolyLine(lines[k], 1.0, 1.0, 0.0, 0.0);
		}

		// draw points
		for (int k = 0; k < points.size(); k++)
		{
			icVector3 point = points[k];
			drawDot(point.x, point.y, point.z);
		}
		break;
	}
	break;

	case 5:	// IBFV vector field display
	{
		displayIBFV();
		glutPostRedisplay();
	}
	break;
	case 7: //color and height 
	case 8:
	{
		for (int k = 0; k < lines.size(); k++)
		{
			double s = lines[k][0].start.z;
			double val_1 = (s - poly->min_scalar) / (poly->max_scalar - poly->min_scalar);
			double val_2 = (poly->max_scalar - s) / (poly->max_scalar - poly->min_scalar);

			double R = ((double)1.0 * val_1);
			double G = ((double)1.0 * val_2);
			double B = 1.0;
			drawPolyLine(lines[k], 1.0, R, G, B);
		}
		for (int k = 0; k < points.size(); k++) {
			icVector3 point = points.at(k);
			drawDot(point.entry[0], point.entry[1], point.entry[2]);

		}
	}
	break;

	default:
	{
		// don't draw anything
	}

	}
}
