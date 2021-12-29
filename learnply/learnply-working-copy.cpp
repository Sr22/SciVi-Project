#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <tuple>
#include <iostream>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "trackball.h"
#include "tmatrix.h"

#define PI 3.14159265

#include "drawUtil.h"

//TODO:
// Get quads
// Get critical points
// Classify critical points
// Make lines
Polyhedron* poly;
std::vector<PolyLine> morse_smale_lines;
std::vector<PolyLine> lines;
std::vector<icVector3> points;


//CUSTOM VARS HERE:
double const TINY_RADIUS = 0.005;
int const CIRCLE_POINTS = 5;//TODO: make bigger
double const DEGREES_TO_ADD = (360 / CIRCLE_POINTS);
double const THRESHOLD = 0.25;
double const STEP_SIZE = 0.05;
std::vector<icVector3> mins;
std::vector<icVector3> max;
std::vector<icVector3> saddles;

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

double min_scalar = NULL;
double max_scalar = NULL;
double max_bounds = NULL;
double c1[3] = { 1, 0, 0 };
double c2[3] = { 0, 1, .45 };
double HEIGHT_AMPLIFY = 0.5;

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

// Added functions
icVector2 normalized_grad_at_point(Quad* quad, double x, double y);
float bilinear_interpolate_at_point(Quad* quad, float x, float y);
Quad* get_quad_by_position(float x, float y);
Vertex* find_with_coords(double x, double y, int vert_count, Vertex** vlist);
void find_critical_points();
double det_hessian_matrix(Quad* quad);
void init_verts();
void morse_smale();
void update_scalars();
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

	//Mess it upupdate_scalars();

	//Give verts height and color.
	init_verts();

	//Call find critical points
	find_critical_points();

	//Morse Smale
	morse_smale();

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

void update_scalars() {
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vert = poly->vlist[i];
		double x = vert->x;
		double y = vert->y;
		vert->scalar = (x + .1) * (y + .1) - pow((x + .15), 4) - pow((x + .2), 6) * pow((y + 0.08), 3);
	}
}

/******************************************************************************
Find Verts
******************************************************************************/
void init_verts() {

	icVector3 center = poly->center;
	max_bounds = abs(center.array[0] - poly->vlist[0]->x);

	//First compute the min_scalar/max_scalar vals
	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vertex = poly->vlist[i];
		if (vertex->scalar < min_scalar || min_scalar == NULL) {
			min_scalar = vertex->scalar;
		}
		if (vertex->scalar > max_scalar || max_scalar == NULL) {
			max_scalar = vertex->scalar;
		}
		if (abs(center.array[0] - vertex->x) > max_bounds) {
			max_bounds = abs(center.array[0] - vertex->x);
		}
		if (abs(center.array[1] - vertex->y) > max_bounds) {
			max_bounds = abs(center.array[1] - vertex->y);
		}
	}

	for (int i = 0; i < poly->nverts; i++) {
		Vertex* vertex = poly->vlist[i];
		double gray_val = (vertex->scalar - min_scalar) / (max_scalar - min_scalar);
		double other_val = (max_scalar - vertex->scalar) / (max_scalar - min_scalar);
		vertex->R = gray_val * c1[0] + other_val * c2[0];
		vertex->G = gray_val * c1[1] + other_val * c2[1];
		vertex->B = gray_val * c1[2] + other_val * c2[2];
		double z_val = (vertex->scalar - min_scalar) / (max_scalar - min_scalar);
		vertex->scalar = z_val;
		vertex->z = z_val * max_bounds * HEIGHT_AMPLIFY;
	}
}


/******************************************************************************
Find Hessian Matrix //TODO: write definition
******************************************************************************/
double det_hessian_matrix(Quad* quad)
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

		double denominator = ((x2 - x1) * (y2 - y1));
		double fxy = (fx1y1/denominator) - (fx2y1/denominator) + (fx1y2 / denominator) - (fx2y2/denominator);
		double det = -(fxy) * (fxy);
		return det;
}
/******************************************************************************
Find Critical Points 
******************************************************************************/
void find_critical_points() 
{

	//Loop through each quad. For each quad,
	//a) Find when partial f/partial x = partial f / partial y = 0 (Peter / we did). 
	//If this point (x0, y0) is not in the quad we don't care
	//if it is, classify it and then add it to our list with its classification.
	//Maxima
	//Minima
	//Saddle
	for (int i = 0; i < poly->nquads; i++)
	{
		Quad* quad = poly->qlist[i];
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

		double x0 = (x2 * fx1y1 - x1 * fx2y1 - x2 * fx1y2 + x1 * fx2y2) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);
		double y0 = (y2 * fx1y1 - y2 * fx2y1 - y1 * fx1y2 + y1 * fx2y2) / (fx1y1 - fx2y1 - fx1y2 + fx2y2);
		
		
		double fx0y0 = bilinear_interpolate_at_point(quad, x0, y0) ;
		
		if (x1 <= x0 && x0 <= x2 && y1 <= y0 && y0 <= y2) {
			double det = det_hessian_matrix(quad);
			if (det < 0) {
				 saddles.push_back(icVector3(x0, y0, fx0y0* max_bounds * HEIGHT_AMPLIFY ) ); 
			} else {
				double x = x0 + TINY_RADIUS;
				double y = y0;
				double theta = 0;
				bool any_greater = false;
				bool any_less = false;
				for (int i = 0; i < CIRCLE_POINTS; i ++) {	
					double fxy = bilinear_interpolate_at_point(quad, x, y);
					if (fxy > fx0y0) {
						any_greater = true;
					} else if (fxy < fx0y0) {
						any_less = true;
					}
					theta = theta + DEGREES_TO_ADD;
					x = TINY_RADIUS * cos(theta * PI / 180.0) + (x0);
					y = TINY_RADIUS * sin(theta * PI / 180.0) + (y0);
				}

				if (any_greater && !any_less) {
					//minimum
					mins.push_back(icVector3(x0, y0, fx0y0 * max_bounds * HEIGHT_AMPLIFY));

				} else if (any_less && !any_greater) {
					//maximum
					max.push_back(icVector3(x0, y0, fx0y0* max_bounds * HEIGHT_AMPLIFY ));
				} else {
					//Saddle
					saddles.push_back(icVector3(x0, y0, fx0y0* max_bounds * HEIGHT_AMPLIFY ));
				}
			}
		}
	}
	for (int i = 0; i < poly->nverts; i++)
	{
		//Loop through all verts
		Vertex* vert = poly->vlist[i];
		std::vector<Vertex*> neighbors;
		int neighbor_count = 0;
		//For each vert, loop through its edges to find all of its neighbors. Store those in an array of size 4
		for (int j = 0; j < vert->nedges; j++) {
			Edge* edge = vert->edges[j];
			//if this edge contains vert, add the other vertex to neighbors
			if (edge->verts[0] == vert) {
				neighbors.push_back(edge->verts[1]);
				neighbor_count ++;
			}
			else if (edge->verts[1] == vert)  {
				neighbors.push_back(edge->verts[0]);
				neighbor_count++;
			}
		}
		if (neighbor_count == 4) {
			double fxy = vert->scalar;
			Vertex* neighbor_1 = neighbors.at(0);
			double f_n1 = neighbor_1->scalar;
			Vertex* neighbor_2 = neighbors.at(1);
			double f_n2 = neighbor_2->scalar;
			Vertex* neighbor_3 = neighbors.at(2);
			double f_n3 = neighbor_3->scalar;
			Vertex* neighbor_4 = neighbors.at(3);
			double f_n4 = neighbor_4->scalar;
			if (f_n1 < fxy && f_n2 < fxy && f_n3 < fxy && f_n4 < fxy) {
				// then add (x, y, fxy) to max list
				max.push_back(icVector3(vert->x, vert->y, fxy * max_bounds * HEIGHT_AMPLIFY));
			} else if (f_n1 > fxy && f_n2 > fxy && f_n3 > fxy && f_n4 > fxy) {
				//then add (x, y, fxy) to min list
				mins.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
			} else if (
				f_n1 > fxy && f_n2 > fxy && f_n3 < fxy && f_n4 < fxy
			) {
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_2 || edge->verts[1] == neighbor_2) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
				}
			} else if (
				f_n1 > fxy && f_n2 < fxy && f_n3 > fxy && f_n4 < fxy
			) {
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_3 || edge->verts[1] == neighbor_3) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy * max_bounds * HEIGHT_AMPLIFY));
				}
				//If n1 and n3 not adjacent, add to saddles list
			} else if (
				f_n1 > fxy && f_n2 < fxy && f_n3 < fxy && f_n4 > fxy
			) {
				//If n1 and n4 not adjacent, add to saddles list
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_4 || edge->verts[1] == neighbor_4) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
				}
			}  else if (
				f_n1 < fxy && f_n2 < fxy && f_n3 > fxy && f_n4 > fxy
			) {
				//If n1 and n2 not adj, saddle
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_2 || edge->verts[1] == neighbor_2) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
				}
			} else if (
				f_n1 < fxy && f_n2 > fxy && f_n3 < fxy && f_n4 > fxy
			) {
				//If n1 and n3 not adj, saddle
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_3 || edge->verts[1] == neighbor_3) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
				}
			} else if (
				f_n1 < fxy && f_n2 > fxy && f_n3 > fxy && f_n4 < fxy
			) {
				//(n1, n4) not adj=>saddle
				bool adj = false;
				for (int k = 0; k < neighbor_1->nedges; k++) {
					Edge* edge = neighbor_1->edges[k];
					if (edge->verts[0] == neighbor_4 || edge->verts[1] == neighbor_4) {
						adj = true;
					}
				}
				if (!adj) {
					//Add to saddles list.
					saddles.push_back(icVector3(vert->x, vert->y, fxy* max_bounds* HEIGHT_AMPLIFY));
				}
			} 
		}
		//If there are 4 neighbors,
			//Using that array, check if:
			//A) all are less. if so, add to max crit point list
			//B) all are greater. If so, add to min crit point list
			//C) 2 are greater and 2 are less, and the greater ones are NOT adjacent (Have to loop through edges again probably is easiest)
		//Else, (TODO if time)
	}
}
/******************************************************************************
Morse-Smale
******************************************************************************/
void morse_smale()
{
 //Loop through all min vertices
	//For each CIRCLE_POINT, draw an integral line from that circle point using a polyline
	//by taking STEP_SIZE size steps in the direction of gradient until no longer increasing.
	//Then check if we are anywhere near any of our saddle points. If we are, add the polyline to our morsesmale lines collection
	
	float x;
	float y;

	for (int i = 0; i < mins.size(); i++) {
		std::vector<PolyLine> polylines;
		Quad* quad_xy;
		Quad* quad_next_xy;
		double x0 = mins[i].x;
		double y0 = mins[i].y;
		
		x = x0 + TINY_RADIUS;
		y = y0;
		double theta = 0;
		for (int j = 0; j < CIRCLE_POINTS; j ++) {	
		    quad_xy = get_quad_by_position(x, y);
			double fxy = bilinear_interpolate_at_point(quad_xy, x, y);
			icVector2 grad = normalized_grad_at_point(quad_xy, x, y);

			double next_x = x + (grad.x * STEP_SIZE);
			double next_y = y + (grad.y * STEP_SIZE);
			quad_next_xy = get_quad_by_position(next_x, next_y);
			double f_next_xy = bilinear_interpolate_at_point(quad_next_xy, next_x, next_y);
			PolyLine pl;
			while (f_next_xy > fxy) {
				//Add to the polyline
				pl.push_back(LineSegment(x, y, fxy, next_x, next_y, f_next_xy));
				//Move to next step
				x = next_x;
				y = next_y;
				fxy = f_next_xy;
				grad = normalized_grad_at_point(quad_xy, x, y);
				next_x = x + grad.x * STEP_SIZE;
				next_y = y + grad.y * STEP_SIZE;
				quad_next_xy = get_quad_by_position(next_x, next_y);
				f_next_xy = bilinear_interpolate_at_point(quad_next_xy, next_x, next_y);
			}

			double last_x = x + (grad.x * STEP_SIZE) / 2;
			double last_y = y + (grad.y * STEP_SIZE) / 2;
			Quad* last_quad = get_quad_by_position(last_x, last_y);
			double last_fxy = bilinear_interpolate_at_point(last_quad, last_x, last_y);
			pl.push_back(LineSegment(x, y, fxy, last_x, last_y, last_fxy));
			polylines.push_back(pl);

		    theta = theta + DEGREES_TO_ADD;
			double val = cos(theta * PI / 180);
			x = TINY_RADIUS * cos(theta* PI / 180.0) + (x0);
			y = TINY_RADIUS * sin(theta* PI / 180.0) + (y0);
		}

		//For each saddle check if any of our last step is within radius and closest to any saddle points. If it is, then add to morse smale lines
		for (int k = 0; k < saddles.size(); k++) {
			double closest_polyline_distance = THRESHOLD;
			icVector3 saddle_point = saddles[k];
			PolyLine line_to_add;
			bool add = false;
			for (int j = 0; j < CIRCLE_POINTS; j++) {
					PolyLine pl = polylines[j];
					//If the line has any point "close enough" to a saddle point, add it.
					for (int b = 0; b < pl.size(); b++) {
						LineSegment ls = pl[b];
						icVector3 end = ls.end;
						double dist = sqrt(pow((saddle_point.x - end.x), 2) + pow((saddle_point.y - end.y), 2));
						if (dist < closest_polyline_distance) {
							line_to_add = pl;
							add = true;
						}
					}
					
					
			}
			if (add == true) {
				morse_smale_lines.push_back(line_to_add);
			}
		}
	}

		for (int i = 0; i < max.size(); i++) {
				std::vector<PolyLine> polylines;
				Quad* quad_xy;
				Quad* quad_next_xy;
				double x0 = max[i].x;
				double y0 = max[i].y;
		
				x = x0 + TINY_RADIUS;
				y = y0;
				double theta = 0;
				for (int j = 0; j < CIRCLE_POINTS; j ++) {	
					quad_xy = get_quad_by_position(x, y);
					double fxy = bilinear_interpolate_at_point(quad_xy, x, y);
					icVector2 grad = normalized_grad_at_point(quad_xy, x, y);

					double next_x = x - (grad.x * STEP_SIZE);
					double next_y = y - (grad.y * STEP_SIZE);
					quad_next_xy = get_quad_by_position(next_x, next_y);
					double f_next_xy = bilinear_interpolate_at_point(quad_next_xy, next_x, next_y);
					PolyLine pl;
					while (f_next_xy < fxy) {
						//Add to the polyline
						pl.push_back(LineSegment(x, y, fxy, next_x, next_y, f_next_xy));
						//Move to next step
						x = next_x;
						y = next_y;
						fxy = f_next_xy;
						grad = normalized_grad_at_point(quad_xy, x, y);
						next_x = x + grad.x * STEP_SIZE;
						next_y = y + grad.y * STEP_SIZE;
						quad_next_xy = get_quad_by_position(next_x, next_y);
						f_next_xy = bilinear_interpolate_at_point(quad_next_xy, next_x, next_y);
					}

					double last_x = x + (grad.x * STEP_SIZE) / 2;
					double last_y = y + (grad.y * STEP_SIZE) / 2;
					Quad* last_quad = get_quad_by_position(last_x, last_y);
					double last_fxy = bilinear_interpolate_at_point(last_quad, last_x, last_y);
					pl.push_back(LineSegment(x, y, fxy, last_x, last_y, last_fxy));
					polylines.push_back(pl);

					theta = theta + DEGREES_TO_ADD;
					double val = cos(theta * PI / 180);
					x = TINY_RADIUS * cos(theta* PI / 180.0) + (x0);
					y = TINY_RADIUS * sin(theta* PI / 180.0) + (y0);
				}

				//For each saddle check if any of our last step is within radius and closest to any saddle points. If it is, then add to morse smale lines
				for (int k = 0; k < saddles.size(); k++) {
					double closest_polyline_distance = THRESHOLD;
					icVector3 saddle_point = saddles[k];
					PolyLine line_to_add;
					bool add = false;
					for (int j = 0; j < CIRCLE_POINTS; j++) {
							PolyLine pl = polylines[j];
							//If the line has any point "close enough" to a saddle point, add it.
							for (int b = 0; b < pl.size(); b++) {
								LineSegment ls = pl[b];
								icVector3 end = ls.end;
								double dist = sqrt(pow((saddle_point.x - end.x), 2) + pow((saddle_point.y - end.y), 2));
								if (dist < closest_polyline_distance) {
									line_to_add = pl;
									add = true;
								}
							}
					
					
					}
					if (add == true) {
						morse_smale_lines.push_back(line_to_add);
					}
				}
			}	

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
icVector2 normalized_grad_at_point(Quad* quad, double x0, double y0)
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
				double partial_p1_partial_x = (-1 / (x2 - x1)) * ((y2 - y0) / (y2 - y1));
				double partial_p2_partial_x = (1 / (x2 - x1)) * ((y2 - y0) / (y2 - y1));
				double partial_p3_partial_x = (-1 / (x2 - x1)) * ((y0 - y1) / (y2 - y1));
				double partial_p4_partial_x = (1 / (x2 - x1)) * ((y0 - y1) / (y2 - y1));
				double grad_x = (partial_p1_partial_x * fx1y1) + (partial_p2_partial_x * fx2y1) + (partial_p3_partial_x * fx1y2) + (partial_p4_partial_x * fx2y2);
				double partial_p1_partial_y = (-1 / (y2 - y1)) * ((x2 - x0) / (x2 - x1));
				double partial_p2_partial_y = (-1 / (y2 - y1)) * ((x0 - x1) * (x2 - x1));
				double partial_p3_partial_y = (1 / (y2 - y1)) * ((x2 - x0) / (x2 - x1));
				double partial_p4_partial_y = (1 / (y2 - y1)) * ((x0 - x1) / (x2 - x1));
				double grad_y = (partial_p1_partial_y * fx1y1) + (partial_p2_partial_y * fx2y1) + (partial_p3_partial_y * fx1y2) + (partial_p4_partial_y * fx2y2);
				icVector2 grad = icVector2(grad_x, grad_y);
				normalize(grad);
				return grad;
			}
		}
		std::cout<<"OUT OF BOUNDS"<<std::endl;
		return icVector2(0.0, 0.0);
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
				return (p1 * fx1y1 + p2 * fx2y1 + p3 * fx1y2 + p4 * fx2y2);
			}
		} else {
			return NULL;
		}
}

/******************************************************************************
Get Quad by Position 
******************************d************************************************/
Quad* get_quad_by_position(float x, float y)
{
	//TODO: Optimize if time.
	//loop through all quads
	//Return quad that works
	for (int i = 0; i < poly->nquads; i++) {
		Quad* quad = poly->qlist[i];
		Vertex* v1 = (*quad).verts[0];
		Vertex* v2 = quad->verts[1];
		Vertex* v3 = quad->verts[2];
		Vertex* v4 = quad->verts[3];

		double x1 = std::min(std::min(v3->x, v1->x), std::min(v2->x, v4->x));
		double x2 = std::max(std::max(v3->x, v1->x), std::max(v2->x, v4->x));
		double y1 = std::min(std::min(v3->y, v1->y), std::min(v2->y, v4->y));
		double y2 = std::max(std::max(v3->y, v1->y), std::max(v2->y, v4->y));

		if (x1 <= x && x <= x2 && y1 <= y && y <= y2) {
			return quad;
		}
	}
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
	//lines.clear();
	//points.clear();

	switch (key) { 
	case 27:	// set excape key to exit program
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '1':	//Critical points with their classifying colors
		display_mode = 1;
		glutPostRedisplay();
		break;

	case '2':	// Critical points, classifying colors, Morse-Smale lines, height
		display_mode = 2;
		glutPostRedisplay();
		break;

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

	switch (display_mode) {
	case 1:	// just critical points and colors
	{
		glDisable(GL_LIGHTING);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glColor3f(temp_v->R, temp_v->G, temp_v->B);
				glVertex3d(temp_v->x, temp_v->y, 0);
			}
			glEnd();
		}
		// draw points
		for (int k = 0; k < saddles.size(); k++)
		{
			icVector3 point = saddles[k];
			drawDot(point.x, point.y, 0, 0.15, 0, 0, 1);
		}
		// draw points
		for (int k = 0; k < mins.size(); k++)
		{
			icVector3 point = mins[k];
			drawDot(point.x, point.y, 0, .15, 1, 0, 0);
		}
		// draw points
		for (int k = 0; k < max.size(); k++)
		{
			icVector3 point = max[k];
			drawDot(point.x, point.y, 0, 0.15, 0, 1, 0);
		}
	}
	break;

	case 2:	//morse smale + points + height
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_SMOOTH);
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
		// draw points
		for (int k = 0; k < saddles.size(); k++)
		{
			icVector3 point = saddles[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 0, 1);
		}
		// draw points
		for (int k = 0; k < mins.size(); k++)
		{
			icVector3 point = mins[k];
			drawDot(point.x, point.y, point.z, .15, 1, 0, 0);
		}
		// draw points
		for (int k = 0; k < max.size(); k++)
		{
			icVector3 point = max[k];
			drawDot(point.x, point.y, point.z, 0.15, 0, 1, 0);
		}
		// draw lines
		for (int k = 0; k < morse_smale_lines.size(); k++)
		{
			drawPolyLine(morse_smale_lines[k], 1.0, 1.0, 0.0, 0.0);
		}
	}
	break;
	default:
	{
		// don't draw anything
	}
	}
}
