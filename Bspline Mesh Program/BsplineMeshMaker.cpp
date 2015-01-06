// Lab5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

/* E. Angel, Interactive Computer Graphics */
/* A Top-Down Approach with OpenGL, Third Edition */
/* Addison-Wesley Longman, 2003 */

/*
**
**  Simple interactive curve drawing program.
**
**  The following keyboard commands are used to control the
**  program:
**
**    q - Quit the program
**    c - Clear the screen
**    e - Erase the curves
**    b - Draw Bezier curves
**
*/

/* All control points are converted to Bezier
control point to allow use of OpenGL evaluators */

#include <GL/glut.h>
#include <math.h>

typedef enum
{
	BEZIER,
} curveType;

void keyboard(unsigned char key, int x, int y);

/* Colors to draw them in */
GLfloat colors[][3] =
{
	{ 1.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 1.0 }
};


#define MAX_CPTS  25            /* Fixed maximum number of control points */
#define MAX_KNOTS 30
#define MAX_VERTICES 60000		
#define MAX_TRIANGLES 60000
#define MAX_SIZE 256
GLfloat cpts[MAX_CPTS][3];
GLfloat vertex[MAX_VERTICES][3];
GLfloat vertex_normal[MAX_VERTICES][3];
int triangle[MAX_TRIANGLES][3];
GLfloat tri_normal[MAX_TRIANGLES][3];
GLfloat texture[3 * MAX_SIZE*MAX_SIZE];	/* Texture pattern. */
int ncpts = 0;

static int width = 500, height = 500;           /* Window width and height */

/*Flags*/

int ctrlpt_flag = 1;
int ctrlpolygon_ON = 1;
int Bspline_ON = 1;
int move_flag = 0; /*If move flag is set to 1, then the selected point can move*/
int mouse_mode_flag = 0; /* 0 is insertion, 1 is selection*/
int selectp = -1; /*Will be used to select the point*/
int surface_flag = 1; /*Default is set to off, if on, it will sweep the curve to create a surface*/
int selectp_flag = 0; /*Default, there is no selected point (0), if set to 1, it will be used in motion*/

float knot[MAX_KNOTS]; /*Bspline knot array*/

/*B points*/
int numBpts = -1;
float Bpts[MAX_KNOTS * 50][3];

FILE *savefile = NULL;/*file pointer used for saving*/
FILE *outfile = NULL; /*file pointer used for saving mesh*/
char *fileName = "Lab5.txt"; /*filename for the savefile*/

/*Constants*/
float zmin = -5.0;
float zmax = 5.0;
float zd = 0.25;
int zflag = -1; /*If zero then it changes zmin, if 1 then it changes zmax */
float increment = 0.5;

/*Surface variables*/
int num_triangles; int num_vertices;
int type_of_surface = -1; /*If set to 0, it is wireframe, if it is set to 1, it is shaded surface, if it is set to 2, it is a textured surface*/
float radian, convert, costheta, sintheta;
int current, offset, last_tri_strip, prev_vertex, theta;
float surf_x_min, surf_x_max, surf_y_min, surf_y_max, surf_z_min, surf_z_max, surf_x_center, surf_y_center, surf_z_center;

/*Texture_Variables*/
GLubyte image[64][64][3];
GLfloat texCoord[MAX_KNOTS * 50][2];


/* Draw the indicated curves using the current control points. */
static void drawCurves(curveType type)
{
	int i;
	int step;
	GLfloat newcpts[4][3];

	float m[4][4];

	/* Set the step size. */
	step = 3;

	glColor3fv(colors[type]);

	/* Draw the curves */
	i = 0;
	while (i + 3 < ncpts)
	{
		glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 4, &cpts[i][0]);
		glMapGrid1f(30, 0.0, 1.0);
		glEvalMesh1(GL_LINE, 0, 30);

		/* Advance to the next segment */
		i += step;
	}
	glFlush();
}

float CoxdeBoor(int i, int p, float t)
{
	float left, right;
	if (p == 1) {
		if (knot[i] < knot[i + 1] && knot[i] <= t && t < knot[i + 1])
			return(1.0);
		else
			return(0.0);
	}
	else {
		if (knot[i + p - 1] - knot[i] != 0.0)
			left = CoxdeBoor(i, p - 1, t)*(t - knot[i]) /
			(knot[i + p - 1] - knot[i]);
		else
			left = 0.0;
		if (knot[i + p] - knot[i + 1] != 0.0)
			right = CoxdeBoor(i + 1, p - 1, t)*(knot[i + p] - t) /
			(knot[i + p] - knot[i + 1]);
		else
			right = 0.0;
		return(left + right);
	}
}


static void computeBspline()
{
	int i, j;
	float t;
	float B0, B1, B2, B3;
	int m = ncpts - 1;
	if (m < 3)
	{
		return;
	}
	for (i = 0; i <= 3; i++)
		knot[i] = 0.0;
	for (i = 4; i <= m; i++)
		knot[i] = i - 3.0;
	for (i = m + 1; i <= m + 4; i++)
		knot[i] = m - 2.0;
	int num_knots = m + 4.0;

	numBpts = 0;
	Bpts[0][0] = cpts[0][0];
	Bpts[0][1] = cpts[0][1];

	for (i = 3; i < num_knots - 3; i++) {
		for (t = knot[i]; t < knot[i + 1]; t += 0.2) {
			B0 = CoxdeBoor(i, 4, t);
			B1 = CoxdeBoor(i - 1, 4, t);
			B2 = CoxdeBoor(i - 2, 4, t);
			B3 = CoxdeBoor(i - 3, 4, t);
			numBpts++;
			Bpts[numBpts][0] =
				cpts[i][0] * B0 +
				cpts[i - 1][0] * B1 +
				cpts[i - 2][0] * B2 +
				cpts[i - 3][0] * B3;
			Bpts[numBpts][1] =
				cpts[i][1] * B0 +
				cpts[i - 1][1] * B1 +
				cpts[i - 2][1] * B2 +
				cpts[i - 3][1] * B3;
		}
	}
	numBpts++;
	Bpts[numBpts][0] = cpts[ncpts - 1][0];
	Bpts[numBpts][1] = cpts[ncpts - 1][1];

}
static void drawPolygon()
{
	int i = 0;
	if (ctrlpolygon_ON != 1)
	{
		return;
	}
	else
	{
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < ncpts; i++)
		{
			glVertex2f(cpts[i][0], cpts[i][1]);
		}
		glEnd();
	}
}

void Accumulate(int triangle_num, int vertex1, int vertex2, int vertex3)
{
	vertex_normal[vertex1][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex1][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex1][2] += tri_normal[triangle_num][2];
	vertex_normal[vertex2][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex2][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex2][2] += tri_normal[triangle_num][2];
	vertex_normal[vertex3][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex3][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex3][2] += tri_normal[triangle_num][2];
}


/* Compute the surface normal for each triangle. */
void Compute_Triangle_Normals()
{
	int i, counter_clock_wise = 0; /* for one sided. */
	int vertex1, vertex2, vertex3;
	float ax, ay, az, bx, by, bz, trianglenorm;
	for (i = 0; i <= num_triangles; i++) {
		vertex1 = triangle[i][0];
		vertex2 = triangle[i][1];
		vertex3 = triangle[i][2];
		if (counter_clock_wise) {
			ax = vertex[vertex2][0] - vertex[vertex1][0];
			ay = vertex[vertex2][1] - vertex[vertex1][1];
			az = vertex[vertex2][2] - vertex[vertex1][2];

			bx = vertex[vertex3][0] - vertex[vertex1][0];
			by = vertex[vertex3][1] - vertex[vertex1][1];
			bz = vertex[vertex3][2] - vertex[vertex1][2];
		}
		else {
			ax = vertex[vertex3][0] - vertex[vertex1][0];
			ay = vertex[vertex3][1] - vertex[vertex1][1];
			az = vertex[vertex3][2] - vertex[vertex1][2];

			bx = vertex[vertex2][0] - vertex[vertex1][0];
			by = vertex[vertex2][1] - vertex[vertex1][1];
			bz = vertex[vertex2][2] - vertex[vertex1][2];
		}
		/* Compute the cross product. */
		tri_normal[i][0] = ay*bz - by*az;
		tri_normal[i][1] = bx*az - ax*bz;
		tri_normal[i][2] = ax*by - bx*ay;
		trianglenorm = sqrt(tri_normal[i][0] * tri_normal[i][0] +
			tri_normal[i][1] * tri_normal[i][1] +
			tri_normal[i][2] * tri_normal[i][2]);
		if (trianglenorm != 0.0) {
			tri_normal[i][0] = tri_normal[i][0] / trianglenorm;
			tri_normal[i][1] = tri_normal[i][1] / trianglenorm;
			tri_normal[i][2] = tri_normal[i][2] / trianglenorm;
		}
	}
}
/* Compute the surface normal at each vertex. */
void Compute_Vertex_Normals()
{
	int i;
	float vertexnorm;
	for (i = 0; i <= num_vertices; i++) {
		vertex_normal[i][0] = 0.0;
		vertex_normal[i][1] = 0.0;
		vertex_normal[i][2] = 0.0;
	}
	/* Add the normal at each triangle to the triangle's vertices. */
	for (i = 0; i <= num_triangles; i++)
		Accumulate(i, triangle[i][0], triangle[i][1], triangle[i][2]);
	/* Normalize the normal at each vertex. */
	for (i = 0; i <= num_vertices; i++) {
		vertexnorm = sqrt(vertex_normal[i][0] * vertex_normal[i][0] +
			vertex_normal[i][1] * vertex_normal[i][1] +
			vertex_normal[i][2] * vertex_normal[i][2]);
		if (vertexnorm != 0.0) {
			vertex_normal[i][0] /= vertexnorm;
			vertex_normal[i][1] /= vertexnorm;
			vertex_normal[i][2] /= vertexnorm;
		}
	}
}

/*Function to Sweep Surface*/
void sweep_surface(void)
{
	if (numBpts <= -1)
	{
		printf("Not enough curve points");
		return;
	}
	radian = 0.0;
	convert = 3.1415629 / 180.0;
	costheta = cos(radian);
	sintheta = sin(radian);
	/* The B-Spline curve at theta=0 deg. */
	current = 0;
	num_vertices = -1;
	num_triangles = -1;
	int i;
	int theta_incr = 30;
	float x_curve_rotate[2][MAX_KNOTS * 50];
	float y_curve_rotate[2][MAX_KNOTS * 50];
	float z_curve_rotate[2][MAX_KNOTS * 50];
	for (i = 0; i <= numBpts; i++) {
		z_curve_rotate[current][i] = 0.0;
		if (num_vertices >= MAX_VERTICES) {
			printf("Error: Exceeded MAX_VERTICES\n");
			exit(0);
		}
		num_vertices++;
		vertex[num_vertices][0] = Bpts[i][0];
		vertex[num_vertices][1] = Bpts[i][1];
		vertex[num_vertices][2] = Bpts[i][2];
	}
	offset = num_vertices + 1;
	/* Rotate the curve. */
	for (theta = theta_incr; theta < 360.0; theta += theta_incr) {
		radian = theta*convert;
		costheta = cos(radian);
		sintheta = sin(radian);
		current = !current;
		for (i = 0; i <= numBpts; i++) {
			x_curve_rotate[current][i] =
				Bpts[i][0] * costheta;
			z_curve_rotate[current][i] =
				-Bpts[i][0] * sintheta;
			num_vertices++;
			if (num_vertices >= MAX_VERTICES) {
				printf("Error: Exceeded MAX_VERTICES\n");
				exit(0);
			}
			vertex[num_vertices][0] = x_curve_rotate[current][i];
			vertex[num_vertices][1] = Bpts[i][1];
			vertex[num_vertices][2] = z_curve_rotate[current][i];
			if (i != 0) {
				num_triangles++;
				if (num_triangles >= MAX_TRIANGLES) {
					printf("Error: Exceeded MAX_TRIANGLES\n"); exit(1);
				}
				triangle[num_triangles][0] = num_vertices;
				triangle[num_triangles][1] =
					num_vertices - offset - 1;
				triangle[num_triangles][2] = num_vertices - 1;
				num_triangles++;
				if (num_triangles >= MAX_TRIANGLES) {
					printf("Error: Exceeded MAX_TRIANGLES\n");
					exit(2);
				}
				triangle[num_triangles][0] = num_vertices;
				triangle[num_triangles][1] =
					num_vertices - offset;
				triangle[num_triangles][2] =
					num_vertices - offset - 1;
			}
		}
	}
	/* The last strip of the B-Spline surface */
	last_tri_strip = num_triangles + 1;
	prev_vertex = num_vertices - offset + 1;
	for (i = 1; i <= numBpts; i++) {
		prev_vertex++;
		num_triangles++;
		if (num_triangles >= MAX_TRIANGLES) {
			printf("Error: Exceeded MAX_TRIANGLES\n");
			exit(3);
		}
		triangle[num_triangles][0] = i;
		triangle[num_triangles][1] = prev_vertex - 1;
		triangle[num_triangles][2] = i - 1;
		num_triangles++;
		if (num_triangles >= MAX_TRIANGLES) {
			printf("Error: Exceeded MAX_TRIANGLES\n");
			exit(4);
		}
		triangle[num_triangles][0] = i;
		triangle[num_triangles][1] = prev_vertex;
		triangle[num_triangles][2] = prev_vertex - 1;
	}
	Compute_Triangle_Normals();
	Compute_Vertex_Normals();
}

void shade_surface()
{
	printf("shading surface\n");
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_diffuse[] = { 0.8, 0.0, 0.0, 1.0 };
	GLfloat mat_ambient[] = { 0.2, 0.0, 0.0, 1.0 };
	GLfloat mat_shininess = 100.0;

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialf(GL_FRONT, GL_SHININESS, mat_shininess);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
}

void read_text_pattern(int maxs, int maxt)
{
	FILE *infile;
	int s, t, index;
	float rval, gval, bval;
	int textdataR[MAX_SIZE*MAX_SIZE];
	int textdataG[MAX_SIZE*MAX_SIZE];
	int textdataB[MAX_SIZE*MAX_SIZE];

	infile = fopen("texture.bin", "rb");
	fread(textdataR, sizeof(int), maxs*maxt, infile);
	fread(textdataG, sizeof(int), maxs*maxt, infile);
	fread(textdataB, sizeof(int), maxs*maxt, infile);
	fclose(infile);

	for (t = 0; t < maxt; t++) {
		for (s = 0; s < maxs; s++) {
			rval = (float)textdataR[s + maxs * t] / 255.0;
			gval = (float)textdataG[s + maxs * t] / 255.0;
			bval = (float)textdataB[s + maxs * t] / 255.0;

			index = 3 * (s + maxs * t);
			texture[index] = rval;
			texture[index + 1] = gval;
			texture[index + 2] = bval;
		}
	}
	return;
}

void texture_surface()
{
	// create a texture 

	read_text_pattern(MAX_SIZE, MAX_SIZE);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, MAX_SIZE, MAX_SIZE, 0, GL_RGB, GL_FLOAT, texture);

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	// Blend the texture pattern with the shaded object.
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// Replace the shaded object's surface with the texture pattern
	//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	// Enable depth test and texture map
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
}

/* This routine displays the control points */
static void display(void)
{
	int i, v1, v2, v3;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor3f(0.0, 0.0, 0.0);
	glPointSize(5.0);
	if (ctrlpt_flag == 1)
	{
		glBegin(GL_POINTS);
		for (i = 0; i < ncpts; i++)
			glVertex3fv(cpts[i]);
		glEnd();
	}
	drawPolygon();
	if (Bspline_ON == 1)
	{
		if (ncpts > 3)
		{
			computeBspline();
			glBegin(GL_LINE_STRIP);
			for (i = 0; i <= numBpts; i++)
				glVertex2fv(Bpts[i]);
			glEnd();
		}
	}
	/*Redraw the point to highlight*/
	if (mouse_mode_flag == 1)
	{
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		glVertex3f(cpts[selectp][0], cpts[selectp][1], 0.0);
		glEnd();
	}
	if (surface_flag == 1)
	{
		if (type_of_surface == 0) /*Wireframe Surface*/
		{
			computeBspline();
			sweep_surface();
			for (i = 0; i <= num_triangles; i++)
			{
				v1 = triangle[i][0];
				v2 = triangle[i][1];
				v3 = triangle[i][2];
				glBegin(GL_LINE_LOOP);
				{
					glVertex3fv(vertex[v1]);
					glVertex3fv(vertex[v2]);
					glVertex3fv(vertex[v3]);
				}
				glEnd();
			}
		}
		else if (type_of_surface == 1) /*Shaded Surface*/
		{
			shade_surface();
			computeBspline();
			sweep_surface();
			glBegin(GL_TRIANGLES);
			for (i = 0; i <= num_triangles; i++)
			{
				glColor3f(1.0, 0.0, 0.0);
				glVertex3fv(vertex[triangle[i][0]]);
				glVertex3fv(vertex[triangle[i][1]]);
				glVertex3fv(vertex[triangle[i][2]]);
			}
			glEnd();
		}
		else if (type_of_surface == 2) /*Textured Surface*/
		{
			computeBspline();
			sweep_surface();

			glBegin(GL_TRIANGLES);
			for (i = 0; i <= num_triangles; i++)
			{
				v1 = triangle[i][0];
				v2 = triangle[i][1];
				v3 = triangle[i][2];

				glNormal3fv(vertex_normal[triangle[i][0]]);
				glTexCoord2fv(texCoord[v1]);
				glVertex3fv(vertex[triangle[i][0]]);

				glNormal3fv(vertex_normal[triangle[i][1]]);
				glTexCoord2fv(texCoord[v2]);
				glVertex3fv(vertex[triangle[i][1]]);

				glNormal3fv(vertex_normal[triangle[i][2]]);
				glTexCoord2fv(texCoord[v3]);
				glVertex3fv(vertex[triangle[i][2]]);
			}
			glEnd();
		}
	}
	glutSwapBuffers();
	glFlush();
}

float distanceform(float x1, float x2, float y1, float y2)
{
	float distancesquared = (((x1 - x2) * (x1 - x2)) + ((y1 - y2) - (y1 - y2)));
	return sqrtf(distancesquared);
}

/* This routine inputs new control points */
static void mouse(int button, int state, int x, int y)
{
	int i;
	float distance, distance_comp;
	float wx, wy;

	/* Translate back to our coordinate system */

	wx = (20.0 * x) / (float)(width - 1) - 10.0;
	wy = (20.0 * (height - 1 - y)) / (float)(height - 1) - 10.0;
	printf("%f, %f\n", wx, wy);
	if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN)
		return;
	if (ctrlpt_flag == 1)
	{
		if (mouse_mode_flag == 0) /*insertion mode*/
		{
			/* See if we have room for any more control points */
			if (ncpts == MAX_CPTS)
				return;

			/* Save the point */
			cpts[ncpts][0] = wx;
			cpts[ncpts][1] = wy;
			cpts[ncpts][2] = 0.0;
			ncpts++;

			/* Draw the point */
			glColor3f(0.0, 0.0, 0.0);
			glPointSize(5.0);
			glBegin(GL_POINTS);
			glVertex3f(wx, wy, 0.0);
			glEnd();
		}
	}
	if (mouse_mode_flag == 1) /*selection mode*/
	{
		/*find point that is the shortest distance from the mouse*/
		distance = distanceform(cpts[0][0], wx, cpts[0][1], wy);
		selectp = 0;
		for (i = 1; i < ncpts; i++)
		{
			distance_comp = distanceform(cpts[i][0], wx, cpts[i][1], wy);
			printf("distance:%f, distance_comp:%f\n", distance, distance_comp);
			if (distance > distance_comp)
			{
				distance = distance_comp;
				selectp = i;
				printf("selectp: %d\n", selectp);
			}
		}

		/*Redraw the point to highlight*/
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		glVertex3f(cpts[selectp][0], cpts[selectp][1], 0.0);
		glEnd();

		if (move_flag == 1)
		{
			selectp_flag = 1;
		}

	}
	glFlush();
}

void deletepoint()
{
	int i;
	for (i = selectp; i < ncpts - 1; i++)
	{
		cpts[i][0] = cpts[i + 1][0];
		cpts[i][1] = cpts[i + 1][1];
		cpts[i][2] = 0.0;
	}
	--ncpts;
	display();
}

void MouseMotion(int x, int y)
{
	float wx, wy;
	wx = (20.0 * x) / (float)(width - 1) - 10.0;
	wy = (20.0 * (height - 1 - y)) / (float)(height - 1) - 10.0;
	if (selectp_flag)
	{
		cpts[selectp][0] = wx;
		cpts[selectp][1] = wy;
		glutPostRedisplay();
	}
	else
		return;
}

void save()
{
	savefile = fopen(fileName, "w");
	if (savefile == NULL)
	{
		printf("Warning: Could not open %s\n", fileName);
	}
	else
	{
		for (int i = 0; i < ncpts; i++)
		{
			fprintf(savefile, " %f %f %f ", cpts[i][0], cpts[i][1], cpts[i][2]);
		}
	}
}

float findMin(int position)
{
	int i;
	float min = vertex[0][position];
	for (i = 1; i < 60000; i++)
	{
		if (min > vertex[i][position])
			min = vertex[i][position];
	}
	return min;
}

float findMax(int position)
{
	int i;
	float max = vertex[0][position];
	for (i = 1; i < 60000; i++)
	{
		if (max < vertex[i][position])
			max = vertex[i][position];
	}
	return max;
}
void savemesh()
{
	surf_x_min = findMin(0);
	surf_x_max = findMax(0);
	surf_y_min = findMin(1);
	surf_y_max = findMax(1);
	surf_z_min = findMin(2);
	surf_z_max = findMax(2);
	surf_x_center = (surf_x_min + surf_x_max) / 2;
	surf_y_center = (surf_y_min + surf_y_max) / 2;
	surf_z_center = (surf_z_min + surf_z_max) / 2;

	if ((outfile = fopen("bsplinemesh.jv", "wb")) == NULL) {
		printf("Could not open bsplinemesh.dk for write\n");
		return;
	}

	// Save the center and the bounding box of the mesh.

	fwrite(&surf_x_center, sizeof(float), 1, outfile);
	fwrite(&surf_y_center, sizeof(float), 1, outfile);
	fwrite(&surf_z_center, sizeof(float), 1, outfile);
	fwrite(&surf_x_min, sizeof(float), 1, outfile);
	fwrite(&surf_x_max, sizeof(float), 1, outfile);
	fwrite(&surf_y_min, sizeof(float), 1, outfile);
	fwrite(&surf_y_max, sizeof(float), 1, outfile);
	fwrite(&surf_z_min, sizeof(float), 1, outfile);
	fwrite(&surf_z_max, sizeof(float), 1, outfile);

	fwrite(&num_triangles, sizeof(int), 1, outfile);
	fwrite(triangle, sizeof(int), (num_triangles + 1) * 3, outfile);
	fwrite(tri_normal, sizeof(float), (num_triangles + 1) * 3, outfile);

	fwrite(&num_vertices, sizeof(int), 1, outfile);
	fwrite(vertex, sizeof(float), (num_triangles + 1) * 3, outfile);
	fclose(outfile);

}

void retrieve()
{
	int i = 0;
	int j = 0;
	savefile = fopen(fileName, "r");
	if (savefile == NULL)
	{
		printf("Warning: Could not open %s\n", fileName);
	}
	else
	{
		while (fscanf(savefile, "%f %f %f", &cpts[i][0], &cpts[i][1], &cpts[i][2]) != EOF)
			i++;
		ncpts = i;
		for (j = 0; j < i; j++);
		{
			printf("%f, %f, %f \n", cpts[j][0], cpts[j][1], cpts[j][2]);
		}
	}
}

/* This routine handles window resizes */
void reshape(int w, int h)
{
	width = w;
	height = h;

	/* Set the transformations */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
	glMatrixMode(GL_MODELVIEW);
	glViewport(0, 0, w, h);
}

void init()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1.3, 1.3, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	return;
}

void menu(int id)
{
	if (id == 1)
	{
		ctrlpt_flag = !ctrlpt_flag;
		glutPostRedisplay();
	}
	if (id == 2)
	{
		ctrlpolygon_ON = !ctrlpolygon_ON;
		glutPostRedisplay();
	}
	if (id == 3)
	{
		Bspline_ON = !Bspline_ON;
		glutPostRedisplay();
	}
	if (id == 4)
	{
		mouse_mode_flag = !mouse_mode_flag;
		glutPostRedisplay();
	}
	if (id == 5)
	{
		move_flag = !move_flag;
		glutPostRedisplay();
	}
	if (id == 6)
	{
		deletepoint();
		glutPostRedisplay();
	}
	if (id == 7)
	{
		type_of_surface = 0;
		glutPostRedisplay();
	}
	if (id == 8)
	{
		type_of_surface = 1;
		glutPostRedisplay();
	}
	if (id == 9)
	{
		type_of_surface = 2;
		glutPostRedisplay();
	}
	if (id == 10)
	{
		ncpts = 0;
		glutPostRedisplay();
	}
	if (id == 11);
	if (id == 12) save();
	if (id == 13)
	{
		retrieve();
		glutPostRedisplay();
	}
	if (id == 14) savemesh();
	if (id == 15) exit(0);
}

void main(int argc, char **argv)
{
	/*
	int i, j, r, c;
	for (i = 0; i<64; i++)
	{
		for (j = 0; j<64; j++)
		{
			c = ((((i & 0x8) == 0) ^ ((j & 0x8)) == 0)) * 255;
			image[i][j][0] = (GLubyte)c;
			image[i][j][1] = (GLubyte)c;
			image[i][j][2] = (GLubyte)c;
		}
	}
	*/
	/* Intialize the program */
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("curves");

	/* Register the callbacks */
	glutDisplayFunc(display);
	glutMouseFunc(mouse);

	glutCreateMenu(menu);
	glutAddMenuEntry("Control Point Switch", 1);
	glutAddMenuEntry("Control Polygon Switch", 2);
	glutAddMenuEntry("Bspline Curve Switch", 3);
	glutAddMenuEntry("Select/Insert Switch", 4);
	glutAddMenuEntry("Move the selected point", 5);
	glutAddMenuEntry("Delete the selected point", 6);
	glutAddMenuEntry("Wireframe Surface", 7);
	glutAddMenuEntry("Shade Surface", 8);
	glutAddMenuEntry("Texture Surface", 9);
	glutAddMenuEntry("Clear the Screen", 10);
	glutAddMenuEntry("Erase the Curve", 11);
	glutAddMenuEntry("Save Control Points", 12);
	glutAddMenuEntry("Retrieve Control Points", 13);
	glutAddMenuEntry("Save Mesh", 14);
	glutAddMenuEntry("Quit", 15);
	glutAttachMenu(GLUT_RIGHT_BUTTON);


	glutMotionFunc(MouseMotion);
	glutReshapeFunc(reshape);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_MAP1_VERTEX_3);
	init();
	glutMainLoop();
}