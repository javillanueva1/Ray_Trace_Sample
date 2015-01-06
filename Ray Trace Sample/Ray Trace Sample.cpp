/*******************************************************************************

Purpose:	This program renders three spheres and three planes using
non - recursive ray tracing.

The purpose of this program is to demonstrate a basic ray tracer.

Joseph Villanueva
Fall 2014
Ray Trace Sample

*******************************************************************************/

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#include <time.h>
#include <stdlib.h>

#define DEBUG   0

/*  Define some constants.  */

#define	PI	3.1415926
#define PIHALF PI/2.0
#define TWOPI PI*2.0


/* Table Size for the noise table. */

#define TS 65
#define TS_1 64

float max_marble_s = 18.0, max_marble_t = 18.0;
float noise_table[TS][TS][TS];

/* Max image size allowed. */

#define MAX_SIZE 512

/*  Define some structures.  */

struct	points	{
	float   x, y, z;
};

typedef struct	rgb_struct	{
	float   r, g, b;
} rgb;

/*  Viewing parameters.  */

struct	points	from, at, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz;
float	viewangle, angle, tanv2;
float	xinterval, yinterval;

/*  Illumination parameters.  */

float	lsx, lsy, lsz;
rgb	il, ia;
rgb	ka1, kd1, ks1;
rgb	ka2, kd2, ks2;
rgb	ka3, kd3, ks3;
rgb	ka4, kd4, ks4;
rgb	ka5, kd5, ks5;
rgb	ka6, kd6, ks6;
rgb ka7, kd7, ks7;
rgb	tka1, tkd1, tka2, tkd2;
rgb tka3, tkd3, tka4, tkd4;
rgb tka5, tkd5, tka6, tkd6;
rgb tka7, tkd7;

int	phong1, phong2, phong3, phong4, phong5, phong6, phong7;

/*  Image parameters.  */

int		xmax_pixel, ymax_pixel;

/* Image buffer.  A more efficient approach is to use one single array (texture_RGB)
rather than using the three rays.  */

float *texture_R;
float *texture_G;
float *texture_B;

/*  Object parameters.  */

float	xc1, yc1, zc1, r1, xc2, yc2, zc2, r2, xc3, yc3, zc3, r3;
float	a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
float	px1_min, px1_max, py1_min, py1_max, pz1_min, pz1_max;
float	px2_min, px2_max, py2_min, py2_max, pz2_min, pz2_max;
float	px3_min, px3_max, py3_min, py3_max, pz3_min, pz3_max;

/* Image output file. */

FILE *outpfile;

/* Mesh input file*/

FILE *meshfile;

/*Project 3 Required Compatibility Variables*/

#define MAX_VERTICES 60000		
#define MAX_TRIANGLES 60000
GLfloat vertex[MAX_VERTICES][3];
GLfloat vertex_normal[MAX_VERTICES][3];
int triangle[MAX_TRIANGLES][3];
GLfloat tri_normal[MAX_TRIANGLES][3];

float surf_x_min;
float surf_x_max;
float surf_y_min;
float surf_y_max;
float surf_z_min;
float surf_z_max;
float surf_x_center;
float surf_y_center;
float surf_z_center;
float surf_radius;
float dy, dz;

int num_triangles;
int num_vertices;

//For recursion
#define MAX_REC_NUM 2

/*******************************************************************************

Title:	Read_Information

Purpose:	This function reads in the information about the objects (three spheres
and three planes) to be rendered.  The information is assumed to be
stored in an ASCII text file called "ray.dat".

*******************************************************************************/

int Read_Information()
{
	char	str[132];
	char	filename[50];
	FILE	*inpfile;

	if ((inpfile = fopen("ray.dat", "r")) == NULL) {
		printf("ERROR: Could not open ray.dat for read!\n");
		return(0);
	}

	/*  Read in viewing information.  */

	fscanf(inpfile, "%s %f %f %f\n", str, &from.x, &from.y, &from.z);
	fscanf(inpfile, "%s %f %f %f\n", str, &at.x, &at.y, &at.z);
	fscanf(inpfile, "%s %f %f %f\n", str, &up.x, &up.y, &up.z);
	fscanf(inpfile, "%s %f\n", str, &viewangle);
	angle = viewangle * PI / 180.0;
	tanv2 = tan(angle / 2.0);
	printf("tanv2 = %f\n", tanv2);

	/* Read in the viewport information. */

	fscanf(inpfile, "%s %f %f\n", str, &VXL, &VXR);
	fscanf(inpfile, "%s %f %f\n", str, &VYB, &VYT);

	/* Read in the light vector (L).  Note: To specify a light source position, the light
	vector will need to be computed in Compute_Color().   */

	fscanf(inpfile, "%s %f %f %f\n", str, &lsx, &lsy, &lsz);

	printf("From: %f %f %f\n", from.x, from.y, from.z);
	printf("At: %f %f %f\n", at.x, at.y, at.z);
	printf("Up: %f %f %f  View Angle: %f \n", up.x, up.y, up.z, viewangle);

	/*  Read in spheres' information.  */

	fscanf(inpfile, "%s %f %f %f %f\n", str, &xc1, &yc1, &zc1, &r1);
	fscanf(inpfile, "%s %f %f %f %f\n", str, &xc2, &yc2, &zc2, &r2);
	fscanf(inpfile, "%s %f %f %f %f\n", str, &xc3, &yc3, &zc3, &r3);

	/*  Read in checker boards' (planes') information.  */

	fscanf(inpfile, "%s %f %f %f %f\n", str, &a1, &b1, &c1, &d1);
	fscanf(inpfile, "%s %f %f %f %f\n", str, &a2, &b2, &c2, &d2);
	fscanf(inpfile, "%s %f %f %f %f\n", str, &a3, &b3, &c3, &d3);

	/*  Read in the boundaries of the planes.  */

	fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px1_min, &py1_min,
		&pz1_min, &px1_max, &py1_max, &pz1_max);
	fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px2_min, &py2_min,
		&pz2_min, &px2_max, &py2_max, &pz2_max);
	fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px3_min, &py3_min,
		&pz3_min, &px3_max, &py3_max, &pz3_max);

	/*  Read in the image size.  */

	fscanf(inpfile, "%s %d %d\n", str, &xmax_pixel, &ymax_pixel);

	/* Make sure the image size does not exceed MAX_SIZE x MAX_SIZE.  */

	if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE) {
		printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
		printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
		xmax_pixel = MAX_SIZE - 1;
		ymax_pixel = MAX_SIZE - 1;
	}
	/*Read in the the mesh data file*/
	fscanf(inpfile, "%s %s\n", str, &filename);
	if ((meshfile = fopen(filename, "rb")) == NULL) {
		printf("Could not open bsplinemesh.jv for read\n");
		num_triangles = -1;
		num_vertices = -1;
	}
	else {
		printf("Read in %s\n", filename);
		// Read the triangle mesh data

		fread(&surf_x_center, sizeof(float), 1, meshfile);
		fread(&surf_y_center, sizeof(float), 1, meshfile);
		fread(&surf_z_center, sizeof(float), 1, meshfile);
		fread(&surf_x_min, sizeof(float), 1, meshfile);
		fread(&surf_x_max, sizeof(float), 1, meshfile);
		fread(&surf_y_min, sizeof(float), 1, meshfile);
		fread(&surf_y_max, sizeof(float), 1, meshfile);
		fread(&surf_z_min, sizeof(float), 1, meshfile);
		fread(&surf_z_max, sizeof(float), 1, meshfile);

		fread(&num_triangles, sizeof(int), 1, meshfile);
		fread(triangle, sizeof(int), (num_triangles + 1) * 3, meshfile);
		fread(tri_normal, sizeof(float), (num_triangles + 1) * 3, meshfile);

		fread(&num_vertices, sizeof(int), 1, meshfile);
		fread(vertex, sizeof(float), (num_vertices + 1) * 3, meshfile);
		fclose(meshfile);


		printf("Number of triangles in surface mesh: %d\n", num_triangles);
		printf("Number of vertices in surface mesh: %d\n", num_vertices);

		printf("Bounding box: (%.2f, %.2f) (%.2f, %.2f) (%.2f, %.2f)\n",
			surf_x_min, surf_x_max, surf_y_min, surf_y_max, surf_z_min, surf_z_max);
		printf("Center: (%.2f %.2f %.2f)\n", surf_x_center, surf_y_center, surf_z_center);

		// compute the radius of the bounding sphere

		surf_radius = surf_x_max - surf_x_min;
		dy = surf_y_max - surf_y_min;
		if (dy > surf_radius) surf_radius = dy;
		dz = surf_z_max - surf_z_min;
		if (dz > surf_radius) surf_radius = dz;
		printf("Bounding sphere radius : %.3f\n", surf_radius);
	}

	fclose(inpfile);

	/*  Open an output file to store the intensity values of the output image.  */

	if ((outpfile = fopen("image.out", "wb")) == NULL) {
		printf("ERROR:  cannot open image.out for write.\n");
		return(0);
	}

	/*  Allocate memory for the image buffer.  */

	texture_R = new float[xmax_pixel * ymax_pixel];
	texture_G = new float[xmax_pixel * ymax_pixel];
	texture_B = new float[xmax_pixel * ymax_pixel];
	printf("image_buf allocated.  Image size %d x %d\n", xmax_pixel, ymax_pixel);

	return(1);
}



/*******************************************************************************

Title:	Normalize

Purpose:	This function normalizes the given vector.

*******************************************************************************/

void Normalize(float *x, float *y, float *z)
{
	float	norm;

	norm = sqrt(*x * *x + *y * *y + *z * *z);
	if (norm != 0.0) {
		*x = *x / norm;
		*y = *y / norm;
		*z = *z / norm;
	}
}

/*******************************************************************************

Title:	Power

Purpose:	This function computes the power of the given base and
exponent.

*******************************************************************************/

float 	Power(float base, int exp)
{
	int	i;
	float	value;

	value = 1.0;
	for (i = 1; i <= exp; i++)
		value *= base;

	return(value);
}


/*******************************************************************************

Title:	Compute_M

Purpose:	This function computes the transformation matrix to be used
in the perspective viewing model.

*******************************************************************************/

void Compute_M()
{

	/*  Compute the line-of-sight vector, c.  */

	cx = at.x - from.x;
	cy = at.y - from.y;
	cz = at.z - from.z;
	Normalize(&cx, &cy, &cz);

	/*  Compute the cross product of vector c and the up vector.  */

	ax = cy*up.z - up.y*cz;
	ay = up.x*cz - cx*up.z;
	az = cx*up.y - up.x*cy;
	Normalize(&ax, &ay, &az);

	/*  Compute the cross product of vector a and c.  */

	bx = ay*cz - cy*az;
	by = cx*az - ax*cz;
	bz = ax*cy - cx*ay;
}

/*******************************************************************************

The following are functions which can create WoodGrain, Marble and BumpMap
textures.

*******************************************************************************/

/* Compute Wood grain texture. */

void woodgrain(float u, float v, float w, float *p_color)
{
	float r, grain;

	r = sqrt(u*u + v*v);
	if (w == 0.0)
		angle = PIHALF;
	else
		angle = atan2(u, w);

	r += 2 * sin(20.0*angle + v / 5);
	grain = (int)r % 2;
	if (grain < 0.1) {
		p_color[0] = 0.6;
		p_color[1] = 0.6;
		p_color[2] = 0.4;
	}
	else {
		p_color[0] = 0.2;
		p_color[1] = 0.2;
		p_color[2] = 0.1;
	}
}

/* Compute partials of the noise function. */

float calc_noise(float iu, float iv, int direction)
{
	int i, j, x, y, left, right;
	float noise, u, v, w, u1, v1, w1;

	i = (int)iu;
	j = (int)iv;
	x = i % TS_1;
	y = j % TS_1;

	if (direction == 1) {
		if (x <= 0)
			left = 0;
		else
			left = x - 1;

		if (x >= TS_1)
			right = TS_1;
		else
			right = x + 1;

		noise = (noise_table[right][y][0] - noise_table[left][y][0]) / 2.0;
	}
	else {
		if (y <= 0)
			left = 0;
		else
			left = y - 1;

		if (y >= TS_1)
			right = TS_1;
		else
			right = y + 1;

		noise = (noise_table[x][right][0] - noise_table[x][left][0]) / 2.0;
	}
	return(noise);
}


/* Perturb the given normal vector based on bump mapping.  */

void Bump_Map(int noise_method, float x, float y, float z, float xc, float yc, float zc, float r,
	float *nx, float *ny, float *nz)
{
	float xp, yp, zp, iu, iv, xu, yu, zu, xv, yv, zv;
	float fu, fv, a, dx, dy, dz, u, v, nnx, nny, nnz, norm_n;

	/* Translate to the orgin. */

	xp = (x - xc) / r;
	yp = (y - yc) / r;
	zp = (z - zc) / r;

	/* Convert to (u, v) coordinates. */

	v = asin(zp);
	u = atan2(yp, xp);

	/* convert to integer (iu, iv) coordinates. */

	iu = (u + PI) / TWOPI * TS_1;
	iv = (v + PIHALF) / PI * TS_1;

	/* Get the partials. */

	xu = -r * cos(v) * sin(u);
	xv = -r * sin(v) * cos(u);
	yu = r * cos(v) * cos(u);
	yv = -r * sin(v) * sin(u);
	zu = 0.0;
	zv = r * cos(v);

	/* Calculate the perturbations. */

	fu = calc_noise(iu, iv, 1);
	fv = calc_noise(iu, iv, 2);

	/* Compute D. */

	if (noise_method == 1) {
		dx = fu*xu + fv*xv;
		dy = fu*yu + fv*yv;
		dz = fu*zu + fv*zv;
	}
	else {
		nnx = *ny;
		nny = *ny;
		nnz = *nz;
		norm_n = sqrt(nnx*nnx + nny*nny + nnz*nnz);
		Normalize(&nnx, &nny, &nnz);

		/* Compute the cross product of Pu x N. */

		dx = fv*(yu*nnz - nny*zu);
		dy = fv*(nnx*zu - xu*nnz);
		dz = fv*(xu*nny - nnx*yu);

		/* Compute the cross product of Pv x N.  */

		dx += fu*(yv*nnz - nny*zv);
		dy += fu*(nnx*zv - xv*nnz);
		dz += fu*(xv*nny - nnx*yv);
	}

	/* Normalize and scale D. */

	Normalize(&dx, &dy, &dz);

	a = sqrt(fu*fu + fv*fv);
	if (noise_method == 1) {
		dx *= a;
		dy *= a;
		dz *= a;
	}
	else {
		dx *= a * norm_n;
		dy *= a * norm_n;
		dz *= a * norm_n;
	}

	*nx += dx;
	*ny += dy;
	*nz += dz;
}


/* Create a noise table.  */

void create_table()
{
	int i, j, k, ii, jj, kk;
	int val, minval = 9999.0, maxval = -9999.0;
	FILE *outfile;
	srand(time(NULL));
	rand();
	for (i = 0; i<TS_1; i++) {
		for (j = 0; j<TS_1; j++) {
			for (k = 0; k<TS_1; k++) {
				val = rand() % 255;
				noise_table[i][j][k] = val / 256.0;
			}
		}
	}

	for (i = 0; i<TS_1; i++) {
		ii = (i == TS_1) ? 0 : 1;
		for (j = 0; j<TS_1; j++) {
			jj = (j == TS_1) ? 0 : 1;
			for (k = 0; k<TS_1; k++) {
				kk = (k = TS_1) ? 0 : 1;
				noise_table[i][j][k] = noise_table[ii][jj][kk];
			}
		}
	}
}

/* Use trilinear interpolation to calculate noise at the given point.  */

float calc_noise(float px, float py, float pz)
{
	int	x, y, z;
	float noise, u, v, w, u1, v1, w1;
	int npx, npy, npz;

	npx = (int)px;
	npy = (int)py;
	npz = (int)pz;

	x = npx % TS;
	y = npy % TS;
	z = npz % TS;

	u = px - npx;
	v = py - npy;
	w = pz - npz;
	u1 = 1.0 - u;
	v1 = 1.0 - v;
	w1 = 1.0 - w;

	noise = u1*v1*w1 * noise_table[x][y][z] +
		u1*v1*w * noise_table[x][y][z + 1] +
		u1*v*w1 * noise_table[x][y + 1][z] +
		u1*v*w * noise_table[x][y + 1][z + 1] +
		u*v1*w1 * noise_table[x + 1][y][z] +
		u*v1*w * noise_table[x + 1][y][z + 1] +
		u*v*w1 * noise_table[x + 1][y + 1][z] +
		u*v*w * noise_table[x + 1][y + 1][z + 1];
	return(noise);
}

/* Evaluate the turbulence function based on the given point.  */

float turbulence(float px, float py, float pz, float pixel_size)
{

	float t, scale, npx, npy, npz;

	t = 0;
	npx = px;
	npy = py;
	npz = pz;
	for (scale = 1.0; scale>pixel_size; scale /= 2.0) {
		npx /= scale;
		npy /= scale;
		npz /= scale;

		t += calc_noise(npx, npy, npz)*scale;
	}
	return(t);
}

/* Map a scalar to a color.  */

void get_marble_color(float value, float *rgb)
{
	float x;

#if 1
	x = sqrt(value + 1.0)*0.7071;
	rgb[1] = 0.3 + 0.8*x;
	x = sqrt(x);
	rgb[0] = 0.3 + 0.6*x;
	rgb[2] = 0.6 * 0.4*x;
#endif

#if 0
	/*  Yellow */

	x = sqrt(value + 1.0);
	rgb[1] = 0.3 + 0.7 * x;
	x = sqrt(x);
	rgb[0] = 0.3 + 0.5*x;
	rgb[2] = 0.6 * 0.4*x;
#endif


#if 0

	x = sqrt(value + 1.0);
	//rgb[1] = 0.1 + 0.3 * x;
	rgb[2] = 0.6 * 0.4*x;
	x = sqrt(x);
	rgb[0] = 0.3 + 0.6*x;
	//rgb[2] = 0.7 * 0.2*x;
	rgb[1] = 0.3 + 0.8 * x;
#endif
}

void marble_y(float px, float py, float pz, float pixel_size, float *rgb)
{
	float value, t;

	t = 3.0*turbulence(px, py, pz, pixel_size);
	value = py + t;
	get_marble_color((float)sin(value*3.1415926), rgb);
}

/* Create a turbulence noise.  */

void marble_texture(float s, float t, float *p_color)
{
	float u, v;
	float scale_factor = 5.0;
	float pixel_factor = 175.0;
	float px, py, pz;

	pz = 0.0;
	u = s / px1_max;
	v = t / py1_max;
	px = u*scale_factor;
	py = v*scale_factor;
	marble_y(px, py, pz, 1.0 / pixel_factor, p_color);
}

/*******************************************************************************

Title:	Setup_Parameters

Purpose:	This function sets up the necessary parameters for
performing the ray trace.  It first computes the
transformation matrix for the perspective viewing model, then
sets up the default illumination parameters.

*******************************************************************************/

void Setup_Parameters()
{

	/*  Compute the transformation matrix for converting world coordinates to eye
	coordinates.  */

	Compute_M();

	/*  Normalized the given directional light vector.

	Note:  DO NOT normalize this vector, if the position of the light source is given.
	The light vector (L) would need to computed for each intersction point,
	then normalized in Compute_Color().   */

	Normalize(&lsx, &lsy, &lsz);
	printf("light position %f %f %f\n", lsx, lsy, lsz);

	/*  Set up the conversion factors for converting from pixel coordinates to
	view port coordinates.  */

	xinterval = (VXR - VXL) / xmax_pixel;
	yinterval = (VYT - VYB) / ymax_pixel;

	/*  Set up default illumination (Phong lighting) parameters.  */

	il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
	ia.r = 5.0;	ia.g = 5.0;	ia.b = 5.0;

	/*  Phone lighting parameters for the three spheres.  */

	ka1.r = 0.3;	ka1.g = 0.0;	ka1.b = 0.0;
	kd1.r = 0.7;	kd1.g = 0.0;	kd1.b = 0.0;
	ks1.r = 1.0;	ks1.g = 1.0;	ks1.b = 1.0;
	tka1.r = 0.2;	tka1.g = 0.0;	tka1.b = 0.0;
	tkd1.r = 0.2;	tkd1.g = 0.0;	tkd1.b = 0.0;
	phong1 = 60;

	ka2.r = 0.0;	ka2.g = 0.3;	ka2.b = 0.0;
	kd2.r = 0.0;	kd2.g = 0.7;	kd2.b = 0.0;
	ks2.r = 1.0;	ks2.g = 1.0;	ks2.b = 1.0;
	tka2.r = 0.0;	tka2.g = 0.2;	tka2.b = 0.0;
	tkd2.r = 0.0;	tkd2.g = 0.2;	tkd2.b = 0.0;
	phong2 = 90;

	ka3.r = 0.0;	ka3.g = 0.0;	ka3.b = 0.3;
	kd3.r = 0.0;	kd3.g = 0.0;	kd3.b = 0.7;
	ks3.r = 1.0;	ks3.g = 1.0;	ks3.b = 1.0;
	tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.2;
	tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.2;
	phong3 = 120;

	/*  Phone lighting parameters for the three planes (not shown).  */

	ka4.r = 0.1;	ka4.g = 0.1;	ka4.b = 0.0;
	kd4.r = 0.7;	kd4.g = 0.7;	kd4.b = 0.0;
	ks4.r = 1.0;	ks4.g = 1.0;	ks4.b = 1.0;
	tka4.r = 0.1;	tka4.g = 0.0;	tka4.b = 0.0;
	tkd4.r = 0.7;	tkd4.g = 0.0;	tkd4.b = 0.0;
	phong4 = 120;

	ka5.r = 0.1;	ka5.g = 0.0;	ka5.b = 0.1;
	kd5.r = 0.7;	kd5.g = 0.0;	kd5.b = 0.7;
	ks5.r = 1.0;	ks5.g = 1.0;	ks5.b = 1.0;
	tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.1;
	tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;
	phong5 = 120;

	ka6.r = 0.1;	ka6.g = 0.1;	ka6.b = 0.1;
	kd6.r = 0.8;	kd6.g = 0.6;	kd6.b = 0.2;
	ks6.r = 1.0;	ks6.g = 1.0;	ks6.b = 1.0;
	tka6.r = 0.0;	tka6.g = 0.1;	tka6.b = 0.0;
	tkd6.r = 0.0;	tkd6.g = 0.7;	tkd6.b = 0.0;
	phong6 = 120;

	/*Phong lighting parameters for the triangle mesh*/
	ka7.r = 0.0;	ka7.g = 0.0;	ka7.b = 1.0;
	kd7.r = 0.0;	kd7.g = 0.0;	kd7.b = 1.0;
	ks7.r = 0.0;	ks7.g = 0.0;	ks7.b = 1.0;
	tka7.r = 1.0;	tka7.g = 1.0;	tka7.b = 0.0;
	tkd7.r = 1.0;	tka7.g = 1.0;	tkd7.b = 0.0;
	phong7 = 80;

}




/*******************************************************************************

Title:	Check_Sphere

Purpose:	This function determines if the give ray intercepts the given
sphere.

*******************************************************************************/

void Check_Sphere(float px, float py, float pz, float dx, float dy, float dz, float xc, float yc, float zc, float r, float *t1, float *t2)
{
	float	a, b, c, xdiff, ydiff, zdiff, discr;

	xdiff = px - xc;
	ydiff = py - yc;
	zdiff = pz - zc;
	a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
	b = 2.0*(dx*xdiff + dy*ydiff + dz*zdiff);
	c = dx*dx + dy*dy + dz*dz;

	/*  Check if there are any intersections.  */

	discr = b*b - 4.0*a*c;
	if (discr < 0.0) {
		*t1 = -1.0;
		*t2 = -1.0;
	}
	else if (discr == 0.0) {
		*t1 = -b / (2.0*c);
		*t2 = -1.0;
	}
	else {
		discr = sqrt(discr);
		*t1 = (-b + discr) / (2.0*c);
		*t2 = (-b - discr) / (2.0*c);
	}
}


/*******************************************************************************

Title:	Check_Plane

Purpose:	This function checks if the given ray intercepts the given
plane.

*******************************************************************************/

void Check_Plane(float px, float py, float pz, float dx, float dy, float dz, float a, float b, float c, float d, float *t1)
{
	*t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}

/********************************************************************************

Title: Dot

Purpose:	This function computes the dot product
********************************************************************************/
// Compute the dot project of v0 and v1

float dot(float v0[3], float v1[3])
{
	return(v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]);
}

/*******************************************************************************

Title: Check_Triangle

Purpose: This function will, given a triangle from the mesh, determine if an
intersection point is inside the triangle

Note: Given the triangle number (triNumber), the point (px,py,pz),
the ray described by (dx,dy,dz), and the equation of the plane described by (a,b,c,d).
Return t1.  Set t1 = -1 if the ray does not intersect the plane or
the intersection point is not inside the triangle triNumber.
*******************************************************************************/
void Check_Triangle(int triNumber, float px, float py, float pz, float dx, float dy, float dz, float a, float b, float c, float d, float *t1)
{
	int	 vert1, vert2, vert3;
	float dot00, dot01, dot02, dot11, dot12, invDenom;
	float t, ipx, ipy, ipz, u, v, p[3][3], v0[3], v1[3], v2[3];

	t = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
	if (t > 0.00001) {

		// Check if the intersection point is inside the triangle.

		ipx = px + t*dx;
		ipy = py + t*dy;
		ipz = pz + t*dz;

		vert1 = triangle[triNumber][0];
		vert2 = triangle[triNumber][1];
		vert3 = triangle[triNumber][2];

		p[0][0] = vertex[vert1][0];
		p[0][1] = vertex[vert1][1];
		p[0][2] = vertex[vert1][2];

		p[1][0] = vertex[vert2][0];
		p[1][1] = vertex[vert2][1];
		p[1][2] = vertex[vert2][2];

		p[2][0] = vertex[vert3][0];
		p[2][1] = vertex[vert3][1];
		p[2][2] = vertex[vert3][2];

		v0[0] = p[2][0] - p[0][0];
		v0[1] = p[2][1] - p[0][1];
		v0[2] = p[2][2] - p[0][2];

		v1[0] = p[1][0] - p[0][0];
		v1[1] = p[1][1] - p[0][1];
		v1[2] = p[1][2] - p[0][2];

		v2[0] = ipx - p[0][0];
		v2[1] = ipy - p[0][1];
		v2[2] = ipz - p[0][2];

		dot00 = dot(v0, v0);
		dot01 = dot(v0, v1);
		dot02 = dot(v0, v2);
		dot11 = dot(v1, v1);
		dot12 = dot(v1, v2);

		invDenom = 1.0 / (dot00*dot11 - dot01*dot01);

		u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		if (u >= 0.0 && v >= 0.0 && u + v <= 1.0)
			*t1 = t;
		else
			*t1 = -1;
	}
	else
		*t1 = -1;
}
/*******************************************************************************

Title:	Check_Shadow

Purpose:	This funciton will check the shadow ray, a ray from the shadow to the
eye source, to see if this ray intersects any other object. It will return shadow
hit or not hit, which will then be used in Compute_Color(...). In Compute_Color(...)
if shadow is hit, only use ambient lighting to determine color.

Note: Given the point (px, py, pz) and the ray (dx, dy, dz), Check other objects.
If other objects intersect, return 1, otherwise return 0.
********************************************************************************/

int Check_Shadow(float px, float py, float pz, float nx, float ny, float nz)
{
	//Compute Shadow Ray
	float dirx, diry, dirz;
	dirx = lsx - px;
	diry = lsy - py;
	dirz = lsz - pz;
	Normalize(&dirx, &diry, &dirz);

	//Generate Ray
	float dx, dy, dz;
	dx = nx*dirx;
	dy = ny*diry;
	dz = nz*dirz;

	float t_min, t1, t2;
	t_min = 999.0;
	/*  Check if the current ray intercepts sphere #1.  */

	Check_Sphere(px, py, pz, dx, dy, dz, xc1,
		yc1, zc1, r1, &t1, &t2);

	if (t1 > 0.0 && t1 < t_min || t2 > 0.0 && t2<t_min) {
		return 1;
	}

	/*  Check if the current ray intercepts sphere #2.  */

	Check_Sphere(px, py, pz, dx, dy, dz, xc2,
		yc2, zc2, r2, &t1, &t2);

	if ((t1 > 0.0 && t1<t_min) || (t2 > 0.0 && t2<t_min)) {
		return 1;
	}

	/*  Check if the current ray intercepts sphere #3.  */

	Check_Sphere(px, py, pz, dx, dy, dz, xc3,
		yc3, zc3, r3, &t1, &t2);

	if ((t1 > 0.0 && t1<t_min) || (t2 > 0.0 && t2 < t_min)) {
		return 1;
	}

	/*  Check if the current ray intercepts plane #1.  */

	Check_Plane(px, py, pz, dx, dy, dz, a1, b1, c1, d1, &t1);

	if (t1 > 0.0 && t1<t_min) {
		/*  Check if the intersection point is inside the min/max values. */
		return 1;
	}

	/*  Check if the current ray intercepts plane #2.  */

	Check_Plane(px, py, pz, dx, dy, dz, a2, b2, c2, d2, &t1);

	if (t1 > 0.0 && t1<t_min) {
		return 1;
	}

	/*  Check if the current ray intercepts plane #3.  */

	Check_Plane(px, py, pz, dx, dy, dz, a3, b3, c3, d3, &t1);

	if (t1 > 0.0 && t1 < t_min) {
		return 1;
	}

	Check_Sphere(from.x, from.y, from.z, dx, dy, dz, surf_x_center, surf_y_center, surf_z_center, surf_radius, &t1, &t2);
	if (t1 > 0.0 && t1 < t_min || t2 > 0.0 && t2 < t_min) {
		for (int i = 0; i<num_triangles; i++) {
			int v0 = triangle[i][0];

			float a = tri_normal[i][0];
			float b = tri_normal[i][1];
			float c = tri_normal[i][2];


			float d0 = -a * vertex[v0][0] - b * vertex[v0][1] - c * vertex[v0][2];
			// Check if the ray intersects the plane containing the triangle.

			Check_Triangle(i, from.x, from.y, from.z, dx, dy, dz, a, b, c, d0, &t1);
			if (t1 > 0.0 && t1 < t_min) {

				return 1;
				//printf("intersect triangle %d   t_min = %f\n", tri_num, t_min);
			}
		}
	}
	return 0;
}
/*******************************************************************************

Title:	Compute_Intersection

Purpose:	This function computes the intersection of ray with an
object.  The intersection point is given by a parametric value
t, where ray = p + d*t, d = the direction of the ray, and p is
the starting point of the ray.

*******************************************************************************/

void Compute_Intersection(float px, float py, float pz, float dx, float dy, float dz, float t, float *newx, float *newy, float *newz)
{
	*newx = px + t*dx;
	*newy = py + t*dy;
	*newz = pz + t*dz;
}


/*******************************************************************************

Title:	Compute_Color

Purpose:	This function computes the intensity of the color for the
given location based on the Phong lighting model.

*******************************************************************************/

void Compute_Color(int shadow_flag, float ipx, float ipy, float  ipz, float  nx, float  ny, float  nz,
	rgb ia, rgb ka, rgb kd, rgb ks, int n, float *r, float *g, float *b)
{
	float	vx, vy, vz;
	float	rx, ry, rz;
	float	lightx, lighty, lightz;
	float	ndotl, vdotr, cosalphapower;

	/*  Compute the view vector.  */

	vx = from.x - ipx;
	vy = from.y - ipy;
	vz = from.z - ipz;
	Normalize(&vx, &vy, &vz);

	/* Compute the light(L) vector */
	lightx = abs(ipx - lsx);
	lighty = abs(ipy - lsy);
	lightz = abs(ipz - lsz);

	Normalize(&lightx, &lighty, &lightz);

	/*  Compute the R (reflection) vector.  */

	ndotl = nx*lightx + ny*lighty + nz*lightz;
	rx = 2.0*ndotl*nx - lsx;
	ry = 2.0*ndotl*ny - lsy;
	rz = 2.0*ndotl*nz - lsz;

	/* Compute the V (view) vector. */

	vdotr = vx*rx + vy*ry + vz*rz;

	/* Compute Ia * Ka.  */

	*r = ia.r * ka.r;
	*g = ia.g * ka.g;
	*b = ia.b * ka.b;

	/* Compute diffuse reflection. */

	if (ndotl >= 0.0 && shadow_flag == 0) {

		/*  diffuse reflection = kd * N dot L * Il  */

		*r = *r + kd.r*ndotl*il.r;
		*g = *g + kd.g*ndotl*il.g;
		*b = *b + kd.b*ndotl*il.b;

		if (vdotr >= 0.0) {

			/*  specular reflection = ks * cos(alpha)**K^n * Il */

			cosalphapower = Power(vdotr, n);

			*r = *r + ks.r*cosalphapower*il.r;
			*g = *g + ks.g*cosalphapower*il.g;
			*b = *b + ks.b*cosalphapower*il.b;
		}
	}

	/*  Make sure that the color is within range.  */

	if (*r > 1.0) *r = 1.0;
	if (*g > 1.0) *g = 1.0;
	if (*b > 1.0) *b = 1.0;
}

/*******************************************************************************

Title: Check_intersection

Purpose: This will take in a ray (dx, dy, dz) and return numbers by putting them
into pointers

*******************************************************************************/
void check_intersection(float *t_min, float *ipx, float *ipy, float *ipz, int *obj_num, int *tri_num, float dx, float dy, float dz){
	float	t1, t2;
	Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc1,
		yc1, zc1, r1, &t1, &t2);

	if (t1 >= 0.0) {
		*t_min = t1;
		*obj_num = 1;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);
	}

	if (t2 >= 0.0 && t2<*t_min) {
		*t_min = t2;
		*obj_num = 1;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t2, ipx, ipy, ipz);
	}

	Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc2,
		yc2, zc2, r2, &t1, &t2);

	if (t1 >= 0.0 && t1<*t_min) {
		*t_min = t1;
		*obj_num = 2;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);
	}

	if (t2 >= 0.0 && t2<*t_min) {
		*t_min = t2;
		*obj_num = 2;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t2, ipx, ipy, ipz);
	}


	Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc3,
		yc3, zc3, r3, &t1, &t2);

	if (t1 >= 0.0 && t1<*t_min) {
		*t_min = t1;
		*obj_num = 3;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);
	}

	if (t2 >= 0.0 && t2<*t_min) {
		*t_min = t2;
		*obj_num = 3;
		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t2, ipx, ipy, ipz);
	}

	Check_Plane(from.x, from.y, from.z, dx, dy, dz, a1, b1, c1, d1, &t1);

	if (t1 >= 0.0 && t1<*t_min) {

		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);

		if (*ipx >= px1_min && *ipx <= px1_max &&
			*ipy >= py1_min  && *ipy <= py1_max &&
			*ipz >= pz1_min && *ipz <= pz1_max) {

			*t_min = t1;
			*obj_num = 4;
		}
	}

	Check_Plane(from.x, from.y, from.z, dx, dy, dz, a2, b2, c2, d2, &t1);

	if (t1 >= 0.0 && t1<*t_min) {
		/*  Check if the intersection point is inside the min/max values. */

		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);

		if (*ipx >= px2_min && *ipx <= px2_max &&
			*ipy >= py2_min  && *ipy <= py2_max &&
			*ipz >= pz2_min && *ipz <= pz2_max) {

			*t_min = t1;
			*obj_num = 5;
		}
	}

	/*  Check if the current ray intercepts plane #3.  */

	Check_Plane(from.x, from.y, from.z, dx, dy, dz, a3, b3, c3, d3, &t1);

	if (t1 >= 0.0 && t1 < *t_min) {
		/*  Check if the intersection point is inside the min/max values. */

		Compute_Intersection(from.x, from.y, from.z,
			dx, dy, dz, t1, ipx, ipy, ipz);

		if (*ipx >= px3_min && *ipx <= px3_max &&
			*ipy >= py3_min  && *ipy <= py3_max &&
			*ipz >= pz3_min && *ipz <= pz3_max) {

			*t_min = t1;
			*obj_num = 6;
		}
	}
	Check_Sphere(from.x, from.y, from.z, dx, dy, dz, surf_x_center, surf_y_center, surf_z_center, surf_radius, &t1, &t2);
	if (t1 > 0.0 && t1 < *t_min || t2 > 0.0 && t2 < *t_min) {
		for (int i = 0; i<num_triangles; i++) {
			int v0 = triangle[i][0];

			float a = tri_normal[i][0];
			float b = tri_normal[i][1];
			float c = tri_normal[i][2];


			float d0 = -a * vertex[v0][0] - b * vertex[v0][1] - c * vertex[v0][2];
			// Check if the ray intersects the plane containing the triangle.

			Check_Triangle(i, from.x, from.y, from.z, dx, dy, dz, a, b, c, d0, &t1);
			if (t1 > 0.0 && t1 < *t_min) {

				Compute_Intersection(from.x, from.y, from.z,
					dx, dy, dz, t1, ipx, ipy, ipz);

				*t_min = t1;
				*tri_num = i;
				*obj_num = 7;
				//printf("intersect triangle %d   t_min = %f\n", tri_num, t_min);
			}
		}
	}

}

/*******************************************************************************

Title: Recurse

Purpose: This function will perform the recursive part of the ray tracer
*******************************************************************************/
void Recurse(int depth, float *r, float *g, float *b, float dx, float dy, float dz){
	int	    xp, yp, obj_num, shadow_flag;
	int	    texture;
	float	t_min, t1, t2, ipx, ipy, ipz;
	float   nx, ny, nz;
	float newr;
	float newg;
	float newb;
	int rflag = 0;
	float color[3];
	int	tri_num;
	newr = 0.0;
	newg = 0.0;
	newb = 0.0;
	if (depth == -1){
		return;
	}
	else {
		t_min = 999.0;
		obj_num = 0;
		texture = 0;
		check_intersection(&t_min, &ipx, &ipy, &ipz, &obj_num, &tri_num, dx, dy, dz);

		switch (obj_num) {

		case 0: *r = 0.0;
			*g = 0.7;
			*b = 0.8;
			return;
			break;

			/*  The current ray intercept sphere #1.  */

		case 1:
			nx = ipx - xc1;
			ny = ipy - yc1;
			nz = ipz - zc1;
			woodgrain(ipx, ipy, ipz, color);
			*r = color[0];
			*g = color[1];
			*b = color[2];
			return;
			break;

			/*  The current ray intercepts sphere #2.  */

		case 2:
			nx = ipx - xc2;
			ny = ipy - yc2;
			nz = ipz - zc2;
			Normalize(&nx, &ny, &nz);
			shadow_flag = 0;// Check_Shadow(ipx, ipy, ipz, nx, ny, nz);

			Compute_Color(shadow_flag, ipx, ipy, ipz, normaliznx, ny, nz, ia, ka2, kd2, ks2, phong2, &newr, &newg, &newb);
			rflag = 1;
			break;

			/*  The current ray intercepts sphere #3.  */

		case 3:
			nx = ipx - xc3;
			ny = ipy - yc3;
			nz = ipz - zc3;
			Normalize(&nx, &ny, &nz);
			Bump_Map(1, ipx, ipy, ipz, xc3, yc3, zc3, r3, &nx, &ny, &nz);
			shadow_flag = 0;// Check_Shadow(ipx, ipy, ipz, nx, ny, nz);

			texture = 1;

			Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &newr, &newg, &newb);

			break;

			/*  The current ray intercepts wood board #1.  */

		case 4:
			nx = a1;
			ny = b1;
			nz = c1;
			shadow_flag = Check_Shadow(ipx, ipy, ipz, nx, ny, nz);
			// Use wood grain as the color
			woodgrain(ipx, ipy, ipz, color);
			*r = color[0];
			*g = color[1];
			*b = color[2];
			return;
			break;

			/*  The current ray intercepts checker board #2.  */

		case 5:
			nx = a2;
			ny = b2;
			nz = c2;
			shadow_flag = Check_Shadow(ipx, ipy, ipz, nx, ny, nz);

			// Use marble texture as the color
			marble_texture(ipy, ipz, color);

			*r = color[0];
			*g = color[1];
			*b = color[2];
			return;

			break;

			/*  The current ray intercepts plane #3.  */

		case 6:
			nx = a3;
			ny = b3;
			nz = c3;

			shadow_flag = Check_Shadow(ipx, ipy, ipz, nx, ny, nz);
			if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
				if ((ipz >= 2.0 && ipz < 4.0) || (ipz >= 6.0))
					texture = 1;
				else
					texture = 0;
			}
			else {
				if ((ipz < 2.0) || (ipz >= 4.0 && ipz < 6.0))
					texture = 1;
				else
					texture = 0;
			}
			if (texture == 1)
				Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &newr, &newg, &newb);
			else
				Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &newr, &newg, &newb);
			break;

		case 7:
			nx = tri_normal[tri_num][0];
			ny = tri_normal[tri_num][1];
			nz = tri_normal[tri_num][2];
			shadow_flag = 0;// Check_Shadow(ipx, ipy, ipz, nx, ny, nz);

			texture = 0;
			Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka7, kd7, ks7, phong7, &newr, &newg, &newb);
			break;
		}
	}

	if (rflag){
		float rrx, rry, rrz;
		float ndoti;
		float drx, dry, drz;
		*r = newr*Power(0.1, (1 - depth)) + *r;
		*g = newg*Power(0.1, (1 - depth)) + *g;
		*b = newb*Power(0.1, (1 - depth)) + *b;


		rrx = abs(ipx - nx);
		rry = abs(ipy - ny);
		rrz = abs(ipz - nz);


		ndoti = nx*rrx + ny*rry + nz*rrz;

		drx = ipx - 2 * ndoti*ipx;
		dry = ipy - 2 * ndoti*ipy;
		drz = ipz - 2 * ndoti*ipz;

		Recurse(depth - 1, r, g, b, drx, dry, drz);

	}
	else {
		*r = newr;
		*b = newb;
		*g = newg;
		return;
	}
}

/*******************************************************************************

Title:	Ray_Trace

Purpose:	This function performs simple ray tracing.  This is a non-recursive
ray tracing without any reflection and refraction rays.

*******************************************************************************/

void Ray_Trace() /*should accept ray and depth for recursive purposes*/
{
	printf("ray tracing...\n");
	int	    xp, yp, obj_num, shadow_flag;
	int		rflag; //If 1, surface is reflective, if 0 then not reflective
	int	    texture, buf_ptr;
	int     num_image = 1;
	float	xv, yv, dx, dy, dz, nx, ny, nz;
	float	t_min, t1, t2, ipx, ipy, ipz;
	float	r, g, b;
	float   u, v;
	int		tri_num;
	float	ndotl, rx, ry, rz;

	/*  Generate a ray for each pixel in the desired image.  */

	buf_ptr = 0;
	for (xp = 0; xp < xmax_pixel; xp++) {
		u = (float)xp / xmax_pixel;

		for (yp = 0; yp < ymax_pixel; yp++) {
			v = (float)yp / ymax_pixel;

			/*  Compute the corresponding view port coordinates.  */

			xv = VXL + xp * xinterval;
			yv = VYB + yp * yinterval;

			/*  Compute the direction of the current ray from the "From" point to the
			current position on the image.  */
			dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
			dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
			dz = az*xv*tanv2 + bz*yv*tanv2 + cz;

			/*Begin depth trace at depth 1*/
			r = 0.0; g = 0.0; b = 0.0;
			Recurse(1, &r, &g, &b, dx, dy, dz);

			/* Save the computed color intensity to the image buffer. */
			//if (depth == MAX_REC_NUM)
			texture_R[xp + xmax_pixel * yp] = r;
			texture_G[xp + xmax_pixel * yp] = g;
			texture_B[xp + xmax_pixel * yp] = b;
		}
	}
	/*  Write the image to the output file.  */
	printf("Writing to image...\n");
	fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
	fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
	fwrite(texture_R, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_G, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_B, sizeof(float), xmax_pixel*ymax_pixel, outpfile);

	fclose(outpfile);
}


/* Initialize the projection matrix.  */

void myinit(void)
{
	/* attributes */

	glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */

	/* set up viewing */
	/* 512 x 512 window with origin lower left */

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 512.0, 0.0, 512.0);
	glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image.   A more efficient method is
to use glDrawPixels(). */

void display(void)
{
	int s, t;
	float  r, g, b;

	glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */

	for (t = 0; t < ymax_pixel; t++) {
		for (s = 0; s < xmax_pixel; s++) {

			r = texture_R[s + xmax_pixel * t];
			g = texture_G[s + xmax_pixel * t];
			b = texture_B[s + xmax_pixel * t];

			glColor3f(r, g, b);
			glBegin(GL_POINTS);
			glVertex2f(s, t);
			glEnd();
		}
	}

	glFlush(); /* clear buffers */
}

/*  Main routine.  */

int main(int argc, char**argv)
{

	Read_Information();
	Setup_Parameters();
	create_table();
	Ray_Trace();

	/* Standard GLUT initialization */

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); /* default, not needed */
	glutInitWindowSize(500, 500); /* 500 x 500 pixel window */
	glutInitWindowPosition(0, 0); /* place window top left on display */
	glutCreateWindow("Ray Trace"); /* window title */
	glutDisplayFunc(display); /* display callback invoked when window opened */

	myinit(); /* set attributes */
	glutMainLoop(); /* enter event loop */

	return(0);
}
