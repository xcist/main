// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

typedef long INDEX;
typedef short FLAG;
typedef int INTEGER;
typedef char * STRING;
typedef short DEGREE;
typedef double PARAMETER;

#define NR_END 1
#define SWAP(a,b) itemp=(a); (a) = (b); (b) = itemp;
#define NSTACK 50
#define FREE_ARG char*


typedef struct point
{
  float x;
  float y;
  float z;
} POINT;

typedef struct list_node{
  POINT start;
  POINT end;
  float radius;
  int ID;
  float minz, maxz, miny, maxy, minx, maxx;
  int visible;
} listelement;

/*--Definition and functions of a Binary Tree --*/
typedef struct bvh_node {
  int num;
  int *patches;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  struct bvh_node *right, *left;
} bvh_element;

typedef struct bezier_patch
{
  double cpoints[4][4][3];
  double xform[4][3];
  double xpoints[4][4][3];
  double minz, maxz, miny, maxy, minx, maxx;
  int slice_minz, slice_maxz;
} BEZIER_PATCH;

typedef struct bezier_model
{
  BEZIER_PATCH *patches;
  int num_patches;
} BEZIER_MODEL;


typedef struct qpoints
{
  INDEX n;
  POINT *Qw;
} QPOINTS;

typedef struct cpoint
{
  float  x,
         y,
         z, 
         w;
} CPOINT;

typedef struct cpolygon
{
  INDEX n;
  CPOINT *Pw;
} CPOLYGON;

typedef struct knotvector
{
    INDEX m;
    float *U;
} KNOTVECTOR;

typedef struct curve
{
  CPOLYGON pol;
  DEGREE p;
  KNOTVECTOR knt;
  float *uk;
} CURVE;

typedef struct cnet
{
  INDEX n,
        m;
  CPOINT **Pw;
} CNET;

typedef struct qnet
{
  INDEX n,
        m;
  POINT **Qw;
} QNET;

typedef struct surface
{
  CNET net;
  DEGREE p,
         q;
  KNOTVECTOR knu, 
             knv;
  float min_x, max_x, min_y, max_y, min_z, max_z;
  float centerx, centery, centerz;
  int MU_ID;
} SURFACE;

typedef struct xpoint
{
  double x;
  int organ_id;
  float costheta;
} XPOINT;

typedef struct xp_array
{
  XPOINT xp[150];
  int length;
} XP_ARRAY;

typedef struct triangle
{
  POINT vertex[3];
  float minx, maxx, miny, maxy, minz, maxz;
} TRIANGLE;

typedef struct tri_model
{
  TRIANGLE *tris;
  int num_tris;
  float min_x, max_x, min_y, max_y, min_z, max_z;
  int MU_ID;
  float density;
} TRI_MODEL;


/*NURBS_BEZ.H*/
#include <stdio.h>
#include <assert.h>

#define MAX_ORDER 20

typedef unsigned int boolean;

/* Declaration of a Nurb surface patch */
typedef struct {
  int             numU, numV;   /* #points in U and V directions */
  int             ordU, ordV;   /* order of the spline in U and V */
  float          *kU, *kV;      /* knot vectors */
                                /* length(kU) = [0...numU - 1 + ordU] */
                                /* length(kV) = [0...numV - 1 + ordV] */
  CPOINT         **points;      /* [0..numV - 1][0..numU -1] array */
} patch;

/* Structure to hold the statistics */
typedef struct stat_s {
  int count;                    /* number of patches */
  int nbezs;                    /* #Bezier patches generated */
  int tot_nbezs;                /* Total #Bez patches for the whole */
  int tot_trimpts1;             /* Max PW Trim Curve Size */
  int tot_trimpts2;             /* Max SPLINE Trim Curve Size */
                                /* file  */
} statistics_t;


/*Added to try new way*/
typedef struct line_seg
{
  double x1, x2;
  int organ_id;
} LINE_SEG;
           
/*Added to try new way*/
typedef struct seg_array
{
  int length;
  LINE_SEG sp[5000];
} SEG_ARRAY;


typedef struct priority
{
  double x;
  int organ_id;
  int next;
} PRIORITY;

typedef struct priority_array
{
  int length;
  PRIORITY pp[1000];
} PRIORITY_ARRAY;

typedef struct hit
{
  double x, x2;
  int organ_id;
  int sense;
} HIT;

typedef struct hit_array
{
  int length;
  HIT sp[10000];
} HIT_ARRAY;

