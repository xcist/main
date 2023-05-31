// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __YYcode_h__
#define __YYcode_h__

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define VERY_BIG 1e300
#define MATERIAL_INDEX_BASIS 1  //can be 1 or zero
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <pthread.h>
#include <unistd.h>

// Global Variable Definitions /////////////

int COMPARISON_NUMBER = 0;

//Parameters for Accurate Detector Model, Mingye
int Accurate_Detector_Model_is_ON = 0;
int n_col_oversample = 1;
int n_row_oversample = 1;
int n_col_oversample_add_xtalk = 1;
int n_row_oversample_add_xtalk = 1;

struct module_info
{
  double *Height;
  double *Width;
  int *Pix;
  double *Coords;
  int *Sub;
  double *Sampling;
  double *Weight;
  double *sourceWeights;
  int nSubSources;
  int maxPixPerModule;
  int maxSubDets;
  int moduleOverlapType;
};

struct module_info modules = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0};

struct pnt
{
  double x;
  double y;
  double tan_angle;
  int index;
};

struct bounding_info
{
  int *vertexStartIndex;
  double *vertexLocations;
};

struct bounding_info bounding = {NULL,NULL};

struct phantom_info
{
  int numObjects;
  int *objectType;
  double *objectCenter;
  double *shape;
  double *Qmatrix;
  int *clipStartIndex;
  double *clipNormalVector;
  double *clipDistance;
  int *nClipPlanes;
  double *density;
  int *materialInd;
  double xbounds[2];
};

struct phantom_info phantom = {0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,{0,0}};

struct material_info
{
  int materialCount;
  int eBinCount;
  double *muTable;
};

struct material_info materials = {0,0,NULL};

struct projector_args
{
  double *Paras; //parameters for Accurate Detector Model, Mingye
  double subviewWeight;
  double *thisView;
  double *sourcePoints;
  int nSubSources;
  double *srcHullPoints;
  int nSrcHullPoints;
  int *firstDetIndex;
  int nModulesIn;
  int *modTypeInds;
  double *Up;
  double *Right;
  double *Center;
  int UNUSED_tvLength;  
};

struct height_lims
{
  double min;
  double max;
};

struct box_lims
{
  double minR;
  double maxR;
  double minU;
  double maxU;
};

double limitLoadAverage = 0;
int nextModuleInQ;
int modulesInQ;
pthread_mutex_t QLock = PTHREAD_MUTEX_INITIALIZER;
int thread_count;
pthread_t *t_id = NULL;

static int debug_flag = 1;

struct indexed_list
{
  double value;
  double originalIndex;
};

//////////////////////////////////////////////

static void dbug(int level,const char *s,...);

double solve_cubic(double *a);
    
void cross(double *vec1, double *vec2, double *vecOut);

double magnitude(double r, double i);

void sqrtm(double in_r,double in_i,double *out_r,double *out_i);

void complex_multiply(double r1, double i1, double r2, double i2,double *real_out, double *imag_out);

void solve_cubic_all(double *a,double *zr,double *zi);

int compare_doubles (const void *a, const void *b);

int compIL (const void *a, const void *b);

int solve_quartic2(double *a,double *z);

int solve_quartic(double *a,double *z);

double quadratic_form(double *vec1,double *matrix, double *vec2);

int quartic_intersect(double *a0,double *alpha,double *tc,int obj);

int clip_all(double *a,double *alpha,double rayLength,double *tc2,int out,double *st_list,double *en_list,double *den_list,int *pri_list,int *mat_list,int num_int,int i);

int quadratic_intersect(double *a0,double *alpha,int pars11,double *tc2,int obj);

int* my_memcpyi(int *src, int *dest, int bytes);

double* my_memcpyd(double *src, double *dest, int bytes);

void set_module_info(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int moduleOverlapType, double *sourceWeights, int nSubSources);

void set_bounding_info(int numObjs, int *vertexStartInd, double *vertLocs, int numVerts);

void set_phantom_info(int numObjs, int *objType, int *clipStInd, int *nPlanes, int *matInd, double *objCent, double *shp, double *Qmat, double *clipNormVec, double *clipDist, double *den, int totalNumPlanes);

void set_material_info(int materialCount, int eBinCount, double *muTable);

void intersections_full_list(int *objlist, int lenObjList,double *a,double *alpha,double rayLength,double *t_ends,int *matls,double *dens,int *num_segs);

void make_vol(float *volume, int Nx, double xoff, double dx, int Ny, double yoff, double dy, int Nz, double zoff, double dz, int oversampling);

void intersections(int *objlist,int lenObjList,double *a,double *alpha,double rayLength,double *mb);

int compare_pts (const void *a, const void *b);

int compute_convex_hull_2d(double* x,double* y,int pts,int* list);

void crop_polygon_FMtest(double *x, double *y, int *indexList, int *listLength, double height, double width, int projPoints);

void crop_polygon(double *x, double *y, int *indexList, int *listLength, double height, double width, int projPoints);

void store_height_lims(double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries);

void store_half_planes(double *pixel_R, double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries);

void store(double *pixel_R, double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries);

void compute_object_projections(double *srcHullPoints, int nSrcHullPoints, int moduleTypeIndex, double *up, double *right, double *center, void *boundaries);

int any_objects_1(int moduleNumber, void *boundaries);

int any_objects_2(int moduleNumber, void *boundaries);

void build_object_list1(double *pix_vlims, int *objectList, int *n_objlist, int moduleNumber, double v, void *boundaries);

void build_object_list2(double *pix_vlims, double *pix_ulims, int *objectList, int *n_objlist, int moduleNumber, double *uv, void *boundaries);

void build_object_list3();

void intersections_loop(double *detCenter, double *right, double *up, double *sampling, int nSubDets, double *sourcePoints, int *objectList, int nListObjects, double *thisView, int detIndex, double subviewWeight, double *detWeights);

void set_Accurate_Detector_Model(double *Paras); //Set parameters for Accurate Detector Model, Mingye

void Projector(double *Paras, double subviewWeight, double *thisView, double *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength);

void *projector_wrapper(void *pointerIn);

void Projector_threaded(double *Paras, double subviewWeight, double *thisView, double *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength, int numThreads, double maxLoadAverage);

void diff(double *d1,double *d2,int len);

void CO(double *c_values);

void COf(float *c_values);

void COi(int *c_values);

void COs(short *c_values);

#endif

