// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <pthread.h>
#include "ct_nurbs.h"
#define MATERIAL_INDEX_BASIS 0 // can be 1 or zero

#include <unistd.h>
#include <ctype.h>

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])
#define MAXIMUM_NUMBER_OF_MATERIALS 50
#define EPSILON 0.000001
#define EPSILON2 0.0001

// What is taking the bulk of the time?
//  For a view of 2.7 seconds, I found that:
//     - The three main functions in Intersect_bez (Test_patch, subdivide_patch, and Test_extents2) account for 1.6 seconds (59%)
//             this was measured by adding Intersect_dum, which added 1.6 seconds to the time
//             much of this was found to be in computing the x, y, and z patch bounds in Test_extents2
//     - The Segm_inside stuff now takes very little time since it is called rarely (0.3 seconds) (11%)
//             this was determined by calling Segm_inside function twice each time
//     - The final thisView calculations (including Break_segment2 and Calc_line_int2) take 0.1 seconds (4%)
//             this was measured by commenting out this section (the time dropped from 2.7 to 2.6)
//             This is much better than it was before with the original Break_segment() and Calc_line_int()
//     - The other 0.7 seconds (26%) must be spent somewhere (likely the remaining part of Intersect_bez accounts for a good part of this)

static double mdps[362][3]; // maximally dispersed points on a sphere
static long pval_cnt = 0, pinv_cnt = 0;
static double tol1 = 0.001, tol2 = 0.0001, pad = 0.02;
// #define TOL1 0.001    //original value: 0.001
// #define TOL2 0.0001    //original value: 0.0001
// #define PAD 0.02    //original value: 0.5
//  Based on some quick experiments, PAD should be at least 10x the larger of TOL1 and TOL2
//  a lower number makes it possible to get large errors since rays can pass through gaps between the bezier subpatches
//  a higher number increases the amount of speckle for near-tangent rays

int *my_memcpyi(int *src, int *dest, int bytes);
double *my_memcpyd(double *src, double *dest, int bytes);

struct module_info_NCAT
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
};

struct projector_args
{
  double subviewWeight;
  double *thisView;
  float *sourcePoints;
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

// static int max_split_count = -1;
static int nextModuleInQ;
static int modulesInQ;
static pthread_mutex_t QLock = PTHREAD_MUTEX_INITIALIZER;
static int thread_count;
static pthread_t *t_id = NULL;
static int debug_flag = 0;
// static int d_flag = 0;
// static int d_flag2 = 0;
struct module_info_NCAT modules_NCAT = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0};
static int n_energies = 0;
static int n_materials = 0;
static int NUM_MODELS = 0;
static int NUM_NRB = 0;
static int NUM_POLY = 0;
static FILE *fpo = NULL;
SURFACE *nrb_model = NULL;
BEZIER_MODEL *bez_model = NULL;
TRI_MODEL *tri_model = NULL;
bvh_element **treepointer_nrb = NULL;                    /*Holds bounding volume hierarchies for organs, for nurbs*/
bvh_element **treepointer_tri = NULL;                    /*Holds bounding volume hierarchies for organs, for polygon*/
static float mu_table[MAXIMUM_NUMBER_OF_MATERIALS][300]; /*FIXME: harcoded number of materials and energy bins*/
static float **mu_table_tri = NULL;
static int use_tri_model = 0;
// static int s_stack[10000];
// static int s_depth = 0;

static int max_num_models = 5000;

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#define NR_END 1
#define FREE_ARG char *
void Subdivide_patch(double patch[4][4][3], double ul_patch[4][4][3], double ur_patch[4][4][3], double dl_patch[4][4][3], double dr_patch[4][4][3]);

static void dbug(int level, const char *s, ...)

{
  if (debug_flag >= level)
  {
    va_list ap;
    va_start(ap, s);
    vprintf(s, ap);
    va_end(ap);
    fflush(stdout);
  }
}

int comp_lines(const void *a, const void *b)

{
  const LINE_SEG da = *(const LINE_SEG *)a;
  const LINE_SEG db = *(const LINE_SEG *)b;

  // return (da.x2 > db.x2) - (da.x2 < db.x2);  // why did this fail?
  return (da.x2 + da.x1 > db.x2 + db.x1) - (da.x2 + da.x1 < db.x2 + db.x1);
}

int comp_intersections(const void *a, const void *b)

{
  const HIT da = *(const HIT *)a;
  const HIT db = *(const HIT *)b;

  return (da.x > db.x) - (da.x < db.x);
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

listelement *list_vector(long nl, long nh)
{
  listelement *v;
  long i;
  v = (listelement *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(listelement)));
  if (!v)
    nrerror("allocation failure in listelement vector()");
  return v - nl + NR_END;
}

void free_list_vector(listelement *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
  nh = nh;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  float ***t;

  t = (float ***)malloc((size_t)((nrow + NR_END) * sizeof(float **)));
  if (!t)
  {
    printf("allocation failure 1 in f3tensor()");
    exit(1);
  }
  t += NR_END;
  t -= nrl;

  t[nrl] = (float **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float *)));
  if (!t[nrl])
  {
    printf("allocation failure 2 in f3tensor()");
    exit(1);
  }
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  t[nrl][ncl] = (float *)malloc((size_t)((nrow * ncol * ndep + NR_END) *
                                         sizeof(float)));
  if (!t[nrl][ncl])
  {
    printf("allocation failure 3 in f3tensor()");
    exit(1);
  }
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  return t;
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
  nrh = nrh;
  nch = nch;
  ndh = ndh;
}

static float *vector(long nl, long nh) // JDP: made static
/* allocate & initialize a float vector with subscript range v[nl..nh] */
{
  float *v;
  long i;
  v = (float *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v)
    nrerror("allocation failure in vector()");
  for (i = nl; i <= nh; i++)
    v[i] = 0.0;
  return v - nl + NR_END;
}

static void free_vector(float *v, long nl, long nh) // JDP: made static
/* free a float vector allocated with vector() */
{
  free((FREE_ARG)(v + nl - NR_END));
  nh = nh;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl])
    nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
  nch = nch;
  nrh = nrh;
}

/* Allocation for structures */
TRIANGLE *tri_vector(long nl, long nh)
{
  TRIANGLE *v;

  v = (TRIANGLE *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(TRIANGLE)));
  if (!v)
  {
    printf("\nallocation failure in tri_vector()");
    exit(1);
  }

  return v - nl + NR_END;
}

void free_tri_vector(TRIANGLE *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
  nh = nh;
}

/*----------NURBS.C-------------*/
#define TINY 1.0e-20;
#define NR_END 1
#define FREE_ARG char *

BEZIER_PATCH *bp_vector(long nl, long nh)
{
  BEZIER_PATCH *v;

  v = (BEZIER_PATCH *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(BEZIER_PATCH)));
  if (!v)
  {
    printf("\nallocation error in bp_vector");
    exit(1);
  }
  return v - nl + NR_END;
}

void free_bpvector(BEZIER_PATCH *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

int *ivector(long nl, long nh)
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v)
  {
    printf("\nallocation error in ivector");
    exit(1);
  }
  return v - nl + NR_END;
}

void free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

CPOINT *cp_vector(long nl, long nh)
{
  CPOINT *v;

  v = (CPOINT *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(CPOINT)));
  if (!v)
  {
    printf("\nallocation error in cp_vector");
    exit(1);
  }
  return v - nl + NR_END;
}

POINT *p_vector(long nl, long nh)
{
  POINT *v;

  v = (POINT *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(POINT)));
  if (!v)
  {
    printf("\nallocation error in p_vector");
    exit(1);
  }
  return v - nl + NR_END;
}

void free_cpvector(CPOINT *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_pvector(POINT *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

POINT **p_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  POINT **m;

  m = (POINT **)malloc((size_t)((nrow + NR_END) * sizeof(POINT *)));

  if (!m)
  {
    printf("/n allocation error in p_matrix");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  m[nrl] = (POINT *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(POINT)));
  if (!m[nrl])
  {
    printf("/n allocation error in p_matrix");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  return m;
}

POINT ***p_3d(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  POINT ***t;

  t = (POINT ***)malloc((size_t)((nrow + NR_END) * sizeof(POINT **)));
  if (!t)
  {
    printf("allocation failure 1 in f3tensor()");
    exit(1);
  }
  t += NR_END;
  t -= nrl;

  t[nrl] = (POINT **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(POINT *)));
  if (!t[nrl])
  {
    printf("allocation failure 2 in f3tensor()");
    exit(1);
  }
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  t[nrl][ncl] = (POINT *)malloc((size_t)((nrow * ncol * ndep + NR_END) *
                                         sizeof(POINT)));
  if (!t[nrl][ncl])
  {
    printf("allocation failure 3 in f3tensor()");
    exit(1);
  }
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  return t;
}

void free_p_3d(POINT ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
  nrh = nrh;
  nch = nch;
  ndh = ndh;
}

CPOINT ***c_3d(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  CPOINT ***t;

  t = (CPOINT ***)malloc((size_t)((nrow + NR_END) * sizeof(CPOINT **)));
  if (!t)
  {
    printf("allocation failure 1 in f3tensor()");
    exit(1);
  }
  t += NR_END;
  t -= nrl;

  t[nrl] = (CPOINT **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(CPOINT *)));
  if (!t[nrl])
  {
    printf("allocation failure 2 in f3tensor()");
    exit(1);
  }
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  t[nrl][ncl] = (CPOINT *)malloc((size_t)((nrow * ncol * ndep + NR_END) *
                                          sizeof(CPOINT)));
  if (!t[nrl][ncl])
  {
    printf("allocation failure 3 in f3tensor()");
    exit(1);
  }
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  return t;
}

void free_c_3d(CPOINT ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
  nrh = nrh;
  nch = nch;
  ndh = ndh;
}

CURVE *curve_vector(long nl, long nh)
{
  CURVE *v;

  v = (CURVE *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(CURVE)));
  if (!v)
  {
    printf("\nallocation error in curve_vector");
    exit(1);
  }
  return v - nl + NR_END;
}

void free_curve_vector(CURVE *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
  nh = nh;
}

CURVE **curve_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  CURVE **m;

  m = (CURVE **)malloc((size_t)((nrow + NR_END) * sizeof(CURVE *)));

  if (!m)
  {
    printf("/n allocation error in curve_matrix");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  m[nrl] = (CURVE *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(CURVE)));
  if (!m[nrl])
  {
    printf("/n allocation error in curve_matrix");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  return m;
}

void free_curve_matrix(CURVE **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

CURVE ***curve_3d(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  CURVE ***t;

  t = (CURVE ***)malloc((size_t)((nrow + NR_END) * sizeof(CURVE **)));
  if (!t)
  {
    printf("allocation failure 1 in f3tensor()");
    exit(1);
  }
  t += NR_END;
  t -= nrl;

  t[nrl] = (CURVE **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(CURVE *)));
  if (!t[nrl])
  {
    printf("allocation failure 2 in f3tensor()");
    exit(1);
  }
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  t[nrl][ncl] = (CURVE *)malloc((size_t)((nrow * ncol * ndep + NR_END) *
                                         sizeof(CURVE)));
  if (!t[nrl][ncl])
  {
    printf("allocation failure 3 in f3tensor()");
    exit(1);
  }
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }
  return t;
}

void free_curve_3d(CURVE ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
  free((FREE_ARG)(t[nrl] + ncl - NR_END));
  free((FREE_ARG)(t + nrl - NR_END));
  nrh = nrh;
  nch = nch;
  ndh = ndh;
}

void free_pmatrix(POINT **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

CPOINT **cp_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  CPOINT **m;

  m = (CPOINT **)malloc((size_t)((nrow + NR_END) * sizeof(CPOINT *)));

  if (!m)
  {
    printf("/n allocation error in cp_matrix");
    exit(1);
  }
  m += NR_END;
  m -= nrl;

  m[nrl] = (CPOINT *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(CPOINT)));
  if (!m[nrl])
  {
    printf("/n allocation error in cp_matrix");
    exit(1);
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;
  return m;
}

void free_cpmatrix(CPOINT **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

void ludcmp(float **a, int n, int *indx)
{
  int i, imax, j, k;
  float big, dum, sum, temp;
  float *vv;

  vv = vector(1, n);
  for (i = 1; i <= n; i++)
  {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i][j])) > big)
        big = temp;
    if (big == 0.0)

    {
      printf("Singular matrix in routine LUDCMP");
      exit(1);
    }
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++)
  {
    for (i = 1; i < j; i++)
    {
      sum = a[i][j];
      for (k = 1; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++)
    {
      sum = a[i][j];
      for (k = 1; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 1; k <= n; k++)
      {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n)
    {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i <= n; i++)
        a[i][j] *= dum;
    }
  }
  free_vector(vv, 1, n);
}

#undef TINY

void lubksb(float **a, int n, int *indx, float b[])
{
  int i, ii = 0, ip, j;
  float sum;

  for (i = 1; i <= n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i][j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n; i >= 1; i--)
  {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

void MakeCurve(CURVE *C, INDEX n, INDEX m, DEGREE p)
{
  int i;

  C->pol.n = n;
  C->pol.Pw = cp_vector(0, n);
  for (i = 0; i <= n; i++)
  {
    C->pol.Pw[i].x = 0.0;
    C->pol.Pw[i].y = 0.0;
    C->pol.Pw[i].z = 0.0;
    C->pol.Pw[i].w = 0.0;
  }
  C->p = p;
  C->knt.m = m;
  C->knt.U = vector(0, m);
  for (i = 0; i <= m; i++)
    C->knt.U[i] = 0.0;

  C->uk = vector(0, n - 1);
  for (i = 0; i < n; i++)
    C->uk[i] = 0.0;
}

void BasisFuns(INDEX i, float u, DEGREE p, KNOTVECTOR knot_v, float *N)
{
  INDEX j, r;
  float *left, *right;
  float saved, temp;

  left = vector(1, p);
  right = vector(1, p);
  N[0] = 1.0;

  for (j = 1; j <= p; j++)
  {
    left[j] = u - knot_v.U[i + 1 - j];
    right[j] = knot_v.U[i + j] - u;
    saved = 0.0;
    for (r = 0; r < j; r++)
    {
      temp = N[r] / (right[r + 1] + left[j - r]);
      N[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    N[j] = saved;
  }
  free_vector(left, 1, p);
  free_vector(right, 1, p);
}

int FindSpan(INDEX n, DEGREE p, float u, KNOTVECTOR knot_v)
{
  int low, high, mid;

  if ((int)(u) == 1)
    return (knot_v.m - p - 1);

  low = 0;
  high = n + 1;
  mid = (low + high) / 2;

  while (u < knot_v.U[mid] || u >= knot_v.U[mid + 1])
  {
    if (u < knot_v.U[mid])
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2;
  }
  return (mid);
}

POINT CurvePoint(INDEX n, DEGREE p, KNOTVECTOR knot_v, CPOINT *Pw, float u)
{
  CPOINT Cw;
  POINT C;
  float *N;
  int i, j, span;

  N = vector(0, p);
  for (i = 0; i <= p; i++)
    N[i] = 0;

  span = FindSpan(n, p, u, knot_v);
  BasisFuns(span, u, p, knot_v, N);

  Cw.x = 0.0;
  Cw.y = 0.0;
  Cw.z = 0.0;
  Cw.w = 0.0;

  for (j = 0; j <= p; j++)
  {
    Cw.x = Cw.x + N[j] * Pw[span - p + j].x;
    Cw.y = Cw.y + N[j] * Pw[span - p + j].y;
    Cw.z = Cw.z + N[j] * Pw[span - p + j].z;
    Cw.w = Cw.w + N[j] * Pw[span - p + j].w;
  }

  if (Cw.w != 0)
  {
    C.x = Cw.x / Cw.w;
    C.y = Cw.y / Cw.w;
    C.z = Cw.z / Cw.w;
  }
  else
  {
    C.x = Cw.x;
    C.y = Cw.y;
    C.z = Cw.z;
  }

  free_vector(N, 0, p);
  return C;
}

POINT SurfacePoint(INDEX n, DEGREE p, KNOTVECTOR u_knot, INDEX m, DEGREE q, KNOTVECTOR v_knot, CPOINT **Pw, float u, float v)
{
  CPOINT Sw, *temp;
  POINT S;
  float *Nu, *Nv;
  int l, k, uspan, vspan;

  temp = cp_vector(0, q);

  Nu = vector(0, p);
  Nv = vector(0, q);

  for (l = 0; l <= p; l++)
    Nu[l] = 0;

  for (l = 0; l <= q; l++)
    Nv[l] = 0;

  uspan = FindSpan(n, p, u, u_knot);
  BasisFuns(uspan, u, p, u_knot, Nu);

  vspan = FindSpan(m, q, v, v_knot);
  BasisFuns(vspan, v, q, v_knot, Nv);

  for (l = 0; l <= q; l++)
  {
    temp[l].x = 0.0;
    temp[l].y = 0.0;
    temp[l].z = 0.0;
    temp[l].w = 0.0;

    for (k = 0; k <= p; k++)
    {
      temp[l].x = temp[l].x + Nu[k] * Pw[uspan - p + k][vspan - q + l].x;
      temp[l].y = temp[l].y + Nu[k] * Pw[uspan - p + k][vspan - q + l].y;
      temp[l].z = temp[l].z + Nu[k] * Pw[uspan - p + k][vspan - q + l].z;
      temp[l].w = temp[l].w + Nu[k] * Pw[uspan - p + k][vspan - q + l].w;
    }
  }

  Sw.x = 0.0;
  Sw.y = 0.0;
  Sw.z = 0.0;
  Sw.w = 0.0;

  for (l = 0; l <= q; l++)
  {
    Sw.x = Sw.x + Nv[l] * temp[l].x;
    Sw.y = Sw.y + Nv[l] * temp[l].y;
    Sw.z = Sw.z + Nv[l] * temp[l].z;
    Sw.w = Sw.w + Nv[l] * temp[l].w;
  }

  if (Sw.w != 0)
  {
    S.x = Sw.x / Sw.w;
    S.y = Sw.y / Sw.w;
    S.z = Sw.z / Sw.w;
  }
  else
  {
    S.x = Sw.x;
    S.y = Sw.y;
    S.z = Sw.z;
  }

  free_cpvector(temp, 0, q);
  free_vector(Nu, 0, p);
  free_vector(Nv, 0, q);

  return S;
}

float distance3d(POINT p1, POINT p2)
{
  float tempx, tempy, tempz;

  tempx = p1.x - p2.x;
  tempy = p1.y - p2.y;
  tempz = p1.z - p2.z;

  return (sqrt(tempx * tempx + tempy * tempy + tempz * tempz));
}

void Calc_UniformKnotVector(INDEX n, DEGREE p, KNOTVECTOR knot_v, float *uk)
{
  INDEX i, j, m;
  float d = 0.0;
  float sum = 0.0;

  m = n + p + 1;

  uk[0] = 0;
  uk[n] = 1;

  for (i = 1; i < n; i++)
    uk[i] = (float)i / (float)n;

  /* Calculate Knot Vector */
  for (i = 0; i <= p; i++)
    knot_v.U[i] = 0;
  for (i = m - p; i <= m; i++)
    knot_v.U[i] = 1;
  for (j = 1; j <= n - p; j++)
  {
    sum = 0.0;
    for (i = j; i <= j + p - 1; i++)
      sum += uk[i];
    knot_v.U[j + p] = sum / p;
  }

  knot_v.m = m;
}

void Calc_KnotVector(INDEX n, DEGREE p, POINT *Qw, KNOTVECTOR knot_v, float *uk)
{
  INDEX i, j, m;
  float d = 0.0;
  float sum = 0.0;

  m = n + p + 1;

  uk[0] = 0;
  uk[n] = 1;
  for (i = 1; i <= n; i++)
    d += distance3d(Qw[i], Qw[i - 1]);

  if (d == 0)
  {
    for (i = 1; i < n; i++)
      uk[i] = (float)i / (float)n;
  }

  else
  {
    for (i = 1; i < n; i++)
      uk[i] = uk[i - 1] + (distance3d(Qw[i], Qw[i - 1]) / d);
  }

  /* Calculate Knot Vector */
  for (i = 0; i <= p; i++)
    knot_v.U[i] = 0;
  for (i = m - p; i <= m; i++)
    knot_v.U[i] = 1;
  for (j = 1; j <= n - p; j++)
  {
    sum = 0.0;
    for (i = j; i <= j + p - 1; i++)
      sum += uk[i];
    knot_v.U[j + p] = sum / p;
  }

  knot_v.m = m;
}

void GlobalCurveInterp(INDEX n, POINT *Qw, int r, DEGREE p, KNOTVECTOR knot_v, float *uk, CPOINT *Pw)
{
  float **mA, **mA2;
  float *rhs, *ipointer;
  int *indx;
  float d = 0.0;
  float sum = 0.0;
  INDEX i, j;
  int span;

  mA = matrix(0, n, 0, n);
  mA2 = matrix(1, n + 1, 1, n + 1);
  indx = ivector(1, n + 1);
  rhs = vector(1, n + 1);

  /*Initialize mA to zero*/
  for (i = 0; i <= n; i++)
    for (j = 0; j <= n; j++)
      mA[i][j] = 0.0;

  for (i = 0; i <= n; i++)
  {
    span = FindSpan(n, p, uk[i], knot_v);
    ipointer = &mA[i][span - p];
    BasisFuns(span, uk[i], p, knot_v, ipointer);
  }

  for (i = 0; i <= n; i++)
  {
    for (j = 0; j <= n; j++)
    {
      mA2[i + 1][j + 1] = mA[i][j];
    }
  }

  ludcmp(mA2, n + 1, indx);

  for (i = 0; i < r; i++)
  {
    for (j = 0; j <= n; j++)
    {
      if (i == 0)
        rhs[j + 1] = Qw[j].x;
      else if (i == 1)
        rhs[j + 1] = Qw[j].y;
      else if (i == 2)
        rhs[j + 1] = Qw[j].z;
    }
    lubksb(mA2, n + 1, indx, rhs);
    for (j = 0; j <= n; j++)
    {
      if (i == 0)
        Pw[j].x = rhs[j + 1];
      if (i == 1)
        Pw[j].y = rhs[j + 1];
      if (i == 2)
        Pw[j].z = rhs[j + 1];
    }
  }

  free_ivector(indx, 1, n + 1);
  free_vector(rhs, 1, n + 1);
  free_matrix(mA, 0, n, 0, n);
  free_matrix(mA2, 1, n + 1, 1, n + 1);
}

void SurfMeshParams(INDEX n, INDEX m, QNET Q, float *uk, float *vl)
{
  INDEX num, k, l, max;
  float total, *cds, d;

  if (m >= n)
    max = m;
  else
    max = n;

  cds = vector(1, max + 1);
  num = m + 1;

  uk[0] = 0.0;
  uk[n] = 1.0;

  for (k = 1; k < n; k++)
    uk[k] = 0.0;

  for (l = 0; l <= m; l++)
  {
    total = 0.0;

    for (k = 1; k <= n; k++)
    {
      cds[k] = distance3d(Q.Qw[k][l], Q.Qw[k - 1][l]);
      total = total + cds[k];
    }

    if (total == 0.0)
      num = num - 1;
    else
    {
      d = 0.0;
      for (k = 1; k < n; k++)
      {
        d = d + cds[k];
        uk[k] = uk[k] + d / total;
      }
    }
  }

  if (num == 0)
  {
    printf("\n Error in SurfMeshParams!!!");
    exit(1);
  }

  for (k = 0; k < n; k++)
    uk[k] = uk[k] / num;

  num = n + 1;

  vl[0] = 0.0;
  vl[m] = 1.0;
  for (l = 1; l < m; l++)
    vl[l] = 0.0;
  for (k = 0; k <= n; k++)
  {
    total = 0.0;
    for (l = 1; l <= m; l++)
    {
      cds[l] = distance3d(Q.Qw[k][l], Q.Qw[k][l - 1]);
      total = total + cds[l];
    }
    if (total == 0.0)
      num = num - 1;
    else
    {
      d = 0.0;
      for (l = 1; l < m; l++)
      {
        d = d + cds[l];
        vl[l] = vl[l] + d / total;
      }
    }
  }

  if (num == 0)
  {
    printf("\n Error in SurfMeshParams!!!");
    exit(1);
  }

  for (l = 0; l < m; l++)
    vl[l] = vl[l] / num;

  free_vector(cds, 1, max + 1);
}

void GlobalSurfInterp(INDEX n, INDEX m, QNET Q, DEGREE p, DEGREE q, KNOTVECTOR u_knot, KNOTVECTOR v_knot, CNET P)
{
  INDEX i, j, l, k, max, uknots, vknots;
  float sum;
  QPOINTS QP;
  CNET R;
  float *uk, *vl;
  CPOINT **temp, *Rpw;

  temp = cp_matrix(0, n, 0, m);
  uk = vector(0, n + 1);
  vl = vector(0, m + 1);

  if (m >= n)
    max = m;
  else
    max = n;

  Rpw = cp_vector(0, max);
  QP.Qw = p_vector(0, max);

  for (i = 0; i <= n; i++)
    for (j = 0; j <= m; j++)
    {
      temp[i][j].x = 0;
      temp[i][j].y = 0;
      temp[i][j].z = 0;
      temp[i][j].w = 0;
    }

  R.Pw = temp;
  SurfMeshParams(n, m, Q, uk, vl);
  uknots = n + p + 1;
  vknots = m + q + 1;

  /* Calculate Knot Vector U */
  for (i = 0; i <= p; i++)
    u_knot.U[i] = 0;
  for (i = uknots - p; i <= uknots; i++)
    u_knot.U[i] = 1;
  for (j = 1; j <= n - p; j++)
  {
    sum = 0.0;
    for (i = j; i <= j + p - 1; i++)
      sum += uk[i];
    u_knot.U[j + p] = sum / p;
  }

  /* Calculate Knot Vector V */
  for (i = 0; i <= q; i++)
    v_knot.U[i] = 0;
  for (i = vknots - q; i <= vknots; i++)
    v_knot.U[i] = 1;
  for (j = 1; j <= m - q; j++)
  {
    sum = 0.0;
    for (i = j; i <= j + q - 1; i++)
      sum += vl[i];
    v_knot.U[j + q] = sum / q;
  }

  for (l = 0; l <= m; l++)
  {
    for (k = 0; k <= n; k++)
    {
      QP.Qw[k].x = Q.Qw[k][l].x;
      QP.Qw[k].y = Q.Qw[k][l].y;
      QP.Qw[k].z = Q.Qw[k][l].z;
    }
    GlobalCurveInterp(n, QP.Qw, 3, p, u_knot, uk, Rpw);
    for (i = 0; i <= n; i++)
    {
      R.Pw[i][l].x = Rpw[i].x;
      R.Pw[i][l].y = Rpw[i].y;
      R.Pw[i][l].z = Rpw[i].z;
      R.Pw[i][l].w = Rpw[i].w;
    }
  }

  for (i = 0; i <= n; i++)
  {
    for (l = 0; l <= m; l++)
    {
      QP.Qw[l].x = R.Pw[i][l].x;
      QP.Qw[l].y = R.Pw[i][l].y;
      QP.Qw[l].z = R.Pw[i][l].z;
    }
    GlobalCurveInterp(m, QP.Qw, 3, q, v_knot, vl, Rpw);
    for (l = 0; l <= m; l++)
    {
      P.Pw[i][l].x = Rpw[l].x;
      P.Pw[i][l].y = Rpw[l].y;
      P.Pw[i][l].z = Rpw[l].z;
      P.Pw[i][l].w = Rpw[l].w;
    }
  }

  free_cpmatrix(temp, 0, n, 0, m);

  free_cpvector(Rpw, 0, max);
  free_pvector(QP.Qw, 0, max);
}

/*-----------------------------------Bezier Patch routines-------------------------------------------*/
/* The following routines are used in converting a NURBS surface into Bezier patches                 */
/*---------------------------------------------------------------------------------------------------*/
int get_breakpoint(int len_kU, float *kU, float w)
{
  register int i;
  i = 0;
  while ((i < len_kU) && (kU[i] <= w))
    i++;
  return (i - 1);
}

static struct knotmultCnt_s
{
  struct knotmult_s
  {
    int mu;    /* multiplicity */
    float val; /* value */
  } *umult, *vmult;
  int numU, numV; /* # distinct knots */
} knots;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
void refine_patch(patch *src, patch *dest)
/*  This routine refines the Bezier patches created for a NURBS surface */
{
  register int i, j, s, i2;
  int k, last, r, li, delta;
  float *kU, *kW, *kV;
  CPOINT C[4], *D;
  float lenU, lenV, omega;

  /* Refine along U-direction first */
  if (dest->numU > src->numU) /* there are new breakpoints */
  {
    k = src->ordU;
    kU = src->kU;
    kW = dest->kU;
    lenU = src->numU + src->ordU; /* #knots in U-direction */

    for (i = 0; i < src->numV; i++) /* for each V-vertex row do */
      for (j = 0; j < dest->numU; j++)
      { /* for each w[j], calculate the new */
        /* vertex  */
        delta = get_breakpoint(lenU, kU, kW[j]);

        for (s = 0; s <= MIN(k - 1, delta); s++)
          C[s] = src->points[i][delta - s]; /* Initialize C-array */

        /* Calculate the new vertex W[j] for this value of kW[j] */
        for (r = k; r > 1; r--)
        {
          li = delta;
          last = MIN(r - 2, delta);
          for (s = 0; s < last; s++)
          {
            omega = (kW[j + r - 1] - kU[li]) / (kU[li + r - 1] - kU[li]);
            C[s].x = omega * C[s].x + (1 - omega) * C[s + 1].x;
            C[s].y = omega * C[s].y + (1 - omega) * C[s + 1].y;
            C[s].z = omega * C[s].z + (1 - omega) * C[s + 1].z;
            C[s].w = omega * C[s].w + (1 - omega) * C[s + 1].w;
            li--;
          }
          omega = (kW[j + r - 1] - kU[li]) / (kU[li + r - 1] - kU[li]);
          if (last < (r - 2))
          {
            C[s].x = omega * C[s].x;
            C[s].y = omega * C[s].y;
            C[s].z = omega * C[s].z;
            C[s].w = omega * C[s].w;
          }
          else
          {
            C[s].x = omega * C[s].x + (1 - omega) * C[s + 1].x;
            C[s].y = omega * C[s].y + (1 - omega) * C[s + 1].y;
            C[s].z = omega * C[s].z + (1 - omega) * C[s + 1].z;
            C[s].w = omega * C[s].w + (1 - omega) * C[s + 1].w;
          }
        }
        /* Now, the value of the new vertex is available in C[0] */
        dest->points[i][j] = C[0];
      }
  }
  else /* no new breakpoints in U-direction */
    for (i = 0; i < src->numV; i++)
      memcpy(dest->points[i], src->points[i], src->numU * sizeof(CPOINT));

  /* Now perform refinement along the V-direction */
  if (dest->numV > src->numV) /* there are new breakpoints */
  {
    k = src->ordV;
    kV = src->kV;
    kW = dest->kV;
    lenV = src->numV + src->ordV; /* #knots in V-direction */

    D = (CPOINT *)malloc(src->numV * (sizeof(CPOINT)));

    for (i = 0; i < dest->numU; i++) /* for each U-vertex col do */
    {
      for (i2 = 0; i2 < src->numV; i2++) /* Make a copy of the */
        D[i2] = dest->points[i2][i];     /* i'th column */

      for (j = 0; j < dest->numV; j++)
      { /* for each w[j], calculate the new */
        /* vertex  */
        delta = get_breakpoint(lenV, kV, kW[j]);

        for (s = 0; s <= MIN(k - 1, delta); s++)
          C[s] = D[delta - s]; /* Initialize C-array */

        /* Calculate the new vertex W[j] for this value of kW[j] */
        for (r = k; r > 1; r--)
        {
          li = delta;
          last = MIN(r - 2, delta);
          for (s = 0; s < last; s++)
          {
            omega = (kW[j + r - 1] - kV[li]) / (kV[li + r - 1] - kV[li]);
            C[s].x = omega * C[s].x + (1 - omega) * C[s + 1].x;
            C[s].y = omega * C[s].y + (1 - omega) * C[s + 1].y;
            C[s].z = omega * C[s].z + (1 - omega) * C[s + 1].z;
            C[s].w = omega * C[s].w + (1 - omega) * C[s + 1].w;
            li--;
          }
          omega = (kW[j + r - 1] - kV[li]) / (kV[li + r - 1] - kV[li]);
          if (last < (r - 2))
          {
            C[s].x = omega * C[s].x;
            C[s].y = omega * C[s].y;
            C[s].z = omega * C[s].z;
            C[s].w = omega * C[s].w;
          }
          else
          {
            C[s].x = omega * C[s].x + (1 - omega) * C[s + 1].x;
            C[s].y = omega * C[s].y + (1 - omega) * C[s + 1].y;
            C[s].z = omega * C[s].z + (1 - omega) * C[s + 1].z;
            C[s].w = omega * C[s].w + (1 - omega) * C[s + 1].w;
          }
        }
        /* Now, the value of the new vertex is available in C[0] */
        dest->points[j][i] = C[0];
      }
    }
    free(D);
  }
  free(knots.umult);
  free(knots.vmult);
}

#define KNOTEPS 1e-15
#define EPSEQ(a, b) (fabs((float)(a - b)) < KNOTEPS)
int get_knot_multiplicities(patch *src)
/* This function checks the multiplicity for each knot of a NURBS surface */
/* To convert a NURBS surface into Bezier patches, each knot must have a multiplicity of 4*/
{
  register int i;
  int curr, lenU, lenV, ncurr;
  float u, v;

  i = curr = 0;
  lenU = src->ordU + src->numU;
  knots.umult = (struct knotmult_s *)malloc(lenU * sizeof(struct knotmult_s));

  while (i < lenU) /* traverse all the U-knots */
  {
    u = src->kU[i];
    for (ncurr = 0; (i < lenU) && EPSEQ(src->kU[i], u); i++, ncurr++)
      ;
    knots.umult[curr].val = u;
    if (ncurr > src->ordU)
      ncurr = src->ordU;
    knots.umult[curr].mu = ncurr;
    curr++;
  }
  knots.numU = curr - 1;

  i = 0;
  curr = 0;
  lenV = src->ordV + src->numV;
  knots.vmult = (struct knotmult_s *)malloc(lenV * sizeof(struct knotmult_s));

  while (i < lenV) /* traverse all the V-knots */
  {
    v = src->kV[i];
    for (ncurr = 0; (i < lenV) && EPSEQ(src->kV[i], v); i++, ncurr++)
      ;
    knots.vmult[curr].val = v;
    if (ncurr > src->ordV)
      ncurr = src->ordV;
    knots.vmult[curr].mu = ncurr;
    curr++;
  }
  knots.numV = curr - 1;
  return 0;
}

int alloc_patch(patch *p)
/* Allocate space for storing knot arrays and the vertices */
/* constituting the patch */
{
  int i;

  if (((p->kU = (float *)malloc((p->numU + p->ordU) * sizeof(float))) == NULL) ||
      ((p->kV = (float *)malloc((p->numV + p->ordV) * sizeof(float))) == NULL) ||
      ((p->points = (CPOINT **)malloc(p->numV * sizeof(CPOINT *))) ==
       NULL))
  {
    perror("malloc");
    exit(1);
  }

  for (i = 0; i < p->numV; i++)
    if ((p->points[i] = (CPOINT *)malloc(p->numU * sizeof(CPOINT))) == NULL)
    {
      perror("malloc");
      exit(1);
    }
  return 0;
}

int free_patch(patch *p)
/* Free the space allocated for this patch since it is no longer */
/* required */
{
  int i;

  free(p->kU);
  free(p->kV);
  for (i = 0; i < p->numV; i++)
    free(p->points[i]);

  free(p->points);
  return 0;
}

#define MAXU 1000
#define MAXV 1000
#define MAXT 1000
int insert_multiple_knots(patch *src, patch *dest)
/* Function inserts knots into a NURBS surface: Knots are inserted into a NURBS surface in order to */
/* convert it into Bezier patches */
{
  register int i, j;
  int ucurr, vcurr, kcurr;
  float tmpkU[MAXU], tmpkV[MAXV];

  assert(src->ordU <= MAX_ORDER);
  assert(src->ordV <= MAX_ORDER);

  *dest = *src;
  get_knot_multiplicities(src);

  dest->kU = tmpkU;
  dest->kV = tmpkV;

  kcurr = ucurr = i = 0;

  while (i + knots.umult[kcurr].mu < dest->ordU)
  {
    for (j = 0; (j < knots.umult[kcurr].mu); j++, i++)
      dest->kU[ucurr++] = knots.umult[kcurr].val;
    kcurr++;
  }
  while (i < src->numU + 1)
  {
    for (j = 0; j < dest->ordU; j++)
      dest->kU[ucurr++] = knots.umult[kcurr].val;

    i += knots.umult[kcurr].mu;
    kcurr++;
  }
  while (i < src->numU + src->ordU)
    dest->kU[ucurr++] = src->kU[i++];

  /* Repeat the same process for V */
  kcurr = vcurr = i = 0;

  while (i + knots.vmult[kcurr].mu < dest->ordV)
  {
    for (j = 0; (j < knots.vmult[kcurr].mu); j++, i++)
      dest->kV[vcurr++] = knots.vmult[kcurr].val;
    kcurr++;
  }
  while (i < src->numV + 1)
  {
    for (j = 0; j < dest->ordV; j++)
      dest->kV[vcurr++] = knots.vmult[kcurr].val;

    i += knots.vmult[kcurr].mu;
    kcurr++;
  }
  while (i < src->numV + src->ordV)
    dest->kV[vcurr++] = src->kV[i++];

  /* Now allocate space for new patch to hold the knot arrays and the */
  /* vertices */
  dest->numU = ucurr - dest->ordU;
  dest->numV = vcurr - dest->ordV;
  alloc_patch(dest);
  memcpy(dest->kU, tmpkU, ucurr * sizeof(float));
  memcpy(dest->kV, tmpkV, vcurr * sizeof(float));
  return 0;
}

int setup_initial_patch(patch *p, SURFACE *nrb_model)
/* This routine sets up the initial Bezier patch (p) for the NURBS surface (nrb_model) */
{
  register int i, j;
  int lenU, lenV;
  float *tmpkU;

  p->ordU = 3;
  p->ordV = 3;
  p->ordU++;
  p->ordV++;

  lenU = nrb_model->net.n + 4;
  p->numU = lenU - p->ordU;
  assert((tmpkU = (float *)malloc(lenU * sizeof(float))) != NULL);
  for (i = 0; i < lenU; i++) /* Read in the U-knots */
    tmpkU[i] = nrb_model->knu.U[i];

  lenV = nrb_model->net.m + 4;
  p->numV = lenV - p->ordV;
  alloc_patch(p); /* Allocate space for the tables */
  memcpy(p->kU, tmpkU, lenU * sizeof(float));
  free(tmpkU);
  for (i = 0; i < lenV; i++) /* Read in the V-knots */
    p->kV[i] = nrb_model->knv.U[i];

  lenU = nrb_model->net.n;
  lenV = nrb_model->net.m;

  assert(lenU == p->numU && lenV == p->numV);

  for (i = 0; i < p->numU; i++) /* read rational vertices */
    for (j = 0; j < p->numV; j++)
    {
      p->points[j][i].x = nrb_model->net.Pw[i][j].x;
      p->points[j][i].y = nrb_model->net.Pw[i][j].y;
      p->points[j][i].z = nrb_model->net.Pw[i][j].z;
      p->points[j][i].w = 1;
    }
  return 1;
}

void apply_rotation(double xform[4][3], double vin[3], double vout[3])

{
  int i;
  for (i = 0; i < 3; i++)
    vout[i] = DOT(vin, xform[i]);
}

void apply_rotation_f(double xform[4][3], float vin[3], float vout[3])

{
  int i;
  for (i = 0; i < 3; i++)
    vout[i] = DOT(vin, xform[i]);
}

void apply_xform_f(double xform[4][3], float vin[3], float vout[3])

{
  int i;
  for (i = 0; i < 3; i++)
    vout[i] = DOT(vin, xform[i]) + xform[3][i];
}

void apply_xform(double xform[4][3], double vin[3], double vout[3])

{
  int i;
  for (i = 0; i < 3; i++)
    vout[i] = DOT(vin, xform[i]) + xform[3][i];
}

void invert_xform(double xform[4][3], double xform_inv[4][3])

{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      xform_inv[i][j] = xform[j][i];
  apply_rotation(xform_inv, xform[3], xform_inv[3]);
  for (i = 0; i < 3; i++)
    xform_inv[3][i] = -xform_inv[3][i];
}

void cross_product(double *u, double *v, double *result)
{
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = -(u[0] * v[2] - u[2] * v[0]);
  result[2] = u[0] * v[1] - u[1] * v[0];
}

void get_patch_xform(double patch[4][4][3], double xform[4][3], double patch_xf[4][4][3])

{
  double axes[3][3], nrm[3], dot, center[3];
  int i, j, k;
  // First, we determine which axis on the patch is bigger... this becomes the main surface axis
  nrm[0] = 0;
  nrm[1] = 0;
  for (i = 0; i < 3; i++)
  {
    axes[0][i] = patch[0][0][i] + patch[0][3][i] - patch[3][0][i] - patch[3][3][i];
    nrm[0] += axes[0][i] * axes[0][i];
    axes[1][i] = patch[0][0][i] + patch[3][0][i] - patch[0][3][i] - patch[3][3][i];
    nrm[1] += axes[1][i] * axes[1][i];
  }
  if (nrm[1] > nrm[0])
  {
    // we swap them so the bigger one is first if needed
    for (i = 0; i < 3; i++)
    {
      axes[2][i] = axes[1][i];
      axes[1][i] = axes[0][i];
      axes[0][i] = axes[2][i];
    }
    nrm[0] = sqrt(nrm[1]);
  }
  else
    nrm[0] = sqrt(nrm[0]);
  // We ensure the other axis is orthogonal to the main one
  for (i = 0; i < 3; i++)
    axes[0][i] /= nrm[0]; // normalize main
  for (i = 0; i < 3; i++)
    dbug(3, "%1.10f --------------------*******************************************************************\n", axes[0][i]);
  // for(i=0;i<3;i++) printf("%1.15f\n",axes[0][i]);
  dot = DOT(axes[0], axes[1]);
  for (i = 0; i < 3; i++)
    axes[1][i] -= dot * axes[0][i]; // orthogonalize
  nrm[1] = 0;
  for (i = 0; i < 3; i++)
    nrm[1] += axes[1][i] * axes[1][i];
  nrm[1] = sqrt(nrm[1]);
  for (i = 0; i < 3; i++)
    axes[1][i] /= nrm[1]; // normalize 2nd
  // for(i=0;i<3;i++) printf("%1.15f\n",axes[1][i]);
  // We then compute the cross product to give the third axis
  cross_product(axes[0], axes[1], axes[2]);
  // We form the four key vectors of a 4x4 transform matrix, which includes rotation (3 vecs) and translation (1 vec)
  for (i = 0; i < 3; i++)
    center[i] = 0;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 3; k++)
        center[k] += patch[i][j][k];
  for (i = 0; i < 3; i++)
    center[i] /= 16;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
      xform[i][j] = axes[i][j];
    xform[3][i] = -DOT(axes[i], center);
  }
  // We transform all points into the new coordinate frame
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      apply_xform(xform, patch[i][j], patch_xf[i][j]);
}

int create_bezier_patches(patch *p, BEZIER_PATCH *patches)
/* Routine decomposes a NURBS surface (initial patch p) into many Bezier patches (patches) */
/* This is done by inserting knots into the NURBS surface until each knot has a multiplicity of 4 */
{
  register int i, j, k, l, l1, k1;
  int umark, vmark;
  int numpt1, numpt2;
  float uval, vval;
  extern statistics_t pstat;
  float minz, maxz, miny, maxy, minx, maxx;

  pstat.nbezs = 0;

  uval = p->kU[p->ordU - 1];
  umark = p->ordU - 1;
  while ((umark > 0) && (p->kU[umark] == uval))
    umark--;
  if (p->kU[umark] < uval)
    umark++;

  vval = p->kV[p->ordV - 1];
  vmark = p->ordV - 1;
  while ((vmark > 0) && (p->kV[vmark] == vval))
    vmark--;
  if (p->kV[vmark] < vval)
    vmark++;

  numpt1 = (p->numU - umark) / p->ordU;
  numpt2 = (p->numV - vmark) / p->ordV;

  for (l1 = 0, l = umark; l <= (p->numU - p->ordU); l += p->ordU, l1++)
  {
    for (k1 = 0, k = vmark; k <= (p->numV - p->ordV); k += p->ordV, k1++)
    {
      minz = p->points[k][l].z;
      maxz = minz;
      miny = p->points[k][l].y;
      maxy = miny;
      minx = p->points[k][l].x;
      maxx = minx;
      for (j = 0; j < p->ordU; j++)
        for (i = 0; i < p->ordV; i++)
        {
          patches[pstat.nbezs].cpoints[j][i][0] = p->points[k + i][l + j].x;
          patches[pstat.nbezs].cpoints[j][i][1] = p->points[k + i][l + j].y;
          patches[pstat.nbezs].cpoints[j][i][2] = p->points[k + i][l + j].z;

          if (p->points[k + i][l + j].z > maxz)
            maxz = p->points[k + i][l + j].z;
          if (p->points[k + i][l + j].z < minz)
            minz = p->points[k + i][l + j].z;
          if (p->points[k + i][l + j].y > maxy)
            maxy = p->points[k + i][l + j].y;
          if (p->points[k + i][l + j].y < miny)
            miny = p->points[k + i][l + j].y;
          if (p->points[k + i][l + j].x > maxx)
            maxx = p->points[k + i][l + j].x;
          if (p->points[k + i][l + j].x < minx)
            minx = p->points[k + i][l + j].x;
        }

      numpt1 = 0;
      numpt2 = 0;
      patches[pstat.nbezs].maxx = maxx;
      patches[pstat.nbezs].minx = minx;
      patches[pstat.nbezs].maxy = maxy;
      patches[pstat.nbezs].miny = miny;
      patches[pstat.nbezs].maxz = maxz;
      patches[pstat.nbezs].minz = minz;

      get_patch_xform(patches[pstat.nbezs].cpoints, patches[pstat.nbezs].xform, patches[pstat.nbezs].xpoints);

      pstat.nbezs++;
      pstat.tot_trimpts1 += numpt1;
      pstat.tot_trimpts2 += numpt2;
    }
  }
  return 0;
}

statistics_t pstat;

void SETUP_BEZIER_MODEL(SURFACE surf, BEZIER_MODEL *bez_model)
{
  bez_model->num_patches = (surf.net.n - 3) * (surf.net.m - 3);
  bez_model->patches = bp_vector(0, bez_model->num_patches);
}

void Plane_eqn(double *p1, double *p2, double *p3, double *p4, double *A, double *B, double *C, double *D)
{

  // In sept 2014, this function was updated so that it provides a better approximation to the surface
  // We deal with the problem of 0 area triangles when one side of the patch goes to length=0 by combining normal vectors computed from two triangular patches that span the 4 corners
  double v1[3], v2[3];
  double result[2][3], nrm[2], sum_nrm;
  int i;

  v1[0] = p2[0] - p1[0];
  v1[1] = p2[1] - p1[1];
  v1[2] = p2[2] - p1[2];

  v2[0] = p3[0] - p1[0];
  v2[1] = p3[1] - p1[1];
  v2[2] = p3[2] - p1[2];

  cross_product(v1, v2, result[0]);
  nrm[0] = DOT(result[0], result[0]);

  v1[0] = p3[0] - p4[0];
  v1[1] = p3[1] - p4[1];
  v1[2] = p3[2] - p4[2];

  v2[0] = p2[0] - p4[0];
  v2[1] = p2[1] - p4[1];
  v2[2] = p2[2] - p4[2];

  cross_product(v1, v2, result[1]);
  nrm[1] = DOT(result[1], result[1]);

  sum_nrm = nrm[0] + nrm[1];

  if(sum_nrm == 0.0) {
	  *A = 0; *B = 0; *C = 0; *D = 0;
	  // Fixed the issue, Mingye Wu, 7/15/2024
	  //dbug(-1,"*************************************XCAT*************************ERROR************************************");
	  return;
  }
  
  *A = (nrm[0] * result[0][0] + nrm[1] * result[1][0]) / sum_nrm;
  *B = (nrm[0] * result[0][1] + nrm[1] * result[1][1]) / sum_nrm;
  *C = (nrm[0] * result[0][2] + nrm[1] * result[1][2]) / sum_nrm;

  *D = -(*A) * p1[0] - (*B) * p1[1] - (*C) * p1[2];
}

void random_unit_vector(double *vec)

// note: the result is not uniformly distributed in S2

{
  double tmp;
  vec[0] = ((double)(rand() % 1000)) - 499.5;
  vec[1] = ((double)(rand() % 1000)) - 499.5;
  vec[2] = ((double)(rand() % 1000)) - 499.5;
  tmp = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  vec[0] /= tmp;
  vec[1] /= tmp;
  vec[2] /= tmp;
}

void random_unit_vector_f(float *vec)

// note: the result is not uniformly distributed in S2

{
  float tmp;
  vec[0] = ((float)(rand() % 1000)) - 499.5;
  vec[1] = ((float)(rand() % 1000)) - 499.5;
  vec[2] = ((float)(rand() % 1000)) - 499.5;
  tmp = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  vec[0] /= tmp;
  vec[1] /= tmp;
  vec[2] /= tmp;
}

int get_nrb_orientation(SURFACE *nrb_model)

//  THIS FUNCTION IS UNFINISHED AND UNUSED FOR NOW... SOME NRB SURFACES FOLD ON THEMSELVES, SO THE ORIENTATION IS NOT RELEVANT

{
  int i, j, ii, jj, flag = 1, cnt = 0;
  double vec[3], dot, tmp, out_point[3], mxd, mnd;
  double loop_points[5][3];
  double A, B, C, D;
  const int loop_ind1[] = {0, -1, 0, 1, 0};
  const int loop_ind2[] = {0, 0, -1, 0, 1};
  while (flag)
  {
    ii = 0;
    jj = 0;
    while ((ii == 0) || (jj == 0) || (ii == nrb_model->net.n - 1) || (jj == nrb_model->net.m - 1))
    {
      random_unit_vector(vec);
      dot = -10000;
      for (i = 0; i < nrb_model->net.n; i++)
        for (j = 0; j < nrb_model->net.m; j++)
        {
          tmp =
              vec[0] * nrb_model->net.Pw[i][j].x +
              vec[1] * nrb_model->net.Pw[i][j].y +
              vec[2] * nrb_model->net.Pw[i][j].z;
          if (tmp > dot)
          {
            dot = tmp;
            ii = i;
            jj = j;
          }
        }
    }
    // here is a point that is outside the surface
    out_point[0] = nrb_model->net.Pw[ii][jj].x + vec[0];
    out_point[1] = nrb_model->net.Pw[ii][jj].y + vec[1];
    out_point[2] = nrb_model->net.Pw[ii][jj].z + vec[2];

    for (i = 0; i < 5; i++)
    {
      loop_points[i][0] = nrb_model->net.Pw[ii + loop_ind1[i]][jj + loop_ind2[i]].x;
      loop_points[i][1] = nrb_model->net.Pw[ii + loop_ind1[i]][jj + loop_ind2[i]].y;
      loop_points[i][2] = nrb_model->net.Pw[ii + loop_ind1[i]][jj + loop_ind2[i]].z;
    }

    mxd = -100000;
    mnd = 100000;
    for (i = 0; i < 4; i++)
    {
      Plane_eqn(loop_points[0], loop_points[0], loop_points[1 + i], loop_points[(1 + i) % 4 + 1], &A, &B, &C, &D);
      dot = A * out_point[0] + B * out_point[1] + C * out_point[2] + D;
      if (dot > mxd)
        mxd = dot;
      if (dot < mnd)
        mnd = dot;
    }
    flag = (mxd * mnd <= 0);
    if (flag)
    {
      dbug(-1, "\n\n pts=[");
      for (i = 0; i < 5; i++)
        dbug(-1, "%1.14f %1.14f %1.14f;\r\n", loop_points[i][0], loop_points[i][1], loop_points[i][2]);
      dbug(-1, "];\nop=[%1.12f %1.12f %1.12f];\nclf;hold off;plot3([pts(1,1) op(1)], [pts(1,2) op(2)], [pts(1,3) op(3)],'r');\n",
           out_point[0], out_point[1], out_point[2]);
      dbug(-1, "hold on;ind=[0 1 2 0 2 3 0 3 4 0 4 1 0]+1;plot3(pts(ind,1), pts(ind,2), pts(ind,3));\n");
      dbug(-1, "plot3(op(1),op(2),op(3),'k*');\n");
      cnt++;
    }
    if (cnt > 20)
      exit(1);
  }
  return flag;
}

void write_nrb(patch *p)

{
  int i, j;
  fpo = fopen("output", "w");
  fprintf(fpo, "%d %d\n", p->numV, p->numU);
  for (i = 0; i < p->numV; i++)
    for (j = 0; j < p->numU; j++)
      fprintf(fpo, "%1.12f %1.12f %1.12f\n", p->points[i][j].x, p->points[i][j].y, p->points[i][j].z);
  fclose(fpo);
  exit(1);
}

void SPLINE2BEZ(SURFACE *nrb_model, BEZIER_MODEL *bez_model)
/* This is the main routine for converting a NURBS surface (nrb_model) into a Bezier model (bez_model) */
{
  patch src, dest;
  // int right_hand_out;

  pstat.count = 0;
  pstat.tot_nbezs = 0;
  pstat.tot_trimpts1 = 0;
  pstat.tot_trimpts2 = 0;

  // The original plan here was to ensure that the surface orientation of all nrb surfaces matched so that I could tell whether
  // any given ray was going into or out of the object by looking at the ray-plane orientation, but it turns out that some nurb
  // surfaces in XCAT fold on themselves (e.g., the inner surface of rrib8), so this technique will not work.
  //  right_hand_out = get_nrb_orientation(nrb_model);

  setup_initial_patch(&src, nrb_model);
  insert_multiple_knots(&src, &dest);
  refine_patch(&src, &dest);
  // if (debug_flag == 1) write_nrb(&dest);
  create_bezier_patches(&dest, bez_model->patches);

  free_patch(&dest);
  free_patch(&src);
}

int get_refinement_depth(BEZIER_MODEL *bez_model)

{
  return 1;
}

void tri_bbox(TRIANGLE *tri)

{
  tri->minx = MIN(tri->vertex[0].x, MIN(tri->vertex[1].x, tri->vertex[2].x));
  tri->maxx = MAX(tri->vertex[0].x, MAX(tri->vertex[1].x, tri->vertex[2].x));
  tri->miny = MIN(tri->vertex[0].y, MIN(tri->vertex[1].y, tri->vertex[2].y));
  tri->maxy = MAX(tri->vertex[0].y, MAX(tri->vertex[1].y, tri->vertex[2].y));
  tri->minz = MIN(tri->vertex[0].z, MIN(tri->vertex[1].z, tri->vertex[2].z));
  tri->maxz = MAX(tri->vertex[0].z, MAX(tri->vertex[1].z, tri->vertex[2].z));
}

void Calc_extents_tri(TRI_MODEL *tri_model)

{
  int i, j;
  float tx, ty, tz;

  /*Initialize*/
  tri_model->max_x = tri_model->tris[0].vertex[0].x;
  tri_model->max_y = tri_model->tris[0].vertex[0].y;
  tri_model->max_z = tri_model->tris[0].vertex[0].z;

  tri_model->min_x = tri_model->tris[0].vertex[0].x;
  tri_model->min_y = tri_model->tris[0].vertex[0].y;
  tri_model->min_z = tri_model->tris[0].vertex[0].z;

  for (i = 0; i < tri_model->num_tris; i++)
  {
    for (j = 0; j < 3; j++)
    {
      tx = tri_model->tris[i].vertex[j].x;
      ty = tri_model->tris[i].vertex[j].y;
      tz = tri_model->tris[i].vertex[j].z;

      /*Find the maximums*/
      if (tx > tri_model->max_x)
        tri_model->max_x = tx;
      if (ty > tri_model->max_y)
        tri_model->max_y = ty;
      if (tz > tri_model->max_z)
        tri_model->max_z = tz;

      /*Find the minimums*/
      if (tx < tri_model->min_x)
        tri_model->min_x = tx;
      if (ty < tri_model->min_y)
        tri_model->min_y = ty;
      if (tz < tri_model->min_z)
        tri_model->min_z = tz;
    }
  }
}

void Add_polygon(TRI_MODEL *tri_model, int ID, float density, float *verts, int num_tris)

{
  int tnum;

  tri_model->MU_ID = ID; // Material integer ID for the object, ID's I use are in ct_global_vars.h file
  tri_model->density = density;
  tri_model->num_tris = num_tris;
  dbug(1, "Triangles: %d\r\n", tri_model->num_tris);
  tri_model->tris = tri_vector(0, tri_model->num_tris);
  for (tnum = 0; tnum < num_tris; tnum++)
  {
    tri_model->tris[tnum].vertex[0].x = *verts++;
    tri_model->tris[tnum].vertex[0].y = *verts++;
    tri_model->tris[tnum].vertex[0].z = *verts++;
    tri_model->tris[tnum].vertex[1].x = *verts++;
    tri_model->tris[tnum].vertex[1].y = *verts++;
    tri_model->tris[tnum].vertex[1].z = *verts++;
    tri_model->tris[tnum].vertex[2].x = *verts++;
    tri_model->tris[tnum].vertex[2].y = *verts++;
    tri_model->tris[tnum].vertex[2].z = *verts++;
    tri_bbox(&tri_model->tris[tnum]);
  }
}

void add_triangles(double patch[4][4][3], TRI_MODEL *tri_model, int depth, int *t_num)

{
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];
  int tn1, tn2;

  if (depth > 0)
  {
    Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);
    add_triangles(ul_patch, tri_model, depth - 1, t_num);
    add_triangles(ur_patch, tri_model, depth - 1, t_num);
    add_triangles(dl_patch, tri_model, depth - 1, t_num);
    add_triangles(dr_patch, tri_model, depth - 1, t_num);
    return;
  }

  tn1 = *t_num;
  tn2 = (*t_num) + 1;
  t_num[0] += 2;

  tri_model->tris[tn1].vertex[0].x = patch[0][0][0];
  tri_model->tris[tn1].vertex[0].y = patch[0][0][1];
  tri_model->tris[tn1].vertex[0].z = patch[0][0][2];

  tri_model->tris[tn1].vertex[1].x = patch[0][3][0];
  tri_model->tris[tn1].vertex[1].y = patch[0][3][1];
  tri_model->tris[tn1].vertex[1].z = patch[0][3][2];

  tri_model->tris[tn1].vertex[2].x = patch[3][3][0];
  tri_model->tris[tn1].vertex[2].y = patch[3][3][1];
  tri_model->tris[tn1].vertex[2].z = patch[3][3][2];

  tri_model->tris[tn2].vertex[0].x = patch[3][3][0];
  tri_model->tris[tn2].vertex[0].y = patch[3][3][1];
  tri_model->tris[tn2].vertex[0].z = patch[3][3][2];

  tri_model->tris[tn2].vertex[1].x = patch[3][0][0];
  tri_model->tris[tn2].vertex[1].y = patch[3][0][1];
  tri_model->tris[tn2].vertex[1].z = patch[3][0][2];

  tri_model->tris[tn2].vertex[2].x = patch[0][0][0];
  tri_model->tris[tn2].vertex[2].y = patch[0][0][1];
  tri_model->tris[tn2].vertex[2].z = patch[0][0][2];

  tri_bbox(&tri_model->tris[tn1]);
  tri_bbox(&tri_model->tris[tn2]);
}

void BEZ2TRI(BEZIER_MODEL *bez_model, TRI_MODEL *tri_model)

{
  int i, depth, triangles_per_patch, t_num = 0;

  depth = get_refinement_depth(bez_model);

  triangles_per_patch = 2;
  for (i = 0; i < depth; i++)
    triangles_per_patch *= 4;

  tri_model->num_tris = triangles_per_patch * bez_model->num_patches;
  dbug(1, "Triangles: %d\r\n", tri_model->num_tris);
  tri_model->tris = tri_vector(0, tri_model->num_tris);

  for (i = 0; i < bez_model->num_patches; i++)
  {
    add_triangles(bez_model->patches[i].cpoints, tri_model, depth, &t_num);
  }
}

int AddItem(bvh_element **treepointer, int num, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, int *patches)
{
  int i;

  if ((*treepointer = (bvh_element *)malloc(sizeof(bvh_element))) == NULL)
    return 0;

  (*treepointer)->left = NULL;
  (*treepointer)->right = NULL;

  (*treepointer)->patches = ivector(0, num - 1);

  for (i = 0; i < num; i++)
    (*treepointer)->patches[i] = patches[i];

  (*treepointer)->num = num;
  (*treepointer)->xmin = xmin - EPSILON2;
  (*treepointer)->ymin = ymin - EPSILON2;
  (*treepointer)->zmin = zmin - EPSILON2;
  (*treepointer)->xmax = xmax + EPSILON2;
  (*treepointer)->ymax = ymax + EPSILON2;
  (*treepointer)->zmax = zmax + EPSILON2;
  return 1;
}

/* This function prints the elements in the queue */
void PrintTree(bvh_element *treepointer)
{
  int i;

  /* Traverse the left branch */
  if (treepointer->left != NULL)
    PrintTree(treepointer->left);

  /* Traverse the right branch */
  if (treepointer->right != NULL)
    PrintTree(treepointer->right);

  /* Print the current node */

  for (i = 0; i < treepointer->num; i++)
  {
    ;
  }
}

void CalcBVH(BEZIER_MODEL bez_model, int num, int *patches, float *xmin, float *xmax, float *ymin, float *ymax, float *zmin, float *zmax)
{
  int i;

  *xmin = 10000;
  *xmax = -10000;
  *ymin = 10000;
  *ymax = -10000;
  *zmin = 10000;
  *zmax = -10000;
  for (i = 0; i < num; i++)
  {
    if (bez_model.patches[patches[i]].minx < *xmin)
      *xmin = bez_model.patches[patches[i]].minx;
    if (bez_model.patches[patches[i]].miny < *ymin)
      *ymin = bez_model.patches[patches[i]].miny;
    if (bez_model.patches[patches[i]].minz < *zmin)
      *zmin = bez_model.patches[patches[i]].minz;

    if (bez_model.patches[patches[i]].maxx > *xmax)
      *xmax = bez_model.patches[patches[i]].maxx;
    if (bez_model.patches[patches[i]].maxy > *ymax)
      *ymax = bez_model.patches[patches[i]].maxy;
    if (bez_model.patches[patches[i]].maxz > *zmax)
      *zmax = bez_model.patches[patches[i]].maxz;
  }
}

void CreateBVH(bvh_element *treepointer, BEZIER_MODEL bez_model, int num, int *patches)
{
  int i;
  float minx, maxx;
  float miny, maxy;
  float minz, maxz;

  float dx, dy, dz;
  int axis; /*0 = x, 1 = y, 2 = z */
  int count1, count2;
  int count1x, count2x;
  int count1y, count2y;
  int count1z, count2z;
  int *right, *left;

  if (num == 1)
    return;

  right = ivector(0, num);
  left = ivector(0, num);

  minx = 10000;
  miny = 10000;
  minz = 10000;
  maxx = -10000;
  maxy = -10000;
  maxz = -10000;
  for (i = 0; i < num; i++)
  {
    if (bez_model.patches[patches[i]].minx < minx)
      minx = bez_model.patches[patches[i]].minx;
    if (bez_model.patches[patches[i]].miny < miny)
      miny = bez_model.patches[patches[i]].miny;
    if (bez_model.patches[patches[i]].minz < minz)
      minz = bez_model.patches[patches[i]].minz;

    if (bez_model.patches[patches[i]].maxx > maxx)
      maxx = bez_model.patches[patches[i]].maxx;
    if (bez_model.patches[patches[i]].maxy > maxy)
      maxy = bez_model.patches[patches[i]].maxy;
    if (bez_model.patches[patches[i]].maxz > maxz)
      maxz = bez_model.patches[patches[i]].maxz;
  }

  count1x = 0;
  count2x = 0;
  count1y = 0;
  count2y = 0;
  count1z = 0;
  count2z = 0;
  for (i = 0; i < num; i++)
  {
    if (bez_model.patches[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      count1x++;
    else
      count2x++;

    if (bez_model.patches[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      count1y++;
    else
      count2y++;

    if (bez_model.patches[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      count1z++;
    else
      count2z++;
  }

  dx = fabs(0.5 - (float)count1x / (count1x + count2x));
  dy = fabs(0.5 - (float)count1y / (count1y + count2y));
  dz = fabs(0.5 - (float)count1z / (count1z + count2z));

  axis = 0;
  if (dy < dx)
    axis = 1;
  if (dz < dy)
    axis = 2;

  if ((count1x == 0 || count2x == 0) && (count1y == 0 || count2y == 0) && (count1z == 0 || count2z == 0))
  {
    free_ivector(right, 0, num);
    free_ivector(left, 0, num);
    return;
  }

  if (axis == 0)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (bez_model.patches[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 1)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (bez_model.patches[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 2)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (bez_model.patches[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }

  if (count1 > 0)
  {
    CalcBVH(bez_model, count1, right, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->right, count1, minx, maxx, miny, maxy, minz, maxz, right);
    free_ivector(right, 0, num);
    CreateBVH(treepointer->right, bez_model, count1, treepointer->right->patches);
  }
  else
    free_ivector(right, 0, num);

  if (count2 > 0)
  {
    CalcBVH(bez_model, count2, left, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->left, count2, minx, maxx, miny, maxy, minz, maxz, left);
    free_ivector(left, 0, num);
    CreateBVH(treepointer->left, bez_model, count2, treepointer->left->patches);
  }
  else
    free_ivector(left, 0, num);
}

void FreeItem_BVH(bvh_element *treepointer)

{
  if (treepointer->left != NULL)
    FreeItem_BVH(treepointer->left);
  if (treepointer->right != NULL)
    FreeItem_BVH(treepointer->right);
  free_ivector(treepointer->patches, 0, treepointer->num - 1);
  free(treepointer);
  treepointer = NULL;
}

void Create_Bounding_Box(SURFACE surf, BEZIER_MODEL bez_model, int ID)
{
  int *patches, i;

  patches = ivector(0, bez_model.num_patches);
  for (i = 0; i < bez_model.num_patches; i++)
    patches[i] = i;

  if (treepointer_nrb[ID] != NULL)
    FreeItem_BVH(treepointer_nrb[ID]);

  AddItem(&treepointer_nrb[ID], bez_model.num_patches, surf.min_x, surf.max_x, surf.min_y, surf.max_y, surf.min_z, surf.max_z, patches);
  CreateBVH(treepointer_nrb[ID], bez_model, bez_model.num_patches, patches);

  free_ivector(patches, 0, bez_model.num_patches);
}

void CalcBVHCyl(listelement *listpointer, int num, int *patches, float *xmin, float *xmax, float *ymin, float *ymax, float *zmin, float *zmax)
{
  int i;

  *xmin = 10000;
  *xmax = -10000;
  *ymin = 10000;
  *ymax = -10000;
  *zmin = 10000;
  *zmax = -10000;
  for (i = 0; i < num; i++)
  {
    if (listpointer[patches[i]].minx < *xmin)
      *xmin = listpointer[patches[i]].minx;
    if (listpointer[patches[i]].miny < *ymin)
      *ymin = listpointer[patches[i]].miny;
    if (listpointer[patches[i]].minz < *zmin)
      *zmin = listpointer[patches[i]].minz;

    if (listpointer[patches[i]].maxx > *xmax)
      *xmax = listpointer[patches[i]].maxx;
    if (listpointer[patches[i]].maxy > *ymax)
      *ymax = listpointer[patches[i]].maxy;
    if (listpointer[patches[i]].maxz > *zmax)
      *zmax = listpointer[patches[i]].maxz;
  }
}

void CreateBVHCyl(bvh_element *treepointer, listelement *listpointer, int num, int *patches)
{
  int i;
  float minx, maxx;
  float miny, maxy;
  float minz, maxz;

  float dx, dy, dz;
  int axis; /*0 = x, 1 = y, 2 = z */
  int count1, count2;
  int count1x, count2x;
  int count1y, count2y;
  int count1z, count2z;
  int *right, *left;

  if (num == 1)
    return;

  right = ivector(0, num);
  left = ivector(0, num);

  minx = 10000;
  miny = 10000;
  minz = 10000;
  maxx = -10000;
  maxy = -10000;
  maxz = -10000;
  for (i = 0; i < num; i++)
  {
    if (listpointer[patches[i]].minx < minx)
      minx = listpointer[patches[i]].minx;
    if (listpointer[patches[i]].miny < miny)
      miny = listpointer[patches[i]].miny;
    if (listpointer[patches[i]].minz < minz)
      minz = listpointer[patches[i]].minz;

    if (listpointer[patches[i]].maxx > maxx)
      maxx = listpointer[patches[i]].maxx;
    if (listpointer[patches[i]].maxy > maxy)
      maxy = listpointer[patches[i]].maxy;
    if (listpointer[patches[i]].maxz > maxz)
      maxz = listpointer[patches[i]].maxz;
  }

  count1x = 0;
  count2x = 0;
  count1y = 0;
  count2y = 0;
  count1z = 0;
  count2z = 0;
  for (i = 0; i < num; i++)
  {
    if (listpointer[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      count1x++;
    else
      count2x++;

    if (listpointer[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      count1y++;
    else
      count2y++;

    if (listpointer[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      count1z++;
    else
      count2z++;
  }

  dx = fabs(0.5 - (float)count1x / (count1x + count2x));
  dy = fabs(0.5 - (float)count1y / (count1y + count2y));
  dz = fabs(0.5 - (float)count1z / (count1z + count2z));

  axis = 0;
  if (dy < dx)
    axis = 1;
  if (dz < dy)
    axis = 2;

  if ((count1x == 0 || count2x == 0) && (count1y == 0 || count2y == 0) && (count1z == 0 || count2z == 0))
  {
    free_ivector(right, 0, num);
    free_ivector(left, 0, num);
    return;
  }

  if (axis == 0)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (listpointer[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 1)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (listpointer[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 2)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (listpointer[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }

  if (count1 > 0)
  {
    CalcBVHCyl(listpointer, count1, right, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->right, count1, minx, maxx, miny, maxy, minz, maxz, right);
    free_ivector(right, 0, num);
    CreateBVHCyl(treepointer->right, listpointer, count1, treepointer->right->patches);
  }
  else
    free_ivector(right, 0, num);

  if (count2 > 0)
  {
    CalcBVHCyl(listpointer, count2, left, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->left, count2, minx, maxx, miny, maxy, minz, maxz, left);
    free_ivector(left, 0, num);
    CreateBVHCyl(treepointer->left, listpointer, count2, treepointer->left->patches);
  }
  else
    free_ivector(left, 0, num);
}

void Create_Bounding_Box_Cyl(SURFACE surf, listelement *listpointer, int begin_num, int end_num, int ID)
{
  int *patches, i;
  int num;

  num = end_num - begin_num + 1;
  patches = ivector(0, num);
  for (i = 0; i < num; i++)
    patches[i] = i + begin_num;

  treepointer_nrb[ID] = NULL;
  AddItem(&treepointer_nrb[ID], num, surf.min_x, surf.max_x, surf.min_y, surf.max_y, surf.min_z, surf.max_z, patches);

  CreateBVHCyl(treepointer_nrb[ID], listpointer, num, patches);

  free_ivector(patches, 0, num);
}

void Create_Bounding_Box_Cyl2(listelement *listpointer, int begin_num, int end_num)
{
  int *patches, i;
  int num;
  int ID = 0;

  float minx, maxx, miny, maxy, minz, maxz;
  num = end_num - begin_num + 1;
  patches = ivector(0, num);
  for (i = 0; i < num; i++)
    patches[i] = i + begin_num;

  minx = 10000;
  miny = 10000;
  minz = 10000;
  maxx = -10000;
  maxy = -10000;
  maxz = -10000;
  for (i = 0; i < num; i++)
  {
    if (listpointer[patches[i]].minx < minx)
      minx = listpointer[patches[i]].minx;
    if (listpointer[patches[i]].miny < miny)
      miny = listpointer[patches[i]].miny;
    if (listpointer[patches[i]].minz < minz)
      minz = listpointer[patches[i]].minz;

    if (listpointer[patches[i]].maxx > maxx)
      maxx = listpointer[patches[i]].maxx;
    if (listpointer[patches[i]].maxy > maxy)
      maxy = listpointer[patches[i]].maxy;
    if (listpointer[patches[i]].maxz > maxz)
      maxz = listpointer[patches[i]].maxz;
  }

  treepointer_nrb[ID] = NULL;
  AddItem(&treepointer_nrb[ID], num, minx, maxx, miny, maxy, minz, maxz, patches);

  CreateBVHCyl(treepointer_nrb[ID], listpointer, num, patches);

  free_ivector(patches, 0, num);
}

void CalcBVHTri(TRI_MODEL tmodel, int num, int *patches, float *xmin, float *xmax, float *ymin, float *ymax, float *zmin, float *zmax)
{
  int i;

  *xmin = 10000;
  *xmax = -10000;
  *ymin = 10000;
  *ymax = -10000;
  *zmin = 10000;
  *zmax = -10000;
  for (i = 0; i < num; i++)
  {
    if (tmodel.tris[patches[i]].minx < *xmin)
      *xmin = tmodel.tris[patches[i]].minx;
    if (tmodel.tris[patches[i]].miny < *ymin)
      *ymin = tmodel.tris[patches[i]].miny;
    if (tmodel.tris[patches[i]].minz < *zmin)
      *zmin = tmodel.tris[patches[i]].minz;

    if (tmodel.tris[patches[i]].maxx > *xmax)
      *xmax = tmodel.tris[patches[i]].maxx;
    if (tmodel.tris[patches[i]].maxy > *ymax)
      *ymax = tmodel.tris[patches[i]].maxy;
    if (tmodel.tris[patches[i]].maxz > *zmax)
      *zmax = tmodel.tris[patches[i]].maxz;
  }
}

void CreateBVHTri(bvh_element *treepointer, TRI_MODEL tmodel, int num, int *patches)
{
  int i, j;
  float minx, maxx;
  float miny, maxy;
  float minz, maxz;

  float dx, dy, dz;
  int axis; /*0 = x, 1 = y, 2 = z */
  int count1, count2;
  int count1x, count2x;
  int count1y, count2y;
  int count1z, count2z;
  int *right, *left;

  if (num == 1)
    return;

  right = ivector(0, num);
  left = ivector(0, num);

  minx = 10000;
  miny = 10000;
  minz = 10000;
  maxx = -10000;
  maxy = -10000;
  maxz = -10000;
  for (i = 0; i < num; i++)
  {
    for (j = 0; j < 3; j++)
    {
      if (tmodel.tris[patches[i]].minx < minx)
        minx = tmodel.tris[patches[i]].minx;
      if (tmodel.tris[patches[i]].miny < miny)
        miny = tmodel.tris[patches[i]].miny;
      if (tmodel.tris[patches[i]].minz < minz)
        minz = tmodel.tris[patches[i]].minz;

      if (tmodel.tris[patches[i]].maxx > maxx)
        maxx = tmodel.tris[patches[i]].maxx;
      if (tmodel.tris[patches[i]].maxy > maxy)
        maxy = tmodel.tris[patches[i]].maxy;
      if (tmodel.tris[patches[i]].maxz > maxz)
        maxz = tmodel.tris[patches[i]].maxz;
    }
  }

  count1x = 0;
  count2x = 0;
  count1y = 0;
  count2y = 0;
  count1z = 0;
  count2z = 0;
  for (i = 0; i < num; i++)
  {
    if (tmodel.tris[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      count1x++;
    else
      count2x++;

    if (tmodel.tris[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      count1y++;
    else
      count2y++;

    if (tmodel.tris[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      count1z++;
    else
      count2z++;
  }

  dx = fabs(0.5 - (float)count1x / (count1x + count2x));
  dy = fabs(0.5 - (float)count1y / (count1y + count2y));
  dz = fabs(0.5 - (float)count1z / (count1z + count2z));

  axis = 0;
  if (dy < dx)
    axis = 1;
  if (dz < dy)
    axis = 2;

  if ((count1x == 0 || count2x == 0) && (count1y == 0 || count2y == 0) && (count1z == 0 || count2z == 0))
  {
    free_ivector(right, 0, num);
    free_ivector(left, 0, num);
    return;
  }

  if (axis == 0)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (tmodel.tris[patches[i]].maxx < minx + (maxx - minx) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 1)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (tmodel.tris[patches[i]].maxy < miny + (maxy - miny) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }
  else if (axis == 2)
  {
    count1 = 0;
    count2 = 0;
    for (i = 0; i < num; i++)
    {
      if (tmodel.tris[patches[i]].maxz < minz + (maxz - minz) / 2.0)
      {
        right[count1] = patches[i];
        count1++;
      }
      else
      {
        left[count2] = patches[i];
        count2++;
      }
    }
  }

  if (count1 > 0)
  {
    CalcBVHTri(tmodel, count1, right, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->right, count1, minx, maxx, miny, maxy, minz, maxz, right);
    free_ivector(right, 0, num);
    CreateBVHTri(treepointer->right, tmodel, count1, treepointer->right->patches);
  }
  else
    free_ivector(right, 0, num);

  if (count2 > 0)
  {
    CalcBVHTri(tmodel, count2, left, &minx, &maxx, &miny, &maxy, &minz, &maxz);
    AddItem(&treepointer->left, count2, minx, maxx, miny, maxy, minz, maxz, left);
    free_ivector(left, 0, num);
    CreateBVHTri(treepointer->left, tmodel, count2, treepointer->left->patches);
  }
  else
    free_ivector(left, 0, num);
}

void Create_Bounding_Box_TRI(TRI_MODEL tmodel, int ID)
{
  int *patches, i;
  int num;

  num = tmodel.num_tris;

  patches = ivector(0, num);
  for (i = 0; i < num; i++)
    patches[i] = i;

  if (treepointer_tri[ID] != NULL)
  {
    dbug(1, "Free TP: %d\n", ID);
    FreeItem_BVH(treepointer_tri[ID]);
  }

  treepointer_tri[ID] = NULL;

  CalcBVHTri(tmodel, num, patches, &tmodel.min_x, &tmodel.max_x, &tmodel.min_y, &tmodel.max_y, &tmodel.min_z, &tmodel.max_z);
  dbug(1, "Create TP: %d\n", ID);
  AddItem(&treepointer_tri[ID], num, tmodel.min_x, tmodel.max_x, tmodel.min_y, tmodel.max_y, tmodel.min_z, tmodel.max_z, patches);
  CreateBVHTri(treepointer_tri[ID], tmodel, num, patches);
  free_ivector(patches, 0, num);
}
/*------------------------------------End Tree Definition--------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Abort(Mesg)                \
  {                                \
    fprintf(stderr, "%s\n", Mesg); \
    exit(1);                       \
  }
#define p_degree 3
#define q_degree 3

/*SUBROUTINE SAVE_TO_FILE(array,xdim,ydim,zdim,name,name_length)

  **************************************************************

  This subroutine saves images as binary files with NO HEADER.
  Each voxel value in output image is a
  4 BYTE FLOATING POINT number.

  **************************************************************
  This is simply a subroutine which saves the phantom to a file.
  This subroutine can be modified by the user to save the
  phantom any whatever format the user wants. The purpose of
  this seperate subroutine is to allow various users to
  modify the file format of the saved phantom without having
  to modify any of the other files containing the phantom code.
  **************************************************************

  array:           contains the data to be saved
  xdim,ydim,zdim:  dimensions of data
  name:            name of file
  --------------------------------------------------------------
*/
typedef float FLOATTYPE;

/* Not currently used, but here if you need them: */
/* typedef unsigned char BYTETYPE;    value range: 0 to 255 */
/* typedef short INTTYPE;             value range: -32768 to 32767 */

int flag = 0;
/*-------------------------------------------------------------*/
void SAVE_TO_FILE(float *array, int xdim, int ydim, int zdim, char *name)
{
  int n, TotalPix;
  char new_name[64];
  FILE *fp_out;
  /*
  ** ADD EXTENSION TO FILENAME:
  */
  strcpy(new_name, name);
  strcat(new_name, ".bin"); /* add extension to file name*/
  unlink(new_name);

  if (!flag)
  {
    printf("\nPhantoms saved as %i x %i x %i raw images (32 bit float) with no header.\n", xdim, ydim, zdim);
    flag = 1;
  }

  /*
  ** WRITE TO FILE:
  */
  TotalPix = xdim * ydim * zdim; /* total number of pixels in 3D image */
  if ((fp_out = fopen(new_name, "wb")) == NULL)
    Abort("Cannot open raw output file");

  n = fwrite(array, sizeof(FLOATTYPE), TotalPix, fp_out);
  if (n != TotalPix)
  {
    printf("Error : fwrite return %d\n", n);
    Abort("Failure writing pixels to output image\n");
    unlink(new_name);
  }
  fclose(fp_out);
}

/*----------------------------------------------------------------*/
void Calc_extents(SURFACE *nrb_model)
/*----------------------------------------------------------------*/
/*------------------------------------------------------------------
**  This subroutine is used to set the maximum and minimum dimensions
**  of a given NURBS model.  These dimensions are used to pre-render
**  the NURBS surface into its own image before combining the image
**  into the final output image.
**------------------------------------------------------------------
*/
{
  int i, j;
  float tx, ty, tz;

  /*Initialize the maximums and minimums to the first control point*/
  nrb_model->max_x = nrb_model->net.Pw[0][0].x;
  nrb_model->max_y = nrb_model->net.Pw[0][0].y;
  nrb_model->max_z = nrb_model->net.Pw[0][0].z;

  nrb_model->min_x = nrb_model->net.Pw[0][0].x;
  nrb_model->min_y = nrb_model->net.Pw[0][0].y;
  nrb_model->min_z = nrb_model->net.Pw[0][0].z;

  /*Cycle through the model's control points*/
  for (i = 0; i < nrb_model->net.m; i++)
  {
    for (j = 0; j < nrb_model->net.n; j++)
    {
      tx = nrb_model->net.Pw[j][i].x;
      ty = nrb_model->net.Pw[j][i].y;
      tz = nrb_model->net.Pw[j][i].z;

      /*Find the maximums*/
      if (tx > nrb_model->max_x)
        nrb_model->max_x = tx;
      if (ty > nrb_model->max_y)
        nrb_model->max_y = ty;
      if (tz > nrb_model->max_z)
        nrb_model->max_z = tz;

      /*Find the minimums*/
      if (tx < nrb_model->min_x)
        nrb_model->min_x = tx;
      if (ty < nrb_model->min_y)
        nrb_model->min_y = ty;
      if (tz < nrb_model->min_z)
        nrb_model->min_z = tz;
    }
  }
}

/*--------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
void Allocate_NURBS(SURFACE *nrb_model, int n, int m)
/*------------------------------------------------------------------------------
**  This subroutine allocates memory for the specified NURBS surface
**    m,n:      number of control points in the v and u directions of the surface
**-------------------------------------------------------------------------------
*/
{
  int i, j;

  /*Setup NURB surfaces of the heart*/
  nrb_model->net.Pw = cp_matrix(0, n - 1, 0, m - 1);
  nrb_model->net.n = n;
  nrb_model->net.m = m;
  nrb_model->knu.U = vector(0, n + p_degree);
  nrb_model->knu.m = n + p_degree;
  nrb_model->knv.U = vector(0, m + q_degree);
  nrb_model->knv.m = m + q_degree;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      nrb_model->net.Pw[j][i].x = 0.0;
      nrb_model->net.Pw[j][i].y = 0.0;
      nrb_model->net.Pw[j][i].z = 0.0;
      nrb_model->net.Pw[j][i].w = 0.0;
    }
}
/*-------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
void Free_NURBS(SURFACE *nrb_model)
/*------------------------------------------------------------------------------
**  This subroutine frees the memory for the specified NURBS surface
**-------------------------------------------------------------------------------
*/
{
  free_cpmatrix(nrb_model->net.Pw, 0, nrb_model->net.n - 1, 0, nrb_model->net.m - 1);
  free_vector(nrb_model->knu.U, 0, nrb_model->net.n + p_degree);
  free_vector(nrb_model->knv.U, 0, nrb_model->net.m + q_degree);
}
/*-------------------------------------------------------------------------------*/

void print_bvh(bvh_element *treepointer, int depth)

{
  int i;
  if (treepointer->right != NULL)
    print_bvh(treepointer->right, depth + 1);
  if (treepointer->left != NULL)
    print_bvh(treepointer->left, depth + 1);

  if ((treepointer->right == NULL) && (treepointer->left == NULL))
    for (i = 0; i < treepointer->num; i++)
      dbug(3, "patches[%d]: %d\r\n", i, treepointer->patches[i]);
  dbug(3, "depth: %d\r\n", depth);
  dbug(3, "x: %f %f \r\ny: %f %f \r\nz: %f %f \r\n\n", treepointer->xmin, treepointer->xmax, treepointer->ymin, treepointer->ymax, treepointer->zmin, treepointer->zmax);
}

void print_vec(float *vec)

{
  printf("print_vec:  [%f %f %f]\r\n", vec[0], vec[1], vec[2]);
}

void print_verts(TRIANGLE *tris, int num_tris)

{
  int i;

  for (i = 0; i < num_tris; i++)
  {
    printf("Triangle #%d\r\n", i);
    printf(" Vert 1: %f %f %f\r\n", tris[i].vertex[0].x, tris[i].vertex[0].y, tris[i].vertex[0].z);
    printf(" Vert 2: %f %f %f\r\n", tris[i].vertex[1].x, tris[i].vertex[1].y, tris[i].vertex[1].z);
    printf(" Vert 3: %f %f %f\r\n\n", tris[i].vertex[2].x, tris[i].vertex[2].y, tris[i].vertex[2].z);
  }
}

void print_poly(int model_num)

{
  TRI_MODEL *tmodel = &tri_model[model_num];

  printf("Triangles:");
  print_verts(tmodel->tris, tmodel->num_tris);
  printf("num_triangles: %d\r\n", tmodel->num_tris);
  printf("minx: %f\r\n", tmodel->min_x);
  printf("maxx: %f\r\n", tmodel->max_x);
  printf("miny: %f\r\n", tmodel->min_y);
  printf("maxy: %f\r\n", tmodel->max_y);
  printf("minz: %f\r\n", tmodel->min_z);
  printf("maxz: %f\r\n", tmodel->max_z);
  print_bvh(treepointer_tri[model_num], 0);
}

DLLEXPORT void clear_polygonalized_phantom(int num_polygons)

{
  int i;

  /*
    JDP : Here we need to make room for max_num_models polygons (currently set to 5000).
    Otherwise, you can run once with 10 polygons and then when you run with 20 later, it won't have room --> possible crash.
  */

  if (treepointer_tri == NULL)
  {
    dbug(1, " Allocating space for bvh trees\n");
    treepointer_tri = (bvh_element **)malloc(max_num_models * sizeof(bvh_element *));
    for (i = 0; i < max_num_models; i++)
      treepointer_tri[i] = NULL;
  }

  if (tri_model == NULL)
  {
    dbug(1, " Allocating space for tri_models\n");
    tri_model = (TRI_MODEL *)malloc(max_num_models * sizeof(TRI_MODEL));
    for (i = 0; i < max_num_models; i++)
      tri_model[i].tris = NULL;
  }
  else
  {
    dbug(1, " Clearing all tri_models: ");
    for (i = 0; i < max_num_models; i++)
    {
      if (tri_model[i].tris != NULL)
      {
        dbug(1, "%d ", i);
        free_tri_vector(tri_model[i].tris, 0, tri_model[i].num_tris);
        tri_model[i].tris = NULL;
      }
    }
    dbug(1, "\n");
  }

  NUM_POLY = 0;
}

DLLEXPORT void pass_polygon_to_c(float *vertices, int num_triangles, float density, int ID)
/*  This routine passes a polygonalized object to the projector */
{
  int model_num;

  model_num = NUM_POLY;
  // dbug(1,"\rmodel: %d\n",model_num);
  Add_polygon(&tri_model[model_num], ID, density, vertices, num_triangles);
  Calc_extents_tri(&tri_model[model_num]);
  Create_Bounding_Box_TRI(tri_model[model_num], model_num);
  dbug(2, "f\n");
  NUM_POLY++;
  // print_poly(model_num);
}

int isEmptyString(char *line)
{
  while (*line)
  {
    if (!isspace(*line++))
    {
      return 0;
    }
  }
  return 1;
}

DLLEXPORT void Parse_Phantom(char *filename, int *materials, float *coord_origin_offset, float scale) //  JDP
/*  This routine parses a NURBS phantom file output from the NCAT */
{
  FILE *fp;
  int i, j, k;
  int n, m;
  float temp, tx, ty, tz;
  char line[160];
  int tmp = 0;
  int ID;

  use_tri_model = 0;
  if ((fp = fopen(filename, "r")) == NULL)
    Abort("Can not open anatomy nurbs datafile");

  dbug(0, "\n\rStarting to parse nCAT phantom.\r\n");

  for (i = 0; i < 3; i++)
    dbug(1, "\n\roffset[%d]:  %f \r\n", i, coord_origin_offset[i]);

  // JDP

  if (treepointer_nrb == NULL)
  {
    treepointer_nrb = (bvh_element **)malloc(max_num_models * sizeof(bvh_element *));
    for (i = 0; i < max_num_models; i++)
      treepointer_nrb[i] = NULL;
  }

  if (nrb_model == NULL)
  {
    nrb_model = (SURFACE *)malloc(max_num_models * sizeof(SURFACE));
    for (i = 0; i < max_num_models; i++)
    {
      nrb_model[i].net.Pw = NULL;
      nrb_model[i].knu.U = NULL;
      nrb_model[i].knv.U = NULL;
    }
  }
  else
  {
    for (i = 0; i < max_num_models; i++)
    {
      if (nrb_model[i].net.Pw != NULL)
      {
        Free_NURBS(&nrb_model[i]);
        nrb_model[i].net.Pw = NULL;
      }
    }
  }

  if (bez_model == NULL)
  {
    bez_model = (BEZIER_MODEL *)malloc(max_num_models * sizeof(BEZIER_MODEL));
    for (i = 0; i < max_num_models; i++)
    {
      bez_model[i].patches = NULL;
    }
  }
  else
  {
    for (i = 0; i < max_num_models; i++)
    {
      if (bez_model[i].patches != NULL)
      {
        free_bpvector(bez_model[i].patches, 0, bez_model[i].num_patches);
        bez_model[i].patches = NULL;
      }
    }
  }

  if (tri_model == NULL)
  {
    tri_model = (TRI_MODEL *)malloc(max_num_models * sizeof(TRI_MODEL));
    for (i = 0; i < max_num_models; i++)
    {
      tri_model[i].tris = NULL;
    }
  }
  else
  {
    for (i = 0; i < max_num_models; i++)
    {
      if (tri_model[i].tris != NULL)
      {
        free_tri_vector(tri_model[i].tris, 0, tri_model[i].num_tris);
        tri_model[i].tris = NULL;
      }
    }
  }
  // end JDP

  //  tmp = fscanf(fp, "%s", line);  //Title line
  do
  {
    char *p = fgets(line, sizeof(line) - 1, fp); // Title line
    if (!p)
    {
      tmp = EOF;
      break;
    }
  } while (line[0] == '%' || isEmptyString(line));
  k = 0;
  while (tmp != EOF)
  {
    dbug(1, "\rmodel: %d\n", k);
    if (k > max_num_models)
    {
      dbug(0, "Error: Maximum number of nurbs models (%d) exceeded.\r\nSee Parse_Phantom routine in nCAT_main.c\r\n", max_num_models);
      exit(1);
    }
    fscanf(fp, "%i", &ID);
    nrb_model[k].MU_ID = ID; // Material integer ID for the object, ID's I use are in ct_global_vars.h file
    materials[ID] = 1;       // Place a 1 in the materials array to let it know material[ID] is included

    /* Read in M and N parameters */
    fscanf(fp, "%i", &m);
    fscanf(fp, "%s", line);
    fscanf(fp, "%i", &n);
    fscanf(fp, "%s", line);

    /* Setup the surface */
    Allocate_NURBS(&nrb_model[k], n, m);

    /* Read in U Knot Vector */
    fscanf(fp, "%s", line);
    fscanf(fp, "%s", line);
    fscanf(fp, "%s", line);
    for (i = 0; i <= n + p_degree; i++)
    {
      fscanf(fp, "%f", &temp);
      nrb_model[k].knu.U[i] = temp;
    }

    /* Read in V Knot Vector */
    fscanf(fp, "%s", line);
    fscanf(fp, "%s", line);
    fscanf(fp, "%s", line);
    for (i = 0; i <= m + q_degree; i++)
    {
      fscanf(fp, "%f", &temp);
      nrb_model[k].knv.U[i] = temp;
    }

    /*Readin Control Points*/
    fscanf(fp, "%s", line);
    fscanf(fp, "%s", line);
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
      {
        fscanf(fp, "%f %f %f", &tx, &ty, &tz);
        nrb_model[k].net.Pw[j][i].x = scale * tx + coord_origin_offset[0]; // FIXME: any way to get these from nrb file or ncat cfg?
        nrb_model[k].net.Pw[j][i].y = scale * ty + coord_origin_offset[1]; // (origin_offsets)
        nrb_model[k].net.Pw[j][i].z = tz + coord_origin_offset[2];
        nrb_model[k].net.Pw[j][i].w = 0.0;
      }

    dbug(2, "e\n");
    SETUP_BEZIER_MODEL(nrb_model[k], &bez_model[k]);
    // dbug(-1,"obj_num:%d",k);
    // if (k==161) debug_flag = 1;
    SPLINE2BEZ(&nrb_model[k], &bez_model[k]);
    // if (k==161) debug_flag = 0;
    Calc_extents(&nrb_model[k]);
    // BEZ2TRI(&bez_model[k],&tri_model[k]);
    // Create_Bounding_Box_TRI(tri_model[k], k);
    Create_Bounding_Box(nrb_model[k], bez_model[k], k);
    dbug(2, "f\n");

    k++;
    tmp = fscanf(fp, "%s", line);
  }
  NUM_NRB = k;
  fclose(fp);
  NUM_MODELS = NUM_NRB;
  
  dbug(0, "Done parsing nCAT phantom.\r\n");
}

int Intersect_segments(LINE_SEG A, LINE_SEG B) /*A = higher priority segment*/
{
  int int_type;

  if (B.x1 <= A.x1 && B.x2 <= A.x2 && B.x2 >= A.x1)
    int_type = 1;
  else if (B.x1 < A.x1 && B.x2 > A.x2)
    int_type = 2;
  else if (B.x1 >= A.x1 && B.x2 >= A.x2 && B.x1 <= A.x2)
    int_type = 3;
  else if (B.x1 >= A.x1 && B.x2 <= A.x2)
    int_type = 4;
  else
    int_type = 0;

  return int_type;
}

void Break_segment(SEG_ARRAY lines, LINE_SEG A, SEG_ARRAY *segments, int line_loc)
{
  int i, flag;
  LINE_SEG C;

  i = line_loc + 1;
  flag = 0;

  while (flag == 0 && i < lines.length)
  {
    flag = Intersect_segments(lines.sp[i], A);
    i++;
  }

  if (flag > 0)
    switch (flag)
    {
    case 1:
      C.x1 = A.x1;
      C.x2 = lines.sp[i - 1].x1;
      C.organ_id = A.organ_id;
      Break_segment(lines, C, segments, i - 1);
      break;
    case 2:
      C.x1 = A.x1;
      C.x2 = lines.sp[i - 1].x1;
      C.organ_id = A.organ_id;
      Break_segment(lines, C, segments, i - 1);

      C.x1 = lines.sp[i - 1].x2;
      C.x2 = A.x2;
      C.organ_id = A.organ_id;
      Break_segment(lines, C, segments, i - 1);
      break;
    case 3:
      C.x1 = lines.sp[i - 1].x2;
      C.x2 = A.x2;
      C.organ_id = A.organ_id;
      Break_segment(lines, C, segments, i - 1);
      break;
    case 4:
      C.x1 = 0;
      C.x2 = 0;
      C.organ_id = A.organ_id;
      break;
    default:
      break;
    };

  if (flag == 0)
  {
    segments->sp[segments->length].x1 = A.x1;
    segments->sp[segments->length].x2 = A.x2;
    segments->sp[segments->length].organ_id = A.organ_id;
    segments->length += 1;
  }
}

void Break_segment2(SEG_ARRAY lines, SEG_ARRAY *segments)

{
  /*
    The purpose of this function is to replace the functionality of Break_segment, but be much faster.
    One difference is that Break_segment does not produce a segment for the last line, whereas this one does.
    As a result, we need to change Calc_line_int to the simpler Calc_line_int2 when using Break_segment2.

    This function is many many times faster than Break_segment.  The way it works is the following:
    1) All intersections are listed separately (rather than in pairs)
    2) The intersections are sorted by position along the ray
    3) There is a loop over the intersections. During the loop, we keep track (in a prioritized list... with the priority defined by the organ_id) of all the relevant segments which we are inside as we march along.  A segment is considered relevant if it will be the highest priority segment for at least a portion of its length.  During the loop:
        - If the intersection is at the start of the segment (sense==0), we first find the position of the segment in the priority list.  If the new segment ends after the one with just higher priority, we add it to the list.  We also remove any existing segments from the list if they are completely covered by the new segment (the new segment has higher priority and a higher endpoint location).  Finally, if the new segment has the highest priority of any in the list, we complete any segment that is in progress and start a new one at the current location with the new priority.
        - If the intersection is at the end of the segment (sense == 1), we ignore it unless it is the top priority segment.  If it is the top priority segment, we finish the segment, remove it from the list, and start a new segment at the next highest priority.  if there is no other segment in the list, we find ourselves outside the body and we keep track of this with a flag (in_body=0) so that we will know we don't need to end a segment the next time we start one.
  */

  HIT_ARRAY critical_pts;
  PRIORITY_ARRAY priority_list;
  int i, j, k, this, old_next, in_body = 0;
  HIT intersection;

  // store the intersections individually instead of in pairs, with a flag to indicate if its a start or an end

  //  dbug(1,"length(lines): %d\n",lines.length);
  if (0)
    for (i = 0; i < lines.length; i++)
    {
      dbug(-1, "lines.sp[%d].x1: %f\n", i, lines.sp[i].x1);
      dbug(-1, "lines.sp[%d].x2: %f\n", i, lines.sp[i].x2);
      dbug(-1, "lines.sp[%d].organ_id: %d\n\n", i, lines.sp[i].organ_id);
    }
  critical_pts.length = lines.length * 2;
  for (i = 0; i < lines.length; i++)
  {
    critical_pts.sp[i * 2].x = lines.sp[i].x1;
    critical_pts.sp[i * 2].x2 = lines.sp[i].x2;
    critical_pts.sp[i * 2].sense = 0; // indicates start of a segment
    critical_pts.sp[i * 2].organ_id = lines.sp[i].organ_id;

    critical_pts.sp[i * 2 + 1].x = lines.sp[i].x2;
    critical_pts.sp[i * 2 + 1].sense = 1; // indicates end of a segment
    critical_pts.sp[i * 2 + 1].organ_id = lines.sp[i].organ_id;
  }
  // sort them
  qsort(critical_pts.sp, critical_pts.length, sizeof(HIT), comp_intersections);
  if (0)
    for (i = 0; i < critical_pts.length; i++)
    {
      printf("critical_pts.sp[%d].x: %f\n", i, critical_pts.sp[i].x);
      printf("critical_pts.sp[%d].x2: %f\n", i, critical_pts.sp[i].x2);
      printf("critical_pts.sp[%d].sense: %d\n", i, critical_pts.sp[i].sense);
      printf("critical_pts.sp[%d].organ_id: %d\n\n", i, critical_pts.sp[i].organ_id);
    }
  // loop through the segments
  priority_list.length = 1;
  priority_list.pp[0].next = -1;
  priority_list.pp[0].x = -1;
  priority_list.pp[0].organ_id = NUM_MODELS + 1;
  for (i = 0; i < critical_pts.length; i++)
  {
    intersection = critical_pts.sp[i];
    if (critical_pts.sp[i].sense == 0)
    { // start of a segment
      this = 0;
      while ((priority_list.pp[this].next != -1) && (priority_list.pp[priority_list.pp[this].next].organ_id > intersection.organ_id)) {
        this = priority_list.pp[this].next;
        if (this == 0)  break;  // To avoid infinite loop, Mingye Wu, 7/15/2024
      }
      // dbug(1," this:%d i:%d\n",this,i);
      if (this == 0) // we  start a new segment (finishing the existing one, if needed), since this is highest priority
      {
        if (in_body)
          segments->sp[segments->length - 1].x2 = intersection.x;

        segments->sp[segments->length].x1 = intersection.x; // x2 will be filled in later
        segments->sp[segments->length].organ_id = intersection.organ_id;
        segments->length++;
        in_body = 1;
        // dbug(1,"A:%d %d\n",i, intersection.organ_id);
      }
      // should now be pointing to the one I want to put it after.
      if (priority_list.pp[this].x < intersection.x2) /* Note: if these are equal, we don't bother inserting anything... we want lower priority lines
               to have endpoints that are strictly greater than all higher priority lines wince the sorting
               doesn't put these in a predictable order */
      {
        old_next = priority_list.pp[this].next;
        priority_list.pp[this].next = priority_list.length;
        // dbug(1," //// %d %d ////\n", this, priority_list.pp[this].next);
        this = priority_list.length;

        // add a new node
        priority_list.length++;
        priority_list.pp[this].organ_id = intersection.organ_id;
        priority_list.pp[this].x = intersection.x2;

        // remove any obsolete ones by updating old_next if needed
        while ((old_next != -1) && (priority_list.pp[old_next].x <= intersection.x2)) /* Note: if they are equal, we make sure to get rid of the lower priority one.
                             Again, we want to keep lower priority ones only if they have endpoints that
                             are strictly greater than the higher priority one.*/
          old_next = priority_list.pp[old_next].next;
        priority_list.pp[this].next = old_next;
        // dbug(1," -- %d %d --\n", this, priority_list.pp[this].next);
      }
    }
    else
    { // end of a segment
      // dbug(1," -----------------------------------\n");
      // dbug(1," -- %d %d ---------------------------------\n", priority_list.pp[priority_list.pp[0].next].organ_id, priority_list.pp[0].next);
      // dbug(1," -- %d ---------------------------------\n", intersection.organ_id);
      // dbug(1," -----------------------------------\n");
      if (priority_list.pp[priority_list.pp[0].next].organ_id == intersection.organ_id) // we finish the segment, then start a new segment
      {
        // first remove from list
        priority_list.pp[0].next = priority_list.pp[priority_list.pp[0].next].next;
        // dbug(1," ** %d **\n",  priority_list.pp[0].next);
        if (priority_list.pp[0].next == -1)
          in_body = 0;

        segments->sp[segments->length - 1].x2 = intersection.x;
        if (in_body)
        {
          segments->sp[segments->length].x1 = intersection.x; // x2 will be filled in later
          segments->sp[segments->length].organ_id = priority_list.pp[priority_list.pp[0].next].organ_id;
          segments->length++;
          // dbug(1,"B:%d %d\n",i,priority_list.pp[priority_list.pp[0].next].organ_id);
        }
      }
    }
  }

  if (0)
    for (i = 0; i < segments->length; i++)
    {
      printf("segments->sp[%d].x1: %f\n", i, segments->sp[i].x1);
      printf("segments->sp[%d].x2: %f\n", i, segments->sp[i].x2);
      printf("segments->sp[%d].organ_id: %d\n\n", i, segments->sp[i].organ_id);
    }
}

void Calc_line_int(SEG_ARRAY lines, SEG_ARRAY segments, float mu_table[MAXIMUM_NUMBER_OF_MATERIALS][300], int energy_ID, float *line_int)
{

  // THIS FUNCTION IS NOW OBSOLETE

  int i;
  float dist;
  float atten;

  *line_int = 0.0;
  if (lines.length != 0)
  {
    for (i = 0; i < segments.length; i++)
    {
      dist = segments.sp[i].x2 - segments.sp[i].x1;
      atten = mu_table[nrb_model[segments.sp[i].organ_id].MU_ID][energy_ID];
      *line_int += (atten * dist);
      // dbug(1,"dist[%d]: %1.14f  atten: %1.14f\r\n",i,dist,atten);
    }

    dist = lines.sp[lines.length - 1].x2 - lines.sp[lines.length - 1].x1;
    dbug(1, "dist[X]: %1.14f\r\n", dist);
    atten = mu_table[nrb_model[lines.sp[lines.length - 1].organ_id].MU_ID][energy_ID];
    *line_int += (atten * dist);
  }
}

void Calc_line_int_tri(SEG_ARRAY lines, SEG_ARRAY segments, float **mu_table, int energy_ID, float *line_int)
{

  // THIS FUNCTION IS NOW OBSOLETE

  int i;
  float dist;
  float atten;

  *line_int = 0.0;
  if (lines.length != 0)
  {
    for (i = 0; i < segments.length; i++)
    {
      dist = segments.sp[i].x2 - segments.sp[i].x1;
      atten = mu_table[tri_model[segments.sp[i].organ_id].MU_ID][energy_ID] * tri_model[segments.sp[i].organ_id].density;
      *line_int += (atten * dist);
    }

    dist = lines.sp[lines.length - 1].x2 - lines.sp[lines.length - 1].x1;
    atten = mu_table[tri_model[lines.sp[lines.length - 1].organ_id].MU_ID][energy_ID] * tri_model[lines.sp[lines.length - 1].organ_id].density;
    *line_int += (atten * dist);
    //    dbug(2,"line_int: %f atten: %f dist: %f ind: %d mu: %f\r\n",*line_int,atten,dist,tri_model[lines.sp[lines.length-1].organ_id].MU_ID,mu_table[tri_model[lines.sp[lines.length-1].organ_id].MU_ID][energy_ID]);
  }
}

void Calc_line_int2(SEG_ARRAY segments, float mu_table[MAXIMUM_NUMBER_OF_MATERIALS][300], int energy_ID, float *line_int)
{
  int i;
  float dist;
  float atten;

  *line_int = 0.0;
  for (i = 0; i < segments.length; i++)
  {
    // dbug(1,"%d %d %d\r\n",i, segments.sp[i].organ_id, nrb_model[segments.sp[i].organ_id].MU_ID);
    // dbug(0,"%d %d\r\n",i, segments.sp[i].organ_id);
    // dbug(0,"%d\r\n",nrb_model[segments.sp[i].organ_id].MU_ID);
    dist = segments.sp[i].x2 - segments.sp[i].x1;
    atten = mu_table[nrb_model[segments.sp[i].organ_id].MU_ID][energy_ID];
    *line_int += (atten * dist);
    // dbug(1,"dist[%d]: %1.14f  atten: %1.14f\r\n",i,dist,atten);
  }
}

void Calc_line_int2_tri(SEG_ARRAY segments, float **mu_table, int energy_ID, float *line_int)
{
  int i;
  float dist;
  float atten;

  *line_int = 0.0;
  for (i = 0; i < segments.length; i++)
  {
    dist = segments.sp[i].x2 - segments.sp[i].x1;
    atten = mu_table[tri_model[segments.sp[i].organ_id].MU_ID][energy_ID] * tri_model[segments.sp[i].organ_id].density;
    *line_int += (atten * dist);
  }
}

void hull_split_u(double P[4][4][3], double Q[4][4][3], double R[4][4][3])
{
  int i, iv;

  for (iv = 3; iv >= 0; iv--)
  {
    for (i = 2; i >= 0; i--)
    {
      Q[0][iv][i] = P[0][iv][i];
      Q[1][iv][i] = (P[0][iv][i] + P[1][iv][i]) / 2.0;
      Q[2][iv][i] = Q[1][iv][i] / 2.0 + (P[1][iv][i] + P[2][iv][i]) / 4.0;

      R[3][iv][i] = P[3][iv][i];
      R[2][iv][i] = (P[2][iv][i] + P[3][iv][i]) / 2.0;
      R[1][iv][i] = R[2][iv][i] / 2.0 + (P[1][iv][i] + P[2][iv][i]) / 4.0;

      Q[3][iv][i] = (Q[2][iv][i] + R[1][iv][i]) / 2.0;
      R[0][iv][i] = Q[3][iv][i];
    }
  }
}

void hull_split_v(double P[4][4][3], double Q[4][4][3], double R[4][4][3])
{
  int i, iu;

  for (i = 2; i >= 0; i--)
  {

    for (iu = 3; iu >= 0; iu--)
    {
      Q[iu][0][i] = P[iu][0][i];
      Q[iu][1][i] = (P[iu][0][i] + P[iu][1][i]) / 2.0;
      Q[iu][2][i] = Q[iu][1][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      R[iu][3][i] = P[iu][3][i];
      R[iu][2][i] = (P[iu][2][i] + P[iu][3][i]) / 2.0;
      R[iu][1][i] = R[iu][2][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      Q[iu][3][i] = (Q[iu][2][i] + R[iu][1][i]) / 2.0;
      R[iu][0][i] = Q[iu][3][i];
    }
  }
}

void find_bounds(double cpoints[4][4][3], double extent[3][2])

{
  int i, j;
  double *tmp;
  for (i = 0; i < 3; i++)
  {
    extent[i][0] = 100000;
    extent[i][1] = -100000;
  }
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      tmp = cpoints[i][j];

      if (tmp[0] < extent[0][0])
        extent[0][0] = tmp[0];
      if (tmp[0] > extent[0][1])
        extent[0][1] = tmp[0];

      if (tmp[1] < extent[1][0])
        extent[1][0] = tmp[1];
      if (tmp[1] > extent[1][1])
        extent[1][1] = tmp[1];

      if (tmp[2] < extent[2][0])
        extent[2][0] = tmp[2];
      if (tmp[2] > extent[2][1])
        extent[2][1] = tmp[2];
    }
}

void hull_split_v_bnds(double P[4][4][3], double Q[4][4][3], double R[4][4][3], double Q_ext[3][2], double R_ext[3][2])
{

  // THIS FUNCTION IS NO LONGER USED

  int i, iu;

  for (i = 2; i >= 0; i--)
  {
    Q_ext[i][0] = 100000;
    Q_ext[i][1] = -100000;

    R_ext[i][0] = 100000;
    R_ext[i][1] = -100000;

    for (iu = 3; iu >= 0; iu--)
    {
      Q[iu][0][i] = P[iu][0][i];
      Q[iu][1][i] = (P[iu][0][i] + P[iu][1][i]) / 2.0;
      Q[iu][2][i] = Q[iu][1][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      R[iu][3][i] = P[iu][3][i];
      R[iu][2][i] = (P[iu][2][i] + P[iu][3][i]) / 2.0;
      R[iu][1][i] = R[iu][2][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      Q[iu][3][i] = (Q[iu][2][i] + R[iu][1][i]) / 2.0;
      R[iu][0][i] = Q[iu][3][i];

      if (Q[iu][0][i] < Q_ext[i][0])
        Q_ext[i][0] = Q[iu][0][i];
      if (Q[iu][1][i] < Q_ext[i][0])
        Q_ext[i][0] = Q[iu][1][i];
      if (Q[iu][2][i] < Q_ext[i][0])
        Q_ext[i][0] = Q[iu][2][i];
      if (Q[iu][3][i] < Q_ext[i][0])
        Q_ext[i][0] = Q[iu][3][i];

      if (Q[iu][0][i] > Q_ext[i][1])
        Q_ext[i][1] = Q[iu][0][i];
      if (Q[iu][1][i] > Q_ext[i][1])
        Q_ext[i][1] = Q[iu][1][i];
      if (Q[iu][2][i] > Q_ext[i][1])
        Q_ext[i][1] = Q[iu][2][i];
      if (Q[iu][3][i] > Q_ext[i][1])
        Q_ext[i][1] = Q[iu][3][i];

      if (R[iu][0][i] < R_ext[i][0])
        R_ext[i][0] = R[iu][0][i];
      if (R[iu][1][i] < R_ext[i][0])
        R_ext[i][0] = R[iu][1][i];
      if (R[iu][2][i] < R_ext[i][0])
        R_ext[i][0] = R[iu][2][i];
      if (R[iu][3][i] < R_ext[i][0])
        R_ext[i][0] = R[iu][3][i];

      if (R[iu][0][i] > R_ext[i][1])
        R_ext[i][1] = R[iu][0][i];
      if (R[iu][1][i] > R_ext[i][1])
        R_ext[i][1] = R[iu][1][i];
      if (R[iu][2][i] > R_ext[i][1])
        R_ext[i][1] = R[iu][2][i];
      if (R[iu][3][i] > R_ext[i][1])
        R_ext[i][1] = R[iu][3][i];
    }
  }
}

void print_anisotropy(double patch[4][4][3])
{
  double hvec_nrm2 = 0, vvec_nrm2 = 0, tmp, anisotropy;
  int i;

  for (i = 0; i < 3; i++)
  {
    tmp = patch[1][0][i] + patch[2][0][i] - patch[1][3][i] - patch[2][3][i];
    hvec_nrm2 += tmp * tmp;
    tmp = patch[0][1][i] + patch[0][2][i] - patch[3][1][i] - patch[3][2][i];
    vvec_nrm2 += tmp * tmp;
  }
  anisotropy = sqrt(hvec_nrm2 / vvec_nrm2);
  printf("Anisotropy: %1.10f  [100]: %15.10f\n", anisotropy, patch[1][0][0]);
}

void Subdivide_patch(double patch[4][4][3], double ul_patch[4][4][3], double ur_patch[4][4][3], double dl_patch[4][4][3], double dr_patch[4][4][3])
{
  /*
    As we subdivide the patch, we try to force it to be more square by dividing in four pieces along the same direction if that direction is much longer than the other (based on the computed anisotropy).  Also, we no longer compute the bounds as this is slow.  They are computed only when needed.  Sometimes only the x direction bounds need to be computed (or x and y), rather than all three directions in order to prove that the patch does not intersect the ray.
  */

  double up_patch[4][4][3], down_patch[4][4][3];
  double hvec_nrm2 = 0, vvec_nrm2 = 0, tmp, anisotropy;
  int i;

  for (i = 0; i < 3; i++)
  {
    tmp = patch[1][0][i] + patch[2][0][i] - patch[1][3][i] - patch[2][3][i];
    hvec_nrm2 += tmp * tmp;
    tmp = patch[0][1][i] + patch[0][2][i] - patch[3][1][i] - patch[3][2][i];
    vvec_nrm2 += tmp * tmp;
  }
  anisotropy = sqrt(hvec_nrm2 / vvec_nrm2);

  if (anisotropy > 2)
  {
    hull_split_v(patch, up_patch, down_patch);
    hull_split_v(up_patch, ul_patch, ur_patch);
    hull_split_v(down_patch, dl_patch, dr_patch);
  }
  else if (anisotropy < 0.5)
  {
    hull_split_u(patch, up_patch, down_patch);
    hull_split_u(up_patch, ul_patch, ur_patch);
    hull_split_u(down_patch, dl_patch, dr_patch);
  }
  else
  {
    hull_split_u(patch, up_patch, down_patch);
    hull_split_v(up_patch, ul_patch, ur_patch);
    hull_split_v(down_patch, dl_patch, dr_patch);
  }
}

/*------------------------------------------------------------------------------------------*/
int Test_extents_surface(int model_num, float line_origin[3], float line_vector[3], float line_vector_inv[3])
/*------------------------------------------------------------------------------------------*/
/* This subroutine tests the bounding box of a surface to see if the projection ray will    */
/* intersect it										    */
/*------------------------------------------------------------------------------------------*/
{
  float xl, yl, zl;
  float xh, yh, zh;
  int i, j;
  float x0, y0, z0;
  float xd, yd, zd;
  float t1, t2, tnear = -100000, tfar = 100000;
  float temp;
  SURFACE *smodel;
  TRI_MODEL *tmodel;

  if (!use_tri_model)
  {
    smodel = (SURFACE *)&nrb_model[model_num];
    xl = smodel->min_x;
    xh = smodel->max_x;

    yl = smodel->min_y;
    yh = smodel->max_y;

    zl = smodel->min_z;
    zh = smodel->max_z;
  }
  else
  {
    tmodel = (TRI_MODEL *)&tri_model[model_num];
    xl = tmodel->min_x;
    xh = tmodel->max_x;

    yl = tmodel->min_y;
    yh = tmodel->max_y;

    zl = tmodel->min_z;
    zh = tmodel->max_z;
  }

  x0 = line_origin[0];
  y0 = line_origin[1];
  z0 = line_origin[2];
  xd = line_vector[0];
  yd = line_vector[1];
  zd = line_vector[2];

  /* X PLANES */
  if (xd == 0)
  {
    if (x0 < xl || x0 > xh)
      return 0;
  }
  else
  {
    t1 = (xl - x0) * line_vector_inv[0];
    t2 = (xh - x0) * line_vector_inv[0];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    if (tfar < 0)
      return 0;
  }

  /* Y PLANES */
  if (yd == 0)
  {
    if (y0 < yl || y0 > yh)
      return 0;
  }
  else
  {
    t1 = (yl - y0) * line_vector_inv[1];
    t2 = (yh - y0) * line_vector_inv[1];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    if (tfar < 0)
      return 0;
  }

  /* Z PLANES */
  if (zd == 0)
  {
    if (z0 < zl || z0 > zh)
      return 0;
  }
  else
  {
    t1 = (zl - z0) * line_vector_inv[2];
    t2 = (zh - z0) * line_vector_inv[2];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;

    if (tfar < 0)
      return 0;
  }

  return 1;
}

void vec_inv(float vec[3], float vec_inv[3])

{
  float rayLength;

  // normalize
  rayLength = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  vec[0] /= rayLength;
  vec[1] /= rayLength;
  vec[2] /= rayLength;

  vec_inv[0] = 0;
  vec_inv[1] = 0;
  vec_inv[2] = 0;
  if (vec[0] != 0)
    vec_inv[0] = 1 / vec[0];
  if (vec[1] != 0)
    vec_inv[1] = 1 / vec[1];
  if (vec[2] != 0)
    vec_inv[2] = 1 / vec[2];
}

void patch_in_matlab(double patch[4][4][3])

{
  int ii, jj, kk;
  FILE *file_desc;
  //  dbug(-1,"\n\nfigure(1);clf;\n");
  // file_desc = fpo;  // to write to a file, it needs to be opened before this.
  file_desc = stdout;
  fprintf(file_desc, "hold on;\n\n");
  for (ii = 0; ii < 4; ii++)
    for (jj = 0; jj < 4; jj++)
      for (kk = 0; kk < 3; kk++)
      {
        fprintf(file_desc, "patch(%d,%d,%d) = %1.6f;", ii + 1, jj + 1, kk + 1, patch[ii][jj][kk]);
        if (kk == 2)
          fprintf(file_desc, "\n");
      }
  fprintf(file_desc, "plot3(patch(:,:,1),patch(:,:,2),patch(:,:,3),'k');\n");
  fprintf(file_desc, "plot3(patch(:,:,1)',patch(:,:,2)',patch(:,:,3)','k');\n\n\n");
}

void patch_in_matlab_xf(double patch[4][4][3], double xform[4][3])

{
  double xform_inv[4][3];
  double patch_xf[4][4][3], mind, tmp, tmp2;
  double bbox[8][3];
  int i, j;
  // first invert the xform, then apply it and pass the result to the standard patch_in_matlab
  invert_xform(xform, xform_inv);
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      apply_xform(xform_inv, patch[i][j], patch_xf[i][j]);

  // only print it if its close to the point 12.379941 -6.546379  3.662493;
  mind = 100;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      tmp2 = 0;
      tmp = fabs(patch_xf[i][j][0] - 12.5);
      tmp2 += tmp * tmp;
      tmp = fabs(patch_xf[i][j][1] + 6.68);
      tmp2 += tmp * tmp;
      tmp = fabs(patch_xf[i][j][2] - 4);
      tmp2 += tmp * tmp;
      tmp2 = sqrt(tmp2);
      mind = MIN(mind, tmp2);
    }

  mind = 0; // tmp disable check

  if (mind < 1)
    patch_in_matlab(patch_xf);
}

/*-------------------------------------------------------------------------*/
int Test_patch(double patch[4][4][3], double tol, float line_vec[3], float *costheta)
/*-------------------------------------------------------------------------*/
/* This function tests a Bezier patch to see if it is flat and can be      */
/* approximated by a rectangle                                             */
/* The new flatness criteria takes into account the angle between the      */
/* the ray and the patch, since small errors can be magnified if this      */
/* angle is small                                                          */
/*-------------------------------------------------------------------------*/
{
  int i, j;
  double A, B, C, D;
  double denom;
  double t_error, hvec_nrm2 = 0, vvec_nrm2 = 0;
  double max_error, anisotropy, tmp;
  double line_vec_nrm;

  // First, we compute the normal and return 1 if it has small enough magnitude (indicating the patch is small enough)
  Plane_eqn(patch[0][0], patch[0][3], patch[3][0], patch[3][3], &A, &B, &C, &D);
  denom = sqrt(A * A + B * B + C * C);
  if (denom < 4 * tol * tol)
    return 1; // The patch is small enough

  // Next, we check to see if the flatness tolerance is met.

  // line_vec_nrm=0;
  // for(i=0;i<3;i++)
  //   line_vec_nrm += line_vec[i]*line_vec[i];
  line_vec_nrm = sqrt(DOT(line_vec, line_vec));
  costheta[0] = (A * line_vec[0] + B * line_vec[1] + C * line_vec[2]) / (line_vec_nrm * denom);
  // dbug(2,"costheta: %1.12f\r\n",costheta[0]);
  max_error = fabs(tol * costheta[0]);

  // Check points, returning 0 if there is one that is one out of bounds
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      t_error = fabs((A * patch[i][j][0] + B * patch[i][j][1] + C * patch[i][j][2] + D) / denom);
      if (t_error > max_error)
        return 0; // Still need to divide this patch up
    }

  // Lastly, we check the anisotropy and return 0 if it is extreme
  for (i = 0; i < 3; i++)
  {
    tmp = patch[1][0][i] + patch[2][0][i] - patch[1][3][i] - patch[2][3][i];
    hvec_nrm2 += tmp * tmp;
    tmp = patch[0][1][i] + patch[0][2][i] - patch[3][1][i] - patch[3][2][i];
    vvec_nrm2 += tmp * tmp;
  }
  anisotropy = sqrt(hvec_nrm2 / vvec_nrm2);
  // dbug(3,"aniso-------------%15.10f--------------------------------\n", anisotropy);
  if ((anisotropy > 4) || (anisotropy < 1 / 4))
    return 0;

  return 1; // If we make it here, all errors are small... the patch is good
}

int Test_extents2(double patch[4][4][3], float line_origin[3], float line_vector[3], float line_vector_inv[3])
{
  double xl, yl, zl;
  double xh, yh, zh;
  int i, j;
  double x0, y0, z0;
  double xd, yd, zd;
  double t1, t2, tnear = 0, tfar = 100000;
  double tmp;

  // Speed of this routine is critical.   The following optimizations have been made:
  // tnear is initialized to 0 instead of -100000 so that we don't have to check the sign of tfar all the time (this is done implicitly by comparing to tnear, which we are already doing since tnear is always at least 0.
  // We no longer need to sort t1 and t2.  Instead, we ensure that t2 is bigger by computing it based on the sign of line_vector_inv
  // We do not compute the bounds until (unless) they are needed.  The x bounds are computed and then tested.  Only if we have not returned do we then compute the ybounds, etc.
  // further speed improvements could probably be made by vectorizing some operations (min, max), but this has not been done, and I'm not sure if the added overhead would make it worth it

  x0 = line_origin[0];
  y0 = line_origin[1];
  z0 = line_origin[2];
  xd = line_vector[0];
  yd = line_vector[1];
  zd = line_vector[2];

  /* X PLANES */
  xl = 100000;
  xh = -100000;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      tmp = patch[i][j][0];

      if (tmp < xl)
        xl = tmp;
      if (tmp > xh)
        xh = tmp;
    }

  if (xd == 0)
  {
    if (x0 < xl || x0 > xh)
      return 0;
    //            {if (d_flag2) dbug(3,"        -----a----  low: %f   test: %f    high: %f \r\n",xl,x0,xh);return 0;}
  }
  else
  {
    if (line_vector_inv[0] > 0)
    {
      t1 = (xl - x0) * line_vector_inv[0];
      t2 = (xh - x0) * line_vector_inv[0];
    }
    else
    {
      t1 = (xh - x0) * line_vector_inv[0];
      t2 = (xl - x0) * line_vector_inv[0];
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
  }

  /* Y PLANES */
  yl = 100000;
  yh = -100000;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      tmp = patch[i][j][1];

      if (tmp < yl)
        yl = tmp;
      if (tmp > yh)
        yh = tmp;
    }

  if (yd == 0)
  {
    if (y0 < yl || y0 > yh)
      return 0;
    //      {if (d_flag2) dbug(3,"        -----d----   \r\n");return 0;}
  }
  else
  {
    if (line_vector_inv[1] > 0)
    {
      t1 = (yl - y0) * line_vector_inv[1];
      t2 = (yh - y0) * line_vector_inv[1];
    }
    else
    {
      t1 = (yh - y0) * line_vector_inv[1];
      t2 = (yl - y0) * line_vector_inv[1];
    }

    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    //      {if (d_flag2) dbug(3,"        -----e----   \r\n");return 0;}
    //      {if (d_flag2) dbug(3,"        -----f----   \r\n");return 0;}
  }

  /* Z PLANES */
  zl = 100000;
  zh = -100000;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      tmp = patch[i][j][2];

      if (tmp < zl)
        zl = tmp;
      if (tmp > zh)
        zh = tmp;
    }

  if (zd == 0)
  {
    if (z0 < zl || z0 > zh)
      return 0;
    //      {if (d_flag2) dbug(3,"        -----g----   \r\n");return 0;}
  }
  else
  {
    if (line_vector_inv[2] > 0)
    {
      t1 = (zl - z0) * line_vector_inv[2];
      t2 = (zh - z0) * line_vector_inv[2];
    }
    else
    {
      t1 = (zh - z0) * line_vector_inv[2];
      t2 = (zl - z0) * line_vector_inv[2];
    }

    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    //      {if (d_flag2) dbug(3,"        -----h----   \r\n");return 0;}
    //            {if (d_flag2) dbug(3,"        -----i----   \r\n");return 0;}
  }

  return 1;
}

int Test_extents(float xl, float xh, float yl, float yh, float zl, float zh, float line_origin[3], float line_vector[3], float line_vector_inv[3])
{
  int i, j;

  float x0, y0, z0;
  float xd, yd, zd;

  float t1, t2, tnear = 0, tfar = 100000;

  x0 = line_origin[0];
  y0 = line_origin[1];
  z0 = line_origin[2];
  xd = line_vector[0];
  yd = line_vector[1];
  zd = line_vector[2];

  /* X PLANES */
  if (xd == 0)
  {
    if (x0 < xl || x0 > xh)
      return 0;
  }
  else
  {
    if (line_vector_inv[0] > 0)
    {
      t1 = (xl - x0) * line_vector_inv[0];
      t2 = (xh - x0) * line_vector_inv[0];
    }
    else
    {
      t1 = (xh - x0) * line_vector_inv[0];
      t2 = (xl - x0) * line_vector_inv[0];
    }

    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
  }

  /* Y PLANES */
  if (yd == 0)
  {
    if (y0 < yl || y0 > yh)
      return 0;
  }
  else
  {
    if (line_vector_inv[1] > 0)
    {
      t1 = (yl - y0) * line_vector_inv[1];
      t2 = (yh - y0) * line_vector_inv[1];
    }
    else
    {
      t1 = (yh - y0) * line_vector_inv[1];
      t2 = (yl - y0) * line_vector_inv[1];
    }

    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
  }

  /* Z PLANES */
  if (zd == 0)
  {
    if (z0 < zl || z0 > zh)
      return 0;
  }
  else
  {
    if (line_vector_inv[2] > 0)
    {
      t1 = (zl - z0) * line_vector_inv[2];
      t2 = (zh - z0) * line_vector_inv[2];
    }
    else
    {
      t1 = (zh - z0) * line_vector_inv[2];
      t2 = (zl - z0) * line_vector_inv[2];
    }

    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
  }

  return 1;
}

/*------------------------------------------------------------------------------------------*/
int Test_extents_TriModel(TRI_MODEL *tmodel, float line_origin[3], float line_vector[3], float line_vector_inv[3])
/*------------------------------------------------------------------------------------------*/
/* This subroutine tests the bounding box of a triangle model to see if the projection ray will    */
/* intersect it										    */
/*------------------------------------------------------------------------------------------*/
{
  float xl, yl, zl;
  float xh, yh, zh;
  int i, j;
  float x0, y0, z0;
  float xd, yd, zd;
  float t1, t2, tnear = -100000, tfar = 100000;
  float temp;

  xl = tmodel->min_x;
  xh = tmodel->max_x;

  yl = tmodel->min_y;
  yh = tmodel->max_y;

  zl = tmodel->min_z;
  zh = tmodel->max_z;

  x0 = line_origin[0];
  y0 = line_origin[1];
  z0 = line_origin[2];
  xd = line_vector[0];
  yd = line_vector[1];
  zd = line_vector[2];

  /* X PLANES */
  if (xd == 0)
  {
    if (x0 < xl || x0 > xh)
      return 0;
  }
  else
  {
    t1 = (xl - x0) * line_vector_inv[0];
    t2 = (xh - x0) * line_vector_inv[0];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    if (tfar < 0)
      return 0;
  }

  /* Y PLANES */
  if (yd == 0)
  {
    if (y0 < yl || y0 > yh)
      return 0;
  }
  else
  {
    t1 = (yl - y0) * line_vector_inv[1];
    t2 = (yh - y0) * line_vector_inv[1];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;
    if (tfar < 0)
      return 0;
  }

  /* Z PLANES */
  if (zd == 0)
  {
    if (z0 < zl || z0 > zh)
      return 0;
  }
  else
  {
    t1 = (zl - z0) * line_vector_inv[2];
    t2 = (zh - z0) * line_vector_inv[2];

    if (t1 > t2)
    {
      temp = t1;
      t1 = t2;
      t2 = temp;
    }
    if (t1 > tnear)
      tnear = t1;
    if (t2 < tfar)
      tfar = t2;

    if (tnear > tfar)
      return 0;

    if (tfar < 0)
      return 0;
  }

  return 1;
}

void Check_difference(double xint, double *min_error, XP_ARRAY *int_points, int *which_one)
{
  int i;
  double tmp_error;

  *min_error = 100.0;

  //   dbug(4,"cd0\n");
  for (i = 0; i < int_points->length; i++)
  {
    tmp_error = fabs(xint - int_points->xp[i].x);
    //    dbug(3,"           %f       ",tmp_error);
    if (tmp_error < *min_error)
    {
      *min_error = tmp_error;
      *which_one = i;
    }
  }
  //   dbug(4,"cd1\n");
}

int Check_IntPoint(double p_ext[3][2], float px, float py, float pz, double my_tol)
{
  float minx, maxx, miny, maxy, minz, maxz;

  minx = (float)p_ext[0][0] - my_tol;
  miny = (float)p_ext[1][0] - my_tol;
  minz = (float)p_ext[2][0] - my_tol;
  maxx = (float)p_ext[0][1] + my_tol;
  maxy = (float)p_ext[1][1] + my_tol;
  maxz = (float)p_ext[2][1] + my_tol;

  if (px >= minx && px <= maxx && py >= miny && py <= maxy && pz >= minz && pz <= maxz)
    return 1;
  else
    return 0;
}

/*
void clear_stack()

{
  s_depth=0;
}

void push_stack(int x)

{
  s_stack[s_depth]=x;
  s_depth++;
}

int pop_stack()

{
  s_depth--;
  return s_stack[s_depth];
}

void dbug_stack()

{
  int i;
  dbug(3,"stack: ");
  for(i=0;i<s_depth;i++)
    dbug(3,"%d ",s_stack[i]);
  dbug(3,"\r\n");
}

void set_d_flag2()

{
  d_flag2 = 0;
  if ((s_depth>2)&&(s_stack[0]==1)&&(s_stack[1]==3)&&(s_stack[2]==3))
    d_flag2 = 1;
}
*/

void normalize(double v[3])

{
  double nrm;
  int i;
  nrm = sqrt(DOT(v, v));

  for (i = 0; i < 3; i++)
    v[i] /= nrm;
}

void orthogonalize_2nd(double v1[3], double v2[3], double out[3])

{
  int i;

  for (i = 0; i < 3; i++)
    out[i] = v2[i] - DOT(v1, v2) * v1[i];
  normalize(out);
}

void matrix_mult(double xform[3][3], double vin[3], double vout[3])

{
  int i;
  for (i = 0; i < 3; i++)
    vout[i] = DOT(vin, xform[i]);
}

int parallelepiped_intersect(double patch[4][4][3], float line_origin[3], float line_vector[3])

{
  // these parameters (axes, bounds) describe the parallelepiped... since the parallelepiped is convex and all the points in the patch are
  // inside the parallelepiped, the patch surface described by the patch must also be contained in the parallelepiped
  double axes[3][3], bounds[3][2];
  double patch_axes[2][3];
  double dot, xint[2], or_dot_ax[3], vec_dot_ax[3];
  double origin[3], vec[3];
  int i, j, k;

  // compute patch axes and vectors orthogonal to patch axes
  for (i = 0; i < 3; i++)
  {
    patch_axes[0][i] = patch[1][0][i] + patch[2][0][i] - patch[1][3][i] - patch[2][3][i];
    patch_axes[1][i] = patch[0][1][i] + patch[0][2][i] - patch[3][1][i] - patch[3][2][i];
    origin[i] = (double)line_origin[i];
    vec[i] = (double)line_vector[i];
  }
  normalize(patch_axes[0]);
  normalize(patch_axes[1]);
  orthogonalize_2nd(patch_axes[0], patch_axes[1], axes[1]);
  orthogonalize_2nd(patch_axes[1], patch_axes[0], axes[0]);

  // The 3rd axis is computed as the cross product of the first two
  cross_product(axes[0], axes[1], axes[2]);
  normalize(axes[2]);

  // orient the axes so that DOT(vec,axes) is not negative
  matrix_mult(axes, vec, vec_dot_ax);
  for (i = 0; i < 3; i++)
  {
    if (vec_dot_ax[i] < 0)
    {
      vec_dot_ax[i] = -vec_dot_ax[i];
      for (j = 0; j < 3; j++)
        axes[i][j] = -axes[i][j];
    }
  }
  matrix_mult(axes, origin, or_dot_ax);

  // Now we compute bounds
  for (i = 0; i < 3; i++)
  {
    bounds[i][0] = 1e10;
    bounds[i][1] = -1e10;
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
      {
        dot = patch[j][k][0] * axes[i][0] + patch[j][k][1] * axes[i][1] + patch[j][k][2] * axes[i][2];
        bounds[i][0] = MIN(bounds[i][0], dot);
        bounds[i][1] = MAX(bounds[i][1], dot);
      }
  }

  xint[0] = 0;
  xint[1] = 1e10;
  for (i = 0; i < 3; i++)
  {
    if (vec_dot_ax[i] < 1e-9) // line parallel to this set of planes... just need to determine if the origin is in-bounds, if not return 0
    {
      dbug(2, "PARALLEL!  %1.12f\n", vec_dot_ax[i]);
      if ((bounds[i][0] > or_dot_ax[i]) || (bounds[i][1] < or_dot_ax[i]))
        return 0;
    }
    else // line not parallel to this set of planes
    {
      xint[0] = MAX(xint[0], (bounds[i][0] - or_dot_ax[i]) / vec_dot_ax[i]);
      xint[1] = MIN(xint[1], (bounds[i][1] - or_dot_ax[i]) / vec_dot_ax[i]);
      dbug(2, "i: %d     xi[0]: %15.10f  xi[1]: %15.10f\n", i, xint[0], xint[1]);
      if (xint[0] > xint[1])
        return 0;
    }
  }
  dbug(2, "                                                                                           hit\n");
  return 1;
}

void Intersect_dum(double patch[4][4][3], int organ_id, XP_ARRAY *int_points, double tol, float line_origin[3], float line_vector[3], float line_vector_inv[3], double xform[4][3])
/* ------------------------------------------------------------------------------------ */
/* This routine was created for timing purposes only... it matched Intersect_bez, but   */
/* only included the main routines needed to recurse, not the stuff done at the leaves  */
/* ------------------------------------------------------------------------------------ */
{
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];
  double A, B, C, D;
  double numerator, denominator, denom;
  double xint, min_error, midpt[3];
  int count, count2;
  int i, j, k;
  int flag;
  float px, py, pz;
  float costheta;

  // dbug(3,"bez0             len:%d\n",int_points->length);
  if (Test_patch(patch, tol, line_vector, &costheta))
  {
    dbug(3, "test");
  }
  else
  {
    // dbug(3,"subd\n");

    Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);

    if (Test_extents2(ul_patch, line_origin, line_vector, line_vector_inv))
      Intersect_dum(ul_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(ul_patch, xform);
    if (Test_extents2(ur_patch, line_origin, line_vector, line_vector_inv))
      Intersect_dum(ur_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(ur_patch, xform);
    if (Test_extents2(dl_patch, line_origin, line_vector, line_vector_inv))
      Intersect_dum(dl_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(dl_patch, xform);
    if (Test_extents2(dr_patch, line_origin, line_vector, line_vector_inv))
      Intersect_dum(dr_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(dr_patch, xform);
    // dbug(3,"bez-done-----------------\n");
  }
}

void Intersect_bez(double patch[4][4][3], int organ_id, XP_ARRAY *int_points, double tol, float line_origin[3], float line_vector[3], float line_vector_inv[3], double xform[4][3])
/* ------------------------------------------------------------------------------------ */
/* This routine iteratively breaks down surface patches until it finds one that is flat */
/* and is intersected by the projection ray.  Once it finds a flat patch, it calculates */
/* the intersection point                                                               */
/* ------------------------------------------------------------------------------------ */
{
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];
  double A, B, C, D;
  double numerator, denominator, denom;
  double xint, min_error, midpt[3];
  int count, count2;
  int i, j, k, which_one;
  int flag;
  float px, py, pz;
  float costheta;

  // dbug(3,"bez0             len:%d\n",int_points->length);
  if (Test_patch(patch, tol, line_vector, &costheta))
  {
    if (debug_flag == 3)
      patch_in_matlab_xf(patch, xform);
    // dbug(3,"                                                                              FLAT  \r\n");
    Plane_eqn(patch[0][0], patch[0][3], patch[3][0], patch[3][3], &A, &B, &C, &D);

    // refine D to get the plane to fit the patch even more closely
    for (i = 0; i < 3; i++)
    {
      midpt[i] = 0;
      for (j = 0; j < 4; j++)
        for (k = 0; k < 4; k++)
          midpt[i] += patch[j][k][i];
      midpt[i] *= (1.0 / 16.0);
    }
    D = -A * midpt[0] - B * midpt[1] - C * midpt[2];

    numerator = A * line_origin[0] + B * line_origin[1] + C * line_origin[2] + D;
    denominator = -A * line_vector[0] - B * line_vector[1] - C * line_vector[2];

    if (denominator != 0)
      xint = numerator / denominator;
    else
      xint = 0;

    px = line_origin[0] + xint * line_vector[0];
    py = line_origin[1] + xint * line_vector[1];
    pz = line_origin[2] + xint * line_vector[2];
    // if (debug_flag >= 2) patch_in_matlab(patch);

    if (parallelepiped_intersect(patch, line_origin, line_vector))
    {
      // dbug(3,"                                                                              INSIDE  \r\n");
      if (int_points->length == 0 && xint > 0)
      {
        int_points->length = 1;
        int_points->xp[0].x = xint;
        int_points->xp[0].organ_id = organ_id;
        int_points->xp[0].costheta = costheta;
        // dbug(2,"               costheta: %1.12f\r\n",costheta);
        dbug(2, "                 numerator:   %1.14f   xint: %1.14f\r\n", numerator, xint);
      }
      else if (xint > 0)
      {
        // Here, we check to make sure that this new intersection is not essentially the same (within pad*2) as one that we already computed.
        //  If it is the same as one we already have, we don't add it to the list (don't want duplicates).
        Check_difference(xint, &min_error, int_points, &which_one);
        if (min_error >= pad * 2.0)
        {
          count = 0;
          flag = 1;
          while (count < int_points->length && flag)
            if (xint > int_points->xp[count].x)
              count++;
            else
              flag = 0;

          for (count2 = int_points->length; count2 > count; count2--)
          {
            if (count2 >= 150 || count2 == 0)
            {
              dbug(-1, "\nError... count2 has unexpected value: %i\n", count2);
              exit(1);
            }

            int_points->xp[count2].x = int_points->xp[count2 - 1].x;
            int_points->xp[count2].organ_id = int_points->xp[count2 - 1].organ_id;
            int_points->xp[count2].costheta = int_points->xp[count2 - 1].costheta;
          }

          dbug(2, "                 numerator:   %1.14f   xint: %1.14f\r\n", numerator, xint);
          int_points->xp[count].x = xint;
          int_points->xp[count].organ_id = organ_id;
          int_points->xp[count].costheta = costheta;
          int_points->length++;
        }
        else
        {
          // If we find another intersection nearby, we mark it as having a tangential intersection (even though it may not) as a means of ensuring that we don't rely on this intersection for determining whether it is going in or coming out of the object.
          int_points->xp[which_one].costheta = 0;
        }
      }
    }
  }
  else
  {
    // dbug(3,"subd\n");

    Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);

    if (Test_extents2(ul_patch, line_origin, line_vector, line_vector_inv))
      Intersect_bez(ul_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(ul_patch, xform);
    if (Test_extents2(ur_patch, line_origin, line_vector, line_vector_inv))
      Intersect_bez(ur_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(ur_patch, xform);
    if (Test_extents2(dl_patch, line_origin, line_vector, line_vector_inv))
      Intersect_bez(dl_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(dl_patch, xform);
    if (Test_extents2(dr_patch, line_origin, line_vector, line_vector_inv))
      Intersect_bez(dr_patch, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform);
    // else if (debug_flag==3) patch_in_matlab_xf(dr_patch, xform);
    // dbug(3,"bez-done-----------------\n");
  }
}

int isAwayFromPatch(double patch[4][4][3], float position[3], double my_tol)
{
  // When this function returns a 1, it means that the input point is at least "my_tol" away from all points on the patch
  // When this function returns a 0, it means that the input point is at most ~3*my_tol away from all points on the patch
  // Note that when the distance is between my_tol and ~3*my_tol, the output can be either 0 or 1!
  int i, tinyBox = 1;
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];
  double p_ext[3][2];

  find_bounds(patch, p_ext);

  if (!Check_IntPoint(p_ext, position[0], position[1], position[2], my_tol))
    return 1; // CASE 1: if its not in the patch bounding box (paded by my_tol), then we say it is definitely "AWAY"
  else
  {

    // check to see if the box is tiny
    for (i = 0; i < 3; i++)
    {
      if ((p_ext[i][1] - p_ext[i][0]) > my_tol)
      {
        tinyBox = 0;
        break;
      }
    }

    if (tinyBox)
      return 0; // CASE 2: if it is in the bounding box and the box is tiny, it has to be close to the surface
    else
    {
      // If neither CASE 1 or CASE 2, we recurse until we satisfy one of the two conditions.
      Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);

      return (isAwayFromPatch(ul_patch, position, my_tol) && isAwayFromPatch(ur_patch, position, my_tol) && isAwayFromPatch(dl_patch, position, my_tol) && isAwayFromPatch(dr_patch, position, my_tol));
    }
  }
}

int isAwayFromBez(bvh_element *treepointer, BEZIER_MODEL *bez_model, float position[3], double my_tol)
// When this function returns a 1, it means that the input point is at least "my_tol" away from all points on the surface of the bez_model.
// When this function returns a 0, it means that the input point is at most ~3*my_tol away from all points on the surface of the bez_model.
// Note that when the distance is between my_tol and ~3*my_tol, the output can be either 0 or 1!
{
  int i;
  double p_ext[3][2];

  p_ext[0][0] = treepointer->xmin;
  p_ext[0][1] = treepointer->xmax;
  p_ext[1][0] = treepointer->ymin;
  p_ext[1][1] = treepointer->ymax;
  p_ext[2][0] = treepointer->zmin;
  p_ext[2][1] = treepointer->zmax;

  if (Check_IntPoint(p_ext, position[0], position[1], position[2], my_tol))
  {
    if (treepointer->right != NULL)
      if (!isAwayFromBez(treepointer->right, bez_model, position, my_tol))
        return 0;
    if (treepointer->left != NULL)
      if (!isAwayFromBez(treepointer->left, bez_model, position, my_tol))
        return 0;
    if (treepointer->right == NULL && treepointer->left == NULL)
      for (i = 0; i < treepointer->num; i++)
      {
        if (!isAwayFromPatch(bez_model->patches[treepointer->patches[i]].cpoints, position, my_tol))
          return 0;
      }
  }
  return 1;
}

void Find_Intersections(BEZIER_MODEL *bez_model, int organ_id, float line_origin[3], float line_vector[3], float line_vector_inv[3], XP_ARRAY *int_points, double tol)
{

  // THIS FUNCTION IS OBSOLETE

  int i;
  double p_ext[3][2];

  for (i = 0; i < bez_model->num_patches; i++)
  {
    if (Test_extents(bez_model->patches[i].minx, bez_model->patches[i].maxx,
                     bez_model->patches[i].miny, bez_model->patches[i].maxy,
                     bez_model->patches[i].minz, bez_model->patches[i].maxz,
                     line_origin, line_vector, line_vector_inv))
    {
      find_bounds(bez_model->patches[i].cpoints, p_ext);
      dbug(2, "patch: %d  len: %d \r\n", i, int_points->length);
      // if (debug_flag >= 2) patch_in_matlab(bez_model->patches[i].cpoints);
      // Intersect_bez(bez_model->patches[i].cpoints, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, p_ext);
    }
  }
}

void Find_Intersections2(bvh_element *treepointer, BEZIER_MODEL *bez_model, int organ_id, float line_origin[3], float line_vector[3], float line_vector_inv[3], XP_ARRAY *int_points, double tol)
{
  int i, hitsBox, ilo, xform_patches, j, k;
  //  double patch_xf[4][4][3], xform[4][3];
  double xform2[4][3];
  double(*xform)[3];
  double(*patch_xf)[4][3];
  float line_vector_xf[3], line_vector_inv_xf[3], line_origin_xf[3];

  // dbug(2,"fin  %d\n", int_points->length);

  hitsBox = Test_extents(treepointer->xmin, treepointer->xmax, treepointer->ymin, treepointer->ymax,
                         treepointer->zmin, treepointer->zmax, line_origin, line_vector, line_vector_inv);

  if (hitsBox)
  {
    if (treepointer->right != NULL)
      Find_Intersections2(treepointer->right, bez_model, organ_id, line_origin, line_vector, line_vector_inv, int_points, tol);

    if (treepointer->left != NULL)
      Find_Intersections2(treepointer->left, bez_model, organ_id, line_origin, line_vector, line_vector_inv, int_points, tol);

    if (treepointer->right == NULL && treepointer->left == NULL)
      for (i = 0; i < treepointer->num; i++)
      {
        dbug(2, "patch: %d  len: %d \r\n", treepointer->patches[i], int_points->length);
        ilo = int_points->length;
        // need to xform into patch coordinates here!... can precompute xform, as well as xformed patch coordinates in xformed coords
        //  All that remains is to compute the xformed line_vector, line_vector_inv, and line_origin

        xform_patches = 1;
        if (xform_patches)
        {
          // if ((debug_flag==2)&&(treepointer->patches[i]==1692)) debug_flag = 3;

          // This line is now precomputed and the results stored for every patch to save some runtime
          // get_patch_xform(bez_model->patches[treepointer->patches[i]].cpoints, xform, patch_xf);  // this line takes ~15% now

          xform = bez_model->patches[treepointer->patches[i]].xform;
          patch_xf = bez_model->patches[treepointer->patches[i]].xpoints;

          apply_xform_f(xform, line_origin, line_origin_xf);
          apply_rotation_f(xform, line_vector, line_vector_xf);
          vec_inv(line_vector_xf, line_vector_inv_xf);
          // if (debug_flag==3) patch_in_matlab(bez_model->patches[treepointer->patches[i]].cpoints);
          // if (debug_flag==3) patch_in_matlab(patch_xf);
          // Intersect_dum(patch_xf, organ_id, int_points, tol, line_origin_xf, line_vector_xf, line_vector_inv_xf, xform);
          // int_points->length = ilo;
          Intersect_bez(patch_xf, organ_id, int_points, tol, line_origin_xf, line_vector_xf, line_vector_inv_xf, xform);
          // debug_flag = MIN(2,debug_flag);
        }
        else
        {
          // identity xform
          for (j = 0; j < 4; j++)
            for (k = 0; k < 3; k++)
              if (j == k)
                xform2[j][k] = 1;
              else
                xform2[j][k] = 0;
          Intersect_bez(bez_model->patches[treepointer->patches[i]].cpoints, organ_id, int_points, tol, line_origin, line_vector, line_vector_inv, xform2);
        }
        // if (debug_flag == 3) debug_flag = 2;
        if (int_points->length > ilo)
        {
          dbug(2, "patchhit: %d  int_points->length: %d\r\n", treepointer->patches[i], int_points->length);
          if (debug_flag >= 2)
            patch_in_matlab(bez_model->patches[treepointer->patches[i]].cpoints);
        }
      }
    // dbug(2,"fend1\n");
  }
  else
  {
    // dbug(2,"fend2\n");
    return;
  }
}

void print_xform(double xform[4][3])

{
  int i, j;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 3; j++)
      printf("xform[%d][%d]: %10.5f\n", i, j, xform[i][j]);
    printf("\n");
  }
}

int get_patch_xform_unit_test()

{
  double patch[4][4][3], xform[4][3], patch_xf[4][4][3];
  int i, j, k;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
    {
      patch[i][j][0] = i;
      patch[i][j][1] = j;
      patch[i][j][2] = 0.1 * (i - 1.5) * (i - 1.5);
    }
  get_patch_xform(patch, xform, patch_xf);
  print_xform(xform);
  return 0;
}

double findFarPoint(float origin[3], float vector[3], float x1, float x2, int organ_id, double my_tol)

{
  // we pick points between x1 and x2 (starting with the midpoint) until we find one that is far from the surface
  int cnt;
  float x, r, w;
  float vec[3], tol_decay = 0.85;
  double golden_ratio1;

  golden_ratio1 = (sqrt(5.0) - 1.0) / 2.0;
  x = (x1 + x2) / 2;
  r = 0.5;
  cnt = 0;
  while (cnt < 50) // usually will return on first iteration (very rarely will take more than 2 or 3)
  {
    vec[0] = origin[0] + x * vector[0];
    vec[1] = origin[1] + x * vector[1];
    vec[2] = origin[2] + x * vector[2];

    if (isAwayFromBez(treepointer_nrb[organ_id], &bez_model[organ_id], vec, my_tol))
      return x; // usually done 1st time here

    // otherwise, start testing various points along the segment (not too close to the ends)
    r += golden_ratio1;
    if (r > 1)
      r -= 1;
    dbug(2, "r: %1.10f\n", r);
    w = 0.1 + 0.8 * r / 10000.0;
    x = x1 + (x2 - x1) * w;

    cnt++;
    if (cnt > 10)
      my_tol = my_tol * tol_decay; // if we are getting desperate, we start to reduce the tolerance
  }
  return (x1 + x2) / 2;
}

int Segm_Inside_Object(float origin[3], float vector[3], float x1, float x2, int organ_id)

{
  /*
    This function tests to see whether a given segment is inside or outside the object.  If it determines it is inside, it returns a 1
    If it determines it is outside, it returns 0.  In very rare cases, the segment is so close to the boundary, it may not be able to tell.

    The main idea is to cast rays in many directions from a point on the segment that is as far as possible from the surface boundary
    until no intersection is either (i) tangential to the surface boundary or (ii) close to another intersection.
    The number of intersections is then counted and if it is odd, the segment is deemed to be inside the object surface.
    Unfortunately, there are some rare cases that cause this approach to fail (like saddle-shaped patches), so we repeat this procedure
    until we have 3 rays that all agree.  This takes a lot of time, but fortunately this function is only used for the less than 1% of segments
    that are on rays where the "insideness" cannot be determined through other (simpler) means.
  */

  float line_o[3], line_v[3], line_v_inv[3], x_far;
  XP_ARRAY check;
  int have_problem = 1, i, j, cnt, done = 0, inside_cnt = 0, outside_cnt = 0;
  int cnt_lim = 362; // This should not exceed the number of points in mdps!
  int quorum_cnt = 3;

  x_far = findFarPoint(origin, vector, x1, x2, organ_id, tol1);
  dbug(2, "x_far1: %1.10f\n", x_far);
  line_o[0] = origin[0] + vector[0] * x_far;
  line_o[1] = origin[1] + vector[1] * x_far;
  line_o[2] = origin[2] + vector[2] * x_far;

  // We cast rays in random directions: FIXME: make this deterministic (maximally dispersed points on a sphere)... perhaps precompute directions and store... also add comments to this function
  cnt = 0;
  while (!done)
  {
    //    random_unit_vector_f(line_v);  //the result has no 0 components, so the inverse is easy to compute
    for (i = 0; i < 3; i++)
      line_v[i] = mdps[cnt][i];
    cnt++;
    line_v_inv[0] = 1 / line_v[0];
    line_v_inv[1] = 1 / line_v[1];
    line_v_inv[2] = 1 / line_v[2];
    check.length = 0;
    Find_Intersections2(treepointer_nrb[organ_id], &bez_model[organ_id], organ_id, line_o, line_v, line_v_inv, &check, tol2);

    for (i = 0; i < 3; i++)
      dbug(1, "%1.10f\n", line_v[i]);
    have_problem = 0;
    for (i = 0; i < check.length; i++)
    {
      dbug(2, "x_far7: %1.10f\n", x_far);
      dbug(2, "check.xp[%d].x: %1.10f\n", i, check.xp[i].x);
      dbug(2, "check.xp[%d].costheta: %1.10f\n", i, check.xp[i].costheta);
      if (fabs(check.xp[i].costheta) < 0.01)
      {
        have_problem = 1; // either a tangential intersection or a pair of nearby intersections
        dbug(2, "***************             costheta :%1.12f  *****************", check.xp[i].costheta);
        break;
      }
    }
    if (!have_problem)
    {
      if (check.length % 2)
        inside_cnt++;
      else
        outside_cnt++;
      done = (abs(inside_cnt - outside_cnt) >= quorum_cnt);
    }
    if ((cnt > cnt_lim) && (!done))
    {
      if (inside_cnt > (outside_cnt + 1))
        return 1;
      else if ((inside_cnt + 1) < outside_cnt)
        return 0;
      else
      {
        // This really should almost never happen...
        // Even if it does happen, though, its probably not a big deal, since it would be for a small tangential segment
        if (fabs(x1 - x2) > 1.5)
          dbug(-1, "!!! NURB projector error: can't determine whether this segment is inside the surface... picking randomly !!!  (segment length = %1.6f)\n", x1 - x2);
        return rand() % 2;
      }
    }
  }
  //  if ((cnt > 5)&&(cnt < 100)) dbug(-1,"cnt: %d",cnt);
  return (inside_cnt > outside_cnt);
}

int Check_Y_Boundary(float x, float y, float z, int organ_id)
{

  // THIS FUNCTION IS OBSOLETE!

  int hits = 0, misses = 0;
  float line_o[3], line_v[3], line_v_inv[3], rayLength;
  XP_ARRAY check;
  int cnt;

  line_o[0] = x;
  line_o[1] = y;
  line_o[2] = z;

  dbug(2, "origin_x: %1.14f\r\n", line_o[0]);
  dbug(2, "origin_y: %1.14f\r\n", line_o[1]);
  dbug(2, "origin_z: %1.14f\r\n", line_o[2]);

  dbug(3, "cin\n");
  for (cnt = 0; cnt < 6; cnt++)
  {
    dbug(3, "cnti: %d\n", cnt);
    line_v[0] = 0.0;
    line_v[1] = 0.0;
    line_v[2] = 0.0;
    if (cnt <= 2)
      line_v[cnt % 3] = 1.0;
    else
      line_v[cnt % 3] = -1.0;
    check.length = 0;
    // if ((cnt==3)&&(debug_flag == 2)) debug_flag = 3;
    Find_Intersections2(treepointer_nrb[organ_id], &bez_model[organ_id], organ_id, line_o, line_v, line_v, &check, tol1);
    // if ((cnt==3)&&(debug_flag == 3)) debug_flag = 2;
    dbug(3, "cnto: %d\n", cnt);
    if (check.length != 0)
      hits++;
    else
    {
      misses++;
      dbug(2, "miss!             cnt: %d\r\n", cnt);
    }
    if (misses > 2)
      break;
  }

  dbug(2, "h/m: %d/%d\n", hits, misses);
  if (misses > 2)
    return 0;
  if (hits == 6)
    return 1;
  dbug(3, "c3\n");

  hits = 0;
  misses = 0;
  for (cnt = 0; cnt < 6; cnt++)
  {
    line_v[0] = 0.1;
    line_v[1] = 0.1;
    line_v[2] = 0.1;
    if (cnt <= 2)
      line_v[cnt % 3] = 1.0;
    else
      line_v[cnt % 3] = -1.0;
    vec_inv(line_v, line_v_inv);

    check.length = 0;
    Find_Intersections2(treepointer_nrb[organ_id], &bez_model[organ_id], organ_id, line_o, line_v, line_v_inv, &check, tol2);
    if (check.length != 0)
      hits++;
    else
      misses++;
    dbug(2, "cnt %d\n", cnt);
    if (misses)
      break;
  }

  dbug(2, "h/m  ------  : %d %d\n", hits, misses);
  return (hits == 6);
}

void Fill(XP_ARRAY int_points, int organ_id, float line_origin[3], float line_vector[3], SEG_ARRAY *lines)
{
  int i, inside;
  float mag;
  int pred[150], predValid = 1, ii;

  mag = sqrt(line_vector[0] * line_vector[0] + line_vector[1] * line_vector[1] + line_vector[2] * line_vector[2]);
  line_vector[0] /= mag;
  line_vector[1] /= mag;
  line_vector[2] /= mag;

  // dbug(2,"length: %d\r\n",int_points.length);
  for (i = 0; i < int_points.length; i++)
  {
    dbug(2, "****   int_points.xp[%d].x: %1.14f       costheta: %1.14f\r\n", i, int_points.xp[i].x, int_points.xp[i].costheta);
    if (fabs(int_points.xp[i].costheta) < 0.05)
      predValid = 0;
    else if (predValid)
      if (i == 0)
        pred[0] = 1; // in-to-out
      else
        pred[i] = 1 - pred[i - 1];
  }
  if (pred[int_points.length - 1] == 1)
    predValid = 0;

  // Based on some experiments here (to run this, set to if(1)), we have validity more than 99% of the time (~99.3%)
  // This means we can take our time figuring out the ones that are potentially invalid using a slow version of Segm_Inside_Object
  if (0)
  {
    if (predValid)
      pval_cnt += 1;
    else
      pinv_cnt += 1;
    if (pval_cnt % 100000 == 1)
      dbug(-1, "invalid_fraction: %f  pinv_cnt: %ld  pval_cnt: %ld\n", ((float)pinv_cnt) / ((float)pval_cnt), pinv_cnt, pval_cnt);
  }

  for (i = 0; i < int_points.length - 1; i++)
  {
    if (predValid)      // we have no tangent intersections or very close intersections
      inside = pred[i]; // When we deem the prediction to be valid, we can use the shortcut
    else
      inside = Segm_Inside_Object(line_origin, line_vector, int_points.xp[i].x, int_points.xp[i + 1].x, organ_id); // this is slow, but run rarely
    if (inside)
    {
      dbug(2, " x1         %f\r\n", int_points.xp[i].x);
      dbug(2, " x2         %f\r\n", int_points.xp[i + 1].x);
      dbug(2, " costheta1   %f\r\n", int_points.xp[i].costheta);
      dbug(2, " costheta2   %f\r\n", int_points.xp[i + 1].costheta);
      dbug(2, " organ_id   %d\r\n", organ_id);
      lines->sp[lines->length].x1 = int_points.xp[i].x;
      lines->sp[lines->length].x2 = int_points.xp[i + 1].x;
      lines->sp[lines->length].organ_id = organ_id;
      lines->length += 1;
    }
  }
}

void CROSS(double dest[3], double v1[3], double v2[3])
{
  dest[0] = v1[1] * v2[2] - v1[2] * v2[1];
  dest[1] = v1[2] * v2[0] - v1[0] * v2[2];
  dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void SUB(double dest[3], double v1[3], double v2[3])
{
  dest[0] = v1[0] - v2[0];
  dest[1] = v1[1] - v2[1];
  dest[2] = v1[2] - v2[2];
}

int intersect_triangle(double orig[3], double dir[3],
                       TRIANGLE T,
                       double *t, double *u, double *v)
{
  double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
  double det, inv_det;
  double vert0[3], vert1[3], vert2[3];

  //   dbug(3,"orig : %f %f %f\r\n",orig[0],orig[1],orig[2]);
  // dbug(3,"dir  : %f %f %f\r\n",dir[0],dir[1],dir[2]);
  // print_verts(&T,1);

  vert0[0] = T.vertex[0].x;
  vert0[1] = T.vertex[0].y;
  vert0[2] = T.vertex[0].z;
  vert1[0] = T.vertex[1].x;
  vert1[1] = T.vertex[1].y;
  vert1[2] = T.vertex[1].z;
  vert2[0] = T.vertex[2].x;
  vert2[1] = T.vertex[2].y;
  vert2[2] = T.vertex[2].z;

  if (0)
  {
    printf("\nv1 %f %f %f", vert0[0], vert0[1], vert0[2]);
    printf("\nv2 %f %f %f", vert1[0], vert1[1], vert1[2]);
    printf("\nv3 %f %f %f\n", vert2[0], vert2[1], vert2[2]);
  }

  /* find vectors for two edges sharing vert0 */
  SUB(edge1, vert1, vert0);
  SUB(edge2, vert2, vert0);

  /* begin calculating determinant - also used to calculate U parameter */
  CROSS(pvec, dir, edge2);

  /* if determinant is near zero, ray lies in plane of triangle */
  det = DOT(edge1, pvec);

  if (det > -EPSILON && det < EPSILON)
    return 0;
  inv_det = 1.0 / det;

  /* calculate distance from vert0 to ray origin */
  SUB(tvec, orig, vert0);

  /* calculate U parameter and test bounds */
  *u = DOT(tvec, pvec) * inv_det;
  if (*u < 0.0 || *u > 1.0)
    return 0;

  /* prepare to test V parameter */
  CROSS(qvec, tvec, edge1);

  /* calculate V parameter and test bounds */
  *v = DOT(dir, qvec) * inv_det;
  if (*v < 0.0 || *u + *v > 1.0)
    return 0;

  /* calculate t, ray intersects triangle */
  *t = DOT(edge2, qvec) * inv_det;

  return 1;
}

void Find_Intersections_tri(bvh_element *treepointer, TRI_MODEL tri_model, int organ_id, float line_origin[3], float line_vector[3], float line_vector_inv[3], XP_ARRAY *int_points)
{
  int i;
  double u, v, t;
  double orig[3], dir[3];
  int count, count2, flag;

  orig[0] = line_origin[0];
  orig[1] = line_origin[1];
  orig[2] = line_origin[2];

  dir[0] = line_vector[0];
  dir[1] = line_vector[1];
  dir[2] = line_vector[2];

  if (Test_extents(treepointer->xmin, treepointer->xmax, treepointer->ymin, treepointer->ymax, treepointer->zmin, treepointer->zmax, line_origin, line_vector, line_vector_inv))
  {
    if (treepointer->right != NULL)
      Find_Intersections_tri(treepointer->right, tri_model, organ_id, line_origin, line_vector, line_vector_inv, int_points);

    if (treepointer->left != NULL)
      Find_Intersections_tri(treepointer->left, tri_model, organ_id, line_origin, line_vector, line_vector_inv, int_points);

    if (treepointer->right == NULL && treepointer->left == NULL)
    {
      // dbug(4,"te2\r\n");
      for (i = 0; i < treepointer->num; i++)
      {
        //           dbug(4,"%d\r\n",i);
        if (intersect_triangle(orig, dir, tri_model.tris[treepointer->patches[i]], &t, &u, &v))
        {
          //          dbug(4,"                       ------------------  %d\r\n",i);
          if (t > 0 && int_points->length == 0)
          {
            int_points->length = 1;
            int_points->xp[0].x = t;
            int_points->xp[0].organ_id = organ_id;
          }
          else if (t > 0)
          {
            count = 0;
            flag = 1;
            while (count < int_points->length && flag)
              if (t > int_points->xp[count].x)
                count++;
              else
                flag = 0;

            for (count2 = int_points->length; count2 > count; count2--)
            {
              int_points->xp[count2].x = int_points->xp[count2 - 1].x;
              int_points->xp[count2].organ_id = int_points->xp[count2 - 1].organ_id;
            }

            int_points->xp[count].x = t;
            int_points->xp[count].organ_id = organ_id;
            int_points->length++;
          }
        }
      }
    }
  }
}

int Check_Y_Boundary_tri(bvh_element *treepointer, TRI_MODEL tmodel, float x, float y, float z, int organ_id)
{
  float line_o[3], line_v[3];
  XP_ARRAY check;

  line_o[0] = x;
  line_o[1] = y;
  line_o[2] = z;
  line_v[0] = 0;
  line_v[1] = 1.0;
  line_v[2] = 0.0;

  check.length = 0;
  Find_Intersections_tri(treepointer, tmodel, organ_id, line_o, line_v, line_v, &check);
  if (check.length % 2 == 0)
    return 0;
  else
    return 1;
}

void Fill_tri(XP_ARRAY int_points, TRI_MODEL tmodel, int organ_id, float line_origin[3], float line_vector[3], SEG_ARRAY *lines)
{
  float x, y, z;
  int diff;
  float edge[200][2];
  int edge_counter;
  int count = 0;
  int flag = 0;
  int i, j;
  float mag;

  mag = sqrt(line_vector[0] * line_vector[0] + line_vector[1] * line_vector[1] + line_vector[2] * line_vector[2]);
  line_vector[0] /= mag;
  line_vector[1] /= mag;
  line_vector[2] /= mag;

  // We first remove any segments that are super tiny by replacing the two bounding intersections by a single intersection at the center of the two very close intersections
  // In most cases, this does nothing, but its possible that a line hits exactly at the boundary of two triangles on both its entry and exit from a surface.
  // This is especially common with catvoxel since the rays are orthogonal to the z axis (parallel to voxel surfaces)
  float tol = 1e-4;
  for (i = 0; i < int_points.length - 1; i++)
  {
    if (int_points.xp[i + 1].x - int_points.xp[i].x < tol)
    {
      int_points.xp[i].x = (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0;
      for (j = i + 1; j < int_points.length - 1; j++)
        int_points.xp[j].x = int_points.xp[j + 1].x;
      int_points.length -= 1;
    }
  }

  if (0)
  {
    printf("len: %d\n", int_points.length);
    for (i = 0; i < int_points.length; i++)
    {
      x = line_origin[0] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[0];
      y = line_origin[1] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[1];
      z = line_origin[2] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[2];
      printf("(x,y,z): (%f, %f, %f)\n", x, y, z);
    }
  }

  if (int_points.length % 2 != 0)
  {
    for (i = 0; i < int_points.length - 1; i++)
    {
      x = line_origin[0] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[0];
      y = line_origin[1] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[1];
      z = line_origin[2] + (int_points.xp[i].x + int_points.xp[i + 1].x) / 2.0 * line_vector[2];

      flag = Check_Y_Boundary_tri(treepointer_tri[organ_id], tmodel, x, y, z, organ_id);
      if (flag)
      {
        lines->sp[lines->length].x1 = int_points.xp[i].x;
        lines->sp[lines->length].x2 = int_points.xp[i + 1].x;
        lines->sp[lines->length].organ_id = organ_id;
        lines->length += 1;
        flag = 0;
      }
    }
  }
  else
  {
    for (i = 0; i < int_points.length; i += 2)
    {
      lines->sp[lines->length].x1 = int_points.xp[i].x;
      lines->sp[lines->length].x2 = int_points.xp[i + 1].x;
      lines->sp[lines->length].organ_id = organ_id;
      lines->length += 1;
    }
  }
}

void other_initializations()

{
  // compute a set of directions in S2 that are fairly uniformly spaced.  We start with the verteces of an icosahedron and then fill
  // in the faces with a triangular grid
  int i, j, k, l, n;
  int faces[20][3];
  int edges[30][2];
  double disp[3], dist, tri[3][3], base_pt[3], wt;
  int order = 6; // the number of points for various orders is: 12, 32, 92, 162, 252, 362
  double golden_ratio2;
  double R[3][3], R1[3][3], R2[3][3], ang, c, s, mn = 10;

  golden_ratio2 = (sqrt(5.0) + 1.0) / 2.0;
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
    {
      mdps[j * 2 + k][0] = 0;
      mdps[j * 2 + k][1] = 1 - 2 * j;
      mdps[j * 2 + k][2] = (1 - 2 * k) * golden_ratio2;
    }
  for (j = 4; j < 12; j++)
  {
    mdps[j][0] = mdps[j % 4][(0 + j / 4) % 3];
    mdps[j][1] = mdps[j % 4][(1 + j / 4) % 3];
    mdps[j][2] = mdps[j % 4][(2 + j / 4) % 3];
  }
  n = 0;
  for (i = 0; i < 11; i++)
    for (j = i + 1; j < 12; j++)
    {
      for (l = 0; l < 3; l++)
        disp[l] = mdps[i][l] - mdps[j][l];
      dist = sqrt(DOT(disp, disp));
      if (dist < 2.1)
      {
        edges[n][0] = i;
        edges[n][1] = j;
        n++;
      }
    }
  dbug(1, "EDGES:                    %d                 ----------------\n", n);

  n = 0;
  for (i = 0; i < 10; i++)
    for (j = i + 1; j < 11; j++)
      for (k = j + 1; k < 12; k++)
      {
        for (l = 0; l < 3; l++)
          disp[l] = mdps[i][l] - mdps[j][l];
        dist = sqrt(DOT(disp, disp));
        if (dist < 2.1)
        {
          for (l = 0; l < 3; l++)
            disp[l] = mdps[j][l] - mdps[k][l];
          dist = sqrt(DOT(disp, disp));
          if (dist < 2.1)
          {
            for (l = 0; l < 3; l++)
              disp[l] = mdps[k][l] - mdps[i][l];
            dist = sqrt(DOT(disp, disp));
            if (dist < 2.1)
            {
              faces[n][0] = i;
              faces[n][1] = j;
              faces[n][2] = k;
              n++;
            }
          }
        }
      }
  dbug(1, "FACES:                    %d                 ----------------\n", n);
  k = 12;
  for (n = 0; n < 30; n++)
  {
    for (i = 0; i < 2; i++)
    {
      tri[i][0] = mdps[edges[n][i]][0];
      tri[i][1] = mdps[edges[n][i]][1];
      tri[i][2] = mdps[edges[n][i]][2];
    }
    for (i = 1; i < order; i++)
    {
      wt = ((double)i) / ((double)order);
      for (l = 0; l < 3; l++)
        mdps[k][l] = wt * tri[1][l] + (1 - wt) * tri[0][l];
      k++;
    }
  }
  for (n = 0; n < 20; n++)
  {
    for (i = 0; i < 3; i++)
    {
      tri[i][0] = mdps[faces[n][i]][0];
      tri[i][1] = mdps[faces[n][i]][1];
      tri[i][2] = mdps[faces[n][i]][2];
    }
    for (i = 2; i <= order - 1; i++)
      for (j = 1; j < i; j++)
      {
        wt = ((double)j) / ((double)i);
        for (l = 0; l < 3; l++)
          base_pt[l] = wt * tri[1][l] + (1 - wt) * tri[2][l];
        wt = ((double)(i)) / ((double)(order));
        for (l = 0; l < 3; l++)
          mdps[k][l] = wt * base_pt[l] + (1 - wt) * tri[0][l];
        k++;
      }
  }
  // Rotate by a small angle about all 3 axes so that there are no points that have 0 coordinates
  ang = 0.7446 / order;
  c = cos(ang);
  s = sin(ang);

  R1[0][0] = 1;
  R1[0][1] = 0;
  R1[0][2] = 0;
  R1[1][0] = 0;
  R1[1][1] = c;
  R1[1][2] = -s;
  R1[2][0] = 0;
  R1[2][1] = s;
  R1[2][2] = c;

  R2[0][0] = c;
  R2[0][1] = -s;
  R2[0][2] = 0;
  R2[1][0] = s;
  R2[1][1] = c;
  R2[1][2] = 0;
  R2[2][0] = 0;
  R2[2][1] = 0;
  R2[2][2] = 1;

  for (i = 0; i < 3; i++)
    matrix_mult(R1, R2[i], R[i]);

  R2[0][0] = c;
  R2[0][1] = 0;
  R2[0][2] = s;
  R2[1][0] = 0;
  R2[1][1] = 1;
  R2[1][2] = 0;
  R2[2][0] = -s;
  R2[2][1] = 0;
  R2[2][2] = c;

  for (i = 0; i < 3; i++)
    matrix_mult(R, R2[i], R1[i]);

  for (i = 0; i < 3; i++)
    dbug(1, "%f %f %f\n", R1[i][0], R1[i][1], R1[i][2]);
  for (i = 0; i < k; i++)
  {
    for (j = 0; j < 3; j++)
      base_pt[j] = mdps[i][j];
    matrix_mult(R1, base_pt, mdps[i]);
    normalize(mdps[i]);
    for (j = 0; j < 3; j++)
      if (mn > fabs(mdps[i][j]))
        mn = fabs(mdps[i][j]);
    dbug(2, "p(%d,1) = %f;p(%d,2) = %f;p(%d,3) = %f;\n", i + 1, mdps[i][0], i + 1, mdps[i][1], i + 1, mdps[i][2]);
  }
  dbug(1, "min_coord: %f;\n", mn);
}

DLLEXPORT void set_module_info_NCAT(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int UNUSED)

{
  modules_NCAT.Height = my_memcpyd(Height, modules_NCAT.Height, nModuleTypes * sizeof(double));
  modules_NCAT.Width = my_memcpyd(Width, modules_NCAT.Width, nModuleTypes * sizeof(double));
  modules_NCAT.Pix = my_memcpyi(Pix, modules_NCAT.Pix, nModuleTypes * sizeof(int));
  modules_NCAT.Coords = my_memcpyd(Coords, modules_NCAT.Coords, maxPix * 2 * nModuleTypes * sizeof(double));
  modules_NCAT.Sub = my_memcpyi(Sub, modules_NCAT.Sub, nModuleTypes * sizeof(int));
  dbug(4, "\n\n\nmodule info : %d %d \r\n\n\n\n", maxSubDets, nModuleTypes);
  modules_NCAT.Sampling = my_memcpyd(Sampling, modules_NCAT.Sampling, 2 * maxSubDets * nModuleTypes * sizeof(double));
  modules_NCAT.Weight = my_memcpyd(Weight, modules_NCAT.Weight, maxSubDets * nModuleTypes * sizeof(double));
  modules_NCAT.maxPixPerModule = maxPix;
  modules_NCAT.maxSubDets = maxSubDets;

  other_initializations();
}

DLLEXPORT void set_src_info_NCAT(double *sourceWeights, int nSubSources)

{
  modules_NCAT.sourceWeights = my_memcpyd(sourceWeights, modules_NCAT.sourceWeights, nSubSources * sizeof(double));
  modules_NCAT.nSubSources = nSubSources;
}

DLLEXPORT void set_tolerance_info_NCAT(double t1, double t2, double pd)

{
  // The tolerance should be much (say ~10x) smaller than the thinnest object feature you care about
  tol1 = t1;
  tol2 = t1;      // Set both tolerances the same for now
  pad = t1 * 3.0; // This value seems to eliminate the vast majority of double intersections (even *2 works, using *3 to be safe)
  dbug(1, "%f %f %f\n", tol1, tol2, pad);
  srand(1010);
}

DLLEXPORT void set_material_info_NCAT(int materialCount, int eBinCount, double *muTable) // JDP

{
  int i, j;

  n_energies = eBinCount;
  n_materials = materialCount;
  if (n_materials > MAXIMUM_NUMBER_OF_MATERIALS)
  {
    dbug(-1, " !!! Attempt to use too many materials with nurb projector !!!    exiting...\n");
    exit(1);
  }
  for (i = 0; i < n_energies; i++)
  {
    // For each energy there are mu's for each material
    for (j = 0; j < n_materials; j++)
      mu_table[j][i] = muTable[j + n_materials * i];
  }
}

DLLEXPORT void set_material_info_polygon(int materialCount, int eBinCount, double *muTable) // JDP

{
  int i, j;

  if (mu_table_tri != NULL)
    free_matrix(mu_table_tri, 0, n_materials, 0, n_energies);
  n_energies = eBinCount;
  n_materials = materialCount;
  mu_table_tri = matrix(0, n_materials, 0, n_energies);
  for (i = 0; i < n_energies; i++)
  {
    for (j = 0; j < n_materials; j++)
      mu_table_tri[j][i] = muTable[j + n_materials * i];
  }
}

void intersections_NCAT_all(SEG_ARRAY *segments, float *a, float *alpha, float *alpha_inv, double rayLength)

{
  int i;               // loop counter
  XP_ARRAY int_points; /*Array of intersection points for the projection*/
  SEG_ARRAY lines;     /*Array of line segments*/

  dbug(3, "a[0]: %f\r\na[1]: %f\r\na[2]: %f\r\n \n\r", a[0], a[1], a[2]);
  lines.length = 0;

  // bugfix for nCAT convexity:
  lines.sp[lines.length].x1 = -10000;
  lines.sp[lines.length].x2 = 10000;
  lines.sp[lines.length].organ_id = -1;
  lines.length += 1;

  for (i = 0; i < NUM_MODELS; i++) // Cycle through Models in the Phantom
  {
    // dbug(4,"model: %d\n\r",i);
    if (Test_extents_surface(i, a, alpha, alpha_inv))
    {
      int_points.length = 0;
      if (!use_tri_model)
      {
        /*---Find Intersection points for MODEL #i---*/
        Find_Intersections2(treepointer_nrb[i], &bez_model[i], i, a, alpha,
                            alpha_inv, &int_points, tol2);
        /*---Create line segment from intersection points---*/
        Fill(int_points, i, a, alpha, &lines);
      }
      else
      {
        /*---Find Intersection points for MODEL #i---*/
        Find_Intersections_tri(treepointer_tri[i], tri_model[i], i, a, alpha,
                               alpha_inv, &int_points);
        /*---Create line segment from intersection points---*/
        Fill_tri(int_points, tri_model[i], i, a, alpha, &lines);
      }
    }
  }
  // dbug(4,"ECHO \n\r");

  segments->length = 0;
  if (lines.length != 0)
  {
    // dbug(4,"--------------------------------------------------------------------------------\r\n");
    /*---Break up all line segments into organ paths-----*/
    /*
          for(i = 0; i < lines.length; i++)
            Break_segment(lines, lines.sp[i], segments, i);
    */
    Break_segment2(lines, segments); // Break_segment2 has now replaced Break_segment

    // Note: In Sept 2014, Break_segment2() was created to improve the speed of catsim.  It could probably also be used to improve
    //   the speed of catvoxel by replacing the use of Break_segment() here, but implementing/testing this idea was beyond the scope
    //   of the work done at that time.

    // Update: June 2019: I went ahead and switched this to Break_segment2... seems to work fine
  }
}

DLLEXPORT void make_vol_ncat(float *volume, int Nx, double xoff, double dx, int Ny, double yoff, double dy, int Nz, double zoff, double dz, int oversampling, int UNUSED_num_volumes, int material_volumes)

{
  int i, j, k, l, m, p, energyBinsToSkip = 0, volume_offset, muid, muidp;
  float a[3], alpha[3], alpha_inv[3], *matl_tabl, tmp;
  double rayLength, xbrd;
  SEG_ARRAY segments; /*Array of line segments*/

  other_initializations();
  matl_tabl = malloc(sizeof(float) * n_materials);
  for (i = 0; i < n_materials; i++)
  {
    if (material_volumes)
      matl_tabl[i] = 1.0 / (oversampling * oversampling);
    else if (use_tri_model)
      matl_tabl[i] = mu_table_tri[i][energyBinsToSkip] / (oversampling * oversampling);
    else
      matl_tabl[i] = mu_table[i][energyBinsToSkip] / (oversampling * oversampling);
    dbug(1, "\n\rmat%d:  %1.12lf\n\r", i, matl_tabl[i] * oversampling * oversampling);
  }
  dy = dy / oversampling;
  yoff = oversampling * (yoff - 1) + (oversampling + 1) * 0.5;
  dz = dz / oversampling;
  zoff = oversampling * (zoff - 1) + (oversampling + 1) * 0.5;

  alpha[0] = 1;
  alpha[1] = 0;
  alpha[2] = 0;
  alpha_inv[0] = 1;
  alpha_inv[1] = 0;
  alpha_inv[2] = 0;
  a[0] = -600;      // fixme: hardcoded value... needs to be lower than surface coordinate in x direction
  rayLength = 1200; // fixme: hardcoded value... needs to be big enough to get whole phantom
  for (l = 0; l < Nz; l++)
  {
    dbug(0, "Processing slice: %d     \r", l);
    for (k = l * oversampling; k < l * oversampling + oversampling; k++)
    {
      a[2] = (k + 1 - zoff) * dz;
      for (m = 0; m < Ny; m++)
      {
        for (j = m * oversampling; j < m * oversampling + oversampling; j++)
        {
          a[1] = (j + 1 - yoff) * dy;

          intersections_NCAT_all(&segments, a, alpha, alpha_inv, rayLength);

          // if (segments.length>1) debug_flag = 1; else debug_flag = 0;
          //	  if (debug_flag){
          //   for(i=0;i<segments.length;i++)
          //     {
          //		dbug(0,"segments.sp[%d]:  x1 = %f   x2 = %f  organ_id = %d\n", i, segments.sp[i].x1, segments.sp[i].x2, segments.sp[i].organ_id);
          //     }
          //  }

          if (segments.length > 0)
          {
            // add in a segment to cover the portion of the ray before the first intersection
            segments.sp[segments.length].x1 = -1e9;
            segments.sp[segments.length].x2 = -1e9;
            segments.sp[segments.length].organ_id = -1; // flag for air
            segments.length = segments.length + 1;
            // sort segments by end point
            qsort(segments.sp, segments.length, sizeof(LINE_SEG), comp_lines);
            // put in actual x2 value for added segment (which should now be at the beginning)
            segments.sp[0].x2 = segments.sp[1].x1;

            for (i = 0; i < segments.length; i++)
            {
              if (i < segments.length - 1)
                tmp = segments.sp[i].x2 - segments.sp[i + 1].x1;
              else
                tmp = 0;
              if (segments.sp[i].x1 > segments.sp[i].x2)
                dbug(2, "segment %3d -- min: %10.5f max: %10.5f diff: %10.5f\r\n", i, segments.sp[i].x1, segments.sp[i].x2, tmp);
            }

            // printf("%d %d %d\r\n",j,k,segments.length);
            // printf("%d",j);
            i = 0;
            p = 0;
            xbrd = (i + 0.5 - xoff) * dx - a[0];
            while ((segments.sp[p].x2 < xbrd) && (p < segments.length))
              p++;

            while ((i < Nx) && (p < segments.length))
            {
              if (segments.sp[p].organ_id == -1)
                muid = -1;
              else
              {
                if (use_tri_model)
                  muid = tri_model[segments.sp[p].organ_id].MU_ID;
                else
                  muid = nrb_model[segments.sp[p].organ_id].MU_ID;
              }
              if (segments.sp[p + 1].organ_id == -1)
                muidp = -1;
              else
              {
                if (p + 1 < segments.length)
                  if (use_tri_model)
                    muidp = tri_model[segments.sp[p + 1].organ_id].MU_ID;
                  else
                    muidp = nrb_model[segments.sp[p + 1].organ_id].MU_ID;
              }

              if (segments.sp[p].organ_id != -1)
                volume_offset = Nx * m + Nx * Ny * l + Nx * Ny * Nz * (muid - MATERIAL_INDEX_BASIS) * material_volumes;
              // dbug(1,"muid %d   voffset: %d\n",muid, volume_offset);

              while ((xbrd < segments.sp[p].x2) && (i < Nx))
              {
                if (segments.sp[p].organ_id != -1)
                  volume[i + volume_offset] += matl_tabl[muid - MATERIAL_INDEX_BASIS];
                // dbug(1,"%d ",muid-MATERIAL_INDEX_BASIS);
                // printf("%1.12lf %1.12lf %d\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l],i+Nx*m+Nx*Ny*l);
                xbrd += dx;
                i++;
              }
              // Remove partial volume contributions
              if (xbrd > segments.sp[p].x2)
              {
                i--;
                xbrd -= dx;
                if (segments.sp[p].organ_id != -1)
                  volume[i + volume_offset] -= matl_tabl[muid - MATERIAL_INDEX_BASIS] * (dx + xbrd - segments.sp[p].x2) / dx;
                // printf("B  %1.12lf %1.12lf\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l]);
                if ((p + 1 < segments.length) && (segments.sp[p + 1].organ_id != -1))
                {
                  volume_offset = Nx * m + Nx * Ny * l + Nx * Ny * Nz * (muidp - MATERIAL_INDEX_BASIS) * material_volumes;
                  volume[i + volume_offset] -= matl_tabl[muidp - MATERIAL_INDEX_BASIS] * (segments.sp[p].x2 - xbrd) / dx;
                }
                // printf("C  %1.12lf %1.12lf %d\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l],i+Nx*m+Nx*Ny*l);
              }
              p++;
            }
          }
        }
      }
    }
  }
  free(matl_tabl);
}

void intersections_NCAT(double *detCenter, double *right, double *up, double *sampling, int nSubDets, float *sourcePoints, double *thisView, int detIndex, double subviewWeight, double *detWeights)

{
  int i, j, k, l, m, sl1; // loop counters
  double subDetCenter[3];
  float alpha[3], alpha_inv[3];
  float line_int, line_int1, line_int2;
  double rayLength;
  XP_ARRAY int_points;       /*Array of intersection points for the projection*/
  SEG_ARRAY lines, segments; /*Array of line segments*/

  int myi;
  float dist;
  float atten;

  int skip_call_flag = 1; // set to 1 for speed.

  // dbug(4,"detIndex : %d\n\r",detIndex);
  for (l = 0; l < nSubDets; l++)
  {
    // dbug(1,"subdet: %d\n\r",l);
    for (i = 0; i < 3; i++)
      subDetCenter[i] = detCenter[i] + right[i] * sampling[2 * l] + up[i] * sampling[2 * l + 1];
    for (k = 0; k < modules_NCAT.nSubSources; k++)
    {
      dbug(2, "subsrc: %d\n\r", k);
      for (i = 0; i < 3; i++)
        alpha[i] = (float)subDetCenter[i] - sourcePoints[3 * k + i];
      rayLength = sqrt(alpha[0] * alpha[0] + alpha[1] * alpha[1] + alpha[2] * alpha[2]);
      alpha[0] /= rayLength;
      alpha[1] /= rayLength;
      alpha[2] /= rayLength;

      dbug(2, "sourcePoints[0]: %1.14f\r\n", sourcePoints[3 * k]);
      dbug(2, "sourcePoints[1]: %1.14f\r\n", sourcePoints[3 * k + 1]);
      dbug(2, "sourcePoints[2]: %1.14f\r\n", sourcePoints[3 * k + 2]);
      dbug(2, "alpha[0]: %1.14f\r\n", alpha[0]);
      dbug(2, "alpha[1]: %1.14f\r\n", alpha[1]);
      dbug(2, "alpha[2]: %1.14f\r\n", alpha[2]);

      alpha_inv[0] = 0;
      alpha_inv[1] = 0;
      alpha_inv[2] = 0;
      if (alpha[0] != 0)
        alpha_inv[0] = 1 / alpha[0];
      if (alpha[1] != 0)
        alpha_inv[1] = 1 / alpha[1];
      if (alpha[2] != 0)
        alpha_inv[2] = 1 / alpha[2];

      lines.length = 0;
      for (i = 0; i < NUM_MODELS; i++) // Cycle through Models in the Phantom
      {
        // if ((i==57) && (debug_flag == 1) && (k==0)) debug_flag = 2;
        // dbug(2,"model: %d, %d %d\n\r",i,l,k);
        if (Test_extents_surface(i, &sourcePoints[3 * k], alpha, alpha_inv))
        {
          int_points.length = 0;
          if (!use_tri_model)
          {
            /*---Find Intersection points for MODEL #i---*/
            // if ((i==169) && (debug_flag == 1) && (k==3) && (l==0)) debug_flag = 2;
            Find_Intersections2(treepointer_nrb[i], &bez_model[i], i, &sourcePoints[3 * k], alpha, alpha_inv, &int_points, tol2);
            // if ((i==169) && (debug_flag == 2)) debug_flag = 1;
            // if (int_points.length != 0){
            // dbug(2, "len: %3d   model: %3d\r\n", int_points.length, i);
            //	for(j=0;j<int_points.length;j++)
            //   dbug(2,"int_points.xp[%d]: %1.10f\r\n",j,int_points.xp[j]);
            // }
            /*---Create line segment from intersection points---*/
            Fill(int_points, i, &sourcePoints[3 * k], alpha, &lines);
            // if (int_points.length != 0) dbug(2, "line-len: %3d   model: %3d\r\n", lines.length, i);
          }
          else
          {
            // dbug(2,"start\r\n");
            Find_Intersections_tri(treepointer_tri[i], tri_model[i], i, &sourcePoints[3 * k], alpha,
                                   alpha_inv, &int_points);
            Fill_tri(int_points, tri_model[i], i, &sourcePoints[3 * k], alpha, &lines);
            // dbug(2,"end\r\n");
          }
        }
      }
      // if (debug_flag == 2) debug_flag = 1;

      // Now we do the final calculations for thisView
      if (lines.length != 0)
      {
        dbug(1, "llen: %d \n\r", lines.length);
        // dbug(4,"--------------------------------------------------------------------------------\r\n");
        /*---Break up all line segments into organ paths-----*/
        segments.length = 0;
        Break_segment2(lines, &segments); // Break_segment2 has now replaced Break_segment
        // for(i = 0; i < lines.length-1; i++)
        // Break_segment(lines, lines.sp[i], &segments, i);

        if (0)
        { // use this to sanity test Break_segment2() vs Break_segment()
          sl1 = segments.length - 1;
          Calc_line_int2(segments, mu_table, 0, &line_int1);

          segments.length = 0;
          for (i = 0; i < lines.length - 1; i++)
            Break_segment(lines, lines.sp[i], &segments, i);
          Calc_line_int(lines, segments, mu_table, 0, &line_int2);

          if (fabs(line_int1 - line_int2) > 0.001)
          {
            dbug(-1, "seg_len1: %d\n", sl1);
            dbug(-1, "seg_len2: %d  line_ints: %f %f\n\n", segments.length, line_int1, line_int2);
            exit(1);
          }
        }
        // dbug(4,"sd %d; ss %d\n\r",l,k);
        /*----Cycle through energies in the spectrum and calculate line integral for each material----*/
        for (m = 0; m < n_energies; m++)
        {
          if (use_tri_model)
          {

            if (skip_call_flag == 0)
              Calc_line_int2_tri(segments, mu_table_tri, m, &line_int); // This line is taking the bulk of the time...  weird!
            else
            {
              // we now simply skip the function call as it seems it was the slow part (by inlining the function)
              // void Calc_line_int2_tri(SEG_ARRAY segments, float **mu_table, int energy_ID, float *line_int)
              line_int = 0.0;
              for (myi = 0; myi < segments.length; myi++)
              {
                dist = segments.sp[myi].x2 - segments.sp[myi].x1;
                atten = mu_table_tri[tri_model[segments.sp[myi].organ_id].MU_ID][m] * tri_model[segments.sp[myi].organ_id].density;
                line_int += (atten * dist);
              }
            }
          }
          else
          {
            // Calc_line_int(lines, segments, mu_table, m, &line_int);   //now obsolete
            Calc_line_int2(segments, mu_table, m, &line_int);
          }
          thisView[detIndex * n_energies + m] += exp(-line_int) * subviewWeight * detWeights[l] * modules_NCAT.sourceWeights[k];
        }
      }
      else
      {
        /*----Cycle through energies in the spectrum and calculate line integral for each material----*/
        for (m = 0; m < n_energies; m++)
          thisView[detIndex * n_energies + m] += subviewWeight * detWeights[l] * modules_NCAT.sourceWeights[k];
      }
    }
  }
}

DLLEXPORT void ncat_projector(double subviewWeight, double *thisView, float *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength)

{
  int i, moduleNumber, k, nSubDets, detIndex;
  double *detWeights, *UV, *sampling, detCenter[3], *center, *right, *up;
  int moduleTypeIndex;

  dbug(2, "Number of modules: %d    \r\n", nModulesIn);
  for (moduleNumber = 0; moduleNumber < nModulesIn; moduleNumber++)
  {
    moduleTypeIndex = modTypeInds[moduleNumber];
    detWeights = &modules_NCAT.Weight[moduleTypeIndex * modules_NCAT.maxSubDets];
    UV = &modules_NCAT.Coords[2 * modules_NCAT.maxPixPerModule * moduleTypeIndex];
    sampling = &modules_NCAT.Sampling[2 * modules_NCAT.maxSubDets * moduleTypeIndex];
    nSubDets = modules_NCAT.Sub[moduleTypeIndex];
    center = &Center[3 * moduleNumber];
    right = &Right[3 * moduleNumber];
    up = &Up[3 * moduleNumber];

    dbug(3, "Number of pixels in this module: %d    \r\nmoduleTypeIndex : %d\n\r", modules_NCAT.Pix[moduleTypeIndex], moduleTypeIndex);
    for (k = 0; k < modules_NCAT.Pix[moduleTypeIndex]; k++)
    {

      //        if ((d_flag)&&((k==64+11)||(k==64+11))) debug_flag = 4; else debug_flag = 0;

      // if (k==73) debug_flag = 1;
      dbug(1, "Pixel:%d ", k);
      // if (k==74) debug_flag = 0;
      // if (k==345) exit(1);
      //  Compute pixel center locations
      detIndex = firstDetIndex[moduleNumber] + k;
      for (i = 0; i < 3; i++)
        detCenter[i] = center[i] + UV[k * 2] * right[i] + UV[k * 2 + 1] * up[i];
      // Compute relative intensities
      intersections_NCAT(detCenter, right, up, sampling, nSubDets, sourcePoints, thisView, detIndex, subviewWeight, detWeights);
    } // Loop over pixels (not subpixels)
  }
}

void *ncat_projector_wrapper(void *pointerIn)

{
  struct projector_args *projectorArgs;
  int module;

  projectorArgs = (struct projector_args *)pointerIn;

  while (nextModuleInQ < modulesInQ)
  {
    pthread_mutex_lock(&QLock);
    module = nextModuleInQ;
    nextModuleInQ++;
    pthread_mutex_unlock(&QLock);

    if (module < modulesInQ)
    {
      //          if (module == 3600) d_flag = 1; else d_flag = 0;
      dbug(1, "running module %d\r\n", module);
      ncat_projector(projectorArgs[0].subviewWeight, projectorArgs[0].thisView, projectorArgs[0].sourcePoints, projectorArgs[0].nSubSources, projectorArgs[0].srcHullPoints, projectorArgs[0].nSrcHullPoints, &(projectorArgs[0].firstDetIndex[module]), 1, &(projectorArgs[0].modTypeInds[module]), &(projectorArgs[0].Up[3 * module]), &(projectorArgs[0].Right[3 * module]), &(projectorArgs[0].Center[3 * module]), projectorArgs[0].UNUSED_tvLength);
      dbug(1, "  done w/module %d\r\n", module);
    }
  }
  return NULL; // pthread_exit(NULL);
}

DLLEXPORT void ncat_projector_threaded(double subviewWeight, double *thisView, float *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength, int numThreads, double UNUSED)

{
  int i;
  struct projector_args projectorArgs[1];

  thread_count = numThreads;
  projectorArgs[0].subviewWeight = subviewWeight;
  projectorArgs[0].thisView = thisView;
  projectorArgs[0].sourcePoints = sourcePoints;
  projectorArgs[0].nSubSources = nSubSources;
  projectorArgs[0].srcHullPoints = srcHullPoints;
  projectorArgs[0].nSrcHullPoints = nSrcHullPoints;
  projectorArgs[0].firstDetIndex = firstDetIndex;
  projectorArgs[0].nModulesIn = nModulesIn;
  projectorArgs[0].modTypeInds = modTypeInds;
  projectorArgs[0].Up = Up;
  projectorArgs[0].Right = Right;
  projectorArgs[0].Center = Center;
  projectorArgs[0].UNUSED_tvLength = UNUSED_tvLength;

  nextModuleInQ = 0;
  modulesInQ = nModulesIn;

  // Create the threads
  t_id = malloc(sizeof(pthread_t) * thread_count);
  for (i = 0; i < thread_count; i++)
  {
    pthread_create(&t_id[i], NULL, ncat_projector_wrapper, projectorArgs);
  }
  // Wait for them to complete
  for (i = 0; i < thread_count; i++)
  {
    pthread_join(t_id[i], NULL);
  }
  free(t_id);
}

void set_para_for_Polygon(int flag)
{
	if(flag)
	{
		use_tri_model = 1;
		NUM_MODELS = NUM_POLY;
	}
	else
	{
		use_tri_model = 0;
		NUM_MODELS = NUM_NRB;
	}
}

DLLEXPORT void polygon_projector_one_thread(double subviewWeight, double *thisView, float *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength)

{
  set_para_for_Polygon(1);
  ncat_projector(subviewWeight, thisView, sourcePoints, nSubSources, srcHullPoints, nSrcHullPoints, firstDetIndex, nModulesIn, modTypeInds, Up, Right, Center, UNUSED_tvLength);
  set_para_for_Polygon(0);
}

DLLEXPORT void polygon_projector(double subviewWeight, double *thisView, float *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength, int numThreads, double UNUSED)

{
  set_para_for_Polygon(1);
  ncat_projector_threaded(subviewWeight, thisView, sourcePoints, nSubSources, srcHullPoints, nSrcHullPoints, firstDetIndex, nModulesIn, modTypeInds, Up, Right, Center, UNUSED_tvLength, numThreads, UNUSED);
  set_para_for_Polygon(0);
}

DLLEXPORT void make_vol_ncat_polygon(float *volume, int Nx, double xoff, double dx, int Ny, double yoff, double dy, int Nz, double zoff, double dz, int oversampling, int UNUSED_num_volumes, int material_volumes)
{
  set_para_for_Polygon(1);
  make_vol_ncat(volume, Nx, xoff, dx, Ny, yoff, dy, Nz, zoff, dz, oversampling, UNUSED_num_volumes, material_volumes);
  set_para_for_Polygon(0);
}
