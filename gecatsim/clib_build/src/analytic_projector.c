// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define VERY_BIG 1e300
#define MATERIAL_INDEX_BASIS 1  //can be 1 or zero
#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
//#define LIMIT_LOAD
#endif

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
  double distance;
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
static int nextModuleInQ;
static int modulesInQ;
static pthread_mutex_t QLock = PTHREAD_MUTEX_INITIALIZER;
static int thread_count;
static pthread_t *t_id = NULL;
static int planB_warned = 0;
static int debug_flag = 0;

struct indexed_list
{
  double value;
  double originalIndex;
};

//////////////////////////////////////////////

static void dbug(int level,const char *s,...)

{
  if (debug_flag >= level)
    {
      va_list ap;
      va_start(ap,s);
      vprintf(s,ap);
      va_end(ap);
      fflush(stdout);
    }
}

double solve_cubic(double *a)
     
{
     double a12,Q,Q3,R,D,theta,srD,S,T;

     a[0] = a[0]*(1.0/3.0);
     a12 = a[0]*a[0];
     Q = (1.0/3.0)*a[1]-a12;
     Q3 = Q*Q*Q;
     R = 0.5*(a[0]*a[1]-a[2])-a12*a[0];
     D = Q3+R*R;
     if (D<0) {
       theta = acos(R/sqrt(-Q3)); 
       return 2*sqrt(-Q)*cos(theta*(1.0/3.0))-a[0];
     } else {
       srD = sqrt(D);
       S = cbrt(R+srD);
       T = cbrt(R-srD);
       return S+T-a[0];
     }

}

void cross(double *vec1, double *vec2, double *vecOut)

{
  vecOut[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  vecOut[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  vecOut[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

double magnitude(double r, double i)

{
  return sqrt(r*r+i*i);
}

void sqrtm(double in_r,double in_i,double *out_r,double *out_i)

{
  double ang,mag;
  ang = atan2(-in_i,-in_r)/2.0;
  mag = sqrt(magnitude(in_r,in_i));
  out_r[0] = mag*cos(ang);
  out_i[0] = mag*sin(ang);
}

void complex_multiply(double r1, double i1, double r2, double i2,double *real_out, double *imag_out)

{
  real_out[0] = r1*r2-i1*i2;
  imag_out[0] = r1*i2+r2*i1;
}

void solve_cubic_all(double *a,double *zr,double *zi)

{
  double a0,a1,a2,p,q,tmp,tmpr,tmpi,ur[3],ui[3],mag,ang,rmag2,pi = 3.141592653589793;
  int i;
  //Uses Cardano's method
  
  a2 = a[0];
  a1 = a[1];
  a0 = a[2];
  
  p = a1-a2*a2/3.0;
  q = a0+(2.0*a2*a2*a2-9.0*a2*a1)/27.0;
  
  tmp = q*q/4.0+p*p*p/27.0;
  //printf("P Q TMP: %1.12lf %1.12lf %1.12lf\n\r",p,q,tmp);
  if (q>0.0)
    if(tmp >= 0)
      {ur[0] = cbrt(q/2.0+sqrt(tmp));ui[0] = 0;mag = ur[0];ang = 0;}
    else
      {
	tmpr = q/2.0;tmpi = sqrt(-tmp);
	mag = cbrt(sqrt(tmpr*tmpr+tmpi*tmpi));
	ang = atan2(tmpi,tmpr)/3.0;
	ur[0] = mag*cos(ang);
	ui[0] = mag*sin(ang);
      }
  else
    if(tmp >= 0)
      {ur[0] = cbrt(q/2.0-sqrt(tmp));ui[0] = 0;mag = ur[0];ang = 0;}
    else
      {
	tmpr = q/2.0;tmpi = -sqrt(-tmp);
	mag = cbrt(magnitude(tmpr,tmpi));
	ang = atan2(tmpi,tmpr)/3.0;
	ur[0] = mag*cos(ang);
	ui[0] = mag*sin(ang);
      }
  ang += 2.0*pi/3.0;
  ur[1] = mag*cos(ang);
  ui[1] = mag*sin(ang);

  ang += 2.0*pi/3.0;
  ur[2] = mag*cos(ang);
  ui[2] = mag*sin(ang);

  //  for(i = 0;i<3;i++) printf("%1.12lf %1.12lf %1.12lf\r\n",ur[i],ui[i],magnitude(ui[i],ur[i]));
  
  rmag2 = 1.0/(mag*mag);

  for(i = 0;i<3;i++) {
    zr[i] = p/3.0*(ur[i]*rmag2)-ur[i]-a2/3;
    zi[i] = p/3.0*(-ui[i]*rmag2)-ui[i];
  }
  // for(i = 0;i<3;i++) printf("%1.12lf %1.12lf %1.12lf\r\n",zr[i],zi[i],magnitude(zi[i],zr[i]));
}

int compare_doubles (const void *a, const void *b)

{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}

int compIL (const void *a, const void *b)

{
  const struct indexed_list da = * (const struct indexed_list *) a;
  const struct indexed_list db = * (const struct indexed_list *) b;
  
  return (da.value > db.value) - (da.value < db.value);
}


DLLEXPORT int solve_quartic2(double *a,double *z)

{
// Faucette's method:
//   Faucette, W. M. "A Geometric Interpretation of the Solution of the General Quartic Polynomial." 
//   Amer. Math. Monthly 103, 51-57, 1996. 

  double yr[3],yi[3],a0,a1,a2,a3,p,q,r,cub[3],sr[3],si[3],zr[4],zi[4],errr,erri,err,errt = 0; 
  int num_real,i;
  
  a3 = a[0];
  a2 = a[1];
  a1 = a[2];
  a0 = a[3];
  
  //  First we solve the resolvent cubic:
  p = a2-0.375*a3*a3;
  q = a1-0.5*a2*a3+0.125*a3*a3*a3;
  r = a0-0.25*a1*a3+0.0625*a2*a3*a3-(3.0/256.0)*a3*a3*a3*a3;
  cub[0] = -2.0*p;
  cub[1] = p*p-4.0*r;
  cub[2] = q*q;

  //  for(i = 0;i<3;i++) printf("input: %1.12lf\n\r",cub[i]);

  solve_cubic_all(cub,yr,yi);

  //  for(i = 0;i<3;i++) printf("Y%d :  %1.12lf +i %1.12lf\n\r",i,yr[i],yi[i]);

  for(i = 0;i<3;i++) sqrtm(yr[i],yi[i],&sr[i],&si[i]);

  //for(i = 0;i<3;i++) printf("S%d :  %1.12lf +i %1.12lf\n\r",i,sr[i],si[i]);
  
  zr[0] = (sr[0]+sr[1]-sr[2])*0.5;
  zr[1] = (sr[2]+sr[0]-sr[1])*0.5;
  zr[2] = (sr[2]+sr[1]-sr[0])*0.5;
  zr[3] = -zr[0]-zr[1]-zr[2];

  zi[0] = (si[0]+si[1]-si[2])*0.5;
  zi[1] = (si[2]+si[0]-si[1])*0.5;
  zi[2] = (si[2]+si[1]-si[0])*0.5;
  zi[3] = -zi[0]-zi[1]-zi[2];
  
  for(i = 0;i<4;i++){
    complex_multiply(zr[i],zi[i],zr[i],zi[i],&errr,&erri);
    errr += p;
    complex_multiply(errr,erri,zr[i],zi[i],&errr,&erri);
    errr += q;
    complex_multiply(errr,erri,zr[i],zi[i],&errr,&erri);
    errr += r;
    errr = magnitude(errr,erri);
    errt += errr*errr;
  }
  errt = sqrt(errt);
  
  if (errt>1e-5){
    zr[0] += sr[2];
    zr[1] -= sr[2];
    zr[2] -= sr[2];
    zr[3] += sr[2];
    zi[0] += si[2];
    zi[1] -= si[2];
    zi[2] -= si[2];
    zi[3] += si[2];
  }
  for(i = 0;i<4;i++) zr[i] -= a3/4;
  num_real = 0;
  for(i = 0;i<4;i++) if ((fabs(zi[i])/fabs(zr[i]))<1e-10) {z[num_real] = zr[i];num_real++;}
  qsort(z,num_real,sizeof(double),compare_doubles);
  return num_real;
}

//function z = solve_quartic(a)
//
// -----------------------------------------------------------------------
//
// Routine
//   solve_quartic : Computes the roots of a quartic equation.
//
// Author
//   Jed Pack (GE Global Research)
//   packj@research.ge.com
//
// Inputs
//   a : Vector containing [a_4, a_3, a_2, a_1, a_0] where the equation to 
//         be solved is:
//
//        a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
//
//     Note: If only four parameters are given, it is assumed a_4 is 1.
//
//
// Outputs
//   z : The roots of the equation 
// -----------------------------------------------------------------------


DLLEXPORT int solve_quartic(double *a,double *z)

{
  // Equations found in: 
  //   Weisstein, Eric W. "Quartic Equation." From MathWorld--A Wolfram Web Resource. 
  //   http://mathworld.wolfram.com/QuarticEquation.html 
  
  double tmp,y1,tmp2,R,alpha,beta,D,E,Dt,Et,cub[3];
  int i;
  
  tmp = a[0];
  for(i = 0;i<4;i++) a[i] = a[i+1]/tmp;
  
  tmp = 0.25*a[0]*a[0]-a[1];
  cub[0] = -a[1];
  cub[1] = a[2]*a[0]-4.0*a[3];
  cub[2] = 4.0*a[1]*a[3]-a[2]*a[2]-a[0]*a[0]*a[3];
  y1 = solve_cubic(cub);

  
  //  Next we compute R, D, and E:
  tmp2 = tmp+y1;
  //printf("y1  = %1.12lf\n\rtmp = %1.12lf     %1.12lf\n\r",y1,tmp,fabs(tmp2)/MAX(1.0,fabs(tmp)));
  if (fabs(tmp2)/MAX(1.0,fabs(tmp))>1e-4){
    if (tmp2>0){
      R = sqrt(tmp2);
      alpha = 0.75*a[0]*a[0]-R*R-2.0*a[1];
      beta = 0.25*(4.0*a[0]*a[1]-8.0*a[2]-a[0]*a[0]*a[0])/R;
      Dt = alpha+beta;
      Et = alpha-beta;
      //  printf("Dt = %1.12f\n\rEt = %1.12f\n\r",Dt,Et);

      if (Dt>0){
	D = sqrt(Dt);
	z[1] = -0.25*a[0]+0.5*R+0.5*D;
	z[0] = -0.25*a[0]+0.5*R-0.5*D;
	//printf("                          E = %1.12lf %1.12lf %1.12lf\n\r",D,R,z[0]);
	if (Et>0){
	  E = sqrt(Et);
	  z[3] = -0.25*a[0]-0.5*R+0.5*E;
	  z[2] = -0.25*a[0]-0.5*R-0.5*E;
	  if (z[2]<z[1]){
	    tmp = z[1];z[1] = z[2];z[2] = tmp;
	    if (z[1]<z[0]){tmp = z[0];z[0] = z[1];z[1] = tmp;}
	    if (z[3]<z[2]){tmp = z[2];z[2] = z[3];z[3] = tmp;if (z[2]<z[1]){tmp = z[1];z[1] = z[2];z[2] = tmp;}}
	  }
	  return 4;
	} else return 2;
      } else{
	if (Et>0){
	  E = sqrt(Et);
	  z[1] = -0.25*a[0]-0.5*R+0.5*E;
	  z[0] = -0.25*a[0]-0.5*R-0.5*E;
	  return 2;
	} else return 0;
      }
    } 
    else return 0;
  }
  else {
    //for(i = 0;i<4;i++) printf("input: %1.12lf\n\r",a[i]);
    return solve_quartic2(a,z);
    //    return -1;
  }
}

double quadratic_form(double *vec1,double *matrix, double *vec2)
{
  int i;
  double x = 0;
  for(i = 0;i<3;i++){x += (matrix[i]*vec2[0]+matrix[i+3]*vec2[1]+matrix[i+6]*vec2[2])*vec1[i];}
  return x;
}

/*
int quartic_intersect_C(double *a0,double *Qlp,double *Qrp,double shp,double *alpha,double *tc,int obj)
{
  int i,out;
  double scale,displ,a0t[3],tmp[3],aa,b,c,d,e,f,C[5],*Ql,*Qr,sh;
  
  Ql = &phantom.Qmatrix[obj*18];
  //  printf("%1.12lf %1.12lf",phantom.Qmatrix[obj);
  //  for(i = 1;i<9;i++) printf("L %1.12lf %1.12lf diff: %1.12lf\r\n",Ql[i],Qlp[i],Ql[i]-Qlp[i]);
  Qr = &phantom.Qmatrix[obj*18+9];
  //  for(i = 1;i<9;i++) printf("R %1.12lf %1.12lf diff: %1.12lf\r\n",Qr[i],Qrp[i],Qr[i]-Qrp[i]);
  sh = phantom.shape[obj];
  //  printf("s diff: %1.12lf\r\n",sh-shp);
  
  scale = 1.0;
  displ = 0;
  for(i = 0;i<3;i++) displ += a0[i]*a0[i];
  displ = sqrt(displ);                  // scaling/displacement done only f0r numerical reasons

  for(i = 0;i<3;i++) a0t[i] = (a0[i]+displ*alpha[i])/scale;

  // The quartic formula is of the form    A t^4 + B t^3 + C t^2 + D t + E = 0
  aa = quadratic_form(alpha,Ql,alpha);
  b = 2.0*quadratic_form(a0t,Ql,alpha);
  //c = 1-sh*sh+quadratic_form(a0t,Ql,a0t);
  c = sh+quadratic_form(a0t,Ql,a0t);
  d = 4.0*quadratic_form(alpha,Qr,alpha);
  e = 8.0*quadratic_form(a0t,Qr,alpha);
  f = 4.0*quadratic_form(a0t,Qr,a0t);

  C[0] = aa*aa;
  C[1] = 2*aa*b;
  C[2] = b*b+2*aa*c-d;
  C[3] = 2*b*c-e;
  C[4] = c*c-f;

  out = solve_quartic(C,tc);
  for(i = 0;i<out;i++) tc[i] = tc[i]*scale+displ;
  return out;
}
*/

int quartic_intersect(double *a0,double *alpha,double *tc,int obj)

{
  int i,out;
  double scale,displ,a0t[3],tmp[3],aa,b,c,d,e,f,C[5],*Ql,*Qr,sh;
  
  Ql = &phantom.Qmatrix[obj*18];
  Qr = &phantom.Qmatrix[obj*18+9];
  sh = phantom.shape[obj];
  scale = 1.0;
  displ = 0;
  for(i = 0;i<3;i++) displ += a0[i]*a0[i];
  displ = sqrt(displ);                  // scaling/displacement done only f0r numerical reasons

  for(i = 0;i<3;i++) a0t[i] = (a0[i]+displ*alpha[i])/scale;

  // The quartic formula is of the form    A t^4 + B t^3 + C t^2 + D t + E = 0
  aa = quadratic_form(alpha,Ql,alpha);
  b = 2.0*quadratic_form(a0t,Ql,alpha);
  //c = 1-sh*sh+quadratic_form(a0t,Ql,a0t);
  c = sh+quadratic_form(a0t,Ql,a0t);
  d = 4.0*quadratic_form(alpha,Qr,alpha);
  e = 8.0*quadratic_form(a0t,Qr,alpha);
  f = 4.0*quadratic_form(a0t,Qr,a0t);

  C[0] = aa*aa;
  C[1] = 2*aa*b;
  C[2] = b*b+2*aa*c-d;
  C[3] = 2*b*c-e;
  C[4] = c*c-f;

  out = solve_quartic(C,tc);
  for(i = 0;i<out;i++) tc[i] = tc[i]*scale+displ;
  return out;
}

DLLEXPORT int clip_all(double *a,double *alpha,double rayLength,double *tc2,int out,double *st_list,double *en_list,double *den_list,int *pri_list,int *mat_list,int num_int,int i)

{
  double b[3], s1, s2, tcrit, tmin, tmax, *eta, *s, den;
  int j, cp, firstPlane, materialIndex;

  firstPlane = phantom.clipStartIndex[i];
  eta = &phantom.clipNormalVector[firstPlane*3];
  s = &phantom.clipDistance[firstPlane];
  cp = phantom.nClipPlanes[i];
  den = phantom.density[i];
  materialIndex = phantom.materialInd[i];
  tmin = -VERY_BIG;tmax = VERY_BIG;
  for(j = 0;j<3;j++) b[j] = a[j]+rayLength*alpha[j];
  
  for(j = 0;j<cp;j++){
    s1 = (eta[j*3]*a[0]+eta[j*3+1]*a[1]+eta[j*3+2]*a[2]-s[j]);
    s2 = (eta[j*3]*b[0]+eta[j*3+1]*b[1]+eta[j*3+2]*b[2]-s[j]);
    if(s1*s2<0){
      tcrit = fabs(s1)/(fabs(s1)+fabs(s2))*rayLength;
      if (s1<s2) 
	{if (tcrit<tmax) tmax = tcrit;}
      else 
	{if (tcrit>tmin) tmin = tcrit;}
    }
    else 
    if ((s1+s2)>0){
      tmin = 0;tmax = 0;
      break;
    }
  }
  
  if(out>2) {
    if(tmax<tc2[2])
      out = 2;
    else 
      if(tmin>tc2[1]) {
	tc2[0] = tc2[2];
	tc2[1] = tc2[3];
      }
      else{
	tc2[3] = MIN(tmax,tc2[3]);
	st_list[num_int] = tc2[2];
	en_list[num_int] = tc2[3];
	den_list[num_int] = den;
	pri_list[num_int] = i;
	mat_list[num_int] = materialIndex;
	num_int = num_int+1;
      }
  }
  tc2[1] = MIN(tmax,tc2[1]);
  tc2[0] = MAX(tmin,tc2[0]);
  if (tc2[1]>tc2[0]){
    st_list[num_int] = tc2[0];
    en_list[num_int] = tc2[1];
    den_list[num_int] = den;
    pri_list[num_int] = i;
    mat_list[num_int] = materialIndex;
    num_int = num_int+1;
  }
  return num_int;
}

DLLEXPORT int quadratic_intersect(double *a0,double *alpha,int pars11,double *tc2,int obj)

{
  int out;
  double A,B,C,tmp,*Q,k;
    // The quadratic formula is of the form    A t^2 + B t + C = 0
  Q = &phantom.Qmatrix[obj*18];
  k = phantom.shape[obj];
  C = quadratic_form(a0,Q,a0)-k;
  B = 2.0*quadratic_form(alpha,Q,a0);
  A = quadratic_form(alpha,Q,alpha);
  if ((B*B)>(4.0*A*C)) {                    // If sign of determinant is positive there is an intersection with the
                                            //  non-clipped quadratic
    if(A >= 0.0){
      tmp = sqrt(B*B-4.0*A*C);
      tc2[0] = (-B-tmp)/(2.0*A);
      tc2[1] = (-B+tmp)/(2.0*A);
      out = 2;
    } else {                       // For cones and hyperboloids the segment bounded by the roots may be the part
                                   //  outside the object instead of the part inside the object (i.e., when A<0)
      tmp = sqrt(B*B-4.0*A*C);
      tc2[0] = -VERY_BIG;
      tc2[1] = (-B+tmp)/(2.0*A);
      tc2[2] = (-B-tmp)/(2.0*A);
      tc2[3] = VERY_BIG;
      out = 4;
    }
  }
  else 
    if((pars11 == 5) && (C<0)){
      tc2[0] = -VERY_BIG;                // For hyperboloid of one sheet the entire ray may be inside the quadratic surface
      tc2[1] = VERY_BIG;
      out = 2;
    }
    else out = 0;
  return out;
}

int* my_memcpyi(int *src, int *dest, int bytes)

{
  if (dest != NULL) 
    {
      free(dest);
      dest = NULL;
    }
  dest = (int *) malloc(bytes);
  memcpy(dest, src, bytes);
  return dest;
}

double* my_memcpyd(double *src, double *dest, int bytes)

{
  if (dest != NULL) 
    {
      free(dest);
      dest = NULL;
    }
  dest = (double *) malloc(bytes);
  memcpy(dest, src, bytes);
  return dest;
}

DLLEXPORT void set_src_info(double *sourceWeights, int nSubSources)

{
  modules.sourceWeights = my_memcpyd(sourceWeights, modules.sourceWeights, nSubSources*sizeof(double));
  modules.nSubSources = nSubSources;
}

DLLEXPORT void set_module_info(double *Height, double *Width, int *Pix, double *Coords, int *Sub, double *Sampling, double *Weight, int nModuleTypes, int maxPix, int maxSubDets, int moduleOverlapType)

{
  int i;
  modules.Height = my_memcpyd(Height, modules.Height, nModuleTypes*sizeof(double));
  modules.Width = my_memcpyd(Width, modules.Width, nModuleTypes*sizeof(double));
  for(i = 0; i < nModuleTypes; i++)
    {
      modules.Height[i] = MAX(1e-7,modules.Height[i]);
      modules.Width[i] = MAX(1e-7,modules.Width[i]);
    }
  modules.Pix = my_memcpyi(Pix, modules.Pix, nModuleTypes*sizeof(int));
  modules.Coords = my_memcpyd(Coords, modules.Coords, maxPix*2*nModuleTypes*sizeof(double));
  modules.Sub = my_memcpyi(Sub, modules.Sub, nModuleTypes*sizeof(int));
  modules.Sampling = my_memcpyd(Sampling, modules.Sampling, 2*maxSubDets*nModuleTypes*sizeof(double));
  modules.Weight = my_memcpyd(Weight, modules.Weight, maxSubDets*nModuleTypes*sizeof(double));
  modules.maxPixPerModule = maxPix;
  modules.maxSubDets = maxSubDets;
  modules.moduleOverlapType = moduleOverlapType;
}

DLLEXPORT void set_bounding_info(int numObjs, int *vertexStartInd, double *vertLocs, int numVerts)

{
  int i;

  bounding.vertexStartIndex = my_memcpyi(vertexStartInd, bounding.vertexStartIndex, (numObjs+1)*sizeof(int));
  bounding.vertexLocations = my_memcpyd(vertLocs, bounding.vertexLocations, sizeof(double)*numVerts*3);
  phantom.xbounds[0]=VERY_BIG;
  phantom.xbounds[1]=-VERY_BIG;
  for(i=0;i<numVerts;i++) {
    if (bounding.vertexLocations[i*3]>phantom.xbounds[1]) phantom.xbounds[1]=bounding.vertexLocations[i*3];
    if (bounding.vertexLocations[i*3]<phantom.xbounds[0]) phantom.xbounds[0]=bounding.vertexLocations[i*3];
  }
}

DLLEXPORT void set_phantom_info(int numObjs, int *objType, int *clipStInd, int *nPlanes, int *matInd, double *objCent, double *shp, double *Qmat, double *clipNormVec, double *clipDist, double *den, int totalNumPlanes)

{
  phantom.numObjects = numObjs;

  phantom.objectType = my_memcpyi(objType, phantom.objectType, numObjs*sizeof(int));
  phantom.clipStartIndex = my_memcpyi(clipStInd, phantom.clipStartIndex, (numObjs+1)*sizeof(int));
  phantom.nClipPlanes = my_memcpyi(nPlanes, phantom.nClipPlanes, numObjs*sizeof(int));
  phantom.materialInd = my_memcpyi(matInd, phantom.materialInd, numObjs*sizeof(int));

  phantom.objectCenter = my_memcpyd(objCent, phantom.objectCenter, numObjs*sizeof(double)*3);
  phantom.shape = my_memcpyd(shp, phantom.shape, numObjs*sizeof(double));
  phantom.Qmatrix = my_memcpyd(Qmat, phantom.Qmatrix, numObjs*sizeof(double)*18);
  phantom.clipNormalVector = my_memcpyd(clipNormVec, phantom.clipNormalVector, sizeof(double)*totalNumPlanes*3);
  phantom.clipDistance = my_memcpyd(clipDist, phantom.clipDistance, sizeof(double)*totalNumPlanes);
  phantom.density = my_memcpyd(den, phantom.density, numObjs*sizeof(double));
}

DLLEXPORT void set_material_info(int materialCount, int eBinCount, double *muTable)

{
  materials.materialCount = materialCount;
  materials.eBinCount = eBinCount;
  materials.muTable = my_memcpyd(muTable, materials.muTable, eBinCount*materialCount*sizeof(double));
}

DLLEXPORT void intersections_full_list(int *objlist, int lenObjList,double *a,double *alpha,double rayLength,double *t_ends,int *matls,double *dens,int *num_segs)

{
  int i,n,p11,j,pri_list[lenObjList*2],mat_list[lenObjList*2],num_int = 0,out,*ibuff,sp,ep,lm,ibuffs,p,ip,nn,ii,seg = 0;
  double tc2[4],a0[3],st_list[lenObjList*2],en_list[lenObjList*2],den_list[lenObjList*2],lc,ld;       //fixme: malloc instead of lenObjList
  /*  int *ien, *ist;
      double *el, *sl;*/ // replaced these with structures startList and endList
  struct indexed_list *startList;
  struct indexed_list *endList;

  for(n = 0;n<lenObjList;n++){
    i = objlist[n];         // Loop over objects of interest
    p11 = phantom.objectType[i];
    for(j = 0;j<3;j++) a0[j] = (a[j]-phantom.objectCenter[i*3+j]);              // vector pointing from obj. center to src
    for(j = 0;j<4;j++) tc2[j] = 0.0;
    if((p11 != 3) && (p11 != 7))                  // If not a torus or vessel segment
      out = quadratic_intersect(a0,alpha,p11,tc2,i);
    else                                   // The torus ca$e requires solving a quartic
      out = quartic_intersect(a0,alpha,tc2,i);
    if(out)        // At least part of the ray is inside the quadratic/quartic surface
      num_int = clip_all(a,alpha,rayLength,tc2,out,st_list,en_list,den_list,pri_list,mat_list,num_int,i);
  }                            // Loop over objects is done

  if(num_int){
    /*    ien = malloc(sizeof(double)*num_int);
    ist = malloc(sizeof(double)*num_int);
    el = malloc(sizeof(double)*num_int);
    sl = malloc(sizeof(double)*(num_int+1));*/  // replaced ien and el with endList, etc.
    endList = malloc(sizeof(struct indexed_list)*num_int);
    startList = malloc(sizeof(struct indexed_list)*(num_int+1));
    ibuff = malloc(sizeof(int)*num_int);
    for(i = 0;i<num_int;i++) ibuff[i] = 0;
    if (num_int>1){
      //for(i = 0;i<num_int;i++) {ien[i] = i;ist[i] = i;}
      for(i = 0;i<num_int;i++) 
	{
	  endList[i].value = en_list[i];
	  startList[i].value = st_list[i];
	  endList[i].originalIndex = i;
	  startList[i].originalIndex = i;
	}
      //vector = en_list;//for(i = 0;i<num_int;i++) printf("en_vector: %1.12f\r\n",vector[i]);
      qsort(endList,num_int,sizeof(struct indexed_list),compIL);
      //vector = st_list;//for(i = 0;i<num_int;i++) printf("st_vector: %1.12f\r\n",vector[i]);
      qsort(startList,num_int,sizeof(struct indexed_list),compIL);
      //for(i = 0;i<num_int;i++) 
      //{
      //  el[i] = en_list[ien[i]];sl[i] = st_list[ist[i]];
      //  //printf("sl %1.12lf  el %1.12lf  ist: %d  ien: %d\r\n",sl[i],el[i],ist[i],ien[i]);
      //  //if (i>0) if(sl[i]<sl[i-1]) printf("ERROR!! %1.12lf\r\n",a0[10000]);
      //}
      //printf("b%d\r\n",num_int);
    }
    else {
       /*      sl[0] = st_list[0];
      el[0] = en_list[0];
      ien[0] = 0;
      ist[0] = 0;*/
      startList[0].value = st_list[0];
      endList[0].value = en_list[0];
      startList[0].originalIndex = 0;
      endList[0].originalIndex = 0;
     }
    //printf("c%d\r\n",num_int);
    //    sl[num_int] = el[num_int-1]+1;
    startList[num_int].value = endList[num_int-1].value+1;
    sp = 0;ep = 0;lc = 0;ld = 0;lm = MATERIAL_INDEX_BASIS;ibuffs = 0;   // sp points to the current interval start
                                                               // ep points to the current interval End
							       // lc : last change (transition from one object to another)
							       // ld : last density
							       // lm : last material
    while(ep<num_int){
      //printf("E  %1.12lf  %1.12lf  %1.12lf  %d  %d\r\n",sl[0],el[0],sl[1],sp,ep);
      //if (sl[sp]<el[ep]){
      if (startList[sp].value<endList[ep].value){
	//printf("f%d\r\n",num_int);
	p = startList[sp].originalIndex;
	// push into ibuff
	ip = 0;
	while(ibuffs>ip){
	  if(pri_list[p]>pri_list[ibuff[ip]]) 
	    break;
	  ip++;
	}
	for(i = ibuffs;i>ip;i--) ibuff[i] = ibuff[i-1];
	ibuff[ip] = p;
	ibuffs++;
	if(ip == 0) {
	  //	  mb[lm-MATERIAL_INDEX_BASIS] = mb[lm-MATERIAL_INDEX_BASIS]+ld*(st_list[p]-lc);
	  t_ends[seg] = st_list[p];
	  matls[seg] = lm;
	  dens[seg] = ld;
	  //	  printf("dens[%d]: %1.12lf\r\n",seg,dens[seg]);
	  //	  printf("t_ends[%d]: %1.12lf\r\n",seg,t_ends[seg]);
	  seg++;
	  lc = st_list[p];
	  ld = den_list[p];
	  lm = mat_list[p];
	}
	sp++;
      }
      else {
	p = endList[ep].originalIndex;
	// pull out of ibuff
	ip = 0;
	while(ibuff[ip] != p) ip++;
	for(i = ip;i<(ibuffs-1);i++) ibuff[i] = ibuff[i+1];
	ibuffs--;
	if(ip == 0){
	  //	  mb[lm-MATERIAL_INDEX_BASIS] = mb[lm-MATERIAL_INDEX_BASIS]+ld*(en_list[p]-lc);
	  t_ends[seg] = en_list[p];
	  matls[seg] = lm;
	  dens[seg] = ld;
	  //	  printf("dens[%d]: %1.12lf\r\n",seg,dens[seg]);
	  //	  printf("t_ends[%d]: %1.12lf\r\n",seg,t_ends[seg]);
	  seg++;
	  lc = en_list[p];
	  if(ibuffs) {
	    ld = den_list[ibuff[0]];
	    lm = mat_list[ibuff[0]];
	  }
	  else
	    ld = 0;  // lm shouldn't matter when ld = 0 so we leave it alone
	}
	ep++;
      }
    }
    /*    free(ien);
    free(ist);
    free(el);
    free(sl);*/
    free(startList);
    free(endList);
    free(ibuff);
  }
  num_segs[0] = seg;
}

DLLEXPORT void make_vol(float *volume, int Nx, double xoff, double dx, int Ny, double yoff, double dy, int Nz, double zoff, double dz, int oversampling, int UNUSED_num_volumes, int material_volumes)

{
  int i, j, k, l, m, p, *matls, *objlist, num_segs, energyBinsToSkip = 0, volume_offset;
  double a[3], alpha[3], *t_ends, *dens, rayLength, xbrd, *matl_tabl;

  matl_tabl = malloc(sizeof(double)*materials.materialCount);
  matls = malloc(sizeof(int)*phantom.numObjects*4);
  objlist = malloc(sizeof(int)*phantom.numObjects);
  t_ends = malloc(sizeof(double)*phantom.numObjects*4);
  dens = malloc(sizeof(double)*phantom.numObjects*4);
  for(i = 0;i<materials.materialCount;i++) {
    if (material_volumes) 
      matl_tabl[i] = 1.0/(oversampling*oversampling);
    else
      matl_tabl[i] = materials.muTable[i+materials.materialCount*energyBinsToSkip]/(oversampling*oversampling);
    printf("\n\rmat%d:  %1.12lf\n\r",i,matl_tabl[i]*oversampling*oversampling);
  }
  dy = dy/oversampling;
  yoff = oversampling*(yoff-1)+(oversampling+1)*0.5;
  dz = dz/oversampling;
  zoff = oversampling*(zoff-1)+(oversampling+1)*0.5;
  for(i = 0;i<phantom.numObjects;i++) objlist[i] = i;

  alpha[0] = 1;
  alpha[1] = 0;
  alpha[2] = 0;
  a[0] = phantom.xbounds[0] - 1;
  rayLength = phantom.xbounds[1] + 1 - phantom.xbounds[0];
  dbug(1,"\n\n\rphantom.xbounds[0]:  %1.12lf  \n\r phantom.xbounds[0]:  %1.12lf  \n\rmaterial_volumes: %d",phantom.xbounds[0],phantom.xbounds[1],material_volumes);
  for(l = 0;l<Nz;l++)
    for(k = l*oversampling;k<l*oversampling+oversampling;k++){
      a[2] = (k+1-zoff)*dz;
      for(m = 0;m<Ny;m++)
	for(j = m*oversampling;j<m*oversampling+oversampling;j++){
	  a[1] = (j+1-yoff)*dy;
	  num_segs = 0;
	  intersections_full_list(objlist, phantom.numObjects, a, alpha, rayLength, t_ends, matls, dens, &num_segs);
	  for(i = 0;i<num_segs;i++){
	    t_ends[i] = t_ends[i]+a[0];
	    //printf("%1.12lf\r\n",t_ends[i]);
	  }
	  //printf("%d %d %d\r\n",j,k,num_segs);
	  //printf("%d",j);
	  i = 0;p = 0;
	  xbrd = (i+0.5-xoff)*dx;
	  while((t_ends[p]<xbrd) && (p<num_segs)) p++;
	  while((i<Nx) && (p<num_segs)){
            volume_offset = Nx*m+Nx*Ny*l+Nx*Ny*Nz*(matls[p]-MATERIAL_INDEX_BASIS)*material_volumes;
            //            dbug(1,"%d ",volume_offset);
	    while((xbrd<t_ends[p]) && (i<Nx)) {
	      volume[i+volume_offset] += dens[p]*matl_tabl[matls[p]-MATERIAL_INDEX_BASIS];
              //              dbug(1,"%d ",matls[p]-MATERIAL_INDEX_BASIS);
	      //printf("%1.12lf %1.12lf %d\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l],i+Nx*m+Nx*Ny*l);
	      xbrd += dx;i++;
	    }
            // Remove patrial volume contributions
	    if (xbrd>t_ends[p]){
	      i--;xbrd -= dx;
	      volume[i+volume_offset] -= dens[p]*matl_tabl[matls[p]-MATERIAL_INDEX_BASIS]*(dx+xbrd-t_ends[p])/dx;
	      //printf("B  %1.12lf %1.12lf\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l]);
	      if (p+1<num_segs){
                volume_offset = Nx*m+Nx*Ny*l+Nx*Ny*Nz*(matls[p+1]-MATERIAL_INDEX_BASIS)*material_volumes;
		volume[i+volume_offset] -= dens[p+1]*matl_tabl[matls[p+1]-MATERIAL_INDEX_BASIS]*(t_ends[p]-xbrd)/dx;
              }
	      //printf("C  %1.12lf %1.12lf %d\r\n",xbrd,volume[i+Nx*m+Nx*Ny*l],i+Nx*m+Nx*Ny*l);
	    }
	    p++;
	  }
	}
    }
  free(matl_tabl);
  free(matls);
  free(objlist);
  free(t_ends);
  free(dens);
}

DLLEXPORT void intersections(int *objlist,int lenObjList,double *a,double *alpha,double rayLength,double *mb)

{
  int i,n,p11,j,pri_list[lenObjList*2],mat_list[lenObjList*2],num_int = 0,out,*ibuff,sp,ep,lm,ibuffs,p,ip,nn,ii;
  double tc2[4],a0[3],st_list[lenObjList*2],en_list[lenObjList*2],den_list[lenObjList*2],lc,ld;       //fixme: malloc instead of lenObjList
  /*  int *ien, *ist;
      double *el, *sl;*/ // replaced these with structures startList and endList
  struct indexed_list *startList;
  struct indexed_list *endList;

  //  printf("objlist: %d lenObjList: %d a = [%lf %lf %lf] alpha = [%lf %lf %lf] maxD: %lf mb[0]: %lf\n",objlist[0],lenObjList,a[0],a[1],a[2],alpha[0],alpha[1],alpha[2],maxD,mb[0]);
  for(n = 0;n<lenObjList;n++){
    i = objlist[n];         // Loop over objects of interest
    //    if((i == 4) && (prnt)) pnt = 1; else pnt = 0;
    p11 = phantom.objectType[i];
    for(j = 0;j<3;j++) a0[j] = (a[j]-phantom.objectCenter[i*3+j]);              // vector pointing from obj. center to src
    for(j = 0;j<4;j++) tc2[j] = 0.0;
    if((p11 != 3) && (p11 != 7))                  // If not a torus or vessel segment
      out = quadratic_intersect(a0,alpha,p11,tc2,i);
    else                                   // The torus ca$e requires solving a quartic
      out = quartic_intersect(a0,alpha,tc2,i);
    //if((pnt) && (tc2[0]>100.0))
    //  printf("Intersection points for object 5, module 181:\r\n   %1.12lf   %1.12lf\r\n\r\n",tc2[0],tc2[1]);
    //if (out != 0) {for(ii = 0;ii<out;ii++) printf("%1.12lf  ",tc2[ii]); printf("\r\n");}
    if(out)        // At least part of the ray is inside the quadratic/quartic surface
      num_int = clip_all(a,alpha,rayLength,tc2,out,st_list,en_list,den_list,pri_list,mat_list,num_int,i);
  }                            // Loop over objects is done
  //printf("%d\r\n",num_int);

  if(num_int){
    /*    ien = malloc(sizeof(double)*num_int);
    ist = malloc(sizeof(double)*num_int);
    el = malloc(sizeof(double)*num_int);
    sl = malloc(sizeof(double)*(num_int+1));*/  // replaced ien and el with endList, etc.
    endList = malloc(sizeof(struct indexed_list)*num_int);
    startList = malloc(sizeof(struct indexed_list)*(num_int+1));
    ibuff = malloc(sizeof(int)*num_int);
    for(i = 0;i<num_int;i++) ibuff[i] = 0;
  //printf("a%d\r\n",num_int);
    if (num_int>1){
      //for(i = 0;i<num_int;i++) {ien[i] = i;ist[i] = i;}
      for(i = 0;i<num_int;i++) 
	{
	  endList[i].value = en_list[i];
	  startList[i].value = st_list[i];
	  endList[i].originalIndex = i;
	  startList[i].originalIndex = i;
	}
      //vector = en_list;//for(i = 0;i<num_int;i++) printf("en_vector: %1.12f\r\n",vector[i]);
      qsort(endList,num_int,sizeof(struct indexed_list),compIL);
      //vector = st_list;//for(i = 0;i<num_int;i++) printf("st_vector: %1.12f\r\n",vector[i]);
      qsort(startList,num_int,sizeof(struct indexed_list),compIL);
      //for(i = 0;i<num_int;i++) 
      //{
      //  el[i] = en_list[ien[i]];sl[i] = st_list[ist[i]];
      //  //printf("sl %1.12lf  el %1.12lf  ist: %d  ien: %d\r\n",sl[i],el[i],ist[i],ien[i]);
      //  //if (i>0) if(sl[i]<sl[i-1]) printf("ERROR!! %1.12lf\r\n",a0[10000]);
      //}
      //printf("b%d\r\n",num_int);
    }
    else {
      /*      sl[0] = st_list[0];
      el[0] = en_list[0];
      ien[0] = 0;
      ist[0] = 0;*/
      startList[0].value = st_list[0];
      endList[0].value = en_list[0];
      startList[0].originalIndex = 0;
      endList[0].originalIndex = 0;
    }
    //printf("c%d\r\n",num_int);
    //    sl[num_int] = el[num_int-1]+1;
    startList[num_int].value = endList[num_int-1].value+1;
    sp = 0;ep = 0;lc = 0;ld = 0;lm = MATERIAL_INDEX_BASIS;ibuffs = 0;  // sp points to the current interval start
                                                               // ep points to the current interval End
							       // lc : last change (transition from one object to another)
							       // ld : last density
							       // lm : last material
    //printf("d%d\r\n",num_int);
    while(ep<num_int){
      //printf("E  %1.12lf  %1.12lf  %1.12lf  %d  %d\r\n",sl[0],el[0],sl[1],sp,ep);
      //if (sl[sp]<el[ep]){
      if (startList[sp].value<endList[ep].value){
	//printf("f%d\r\n",num_int);
	p = startList[sp].originalIndex;
	// push into ibuff
	ip = 0;
	while(ibuffs>ip){
	  if(pri_list[p]>pri_list[ibuff[ip]]) 
	    break;
	  ip++;
	}
	for(i = ibuffs;i>ip;i--) ibuff[i] = ibuff[i-1];
	ibuff[ip] = p;
	ibuffs++;
	//for(i = 0;i<ibuffs;i++) printf("Ibuff: %d\n\r",ibuff[i]);
	if(ip == 0) {
	  mb[lm-MATERIAL_INDEX_BASIS] = mb[lm-MATERIAL_INDEX_BASIS]+ld*(st_list[p]-lc);
	  lc = st_list[p];
	  ld = den_list[p];
	  lm = mat_list[p];
	}
	sp++;
      }
      else {
	//printf("g%d\r\n",num_int);
	//printf("g%d%d\r\n",num_int,ibuff[0]);
	p = endList[ep].originalIndex;
	// pull out of ibuff
	ip = 0;
  //printf("h%d %d\r\n",p,ibuffs);
  //for(i = 0;i<ibuffs;i++) printf("ibuff/w: %d\n\r",ibuff[i]);
	while(ibuff[ip] != p) ip++;
  //printf("i%d\r\n",num_int);
	for(i = ip;i<(ibuffs-1);i++) ibuff[i] = ibuff[i+1];
  //printf("j%d\r\n",num_int);
	ibuffs--;
	//for(i = 0;i<ibuffs;i++) printf("ibuff/wo: %d\n\r",ibuff[i]);
  //printf("k%d\r\n",num_int);
	if(ip == 0){
  //printf("l   %1.12lf %d %1.12lf %1.12lf %d %1.12lf\r\n",mb[0],lm,ld,en_list[0],p,lc);
	  mb[lm-MATERIAL_INDEX_BASIS] = mb[lm-MATERIAL_INDEX_BASIS]+ld*(en_list[p]-lc);
	  //	  printf("ld: %1.12lf diff: %1.12lf\r\n",ld,en_list[p]-lc);
  //printf("m%d\r\n",num_int);
	  lc = en_list[p];
  //printf("n%d\r\n",num_int);
	  if(ibuffs) {
  //printf("o%d\r\n",num_int);
	    ld = den_list[ibuff[0]];
  //printf("p%d\r\n",num_int);
	    lm = mat_list[ibuff[0]];
  //printf("q%d\r\n",num_int);
	  }
	  else
	    ld = 0;  // lm shouldn't matter when ld = 0 so we leave it alone
	}
	ep++;
      }
    }
    /*    free(ien);
    free(ist);
    free(el);
    free(sl);*/
    free(startList);
    free(endList);
    free(ibuff);
  }
}

int compare_pts (const void *a, const void *b)

{
  int prelim_return;

  const struct pnt *da = (const struct pnt *) a;
  const struct pnt *db = (const struct pnt *) b;
  
  prelim_return = (da[0].tan_angle > db[0].tan_angle) - (da[0].tan_angle < db[0].tan_angle);

  //  dbug(2,"\n\rcompare_pts: %f %f %d",da[0].tan_angle,db[0].tan_angle,prelim_return);

  if (prelim_return == 0)
    {
      prelim_return = (da[0].distance > db[0].distance) - (da[0].distance < db[0].distance);
      //      dbug(2,"\n\r   compare_pts2: %f %f %d",da[0].distance,db[0].distance,prelim_return);
    }
  return prelim_return;
}

int compute_convex_hull_2d(double* x,double* y,int pts,int* list)

/* we use the graham scan algorithm (http://en.wikipedia.org/wiki/Graham_scan) */

{
  int min_index = 0,i,list_length = 2;
  double min_y = y[0];
  struct pnt *points;
  const int one_based = 0;

  for(i=0;i<pts;i++)
    dbug(2,"x[%d] = %12.4lf    y[%d] = %12.4lf\r\n",i,x[i],i,y[i]);
  //dbug(3,"a\n\r");
  //
  //  First we find the point with lowest y coordinate
  //
  for(i = 1;i<pts;i++)
    if (y[i]<min_y) 
      {min_y = y[i];min_index = i;}
  //dbug(3,"f %d %d\n\n\r",sizeof(struct pnt),pts);
  points = malloc(sizeof(struct pnt)*pts);
  //dbug(3,"g\n\r");
  for(i = 0;i<pts;i++)
    {
      //dbug(3,"e   (%d/%d)\n\r",i,pts-1);
      points[i].x = x[i];
      points[i].y = y[i];
      points[i].index = i;
      points[i].tan_angle = -(x[i] - x[min_index]) / (y[i] - y[min_index]);
      points[i].distance = sqrt((x[i] - x[min_index]) * (x[i] - x[min_index]) + (y[i] - y[min_index]) * (y[i] - y[min_index]));
    }
  points[min_index].tan_angle = -1.0 / 0.0;  //-inf (this ensures we get the right starting point)

  //dbug(3,"m\n\r");

  //
  // Next, we sort by angle, breaking ties by comparing the distance (this is the only step that is O(n log n), the others are O(n))
  //
  qsort(points,pts,sizeof(struct pnt),compare_pts);

  for(i = 0;i<pts;i++)
    dbug(2,"points[%d].x = %12.4lf    points[%d].y = %12.4lf\r\n",i,points[i].x,i,points[i].y);

  //
  //  Finally, we loop through the points discarding one whenever we encounter a left turn instead of a right turn
  //

  list[0] = 0;list[1] = 1;
  
  for(i = 2;i<pts;i++)
    {
      while ((list_length >= 2) && (((points[list[list_length-1]].x - points[list[list_length-2]].x) * (points[i].y - points[list[list_length-2]].y) - (points[list[list_length-1]].y - points[list[list_length-2]].y) * (points[i].x - points[list[list_length-2]].x)) <= 0))
	  list_length--;
      list[list_length] = i;
      list_length++;
    }
  //dbug(3,"q\n\r");

  for(i = 0;i<list_length;i++)
    list[i] = points[list[i]].index + one_based;

  for(i = 0;i<list_length;i++)
    dbug(2,"i: %d    list[i]: %d     points[list[i]].x: %5.4lf    points[list[i]].y: %5.4lf\n\r",i,list[i],points[list[i]].x,points[list[i]].y);// this line may have a bug

  free(points);points = NULL;
  //dbug(3,"z\n\r");
  return list_length;
}

void crop_polygon_FMtest(double *x, double *y, int *indexList, int *listLength, double height, double width, int projPoints)

{
  int i,j,currentIndex;
  //dbug(3,"Point A\n\r");
  double halfHeight = height/2;
  double halfWidth = width/2;
  double *distOutside = NULL;
  int numOutside, lastOutside;
  int *nextIndex = NULL;
  int *prevIndex = NULL;
  int startInside, startOutside;
  int choppedPoints = 0;
  double frac;


  //dbug(3,"Point B\n\r");
  distOutside = malloc(sizeof(double)*(listLength[0]+8)); // the +8 is to allow up to two new points for each module edge
  nextIndex = malloc(sizeof(int)*(listLength[0]+8));
  prevIndex = malloc(sizeof(int)*(listLength[0]+8));

  for(i = 0;i<listLength[0]-1;i++)
    {
      nextIndex[i] = i+1;
      prevIndex[i+1] = i;
    }
  nextIndex[listLength[0]-1] = 0;
  prevIndex[0] = listLength[0]-1;
  currentIndex = 0;

  for(j = 0;j<4;j++)
    {
      //dbug(3,"Point %d, currentIndex: %d\n\r",j,currentIndex);
      switch(j)
	{
	case 0:
	  distOutside[currentIndex] = x[indexList[currentIndex]] - halfWidth;
	  break;
	case 1:
	  distOutside[currentIndex] = - x[indexList[currentIndex]] - halfWidth;
	  break;
	case 2:
	  distOutside[currentIndex] = y[indexList[currentIndex]] - halfHeight;
	  break;
	case 3:
	  distOutside[currentIndex] = - y[indexList[currentIndex]] - halfHeight;
	  break;
	}
      //dbug(3,"Point %d\n\r",j);
      lastOutside = (distOutside[currentIndex] > 0);
      numOutside = lastOutside;
      startInside = -1;
      startOutside = -1;
      currentIndex = nextIndex[currentIndex];
      for(i = 1; i<listLength[0]-choppedPoints; i++, currentIndex = nextIndex[currentIndex])
	{
	  //dbug(3,"Point j: %d  , i: %d\n\r",j,i);
	  switch(j)
	    {
	    case 0:
	      distOutside[currentIndex] = x[indexList[currentIndex]] - halfWidth;
	      break;
	    case 1:
	      distOutside[currentIndex] = - x[indexList[currentIndex]] - halfWidth;
	      break;
	    case 2:
	      distOutside[currentIndex] = y[indexList[currentIndex]] - halfHeight;
	      break;
	    case 3:
	      distOutside[currentIndex] = - y[indexList[currentIndex]] - halfHeight;
	      break;
	    }
	  numOutside += (distOutside[currentIndex] > 0); 
	  if ((lastOutside == 0) && (distOutside[currentIndex] > 0))
	    startOutside = currentIndex;
	  if ((lastOutside == 1) && (distOutside[currentIndex] < 0))
	    startInside = currentIndex;
	  //dbug(2,"lastOutside: %d  currentIndex: %d  nextindex[currentindex]: %d\n\r",lastOutside,currentIndex, nextIndex[currentIndex]);
	  lastOutside = (distOutside[currentIndex] > 0);
	}
      if ((lastOutside == 0) && (distOutside[currentIndex] > 0))
	startOutside = currentIndex;
      if ((lastOutside == 1) && (distOutside[currentIndex] < 0))
	startInside = currentIndex;
      //dbug(3,"startOutside: %d startInside: %d, currentIndex: %d\n\r",startOutside,startInside,currentIndex);
      if (numOutside == listLength[0] - choppedPoints)
	{
	  free(distOutside);distOutside = NULL;
	  free(nextIndex);nextIndex = NULL;
	  free(prevIndex);prevIndex = NULL;
	  //dbug(3,"**********************************free %d\n\r",j);
	  listLength[0] = 0;
	  return;
	}
      //dbug(3,"Point %d\n\r",j);
      choppedPoints += numOutside;
      if (startOutside != -1)
	{
	  //Find new polygon vertex
	  frac = distOutside[startOutside] / (distOutside[startOutside] - distOutside[prevIndex[startOutside]]);
	  switch(j)
	    {
	    case 0:
	      y[projPoints] = y[indexList[startOutside]] + frac * (y[indexList[prevIndex[startOutside]]] - y[indexList[startOutside]]);
	      x[projPoints] = halfWidth;
	      break;
	    case 1:
	      y[projPoints] = y[indexList[startOutside]] + frac * (y[indexList[prevIndex[startOutside]]] - y[indexList[startOutside]]);
	      x[projPoints] = - halfWidth;
	      break;
	    case 2:
	      x[projPoints] = x[indexList[startOutside]] + frac * (x[indexList[prevIndex[startOutside]]] - x[indexList[startOutside]]);
	      y[projPoints] = halfHeight;
	      break;
	    case 3:
	      x[projPoints] = x[indexList[startOutside]] + frac * (x[indexList[prevIndex[startOutside]]] - x[indexList[startOutside]]);
	      y[projPoints] = - halfHeight;
	      break;
	    }
	  //Expand index list and update nextIndex and prevIndex
	  indexList[listLength[0]] = projPoints;
	  nextIndex[prevIndex[startOutside]] = listLength[0];
	  prevIndex[listLength[0]] = prevIndex[startOutside];
	  nextIndex[listLength[0]] = listLength[0] + 1;
	  //Increment counters
	  listLength[0]++;
	  projPoints++;

	  //Find other new polygon vertex
	  frac = distOutside[startInside] / (distOutside[startInside] - distOutside[prevIndex[startInside]]);
	  switch(j)
	    {
	    case 0:
	      y[projPoints] = y[indexList[startInside]] + frac * (y[indexList[prevIndex[startInside]]] - y[indexList[startInside]]);
	      x[projPoints] = halfWidth;
	      break;
	    case 1:
	      y[projPoints] = y[indexList[startInside]] + frac * (y[indexList[prevIndex[startInside]]] - y[indexList[startInside]]);
	      x[projPoints] = - halfWidth;
	      break;
	    case 2:
	      x[projPoints] = x[indexList[startInside]] + frac * (x[indexList[prevIndex[startInside]]] - x[indexList[startInside]]);
	      y[projPoints] = halfHeight;
	      break;
	    case 3:
	      x[projPoints] = x[indexList[startInside]] + frac * (x[indexList[prevIndex[startInside]]] - x[indexList[startInside]]);
	      y[projPoints] = - halfHeight;
	      break;
	    }
	  //Expand index list and update nextIndex and prevIndex
	  indexList[listLength[0]] = projPoints;
	  prevIndex[startInside] = listLength[0];
	  prevIndex[listLength[0]] = listLength[0] - 1;
	  nextIndex[listLength[0]] = startInside;
	  //Increment counters
	  listLength[0]++;
	  projPoints++;
	  currentIndex = startInside;
	}
      //for(i=0;i<listLength[0];i++)
      //  dbug(3,"x:  %1.4f     y:  %1.4f\n\r",x[indexList[i]],y[indexList[i]]);
    }

  listLength[0] = listLength[0] - choppedPoints;
  indexList[0] = currentIndex;
  for(i=1;i<listLength[0];i++)
    {
      currentIndex = nextIndex[currentIndex];
      indexList[i] = currentIndex;
    }
  free(distOutside);distOutside = NULL;
  free(nextIndex);nextIndex = NULL;
  free(prevIndex);prevIndex = NULL;
}

void crop_polygon(double *x, double *y, int *indexList, int *listLength, double height, double width, int projPoints)

{
  int i,j,currentIndex;
  double halfHeight = height/2;
  double halfWidth = width/2;
  double *distOutside = NULL;
  int numOutside, lastOutside;
  int *nextIndex = NULL;
  int *prevIndex = NULL;
  int *tempList = NULL;
  int startInside, startOutside;
  int choppedPoints = 0;
  double frac;
  double maxabsy = 0,maxabsx = 0;
  //  int ci;

  for(i=0;i<listLength[0];i++)
    {
      maxabsy=MAX(maxabsy,fabs(y[indexList[i]]));
      maxabsx=MAX(maxabsx,fabs(x[indexList[i]]));
      //dbug(2,"i:%4d   x[%d] = %12.4lf    y[%d] = %12.4lf\r\n",i,indexList[i],x[indexList[i]],indexList[i],y[indexList[i]]);
    }
  //dbug(3,"mx: %8.4lf my: %8.4lf",maxabsx,maxabsy);
  //dbug(2,"ll:%d\n\r",listLength[0]);
  distOutside = malloc(sizeof(double)*(listLength[0]+8)); // the +8 is to allow up to two new points for each module edge
  nextIndex = malloc(sizeof(int)*(listLength[0]+8));
  prevIndex = malloc(sizeof(int)*(listLength[0]+8));

  for(i = 0;i<listLength[0]-1;i++)
    {
      nextIndex[i] = i+1;
      prevIndex[i+1] = i;
    }
  nextIndex[listLength[0]-1] = 0;
  prevIndex[0] = listLength[0]-1;
  currentIndex = 0;

  for(j = 0;j<4;j++)
    {
      //dbug(2,"-------------------- j :%d\n\r",j);
      switch(j)
	{
	case 0:
	  distOutside[currentIndex] = x[indexList[currentIndex]] - halfWidth;
	  break;
	case 1:
	  distOutside[currentIndex] = - x[indexList[currentIndex]] - halfWidth;
	  break;
	case 2:
	  distOutside[currentIndex] = y[indexList[currentIndex]] - halfHeight;
	  break;
	case 3:
	  distOutside[currentIndex] = - y[indexList[currentIndex]] - halfHeight;
	  break;
	}
      //dbug(2,"currentIndex:%d\n\r",currentIndex);
      //if (x[currentIndex]>0) dbug(3,":");
      //if (y[currentIndex]>0) dbug(3,":");
      //if (distOutside[currentIndex]>0) dbug(3,":");
      //if (currentIndex>0) dbug(3,"*");
      lastOutside = (distOutside[currentIndex] > 0);
      numOutside = lastOutside;
      startInside = -1;
      startOutside = -1;
      currentIndex = nextIndex[currentIndex];
      //dbug(2,"currentIndex:%d\n\r",currentIndex);
      for(i = 1; i<listLength[0]-choppedPoints; i++, currentIndex = nextIndex[currentIndex])
	{
	  //if (choppedPoints>0) dbug(3,"|");
	  switch(j)
	    {
	    case 0:
	      distOutside[currentIndex] = x[indexList[currentIndex]] - halfWidth;
	      break;
	    case 1:
	      distOutside[currentIndex] = - x[indexList[currentIndex]] - halfWidth;
	      break;
	    case 2:
	      distOutside[currentIndex] = y[indexList[currentIndex]] - halfHeight;
	      break;
	    case 3:
	      distOutside[currentIndex] = - y[indexList[currentIndex]] - halfHeight;
	      break;
	    }
	  numOutside += (distOutside[currentIndex] > 0); 
	  if ((lastOutside == 0) && (distOutside[currentIndex] > 0))
	      startOutside = currentIndex;
	  if ((lastOutside == 1) && (distOutside[currentIndex] <= 0))
	      startInside = currentIndex;
	  lastOutside = (distOutside[currentIndex] > 0);
          //dbug(2,"currentIndex:%d  i:%d  distOutside: %1.12lf  outside: %d prevIndex: %d nextIndex: %d\n\r",currentIndex,i,distOutside[currentIndex],distOutside[currentIndex] > 0, prevIndex[currentIndex], nextIndex[currentIndex]);
	}
      if ((lastOutside == 0) && (distOutside[currentIndex] > 0))
	startOutside = currentIndex;
      if ((lastOutside == 1) && (distOutside[currentIndex] <= 0))
	startInside = currentIndex;
      //dbug(2,"currentIndex:%d  i:%d  distOutside: %1.12lf  outside: %d prevIndex: %d nextIndex: %d\n\r",currentIndex,i,distOutside[currentIndex],distOutside[currentIndex] > 0, prevIndex[currentIndex], nextIndex[currentIndex]);
      //dbug(2,"numOutside: %d\n\r",numOutside);
      if (numOutside == listLength[0] - choppedPoints)
	{
	  free(distOutside);distOutside = NULL;
	  free(nextIndex);nextIndex = NULL;
	  free(prevIndex);prevIndex = NULL;
	  listLength[0] = 0;
	  return;
	}
      //      if (numOutside>0) dbug(2,".");
      choppedPoints += numOutside;
      if (startOutside != -1)
	{
	  //Find new polygon vertex
	  frac = distOutside[startOutside] / (distOutside[startOutside] - distOutside[prevIndex[startOutside]]);
	  switch(j)
	    {
	    case 0:
	      y[projPoints] = y[indexList[startOutside]] + frac * (y[indexList[prevIndex[startOutside]]] - y[indexList[startOutside]]);
	      x[projPoints] = halfWidth;
	      break;
	    case 1:
	      y[projPoints] = y[indexList[startOutside]] + frac * (y[indexList[prevIndex[startOutside]]] - y[indexList[startOutside]]);
	      x[projPoints] = - halfWidth;
	      break;
	    case 2:
	      x[projPoints] = x[indexList[startOutside]] + frac * (x[indexList[prevIndex[startOutside]]] - x[indexList[startOutside]]);
	      y[projPoints] = halfHeight;
	      break;
	    case 3:
	      x[projPoints] = x[indexList[startOutside]] + frac * (x[indexList[prevIndex[startOutside]]] - x[indexList[startOutside]]);
	      y[projPoints] = - halfHeight;
	      break;
	    }
	  //Expand index list and update nextIndex and prevIndex
	  indexList[listLength[0]] = projPoints;
	  nextIndex[prevIndex[startOutside]] = listLength[0];
	  prevIndex[listLength[0]] = prevIndex[startOutside];
	  nextIndex[listLength[0]] = listLength[0] + 1;
	  //Increment counters
	  listLength[0]++;
	  //dbug(2,"ll1:%d\n\r",listLength[0]);
	  projPoints++;

	  //Find other new polygon vertex
	  frac = distOutside[startInside] / (distOutside[startInside] - distOutside[prevIndex[startInside]]);
	  switch(j)
	    {
	    case 0:
              //dbug(2,"ll0:%d\n\r",listLength[0]);
	      y[projPoints] = y[indexList[startInside]] + frac * (y[indexList[prevIndex[startInside]]] - y[indexList[startInside]]);
	      x[projPoints] = halfWidth;
              //dbug(2,"ll0:%d\n\r",listLength[0]);
	      break;
	    case 1:
              //dbug(2,"ll1:%d\n\r",listLength[0]);
	      y[projPoints] = y[indexList[startInside]] + frac * (y[indexList[prevIndex[startInside]]] - y[indexList[startInside]]);
	      x[projPoints] = - halfWidth;
              //dbug(2,"ll1:%d\n\r",listLength[0]);
	      break;
	    case 2:
              //dbug(2,"pp %d     si %d\r\n",projPoints,startInside);
	      x[projPoints] = x[indexList[startInside]] + frac * (x[indexList[prevIndex[startInside]]] - x[indexList[startInside]]);
	      y[projPoints] = halfHeight;
              //dbug(2,"ll2:%d\n\r",listLength[0]);
	      break;
	    case 3:
              //dbug(2,"ll3:%d\n\r",listLength[0]);
	      x[projPoints] = x[indexList[startInside]] + frac * (x[indexList[prevIndex[startInside]]] - x[indexList[startInside]]);
	      y[projPoints] = - halfHeight;
              //dbug(2,"ll3:%d\n\r",listLength[0]);
	      break;
	    }
	  //dbug(2,"ll1:%d\n\r",listLength[0]);
	  //Expand index list and update nextIndex and prevIndex
	  indexList[listLength[0]] = projPoints;
	  prevIndex[startInside] = listLength[0];
	  prevIndex[listLength[0]] = listLength[0] - 1;
	  nextIndex[listLength[0]] = startInside;
	  //Increment counters
	  listLength[0]++;
	  dbug(2,"ll2:%d\n\r",listLength[0]);
	  projPoints++;
	  currentIndex = startInside;
	}
      //ci = currentIndex;
      //for(i=0;i<listLength[0]-choppedPoints;i++)
      //  {
      //    dbug(2,"ci:%4d   x[%d] = %12.4lf    y[%d] = %12.4lf\r\n",ci,indexList[ci],x[indexList[ci]],indexList[ci],y[indexList[ci]]);
      //    ci=nextIndex[ci];
      //  }
    }

  listLength[0] = listLength[0] - choppedPoints;
  //dbug(2,"ll3:%d\n\r",listLength[0]);
  tempList = malloc(listLength[0]*sizeof(int));
  tempList[0] = indexList[currentIndex];
  for(i=1;i<listLength[0];i++)
    {
      currentIndex = nextIndex[currentIndex];
      tempList[i] = indexList[currentIndex];
    }
  for(i=0;i<listLength[0];i++)
    {
      indexList[i] = tempList[i];
      //dbug(2,"i: %8d, indexList:%4d  x[]: %4.2lf y[]:%4.2lf \r\n",i,indexList[i],x[indexList[i]],y[indexList[i]]);
    }
  //dbug(2,"prefree\n\r");
  free(distOutside);distOutside = NULL;
  free(nextIndex);nextIndex = NULL;
  free(prevIndex);prevIndex = NULL;
  free(tempList);tempList = NULL;
  //dbug(2,"postfree\n\r");
}

void store_height_lims(double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries)

{
  int i;
  double minHt = pixel_U[indexList[0]];
  double maxHt = pixel_U[indexList[0]];
  struct height_lims *heightLims = (struct height_lims *) boundaries;

  if (listLength == 0)
    {
      heightLims[objectNumber].min = -VERY_BIG;
      heightLims[objectNumber].max = -VERY_BIG;
      return;
    }
  for(i = 1;i<listLength;i++)
    {
      if (maxHt < pixel_U[indexList[i]])
	maxHt = pixel_U[indexList[i]];
      if (minHt > pixel_U[indexList[i]])
	minHt = pixel_U[indexList[i]];
    }
  heightLims[objectNumber].min = minHt;
  heightLims[objectNumber].max = maxHt;
  //dbug(1,"store   --   minHt: %8.4lf     maxHt:  %8.4lf \n\r",minHt,maxHt);
}

void store_box_lims(double *pixel_R, double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries)

{
  int i;
  double minHt = pixel_U[indexList[0]];
  double maxHt = pixel_U[indexList[0]];
  double minWd = pixel_R[indexList[0]];
  double maxWd = pixel_R[indexList[0]];
  struct box_lims *boxLims = (struct box_lims *) boundaries;

  dbug(3,"store_box_lims -- listLength: %d",listLength);
  if (listLength == 0)
    {
      boxLims[objectNumber].minR = -VERY_BIG;
      boxLims[objectNumber].maxR = -VERY_BIG;
      boxLims[objectNumber].minU = -VERY_BIG;
      boxLims[objectNumber].maxU = -VERY_BIG;
      return;
    }
  for(i = 1;i<listLength;i++)
    {
      if (maxHt < pixel_U[indexList[i]])
	maxHt = pixel_U[indexList[i]];
      if (minHt > pixel_U[indexList[i]])
	minHt = pixel_U[indexList[i]];
      if (maxWd < pixel_R[indexList[i]])
	maxWd = pixel_R[indexList[i]];
      if (minWd > pixel_R[indexList[i]])
	minWd = pixel_R[indexList[i]];
    }
  boxLims[objectNumber].minU = minHt;
  boxLims[objectNumber].maxU = maxHt;
  boxLims[objectNumber].minR = minWd;
  boxLims[objectNumber].maxR = maxWd;
}

void store_half_planes(double *pixel_R, double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries)

{
  printf("ERROR: This type of overlap specification is not yet supported!\n\r");
}

void store(double *pixel_R, double *pixel_U, int *indexList, int listLength, int objectNumber, void *boundaries)

{
  switch(modules.moduleOverlapType)
    {
    case 1:
      store_height_lims(pixel_U, indexList, listLength, objectNumber, boundaries);
      break;
    case 2:
      store_box_lims(pixel_R, pixel_U, indexList, listLength, objectNumber, boundaries);
      break;
    case 3:
      store_half_planes(pixel_R, pixel_U, indexList, listLength, objectNumber, boundaries);
      break;
    default:
      printf("\nERROR: Unrecognized overlap specification!\n\r");
      break;
    }
}


void compute_object_projections(double *srcHullPoints, int nSrcHullPoints, int moduleTypeIndex, double *up, double *right, double *center, void *boundaries)

{
  int i, j, k, n;
  int startIndex;
  double *vertex_XYZ, vertex_RUN[3], source_RUN[3];  // RUN : right, up, normal (module coordinates)
  double normal[3];
  double posRelToModule[3];
  double *pixel_R = NULL;
  double *pixel_U = NULL;
  double scale, Ndist;
  int *indexList = NULL;
  int maxVertexPoints = 0, projPoints, listLength;
  const int moduleEdges = 4;
  int *vertexPoints = NULL;
  int planB;

  vertexPoints = malloc(sizeof(int)*phantom.numObjects);
  //dbug(3,"COP1\r\n");
  cross(right,up,normal);

  // Compute the maximum number of vertex points on a bounding polyhedron
  for(i = 0;i<phantom.numObjects;i++)
    {
      vertexPoints[i] = bounding.vertexStartIndex[i+1] - bounding.vertexStartIndex[i];
      if (vertexPoints[i]>maxVertexPoints) maxVertexPoints = vertexPoints[i];
    }
  //dbug(2,"COP1\r\n");
  pixel_R = malloc(sizeof(double)*(maxVertexPoints*nSrcHullPoints+moduleEdges*2));  // Each module edge can introduce 2 new points during the cropping of the convex hull.
  pixel_U = malloc(sizeof(double)*(maxVertexPoints*nSrcHullPoints+moduleEdges*2));
  indexList = malloc(sizeof(int)*(maxVertexPoints*nSrcHullPoints+moduleEdges*2));

  //dbug(3,"COP1\r\n");
  for(i = 0;i<phantom.numObjects;i++)
    {
      //dbug(3,"Loop %d\r\n",i);
      planB = 0;
      startIndex = bounding.vertexStartIndex[i];
      vertex_XYZ = &bounding.vertexLocations[startIndex*3];
      //dbug(2,"\n\n\rvertex_XYZ:\n\r");
      //for(j=0;j<vertexPoints[i];j++)
      //  dbug(2,"%d:  x[] = %12.4lf    y[] = %12.4lf   z[] = %12.4lf\r\n",j,vertex_XYZ[j*3],vertex_XYZ[j*3+1],vertex_XYZ[j*3+2]);
      //dbug(2,"\n\n\rsource_hull_XYZ:\n\r");
      //for(j=0;j<nSrcHullPoints;j++)
      //  dbug(2,"%d:  x[] = %12.4lf    y[] = %12.4lf   z[] = %12.4lf\r\n",j,srcHullPoints[j*3],srcHullPoints[j*3+1],srcHullPoints[j*3+2]);
      projPoints = 0;
      //dbug(3,"loopbeg %d %d\r\n",j,nSrcHullPoints);
      for(j = 0;j<nSrcHullPoints;j++)
	{
	  //dbug(3,"loop %d %d\r\n",j,nSrcHullPoints);
	  for(n=0;n<3;n++)
	    posRelToModule[n] = srcHullPoints[3*j+n] - center[n];
	  source_RUN[0] = right[0]  * posRelToModule[0] + right[1]  * posRelToModule[1] + right[2]  * posRelToModule[2];
	  source_RUN[1] = up[0]     * posRelToModule[0] + up[1]     * posRelToModule[1] + up[2]     * posRelToModule[2];
	  source_RUN[2] = normal[0] * posRelToModule[0] + normal[1] * posRelToModule[1] + normal[2] * posRelToModule[2];
	  // Project object onto module
	  for(k = 0;k<vertexPoints[i];k++)
	    {
	      for(n=0;n<3;n++)
		posRelToModule[n] = vertex_XYZ[3*k+n] - center[n];
	      dbug(2,"posRel_X: %5.4lf  posRel_Y: %5.4lf  posRel_Z: %5.4lf\n\r",posRelToModule[0],posRelToModule[1],posRelToModule[2]);
	      vertex_RUN[0] = right[0]  * posRelToModule[0] + right[1]  * posRelToModule[1] + right[2]  * posRelToModule[2];
	      vertex_RUN[1] = up[0]     * posRelToModule[0] + up[1]     * posRelToModule[1] + up[2]     * posRelToModule[2];
	      vertex_RUN[2] = normal[0] * posRelToModule[0] + normal[1] * posRelToModule[1] + normal[2] * posRelToModule[2];
	      dbug(2,"posRel_R: %5.4lf  posRel_U: %5.4lf  posRel_N: %5.4lf\n\r",vertex_RUN[0],vertex_RUN[1],vertex_RUN[2]);
	      //dbug(2,"posRel_R: %5.4lf  \n\r",vertex_RUN[0]);
              scale = source_RUN[2] / (source_RUN[2] - vertex_RUN[2]);
              if (fabs(source_RUN[2]) - fabs(vertex_RUN[2]) < 1e-7)
                {
                  dbug(2,"center: %5.4lf %5.4lf %5.4lf\n\r",center[0],center[1],center[2]);
                  dbug(2,"src_pt: %5.4lf %5.4lf %5.4lf\n\r",srcHullPoints[3*j],srcHullPoints[3*j+1],srcHullPoints[3*j+2]);
                  dbug(2,"sourceRUN_R: %5.4lf  sourceRUN_U: %5.4lf  sourceRUN_N: %5.4lf\n\r",source_RUN[0],source_RUN[1],source_RUN[2]);
                  dbug(2,"vertexRUN_R: %5.4lf  vertexRUN_U: %5.4lf  vertexRUN_N: %5.4lf\n\r",vertex_RUN[0],vertex_RUN[1],vertex_RUN[2]);
                  dbug(2,"vertexXYZ_X: %5.4lf  vertexXYZ_Y: %5.4lf  vertexXYZ_Z: %5.4lf\n\r",vertex_XYZ[0],vertex_XYZ[1],vertex_XYZ[2]);
		  if (planB_warned == 0)
		    {
		      printf("\rWarning : Bounding polygon vertex is farther from a detector module than the source along the module normal vector direction. This may indicate that the phantom is too big. Please ensure that phantom_scale is not too large. You may also want to check the detector and focal spot geometry.  In any case, this can slow the simulation down if it occurs frequently.\n\r");
		      planB_warned = 1;
		    }
		  //		  printf("|%d,%d,%d,%d|",i,j,k,(int)center[2]);
		  planB = 1;
		  projPoints=4;
		  pixel_R[0] = -VERY_BIG;pixel_U[0] = -VERY_BIG;
		  pixel_R[1] = -VERY_BIG;pixel_U[1] = VERY_BIG;
		  pixel_R[2] = VERY_BIG;pixel_U[2] = VERY_BIG;
		  pixel_R[3] = VERY_BIG;pixel_U[3] = -VERY_BIG;
		  break;
                }
	      pixel_R[projPoints] = source_RUN[0] + (vertex_RUN[0] - source_RUN[0]) * scale;
	      pixel_U[projPoints] = source_RUN[1] + (vertex_RUN[1] - source_RUN[1]) * scale;
	      projPoints++;
	    }
	  if (planB) break;
	}
      //dbug(3,"loopend %d %d\r\n",j,nSrcHullPoints);
      //dbug(3,"loopend %d %d\r\n",j,projPoints);
      listLength = compute_convex_hull_2d(pixel_R, pixel_U, projPoints, indexList);
      //dbug(2,"loopend %d %d\r\n",j,nSrcHullPoints);
      //dbug(3,"loopend2 %d %d\r\n",j,nSrcHullPoints);
      //dbug(3,"listLength: %d \r\n",listLength);
      //if ((moduleTypeIndex == 146)||(moduleTypeIndex == 147)) debug_flag=2; else debug_flag = 1;
      crop_polygon(pixel_R, pixel_U, indexList, &listLength, modules.Height[moduleTypeIndex], modules.Width[moduleTypeIndex], projPoints);
      //dbug(2,"loopend %d %d\r\n",j,nSrcHullPoints);
      //debug_flag=1;
      //if (listLength>0) dbug(2,"listLength: %d\r\n",listLength);
      store(pixel_R, pixel_U, indexList, listLength, i, boundaries);
      //dbug(2,"loopend %d %d\r\n",j,nSrcHullPoints);
      //dbug(3,"loopend4 %d %d\r\n",j,nSrcHullPoints);
    }
  //dbug(2,"free0\r\n");
  free(pixel_R);pixel_R = NULL;
  //dbug(2,"free0\r\n");
  free(pixel_U);pixel_U = NULL;
  //dbug(2,"free0\r\n");
  free(indexList);indexList = NULL;
  //dbug(2,"free0\r\n");
  free(vertexPoints);vertexPoints = NULL;
  //dbug(2,"free0\r\n");
}

int any_objects_1(int moduleNumber, void *boundaries)

{
  int any_objs = 0, i;
  struct height_lims *heightLims = (struct height_lims *) boundaries;

  for(i = 0;i < phantom.numObjects;i++)
    {
      //dbug(3,"Module number: %d minHeight[%d] = %lf   maxHeight[%d] = %lf     Different?  %d",moduleNumber,i,heightLims.min[i],i,heightLims.max[i],(int)(heightLims.min[i] != heightLims.max[i]));
      if (heightLims[i].min != heightLims[i].max)
	{any_objs++;break;}
    }
  return any_objs;
}

int any_objects_2(int moduleNumber, void *boundaries)

{
  int any_objs = 0, i;
  struct box_lims *boxLims = (struct box_lims *) boundaries;

  for(i = 0;i < phantom.numObjects;i++)
    {
      if (boxLims[i].minR != boxLims[i].maxR)
	{any_objs++;break;}
    }
  return any_objs;
}

void build_object_list1(double *pix_vlims, int *objectList, int *n_objlist, int moduleNumber, double v, void *boundaries)

{
  int i;
  struct height_lims *heightLims = (struct height_lims *) boundaries;

  for(i = 0;i<phantom.numObjects;i++)
    if((heightLims[i].min <= v+pix_vlims[1]) && (heightLims[i].max >= v+pix_vlims[0]))
      {objectList[n_objlist[0]] = i;n_objlist[0]++;}
}

void build_object_list2(double *pix_ulims, double *pix_vlims, int *objectList, int *n_objlist, int moduleNumber, double *uv, void *boundaries)

{
  int i;
  struct box_lims *boxLims = (struct box_lims *) boundaries;

  for(i = 0;i<phantom.numObjects;i++)
    {
      dbug(3,"if statement parts: %d %d %d %d\n\r",(int) (boxLims[i].minU <= uv[1]+pix_vlims[1]),(int) (boxLims[i].maxU >= uv[1]+pix_vlims[0]),(int) (boxLims[i].minR <= uv[0]+pix_ulims[1]),(int) (boxLims[i].maxR >= uv[0]+pix_ulims[0]));
      dbug(3,"parts for third one: %1.5lf %1.5lf %1.5lf %1.5lf \r\n ",boxLims[0].minR,uv[0],pix_ulims[1],uv[0]+pix_ulims[1]);
      if((boxLims[i].minU <= uv[1]+pix_vlims[1]) && (boxLims[i].maxU >= uv[1]+pix_vlims[0]) && (boxLims[i].minR <= uv[0]+pix_ulims[1]) && (boxLims[i].maxR >= uv[0]+pix_ulims[0]))
        {objectList[n_objlist[0]] = i;n_objlist[0]++;}
    }
}

void build_object_list3()

{
  printf("ERROR: this moduleOverlapType not yet implemented!\n\n");
}

void intersections_loop(double *detCenter, double *right, double *up, double *sampling, int nSubDets, double *sourcePoints, int *objectList, int nListObjects, double *thisView, int detIndex, double subviewWeight, double *detWeights)

{
  int i, k, l, m, n;  //loop counters
  double subDetCenter[3];
  double alpha[3];
  double rayLength;
  double *materialBuffer = NULL;
  double *pValueSpectrum = NULL;

  materialBuffer = malloc(sizeof(double)*materials.materialCount);
  pValueSpectrum = malloc(sizeof(double)*materials.eBinCount);
  for(l = 0;l<nSubDets;l++)
    {
      for(i = 0;i<3;i++)
	subDetCenter[i] = detCenter[i]+right[i]*sampling[2*l]+up[i]*sampling[2*l+1];
      for(k = 0;k<modules.nSubSources;k++)
	{
	  for(i = 0;i<3;i++)
	    alpha[i] = subDetCenter[i]-sourcePoints[3*k+i];
	  rayLength = sqrt(alpha[0]*alpha[0]+alpha[1]*alpha[1]+alpha[2]*alpha[2]);
	  alpha[0] /= rayLength;
	  alpha[1] /= rayLength;
	  alpha[2] /= rayLength;
	  for(i = 0;i<materials.materialCount;i++) materialBuffer[i] = 0.0;
	  //dbug(3,"sourcePoints:   [%4.4lf %4.4lf %4.4lf]\r\n",sourcePoints[3*k],sourcePoints[3*k+1],sourcePoints[3*k+2]);
	  //dbug(3,"alpha:   [%4.4lf %4.4lf %4.4lf]\r\n",alpha[0],alpha[1],alpha[2]);
	  //dbug(3,"rayLength:   [%4.4lf]\r\n",rayLength);
	  //dbug(3,"phantom.objectCenter:   [%4.4lf %4.4lf %4.4lf]\r\n",phantom.objectCenter[0],phantom.objectCenter[1],phantom.objectCenter[2]);
	  intersections(&objectList[0],nListObjects,&sourcePoints[3*k],alpha,rayLength,&materialBuffer[0]);
	  //if (materialBuffer[0]>0.0) dbug(1,"materialBuffer[0]: %4.4lf\r\n",materialBuffer[0]);
	  for(m = 0;m<materials.eBinCount;m++){
	    pValueSpectrum[m] = 0;
	    for(n = 0;n<materials.materialCount;n++)
	      pValueSpectrum[m] += materialBuffer[n]*materials.muTable[n+materials.materialCount*m];
	  }
	  //if (pValueSpectrum[0]>0.0) dbug(1,"pValueSpectrum[0]: %4.4lf\r\n",pValueSpectrum[0]);

	  //------------------- Accurate Detector Model, Mingye
	  if(Accurate_Detector_Model_is_ON)
	  {
		  for(m = 0;m < materials.eBinCount;m++)
		  {
			  int the_index = (detIndex*n_col_oversample + (int)(l/n_row_oversample_add_xtalk))*materials.eBinCount + m;
			  thisView[the_index] += subviewWeight*detWeights[l]*modules.sourceWeights[k]*exp(-pValueSpectrum[m]);
		  }
	  }
	  //-------------------
	  else
	  {
		  for(m = 0;m < materials.eBinCount;m++)
			thisView[detIndex*materials.eBinCount+m] += subviewWeight*detWeights[l]*modules.sourceWeights[k]*exp(-pValueSpectrum[m]);
	  }
	}
    }
  free(materialBuffer);materialBuffer = NULL;
  free(pValueSpectrum);pValueSpectrum = NULL;
}

//initialize parameters for Accurate Detector Model, Mingye
void set_Accurate_Detector_Model(double *Paras)
{
	if(Paras[1]==1 && Paras[2]==1) //first view and first subview
	{
		Accurate_Detector_Model_is_ON=1;
		n_col_oversample=Paras[3];
		n_row_oversample=Paras[4];
		n_col_oversample_add_xtalk=Paras[6];
		n_row_oversample_add_xtalk=Paras[8];
	}

	if(Paras[5]==1 && Paras[1]<3 && Paras[2]==1) //col_xtalk>0, first subview in the first two views 
		if(modules.Sub[0]>n_col_oversample*n_row_oversample_add_xtalk)
		{
			int n1;
			
			modules.maxSubDets-=2*n_row_oversample_add_xtalk;
			modules.Sub[0]-=2*n_row_oversample_add_xtalk;

			double sum_weights=0;
			for(n1=0;n1<modules.maxSubDets;n1++)
			{
				modules.Weight[n1]=modules.Weight[n1+n_row_oversample_add_xtalk];
				sum_weights+=modules.Weight[n1];
			}
			for(n1=0;n1<modules.maxSubDets;n1++)
				modules.Weight[n1]/=sum_weights;

			for(n1=0;n1<2*modules.maxSubDets;n1++)
				modules.Sampling[n1]=modules.Sampling[n1+2*n_row_oversample_add_xtalk];
		}

}

DLLEXPORT void Projector(double *Paras, double subviewWeight, double *thisView, double *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength)

{
	//------------------- Accurate Detector Model, Mingye
	if(Paras[0]==1)
		set_Accurate_Detector_Model(Paras);
	else if(Paras[0]==0)
		Accurate_Detector_Model_is_ON=0;

	//FILE *fp;
	//fp=fopen("flag_2.txt","a");
	//if(fp)
	//{
	//	fprintf(fp,"%f	%d	%d	%d	%d	%d	%d	%d	%f	%f\n",
	//		Paras[0],Accurate_Detector_Model_is_ON,
	//		n_col_oversample,n_row_oversample,
	//		n_col_oversample_add_xtalk,n_row_oversample_add_xtalk,
	//		modules.Sub[0],modules.maxSubDets,
	//		modules.Weight[0],modules.Sampling[0]);
	//	fclose(fp);
	//}
	//-------------------


  int i, moduleNumber, k, l, m, n, nListObjects, *objectList = NULL, nSubDets, detIndex, timing = 0;
  double *detWeights, *UV, *sampling, detCenter[3], pix_vlims[2], pix_ulims[2], *center, *right, *up;
  double elapsedTime, min_val;
  clock_t clk;
  int moduleTypeIndex;
  void *boundaries;
  int any_objs;

  if (timing != 0){
    clk = clock();
    elapsedTime = (double)( (clk + 0.0) / CLOCKS_PER_SEC);
  }

  objectList = malloc(phantom.numObjects*sizeof(int));

  switch(modules.moduleOverlapType)
    {
    case 1:
      boundaries = malloc(sizeof(struct height_lims)*phantom.numObjects);
      break;
    case 2:
      boundaries = malloc(sizeof(struct box_lims)*phantom.numObjects);
      break;
    case 3:
      //boundaries = malloc(sizeof(struct polygon)*phantom.numObjects);
      break;
    }

  //dbug(2,"Number of modules: %d     modules.moduleOverlapType : %d\r\n",nModulesIn, modules.moduleOverlapType);
  for(moduleNumber = 0;moduleNumber<nModulesIn;moduleNumber++)
    {
      moduleTypeIndex = modTypeInds[moduleNumber];
      detWeights = &modules.Weight[moduleTypeIndex*modules.maxSubDets];//Weight should be transposed before being passed in
      UV = &modules.Coords[2*modules.maxPixPerModule*moduleTypeIndex];
      //if ((moduleNumber%20)==0) 
      //  dbug(1,"moduleTypeIndex: %d, modules.Pix[moduleTypeIndex]: %d\n\r",moduleTypeIndex,modules.Pix[moduleTypeIndex]);
      sampling = &modules.Sampling[2*modules.maxSubDets*moduleTypeIndex];
      //dbug(2,"Module %d  mti:%d\n\r",moduleNumber,moduleTypeIndex);
      nSubDets = modules.Sub[moduleTypeIndex];
      //dbug(2,"Module %db %d\n\r",moduleNumber,moduleTypeIndex);
      pix_vlims[0] = VERY_BIG;pix_vlims[1] = -VERY_BIG;
      pix_ulims[0] = VERY_BIG;pix_ulims[1] = -VERY_BIG;
      for(i = 0;i<modules.Sub[moduleTypeIndex];i++){
	if(sampling[2*i+1]<pix_vlims[0]) pix_vlims[0] = sampling[2*i+1];
	if(sampling[2*i+1]>pix_vlims[1]) pix_vlims[1] = sampling[2*i+1];
	if(sampling[2*i]<pix_ulims[0]) pix_ulims[0] = sampling[2*i];
	if(sampling[2*i]>pix_ulims[1]) pix_ulims[1] = sampling[2*i];
      }
      //dbug(2,"here %db\n\r",moduleNumber);
      center = &Center[3*moduleNumber];
      right = &Right[3*moduleNumber];
      up = &Up[3*moduleNumber];
      //dbug(2,"Module %db\n\r",moduleNumber);

      compute_object_projections(srcHullPoints, nSrcHullPoints, moduleTypeIndex, up, right, center, boundaries);
      //dbug(2,"Module %dc\n\r",moduleNumber);
      switch(modules.moduleOverlapType)
	{
	case 1:
	  any_objs = any_objects_1(moduleNumber,boundaries);
	  break;
	case 2:
	  any_objs = any_objects_2(moduleNumber,boundaries);
	  break;
	default:
	  printf("\nERROR: Unrecognized moduleOverlapType!\n\r");
	}
      if (any_objs) {
	dbug(1,"There are objects that may project onto module %d.    (%d)\n\r",moduleNumber,modules.moduleOverlapType);
        //        dbug(2,"max : %1.12lf   min : %1.12lf\r\n",boundaries[0].max,boundaries[0].min);
	for(k = 0;k<modules.Pix[moduleTypeIndex];k++){
	  
	  // Compute pixel center locations
	  detIndex = firstDetIndex[moduleNumber]+k;
	  for(i = 0;i<3;i++)
	    detCenter[i] = center[i]+UV[k*2]*right[i]+UV[k*2+1]*up[i];
	  
	  // Build object list
	  nListObjects = 0;
	  switch(modules.moduleOverlapType)
	    {
	    case 1:
	      build_object_list1(pix_vlims, objectList, &nListObjects, moduleNumber, UV[k*2+1], boundaries);
	      //if (nListObjects > phantom.numObjects)
	      //dbug(2,"nListObjects = %d, phantom.numObjects = %d\r\n", nListObjects, phantom.numObjects);
	      //if (nListObjects > 0)
	      //dbug(2,"nListObjects = %d, k = %d\r\n", nListObjects, k);
	      break;
	    case 2:
	      build_object_list2(pix_ulims, pix_vlims, objectList, &nListObjects, moduleNumber, &UV[k*2], boundaries);
	      break;
	    case 3:
	      build_object_list3();
	      break;
	    default:
	      printf("\nERROR: Unrecognized moduleOverlapType!\n\r");
	    }
	  
          //dbug(2,"nListObjects : %d.\n\r",nListObjects);
          //dbug(2,"pix_vlims : %1.5lf %1.5lf.\n\r",pix_vlims[0],pix_vlims[1]);
          //dbug(2,"pix_ulims : %1.5lf %1.5lf.\n\r",pix_ulims[0],pix_ulims[1]);
          //dbug(2,"UV : %1.5lf %1.5lf.\n\r",UV[0],UV[1]);
          //struct box_lims *bl = (struct box_lims *) boundaries;
          //dbug(2,"wd lims : %1.5lf %1.5lf.\n\r",bl[0].minR,bl[0].maxR);
          //dbug(2,"ht lims : %1.5lf %1.5lf.\n\r",bl[0].minU,bl[0].maxU);

	  // Compute relative intensities
	  if(nListObjects){
	    intersections_loop(detCenter, right, up, sampling, nSubDets, sourcePoints, objectList, nListObjects, thisView, detIndex, subviewWeight, detWeights);
	  }
	  else {
	    //	  for(l = 0;l<nSubDets;l++)
	    // We are assuming the detector weights sum to one and that the sourceWeights sum to one also.

		  //------------------- Accurate Detector Model, Mingye
		  if(Accurate_Detector_Model_is_ON)
		  {
			  for(m = 0;m<materials.eBinCount;m++)
				  for(l = 0;l<nSubDets;l++)
				  {
					  int the_index=(detIndex*n_col_oversample+(int)(l/n_row_oversample_add_xtalk))*materials.eBinCount+m;
					  thisView[the_index] += subviewWeight*detWeights[l];
				  }
		  }
		  //-------------------
		  else
		  {
			  for(m = 0;m<materials.eBinCount;m++)
				  thisView[detIndex*materials.eBinCount+m] += subviewWeight;
		  }

	    // thisView[detIndex*materials.eBinCount+m] = thisView[detIndex*materials.eBinCount+m]+subviewWeight[n]*detWeights[l]*sourceWeights;
	  }   // else
	}    // Loop over pixels (not subpixels)
      }    // if any_objects
      else {
	//      sum_det_w = 0;                  
	//      for(l = 0;l<nSubDets;l++) sum_det_w += detWeights[l];
	// We are assuming the detector weights sum to one and that the sourceWeights sum to one also.
	for(k = 0;k<modules.Pix[moduleTypeIndex];k++){
	  detIndex = firstDetIndex[moduleNumber]+k;

	  //------------------- Accurate Detector Model, Mingye
	  if(Accurate_Detector_Model_is_ON)
	  {
		  for(m = 0;m<materials.eBinCount;m++)
			  for(l = 0;l<nSubDets;l++)
			  {
				  int the_index=(detIndex*n_col_oversample+(int)(l/n_row_oversample_add_xtalk))*materials.eBinCount+m;
				  thisView[the_index] += subviewWeight*detWeights[l];
			  }
	  }
	  //-------------------
	  else
	  {
		  for(m = 0;m<materials.eBinCount;m++)
			  thisView[detIndex*materials.eBinCount+m] += subviewWeight;
	  }

	}
      }
      if (debug_flag>0)
	{
	  min_val = VERY_BIG;
	  for(k = 0;k<modules.Pix[moduleTypeIndex];k++){
	    detIndex = firstDetIndex[moduleNumber]+k;
	    for(m = 0;m<materials.eBinCount;m++)
	      if (thisView[detIndex*materials.eBinCount+m]<min_val) min_val = thisView[detIndex*materials.eBinCount+m];
	  }
	  dbug(1,"minval: %f (module #%d)\n\r", min_val,moduleNumber);
	}
      //dbug(2,"Module %dd\n\r",moduleNumber);
    }
  
  // Print timing, if requested
  if (timing != 0)
    {
      clk = clock();
      elapsedTime = (double)( (clk + 0.0) / CLOCKS_PER_SEC)-elapsedTime;
      printf("End module loop time = %1.12f\r\n",elapsedTime);
    }

  //dbug(2,"freebefore\n\r",UV[0],UV[1]);
  free(objectList);objectList = NULL;
  free(boundaries);boundaries = NULL;
  //dbug(2,"freeafter\n\r",UV[0],UV[1]);
}

void *projector_wrapper(void *pointerIn)

{
  struct projector_args *projectorArgs;
  FILE *fid;
  int valuesRead;
  int goToSleep;
  int module;
  int notFirstTimeThrough;
  double oneMinuteLoadAverage = 0;
  double loadAdjust;
  pthread_t my_id = pthread_self();
  int threadNumber = -1;
  int i;

  for (i=0;i<thread_count;i++)
    if (pthread_equal(my_id,t_id[i]))
      threadNumber = i;
  loadAdjust = (threadNumber - thread_count/2 + 0.5) * 0.1;
  
  //printf("Start thread %d\r\n", threadNumber);
  projectorArgs = (struct projector_args *) pointerIn;

  while (nextModuleInQ < modulesInQ)
    {
      goToSleep = 0;
#ifdef LIMIT_LOAD
      if (limitLoadAverage != 0)
	{
	  fid = fopen("/proc/loadavg","r");
	  valuesRead = fscanf(fid,"%lf",&oneMinuteLoadAverage);
	  fclose(fid);
	  if (valuesRead == 0)
	    goToSleep = 1;
	  //else
	  //  printf(" One Minute Load Average = %1.2lf   limit: %1.2lf \n\r",oneMinuteLoadAverage,limitLoadAverage);
	  if (oneMinuteLoadAverage > limitLoadAverage + loadAdjust)
	    goToSleep = 1;
	}
      if (goToSleep)
	{dbug(2,"tic\n");_sleep(6);dbug(2,"toc\n");}
      else
	{
#endif
	  //printf("Thread #%d is tyring to get a module.\n\r",threadNumber);
	  pthread_mutex_lock(&QLock);
	  module = nextModuleInQ;
          nextModuleInQ++;
	  pthread_mutex_unlock(&QLock);

	  if (module < modulesInQ)
	    {
	      dbug(1,"Thread #%d acquired module #%d.\n\r",threadNumber,module);
	      //if ((module == 32)) debug_flag = 2; else debug_flag = 1;
	      Projector(projectorArgs[0].Paras, projectorArgs[0].subviewWeight, projectorArgs[0].thisView, projectorArgs[0].sourcePoints, projectorArgs[0].nSubSources, projectorArgs[0].srcHullPoints, projectorArgs[0].nSrcHullPoints, &(projectorArgs[0].firstDetIndex[module]), 1, &(projectorArgs[0].modTypeInds[module]), &(projectorArgs[0].Up[3*module]), &(projectorArgs[0].Right[3*module]), &(projectorArgs[0].Center[3*module]), projectorArgs[0].UNUSED_tvLength);
	      dbug(1,"Thread #%d finished module #%d.\n\r",threadNumber,module);
	    }
#ifdef LIMIT_LOAD 
	}
#endif
    }
}

DLLEXPORT void Projector_threaded(double *Paras, double subviewWeight, double *thisView, double *sourcePoints, int nSubSources, double *srcHullPoints, int nSrcHullPoints, int *firstDetIndex, int nModulesIn, int *modTypeInds, double *Up, double *Right, double *Center, int UNUSED_tvLength, int numThreads, double maxLoadAverage)

{
	//------------------- Accurate Detector Model, Mingye
	if(Paras[0]==1)
	{
		set_Accurate_Detector_Model(Paras);
		Paras[0]=-1;
	}
	else if(Paras[0]==0)
		Accurate_Detector_Model_is_ON=0;

	//FILE *fp;
	//fp=fopen("flag_1.txt","a");
	//if(fp)
	//{
	//	fprintf(fp,"%f	%d	%d	%d	%d	%d	%d	%d	%f	%f\n",
	//		Paras[0],Accurate_Detector_Model_is_ON,
	//		n_col_oversample,n_row_oversample,
	//		n_col_oversample_add_xtalk,n_row_oversample_add_xtalk,
	//		modules.Sub[0],modules.maxSubDets,
	//		modules.Weight[0],modules.Sampling[0]);
	//	fclose(fp);
	//}
	//-------------------


  int i;
  struct projector_args projectorArgs[1];

  //printf("Point T0\n\r");
  thread_count = numThreads;
  limitLoadAverage = maxLoadAverage;
  
  projectorArgs[0].Paras = Paras; //parameters for Accurate Detector Model
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

  //dbug(2,"Point T1\n\r");
  // Create the threads
  t_id = malloc(sizeof(pthread_t)*thread_count);
  for (i=0;i<thread_count;i++)
    {
      pthread_create(&t_id[i],NULL,projector_wrapper,projectorArgs);
      //printf("Created thread %d",i);
    }
  //dbug(2,"Point T2\n\r");
  // Wait for them to complete
  for (i=0;i<thread_count;i++)
    {
      pthread_join(t_id[i],NULL);
      //printf("Finished thread %d",i);
    }
  free(t_id);
}

void diff(double *d1,double *d2,int len)
{
  int i;
  double nrm = 0;
  for(i = 0;i<len;i++)
    nrm += (d1[i]-d2[i])*(d1[i]-d2[i]);
  printf("%.10f  ",nrm);
}


void CO(double *c_values)

{
  int n,i,imax = -1;
  double *dat,mx_diff = 0,nrm_diff = 0,tmp;
  FILE *fid;
  char *fname;

  if (debug_flag){
  COMPARISON_NUMBER++;
  fname = malloc(sizeof(char)*1280);
  sprintf(fname,"CO");
  sprintf(&fname[2],"%d",COMPARISON_NUMBER);
  sprintf(&fname[3],".dat");
  fid = fopen(fname,"r");
  fread(&n,sizeof(int),1,fid);
  printf("\r\n%d\r\n\n\n\n\n\n",n);
  dat = malloc(sizeof(double)*n);
  fread(dat,sizeof(double),n,fid);
  fclose(fid);
  for(i = 0;i<n;i++)
    {
      tmp = fabs(c_values[i]-dat[i]);
      if (tmp>mx_diff) {mx_diff = tmp;imax = i;}
      nrm_diff += tmp*tmp;
    }
  printf("\r\n\n Comparison number %d:\r\n   mx_diff  :  %1.12lf\r\n   (at index #%d)\r\n   nrm_diff :  %1.12lf\r\n\n",COMPARISON_NUMBER,mx_diff,imax,nrm_diff);
  free(dat);
  }
}

void COf(float *c_values)

{
  int n,i,imax = -1;
  float *dat,mx_diff = 0,nrm_diff = 0,tmp;
  FILE *fid;
  char *fname;

  if (debug_flag){
  COMPARISON_NUMBER++;
  fname = malloc(sizeof(char)*1280);
  sprintf(fname,"CO");
  sprintf(&fname[2],"%d",COMPARISON_NUMBER);
  sprintf(&fname[3],".dat");
  fid = fopen(fname,"r");
  fread(&n,sizeof(int),1,fid);
  printf("\r\n%d\r\n\n\n\n\n\n",n);
  dat = malloc(sizeof(float)*n);
  fread(dat,sizeof(float),n,fid);
  fclose(fid);
  for(i = 0;i<n;i++)
    {
      tmp = fabs(c_values[i]-dat[i]);
      if (tmp>mx_diff) {mx_diff = tmp;imax = i;}
      nrm_diff += tmp*tmp;
    }
  printf("\r\n\n Comparison number %d:\r\n   mx_diff  :  %1.12f\r\n   (at index #%d)\r\n   nrm_diff :  %1.12f\r\n\n",COMPARISON_NUMBER,mx_diff,imax,nrm_diff);
  free(dat);
  }
}

void COi(int *c_values)

{
  int n,i,imax = -1;
  int *dat,mx_diff = 0,nrm_diff = 0,tmp;
  FILE *fid;
  char *fname;

  if (debug_flag){
  COMPARISON_NUMBER++;
  fname = malloc(sizeof(char)*1280);
  sprintf(fname,"CO");
  sprintf(&fname[2],"%d",COMPARISON_NUMBER);
  sprintf(&fname[3],".dat");
  fid = fopen(fname,"r");
  fread(&n,sizeof(int),1,fid);
  printf("\r\n%d\r\n\n\n\n\n\n",n);
  dat = malloc(sizeof(int)*n);
  fread(dat,sizeof(int),n,fid);
  fclose(fid);
  if (n<32){
    printf("\r\n\n Comparison number %d:\r\n\n",COMPARISON_NUMBER);
    for(i = 0;i<n;i++)
      {
	tmp = abs(c_values[i]-dat[i]);
	printf("FreeMat: %d   C: %d\r\n",dat[i],c_values[i]);
      }
  }
  else {
    for(i = 0;i<n;i++)
      {
	tmp = abs(c_values[i]-dat[i]);
	if (tmp>mx_diff) {mx_diff = tmp;imax = i;}
	nrm_diff += tmp*tmp;
      }
    printf("\r\n\n Comparison number %d:\r\n   mx_diff  :  %d\r\n   (at index #%d)\r\n   nrm_diff :  %d\r\n\n",COMPARISON_NUMBER,mx_diff,imax,nrm_diff);
  }
  free(dat);
  }
}

void COs(short *c_values)

{
  int n,i,imax = -1;
  short *dat;
  int mx_diff = 0,nrm_diff = 0,tmp;
  FILE *fid;
  char *fname;

  if (debug_flag){
  COMPARISON_NUMBER++;
  fname = malloc(sizeof(char)*1280);
  sprintf(fname,"CO");
  sprintf(&fname[2],"%d",COMPARISON_NUMBER);
  sprintf(&fname[3],".dat");
  fid = fopen(fname,"r");
  fread(&n,sizeof(int),1,fid);
  printf("\r\n%d\r\n\n\n\n\n\n",n);
  dat = malloc(sizeof(short)*n);
  fread(dat,sizeof(short),n,fid);
  fclose(fid);
  for(i = 0;i<n;i++)
    {
      tmp = abs(c_values[i]-dat[i]);
      if (tmp>mx_diff) {mx_diff = tmp;imax = i;}
      nrm_diff += tmp*tmp;
    }
  printf("\r\n\n Comparison number %d:\r\n   mx_diff  :  %d\r\n   (at index #%d)\r\n   nrm_diff :  %d\r\n\n",COMPARISON_NUMBER,mx_diff,imax,nrm_diff);
  free(dat);
  }
}

// gcc -shared -fPIC -O3 -mfpmath=sse -msse -o YYcode64.so YYcode.c
// gcc -shared -fPIC -O3 -mfpmath=sse -msse -m32 -o YYcode.so YYcode.c


