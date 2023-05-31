// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

float* vector(int nl,int nh) {
  float *v;
  
  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  return v-nl;
}

void free_vector(float *v, int nl, int nh) {
  free((char*) (v+nl));
}


void nr_spline(float *x,float *y, int n, float yp1, float ypn, float *y2) {
  int i,k;
  float p,qn,sig,un,*u;
  u=vector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,1,n-1);
}

void nr_splint(float *xa,float *ya,float *y2a,int n, float x, float *y) {
  int klo,khi,k;
  float h,b,a;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) {
    printf("Bad XA input to routine SPLINT");
    assert(0);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#ifdef WIN32
__declspec(dllexport)
#endif
void spline_interpolate(float *xold, float *yold, int nold,
			float *xnew, float *ynew, int nnew) {
  int i;
  float *ytmp = (float*) malloc(sizeof(float)*nold);
  nr_spline(xold-1,yold-1,nold,1E30,1E30,ytmp-1);
  for (i=0;i<nnew;i++) 
    nr_splint(xold-1,yold-1,ytmp-1,nold,xnew[i],ynew+i);
}


