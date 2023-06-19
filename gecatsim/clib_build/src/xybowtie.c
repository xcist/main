// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 *
 * Routine to calculate pathlengths through a closed object defined by a set of (x,y) points
 * that is assumed to have an even number of intersections with any ray omitted from
 * the origin.
 *
 */

int flt_compare(const void *a, const void *b) {
  const float *A, *B;
  A = (const float *) a;
  B = (const float *) b;
  if ((*A) < (*B)) return -1;
  if ((*A) > (*B)) return 1;
  return 0;
}

#ifdef WIN32
__declspec(dllexport)
#endif
void xybowtie(int npoints, float *xvec, float *yvec, float source_x, float source_y, int nfanangles, float *fanangles, float *pathlengths) {
  int i;
  int j;
  float rhsx, rhsy;
  float m00, m01, m10, m11;
  float fan;
  float deltax, deltay;
  float *intersections;
  int nintersections;
  float beta, alpha;

  intersections = (float*) malloc(npoints+1);
  for (i=0;i<nfanangles;i++) {
    /* retrieve the fan angle */
    fan = fanangles[i];
    /* expand into a direction vector, fa=0 --> -y unit vector */
    deltax = sin(fan);
    deltay = -cos(fan);
    m00 = deltax;
    m10 = deltay;
    nintersections = 0;
    for (j=0;j<npoints;j++) {
      float dx, dy;
      dx = xvec[(j+1) % npoints] - xvec[j];
      dy = yvec[(j+1) % npoints] - yvec[j];
      m01 = dx;
      m11 = dy;
      rhsx = source_x - xvec[j];
      rhsy = source_y - yvec[j];
      /* Solve for alpha */
      alpha = (m00*rhsy - m10*rhsx)/(m00*m11-m10*m01);
      if ((alpha >=0) && (alpha < 1)) {
	beta = (m11*rhsx - m01*rhsy)/(m00*m11-m10*m01);
	intersections[nintersections++] = beta;
      }
    }
    /* sort the intersection times */
    if (nintersections & 1) {
      fprintf(stderr,"Warning!  Odd number of intersections encountered!\n");
      nintersections--;
    }
    qsort(intersections,nintersections,sizeof(float),flt_compare);
    pathlengths[i] = 0;
    for (j=0;j<nintersections/2;j++)
      pathlengths[i] += (intersections[2*j+1] - intersections[2*j]);
  }
  free(intersections);
}


