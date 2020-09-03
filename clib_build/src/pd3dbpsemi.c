// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h>

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define SAT(t,tmax) (MAX(MIN(t,tmax-1),0))
#define LOOKUP(u,v) projData[SAT(v,detVcount)*detScount+SAT(u,detScount)]

/*
 * Backproject a view onto a 2D image using a 3D imaging geometry that corresponds
 * to a fan-rebinned view from a cone beam system. =P
 */
#ifdef WIN32
__declspec(dllexport)
#endif
void pd3dbpsemi(float *projData,  /* View data, should be detScount x detVcount */
		float *imgData,   /* Image data - pass by reference, should be pixCount x pixCount */
		int pixCount,     /* Number of image pixels (linear count) */
		float xcenter,    /* Center of image in x coordinate */
		float ycenter,    /* Center of image in y coordinate */
		float pixSize,    /* Size of pixels in mm */
		int detScount,    /* Number of detectors in S direction (parallel beam) */
		int detVcount,    /* Number of detectors in V direction (row direction) */
		float SID,        /* Source to isocenter distance */
		float SDD,        /* Source to detector distance */
		float detVsize,   /* Size of detector pixels */
		float detVcenter, /* Index (fractional) of center row */
		float sMin,       /* Smallest s value */
		float sDelt,      /* Delta s */
		float beta,       /* Rotation angle for this view */
		float zoffset,    /* Distance (in Z) between current slice and view axis */
		float pitch) {    /* Helix pitch, expressed in mm/rotation */
  int row, col;
  float xcoord, ycoord;
  float cosbeta, sinbeta;
  float tcoord, scoord, pixHalf;
  float delbeta, delZ, gamma, Zdet;
  float normalized_pitch;
  float detZmin;
  float ureal, vreal;
  float puv, pu1v, puv1, pu1v1;
  int uind, vind;
  float ueps, veps;
  float pustarv, pustarv1;
  
  cosbeta = cos(beta);
  sinbeta = sin(beta);
  pixHalf = (pixCount-1.0)/2.0;
  normalized_pitch = pitch/(2.0*M_PI);
  detZmin = -detVcenter*detVsize;
  for (row=0;row<pixCount;row++) {
    ycoord = (row-pixHalf)*pixSize + ycenter;
    for (col=0;col<pixCount;col++) {
      xcoord = (col-pixHalf)*pixSize + xcenter;
      /* Map (x,y) location to position along parallel-beam axis */
      scoord = xcoord*cosbeta + ycoord*sinbeta;
      tcoord = -xcoord*sinbeta + ycoord*cosbeta;
      /* 
	 Calculate how much the source had to rotate (relative to current)
         to line up with the current point (along the parallel ray axis)
      */
      delbeta = asin(scoord/SID);
      /* Calculate how much the source moved in z with this rotation */
      delZ = delbeta * normalized_pitch;
      /* Calculate the distance from source to iso as seen along the parallel ray axis */
      gamma = SID*cos(delbeta);
      /* Adjust delZ for the magnification to the detector */
      Zdet = (zoffset+delZ)*SDD/(gamma-tcoord);
      /* map scoord to a u coordinate */
      ureal = (scoord-sMin)/sDelt;
      /* map Zdet to a v coordinate */
      vreal = (Zdet-detZmin)/detVsize;
      uind = (int) floor(ureal);
      vind = (int) floor(vreal);
      /* retrieve the four values */
      puv = LOOKUP(uind,vind);
      pu1v = LOOKUP(uind+1,vind);
      puv1 = LOOKUP(uind,vind+1);
      pu1v1 = LOOKUP(uind+1,vind+1);
      /* linearly interpolate in the u directions */
      ueps = ureal - uind;
      pustarv = puv + (ueps)*(pu1v-puv);
      pustarv1 = puv1 + (ueps)*(pu1v1-puv1);
      veps = vreal - vind;
      /* Update the pixel using bilinear interpolation */
      imgData[row*pixCount+col] +=  pustarv + (veps)*(pustarv1-pustarv);
    }
  }
}
		
		
		
#ifdef WIN32
__declspec(dllexport)
#endif
void pd3dbpsemi_flat(float *projData,  /* View data, should be detScount x detVcount */
		     float *imgData,   /* Image data - pass by reference, should be pixCount x pixCount */
		     int pixCount,     /* Number of image pixels (linear count) */
		     float xcenter,    /* Center of image in x coordinate */
		     float ycenter,    /* Center of image in y coordinate */
		     float pixSize,    /* Size of pixels in mm */
		     int detScount,    /* Number of detectors in S direction (parallel beam) */
		     int detVcount,    /* Number of detectors in V direction (row direction) */
		     float SID,        /* Source to isocenter distance */
		     float SDD,        /* Source to detector distance */
		     float detVsize,   /* Size of detector pixels */
		     float detVcenter, /* Index (fractional) of center row */
		     float sMin,       /* Smallest s value */
		     float sDelt,      /* Delta s */
		     float beta,       /* Rotation angle for this view */
		     float zoffset,    /* Distance (in Z) between current slice and view axis */
		     float pitch) {    /* Helix pitch, expressed in mm/rotation */
  int row, col;
  float xcoord, ycoord;
  float cosbeta, sinbeta;
  float tcoord, scoord, pixHalf;
  float delbeta, delZ, gamma, Zdet;
  float normalized_pitch;
  float detZmin;
  float ureal, vreal;
  float puv, pu1v, puv1, pu1v1;
  int uind, vind;
  float ueps, veps;
  float pustarv, pustarv1;
  float H;
  
  cosbeta = cos(beta);
  sinbeta = sin(beta);
  pixHalf = (pixCount-1.0)/2.0;
  normalized_pitch = pitch/(2.0*M_PI);
  detZmin = -detVcenter*detVsize;
  for (row=0;row<pixCount;row++) {
    ycoord = (row-pixHalf)*pixSize + ycenter;
    for (col=0;col<pixCount;col++) {
      xcoord = (col-pixHalf)*pixSize + xcenter;
      /* Map (x,y) location to position along parallel-beam axis */
      scoord = xcoord*cosbeta + ycoord*sinbeta;
      tcoord = -xcoord*sinbeta + ycoord*cosbeta;
      /* 
	 Calculate how much the source had to rotate (relative to current)
         to line up with the current point (along the parallel ray axis)
      */
      delbeta = asin(scoord/SID);
      /* Calculate how much the source moved in z with this rotation */
      delZ = delbeta * normalized_pitch;
      /* Calculate the distance from source to iso as seen along the parallel ray axis */
      gamma = SID*cos(delbeta);
      H = SDD/cos(delbeta);
      /* Adjust delZ for the magnification to the detector */
      Zdet = (zoffset+delZ)*H/(gamma-tcoord);
      /* map scoord to a u coordinate */
      ureal = (scoord-sMin)/sDelt;
      /* map Zdet to a v coordinate */
      vreal = (Zdet-detZmin)/detVsize;
      uind = (int) floor(ureal);
      vind = (int) floor(vreal);
      /* retrieve the four values */
      puv = LOOKUP(uind,vind);
      pu1v = LOOKUP(uind+1,vind);
      puv1 = LOOKUP(uind,vind+1);
      pu1v1 = LOOKUP(uind+1,vind+1);
      /* linearly interpolate in the u directions */
      ueps = ureal - uind;
      pustarv = puv + (ueps)*(pu1v-puv);
      pustarv1 = puv1 + (ueps)*(pu1v1-puv1);
      veps = vreal - vind;
      /* Update the pixel using bilinear interpolation */
      imgData[row*pixCount+col] +=  pustarv + (veps)*(pustarv1-pustarv);
    }
  }
}
		
		
		


