// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h>

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

static void innerLoop(int N, float initialC, float deltaC,
		      float initialB, float deltaB,
		      float recipAlpha, float detUCenter,
		      float *imgPtr, float *viewPtr,
		      int detUCount) {
  int n;
  float alphaval;
  float alpha_rescaled, alpha_residual;
  int alpha_discrete;
  float weight_val;
  float *dp;
  float c, b;
  float toadd;
  
  dp = imgPtr;
  c = initialC;
  b = initialB;
  
  for (n=0;n<N;n++) {
    alphaval = atan(c/b);
    alpha_rescaled = alphaval*recipAlpha+detUCenter;
    alpha_discrete = ((int) (alpha_rescaled+10)) - 10;
    alpha_residual = alpha_rescaled - alpha_discrete;
    weight_val = 1.0/(b*b+c*c);
    toadd = 0.0;
    if ((alpha_discrete >= 0) && (alpha_discrete < detUCount-1))
      toadd += viewPtr[alpha_discrete]*(1-alpha_residual);
    if ((alpha_discrete >= -1) && (alpha_discrete < detUCount-2))
      toadd += viewPtr[alpha_discrete+1]*alpha_residual;
    dp[n] += toadd*weight_val;
    b += deltaB;
    c += deltaC;
  }
}

#ifdef WIN32
__declspec(dllexport)
#endif
void pd2dbpcurved(float detUCenter, int detUCount, float detUSize,
		  float pixelXSize, float pixelYSize,
		  float pixelXCenter, float pixelYCenter,
		  int pixelXCount, int pixelYCount,
		  float SID, float SDD,
		  int viewCount, float *viewAngles,
		  float *sino, float *image) {
  int i, k;
  float costheta, sintheta;
  float theta;
  float xVal, yVal;
  float recipAlpha;
  float initialB, initialC;
  float deltaB, deltaC;

  /*
    We also need to calculate 1/(detector_size) in inverse radians.
  */
  recipAlpha = 1.0/atan(detUSize/SDD);
  for (i=0;i<viewCount;i++) {
    theta = viewAngles[i];
    /**
       Calculate the cosine and sine of this angle
       for use.
    */
    costheta = cos(theta);
    sintheta = sin(theta);
    for (k=0;k<pixelYCount;k++) {
      /**
	 For this row of the image, we calculate the \f$(x,y)\f$ position
	 of the first pixel in the row.
      */
      xVal = pixelXCenter - (pixelXCount-1)/2.0*pixelXSize;
      yVal = pixelYCenter - ((pixelYCount-1)/2.0-k)*pixelYSize;
      /*
	We map this to a modified coordinate system (explanation forthcoming...)
      */
      initialB = SID - xVal*costheta - yVal*sintheta;
      initialC = -xVal*sintheta + yVal*costheta;
      deltaB = -costheta*pixelXSize;
      deltaC = -sintheta*pixelYSize;
      innerLoop(pixelXCount, initialC, deltaC,
		initialB, deltaB, recipAlpha,
		detUCenter, image + k*pixelXCount, 
		sino + i*detUCount, detUCount);
    }
  }
  for (i=0;i<pixelXCount*pixelYCount;i++)
    image[i] *= 2.0*M_PI/viewCount;
}


