// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h>
#include <stdlib.h> // calloc(), malloc(), free()

void DD2ParProjInner(int direction,
		     int nrcols,
		     int nrrows,
		     float *detDistCopy,
		     float previousDist,
		     float *projCopy,
		     float *imgCopy,
		     float stepSize)

{
  static int colnr;
  static float pixDist;

  colnr=0;
  pixDist=previousDist+stepSize;

  while (colnr <= (nrcols-1))
    { 
      if (pixDist <= *detDistCopy) /* next limit is pixel border */
	{
	  *projCopy += (pixDist-previousDist)*(*imgCopy);
	  imgCopy += nrrows;
	  previousDist=pixDist;
	  pixDist+=stepSize;
	  colnr++;
	}
      else /* next limit is detector border */
	{
	  *projCopy += (*detDistCopy-previousDist)*(*imgCopy);
	  previousDist=*detDistCopy;
	  detDistCopy += direction;
	  projCopy += direction;
	}
    }
  
}

//-----------------------------------------------------------------------------

/*
 * Distance-driven 2D Fan-beam view projector
 */
void DD2ParProjView(int nrdet,
		    int vertical,
		    float *xdis,
		    float *ydis,
		    float *detDist,
		    float *projection,
		    float *newproj,
		    int nrcols,
		    int nrrows,
		    float* originalImgPtr,
		    float* transposeImgPtr,
                    float sinAngle,
                    float cosAngle)
{
  float *detDistCopy, *imgCopy, *imgCCopy, *newprojCopy;
  int dummyint, rownr, distnr, detnr;
  float previousDist, stepSize, dummyfloat;
  float RowShift;

  /*
   * Calculate detector distances
   */
  detDistCopy=detDist;
  if (vertical) /* vertical projection */    
    {
    RowShift = sinAngle/cosAngle;
    *detDistCopy++=1.e12; /* sentinel */
    for (distnr=0 ; distnr <= nrdet ; distnr++)
      {
       *detDistCopy++ = *xdis++ + *ydis++ * RowShift; 
      }
    *detDistCopy=1.e12; /* sentinel */
    imgCopy=originalImgPtr;
    }
  else /* horizontal projection */
    {
    RowShift =  -cosAngle/sinAngle;
    *detDistCopy++=1.e12; /* sentinel */
    for (distnr=0 ; distnr <= (nrdet) ; distnr++)
      {
       *detDistCopy++ = *ydis++ - *xdis++ * RowShift;
      }
    *detDistCopy=1.e12; /* sentinel */
    imgCopy=transposeImgPtr;
    dummyint=nrrows;
    nrrows=nrcols;
    nrcols=dummyint;
    }


  /*
   * Reset newproj
   */
  newprojCopy=newproj+1;
  for (detnr=0 ; detnr <= nrdet-1 ; detnr++)
    {
      *newprojCopy++=0.;
    }

  /*
   * If first detector has lowest distance value
   */  
    if (*(detDist+2) > *(detDist+1))
      {
        *detDist=-1.e12; /* invert sentinel */

        /*
         * For all rows
         */
        for (rownr=0 ; rownr <= (nrrows-1) ; rownr++)
    	  {
 	  /*
           * Initialize image parameters
           */
           previousDist=(-nrcols/2.) - (rownr-(nrrows-1)/2.) * RowShift;
          /*
           * Initialize measurement parameters
           */
           distnr=0;
           detDistCopy=detDist+1;
          while (*detDistCopy <= previousDist)
	    {
             distnr++;
             detDistCopy++;
            } /* now *detDistCopy > previousDist */
	  newprojCopy=newproj+distnr;

	  /*
	   * For all distances
	   */
          imgCCopy = imgCopy;
	  DD2ParProjInner(1, nrcols, nrrows, detDistCopy, 
                          previousDist, newprojCopy, imgCCopy, 1);
          imgCopy++; 

 	}
     }
  /*
   * If last detector has lowest distance value
   */  
   else
     {
     *(detDist+nrdet+2)=-1.e12; /* invert sentinel */
      /*
       * For all rows
       */
      for (rownr=0 ; rownr <= (nrrows-1) ; rownr++)
	{
	  /*
	   * Initialize image parameters
	   */
	  previousDist=(-nrcols/2.) - (rownr-(nrrows-1)/2.) * RowShift;
	  /*
	   * Initialize measurement parameters
	   */
	  distnr=nrdet;
	  detDistCopy=detDist+nrdet+1;
	  while (*detDistCopy <= previousDist)
	    {
	      distnr--;
	      detDistCopy--;
	    } /* now *detDistCopy > previousDist */
	  newprojCopy=newproj+distnr+1;

	  /*
	   * For all distances
	   */
          imgCCopy = imgCopy;
	  DD2ParProjInner(-1, nrcols, nrrows, detDistCopy, 
                          previousDist, newprojCopy, imgCCopy, 1);
          imgCopy++; 
	}
    }
  /*
   * Scale projections
   */
  newprojCopy=newproj+1;
  for (detnr=0 ; detnr <= nrdet-1 ; detnr++)
    {
      if (*newprojCopy != 0)
	{
	  *projection++ = *newprojCopy++;
	} 
      else
	{
	  *newprojCopy++ = 0.;
	  projection++;
	}
    }
}

//-----------------------------------------------------------------------------

/*
 * Distance-driven 2D parallel beam projector
 */
#ifdef WIN32
__declspec(dllexport)
#endif
void DD2ParProj(int nrdet,
		float *xds,
		float *yds,
		float xCor,
		float yCor,
		float *viewangles,
		int nrviews,
		float *sinogram,
		int nrcols,
		int nrrows,
		float *originalImgPtr)
{
  float *xdsCopy, *ydsCopy, *newproj;
  float *xdi, *ydi, *xdiCopy, *ydiCopy;
  float *xdiRot, *ydiRot, *xdiRotCopy, *ydiRotCopy;
  float angle, sinAngle, cosAngle;
  float *sinogramCopy, *distances, *viewanglesCopy;
  float *originalImgPtrCopy, *transposeImgPtr, *transposeImgPtrCopy;
  float *rotateImgPtr, *rotateImgPtrCopy;
  int colnr, rownr, distnr, viewnr, vertical;


  /*
   * Allocate memory for rotated coordinates
   */
  xdi=(float*)malloc((nrdet+1)*sizeof(float));
  ydi=(float*)malloc((nrdet+1)*sizeof(float));
  xdiRot=(float*)malloc((nrdet+1)*sizeof(float));
  ydiRot=(float*)malloc((nrdet+1)*sizeof(float));
  newproj=(float*)calloc((nrdet+2),sizeof(float));

  /*
   * Create transpose image
   */
  transposeImgPtr=(float*)malloc(nrcols*nrrows*sizeof(float));
  transposeImgPtrCopy = transposeImgPtr;
  originalImgPtrCopy  = originalImgPtr;
  for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
    {
    for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
      {
      *transposeImgPtrCopy = *originalImgPtrCopy++;
      transposeImgPtrCopy += nrcols;
      }
    transposeImgPtrCopy += (1-nrrows*nrcols);
    }


  /* 
   * Flip transpose image to get rotated image
   */
  rotateImgPtr=(float*)malloc(nrcols*nrrows*sizeof(float));
  rotateImgPtrCopy     = rotateImgPtr;
  transposeImgPtrCopy  = transposeImgPtr   ;
  rotateImgPtrCopy    += (nrrows-1)*nrcols ;
  for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
    {
    for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
      {
      *rotateImgPtrCopy++ = *transposeImgPtrCopy++;
      }
    rotateImgPtrCopy -= 2*nrcols ;
    }
  free(transposeImgPtr);


  /*
   * Load rotated image in original image
   */
//      rotateImgPtrCopy   = rotateImgPtr;
//      originalImgPtrCopy = originalImgPtr;
//      for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
//        {
//        for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
//          {
//          *originalImgPtrCopy++ = *rotateImgPtrCopy++;
//          }
//        }

  /*
   * Load rotated image in original image
   */
      //  transposeImgPtrCopy   = transposeImgPtr;
//  originalImgPtrCopy = originalImgPtr;
//  for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
//    {
//    for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
//      {
//      *originalImgPtrCopy++ = *transposeImgPtrCopy++;
//      }
//    }


  /*
   * Calculate detector boundaries
   */
  *xdi = 1.5 * *xds - 0.5 * *(xds+1);
  *ydi = 1.5 * *yds - 0.5 * *(yds+1);
  xdsCopy=xds;
  ydsCopy=yds;
  xdiCopy=xdi+1;
  ydiCopy=ydi+1;
  for (distnr=1 ; distnr<=(nrdet-1) ; distnr++)
    {
      *xdiCopy++ = 0.5 * *xdsCopy + 0.5 * *(xdsCopy+1);
      xdsCopy++;
      *ydiCopy++ = 0.5 * *ydsCopy + 0.5 * *(ydsCopy+1);
      ydsCopy++;
    }
  *xdiCopy = 1.5 * *xdsCopy - 0.5 * *(xdsCopy-1);
  *ydiCopy = 1.5 * *ydsCopy - 0.5 * *(ydsCopy-1);

  
  /*
   * Prepare empty array to contain distances
   */
  distances=(float*)malloc((nrdet+3)*sizeof(float)); /* provide 2 spaces
                                                     for sentinels */


  /*
   * Loop over all views
   */
  sinogramCopy=sinogram;
  viewanglesCopy=viewangles;
  angle = *viewanglesCopy++;
  for (viewnr = 0 ; viewnr <= (nrviews-1) ; viewnr++)
    {
      /*
       * Rotate coordinates
       */
      sinAngle=(float)sin(angle);
      cosAngle=(float)cos(angle);
      xdiRotCopy=xdiRot;
      ydiRotCopy=ydiRot;
      xdiCopy=xdi;
      ydiCopy=ydi;
      for (distnr=0 ; distnr<=nrdet ; distnr++)
	{
	  *xdiRotCopy++ = (*xdiCopy-xCor)*cosAngle
	    - (*ydiCopy-yCor)*sinAngle + xCor;
	  *ydiRotCopy++ = (*ydiCopy++ - yCor)*cosAngle
	    + (*xdiCopy++-xCor)*sinAngle + yCor;
	}
      vertical = (fabs(cosAngle) >= fabs(sinAngle)) ;

      /*
       * View projection
       */
      DD2ParProjView(nrdet, vertical, xdiRot, ydiRot,
                     distances, sinogramCopy, newproj, 
                     nrcols, nrrows, originalImgPtr, 
                     rotateImgPtr,sinAngle,cosAngle);
      //                     transposeImgPtr,sinAngle,cosAngle);
      sinogramCopy += nrdet;
      angle = *viewanglesCopy++;
    }

  /*
   * Clean up memory
   */
  free(xdiRot);
  free(ydiRot);
  free(xdi);
  free(ydi);
  free(rotateImgPtr);
  free(distances);
  free(newproj);
}


