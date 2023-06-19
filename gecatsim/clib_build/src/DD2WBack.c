// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h>
#include <stdlib.h>

void DD2FanWBackInner(int direction,
		      int nrcols,
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
	  *imgCopy += stepSize*(pixDist-previousDist)*(*projCopy);
	  imgCopy++;
	  previousDist=pixDist;
	  pixDist+=stepSize;
	  colnr++;
	}
      else /* next limit is detector border */
	{
	  *imgCopy += stepSize*(*detDistCopy-previousDist)*(*projCopy);
	  previousDist=*detDistCopy;
	  detDistCopy += direction;
	  projCopy += direction;
	}
    }
  
}

//-----------------------------------------------------------------------------

/*
 * Distance-driven 2D Fan-beam view WeightedBackprojector
 */
void DD2FanWBackView(int nrdet,
		     int vertical,
		     float x0,
		     float y0,
		     float *xdis,
		     float *ydis,
		     float *detDist,
		     float *projection,
		     float *newproj,
		     int nrcols,
		     int nrrows,
		     float* originalImgPtr,
		     float* transposeImgPtr,
		     float sid, 
		     int mode) // 1=flat, 0=curved
{
  float *detDistCopy, *imgCopy, *newprojCopy;
  int dummyint, rownr, distnr, detnr;
  float previousDist, stepSize, dummyfloat;
  float det_deltax, det_deltay;
  float d2, dl, l2;

  /*
   * Calculate detector distances
   */
  detDistCopy=detDist;
  det_deltax = *(xdis + nrdet) - *xdis;
  det_deltay = *(ydis + nrdet) - *ydis;
  if (vertical) /* vertical projection */
    {
      *detDistCopy++=1.e12; /* sentinel */
      for (distnr=0 ; distnr <= nrdet ; distnr++)
	{
	  *detDistCopy++ = (x0 * *ydis - *xdis++ * y0) / (*ydis - y0);
	  ydis++;
	}
      *detDistCopy=1.e12; /* sentinel */
      imgCopy=originalImgPtr;
    }
  else /* horizontal projection */
    {
      *detDistCopy++=1.e12; /* sentinel */
      for (distnr=0 ; distnr <= (nrdet) ; distnr++)
	{
	  *detDistCopy++ = -(y0 * *xdis - *ydis++ * x0)/(*xdis - x0);
	  xdis++; /* in row direction: i.e. negative y-axis,
			 so after transpose = positive x-axis */
	}
      *detDistCopy=1.e12; /* sentinel */
      imgCopy=transposeImgPtr;
      dummyint=nrrows;
      nrrows=nrcols;
      nrcols=dummyint;
      dummyfloat=y0;
      y0=-x0;
      x0=-dummyfloat;
      dummyfloat=det_deltay;
      det_deltay=-det_deltax;
      det_deltax=-dummyfloat;
    }

  /*
   * projection --> newproj: scale
   */
  newprojCopy=newproj+1;
  detDistCopy=detDist+1;
  d2 = ((det_deltax)*(det_deltax)+(det_deltay)*(det_deltay));
  for (detnr=0 ; detnr <= nrdet-1 ; detnr++)
    {
      if (*projection != 0)
	{
	  previousDist=(*detDistCopy + *(detDistCopy+1))/2.;
          l2 = (y0*y0+(previousDist-x0)*(previousDist-x0));
	  if (mode == 1)
	    {
	      dl = (det_deltax*(previousDist-x0)-det_deltay*y0);
	      *newprojCopy = *projection * sid * sid / l2 / (1.-dl/l2*dl/d2);
	    }
	  else
	    {
	      *newprojCopy = *projection * sid / l2;
	    }
	}
      else
	{
	  *newprojCopy = 0.;
	}
      projection++;
      newprojCopy++;
      detDistCopy++;
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
          stepSize=-y0/(((nrrows-1.)/2.-rownr)-y0);
	  previousDist=x0+(-nrcols/2.-x0)*stepSize;// correct !
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
	  DD2FanWBackInner(1, nrcols, detDistCopy, previousDist,
			   newprojCopy, imgCopy, stepSize);
	  imgCopy += nrcols;
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
          stepSize=-y0/(((nrrows-1.)/2.-rownr)-y0);
	  previousDist=x0+(-nrcols/2.-x0)*stepSize;
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
	  DD2FanWBackInner(-1, nrcols, detDistCopy, previousDist,
			   newprojCopy, imgCopy, stepSize);
	  imgCopy += nrcols;

	}
    }
}

//-----------------------------------------------------------------------------

/*
 * Distance-driven 2D Fan-beam WeightedBackprojector
 */
#ifdef WIN32
__declspec(dllexport)
#endif
void DD2FanWBack(int nrdet,
		  float x0,
		  float y0,
		  float *xds,
		  float *yds,
		  float xCor,
		  float yCor,
		  float *viewangles,
		  int nrviews,
		  float *sinogram,
		  int nrcols,
		  int nrrows,
		  float *originalImgPtr,
		  float sid,
		  int mode)
{
  float *xdsCopy, *ydsCopy, *newproj;
  float *xdi, *ydi, *xdiCopy, *ydiCopy;
  float *xdiRot, *ydiRot, *xdiRotCopy, *ydiRotCopy;
  float angle, sinAngle, cosAngle, x0Rot, y0Rot;
  float *sinogramCopy, *distances, *viewanglesCopy;
  float *originalImgPtrCopy, *transposeImgPtr, *transposeImgPtrCopy;
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
   * Create empty transpose image
   */
  transposeImgPtr=(float*)calloc(nrcols*nrrows,sizeof(float));

  /*
   * Detector interface coordinates xdi and ydi
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
   * Loop for all views
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
      x0Rot = (x0-xCor)*cosAngle - (y0-yCor)*sinAngle + xCor;
      y0Rot = (y0-yCor)*cosAngle + (x0-xCor)*sinAngle + yCor;
      vertical = (fabs(y0Rot-yCor) >= fabs(x0Rot-xCor));
      
      /*
       * View WeightedBackprojection
       */
      DD2FanWBackView(nrdet, vertical, x0Rot, y0Rot,
		      xdiRot, ydiRot,
		      distances, sinogramCopy, newproj, nrcols, nrrows,
		      originalImgPtr, transposeImgPtr, sid, mode);
      sinogramCopy+=nrdet;
      angle = *viewanglesCopy++;
    }

  /*
   * Add transpose image to original image
   */
  transposeImgPtrCopy=transposeImgPtr;
  originalImgPtrCopy=originalImgPtr;
  for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
    {
      for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
	{
	  *originalImgPtrCopy++ += *transposeImgPtrCopy;
	  transposeImgPtrCopy+=nrrows;
	}
      transposeImgPtrCopy += (1-nrrows*nrcols);
    }

  /*
   * Clean up memory
   */
  free(xdiRot);
  free(ydiRot);
  free(xdi);
  free(ydi);
  free(transposeImgPtr);
  free(distances);
  free(newproj);
}


