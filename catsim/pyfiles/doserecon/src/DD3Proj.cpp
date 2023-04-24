// GE Proprietary
// last revision: Bruno De Man, Jun 21, 2004

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()

extern "C"{
#include "DD3.hpp"
}

#define EXPORT

/*
 * DD3 row projector
 */
void DD3ProjRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		int nrdetrows)

{
  int colnr, planenr, jumpcount;
  float previousX, imgZ, previousZ, scale, nextdetZ, dx;
  float *detZcopy, *viewCopy;

  /*
   *  Initialize img X variables
   */
  colnr=0;
  previousX=imgX;
  imgX+=imgXstep;

  /*
   *  Initialize det X variables
   */
  jumpcount=0;
  while (*detX <= previousX)
    {
      detX+=increment;
      jumpcount+=increment;
    } /* now *detX > previousX */
  scales += jumpcount;
  view += jumpcount*(nrdetrows+2);
      
  /*
   *  Loop over all img cols
   */
  while (colnr <= (nrcols-1))
    {
      scale=*scales;

      /*
       *  Initialize img Z variables
       */
      planenr=0;
      previousZ=imgZstart;
      imgZ=imgZstart+imgZstep;

      /*
       *  Initialize det Z variables
       */
      detZcopy=detZ;
      viewCopy=view;
      nextdetZ=z0 + scale * *detZcopy;
      while (nextdetZ <= previousZ)
	{
	  detZcopy++;
	  viewCopy++;
	  nextdetZ=z0 + scale * *detZcopy;
	}

      /*
       *  Update X variables
       */
      if (imgX <= *detX) /* next X boundary is a pixel boundary */
	{
	  dx=imgX-previousX;
	  previousX=imgX;
	  imgX+=imgXstep;
	  colnr++;
	  jumpcount=0;
	}
      else /* next X boundary is a detector boundary */
	{
	  dx=*detX-previousX;
	  previousX=*detX;
	  detX += increment;
	  view += increment * (nrdetrows+2);
	  scales += increment;
	  jumpcount=-nrplanes;
	}

      /*
       *  Loop over all img planes
       */
      while (planenr <= (nrplanes-1))
	{
	  if (imgZ <= nextdetZ) /* next Z boundary is a pixel boundary */
	    {
	      *viewCopy += dx*(imgZ-previousZ)*(*pImg);
	      pImg++;
	      planenr++;
	      previousZ=imgZ;
	      imgZ+=imgZstep;
	    }
	  else /* next Z boundary is a detector boundary */
	    {
	      *viewCopy += dx*(nextdetZ-previousZ)*(*pImg);
	      detZcopy++;
	      viewCopy++;
	      previousZ=nextdetZ;
	      nextdetZ=z0 + scale * *detZcopy;
	    }
	}
      pImg+=jumpcount;
    }
  
}

//-----------------------------------------------------------------------------

/*
 * DD3 view projector
 */
void DD3ProjView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ,
		 float *scales, // empty array (nrdetcols + 2)
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans)    //             a dummy 1 pixel frame

{
  float *detXcopy, *pImg, *viewCopy, *detZcopy, *scalesCopy;
  int xdistnr, nrxdist;
  int dummyint, increment, rownr, detcolnr, detrownr;
  float dummyfloat, imgX, imgXstep, imgZ, imgZstep, invCos;
  float deltaX, deltaZ, detXstep, detZstep;

  /*
   * Calculate detX, and if necessary transpose the geometry
   */
  detXcopy=detX; // size = nrxdist + 2
  scalesCopy=scales+1;
  nrxdist = nrdetcols+1;
  *detXcopy++=1.e12; /* sentinel */
  if (vertical) /* vertical projection */
    {
      for (xdistnr=0 ; xdistnr <= (nrxdist-1) ; xdistnr++)
	{
	  *detXcopy++ = (x0 * *ydi - *xdi * y0) / (*ydi - y0);/* x-intercept */
	  *scalesCopy++ = y0 / (y0-*ydi);
	  ydi++;        
	  xdi++;
	}
      pImg=pOrig;
    }
  else /* horizontal projection */
    {
      for (xdistnr=0 ; xdistnr <= (nrxdist-1) ; xdistnr++)
	{
	  *detXcopy++ = -(y0 * *xdi - *ydi * x0)/(*xdi - x0);
	  /* y-intercept in row direction: i.e. negative y-axis,
	     so after transpose = positive x-axis */
	  *scalesCopy++ = x0 / (x0-*xdi);
	  xdi++; 
	  ydi++;
	}
      pImg=pTrans;
      dummyint=nrrows;
      nrrows=nrcols;
      nrcols=dummyint;
      dummyfloat=y0;
      y0=-x0;
      x0=-dummyfloat;
    }
  *detXcopy=1.e12; /* sentinel */
  *scales=*(scales+1);
  for (scalesCopy=scales+1 ; scalesCopy<=(scales+nrdetcols) ; scalesCopy++)
    {
      *scalesCopy = (*scalesCopy+*(scalesCopy+1))/2.;
    }

  /*
   * Reset view
   */
  for (viewCopy = view ;
       viewCopy <= view+(nrdetcols+2)*(nrdetrows+2)-1 ;
       viewCopy++)
    {
      *viewCopy = 0.;
    }

  /*
   * Initialize detector X variables
   */  
  if (*(detX+2) > *(detX+1))
    {
      increment = 1;
      detXcopy=detX+1;
      viewCopy=view; // pointing to first detector col
      scalesCopy=scales;
    }
  else
    {
      increment = -1;
      detXcopy=detX+nrdetcols+1;
      viewCopy=view+(nrdetcols+1)*(nrdetrows+2);
      // pointing to last detector col
      scalesCopy=scales+nrdetcols+1;
    }

  /*
   * For all rows
   */
  for (rownr=0 ; rownr <= (nrrows-1) ; rownr++)
    {
      /*
       * Initialize image X parameters
       */
      imgXstep=y0/(y0-((nrrows-1.)/2.-rownr));
      imgX=x0-(nrcols/2.+x0)*imgXstep;// correct !

      /*
       * Initialize image Z parameters
       *     (this is the only place were dzdx comes into play)
       */
      imgZstep=imgXstep;//*dzdx;
      imgZ=z0-(nrplanes/2.+z0)*imgZstep;

      /*
       * Call DD3ProjRow
       */
      DD3ProjRow(imgX, imgXstep, nrcols,
		 imgZ, imgZstep, nrplanes,
		 pImg, detXcopy, increment, detZ, scalesCopy,
		 z0, viewCopy, nrdetrows);
      pImg += nrcols*nrplanes;
    }

  /*
   * Scale projections
   */
  viewCopy=view+nrdetrows+2+1; // skip hor. and vert. margin
  detXcopy=detX+1;
  scalesCopy=scales+1;
  for (detcolnr=0 ; detcolnr <= nrdetcols-1 ; detcolnr++)
    {
      deltaX=(*detXcopy + *(detXcopy+1))/2.-x0;
      detXstep=fabs(*(detXcopy+1)- *detXcopy);
      detZcopy=detZ;
      for (detrownr=0 ; detrownr <= nrdetrows-1 ; detrownr++)
	{
	  if (*viewCopy != 0)
	    {
	      deltaZ=(*detZcopy + *(detZcopy+1))/2.* *scalesCopy;
	      detZstep=fabs(*(detZcopy+1)- *detZcopy)* *scalesCopy;
              invCos=sqrt(y0*y0+deltaX*deltaX+dzdx*dzdx*deltaZ*deltaZ)/fabs(y0);
	      *sinogram++ = invCos / (detXstep*detZstep) * *viewCopy++;
	    } // divide by cos(alpha) and divide by det. stepsizes
	  else
	    {
	      viewCopy++;
	      *sinogram++ = 0.;
	    }
	  detZcopy++;
	}
      detXcopy++;
      scalesCopy++;
      viewCopy+=2;
    }

}

//-----------------------------------------------------------------------------

/*
 * DD3 transpose
 */
void DD3Transpose(int nrcols,
		  int nrrows,
		  int nrplanes,
		  float *pOrig,
		  float *pTrans)
{
  int rownr, colnr, planenr;

  for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
    {
      for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
	{
	  for (planenr=0 ; planenr<=(nrplanes-1) ; planenr++)
	    {
	      *pTrans++ = *pOrig++;
	    }
	  pTrans += nrplanes * (nrrows-1);
	}
      pTrans += nrplanes * (1 - nrcols*nrrows);
    }
}

//-----------------------------------------------------------------------------

/*
 * DD3 boundaries
 */
void DD3Boundaries(int nrBoundaries,
		   float *pCenters,
		   float *pBoundaries)
{
  int i;
  if (nrBoundaries >= 3)
    {
    *pBoundaries++ = 1.5 * *pCenters - 0.5 * *(pCenters+1);
    for (i=1 ; i<=(nrBoundaries-2) ; i++)
      {
	*pBoundaries++ = 0.5 * *pCenters + 0.5 * *(pCenters+1);
	pCenters++;
      }
    *pBoundaries = 1.5 * *pCenters - 0.5 * *(pCenters-1);
    }
  else
    {
      *pBoundaries = *pCenters-0.5;
      *(pBoundaries+1) = *pCenters+0.5;
    }
}

//-----------------------------------------------------------------------------

/*
 * DD3 projector
 */
extern "C" {
EXPORT
void DD3Proj(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,         // image
	     int nrrows,         //    does NOT
	     int nrplanes,       //        contain a dummy 1 pixel frame
	     float *pOrig)

{
  int nrxdist, xdistnr, nrzdist, zdistnr, viewnr, vertical;
  float *xdi, *ydi, *xdiRot, *ydiRot, *view, *pTrans, *detX;
  float *xdiCopy, *ydiCopy, *xdiRotCopy, *ydiRotCopy;
  float *detZ, *detZcopy, *detZshift, *detZshiftCopy;
  float *sinogramCopy, *viewanglesCopy, *scales;
  float sinAngle, cosAngle, x0Rot, y0Rot, z0Shift;
 
  /*
   * Allocate memory for original and rotated detector boundaries
   */
  nrxdist = nrdetcols+1;
  nrzdist = nrdetrows+1;
  xdi=(float*)malloc(nrxdist*sizeof(float));
  ydi=(float*)malloc(nrxdist*sizeof(float));
  xdiRot=(float*)malloc(nrxdist*sizeof(float));
  ydiRot=(float*)malloc(nrxdist*sizeof(float));
  detZ=(float*)malloc((nrzdist)*sizeof(float));
  detZshift=(float*)malloc((nrzdist+2)*sizeof(float)); // 1 space for sentinel
  // 2 sentinels ; one for start-loop, one for main-loop

  /*
   * Allocate and reset memory for 1 view with 1 pixel margin all around
   */
  view=(float*)calloc((nrdetcols+2)*(nrdetrows+2),sizeof(float));

  /*
   * Create transpose image
   */
  pTrans=(float*)malloc(nrcols*nrrows*nrplanes*sizeof(float));
  DD3Transpose(nrcols, nrrows, nrplanes, pOrig, pTrans);  

  /*
   * Calculate detector boundaries
   */
  DD3Boundaries(nrxdist, xds, xdi);
  DD3Boundaries(nrxdist, yds, ydi);
  DD3Boundaries(nrzdist, zds, detZ);

  /*
   * Translate detZ and z0 (detZ represents zdi-z0)
   */
  for (zdistnr=0 ; zdistnr <= (nrzdist-1) ; zdistnr++)
    {
      *(detZ+zdistnr)-=(z0);
    }
  z0 -= imgZoffset;

  /*
   * Prepare empty arrays to contain detX and scales
   */
  detX=(float*)malloc((nrxdist+2)*sizeof(float)); /* provide 2 spaces
						     for sentinels */
  scales=(float*)malloc((nrdetcols+2)*sizeof(float)); /* provide 1 pixel margin
							 each side */

  /*
   * Loop for all views
   */
  sinogramCopy=sinogram;
  viewanglesCopy=viewangles;
  for (viewnr = 0 ; viewnr <= (nrviews-1) ; viewnr++)
    {
      /*
       * Rotate xy-coordinates
       */
      sinAngle=(float)sin(*viewanglesCopy);
      cosAngle=(float)cos(*viewanglesCopy++);
      xdiRotCopy=xdiRot;
      ydiRotCopy=ydiRot;
      xdiCopy=xdi;
      ydiCopy=ydi;
      for (xdistnr=0 ; xdistnr<=(nrxdist-1) ; xdistnr++)
	{
	  *xdiRotCopy++ = *xdiCopy * cosAngle - *ydiCopy * sinAngle
	    - imgXoffset;
	  *ydiRotCopy++ = *ydiCopy++ * cosAngle + *xdiCopy++ * sinAngle
	    - imgYoffset;
	}
      x0Rot = x0 * cosAngle - y0 * sinAngle;
      y0Rot = y0 * cosAngle + x0 * sinAngle;
      vertical = (fabs(y0Rot) >= fabs(x0Rot));
      x0Rot -= imgXoffset;
      y0Rot -= imgYoffset;

      /*
       * Shift z coordinates
       */
      detZshiftCopy=detZshift;
      detZcopy=detZ;
      for (zdistnr=0 ; zdistnr<=(nrzdist-1) ; zdistnr++)
      	{
      	  *detZshiftCopy++ = *detZcopy++; // + *zshifts;
      	}
      *(detZshift+nrzdist) = 1.e12; /* sentinel in shifted array */
      *(detZshift+nrzdist+1) = 1.5e12; /* sentinel in shifted array */
      z0Shift = z0 + *zshifts++;
      
      /*
       * View projection
       */
      DD3ProjView(x0Rot, y0Rot, z0Shift, nrdetcols, nrdetrows,
		  vertical, xdiRot, ydiRot, detX, detZshift,
		  scales, dzdx, sinogramCopy, view,
		  nrcols, nrrows, nrplanes,
		  pOrig, pTrans);
      sinogramCopy += nrdetcols * nrdetrows;
    }

  /*
   * Clean up memory
   */
  free(xdiRot);
  free(ydiRot);
  free(xdi);
  free(ydi);
  free(pTrans);
  free(detX);
  free(detZ);
  free(view);
  free(scales);
  free(detZshift);
}

}

