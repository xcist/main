// GE Proprietary
// based on Kai's modificationss
// last revision: Bruno De Man, May 8, 2014

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()
#include <stdio.h>
#include <string.h>

extern "C"{
#include "DD3.hpp"
}

#define EXPORT


#define MIN(a,b) (a < b) ? (a) : (b)
#define MAX(a,b) (a > b) ? (a) : (b)


//Version 2.1 , with boundaries of Z , with X
void DD3ProjRow_mm(float imgX,
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
				   int nrdetrows,
				   int nrdetcols)

{
	int colnr, planenr, jumpcount;
	float *viewCopy;

	float previousX, imgZ, previousZ, scale, nextdetZ, dx,imgXStart=imgX, imgZ_end=imgZstart+(nrplanes-1)*imgZstep;
	float inv_imgZstep=1.0f/imgZstep, inv_imgXstep=1.0f/imgXstep,detX_0=detX[0],detX_end=detX[nrdetcols*increment];

	//kai, index variables
	int detX_ind=0,detZ_ind=0;

	/*
	*  Initialize img X variables
	*/
	colnr=0;
	previousX=imgX;
	imgX+=imgXstep;

	/*
	*  Initialize det X variables
	*/
	int start_imgX_ind=0,end_imgX_ind=nrcols;

	if(detX_end<(imgXStart+(nrcols-1)*imgXstep))
	{
		end_imgX_ind=int((detX_end-imgXStart)*inv_imgXstep)+2; //extra 1 to avoid precision error
		end_imgX_ind=MIN(end_imgX_ind,nrcols);
	}

	if(detX_0>previousX)
	{
		start_imgX_ind=int((detX_0-imgXStart)*inv_imgXstep-1); //extra 1 to avoid precision error
		if(start_imgX_ind>0)
		{
			pImg+=start_imgX_ind*nrplanes;
			colnr=start_imgX_ind;
			previousX=imgXStart+imgXstep*start_imgX_ind;
			imgX=previousX+imgXstep;
		}
	}
	else
	{
		jumpcount=0; //mem starting det chanel which is within the image fov
		while (detX[detX_ind] <= previousX)
		{
			detX_ind+=increment;
			jumpcount+=increment;
		} /* now *detX > previousX */
		view += jumpcount*(nrdetrows+2); //as detrow is the first dimension
	}


	/*
	*  Loop over all img cols
	*/
	while (colnr < end_imgX_ind)
	{
		//scale=*scales;

		scale=scales[detX_ind];

		/*
		*  Initialize img Z variables
		*/
		planenr=0;
		previousZ=imgZstart;
		imgZ=imgZstart+imgZstep;

		/*
		*  Initialize det Z variables
		*/
		detZ_ind=0;
		viewCopy=view;
		nextdetZ=z0 + scale * detZ[0];

		float enddetZ = z0 + scale *detZ[nrdetrows];
		int start_imgZ_ind=0,end_imgZ_ind=nrplanes;
		float *pImgcopy;
		pImgcopy=pImg;

		if(enddetZ<imgZ_end)
		{
			end_imgZ_ind=int((z0+scale*detZ[nrdetrows]-imgZstart)*inv_imgZstep)+2; //extra 1 to avoid precision error
			end_imgZ_ind=MIN(end_imgZ_ind,nrplanes);	
		}

		if(nextdetZ>previousZ )
		{
			start_imgZ_ind=int((nextdetZ-imgZstart)*inv_imgZstep)-1; //extra 1 to avoid precision error
			if(start_imgZ_ind>0)
			{
				pImgcopy+=start_imgZ_ind;
				planenr=start_imgZ_ind;
				previousZ=imgZstart+imgZstep*start_imgZ_ind;
				imgZ=previousZ+imgZstep;	
			}
		}
		else
		{
			while (nextdetZ <= previousZ)
			{
				detZ_ind++;//detZcopy++;
				viewCopy++;
				nextdetZ=z0 + scale * detZ[detZ_ind];
			}
		}


		/*
		*  Update X variables
		*/
		if (imgX <= detX[detX_ind]) /* next X boundary is a pixel boundary */
		{
			dx=imgX-previousX; //
			previousX=imgX;
			imgX+=imgXstep;
			colnr++;
			jumpcount=nrplanes; //jumpcount=0;
		}
		else /* next X boundary is a detector boundary */
		{
			dx=detX[detX_ind]-previousX;
			previousX=detX[detX_ind];
			detX_ind += increment;
			view += increment * (nrdetrows+2); // why not viewCopy ???
			jumpcount=0; //jumpcount=-nrplanes;  
		}

		/*
		*  Loop over all img planes, z is first dimension
		*/
		while (planenr < (end_imgZ_ind))
		{
			if (imgZ <= nextdetZ) /* next Z boundary is a pixel boundary */
			{
				*viewCopy += dx*(imgZ-previousZ)*(*pImgcopy);
				pImgcopy++;
				planenr++;
				previousZ=imgZ;
				imgZ+=imgZstep;
			}
			else /* next Z boundary is a detector boundary */
			{
				*viewCopy += dx*(nextdetZ-previousZ)*(*pImgcopy);
				detZ_ind++;
				viewCopy++;
				previousZ=nextdetZ;
				nextdetZ=z0 + scale * detZ[detZ_ind];
			}
		}
		pImg+=(jumpcount);
	}  
}





//-----------------------------------------------------------------------------

/*
 * DD3 view projector
 */
void DD3ProjView_mm(float x0,
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
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans,//)    //             a dummy 1 pixel frame
		 float vox_xy_size,   // voxel size
		 float vox_z_size)    //field to represent vox sizes

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
	  for (xdistnr=0 ; xdistnr < (nrxdist) ; xdistnr++)
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
	  for (xdistnr=0 ; xdistnr < (nrxdist) ; xdistnr++)
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
  memset(view,0,(nrdetcols+2)*(nrdetrows+2)*sizeof(float));

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
	  scalesCopy=scales+nrdetcols+1;// pointing to last detector col
  }

  /*
   * For all rows
   */
	for (rownr=0 ; rownr < (nrrows) ; rownr++)
	{
		/*
		* Initialize image X parameters
		*/
		float mag_fac=y0/(y0-((nrrows-1)*0.5-rownr)*vox_xy_size);
		imgXstep=mag_fac*vox_xy_size;
		imgX=x0-(nrcols*.5*vox_xy_size+x0)*mag_fac; //!! the row direction is opposite to y direction

		/*
		* Initialize image Z parameters
		*     (this is the only place were dzdx comes into play)
		*/
		imgZstep=mag_fac*vox_z_size;
		imgZ=z0-(nrplanes*0.5*vox_z_size+z0)*mag_fac;


		/*
		* Call DD3ProjRow
		*/
		DD3ProjRow_mm(imgX, imgXstep, nrcols,
			imgZ, imgZstep, nrplanes,
			pImg, detXcopy, increment, detZ, scalesCopy,
			z0, viewCopy, nrdetrows,nrdetcols);

		pImg += nrcols*nrplanes; //original code
	}

  /*
   * Scale projections
   */
  viewCopy=view+nrdetrows+2+1; // skip hor. and vert. margin
  detXcopy=detX+1;
  scalesCopy=scales+1;
  for (detcolnr=0 ; detcolnr <nrdetcols ; detcolnr++)
  {
	  deltaX=(*detXcopy + *(detXcopy+1))*0.5-x0;
	  detXstep=fabs(*(detXcopy+1)- *detXcopy);
	  detZcopy=detZ;
	  for (detrownr=0 ; detrownr <nrdetrows ; detrownr++)
	  {
		  if (*viewCopy != 0)
		  {
			  deltaZ=(*detZcopy + *(detZcopy+1))*0.5 * *scalesCopy;
			  detZstep=fabs(*(detZcopy+1)- *detZcopy)* *scalesCopy;
			  invCos=sqrt(y0*y0+deltaX*deltaX+deltaZ*deltaZ)/fabs(y0)*vox_xy_size; //in mm
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
extern "C"{

	void DD3Proj_mm(float x0,
			float y0,
			float z0,
			int nrdetcols,
			int nrdetrows,
			float *xds,
			float *yds,
			float *zds,
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
			float *pOrig,
			float vox_xy_size, //added fields
			float vox_z_size) //added fields 

	{
	int nrxdist, xdistnr, nrzdist, zdistnr, viewnr, vertical;
	float *xdi, *ydi, *xdiRot, *ydiRot, *view, *pTrans, *detX;
	float *xdiCopy, *ydiCopy, *xdiRotCopy, *ydiRotCopy;
	float *detZ;
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
	detZ=(float*)malloc((nrzdist+1)*sizeof(float));
	detZ[nrzdist] = 1.e12; /* sentinel in detZ array */
	  
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
		detZ[zdistnr]-=(z0);
	}
	z0 -= imgZoffset;

	/*
	* Prepare empty arrays to contain detX and scales
	*/
	detX=(float*)malloc((nrxdist+2)*sizeof(float)); /* provide 2 spaces	for sentinels */
	scales=(float*)malloc((nrdetcols+2)*sizeof(float)); /* provide 1 pixel margin each side */

	printf("2014 Bruno-Kai-Bruno.\n");

	/*
	* Loop for all views
	*/
	sinogramCopy=sinogram;
	viewanglesCopy=viewangles;
	for (viewnr = 0 ; viewnr <nrviews ; viewnr++)
	{
		/*
		* Rotate and translate xy-coordinates
		*/
		sinAngle=(float)sin(*viewanglesCopy); //! sin is more accurate than sinf
		cosAngle=(float)cos(*viewanglesCopy++);//! cos is more accurate than cosf
		xdiRotCopy=xdiRot;
		ydiRotCopy=ydiRot;
		xdiCopy=xdi;
		ydiCopy=ydi;
		for (xdistnr=0 ; xdistnr<nrxdist ; xdistnr++)
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

		z0Shift = z0 + *zshifts++;

		/*
		* View projection
		*/
		DD3ProjView_mm(x0Rot, y0Rot, z0Shift, nrdetcols, nrdetrows,
			vertical, xdiRot, ydiRot, detX, detZ,
			scales, sinogramCopy, view,
			nrcols, nrrows, nrplanes,
			pOrig, pTrans,vox_xy_size,vox_z_size);
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

	}

}
