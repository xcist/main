// GE Proprietary
// last revision: Bruno De Man, May 21, 2014

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()
#include <stdio.h> // printf
#include <string.h>

extern "C"{
#include "DD3.hpp"
}

#define EXPORT


#define MIN(a,b) (a < b) ? (a) : (b)
#define MAX(a,b) (a > b) ? (a) : (b)


void DD3DoseRow(float imgXstart,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *pEnergyMap,
		float *detX, // POINTING TO FIRST REAL DETECTOR BOUNDARY NOW !!!!!!!
		int increment,
		float *detZ, // POINTING TO FIRST REAL DETECTOR BOUNDARY NOW !!!!!!!
		float *scales,
		float z0,
		float *viewEnergyIn,
		float *viewEnergyDep,
		int nrdetrows,
		int nrdetcols,
		float *sils,
		float *pImgFirst,
		float *pImgLast,
		float *energyInFirst,
		float *energyInLast,
		float *energyDepFirst,
		float *energyDepLast,
		float *pEnergyMapFirst,
		float *pEnergyMapLast,
		float *detXFirst,
		float *detXLast,
		float *detZFirst,
		float *detZLast,
		float *silsFirst,
		float *silsLast,
		float *scalesFirst,
		float *scalesLast)

{
  int colnr, planenr, jumpcount, next_is_voxel_boundary;
	float *energyInCopy, *energyDepCopy, *silsCopy;
	float *viewEnergyDep2, *viewEnergyIn2, *sils2;
	float deltaEnergy, inv_det_area;
	float previousX, nextimgX, previousZ, nextimgZ, scale, nextdetZ, dx, imgXlastbutone, imgZlastbutone;
	float inv_imgZstep, inv_imgXstep, detXstart,detXend;
	int detX_ind,detZ_ind;
	int startcol, endcol;
	float detZend;
	int startplane,endplane;
	float *pImgCopy, *pEnergyMapCopy, *pImg2, *pEnergyMap2;
	int detcolnr, detrownr;

	//auxiliary
	inv_imgZstep=1.0f/imgZstep;
	inv_imgXstep=1.0f/imgXstep;

	// copy pointers to temporary pointers
	viewEnergyIn2=viewEnergyIn;
	viewEnergyDep2=viewEnergyDep;
	sils2=sils;
	pImg2=pImg;
	pEnergyMap2=pEnergyMap;

	// img Z variables
	imgZlastbutone=imgZstart+(nrplanes-1)*imgZstep; // LAST BUT ONE img Z boundary

	// img X variables
	colnr=0;
	previousX=imgXstart;
	nextimgX=imgXstart+imgXstep;

	// det X variables
	detX_ind=0;
	detXstart=detX[0];
	detXend=detX[(nrdetcols)*increment];

	//endcol
	endcol=nrcols-1;
	imgXlastbutone=imgXstart+(nrcols-1)*imgXstep; // LAST BUT ONE img X boundary
	if(detXend<imgXlastbutone)
	  {
	    endcol=int((detXend-imgXstart)*inv_imgXstep); //removed +1
	    endcol=MIN(endcol,nrcols-1);
	  }

	//startcol
	startcol=0;
	if(detXstart>previousX)
	  {
	    startcol=int((detXstart-imgXstart)*inv_imgXstep); // removed -1
	    if(startcol>0)
	      {
		pImg2+=startcol*nrplanes;
		pEnergyMap2+=startcol*nrplanes;
		colnr=startcol;
		previousX=imgXstart+imgXstep*startcol;
		nextimgX=previousX+imgXstep;
	      }
	  }
	else // find first det channel that is within the image fov
	  {
	    jumpcount=0; // remember how many detector cells to jump forward
	    while (detX[detX_ind] <= previousX)
	      {
		detX_ind+=increment;
		jumpcount+=increment;
	      } /* now *detX > previousX */
	    viewEnergyIn2 += jumpcount*(nrdetrows+2);
	    viewEnergyDep2 += jumpcount*(nrdetrows+2);
	    sils2 += jumpcount*(nrdetrows+2);
	  }

	//		printf("start-endcol %i %i imgbb %f %f detx %f %f %f \n", startcol, endcol, imgXstart, imgXlastbutone, detX[0], detX[detX_ind], detX[increment*nrdetcols]); getchar();

	/* HIER
	 *  Loop over all img cols
	 */
	while (colnr <= endcol)
	  {
	    scale=scales[detX_ind];

	    // img Z variables
	    planenr=0;
	    previousZ=imgZstart;
	    nextimgZ=imgZstart+imgZstep;

	    // Initialize det Z variables
	    detZ_ind=0;
	    energyInCopy=viewEnergyIn2;
	    energyDepCopy=viewEnergyDep2;
	    silsCopy=sils2;
	    nextdetZ=z0 + scale * detZ[0];

	    detZend = z0 + scale *detZ[nrdetrows];
	    endplane=nrplanes-1;
	    pImgCopy=pImg2;
	    pEnergyMapCopy=pEnergyMap2;

	    if(detZend<imgZlastbutone)
	      {
		endplane=int((detZend-imgZstart)*inv_imgZstep); // removed +1
		endplane=MIN(endplane,nrplanes-1);	
	      }

	    startplane=0;
	    if(nextdetZ>previousZ )
	      {
		startplane=int((nextdetZ-imgZstart)*inv_imgZstep); //removed -1
		if(startplane>0)
		  {
		    pImgCopy+=startplane;
		    pEnergyMapCopy+=startplane;
		    planenr=startplane;
		    previousZ=imgZstart+imgZstep*startplane;
		    nextimgZ=previousZ+imgZstep;	
		  }
	      }
	    else
	      {
		while (nextdetZ <= previousZ)
		  {
		    detZ_ind++;
		    energyInCopy++;
		    energyDepCopy++;
		    silsCopy++;
		    nextdetZ=z0 + scale * detZ[detZ_ind];
		  }
	      }


	    /*
	     *  X OVERLAP KERNEL - PART 1
	     */
	    if (nextimgX <= detX[detX_ind]) /* next X boundary is a voxel boundary */
	      {
		next_is_voxel_boundary=1;
		dx=nextimgX-previousX;
	      }
	    else /* next X boundary is a detector boundary */
	      {
		next_is_voxel_boundary=0;
		dx=detX[detX_ind]-previousX;
	      }

	    /*
	     *  Loop over all img planes (z is first dimension, so this is the inner loop)
	     */
	    while (planenr <= endplane)
	      {
		inv_det_area=1.0/(scale*(detZ[detZ_ind]-detZ[detZ_ind-1])*(detX[detX_ind]-detX[detX_ind-increment]));
		if (nextimgZ <= nextdetZ) /* next Z boundary is a voxel boundary */
		  {
		    if (pEnergyMapCopy < pEnergyMapFirst || pEnergyMapCopy > pEnergyMapLast) {printf("DoseRow1 : pEnergyMap out of bounds\n");}
		    if (energyDepCopy < energyDepFirst || energyDepCopy > energyDepLast) {printf("DoseRow1 : energyDepCopy out of bounds\n");}
		    if (energyInCopy < energyInFirst || energyInCopy > energyInLast) {printf("DoseRow1 : energyInCopy out of bounds\n");}
		    if (silsCopy < silsFirst || silsCopy > silsLast) {printf("DoseRow1 : silsCopy out of bounds\n");}
		    if (pImgCopy < pImgFirst || pImgCopy > pImgLast) {printf("DoseRow1 : pImgCopy out of bounds\n");}

		    deltaEnergy = dx*(nextimgZ-previousZ)*inv_det_area*(*energyInCopy)*(1.0-exp(-(*silsCopy)*(*pImgCopy)));
		    //		    		    		    printf("VXLB col %i, plane %i, dE %f, dx %f, dZ %f, inv %f, EE %f, sils %f, pImg %f\n", colnr, planenr, deltaEnergy, dx, (nextimgZ-previousZ),inv_det_area, *energyInCopy, *silsCopy, *pImgCopy);
		    pImgCopy++;
		    *pEnergyMapCopy++ += deltaEnergy;
		    *energyDepCopy += deltaEnergy;
		    planenr++;
		    previousZ=nextimgZ;
		    nextimgZ+=imgZstep;
		  }
		else /* next Z boundary is a detector boundary */
		  {
		    if (pEnergyMapCopy < pEnergyMapFirst || pEnergyMapCopy > pEnergyMapLast) {printf("DoseRow2 : pEnergyMap out of bounds\n");}
		    if (energyDepCopy < energyDepFirst || energyDepCopy > energyDepLast) {printf("DoseRow2 : energyDepCopy out of bounds\n");}
		    if (energyInCopy < energyInFirst || energyInCopy > energyInLast) {printf("DoseRow2 : energyInCopy out of bounds\n");}
		    if (silsCopy < silsFirst || silsCopy > silsLast) {printf("DoseRow2 : silsCopy out of bounds\n");}
		    if (pImgCopy < pImgFirst || pImgCopy > pImgLast) {printf("DoseRow2 : pImgCopy out of bounds\n");}

		    deltaEnergy = dx*(nextdetZ-previousZ)*inv_det_area*(*energyInCopy)*(1.0-exp(-(*silsCopy)*(*pImgCopy)));
		    //		    		    		    printf("DETB col %i, plane %i, dE %f, dx %f, dZ %f, inv %f, EE %f, sils %f, pImg %f\n", colnr, planenr, deltaEnergy, dx, (nextdetZ-previousZ),inv_det_area, *energyInCopy, *silsCopy, *pImgCopy);
		    *pEnergyMapCopy += deltaEnergy;
		    *energyDepCopy += deltaEnergy;
		    detZ_ind++;
		    energyInCopy++;
		    energyDepCopy++;
		    silsCopy++;
		    previousZ=nextdetZ;
		    nextdetZ=z0 + scale * detZ[detZ_ind];
		  }
	      }

	    /*
	     *  X OVERLAP KERNEL - PART 2
	     */
	    if (next_is_voxel_boundary)
	      {
		previousX=nextimgX;
		nextimgX+=imgXstep;
		colnr++;
		pImg2+=nrplanes;
		pEnergyMap2+=nrplanes;
	      }
	    else /* next X boundary is a detector boundary */
	      {
		previousX=detX[detX_ind];
		detX_ind += increment;
		viewEnergyIn2 += increment * (nrdetrows+2);
		viewEnergyDep2 += increment * (nrdetrows+2);
		sils2 += increment * (nrdetrows+2);
	      }

	  }

	/*
	 * Subtract energyDep from energyIn
	 */
	energyInCopy=viewEnergyIn+(nrdetrows+2)*increment+1; // skip hor. and vert. margin
	energyDepCopy=viewEnergyDep+(nrdetrows+2)*increment+1; // skip hor. and vert. margin
	for (detcolnr=0 ; detcolnr <= (nrdetcols-1) ; detcolnr++)
	  {
	    for (detrownr=0 ; detrownr <= (nrdetrows-1) ; detrownr++)
	      {
		    if (energyInCopy < energyInFirst || energyInCopy > energyInLast) {printf("DoseRow3 : energyInCopy out of bounds\n");}
		    if (energyDepCopy < energyDepFirst || energyDepCopy > energyDepLast) {printf("DoseRow3 : energyDepCopy out of bounds\n");}
		*energyInCopy++ -= *energyDepCopy++;
	      }
	    energyInCopy+=2;
	    energyDepCopy+=2;
	    if (increment == -1)
	      {
		energyInCopy-=2*(nrdetrows+2);
		energyDepCopy-=2*(nrdetrows+2);
	      }
	  }

}





//-----------------------------------------------------------------------------

/*
 * DD3 view projector
 */
void DD3DoseView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 int firstrowfirst,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ, //now size nrzdist + 2
		 float *scales, // empty array (nrdetcols + 2)
		 float *viewEnergyIn, // incoming energy : view with 1 pixel margin all around
		 float *viewEnergyDep, // depositied energy : empty view with 1 pixel margin all around
		 float *sils, // slab-intersection-lengths - empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,
		 int nrplanes,
		 float *pOrig, // size nrplanes x nrcols x nrrows
		 float *pTrans, // size nrplanes x nrcols x nrrows
		 float *pDose, // size nrplanes x nrcols x nrrows
		 float *pDoseTrans, // size nrplanes x nrcols x nrrows
		 float vox_xy_size,
		 float vox_z_size)

{
  float *detXcopy, *pImg, *pEnergyMap, *energyInCopy, *energyDepCopy, *detZcopy, *scalesCopy, *silsCopy;
  int xdistnr, nrxdist, nrzdist;
  int dummyint, increment, rownr, detcolnr, detrownr, firstrow, lastrow, rowinc;
  float dummyfloat, imgXstart, imgXstep, imgZstart, imgZstep;
  float deltaX, deltaZ, mag_fac;

  nrxdist=nrdetcols+1;
  nrzdist=nrdetrows+1;

  /*
   * BOUNDS
   */
  float *pImgFirst, *pImgLast, *energyInFirst, *energyInLast, *energyDepFirst, *energyDepLast,
    *pEnergyMapFirst, *pEnergyMapLast, *detXFirst, *detXLast, *detZFirst, *detZLast,
    *silsFirst, *silsLast, *scalesFirst, *scalesLast;
  energyInFirst=viewEnergyIn;
  energyInLast=viewEnergyIn+(nrdetcols+2)*(nrdetrows+2)-1;
  energyDepFirst=viewEnergyDep;
  energyDepLast=viewEnergyDep+(nrdetcols+2)*(nrdetrows+2)-1;
  detXFirst=detX;
  detXLast=detX+(nrxdist+2)-1;
  detZFirst=detZ;
  detZLast=detZ+(nrzdist+2)-1;
  silsFirst=sils;
  silsLast=sils+(nrdetcols+2)*(nrdetrows+2)-1;
  scalesFirst=scales;
  scalesLast=scales+(nrdetcols+2)-1;

  /*
   * Calculate detX, scales (for detZ) and if necessary transpose the geometry
   */
  detXcopy=detX+1; // size = nrxdist + 2
  scalesCopy=scales+1; // size = nrdetcols + 2
  if (vertical) /* vertical projection */
    {
      for (xdistnr=0 ; xdistnr <= (nrxdist-1) ; xdistnr++) // FOR ALL REAL DET BOUNDARIES
	{
	  if (scalesCopy < scales || scalesCopy > scalesLast){printf("scalesCopy out of bound\n");}
	  if (detXcopy < detX || detXcopy > detXLast){printf("detXCopy out of bound\n");}
	  *detXcopy++ = (x0 * *ydi - *xdi * y0) / (*ydi - y0);/* x-intercept */
	  *scalesCopy++ = y0 / (y0-*ydi); /*   MAG = src-to-x-axis   /   src-to-detectorcell  used to scale det Z coordinates */
	  ydi++;        
	  xdi++;
	}
      pImg=pTrans;
      pEnergyMap=pDoseTrans;

      // BOUNDS
      pImgFirst=pTrans;
      pImgLast=pTrans+(nrcols*nrrows*nrplanes)-1;
      pEnergyMapFirst=pDoseTrans;
      pEnergyMapLast=pDoseTrans+(nrcols*nrrows*nrplanes)-1;

    }
  else /* horizontal projection */
    {
      for (xdistnr=0 ; xdistnr <= (nrxdist-1) ; xdistnr++) // FOR ALL REAL DET BOUNDARIES
	{
	  if (scalesCopy < scales || scalesCopy > scalesLast){printf("scalesCopy out of bound\n");}
	  if (detXcopy < detX || detXcopy > detXLast){printf("detXCopy out of bound\n");}
	  *detXcopy++ = -(y0 * *xdi - *ydi * x0)/(*xdi - x0);
	  /* y-intercept in row direction: i.e. negative y-axis,
	     so after transpose = positive x-axis !!! */
	  *scalesCopy++ = x0 / (x0-*xdi);
	  xdi++; 
	  ydi++;
	}
      pImg=pOrig;
      pEnergyMap=pDose;

      // BOUNDS
      pImgFirst=pOrig;
      pImgLast=pOrig+(nrcols*nrrows*nrplanes)-1;
      pEnergyMapFirst=pDose;
      pEnergyMapLast=pDose+(nrcols*nrrows*nrplanes)-1;

      dummyint=nrrows;
      nrrows=nrcols;
      nrcols=dummyint;
      dummyfloat=y0;
      y0=-x0;
      x0=-dummyfloat;
    }

  /*
   * Scales: shift by 0.5 cells
   */
  *scales=*(scales+1);
  for (scalesCopy=scales+1 ; scalesCopy <= (scales+nrdetcols) ; scalesCopy++)
    {
      *scalesCopy = (*scalesCopy+*(scalesCopy+1))/2.; /* recompute corresponding to cell centers instead of boundaries */
    }

  /*
   * Computer slab intersection lengths (sils)
   */
  silsCopy=sils+(nrdetrows+2)+1; // skip col and row margin
  detXcopy=detX+1; // point to first REAL det col boundary
  scalesCopy=scales+1; // skip col margin
  for (detcolnr=0 ; detcolnr <= (nrdetcols-1) ; detcolnr++)
    {
      deltaX=(*detXcopy + *(detXcopy+1))*0.5-x0;
      detZcopy=detZ+1; // point to first REAL det row boundary
      for (detrownr=0 ; detrownr <= (nrdetrows-1) ; detrownr++)
	{
	  if (silsCopy < silsFirst || silsCopy > silsLast){printf("silsCopy out of bound\n");}
	  if (scalesCopy < scales || scalesCopy > scalesLast){printf("scalesCopy out of bound\n");}
	  if (detXcopy < detXFirst || (detXcopy+1) > detXLast){printf("detXcopy out of bound\n");}
	  if (detZcopy < detZFirst || (detZcopy+1) > detZLast){printf("detZcopy out of bound\n");}
	  deltaZ=(*detZcopy + *(detZcopy+1))*0.5 * (*scalesCopy); //det Z coordinates need to be scaled by col-dependent scale
	  *silsCopy=sqrt(y0*y0+deltaX*deltaX+deltaZ*deltaZ)/fabs(y0)*vox_xy_size/10.0; // intersection length in mm
	  silsCopy++; detZcopy++;
	}
      detXcopy++;
      scalesCopy++;
      silsCopy+=2;
    }

  /*
   * Initialize detector X variables
   */  
  if (*(detX+2) > *(detX+1))
    {
      increment = 1;
      detXcopy=detX+1; // pointing to first REAL detector X boundary
      energyInCopy=viewEnergyIn; // pointing to first MARGIN detector col and row
      energyDepCopy=viewEnergyDep; // pointing to first MARGIN detector col and row
      silsCopy=sils; // pointing to first MARGIN detector col and row
      scalesCopy=scales; // pointing to first MARGIN detector col
      *detX=-1.e12; /* start sentinel */
      *(detX+nrxdist+1)=1.e12; /* end sentinel */
    }
  else
    {
      increment = -1;
      detXcopy=detX+(nrxdist+2)-2; // pointing to last REAL detector x boundary
      energyInCopy=viewEnergyIn+(nrdetcols+1)*(nrdetrows+2); // pointing to last MARGIN detector col, first row
      energyDepCopy=viewEnergyDep+(nrdetcols+1)*(nrdetrows+2); // pointing to last MARGIN detector col, first row
      silsCopy=sils+(nrdetcols+1)*(nrdetrows+2); // pointing to last MARGIN detector col, first row
      scalesCopy=scales+(nrdetcols+2)-1; // pointing to last MARGIN detector col
      *detX=1.e12; /* start sentinel */
      *(detX+nrxdist+1)=-1.e12; /* end sentinel */
    }

  /*
   * For all rows
   */
  if (firstrowfirst)
    {
      firstrow=0;
      lastrow=nrrows-1;
      rowinc=1;
    }
  else
    {
      firstrow=nrrows-1;
      lastrow=0;
      rowinc=-1;
      pImg += (nrrows-1)*nrcols*nrplanes;
      pEnergyMap += (nrrows-1)*nrcols*nrplanes;
    }
  rownr=firstrow;
  for (int ctr=0 ; ctr <= (nrrows-1) ; ctr++)
    {
      rownr+=rowinc;
      /*
       * Reset viewEnergyDep
       */
      memset(viewEnergyDep,0,(nrdetcols+2)*(nrdetrows+2)*sizeof(float));

      /*
       * Initialize image X parameters
       */
      mag_fac=y0/(y0-((nrrows-1)*0.5-rownr)*vox_xy_size); //!! the row direction is opposite to y direction
      imgXstep=mag_fac*vox_xy_size;
      imgXstart=x0-(nrcols*.5*vox_xy_size+x0)*mag_fac;

      /*
       * Initialize image Z parameters
       */
      imgZstep=mag_fac*vox_z_size;
      imgZstart=z0-(nrplanes*0.5*vox_z_size+z0)*mag_fac;

      /*
       * Call DD3ProjRow
       */
      //            printf("row nr %i \n",rownr);
      DD3DoseRow(imgXstart, imgXstep, nrcols,
		 imgZstart, imgZstep, nrplanes,
		 pImg, pEnergyMap, detXcopy, increment, (detZ+1), scalesCopy,
		 z0, energyInCopy, energyDepCopy, nrdetrows,nrdetcols, silsCopy,
		 pImgFirst, pImgLast, energyInFirst,energyInLast,energyDepFirst,energyDepLast,
		 pEnergyMapFirst,pEnergyMapLast, detXFirst, detXLast, detZFirst, detZLast,
		 silsFirst, silsLast, scalesFirst, scalesLast);

      pImg += rowinc*nrcols*nrplanes;
      pEnergyMap += rowinc*nrcols*nrplanes;
    }

}

//-----------------------------------------------------------------------------

/*
 * DD3 add transpose
 */
void DD3AddTranspose(int nrcols,
		     int nrrows,
		     int nrplanes,
		     float *pOrig,
		     float *pTrans)

{
  int rownr, colnr, planenr;

  for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
    {
      for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
	{
	  for (planenr=0 ; planenr<=(nrplanes-1) ; planenr++)
	    {
	      *pTrans++ += *pOrig++; // "+=" iso "=" in DD3Transpose
	    }
	  pTrans += nrplanes * (nrcols-1);
	}
      pTrans += nrplanes * (1 - nrcols*nrrows);
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

  for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
    {
      for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
	{
	  for (planenr=0 ; planenr<=(nrplanes-1) ; planenr++)
	    {
	      *pTrans++ = *pOrig++;
	    }
	  pTrans += nrplanes * (nrcols-1);
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
 * DD3 dose reconstruction
 */
extern "C"{

  void DD3Dose(float x0,
	       float y0,
	       float z0,
	       int nrdetcols,
	       int nrdetrows,
	       float *xds, // size nrdetcols+1
	       float *yds, // size nrdetcols+1
	       float *zds, // size nrdetrows+1
	       float imgXoffset,
	       float imgYoffset,
	       float imgZoffset,
	       float *viewangles, // size nrviews
	       float *zshifts, // size nrviews
	       int nrviews,
	       float *sinogram, // size nrdetrows x nrdetcols x nrviews
	       int nrcols,
	       int nrrows,
	       int nrplanes,
	       float *pOrig, // size nrplanes x nrrows x nrcols
	       float *pDose, // size nrplanes x nrrows x nrcols
	       float vox_xy_size,
	       float vox_z_size)

  {
    int nrxdist, xdistnr, nrzdist, zdistnr, viewnr, vertical, detcolnr, detrownr, firstrowfirst;
    float *xdi, *ydi, *xdiRot, *ydiRot, *viewEnergyIn, *viewEnergyDep, *energyInCopy, *pTrans, *detX;
    float *xdiCopy, *ydiCopy, *xdiRotCopy, *ydiRotCopy;
    float *detZ;
    float *sinogramCopy, *sinogramCopy2, *viewanglesCopy, *scales, *sils, *pDoseTrans;
    float sinAngle, cosAngle, x0Rot, y0Rot, z0Shift;

    /*
     * Parameters
     */
    nrxdist = nrdetcols+1;
    nrzdist = nrdetrows+1;
	 
    /*
     * Allocate memory
     */
    xdi=(float*)malloc(nrxdist*sizeof(float));
    ydi=(float*)malloc(nrxdist*sizeof(float));
    xdiRot=(float*)malloc(nrxdist*sizeof(float));
    ydiRot=(float*)malloc(nrxdist*sizeof(float));
    detZ=(float*)malloc((nrzdist+2)*sizeof(float)); // provide 2 spaces for sentinels
    viewEnergyIn=(float*)calloc((nrdetcols+2)*(nrdetrows+2),sizeof(float)); // 1 cell margin all around
    viewEnergyDep=(float*)calloc((nrdetcols+2)*(nrdetrows+2),sizeof(float)); // 1 cell margin all around
    sils=(float*)calloc((nrdetcols+2)*(nrdetrows+2),sizeof(float)); // 1 cell margin all around
    pTrans=(float*)malloc(nrcols*nrrows*nrplanes*sizeof(float));
    pDoseTrans=(float*)calloc(nrcols*nrrows*nrplanes,sizeof(float));
    detX=(float*)malloc((nrxdist+2)*sizeof(float)); /* provide 2 spaces	for sentinels */
    scales=(float*)malloc((nrdetcols+2)*sizeof(float)); /* provide 1 pixel margin each side */
	  
    /*
     * Create transpose image
     */
    DD3Transpose(nrcols, nrrows, nrplanes, pOrig, pTrans);  

    /*
     * Calculate detector boundaries
     */
    DD3Boundaries(nrxdist, xds, xdi);
    DD3Boundaries(nrxdist, yds, ydi);
    DD3Boundaries(nrzdist, zds, (detZ+1));
    detZ[0] = -1.e12; /* sentinel at the start of detZ array */
    detZ[nrzdist+1] = 1.e12; /* sentinel at the end of detZ array */

    /*
     * Translate detZ and z0 (detZ represents zdi-z0)
     */
    for (zdistnr=0 ; zdistnr <= (nrzdist-1) ; zdistnr++)
      {
	detZ[zdistnr+1]-=(z0);
      }
    z0 -= imgZoffset;

    /*
     * Loop for all views
     */
    sinogramCopy=sinogram;
    sinogramCopy2=sinogram;
    viewanglesCopy=viewangles;
    for (viewnr = 0 ; viewnr <= (nrviews-1) ; viewnr++)
      {
	/*
	 * Rotate around iso and then translate src and det coordinates so they are relative to image center
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
	firstrowfirst = vertical ? (y0Rot>0) : (x0Rot<0);
	x0Rot -= imgXoffset;
	y0Rot -= imgYoffset;
	z0Shift = z0 + *zshifts++;

	/*
	 * Copy sinogram into viewEnergyIn
	 */
	energyInCopy=viewEnergyIn+nrdetrows+2+1; // skip hor. and vert. margin
	for (detcolnr=0 ; detcolnr <= (nrdetcols-1) ; detcolnr++)
	  {
	    for (detrownr=0 ; detrownr <= (nrdetrows-1) ; detrownr++)
	      {
		*energyInCopy++ = *sinogramCopy++;
	      }
	    energyInCopy+=2;
	  }

	/*
	 * View projection
	 */
	DD3DoseView(x0Rot, y0Rot, z0Shift, nrdetcols, nrdetrows,
		    vertical, firstrowfirst, xdiRot, ydiRot, detX, detZ,
		    scales, viewEnergyIn, viewEnergyDep, sils,
		    nrcols, nrrows, nrplanes,
		    pOrig, pTrans, pDose, pDoseTrans, vox_xy_size,vox_z_size);
	//printf("viewnr %i\n",viewnr);


	/*
	 * Copy viewEnergyIn back into sinogram
	 */
	energyInCopy=viewEnergyIn+nrdetrows+2+1; // skip hor. and vert. margin
	for (detcolnr=0 ; detcolnr <= (nrdetcols-1) ; detcolnr++)
	  {
	    for (detrownr=0 ; detrownr <= (nrdetrows-1) ; detrownr++)
	      {
		*sinogramCopy2++ = *energyInCopy++;
	      }
	    energyInCopy+=2;
	  }

      }

    /*
     * Add Dose and DoseTranspose
     */
    DD3AddTranspose(nrrows, nrcols, nrplanes, pDoseTrans, pDose);  

    /*
     * Clean up memory
     */
    free(xdiRot);
    free(ydiRot);
    free(xdi);
    free(ydi);
    free(pTrans);
    free(pDoseTrans);
    free(detX);
    free(detZ);
    free(viewEnergyIn);
    free(viewEnergyDep);
    free(scales);
    free(sils);

  }

}
