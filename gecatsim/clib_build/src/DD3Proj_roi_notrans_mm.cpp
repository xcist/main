// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()
#include <stdio.h>
#include <string.h>

#include "DD3_roi_notrans_mm.hpp"

//bool useUInt16;

//Version 1.0 ,  mask, and no trans
void DD3ProjRow_roi_notrans_mm(float imgX,
							float imgXstep,
							int nrcols,
							float imgZstart,
							float imgZstep,
							int nrplanes,
							void *pImg,
							float *detX,
							int increment,
							float *detZ,
							float *scales,
							float z0,
							float *view,
							int nrdetrows,
							int nrdetcols,
							int trans_rowstep,
							byte* xy_mask)

{
	int colnr, planenr, jumpcount;
	float *viewCopy;

	float previousX, imgZ, previousZ, scale, nextdetZ, dx,
		imgXStart=imgX, imgZ_end=imgZstart+(nrplanes)*imgZstep;
	float inv_imgZstep=1.0f/imgZstep, inv_imgXstep=1.0f/imgXstep,
		detX_0=detX[0],detX_end=detX[nrdetcols*increment];

	//index variables
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
	int start_imgX_ind=0,end_imgX_ind=nrcols; //default scan boundary

	//X end boundary check
	if(detX_end<(imgXStart+(nrcols-1)*imgXstep))
	{
		end_imgX_ind=int((detX_end-imgXStart)*inv_imgXstep)+2; //extra 1 to avoid precision error
		end_imgX_ind=MIN_(end_imgX_ind,nrcols);
	}

	//X start boundary check
	if(detX_0>previousX)
	{
		start_imgX_ind=int((detX_0-imgXStart)*inv_imgXstep-1); //extra 1 to avoid precision error
		if(start_imgX_ind>0)
		{
			//pImg += start_imgX_ind*nrplanes*trans_rowstep;
		    int tmpstep = start_imgX_ind*nrplanes*trans_rowstep;
            if (useUInt16) {
                pImg = static_cast<unsigned short*>(pImg) + tmpstep;
            } else {
                pImg = static_cast<float*>(pImg) + tmpstep;
            }
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
	}//X boundary check


	/*
	*  Loop over img cols [start_x, end_x)
	*/
	while (colnr < end_imgX_ind)
	{
		scale=scales[detX_ind];
		viewCopy=view;

		/*
		*  Update X variables
		*/
		int mask_ind=colnr; //current mask position,it won't be affected by following change
		if (imgX <= detX[detX_ind]) /* next X boundary is a pixel boundary */
		{
			dx=imgX-previousX; //
			previousX=imgX;
			imgX+=imgXstep;			
			colnr++;
			jumpcount=nrplanes*trans_rowstep; //jumpcount=0;in Bruno's code
		}
		else /* next X boundary is a detector boundary */
		{
			dx=detX[detX_ind]-previousX;
			previousX=detX[detX_ind];
			detX_ind += increment;
			view += increment * (nrdetrows+2);
			jumpcount=0; //jumpcount=-nrplanes;in Bruno's code  
		}

		///mask process, roi
		if(xy_mask[mask_ind])
		{
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
			nextdetZ=z0 + scale * detZ[0];

			float enddetZ = z0 + scale *detZ[nrdetrows];
			int start_imgZ_ind=0,end_imgZ_ind=nrplanes;
			void *pImgcopy;
            pImgcopy = pImg;

			//Z end boundary check
			if(enddetZ<imgZ_end)
			{
				end_imgZ_ind=int((z0+scale*detZ[nrdetrows]-imgZstart)*inv_imgZstep)+2; //extra 1 to avoid precision error
				end_imgZ_ind=MIN_(end_imgZ_ind,nrplanes);	
			}

			//Z start boundary check
			if(nextdetZ>previousZ )
			{
				start_imgZ_ind=int((nextdetZ-imgZstart)*inv_imgZstep)-1; //extra 1 to avoid precision error
				if(start_imgZ_ind>0)
				{
					//pImgcopy += start_imgZ_ind;
                    if (useUInt16) {
                        pImgcopy = static_cast<unsigned short*>(pImgcopy) + start_imgZ_ind;
                    } else {
                        pImgcopy = static_cast<float*>(pImgcopy) + start_imgZ_ind;
                    }
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
			}//z boundary check

			/*
			*  Loop over img planes [start_Z, end_Z), z is first dimension
			*/
			while (planenr < (end_imgZ_ind))
			{
				if (imgZ <= nextdetZ) /* next Z boundary is a pixel boundary */
				{
					//*viewCopy += dx*(imgZ-previousZ)*(*pImgcopy);
                    if (useUInt16) {
					    *viewCopy += dx*(imgZ-previousZ)*(*static_cast<unsigned short*>(pImgcopy));
                        pImgcopy = static_cast<unsigned short*>(pImgcopy) + 1;
                    } else {
					    *viewCopy += dx*(imgZ-previousZ)*(*static_cast<float*>(pImgcopy));
                        pImgcopy = static_cast<float*>(pImgcopy) + 1;
                    }
					planenr++;
					previousZ=imgZ;
					imgZ+=imgZstep;
				}
				else /* next Z boundary is a detector boundary */
				{
					//*viewCopy += dx*(nextdetZ-previousZ)*(*pImgcopy);
                    if (useUInt16) {
					    *viewCopy += dx*(nextdetZ-previousZ)*(*static_cast<unsigned short*>(pImgcopy));
                    } else {
					    *viewCopy += dx*(nextdetZ-previousZ)*(*static_cast<float*>(pImgcopy));
                    }
					detZ_ind++;
					viewCopy++;
					previousZ=nextdetZ;
					nextdetZ=z0 + scale * detZ[detZ_ind];
				}
			}//while z-planes
		}//mask
		//pImg+=(jumpcount);
        if (useUInt16) {
            pImg = static_cast<unsigned short*>(pImg) + jumpcount;
        } else {
            pImg = static_cast<float*>(pImg) + jumpcount;
        }
	}//while colnr  
}


//-----------------------------------------------------------------------------

/*
* DD3 view projector
*/
void DD3ProjView_roi_notrans_mm(float x0,
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
							 void *pOrig,     //        contain
							 float vox_xy_size,   // voxel size
							 float vox_z_size,
							 byte* xy_mask,        //mask
							 byte* xy_mask_trans)    //transposed mask

{
	float *detXcopy, *viewCopy, *detZcopy, *scalesCopy;
	int xdistnr, nrxdist;
	int dummyint, increment, rownr, detcolnr, detrownr;
	float dummyfloat, imgX, imgXstep, imgZ, imgZstep, invCos;
	float deltaX, deltaZ, detXstep, detZstep;

    void *pImg;

	int num_trans_plane_step=1,num_trans_rowstep=1;//to compensate for the removal of pTrans

	byte *p_mask;//for roi reconstruction

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
		num_trans_plane_step=nrcols*nrplanes; //for transpose
		num_trans_rowstep=1;
		p_mask=xy_mask;
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
		pImg=pOrig;//pImg=pTrans;
		num_trans_plane_step=nrplanes;
		num_trans_rowstep=nrcols;
		p_mask=xy_mask_trans;

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
		*scalesCopy = (*scalesCopy+*(scalesCopy+1))*0.5f;
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
	for (rownr=0 ; rownr <nrrows ; rownr++)
	{
		/*
		* Initialize image X parameters
		*/
		float mag_fac=y0/(y0-((nrrows-1)*0.5-rownr)*vox_xy_size);
		imgXstep=mag_fac*vox_xy_size;
		//imgX=x0-(nrcols*0.5*vox_xy_size+x0)*mag_fac; //!! the row direction is opposite to y direction
                imgX=-x0/(y0-((nrrows-1)*0.5-rownr)*vox_xy_size) * ((nrrows-1)*0.5-rownr)*vox_xy_size-(nrcols*0.5*vox_xy_size)*mag_fac;

		/*
		* Initialize image Z parameters
		*     (this is the only place were dzdx comes into play)
		*/
		imgZstep=mag_fac*vox_z_size;
		imgZ=z0-(nrplanes*0.5*vox_z_size+z0)*mag_fac;


		/*
		* Call DD3ProjRow
		*/
		DD3ProjRow_roi_notrans_mm(imgX, imgXstep, nrcols,
			imgZ, imgZstep, nrplanes,
			pImg, detXcopy, increment, detZ, scalesCopy,
			z0, viewCopy, nrdetrows,nrdetcols,num_trans_rowstep,p_mask+rownr*nrcols);

		//pImg += nrcols*nrplanes; //original code
		//pImg+=num_trans_plane_step;
        if (useUInt16) {
            pImg = static_cast<unsigned short*>(pImg) + num_trans_plane_step;
        } else {
            pImg = static_cast<float*>(pImg) + num_trans_plane_step;
        }
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


//-----------------------------------------------------------------------------

/*
* DD3 projector
*/
extern "C"{
DLLEXPORT
	void DD3Proj_roi_notrans_mm(float x0,
		float y0,
                                    float z0,  //first src loc in absolute coords
		int nrdetcols,
		int nrdetrows,
                                    float *xds,//det column locs (for central row)
		float *yds,
                                    float *zds,//zlocs for each row
		float imgXoffset,
		float imgYoffset,
		float imgZoffset,
                                    float *viewangles,//pos angle is counter clockwise
                                    float *zshifts,//src z locs for each view
		int nrviews,
                                    float *sinogram,//z is fastest changing, then col, then view
		int nrcols,         // image
		int nrrows,         //    does NOT
		int nrplanes,       //        contain a dummy 1 pixel frame
                                    void *pOrig,//z is fastest changing (start -), then x (start -), then y (start +)
		float vox_xy_size, //added fields
		float vox_z_size,//added fields 
		byte* xy_mask) //added fields for xy_mask [mask,mask_trans] (mask out xy locations such as those outside FOV)

	{
		int nrxdist, xdistnr, nrzdist, zdistnr, viewnr, vertical;
		float *xdi, *ydi, *xdiRot, *ydiRot, *view, *detX;
		float *xdiCopy, *ydiCopy, *xdiRotCopy, *ydiRotCopy;
		float *detZ/*, *detZcopy*/;
		float *sinogramCopy, *viewanglesCopy, *scales;
		float sinAngle, cosAngle, x0Rot, y0Rot, z0Shift;

		/*
		* Allocate memory for original and rotated detector boundaries
		*/

                /*
                printf("hello\n\r\n");

                printf("xds[e]: %f\n\r\n",xds[100]);

                printf("x0: %f\r\n",x0);
                printf("y0: %f\r\n",y0);
                printf("z0: %f\r\n\n",z0);

                printf("nrdetcols: %d\r\n",nrdetcols);
                printf("nrdetrows: %d\r\n\n",nrdetrows);

                printf("xds[0]: %f\n\r",xds[0]);
                printf("yds[0]: %f\n\r",yds[0]);
                printf("zds[0]: %f\n\r\n",zds[0]);

                printf("xds[e]: %f\n\r",xds[nrdetcols-1]);
                printf("yds[e]: %f\n\r",yds[nrdetcols-1]);
                printf("zds[e]: %f\n\r\n",zds[nrdetrows-1]);

                printf("hello\n\r\n");
                printf("hello\n\r\n");
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

                //		printf("Kai's code.\n");

		/*
		* Loop for all views
		*/
		sinogramCopy=sinogram;
		viewanglesCopy=viewangles;
		for (viewnr = 0 ; viewnr <nrviews ; viewnr++)
		{
			/*
			* Rotate xy-coordinates
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
			DD3ProjView_roi_notrans_mm(x0Rot, y0Rot, z0Shift, nrdetcols, nrdetrows,
				vertical, xdiRot, ydiRot, detX, detZ,
				scales, sinogramCopy, view,
				nrcols, nrrows, nrplanes,
				pOrig, vox_xy_size,vox_z_size,xy_mask,(xy_mask+nrrows*nrcols));
			sinogramCopy += nrdetcols * nrdetrows;
		}

		/*
		* Clean up memory
		*/
		free(xdiRot);
		free(ydiRot);
		free(xdi);
		free(ydi);
		free(detX);
		free(detZ);
		free(view);
		free(scales);
	}

}







