// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include <math.h> // fabs(), sqrt(), sin(), cos()
#include <stdlib.h> // malloc(), free()
#include <math.h> // fabs()
/*
 * RD3 intersections (= RD2Intersections)
 *   IN: x0, y0, deltax, deltay, nrcols, nrrows
 *   OUT: startRow, stopRow, intersection, rico
 */
int RD3Intersections(float x0,
		     float y0,
		     float deltax,
		     float deltay,
		     int nrcols,
		     int nrrows,
		     int* pStartRow,
		     int* pStopRow,
		     float* pIntersection,
		     float* pRico)

{
  float interceptCol, interceptCol2, rico, invrico, dummyfloat, intersection;
  int leftCrossRow, rightCrossRow, startRow, stopRow;

  /*
   * Rico: dcol / drow
   */
  rico = -deltax/deltay; 
  rico = (rico+(float)nrcols)-(float)nrcols;

  /*
   * Interceptcols: crossing with rows 0 and (nrrows-1)
   */
  interceptCol=((float)nrcols-1.0f)*0.5f+x0
    -(((float)nrrows-1.0f)*0.5f-y0)*rico;
  interceptCol2=((float)nrcols-1.0f)*0.5f+x0
    -(-((float)nrrows-1.0f)*0.5f-y0)*rico;
  if (interceptCol <= 0 && interceptCol2 <= 0 
      || interceptCol >= (nrcols-1) && interceptCol2 >= (nrcols-1))
    {
      return 0;
    }

  /*
   * Define start and stop rows so that linear interpolation between two cols is
   * always possible
   */

  /*
   * Vertical
   */
  if (fabs(rico) <= 1.e-6)
    {
      rico=0.;
      if (interceptCol <= 0.0 || interceptCol >= ((float)nrcols-1.0))
	{
	  return 0;
	}
      startRow=0;
      stopRow=nrrows-1;
      intersection = interceptCol;
    }
  /*
   * Top left to bottom right
   */
  else if (rico > 0)
    {
      invrico = ((float)1.)/rico;
      /*
       * leftCrossRow (use col=interceptCol+row*rico)
       */
      leftCrossRow = (int)((-interceptCol)*invrico);
      dummyfloat = interceptCol + ((float)leftCrossRow)*rico;
      while (dummyfloat <= 0.0)
	{
	  leftCrossRow++;
	  dummyfloat = interceptCol + ((float)leftCrossRow)*rico;
	}
      if (leftCrossRow > (nrrows-1)) 
	{
	  return 0;
	}
      startRow = (leftCrossRow > 0) ? leftCrossRow : 0;
      intersection = interceptCol + ((float)startRow)*rico;
      /*
       * rightCrossRow (use col=intersection+(row-startRow)*rico)
       */
      rightCrossRow = (int)((nrcols-1.0-intersection)*invrico)+startRow;
      dummyfloat = intersection + ((float)rightCrossRow-startRow)*rico;
      while (dummyfloat >= nrcols-1.0)
	{
	  rightCrossRow--;
	  dummyfloat = intersection + ((float)rightCrossRow-startRow)*rico;
	}
      if (rightCrossRow < 0)
	{
          return 0;
	}
      stopRow = (rightCrossRow < (nrrows-1)) ? rightCrossRow : (nrrows-1);
    }

  /*
   * Top right to bottom left
   */
  else
    {
      invrico = ((float)1.)/rico;
      /*
       * rightCrossRow (use col=interceptCol+row*rico)
       */
      rightCrossRow = (int)((nrcols-1-interceptCol)*invrico);
      dummyfloat = interceptCol + ((float)rightCrossRow)*rico;
      while (dummyfloat >= nrcols-1)
	{
	  rightCrossRow++;
	  dummyfloat = interceptCol + ((float)rightCrossRow)*rico;
	}
      if (rightCrossRow > (nrrows-1)) 
	{
	  return 0;
	}
      startRow = (rightCrossRow > 0) ? rightCrossRow : 0;
      intersection = interceptCol + ((float)startRow)*rico;
      /*
       * leftCrossRow
       */
      leftCrossRow = startRow-(int)(intersection*invrico);
      dummyfloat = intersection + ((float)leftCrossRow-startRow)*rico;
      while (dummyfloat <= 0)
	{
	  leftCrossRow--;
	  dummyfloat = intersection + ((float)leftCrossRow-startRow)*rico;
	}
      if (leftCrossRow < 0)
	{
          return 0;
	}
      stopRow = (leftCrossRow < (nrrows-1)) ? leftCrossRow : (nrrows-1);
    }

  *pStartRow = startRow;
  *pStopRow = stopRow;
  *pIntersection = intersection;
  *pRico = rico;

  //  cout << "x0: " << x0 << " y0: " << y0 
  //       << "deltax: " << deltax << " deltay: " << deltay 
  //       << "startrow: " << startRow << " stoprow: " << stopRow
  //       << "intersection: " << intersection << " rico: " << rico << endl;
  
  return 1;

}


/*
 * RD3 backprojector for 1 detector column
 */
void RD3BackDetcol(float x0,
		   float y0,
		   float z0,
		   float deltax,
		   float deltay,
		   int nrdetrows,
		   float *zds,
                   float dzdx,
		   float *projection,
		   float *scaledProjection,
		   int nrcols,            // image
		   int nrrows,            //    DOES
		   int nrplanes,          //       contain
		   float* originalImgPtr, //          a dummy 1 pixel frame
		   float* transposeImgPtr)
{
  float *imgPtr, *imgPtrCopy, *projectionCopy, *zdsCopy, *scaledProjectionCopy;
  float *pixel, *pixelLower, *pixelUpper;
  float dummyfloat, intersection, rico;
  float weight, zinterA, zinterB, zintersection, sourcePlane, sourceRow;
  float deltaxy, deltaz, planeWeight;
  int nrplanepixels, dummyint, inside, startRow, stopRow, rownr, leftCol;
  int lowerPlane, oldLowerPlane, detrownr;

  /*
   * If backprojection line is horizontal: switch coordinates and use transpose
   * image
   */
  if (fabs(deltax) > fabs(deltay))
    {
      imgPtr=transposeImgPtr;
      dummyint=nrrows;
      nrrows=nrcols;
      nrcols=dummyint;
      dummyfloat=x0;
      x0=-y0;
      y0=-dummyfloat;
      dummyfloat=deltax;
      deltax=-deltay;
      deltay=-dummyfloat;
    }
  /*
   * If backprojection line is vertical: use original image
   */
  else
    {
      imgPtr=originalImgPtr;
    }

  /*
   * Calculate intersections
   */
  inside = RD3Intersections(x0, y0, deltax, deltay, nrcols, nrrows,
			    &startRow, &stopRow, &intersection, &rico);
  if (inside == 0)
    {
      return;
    }

  /*
   * Scale all projection sums with directional cosines
   */
  deltaxy=(float)sqrt(deltax*deltax+deltay*deltay);
  zdsCopy=zds;
  projectionCopy=projection;
  scaledProjectionCopy=scaledProjection;
  for (detrownr=0 ; detrownr <= nrdetrows-1 ; detrownr++)
    {
      deltaz = dzdx * (*zdsCopy++ - z0);
      *scaledProjectionCopy++ = *projectionCopy++
	* (float)(sqrt(deltaxy*deltaxy+deltaz*deltaz)/fabs(deltay));
    }

  /*
   * For all image rows, interpolate and assign values
   */
  imgPtrCopy = imgPtr + startRow*nrcols;

  sourcePlane = ((float)(nrplanes-1))/2. + z0; // assume +delta_z = +delta_plane
  sourceRow = ((float)(nrrows-1))/2. - y0; // assume +delta_y = -delta_row
  nrplanepixels=nrcols*nrrows;
  for (rownr=startRow ; rownr <= stopRow ; rownr++)
    {
      leftCol = (int) intersection; /*FLOOR*/
      pixel = imgPtrCopy + leftCol;
      weight = intersection-(float)leftCol;
      //
      scaledProjectionCopy=scaledProjection;
      zdsCopy=zds;
      zinterB=-(rownr-sourceRow)/deltay;
      zinterA=sourcePlane-zinterB*z0;
      //
      zintersection = zinterA + zinterB * (*zdsCopy++);
      lowerPlane=((int)(zintersection+2.0f))-2;//FLOOR
      oldLowerPlane=lowerPlane;
      pixelLower=pixel+lowerPlane*nrplanepixels;
      pixelUpper=pixelLower+nrplanepixels;
      planeWeight=zintersection-lowerPlane;
      /*
       * Bilinear interpolate values and add
       */
      if (lowerPlane >= 0 && (lowerPlane+1) <= (nrplanes-1))
	{
	  *(pixelUpper+1)+=planeWeight*weight* *scaledProjectionCopy;
	  *pixelUpper+=planeWeight*(1.-weight)* *scaledProjectionCopy;
	  *(pixelLower+1)+=(1.-planeWeight)*weight* *scaledProjectionCopy;
	  *pixelLower+=(1.-planeWeight)*(1.-weight)* (*scaledProjectionCopy++);
	}
      else
	{
	  scaledProjectionCopy++;
	}
      for (detrownr=1 ; detrownr <= (nrdetrows-1) ; detrownr++)
	{
	  zintersection = zinterA + zinterB * (*zdsCopy++);
	  lowerPlane=((int)(zintersection+2.0f))-2;//FLOOR
	  planeWeight=zintersection-lowerPlane;
          pixelLower += (lowerPlane-oldLowerPlane)*nrplanepixels;
	  pixelUpper=pixelLower+nrplanepixels;
	  oldLowerPlane=lowerPlane;
	  /*
	   * Bilinear interpolate values and add
	   */
          if (lowerPlane >= 0 && (lowerPlane+1) <= (nrplanes-1))
	    {
	      *(pixelUpper+1)+=planeWeight*weight* *scaledProjectionCopy;
	      *pixelUpper+=planeWeight*(1.-weight)* *scaledProjectionCopy;
	      *(pixelLower+1)+=(1.-planeWeight)*weight* *scaledProjectionCopy;
	      *pixelLower+=(1.-planeWeight)*(1.-weight)* (*scaledProjectionCopy++);
	    }
	  else
	    {
	      scaledProjectionCopy++;
	    }
	}
      imgPtrCopy += nrcols;
      intersection += rico;
    }

}

//-----------------------------------------------------------------------------

/*
 * RD3 view backprojector
 */
void RD3BackView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int startDetrow,
		 int stopDetrow,
		 float *xds,
		 float *yds,
		 float *zds,
		 float dzdx,
		 float *projection,
		 float *scaledProjection,
		 int nrcols,
		 int nrrows,             // image
		 int nrplanes,           //     DOES
		 float* originalImgPtr,  //        contain
		 float* transposeImgPtr) //             a dummy 1 pixel frame
{
  float deltax, deltay;
  int detcolnr;

  /*
   * Call RD3 detcol backprojector for all detcols
   */
  for (detcolnr=0 ; detcolnr <= (nrdetcols-1) ; detcolnr++)
    {
      deltax=*xds++ - x0;
      deltay=*yds++ - y0;
      RD3BackDetcol(x0, y0, z0, deltax, deltay, stopDetrow-startDetrow+1,
		    zds+startDetrow, dzdx, projection+startDetrow,
		    scaledProjection+startDetrow,
		    nrcols, nrrows, nrplanes,
		    originalImgPtr, transposeImgPtr);
      projection+=nrdetrows;
    }
}

//-----------------------------------------------------------------------------

/*
 * RD3 backprojector
 */
void RD3Back(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float xCor,
	     float yCor,
	     float *viewangles,
	     int nrviews,
             float *zshifts,
	     float *sinogram,
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *resultImgPtr)
{
  float *xdsRot, *ydsRot, *zdsRot, *xdsRotCopy, *ydsRotCopy, *zdsRotCopy;
  float x0Rot, y0Rot, z0Rot;
  float *xdsCopy, *ydsCopy, *zdsCopy;
  float sinAngle, cosAngle;
  float *transposeImgPtr, *transposeImgPtrCopy;
  float *originalImgPtr, *originalImgPtrCopy;
  float *sinogramCopy, *viewanglesCopy, *resultImgPtrCopy;
  float *scaledProjection;
  int colnr, rownr, planenr, viewnr, detcolnr, detrownr;
  float imsize, z1, z2, ymin, ymax, minmag, maxmag, zmin, zmax;
  int startDetrow, stopDetrow;

  /*
   * Allocate memory for rotated coordinates
   */
  xdsRot=(float*)malloc(nrdetcols*sizeof(float));
  ydsRot=(float*)malloc(nrdetcols*sizeof(float));
  zdsRot=(float*)malloc(nrdetrows*sizeof(float));
  scaledProjection=(float*)malloc(nrdetrows*sizeof(float));
  /*
   * Create empty transpose image
   */
  originalImgPtr=(float*)calloc((nrcols+2)*(nrrows+2)*(nrplanes+2),
				sizeof(float));
  transposeImgPtr=(float*)calloc((nrcols+2)*(nrrows+2)*(nrplanes+2),
				 sizeof(float));

  /*
   * Try to drop some of the detector (Nov 21, 2002)
   */
  imsize=nrcols;
  if (nrrows > nrcols) {imsize=nrrows;}
  z2=(nrplanes+2)/2.;
  z1=-z2;

  /*
   * Initialize and loop over all views
   */
  sinogramCopy=sinogram;
  viewanglesCopy=viewangles;
  for (viewnr = 0 ; viewnr <= (nrviews-1) ; viewnr++)
    {
      /*
       * Rotate coordinates
       */
      sinAngle=(float)sin(*viewanglesCopy);
      cosAngle=(float)cos(*viewanglesCopy++);
      xdsRotCopy=xdsRot;
      ydsRotCopy=ydsRot;
      zdsRotCopy=zdsRot;
      xdsCopy=xds;
      ydsCopy=yds;
      zdsCopy=zds;
      for (detcolnr=0 ; detcolnr<=(nrdetcols-1) ; detcolnr++)
	{
	  *xdsRotCopy++ = (*xdsCopy-xCor)*cosAngle
	    - (*ydsCopy-yCor)*sinAngle + xCor;
	  *ydsRotCopy++ = (*ydsCopy++ - yCor)*cosAngle
	    + (*xdsCopy++-xCor)*sinAngle + yCor;
	}
      for (detrownr=0 ; detrownr<=(nrdetrows-1) ; detrownr++)
	{
          *zdsRotCopy++ = *zdsCopy++ + *zshifts;
	}
      x0Rot = (x0-xCor)*cosAngle - (y0-yCor)*sinAngle + xCor;
      y0Rot = (y0-yCor)*cosAngle + (x0-xCor)*sinAngle + yCor;
      z0Rot = z0 + *zshifts++;

      /*
       * Try to drop some of the detector (Nov 21, 2002)
       */
      ymax=(imsize/2.)*(fabs(sinAngle)+fabs(cosAngle));
      ymin=-ymax;
      maxmag=(y0-*yds)/(y0-ymax);
      minmag=(y0-*yds)/(y0-ymin);
      if (z1 > z0Rot) { zmin=z0Rot+(z1-z0Rot)*minmag; }
      else { zmin=z0Rot+(z1-z0Rot)*maxmag; }
      if (z2 > z0Rot) { zmax=z0Rot+(z2-z0Rot)*maxmag; }
      else { zmax=z0Rot+(z2-z0Rot)*minmag; }
      startDetrow=0;
      zdsRotCopy=zdsRot;
      while (*zdsRotCopy < zmin && startDetrow < (nrdetrows-2)) {
	  zdsRotCopy++;
	  startDetrow++; }
      stopDetrow=nrdetrows-1;
      zdsRotCopy=zdsRot+nrdetrows-1;
      while (*zdsRotCopy > zmax && stopDetrow > 1) {
	zdsRotCopy--;
	stopDetrow--; }

      /*
       * Project view
       */
      RD3BackView(x0Rot, y0Rot, z0Rot, nrdetcols, nrdetrows,
		  startDetrow, stopDetrow,
		  xdsRot, ydsRot, zdsRot, dzdx, sinogramCopy,
		  scaledProjection,
		  (nrcols+2), (nrrows+2), (nrplanes+2),
		  originalImgPtr, transposeImgPtr);
      sinogramCopy+=nrdetcols*nrdetrows;
    }

  /*
   * Add transpose and original image to result image
   */
  transposeImgPtrCopy=transposeImgPtr+((nrcols+2)*(nrrows+2))+(nrcols+2)+1;
  resultImgPtrCopy=resultImgPtr;
  originalImgPtrCopy=originalImgPtr+((nrcols+2)*(nrrows+2))+(nrcols+2)+1;
  for (planenr=0 ; planenr<=nrplanes-1 ; planenr++)
    {
      for (rownr=0 ; rownr<=(nrrows-1) ; rownr++)
	{
	  for (colnr=0 ; colnr<=(nrcols-1) ; colnr++)
	    {
	      *resultImgPtrCopy++ = *originalImgPtrCopy++
		+ *transposeImgPtrCopy;
	      transposeImgPtrCopy += nrrows+2;
	    }
          originalImgPtrCopy += 2;
	  transposeImgPtrCopy += 1-nrcols*(nrrows+2);
	}
      originalImgPtrCopy += 2*(nrcols+2);
      transposeImgPtrCopy += (nrrows+2)*(nrcols+2) - nrrows;
    }

  /*
   * Clean up memory
   */
  free(xdsRot);
  free(ydsRot);
  free(zdsRot);
  free(transposeImgPtr);
  free(originalImgPtr);
  free(scaledProjection);
}


