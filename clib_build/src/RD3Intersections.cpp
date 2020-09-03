// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

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


