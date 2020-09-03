// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

// Given a floating point source array of the form
//  view(row,col,plane)
// An integer map of size
//  map_ndx ~ cols x p
// A floating point
//  map_coef ~ cols x p
// Fills an output array via the sum
//  out(row,col) = sum_p view(row,col,map_ndx(row,p))*map_coef(row,p)
#ifdef WIN32
__declspec(dllexport)
#endif
void viewshift(int rows, int cols, int planes, int pmax,
	       float *view, int *ndx, float *coef, float *output) {
  int i, j, p;
  for (i=0;i<cols;i++)
    for (j=0;j<rows;j++) {
      output[j+i*rows] = 0.0f;
      for (p=0;p<pmax;p++)
	output[j+i*rows] += view[j+i*rows+(ndx[j+p*rows]-1)*rows*cols]*coef[j+p*rows];
    }
}


