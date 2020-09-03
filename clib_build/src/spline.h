// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

void nr_spline(float *x,float *y, int n, float yp1, float ypn, float *y2);
void nr_splint(float *xa,float *ya,float *y2a,int n, float x, float *y);
void spline_interpolate(float *xold, float *yold, int nold, float *xnew, float *ynew, int nnew);

