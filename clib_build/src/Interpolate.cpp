// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern "C" {

EXPORT void Interpolate(int size,
		 float *x,
		 float *y,
		 int sizenew,
		 float *xnew,
		 float *ynew)

{
  int i, j;
  float thex;

  if (size == 1){
    for (i=0;i<sizenew;i++){
     ynew[i] = y[0];
    }
  }
  else {
    if (x[1]>x[0]){
      for (i=0;i<sizenew;i++){
	thex=xnew[i];
	for (j=0 ; j < size ; j++){
	  if (x[j] > thex) {
	    break;
	  }
	}
	if (j==size) {j--;}
	if (j==0) {j++;}
	ynew[i] = y[j]+(thex - x[j])*(y[j-1]- y[j])/(x[j-1]- x[j]);
      }
    }
    else {
      for (i=0;i<sizenew;i++){
	thex=xnew[i];
	for (j=size-1 ; j >= 0 ; j--){
	  if (x[j] > thex) {
	    break;
	  }
	}
	if (j<0) {j=0;}
	if (j==(size-1)) {j--;}
	ynew[i] = y[j]+(thex- x[j])*(y[j+1]- y[j])/(x[j+1]- x[j]);
      }
    }  
  }
}
}


