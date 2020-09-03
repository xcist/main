// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Volume.hpp"
#include <stdio.h>
#include <string.h>

void WriteRawVector(std::string fname, Vec t, int len) {
  FILE *fp;
  fp = fopen(fname.c_str(),"wb");
  fwrite(t,sizeof(float),len,fp);
  fclose(fp);
}

void ReadRawVector(std::string fname, Vec t, int len) {
  FILE *fp;
  fp = fopen(fname.c_str(),"rb");
  fread(t,sizeof(float),len,fp);
  fclose(fp);
}

IVec IVecAllocate(int n) {
  int *t = new int[n];
  memset(t,0,sizeof(int)*n);
  return t;
}

Vec VecAllocate(int n) {
  float *t = new float[n];
  memset(t,0,sizeof(float)*n);
  return t;
}

Mat MatrixAllocate(int rows, int cols) {
  float *g = new float[rows*cols];
  float **h = new float*[rows];
  for (int i=0;i<rows;i++)
    h[i] = g + i*cols;
  return h;
}

Mat MatrixAllocateAndZero(int rows, int cols) {
  float *g = new float[rows*cols];
  float **h = new float*[rows];
  for (int i=0;i<rows;i++)
    h[i] = g + i*cols;
  memset(g,0,sizeof(float)*rows*cols);
  return h;
}

Vol VolumeAllocate(int rows,int cols, int slabs) {
  float *g = new float[rows*cols*slabs];
  float **h = new float*[rows*slabs];
  float ***k = new float**[slabs];
  for (int i=0;i<rows*slabs;i++)
    h[i] = g + i*cols;
  for (int i=0;i<slabs;i++)
    k[i] = h + i*rows;
  return k;
}

Volset VolsetAllocate(int rows, int cols, int slabs, int sets) {
  float *g = new float[rows*cols*slabs*sets];
  float **h = new float*[rows*slabs*sets];
  float ***k = new float**[slabs*sets];
  float ****m = new float ***[sets];
  for (int i=0;i<rows*slabs*sets;i++)
    h[i] = g + i*cols;
  for (int i=0;i<slabs*sets;i++)
    k[i] = h + i*rows;
  for (int i=0;i<sets;i++)
    m[i] = k + i*slabs;
  return m;
}

void MatrixFree(Mat t) {
  delete t[0];
  delete t;
}

void VecFree(Vec f) {
  delete f;
}

void VolumeFree(Vol t) {
  delete **t;
  delete *t;
  delete t;
}

void VolsetFree(Volset t) {
  delete ***t;
  delete **t;
  delete *t;
  delete t;
}

void IVecFree(IVec f) {
  delete f;
}


