// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __Volume_hpp__
#define __Volume_hpp__

#include <string>

typedef int* IVec;
typedef float* Vec;
typedef float** Mat;
typedef float*** Vol;
typedef float**** Volset;

void WriteRawVector(std::string fname, Vec t, int len);
void ReadRawVector(std::string fname, Vec t, int len);
IVec IVecAllocate(int len);
Vec VecAllocate(int len);
Mat MatrixAllocate(int rows, int cols);
Mat MatrixAllocateAndZero(int rows, int cols);
Vol VolumeAllocate(int rows,int cols, int slabs);
Volset VolsetAllocate(int rows, int cols, int slabs, int sets);
void MatrixFree(Mat);
void VecFree(Vec);
void IVecFree(IVec);
void VolumeFree(Vol);
void VolsetFree(Volset);
#endif

