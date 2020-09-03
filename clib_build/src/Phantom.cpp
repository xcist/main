// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Phantom.hpp"
#include "Photon.hpp"
#include <math.h>

Phantom::Phantom() {
}

Phantom::~Phantom() {
  if (densityMaps) VolsetFree(densityMaps) ; 
  if (doseMap)     VolumeFree(doseMap)     ; 
}

void Phantom::Initialize(int xycount, int zcount, float voxsize, int numZ, int *zlist) {
  m_xycount = xycount;
  m_zcount = zcount;
  m_voxsize = voxsize;
  NumZ = numZ;
  for (int i=0;i<numZ;i++)
    ZList.push_back(zlist[i]);
  densityMaps = VolsetAllocate(xycount,xycount,zcount,numZ);
  doseMap = VolumeAllocate(xycount,xycount,zcount);
}

void Phantom::Load(int numZ, float *rho) {
  for (int i=0;i<m_xycount*m_xycount*m_zcount;i++) {
    densityMaps[numZ][0][0][i] = rho[i];
  }
}

vector<int> Phantom::GetZList() {
  return ZList;
}

void Phantom::GetDensities(Photon &p, Vec rho) {
  for (int i=0;i<NumZ;i++) 
    rho[i] = densityMaps[i][p.GetP()][p.GetM()][p.GetN()];
}

void Phantom::DepositDose(Photon &p, float Energy) {
  doseMap[p.GetP()][p.GetM()][p.GetN()] += Energy;
}

void Phantom::WriteDoseTable(std::string fname) {
  WriteRawVector(fname+"_dose.dat",doseMap[0][0],m_xycount*m_xycount*m_zcount);
}

// Calculate the integral of the density map for each
// atomic component along the given ray (specified by
// the position and direction vector).
#define TMax(a,b) ((a) > (b) ? (a) : (b))
#define TMin(a,b) ((a) < (b) ? (a) : (b))

#define Swap(a,b) {float t = a; a = b; b = t;}
#define sgn(a) ((a) < 0 ? -1 : 1)
#define ispos(a) ((a) > 0 ? 1 : 0)
#define isneg(a) ((a) < 0 ? 1 : 0)
void Phantom::GetIntegratedDensityLength(Vec3 pos, Vec3 dir, Vec rhoL) {
  for (int i=0;i<NumZ;i++)
    rhoL[i] = 0;
  bool escaped = false;
  float voxsize = GetVoxelSize();
  float voxex_x = GetXExtent();
  float voxex_z = GetZExtent();
  float xycount = GetXYCount();
  float zcount = GetZCount();
  float txmin = ( voxex_x - pos.x)/dir.x;
  float txmax = (-voxex_x - pos.x)/dir.x;
  //  float txmax = (voxex_x+(xycount-1)*voxsize - pos.x)/dir.x;
  if (dir.x < 0)
    Swap(txmin,txmax);
  float tymin = ( voxex_x - pos.y)/dir.y;
  float tymax = (-voxex_x - pos.y)/dir.y;
  //float tymax = (voxex_x+(xycount-1)*voxsize - pos.y)/dir.y;
  if (dir.y < 0)
    Swap(tymin,tymax);
  float tzmin = ( voxex_z - pos.z)/dir.z;
  float tzmax = (-voxex_z - pos.z)/dir.z;
  //float tzmax = (voxex_z+(zcount-1)*voxsize - pos.z)/dir.z;
  if (dir.z < 0)
    Swap(tzmin,tzmax);
  if (dir.x == 0) {
    txmin = -1e10;
    txmax = 1e10;
  }
  if (dir.y == 0) {
    tymin = -1e10;
    tymax = 1e10;
  }
  if (dir.z == 0) {
    tzmin = -1e10;
    tzmax = 1e10;
  }
  float ttotalmin = TMax(0,TMax(TMax(txmin,tymin),tzmin));
  float ttotalmax = TMax(0,TMin(TMin(txmax,tymax),tzmax));
  if (ttotalmax < ttotalmin) {
    escaped = true;
    return;
  }
  pos.x += dir.x*ttotalmin;
  pos.y += dir.y*ttotalmin;
  pos.z += dir.z*ttotalmin;
  int m = (int) floor((pos.x-voxex_x)/voxsize);
  int n = (int) floor((pos.y-voxex_x)/voxsize);
  int p = (int) floor((pos.z-voxex_z)/voxsize);
  if (txmin == ttotalmin && ispos(dir.x) == 1) m = 0; 
  if (tymin == ttotalmin && ispos(dir.y) == 1) n = 0;
  if (tzmin == ttotalmin && ispos(dir.z) == 1) p = 0; 
  if (txmin == ttotalmin && isneg(dir.x) == 1) m = int(xycount) - 1; 
  if (tymin == ttotalmin && isneg(dir.y) == 1) n = int(xycount) - 1;
  if (tzmin == ttotalmin && isneg(dir.z) == 1) p = int( zcount) - 1; 
  int deltam = sgn(dir.x);
  int deltan = sgn(dir.y);
  int deltap = sgn(dir.z);
  float deltatx = fabs(voxsize/dir.x);
  if (dir.x == 0) deltatx = 1e10;
  float deltaty = fabs(voxsize/dir.y);
  if (dir.y == 0) deltaty = 1e10;
  float deltatz = fabs(voxsize/dir.z);  
  if (dir.z == 0) deltatz = 1e10;
  float txnext = (voxex_x+(ispos(dir.x)+m)*voxsize-pos.x)/dir.x;
  if (dir.x == 0) txnext = 1e10;
  float tynext = (voxex_x+(ispos(dir.y)+n)*voxsize-pos.y)/dir.y;
  if (dir.y == 0) tynext = 1e10;
  float tznext = (voxex_z+(ispos(dir.z)+p)*voxsize-pos.z)/dir.z;
  if (dir.z == 0) tznext = 1e10;
  txnext = TMax(0,txnext);
  tynext = TMax(0,tynext);
  tznext = TMax(0,tznext);
  while (!escaped) {
    if ((txnext <= tynext) && (txnext <= tznext)) {
      pos.x += dir.x*txnext;
      pos.y += dir.y*txnext;
      pos.z += dir.z*txnext;
      for (int i=0;i<NumZ;i++)
	rhoL[i] += txnext*densityMaps[i][p][m][n];
      tynext -= txnext;
      tznext -= txnext;
      txnext = deltatx;
      m += deltam;
    } else if ((tynext <= txnext) && (tynext <= tznext)) {
      pos.x += dir.x*tynext;
      pos.y += dir.y*tynext;
      pos.z += dir.z*tynext;
      for (int i=0;i<NumZ;i++)
	rhoL[i] += tynext*densityMaps[i][p][m][n];
      txnext -= tynext;
      tznext -= tynext;
      tynext = deltaty;
      n += deltan;
    } else {
      pos.x += dir.x*tznext;
      pos.y += dir.y*tznext;
      pos.z += dir.z*tznext;
      for (int i=0;i<NumZ;i++)
	rhoL[i] += tznext*densityMaps[i][p][m][n];
      txnext -= tznext;
      tynext -= tznext;
      tznext = deltatz;
      p += deltap;
    }
    if ((m < 0) || (m >= xycount) || (n < 0) || (n >= xycount) || 
	(p < 0) || (p >= zcount))
      escaped = true;
  }
}


