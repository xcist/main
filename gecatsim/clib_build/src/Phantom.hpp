// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __Phantom_hpp__
#define __Phantom_hpp__

#include "Volume.hpp"
#include <vector>
#include "MatVec.h"

using namespace std;

class Photon;

// Describes the phantom
class Phantom {
  int materialCount;
  Volset densityMaps;
  Vol doseMap;
  int m_xycount;
  int m_zcount;
  float m_voxsize;
  vector<int> ZList;
  int NumZ;
public:
  Phantom();
   ~Phantom(); 
 void Initialize(int xycount, int zcount, float voxsize, int numZ, int *zlist);
  void Load(int numZ, float *rho);
  int GetXYCount() {return m_xycount;}
  int GetZCount() {return m_zcount;}
  float GetVoxelSize() {return m_voxsize;}
  float GetXExtent() {return -(m_xycount*m_voxsize)/2.;}
  float GetZExtent() {return -(m_zcount*m_voxsize)/2.;}
  vector<int> GetZList();
  void GetDensities(Photon &p, Vec rho);
  void DepositDose(Photon &p, float Energy);
  void GetIntegratedDensityLength(Vec3 pos, Vec3 dir, Vec rhoL);
  void WriteDoseTable(std::string fname);
  Vol GetDoseMap() {return doseMap;}
};

#endif

