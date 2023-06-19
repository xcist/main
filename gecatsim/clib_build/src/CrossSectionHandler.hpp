// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __CrossSectionHandler_hpp__
#define __CrossSectionHandler_hpp__

#include "Phantom.hpp"
#include "Volume.hpp"
#include <vector>
#include "Photon.hpp"
#include "CrossSection.hpp"

const int MaxZ = 92;

class BaseDiscreteTable {
  Mat mValues;
  float mMinX;
  float mMaxX;
  float mDeltX;
  int NumZ;
  int NumX;
  int MapXToBin(float x);
public:
  void InitializeTable(int numZ, float MinX, float MaxX, float DeltX);
  Mat GetRawTable();
  float GetValue(int Z, float X);
};

class DiscreteTable {
  Mat mValues;
  float mMinX;
  float mMaxX;
  float mDeltX;
  int NumZ;
  int NumX;
  IVec mZmap;
  IVec ZList;
  int MapXToBin(float x);
public:
  void InitializeTable(Phantom *phantom, CrossSection& cs, 
		       float MinX, float MaxX, float DeltX);
  float GetValue(int Z, float X);
};

class CrossSectionHandler {
  Phantom *mPhantom;
  Mat mCSDiscreteZvsEBarns;
  Mat mCSDiscreteZvsEMAC;
  IVec mZmap;
  IVec ZList;
  int NumZ;
  int NumE;
  int MapEnergyToBin(float Energy);
  Vec rhoStore;
  float mMaxEnergy;
  float mMinEnergy;
  float mDeltEnergy;
public:
  void InitializeHandler(Phantom *phantom, CrossSection& cs, 
			 float MinEnergy, float MaxEnergy,
			 float DeltEnergy);
  int GetRandomAtom(Photon &p);
  vector<float> GetZProbabilities(Photon &p);
  float GetIntegratedCrossSectionBarns(int Z, float Energy);
  float GetIntegratedCrossSectionMAC(int Z, float Energy);
  float GetIntegratedCrossSectionMAC(Photon &p);
  float GetIntegratedCrossSectionMuL(float Energy, Vec rhoL);
};

float GetAtomicMass(int Z);

#endif

