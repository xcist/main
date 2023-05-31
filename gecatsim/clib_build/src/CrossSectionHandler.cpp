// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "CrossSectionHandler.hpp"
#include <math.h>
#include <assert.h>

// This is the atomic mass table as extracted from RIVL.
const float AtomicMass[93] = {0,1.00794,4.002602,6.941,9.012182,10.811,12.0107,14.0067,15.9994,18.9984032,20.1797,22.989770,24.3050,26.981538,28.0855,30.973761,32.065,35.453,39.948,39.0983,40.078,44.955910,47.867,50.9415,51.9961,54.938049,55.845,58.933200,58.6934,63.546,65.39,69.723,72.64,74.92160,78.96,79.904,83.80,85.4678,87.62,88.90585,91.224,92.90638,95.94,98.9062546,101.07,102.90550,106.42,107.8682,112.411,114.818,118.710,121.760,127.60,126.90447,131.293,132.90545,137.327,138.9055,140.116,140.90765,144.24,145,150.36,151.964,157.25,158.92534,162.50,164.93032,167.259,168.93421,173.04,174.967,178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.96655,200.59,204.3833,207.2,208.98038,209,210,222,223,226.0254026,227.0277470,232.0381,231.03588,238.02891};


float GetAtomicMass(int Z) {
  return AtomicMass[Z];
}

void BaseDiscreteTable::InitializeTable(int numZ, float MinX, float MaxX, float DeltX) {
  mMaxX = MaxX;
  mMinX = MinX;
  mDeltX = DeltX;
  NumZ = numZ;
  NumX = (int) ((MaxX - MinX)/DeltX+1);
  mValues = MatrixAllocate(NumZ,NumX);
}

Mat BaseDiscreteTable::GetRawTable() {
  return mValues;
}

int BaseDiscreteTable::MapXToBin(float x) {
  return((int) rint((x-mMinX)/mDeltX));
}

float BaseDiscreteTable::GetValue(int Z, float X) {
  int Xbin = MapXToBin(X);
  if (Xbin < 0) return 0;
  assert(Xbin < NumX);
  return mValues[Z][Xbin];
}


void DiscreteTable::InitializeTable(Phantom *phantom, CrossSection& cs, 
				    float MinX, float MaxX, float DeltX) {
  mMaxX = MaxX;
  mMinX = MinX;
  mDeltX = DeltX;
  // Get the list of Z that appear in the phantom
  vector<int> stlZList = phantom->GetZList();
  NumZ = stlZList.size();
  ZList = IVecAllocate(NumZ);
  for (unsigned i=0;i<stlZList.size();i++)
    ZList[i] = stlZList[i];
 // Build the map of Z to index
  mZmap = IVecAllocate(MaxZ);
  for (int i=0;i<NumZ;i++) 
    mZmap[ZList[i]] = i;
  NumX = (int) ((MaxX - MinX)/DeltX+1);
  mValues = MatrixAllocate(NumZ,NumX);
  for (int j=0;j<NumZ;j++) {
    int Z = ZList[j];
    for (int i=0;i<NumX;i++) {
      float X = MinX + i*DeltX;
      mValues[j][i] = cs.GetValue(Z,X);
    }
  }
}

int DiscreteTable::MapXToBin(float x) {
  return((int) rint((x-mMinX)/mDeltX));
}

float DiscreteTable::GetValue(int Z, float X) {
  int Xbin = MapXToBin(X);
  if (Xbin < 0) return Z;
  assert(Xbin < NumX);
  return mValues[mZmap[Z]][Xbin];
}

void CrossSectionHandler::InitializeHandler(Phantom *phantom, CrossSection& cs,
					    float MinEnergy, float MaxEnergy,
					    float DeltEnergy) {
  mMaxEnergy = MaxEnergy;
  mMinEnergy = MinEnergy;
  mDeltEnergy = DeltEnergy;
  mPhantom = phantom;
  // Get the list of Z that appear in the phantom
  vector<int> stlZList = phantom->GetZList();
  NumZ = stlZList.size();
  ZList = IVecAllocate(NumZ);
  for (unsigned i=0;i<stlZList.size();i++)
    ZList[i] = stlZList[i];
  // Build the map of Z to index
  mZmap = IVecAllocate(MaxZ);
  for (int i=0;i<NumZ;i++) 
    mZmap[ZList[i]] = i;
  NumE = (int) ((MaxEnergy - MinEnergy)/DeltEnergy+1);
  mCSDiscreteZvsEBarns = MatrixAllocate(NumZ,NumE);
  mCSDiscreteZvsEMAC = MatrixAllocate(NumZ,NumE);
  for (int j=0;j<NumZ;j++) {
    int Z = ZList[j];
    for (int i=0;i<NumE;i++) {
      float Energy = MinEnergy + i*DeltEnergy;
      mCSDiscreteZvsEBarns[j][i] = cs.GetValue(Z,Energy/1.e6);
      mCSDiscreteZvsEMAC[j][i] = mCSDiscreteZvsEBarns[j][i]*0.6022/AtomicMass[Z]/10.0;
    }
  }
  rhoStore = VecAllocate(NumZ);
}

int CrossSectionHandler::MapEnergyToBin(float Energy) {
  return((int) rint((Energy-mMinEnergy)/mDeltEnergy));
}

float CrossSectionHandler::GetIntegratedCrossSectionBarns(int Z, float Energy) {
  int Ebin = MapEnergyToBin(Energy);
  if (Ebin < 0) return 1e10;
  assert(Ebin < NumE);
  return mCSDiscreteZvsEBarns[mZmap[Z]][Ebin];
}

float CrossSectionHandler::GetIntegratedCrossSectionMAC(int Z, float Energy) {
  int Ebin = MapEnergyToBin(Energy);
  if (Ebin < 0) return 1e10;
  assert(Ebin < NumE);
  return mCSDiscreteZvsEMAC[mZmap[Z]][Ebin];
}

float CrossSectionHandler::GetIntegratedCrossSectionMAC(Photon &p) {
  mPhantom->GetDensities(p,rhoStore);
  int Ebin = MapEnergyToBin(p.GetEnergy());
  assert(Ebin < NumE);
  float CS = 0;
  for (int i=0;i<NumZ;i++) 
    CS += mCSDiscreteZvsEMAC[mZmap[ZList[i]]][Ebin]*rhoStore[i];
  return CS;
}

float CrossSectionHandler::GetIntegratedCrossSectionMuL(float Energy, Vec rhoL) {
  int Ebin = MapEnergyToBin(Energy);
  assert(Ebin < NumE);
  float CS = 0;
  for (int i=0;i<NumZ;i++) 
    CS += mCSDiscreteZvsEMAC[mZmap[ZList[i]]][Ebin]*rhoL[i];
  return CS;
}

vector<float> CrossSectionHandler::GetZProbabilities(Photon &p) {
  vector<float> ret;
  mPhantom->GetDensities(p,rhoStore);
  float totalCS = GetIntegratedCrossSectionMAC(p);
  int Ebin = MapEnergyToBin(p.GetEnergy());
  assert(Ebin < NumE);
  for (int i=0;i<NumZ;i++)
    ret.push_back(mCSDiscreteZvsEMAC[mZmap[ZList[i]]][Ebin]*rhoStore[i]/totalCS);
  return ret;
}

int CrossSectionHandler::GetRandomAtom(Photon &p) {
  return 1;
}


