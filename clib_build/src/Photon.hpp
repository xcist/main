// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __Photon_hpp__
#define __Photon_hpp__

#include "MatVec.h"
#include "Phantom.hpp"

class Photon {
  Vec3 pos;
  Vec3 dir;
  float Energy;
  double Probability;
  int m, n, p;
  float txnext, tynext, tznext;
  float deltatx, deltaty, deltatz;
  int deltam, deltan, deltap;
  float voxsize, voxex_x, voxex_z;
  int xycount, zcount;
  bool escaped;
  bool interacted;
  Phantom* mPhantom;
  float mCosine;
public:
  Photon();
  void UpdateDeltas();
  void AdvanceVoxel();
  inline int GetM() {return m;}
  inline int GetN() {return n;}
  inline int GetP() {return p;}
  void SetPosition(Vec3 ipos);
  void SetDirection(Vec3 idir);
  void SetCosine(float cosine) {mCosine = cosine;};
  float GetCosine() {return mCosine;}
  inline Vec3 GetPosition() {return pos;};
  inline Vec3 GetDirection() {return dir;};
  inline void SetProbability(float t) {Probability = t;};
  inline float GetProbability() {return Probability;};
  inline float GetEnergy() {return Energy;};
  inline void RotateDirectionVector(float theta, float phi) {dir.RotateVector(theta,phi);}
  void Advance(float dist);
  void AdvanceToPhantom(Phantom* phantom);
  float GetStepSize();
  void SetEnergy(float E);
  inline void Escape() {escaped = true;}
  void Print(std::ostream& o);
  bool HasEscaped() {return escaped;}
  bool HasInteracted() {return interacted;}
  void Interact() {interacted = true;}
};

std::ostream& operator<<(std::ostream&o, Photon& v);


#endif

