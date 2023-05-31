// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __Detector_hpp__
#define __Detector_hpp__

#include "MatVec.h"
#include "Volume.hpp"
#include <string>
#include <stdlib.h>
#include "Photon.hpp"
#include <math.h>

class GenericDetector {
protected:
  Vec3 *mDetectorCenters;
  Vec3 *mDetectorNormals;
  Vec mDetectorCellAreas;
  Mat mCoherentScatterSignal;
  Mat mIncoherentScatterSignal;
  Vec mPrimarySignal;
  Vec mMCSignal;
  int mDetectorCellCount;
  int mNumEvents;
  int mEcount;
public:
  void SetEcount(int maxE) {mEcount = maxE;}
  void SetNumEvents(int cnt) {mNumEvents = cnt;}
  int GetCellCount() {return mDetectorCellCount;}
  Mat GetCoherentSignal() {return mCoherentScatterSignal;}
  Mat GetIncoherentSignal() {return mIncoherentScatterSignal;}
  Vec GetPrimarySignal() {return mPrimarySignal;}
  Vec GetMCSignal() {return mMCSignal;}
  Vec3* GetCenters() {return mDetectorCenters;}
  Vec3* GetNormals() {return mDetectorNormals;}
  Vec GetCellAreas() {return mDetectorCellAreas;}
  virtual bool RecordPhotonMC(Photon& p, bool isPrimary = false) = 0;
  virtual float GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir) {return fabs(mDetectorCellAreas[cellNum]*dir.Dot(mDetectorNormals[cellNum]));}
  void WriteData(std::string fname, bool MConly, bool PrimaryValid);
};

class FocallyAlignedXCollimatedDetector : public GenericDetector {
protected:
  float mSDD, mSID, mCellWidth, mCellHeight;
  int mNRows, mNCols;
  int mNRowsDecimated, mNColsDecimated;
  float mCenterRow, mCenterCol;
  float mGridHeight;
  Vec3* mDetectorToCollimatorUnitVec;
  Vec3* mDetectorAlongCollimatorUnitVec;
  bool mDecimated;
  int mDetColDecimationFactor;
  int mDetRowDecimationFactor;
public:
  FocallyAlignedXCollimatedDetector(float SDD, float SID, float cellWidth, 
				    float cellHeight, int nrows, int ncols, 
				    float centerrow, float centercol, 
				    float gridHeight, int coldecimation,
				    int rowdecimation, int maxEinKV);
  virtual ~FocallyAlignedXCollimatedDetector();
  Vec3* GetToCollimatorUnitVecs() {return mDetectorToCollimatorUnitVec;}
  Vec3* GetAlongCollimatorUnitVecs() {return mDetectorAlongCollimatorUnitVec;}
  virtual float GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir);
  bool RecordPhotonMC(Photon& p, bool isPrimary);
  int GetEffectiveRowCount() {return mNRowsDecimated;}
  int GetEffectiveColCount() {return mNColsDecimated;}
  int GetFinalRowCount() {return mNRows;}
  int GetFinalColCount() {return mNCols;}
  int GetRowDecimationFactor() {return mDetRowDecimationFactor;}
  int GetColDecimationFactor() {return mDetColDecimationFactor;}
};

class FocallyAlignedXandZCollimatedDetector : public FocallyAlignedXCollimatedDetector {

  // This type of detector has collimator blades that are focused in x, but not in z
  // The second set of collimator blades are all parallel to the xy-plane.
  // If the cone angle of the system is not too large, this is a good approximation of a collimator that is focussed in both directions

public:
  FocallyAlignedXandZCollimatedDetector(float SDD, float SID, float cellWidth, 
				    float cellHeight, int nrows, int ncols, 
				    float centerrow, float centercol, 
				    float gridHeight, int coldecimation,
				    int rowdecimation, int maxEinKV) : FocallyAlignedXCollimatedDetector(SDD, SID, cellWidth,
													 cellHeight, nrows, ncols,
													 centerrow, centercol,
													 gridHeight, coldecimation,
													 rowdecimation, maxEinKV) {}
  virtual float GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir);
};



// Included flat detector class - Pablo Milioni (GE Healthcare)

class XAlignedZCollimatedDetectorFlat : public GenericDetector {
  
  // FIX ME:
  // This type of detector should have collimator blades that are alligned with the x-axis,
  // The intersection test for the colimator is not ready yet
  // Forced mGridHeight = 0
  
protected:
  float mSDD, mSID, mCellWidth, mCellHeight;
  int mNRows, mNCols;
  int mNRowsDecimated, mNColsDecimated;
  float mCenterRow, mCenterCol;
  float mGridHeight;
  Vec3* mDetectorToCollimatorUnitVec;
  Vec3* mDetectorAlongCollimatorUnitVec;
  bool mDecimated;
  int mDetColDecimationFactor;
  int mDetRowDecimationFactor;
public:
  XAlignedZCollimatedDetectorFlat(float SDD, float SID, float cellWidth, 
				    float cellHeight, int nrows, int ncols, 
				    float centerrow, float centercol, 
				    float gridHeight, int coldecimation,
				    int rowdecimation, int maxEinKV);
  virtual ~XAlignedZCollimatedDetectorFlat();
  Vec3* GetToCollimatorUnitVecs() {return mDetectorToCollimatorUnitVec;}
  Vec3* GetAlongCollimatorUnitVecs() {return mDetectorAlongCollimatorUnitVec;}
  virtual float GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir);
  bool RecordPhotonMC(Photon& p, bool isPrimary);
  int GetEffectiveRowCount() {return mNRowsDecimated;}
  int GetEffectiveColCount() {return mNColsDecimated;}
  int GetFinalRowCount() {return mNRows;}
  int GetFinalColCount() {return mNCols;}
  int GetRowDecimationFactor() {return mDetRowDecimationFactor;}
  int GetColDecimationFactor() {return mDetColDecimationFactor;}
};


#endif

