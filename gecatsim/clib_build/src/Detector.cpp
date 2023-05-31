// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Detector.hpp"
#include <math.h>

void GenericDetector::WriteData(std::string fname, bool MConly, bool PrimaryValid) {
  if (PrimaryValid)
    WriteRawVector(fname+"_primary.dat",mPrimarySignal,mDetectorCellCount);
  else {
    if (!MConly) {
      WriteRawVector(fname+"_coherent.dat",mCoherentScatterSignal[0],mDetectorCellCount*mEcount);
      WriteRawVector(fname+"_incoherent.dat",mIncoherentScatterSignal[0],mDetectorCellCount*mEcount);
    } else {
      WriteRawVector(fname+"_mc.dat",mMCSignal,mDetectorCellCount);
    }
  }
}

FocallyAlignedXCollimatedDetector::FocallyAlignedXCollimatedDetector(float SDD, float SID, 
								     float cellWidth, 
								     float cellHeight,
								     int nrows, int ncols, 
								     float centerrow, 
								     float centercol,
								     float gridHeight,
								     int coldecimation, 
								     int rowdecimation,
								     int maxEinKV) {
  // Count the actual number of rows and columns
  int effrows, effcols;
  effrows = 0;
  for (int i=0;i<nrows;i+=rowdecimation)
    effrows++;
  effcols = 0;
  for (int i=0;i<ncols;i+=coldecimation)
    effcols++;
  mDecimated = ((rowdecimation != 1) || (coldecimation != 1));
  mDetRowDecimationFactor = rowdecimation;
  mDetColDecimationFactor = coldecimation;
  mDetectorCellCount = effrows*effcols;
  mDetectorCenters = new Vec3[mDetectorCellCount];
  mDetectorNormals = new Vec3[mDetectorCellCount];
  mDetectorToCollimatorUnitVec  = new Vec3[mDetectorCellCount];
  mDetectorAlongCollimatorUnitVec = new Vec3[mDetectorCellCount];
  mDetectorCellAreas = VecAllocate(mDetectorCellCount);
  for (int i=0;i<effrows;i++) {
    for (int j=0;j<effcols;j++) {
      float alpha = (j*coldecimation - centercol)*cellWidth/SDD;
      int k=i*effcols+j;
      mDetectorCenters[k].x = sin(alpha)*SDD;
      mDetectorCenters[k].y = SID-cos(alpha)*SDD;
      mDetectorCenters[k].z = (i*rowdecimation - centerrow)*cellHeight;
      mDetectorNormals[k].x = -sin(alpha);
      mDetectorNormals[k].y = cos(alpha);
      mDetectorNormals[k].z = 0;
      mDetectorToCollimatorUnitVec[k].x = cos(alpha);
      mDetectorToCollimatorUnitVec[k].y = sin(alpha);
      mDetectorToCollimatorUnitVec[k].z = 0;
      mDetectorAlongCollimatorUnitVec[k].x = 0;
      mDetectorAlongCollimatorUnitVec[k].y = 0;
      mDetectorAlongCollimatorUnitVec[k].z = 1;
      mDetectorCellAreas[k] = cellWidth*cellHeight;
    }
  }
  mEcount = maxEinKV;
  mCoherentScatterSignal = MatrixAllocateAndZero(maxEinKV,mDetectorCellCount);
  mIncoherentScatterSignal = MatrixAllocateAndZero(maxEinKV,mDetectorCellCount);
  mPrimarySignal = VecAllocate(mDetectorCellCount);
  mMCSignal = VecAllocate(mDetectorCellCount);
  mSDD = SDD;
  mSID = SID;
  mCellWidth = cellWidth;
  mCellHeight = cellHeight;
  mNRows = nrows;
  mNCols = ncols;
  mCenterRow = centerrow;
  mCenterCol = centercol;
  mGridHeight = gridHeight;
  mNRowsDecimated = effrows;
  mNColsDecimated = effcols;
  mCenterRow = centerrow;
  mCenterCol = centercol;
}

FocallyAlignedXCollimatedDetector::~FocallyAlignedXCollimatedDetector() {
  if (mDetectorCenters)                delete[] mDetectorCenters; 
  if (mDetectorNormals)                delete[] mDetectorNormals; 
  if (mDetectorToCollimatorUnitVec)    delete[] mDetectorToCollimatorUnitVec; 
  if (mDetectorAlongCollimatorUnitVec) delete[] mDetectorAlongCollimatorUnitVec; 
  if (mDetectorCellAreas)              VecFree(mDetectorCellAreas); 
  if (mCoherentScatterSignal)          MatrixFree(mCoherentScatterSignal); 
  if (mIncoherentScatterSignal)        MatrixFree(mIncoherentScatterSignal); 
  if (mPrimarySignal)                  VecFree(mPrimarySignal); 
  if (mMCSignal)                       VecFree(mMCSignal); 
}


float FocallyAlignedXCollimatedDetector::GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir) {
  float W = mCellWidth;
  float D = mCellHeight;
  float H = mGridHeight;
  Vec3 d(pos.PointsTo(mDetectorCenters[cellNum]));
  // Inner products of d with the basis vectors
  float N = mDetectorNormals[cellNum].Dot(d);
  float C = mDetectorToCollimatorUnitVec[cellNum].Dot(d);
  float Wprime = W - H*fabs(C/N);
  Wprime = (Wprime < 0) ? 0 : Wprime;
  Wprime = (Wprime > W) ? W : Wprime;
  return fabs(D*Wprime*N);
} 

float FocallyAlignedXandZCollimatedDetector::GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir) {
  float W = mCellWidth;
  float D = mCellHeight;
  float H = mGridHeight;
  Vec3 d(pos.PointsTo(mDetectorCenters[cellNum]));
  // Inner products of d with the basis vectors
  float N = mDetectorNormals[cellNum].Dot(d);
  float C1 = mDetectorToCollimatorUnitVec[cellNum].Dot(d);
  float C2 = mDetectorAlongCollimatorUnitVec[cellNum].Dot(d);

  float Wprime = W - H*fabs(C1/N);
  float Dprime = D - H*fabs(C2/N);

  Wprime = (Wprime < 0) ? 0 : Wprime;
  Dprime = (Dprime < 0) ? 0 : Dprime;
  return fabs(Dprime*Wprime*N);
} 

static void GetPatchCoords(Vec3 x0, Vec3 normal, Vec3 xdir, Vec3 ydir, Vec3 p, Vec3 d, float &px, float &py) {
  Vec3 b(p.x-x0.x,p.y-x0.y,p.z-x0.z);
  float alpha = -normal.Dot(b)/normal.Dot(d);
  px = xdir.Dot(b)+xdir.Dot(d)*alpha;
  py = ydir.Dot(b)+ydir.Dot(d)*alpha;
}

bool FocallyAlignedXCollimatedDetector::RecordPhotonMC(Photon& p, bool isPrimary) {
  // The detector is defined by the points along the source-focused
  // arc -> S + SDD*[sin(alpha),cos(alpha)]
  //     -> [SDD*sin(alpha), SID - SDD*cos(alpha), 0]
  // We want the photon to intersect this arc.  That means we want
  //    p.x + beta*p.dir = [SDD*sin(alpha), SID-SDD*cos(alpha), 0]
  // for some choice of beta, alpha.  To solve this equation, we
  //    px + beta*dirx = SDD*sin(alpha)
  //    py + beta*diry = SID-SDD*cos(alpha)
  //    pz + beta*dirz \in [-deltz,deltz]
  //
  //  So let
  //    tx = px + beta*dirx
  //    ty = SID - py - beta*diry
  //
  //  Then
  //    tx^2+ty^2 = SDD^2
  //
  //  We can root this equation to find the values of beta for
  // which the photon intersects the detector.
  //
  //   beta^2*(dirx^2+diry^2) + beta*(2*px*dirx+2*(SID-py)*(-diry)) + px^2+(SID-py)^2 - SDD^2 = 0
  //
  if (mDecimated) {
    cerr << "Monte carlo mode is not supported with decimated detectors!\n";
    exit(1);
  }
  if (p.GetEnergy()==0) return false; // if the photon has already been absorbed in the phantom
  Vec3 ppos = p.GetPosition();
  Vec3 opos = p.GetPosition();
  Vec3 pdir = p.GetDirection();
  double a, b, c;
  a = pdir.x*pdir.x+pdir.y*pdir.y;
  b = 2*(pdir.x*ppos.x+(mSID-ppos.y)*(-pdir.y));
  c = ppos.x*ppos.x+(mSID-ppos.y)*(mSID-ppos.y)-mSDD*mSDD;
  double det;
  det = b*b-4*a*c;
  if (det <= 0) return false;
  det = sqrt(det);
  double advance;
  if (-b-det >= 0)
    advance = (-b-det)/(2*a);
  else
    advance = (-b+det)/(2*a);
  p.Advance(advance);
  // The photon in question is now on the cylindrical surface of the detector.
  // position of the photon is [SDD*sin(alpha),SID-SDD*cos(alpha),z]
  // We trim on the z height of the detector
  ppos = p.GetPosition();
  // the z-range spanned by the detector is (0:(mNRows-1) - mCenterRow)*mCellHeight
  // or -mCenterRow*mCellHeight to (mNRows-1-mCenterRow)*mCellHeight
  // So to calculate the z bin, we need (zbin - mCenterRow)*mCellHeight = pos.z
  // or equivalently zbin = mCenterRow + pos.z/mCellHeight
  int zbin = (int) rint(mCenterRow+ppos.z/mCellHeight);
  if ((zbin < 0) || (zbin >= mNRows)) return false;
  // Compute the fan angle of the ray - convert to arc-length
  int alpha = (int) rint(atan2(ppos.x,mSID-ppos.y)*mSDD/mCellWidth+mCenterCol);
  if ((alpha>=0) && (alpha < mNCols)) {
    int hit = alpha+zbin*mNCols; 
    // Test for intersection with the leftcollimator blade
    Vec3 x0(mDetectorCenters[hit]);
    Vec3 c(mDetectorToCollimatorUnitVec[hit]);
    Vec3 a(mDetectorAlongCollimatorUnitVec[hit]);
    Vec3 n(mDetectorNormals[hit]);
    float cx_left, cy_left;
    GetPatchCoords(x0.Minus(c.Scaled(mCellWidth/2)),c,n,a,opos,pdir,cx_left,cy_left);
    // Test for intersection with the rightcollimator blade
    float cx_right, cy_right;
    GetPatchCoords(x0.Minus(c.Scaled(-mCellWidth/2)),c,n,a,opos,pdir,cx_right,cy_right);
    bool lefthit, righthit;
    lefthit = (cx_left >= 0) && (cx_left <= mGridHeight);
    righthit = (cx_right >= 0) && (cx_right <= mGridHeight);
    float update = p.GetEnergy()*p.GetProbability()/1.E3/mNumEvents;
    if ((!lefthit) && (!righthit)) {
      if (!isPrimary) // only record events that has been scattered by the phantom
	mMCSignal[hit] += update;
      else
	mPrimarySignal[hit] += update;
      return true;
    }
  }
  return false;
}



// Included flat detector class - Pablo Milioni (GE Healthcare)

XAlignedZCollimatedDetectorFlat::XAlignedZCollimatedDetectorFlat(float SDD, float SID, 
								     float cellWidth, 
								     float cellHeight,
								     int nrows, int ncols, 
								     float centerrow, 
								     float centercol,
								     float gridHeight,
								     int coldecimation, 
								     int rowdecimation,
								     int maxEinKV) {
  // Count the actual number of rows and columns
  int effrows, effcols;
  effrows = 0;
  for (int i=0;i<nrows;i+=rowdecimation)
    effrows++;
  effcols = 0;
  for (int i=0;i<ncols;i+=coldecimation)
    effcols++;
  mDecimated = ((rowdecimation != 1) || (coldecimation != 1));
  mDetRowDecimationFactor = rowdecimation;
  mDetColDecimationFactor = coldecimation;
  mDetectorCellCount = effrows*effcols;
  mDetectorCenters = new Vec3[mDetectorCellCount];
  mDetectorNormals = new Vec3[mDetectorCellCount];
  mDetectorToCollimatorUnitVec  = new Vec3[mDetectorCellCount];
  mDetectorAlongCollimatorUnitVec = new Vec3[mDetectorCellCount];
  mDetectorCellAreas = VecAllocate(mDetectorCellCount);
  
  
  for (int i=0;i<effrows;i++) {
    for (int j=0;j<effcols;j++) {
      int k=i*effcols+j;
      mDetectorCenters[k].x = (j*coldecimation - centercol)*cellWidth;
      mDetectorCenters[k].y = 0;
      mDetectorCenters[k].z = (i*rowdecimation - centerrow)*cellHeight;
      mDetectorNormals[k].x = 0;
      mDetectorNormals[k].y = 1;
      mDetectorNormals[k].z = 0;
      mDetectorToCollimatorUnitVec[k].x = 0;
      mDetectorToCollimatorUnitVec[k].y = 0;
      mDetectorToCollimatorUnitVec[k].z = 1;
      mDetectorAlongCollimatorUnitVec[k].x = 1;
      mDetectorAlongCollimatorUnitVec[k].y = 0;
      mDetectorAlongCollimatorUnitVec[k].z = 0;
      mDetectorCellAreas[k] = cellWidth*cellHeight;
    }
  }
  mEcount = maxEinKV;
  mCoherentScatterSignal = MatrixAllocateAndZero(maxEinKV,mDetectorCellCount);
  mIncoherentScatterSignal = MatrixAllocateAndZero(maxEinKV,mDetectorCellCount);
  mPrimarySignal = VecAllocate(mDetectorCellCount);
  mMCSignal = VecAllocate(mDetectorCellCount);
  mSDD = SDD;
  mSID = SID;
  mCellWidth = cellWidth;
  mCellHeight = cellHeight;
  mNRows = nrows;
  mNCols = ncols;
  mCenterRow = centerrow;
  mCenterCol = centercol;
  mGridHeight = gridHeight;
  mNRowsDecimated = effrows;
  mNColsDecimated = effcols;
  mCenterRow = centerrow;
  mCenterCol = centercol;
}

XAlignedZCollimatedDetectorFlat::~XAlignedZCollimatedDetectorFlat() {
  if (mDetectorCenters)                delete[] mDetectorCenters; 
  if (mDetectorNormals)                delete[] mDetectorNormals; 
  if (mDetectorToCollimatorUnitVec)    delete[] mDetectorToCollimatorUnitVec; 
  if (mDetectorAlongCollimatorUnitVec) delete[] mDetectorAlongCollimatorUnitVec; 
  if (mDetectorCellAreas)              VecFree(mDetectorCellAreas); 
  if (mCoherentScatterSignal)          MatrixFree(mCoherentScatterSignal); 
  if (mIncoherentScatterSignal)        MatrixFree(mIncoherentScatterSignal); 
  if (mPrimarySignal)                  VecFree(mPrimarySignal); 
  if (mMCSignal)                       VecFree(mMCSignal); 
}

//****************************
// FIX ME: The intersection test is not ready yet.
// Collimator blade vector must change with the row index, alligned with the x-axis passing at the source.
// Forced mGridHeight to 0
//****************************
float XAlignedZCollimatedDetectorFlat::GetEffectiveArea(int cellNum, Vec3 pos, Vec3 dir) {
  float W = mCellWidth;
  float D = mCellHeight;
  float H = mGridHeight;
  H = 0;
  Vec3 d(pos.PointsTo(mDetectorCenters[cellNum]));
  // Inner products of d with the basis vectors
  float N = mDetectorNormals[cellNum].Dot(d);
  float C = mDetectorToCollimatorUnitVec[cellNum].Dot(d);
  float Dprime = D - H*fabs(C/N);
  Dprime = (Dprime < 0) ? 0 : Dprime;
  Dprime = (Dprime > D) ? D : Dprime;
  return fabs(W*Dprime*N);
}

bool XAlignedZCollimatedDetectorFlat::RecordPhotonMC(Photon& p, bool isPrimary) {
  // The detector is defined by the points along the flat detector
  // arc -> [x, 0, z]
  // We want the photon to intersect this plane.  That means we want
  //    px + beta*dirx \in [-deltx,deltx]
  //    py + beta*diry = 0
  //    pz + beta*dirz \in [-deltz,deltz]
  //
  //So,   
  //   beta = (SID-SDD-py)/diry;
  //
  
  if (mDecimated) {
    cerr << "Monte carlo mode is not supported with decimated detectors!\n";
    exit(1);
  }
  Vec3 ppos = p.GetPosition();
  Vec3 opos = p.GetPosition();
  Vec3 pdir = p.GetDirection();
//   std::cout << " photon position : <"<<ppos.x<<ppos.y<<ppos.z<<">"<<endl;
  double advance = -ppos.y/pdir.y;
  p.Advance(advance);
  // The photon in question is now on the flat surface of the detector.
  // position of the photon is [x,0,z]
  // We trim on the z height of the detector
  ppos = p.GetPosition();
  // the z-range spanned by the detector is (0:(mNRows-1) - mCenterRow)*mCellHeight
  // or -mCenterRow*mCellHeight to (mNRows-1-mCenterRow)*mCellHeight
  // So to calculate the z bin, we need (zbin - mCenterRow)*mCellHeight = pos.z
  // or equivalently zbin = mCenterRow + pos.z/mCellHeight
  int zbin = (int) rint(mCenterRow+ppos.z/mCellHeight);
  int xbin = (int) rint(mCenterCol+ppos.x/mCellWidth);
  //std::cout << " checked detector index... xbin = "<<xbin<<" zbin = "<<zbin<<"hit = "<<xbin+zbin*mNCols<<endl;
  
  if ((zbin >= 0) && (zbin < mNRows) && (xbin >= 0) && (xbin < mNCols)) {
    int hit = xbin+zbin*mNCols; 
    
    //****************************
    //FIX ME: The intersection test is not ready yet. Please make sure mGridHeight = 0
    //****************************
//     // Test for intersection with the leftcollimator blade
//     Vec3 x0(mDetectorCenters[hit]);
//     Vec3 c(mDetectorToCollimatorUnitVec[hit]);
//     Vec3 a(mDetectorAlongCollimatorUnitVec[hit]);
//     Vec3 n(mDetectorNormals[hit]);
//     float cx_left, cy_left;
//     GetPatchCoords(x0.Minus(c.Scaled(mCellWidth/2)),c,n,a,opos,pdir,cx_left,cy_left);
//     // Test for intersection with the rightcollimator blade
//     float cx_right, cy_right;
//     GetPatchCoords(x0.Minus(c.Scaled(-mCellWidth/2)),c,n,a,opos,pdir,cx_right,cy_right);
//     bool lefthit, righthit;
//     lefthit = (cx_left >= 0) && (cx_left <= mGridHeight);
//     righthit = (cx_right >= 0) && (cx_right <= mGridHeight);
     float update = p.GetEnergy()*p.GetProbability()/1.E3/mNumEvents;
//     if ((!lefthit) && (!righthit)) {
      if (!isPrimary)
	mMCSignal[hit] += update;
      else
	mPrimarySignal[hit] += update;
      return true;
    //}
  }
  return false;
}


