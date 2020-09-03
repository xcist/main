// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Photon.hpp"
#include <math.h>

#define TMax(a,b) ((a) > (b) ? (a) : (b))
#define TMin(a,b) ((a) < (b) ? (a) : (b))

#define Swap(a,b) {float t = a; a = b; b = t;}
#define sgn(a) ((a) < 0 ? -1 : 1)
#define ispos(a) ((a) > 0 ? 1 : 0)
#define isneg(a) ((a) < 0 ? 1 : 0)

void Photon::Print(std::ostream& o) {
  o << "Pos: " << pos << "(" << sqrt(pos.x*pos.x+pos.y*pos.y) << ")  Dir: " << dir << " Energy " << Energy << "\r\n";
}

Photon::Photon() {
  Probability = 1;
  escaped = false;
  interacted = false;
}

void Photon::SetPosition(Vec3 ipos) {
  pos = ipos;
}
void Photon::SetDirection(Vec3 idir) {
  dir = idir;
}

void Photon::SetEnergy(float E) {
  Energy = E;
}

void Photon::Advance(float dist) {
  pos.x += dir.x*dist;
  pos.y += dir.y*dist;
  pos.z += dir.z*dist;
}

void Photon::AdvanceToPhantom(Phantom* phantom) {
  // 
  // Given the phantom's current position at pos, direction dir,
  // we want to advance our position (without changing momentum)
  // so that we are at the face of the voxel cube.  The face
  // of the cube is defined by
  //     abs(z) = zcount*voxsize/2
  //     abs(x) = xycount*voxsize/2
  //     abs(y) = xycount*voxsize/2
  //This advancement
  // is done by computing the intervals:
  //   (pos + txmin*dir)_x = -xycount*voxsize/2
  //   (pos + txmax*dir)_x = xycount*voxsize/2
  // and then taking the intersection of these intervals.
  mPhantom = phantom;
  voxsize = phantom->GetVoxelSize();
  voxex_x = phantom->GetXExtent();
  voxex_z = phantom->GetZExtent();
  xycount = phantom->GetXYCount();
  zcount = phantom->GetZCount();
  float txmin = ( voxex_x - pos.x)/dir.x;
  float txmax = (-voxex_x - pos.x)/dir.x;
  //  float txmax = (voxex_x+(xycount-1)*voxsize - pos.x)/dir.x;
  if (dir.x < 0)
    Swap(txmin,txmax);
  float tymin = ( voxex_x - pos.y)/dir.y;
  float tymax = (-voxex_x - pos.y)/dir.y;
  //  float tymax = (voxex_x+(xycount-1)*voxsize - pos.y)/dir.y;
  if (dir.y < 0)
    Swap(tymin,tymax);
  float tzmin = ( voxex_z - pos.z)/dir.z;
    float tzmax = (-voxex_z - pos.z)/dir.z;
    //  float tzmax = (voxex_z+(zcount-1)*voxsize - pos.z)/dir.z;
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
  Advance(ttotalmin);
  m = (int) floor((pos.x-voxex_x)/voxsize);
  n = (int) floor((pos.y-voxex_x)/voxsize);
  p = (int) floor((pos.z-voxex_z)/voxsize);
  if (txmin == ttotalmin && ispos(dir.x) == 1) m = 0; 
  if (tymin == ttotalmin && ispos(dir.y) == 1) n = 0;
  if (tzmin == ttotalmin && ispos(dir.z) == 1) p = 0; 
  if (txmin == ttotalmin && isneg(dir.x) == 1) m = xycount - 1; 
  if (tymin == ttotalmin && isneg(dir.y) == 1) n = xycount - 1;
  if (tzmin == ttotalmin && isneg(dir.z) == 1) p =  zcount - 1; 
  if (m<0) m = 0;
  if (n<0) n = 0;
  if (p<0) p = 0;
  UpdateDeltas();
}

void Photon::UpdateDeltas() {
  deltam = sgn(dir.x);
  deltan = sgn(dir.y);
  deltap = sgn(dir.z);
  deltatx = fabs(voxsize/dir.x);
  if (dir.x == 0) deltatx = 1e10;
  deltaty = fabs(voxsize/dir.y);
  if (dir.y == 0) deltaty = 1e10;
  deltatz = fabs(voxsize/dir.z);  
  if (dir.z == 0) deltatz = 1e10;
  txnext = (voxex_x+(ispos(dir.x)+m)*voxsize-pos.x)/dir.x;
  if (dir.x == 0) txnext = 1e10;
  tynext = (voxex_x+(ispos(dir.y)+n)*voxsize-pos.y)/dir.y;
  if (dir.y == 0) tynext = 1e10;
  tznext = (voxex_z+(ispos(dir.z)+p)*voxsize-pos.z)/dir.z;
  if (dir.z == 0) tznext = 1e10;
  txnext = TMax(0,txnext);
  tynext = TMax(0,tynext);
  tznext = TMax(0,tznext);
}

void Photon::AdvanceVoxel() {
  if ((txnext <= tynext) && (txnext <= tznext)) {
    Advance(txnext);
    tynext -= txnext;
    tznext -= txnext;
    txnext = deltatx;
    m += deltam;
  } else if ((tynext <= txnext) && (tynext <= tznext)) {
    Advance(tynext);
    txnext -= tynext;
    tznext -= tynext;
    tynext = deltaty;
    n += deltan;
  } else {
    Advance(tznext);
    txnext -= tznext;
    tynext -= tznext;
    tznext = deltatz;
    p += deltap;
  }
  if ((m < 0) || (m >= xycount) || (n < 0) || (n >= xycount) || 
      (p < 0) || (p >= zcount))
    escaped = true;
}


float Photon::GetStepSize() {
  float step = TMin(TMin(txnext,tynext),tznext);
  return step;
}

std::ostream& operator<<(std::ostream&o, Photon& v) {
  v.Print(o);
  return o;
}



