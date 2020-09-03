// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef ClipPlane_h
#define ClipPlane_h

#include <iostream>
#include <vector>
#include "MatVec.h"
#include "Intersection.h"

class ClipPlane {
 public:
  Vec3 normal;
  double intersect;
  ClipPlane();
  ~ClipPlane();
  ClipPlane(Vec3 anormal, double intersect);
  void PrintMe(std::ostream& c);
  void GetIntersection(Vec3 src, Vec3 delta, IntersectionSet &T);
};

typedef std::vector<ClipPlane> ClipPlaneSet;

std::ostream& operator<<(std::ostream& c, ClipPlane& p);
std::ostream& operator<<(std::ostream& c, ClipPlaneSet& p);

#endif

