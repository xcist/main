// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#ifndef Cylinder_h
#define Cylinder_h

#include "BaseObject.h"
#include "MatVec.h"
#include <string.h>

class Cylinder : public BaseObject {
 public:
  // Describe me to stdout
  virtual void PrintMe(std::ostream& c);
  // Calculate the object-path intersection
  virtual void ObjectPathIntersect(Vec3 src, Vec3 det, IntersectionSet& T);
  // Calculate the object-voxel intersection
  virtual IntersectionSet ObjectVoxelIntersect(Vec3 src, Vec3 det);
  // Calculate the bounding sphere for the object
  virtual double CalculateBoundingRadius();
  Cylinder();
  ~Cylinder();
  void WriteVRML(std::ostream& c);
};

#endif

