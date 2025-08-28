// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#ifndef Cube_h
#define Cube_h

#include "BaseObject.h"
#include "MatVec.h"

class Cube : public BaseObject {
 public:
  // Describe me to stdout
  virtual void PrintMe(std::ostream& c);
  // Calculate the object-path intersection
  virtual void ObjectPathIntersect(Vec3 src, Vec3 delta, IntersectionSet& T);
  // Calculate the object-voxel intersection
  virtual IntersectionSet ObjectVoxelIntersect(Vec3 src, Vec3 det);
  // Calculate the bounding sphere for the object
  virtual double CalculateBoundingRadius();
  Cube();
  ~Cube();
  void WriteVRML(std::ostream& c);
 private:
};

#endif

