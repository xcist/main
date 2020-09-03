// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef BaseObject_h
#define BaseObject_h

#include <vector>
#include <iostream>
#include "MatVec.h"
#include "Transformation.h"
#include "ClipPlane.h"
#include "Intersection.h"

class BaseObject {
 protected:
  Transformation xform;
  ClipPlaneSet clips;
  int my_ID;
  int material;
  double density;
  int priority;
  Vec3 boundCenter;
 public:
  double boundRadius;
  float bx, by, bz;
  // Describe me in VRML syntax
  virtual void WriteVRML(std::ostream& c);
  // Describe me to stdout
  virtual void PrintMe(std::ostream& c);
  // Calculate the object-path intersection
  virtual void ObjectPathIntersect(Vec3 src, Vec3 delta, IntersectionSet& T) = 0;
  // Calculate the object-voxel intersection
  virtual IntersectionSet ObjectVoxelIntersect(Vec3 src, Vec3 det) = 0;
  // Returns the bounding radius for this object
  virtual double CalculateBoundingRadius() = 0;
  // Updates our bounding sphere calculation
  void UpdateBoundingSphere();
  // Test for containment of the bounding boxes
  bool Contains(BaseObject* child);
  // Perform the boundary box test
  bool BoundarySphereIntersect(float *src, float *direction, float SDD);
  // Calculate the path intersection (includes clipping plane calculations)
  void PathIntersect(float* src, float* direction, IntersectionSet& T);
  // Calculate a voxel intersection
  IntersectionSet VoxelIntersect(Vec3 src);
  // Set/Get the object material
  int Material();
  void Material(int);
  // Set/Get the density of the object
  void Density(double rho);
  double Density(void);
  // Set/Get the priority of the object
  void Priority(int prio);
  int Priority(void);
  // Set/Get the ID of the object
  void ID(int id);
  int ID(void);
  // Set/Get the objects transform
  Transformation& Transform();
  void Transform(Transformation m);
  // Add a clipping plane (specified by a normal, a minimum distance to the origin, and a point inside the clipping region
  void AddClipPlane(ClipPlane p);
  Vec3 BoundCenter();
  double BoundRadius();
  // Default constructor
  BaseObject();
  virtual ~BaseObject();
};

typedef std::vector<std::string> MaterialTable;

std::ostream& operator<<(std::ostream& c, BaseObject& p);
std::ostream& operator<<(std::ostream& c, MaterialTable& p);
#endif

