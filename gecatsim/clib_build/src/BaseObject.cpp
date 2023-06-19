// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "BaseObject.h"
#include <algorithm>

void BaseObject::WriteVRML(std::ostream& c) {
}

void BaseObject::PrintMe(std::ostream& c) {
  c << "Base object data: " << my_ID << "\n";
  c << "Transform: \n";
  c << xform;
  c << "Clipping planes: \n";
  c << clips;
  c << "Bounding sphere center: \n";
  c << boundCenter;
  c << "\nBounding sphere radius: " << boundRadius << "\n";
  c << "Material: " << material << "\n";
  c << "Density: " << density << "\n";
  c << "Priority: " << priority << "\n";
}

std::ostream& operator<<(std::ostream& c, BaseObject& p) {
  p.PrintMe(c);
  return c;
}

void BaseObject::UpdateBoundingSphere() {
  boundRadius = CalculateBoundingRadius()*xform.GetMaxScale();
  Vec3 iso;
  boundCenter = xform.Forward(iso);
  bx = boundCenter.x;
  by = boundCenter.y;
  bz = boundCenter.z;
}

bool BaseObject::BoundarySphereIntersect(float *src, float *direction, float SDD) {
  bool retval;
  float x1, y1, z1;
  float mindist, ip;
  x1 = src[0] - boundCenter.x;
  y1 = src[1] - boundCenter.y;
  z1 = src[2] - boundCenter.z;
  mindist = x1*x1 + y1*y1 + z1*z1;
  ip = -(direction[0]*x1+direction[1]*y1+direction[2]*z1);
  mindist -= ip*ip;
  retval = ((mindist <= (boundRadius*boundRadius)) && (ip >= 0) && (ip <= SDD));
  return retval; 
}

bool BaseObject::Contains(BaseObject* child) {
  Vec3 tmp;
  tmp = boundCenter.Subtract(child->boundCenter);
  double dist;
  dist = tmp.Norm();
  return ((dist+child->boundRadius) < boundRadius);
}

// Calculate the path intersection (includes clipping plane calculations)
void BaseObject::PathIntersect(float* src, float *direction, IntersectionSet& T) {
  // Transform the source and unit vector...
  Vec3 delta(xform.BackwardDelta(Vec3(direction[0],direction[1],direction[2])));
  Vec3 ct(xform.Backward(Vec3(src[0],src[1],src[2])));
  // Get the object-path intersection
  IntersectionSet S;
  ObjectPathIntersect(ct, delta, S);
  // Loop over the clip planes
  for (unsigned i=0;i<clips.size();i++) {
    clips[i].GetIntersection(Vec3(src[0],src[1],src[2]),Vec3(direction[0],direction[1],direction[2]),S);
  }
  // Add these intersections to the master list
  for (IntersectionSet::iterator j=S.begin(); j!=S.end();j++)
    T.push_back(*j);
}

// Calculate a voxel intersection
IntersectionSet BaseObject::VoxelIntersect(Vec3 src) {
	return IntersectionSet();
}

void BaseObject::Material(int str) {
  material = str;
}

int BaseObject::Material() {
  return material;
}

void BaseObject::Density(double rho) {
  density = rho;
}

double BaseObject::Density() {
  return density;
}

void BaseObject::Priority(int prio) {
  priority = prio;
}

int BaseObject::Priority() {
  return priority;
}

void BaseObject::ID(int id) {
  my_ID = id;
}

int BaseObject::ID() {
  return my_ID;
}

Transformation& BaseObject::Transform() {
  return xform;
}

void BaseObject::Transform(Transformation m) {
  xform = m;
}

void BaseObject::AddClipPlane(ClipPlane p) {
  clips.push_back(p);
}

// Default constructor
BaseObject::BaseObject() {
  boundRadius = 0;
  material = 0;
  density = 1.0;
  priority = 0;
}

BaseObject::~BaseObject() {
}

Vec3 BaseObject::BoundCenter() {
  return boundCenter;
}

double BaseObject::BoundRadius() {
  return boundRadius;
}


