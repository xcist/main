// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#ifndef TreePhantom_h
#define TreePhantom_h

#include "BaseObject.h"
// Routines to manipulate the phantoms
// The linear phantom
typedef std::vector<BaseObject*> LinearPhantom;
std::ostream& operator<<(std::ostream& c, LinearPhantom& p);

class TreePhantom {
 public:
  BaseObject *obj;
  TreePhantom* right;
  TreePhantom* down;
  TreePhantom();
  ~TreePhantom();
  static TreePhantom* BuildTreePhantomFromLinear(LinearPhantom &p);
};

std::ostream& operator<<(std::ostream& c, TreePhantom* p);
void TreePathIntersect(TreePhantom* tree, float* src, float* det, float* dir, 
		       float SDD, IntersectionSet& T);
void TreePathTubeTest(TreePhantom* tree, float radius, float* src, float* dir, 
		      float SDD, LinearPhantom& P);
#endif

