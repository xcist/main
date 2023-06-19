// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "ClipPlane.h"
#include <cmath>

ClipPlane::ClipPlane() {
  intersect = 0;
}

ClipPlane::~ClipPlane() {
}

ClipPlane::ClipPlane(Vec3 anormal, double aintersect) {
  normal = anormal;
  intersect = aintersect;
}

void ClipPlane::PrintMe(std::ostream& c) {
  c << "Clip Plane normal: " << normal << "\n";
  c << "Clip Plane intersection: " << intersect << "\n";
}

std::ostream& operator<<(std::ostream& c, ClipPlane& p) {
  p.PrintMe(c);
  return c;
}

std::ostream& operator<<(std::ostream& c, ClipPlaneSet& p) {
  for (unsigned i=0;i<p.size();i++) {
    p[i].PrintMe(c);
    c << "\n";
  }
  return c;
}

// Returns true if there is an intersection, sets intersection time to the
// time of intersection of the plane and the ray, and sets direction to +/- 1 to
// indicate which half of the real line to include
void ClipPlane::GetIntersection(Vec3 src, Vec3 delta, IntersectionSet &T) {
  // Are we a positive plane or a negative one?
  double t0;
  int direction;
  // Get the sine of the angle between delta and the plane normal
  double sinang;
  sinang = delta.Dot(normal);
  if (fabs(sinang) < 1e-16) {
    // Ray is essentially parallel to the plane... 
    if (src.Dot(normal) > intersect) {
      // Ray is wholly contained in the half space
      t0 = -INFTY;
      direction = +1;
    } else {
      // Ray is wholly discarded by the half space
      t0 = -INFTY;
      direction = -1;
    }
  } else {
    // Ray is not parallel to the plane.  Compute the interesction time
    t0 = (intersect - src.Dot(normal))/sinang;
    if (sinang > 0) direction = +1; else direction = -1;
  }
  // Update the intersection set... we do this by testing each intersection
  // time in the set against our t0.  If we are on the clipped side of the
  // plane, the time is modified.
  IntersectionSet Copy;
  for (IntersectionSet::iterator i=T.begin();i!=T.end();i++) {
    if ((i->time < t0) && (direction == 1)) {
      Intersection n = *i;
      n.time = t0;
      Copy.push_back(n);
    } else if ((i->time > t0) && (direction == -1)) {
      Intersection n = *i;
      n.time = t0;
      Copy.push_back(n);
    } else {
      Copy.push_back(*i);
    }
  }
  T = Copy;
}



