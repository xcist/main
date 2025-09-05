// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#include "Sphere.h"
#include <cmath>

Sphere::Sphere() : BaseObject() {
}

Sphere::~Sphere() {
}

void Sphere::WriteVRML(std::ostream& c) {
  xform.WriteVRMLCode(c);
  c << "Shape { appearance USE App1 geometry Sphere { radius 1.0 } }\n";
  xform.CloseVRMLCode(c);  
}

void Sphere::PrintMe(std::ostream& c) {
  c << "Object Sphere\n";
  BaseObject::PrintMe(c);
}

IntersectionSet Sphere::ObjectVoxelIntersect(Vec3 src, Vec3 det) {
	return IntersectionSet();
}

void Sphere::ObjectPathIntersect(Vec3 ct, Vec3 delta, IntersectionSet& T) {
  double coefA = delta.Dot(delta);
  double coefB = 2.0*(delta.Dot(ct));
  double coefC = ct.Dot(ct)-1.0;
  double discrim = coefB*coefB - 4.0*coefA*coefC;
  if (discrim >= 0.0) {
    T.push_back(Intersection(my_ID,(-coefB - sqrt(discrim))/(2.0*coefA),true,
			     density,priority,material));
    T.push_back(Intersection(my_ID,(-coefB + sqrt(discrim))/(2.0*coefA),false,
			     density,priority,material));
  }
}

double Sphere::CalculateBoundingRadius() {
  return 1.0f;
}


