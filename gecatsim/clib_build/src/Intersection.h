// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef Intersection_h
#define Intersection_h

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <stdio.h>

class Intersection {
 public:
  int ID;
  double time;
  bool enterSense;
  double rho;
  int priority;
  int material;
  Intersection(int id, double a_time, bool enterSense, double a_rho, int a_priority, int a_material);
  Intersection();
  ~Intersection();
  void PrintMe(std::ostream& c) const;
};

class IntersectionSet : public std::vector<Intersection> {
 public:
  void GetHitList(double *ptable);
  void RenderIntersections(int N, float dX, float* IDmap, int ID);
};

std::ostream& operator<<(std::ostream& c, Intersection& p);
std::ostream& operator<<(std::ostream& c, IntersectionSet& p);

#endif

