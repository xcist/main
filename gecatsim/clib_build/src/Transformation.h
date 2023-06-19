// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef Transformation_h
#define Transformation_h

#include <iostream>
#include "MatVec.h"
#include <vector>
#include <string>
#include <stdio.h>

class Transformation {
  Mat4 forward;
  Mat4 backward;
  double maxScale;
  std::vector<std::string> VRMLCode;
 public:
  Transformation();
  ~Transformation();
  void PrintMe(std::ostream& c);
  void Scale(Vec3 gain);
  void Translate(Vec3 delta);
  void RotateX(double ang);
  void RotateY(double ang);
  void RotateZ(double ang);
  Vec3 Forward(Vec3 x);
  Vec3 Backward(Vec3 x);
  Vec3 ForwardDelta(Vec3 x);
  Vec3 BackwardDelta(Vec3 x);
  double GetMaxScale();
  void WriteVRMLCode(std::ostream& c);
  void CloseVRMLCode(std::ostream& c);
};

std::ostream& operator<<(std::ostream& c, Transformation& p);
#endif

