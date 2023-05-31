// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Transformation.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.1415926536
#endif

Transformation::Transformation() {
  forward.SetIdentity();
  backward.SetIdentity();
  maxScale = 1.0;
}

Transformation::~Transformation() {
}

void Transformation::Translate(Vec3 delta) {
  char buffer[4096];
  sprintf(buffer,"Transform { translation %f %f %f\n",delta.x/1000,delta.y/1000,delta.z/1000);
  VRMLCode.push_back(buffer);
  forward = forward.Mat4MultiplyLeft(Mat4::Translate(delta.Negate()));
  backward = backward.Mat4MultiplyRight(Mat4::Translate(delta));
}

void Transformation::RotateX(double ang) {
  char buffer[4096];
  sprintf(buffer,"Transform { rotation 1 0 0 %f\n",ang/180.0*M_PI);
  VRMLCode.push_back(buffer);
  forward = forward.Mat4MultiplyLeft(Mat4::RotateX(ang));
  backward = backward.Mat4MultiplyRight(Mat4::RotateX(-ang));
}

void Transformation::RotateY(double ang) {
  char buffer[4096];
  sprintf(buffer,"Transform { rotation 0 1 0 %f\n",ang/180.0*M_PI);
  VRMLCode.push_back(buffer);
  forward = forward.Mat4MultiplyLeft(Mat4::RotateY(ang));
  backward = backward.Mat4MultiplyRight(Mat4::RotateY(-ang));
}

void Transformation::RotateZ(double ang) {
  char buffer[4096];
  sprintf(buffer,"Transform { rotation 0 0 1 %f\n",ang/180.0*M_PI);
  VRMLCode.push_back(buffer);
  forward = forward.Mat4MultiplyLeft(Mat4::RotateZ(ang));
  backward = backward.Mat4MultiplyRight(Mat4::RotateZ(-ang));
}

void Transformation::Scale(Vec3 gain) {
  char buffer[4096];
  sprintf(buffer,"Transform { scale %f %f %f\n",gain.x/1000,gain.y/1000,gain.z/1000);
  VRMLCode.push_back(buffer);
  forward = forward.Mat4MultiplyLeft(Mat4::Scale(gain));
  backward = backward.Mat4MultiplyRight(Mat4::Scale(gain.Recip()));
  maxScale *= gain.Max();
}

Vec3 Transformation::Forward(Vec3 x) {
  return forward.Vec3Multiply(x);
}

Vec3 Transformation::Backward(Vec3 x) {
  return backward.Vec3Multiply(x);
}

Vec3 Transformation::ForwardDelta(Vec3 x) {
  return forward.Vec3MultiplyZeroConst(x);
}

Vec3 Transformation::BackwardDelta(Vec3 x) {
  return backward.Vec3MultiplyZeroConst(x);
}

double Transformation::GetMaxScale() {
  return maxScale;
}

void Transformation::WriteVRMLCode(std::ostream& c) {
  for (unsigned i=0;i<VRMLCode.size();i++) {
    c << VRMLCode[VRMLCode.size()-1-i];
    c << " children [\n";
  }
}

void Transformation::CloseVRMLCode(std::ostream& c) {
  for (unsigned i=0;i<VRMLCode.size();i++) 
    c << "]}";
  c << "\n";
}

void Transformation::PrintMe(std::ostream& c) {
  c << "   Forward transform:\n";
  c << forward;
  c << "   Backward transform:\n";
  c << backward;
}

std::ostream& operator<<(std::ostream& c, Transformation& p) {
  p.PrintMe(c);
  return c;
}


