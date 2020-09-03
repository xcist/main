// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef MatVec_h
#define MatVec_h

#include <iostream>

#define INFTY 1e30

class Vec3 {
 public:
  double x;
  double y;
  double z;
  Vec3();
  Vec3(double xa, double ya, double za);
  ~Vec3();
  void PrintMe(std::ostream& c);
  double Norm();
  Vec3 Negate();
  Vec3 Recip();
  Vec3 Subtract(const Vec3& b);
  double Dot(const Vec3& b);
  Vec3 Direction(const Vec3& b);
  Vec3 InverseScale(const Vec3& b);
  double Max();
  void Normalize();
  void RotateVector(float theta, float phi);
  Vec3 PointsTo(Vec3 b) {Vec3 ret(b.x-x,b.y-y,b.z-z); ret.Normalize(); return ret;}
  Vec3 Cross(Vec3 b);
  Vec3 Minus(Vec3 b) {return Vec3(x-b.x,y-b.y,z-b.z);}
  Vec3 Plus(Vec3 b) {return Vec3(x+b.x,y+b.y,z+b.z);}
  Vec3 Scaled(float a) {return Vec3(a*x,a*y,a*z);}
};

class Mat4 {
 public:
  double coeff[4][4];
  Mat4();
  ~Mat4();
  void Zero();
  void SetIdentity();
  void PrintMe(std::ostream& c);
  Vec3 Vec3Multiply(Vec3 B);
  Vec3 Vec3MultiplyZeroConst(Vec3 B);
  Mat4 Mat4MultiplyRight(Mat4 B);
  Mat4 Mat4MultiplyLeft(Mat4 B);
  void Decompose(double &scalex,
		 double &scaley,
		 double &scalez,
		 double &rotx,
		 double &roty,
		 double &rotz,
		 double &transx,
		 double &transy,
		 double &transz);
  static Mat4 RotateX(double ang);
  static Mat4 RotateY(double ang);
  static Mat4 RotateZ(double ang);
  static Mat4 Translate(Vec3 delta);
  static Mat4 Scale(Vec3 gain);
};

std::ostream& operator<<(std::ostream& c, Vec3 p);
std::ostream& operator<<(std::ostream& c, Mat4 p);

#endif

