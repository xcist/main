// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "MatVec.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.1415926535
#endif

Vec3::Vec3() {
  x = 0; y = 0; z = 0;
}

Vec3::Vec3(double xa, double ya, double za)  {
  x = xa; 
  y = ya;
  z = za;
}

Vec3::~Vec3() {
}

Vec3 Vec3::Negate() {
  return Vec3(-x,-y,-z);
}

Vec3 Vec3::Recip() {
  if (x == 0 || y == 0 || z == 0)
    throw "zero encountered in Vec3::Recip!";
  Vec3 C;
  C.x = 1.0/x;
  C.y = 1.0/y;
  C.z = 1.0/z;
  return C;
}

Vec3 Vec3::Subtract(const Vec3& b) {
  Vec3 C;
  C.x = x - b.x;
  C.y = y - b.y;
  C.z = z - b.z;
  return C;
}

void Vec3::Normalize() {
  double nrm;
  nrm = Norm();
  if (nrm == 0)
    throw "Encountered zero-length vector in Vec3::Normalize routine - probably a source and detector are co-located?";
  x /= nrm;
  y /= nrm;
  z /= nrm;
}

double Vec3::Dot(const Vec3& b) {
  return x*b.x + y*b.y + z*b.z;
}

void Vec3::PrintMe(std::ostream& c) {
  c << "<" << x << "," << y << "," << z << ">";
}

double Vec3::Norm() {
  return sqrt(x*x+y*y+z*z);
}

Vec3 Vec3::Direction(const Vec3& b) {
  Vec3 C;
  C.x = x - b.x;
  C.y = y - b.y;
  C.z = z - b.z;
  double nrm;
  nrm = C.Norm();
  C.x /= nrm;
  C.y /= nrm;
  C.z /= nrm;
  return C;
}

Vec3 Vec3::InverseScale(const Vec3& b) {
  if (b.x == 0 || b.y == 0 || b.z == 0)
    throw "Encountered zero in Vec3::InverseScale routine - probably an object dimension is zero?";
  Vec3 C;
  C.x = x/b.x;
  C.y = y/b.y;
  C.z = z/b.z;
  return C;
}

double Vec3::Max() {
  return std::max(std::max(x,y),z);
}

Vec3 Vec3::Cross(Vec3 b) {
  Vec3 C;
  C.x = y*b.z - z*b.y;
  C.y = - x*b.z + z*b.x;
  C.z = x*b.y - y*b.x;
  return C;
}

void Vec3::RotateVector(float theta, float phi) {
  //   G4double dirx = sinTheta * std::cos(phi);
  //   G4double diry = sinTheta * std::sin(phi);
  //   G4double dirz = cosTheta ;
  //   G4ThreeVector photonDirection1(dirx,diry,dirz);
  //   photonDirection1.rotateUz(photonDirection0);
  //   aParticleChange.ProposeMomentumDirection(photonDirection1) ;

  double dx = sin(theta)*cos(phi);
  double dy = sin(theta)*sin(phi);
  double dz = cos(theta);

  double u1 = x;
  double u2 = y;
  double u3 = z;

  double up = u1*u1+u2*u2;
  if (up > 0) {
    up = sqrt(up);
    double px = dx;
    double py = dy;
    double pz = dz;
    dx = (u1*u3*px - u2*py)/up + u1*pz;
    dy = (u2*u3*px + u1*py)/up + u2*pz;
    dz =    -up*px +             u3*pz;
  }
  else if (u3 < 0.) { dx = -dx; dz = -dz; }      // phi=0  teta=pi
  else {};
  
  x = dx;
  y = dy;
  z = dz;
}


Mat4::Mat4() {
  Zero();
}

void Mat4::Zero() {
  int i, j;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      coeff[i][j] = 0.0;
}

void Mat4::SetIdentity() {
  Zero();
  for (int i=0;i<4;i++)
    coeff[i][i] = 1.0;
}

Vec3 Mat4::Vec3Multiply(Vec3 B) {
  Vec3 R;
  R.x = coeff[0][0]*B.x + coeff[0][1]*B.y + coeff[0][2]*B.z + coeff[0][3];
  R.y = coeff[1][0]*B.x + coeff[1][1]*B.y + coeff[1][2]*B.z + coeff[1][3];
  R.z = coeff[2][0]*B.x + coeff[2][1]*B.y + coeff[2][2]*B.z + coeff[2][3];
  return R;
}

Vec3 Mat4::Vec3MultiplyZeroConst(Vec3 B) {
  Vec3 R;
  R.x = coeff[0][0]*B.x + coeff[0][1]*B.y + coeff[0][2]*B.z;
  R.y = coeff[1][0]*B.x + coeff[1][1]*B.y + coeff[1][2]*B.z;
  R.z = coeff[2][0]*B.x + coeff[2][1]*B.y + coeff[2][2]*B.z;
  return R;
}

Mat4 Mat4::Mat4MultiplyRight(Mat4 B) {
  Mat4 C;
  for (int i=0;i<4;i++)
    for (int j=0;j<4;j++)
      for (int k=0;k<4;k++)
	C.coeff[i][j] += coeff[i][k]*B.coeff[k][j];
  return C;
}

Mat4 Mat4::Mat4MultiplyLeft(Mat4 B) {
  Mat4 C;
  for (int i=0;i<4;i++)
    for (int j=0;j<4;j++)
      for (int k=0;k<4;k++)
	C.coeff[i][j] += B.coeff[i][k]*coeff[k][j];
  return C;  
}


Mat4 Mat4::RotateX(double ang) {
  ang = -ang*M_PI/180;
  Mat4 C;
  C.SetIdentity();
  C.coeff[1][1] = cos(ang);
  C.coeff[2][2] = cos(ang);
  C.coeff[2][1] = -sin(ang);
  C.coeff[1][2] = sin(ang);
  return C;
}

Mat4 Mat4::RotateY(double ang) {
  Mat4 C;
  ang = -ang*M_PI/180;
  C.SetIdentity();
  C.coeff[0][0] = cos(ang);
  C.coeff[2][2] = cos(ang);
  C.coeff[0][2] = -sin(ang);
  C.coeff[2][0] = sin(ang);
  return C;
}

Mat4 Mat4::RotateZ(double ang) {
  Mat4 C;
  ang = -ang*M_PI/180;
  C.SetIdentity();
  C.coeff[0][0] = cos(ang);
  C.coeff[1][1] = cos(ang);
  C.coeff[0][1] = sin(ang);
  C.coeff[1][0] = -sin(ang);
  return C;
}

Mat4 Mat4::Translate(Vec3 delta) {
  Mat4 C;
  C.SetIdentity();
  C.coeff[0][3] = -delta.x;
  C.coeff[1][3] = -delta.y;
  C.coeff[2][3] = -delta.z;
  return C;
}

Mat4 Mat4::Scale(Vec3 gain) {
  Mat4 C;
  C.SetIdentity();
  C.coeff[0][0] = gain.x;
  C.coeff[1][1] = gain.y;
  C.coeff[2][2] = gain.z;
  return C;
}

Mat4::~Mat4() {
}

void Mat4::PrintMe(std::ostream& c) {
  int i, j;
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      c << coeff[i][j] << "  ";
    c << "\n";
  }
}

std::ostream& operator<<(std::ostream& c, Vec3 p) {
  p.PrintMe(c);
  return c;
}

std::ostream& operator<<(std::ostream& c, Mat4 p) {
  p.PrintMe(c);
  return c;
}


