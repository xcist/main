// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __Objects_hpp__
#define __Objects_hpp__
#include <stdio.h>

#define MAXDIST 4
#define MAXCLIP 4

class vector3{
 public:
  double x;
  double y;
  double z;

  void printMe() { printf("<%f %f %f>",x,y,z); }
};

class ray3 {
public:
  vector3 source;
  vector3 direction;

  void printMe() { 
    printf(" Source: "); 
    source.printMe();
    printf("\n Direction: ");
    direction.printMe();
    printf("\n");
  }
  
};

class boundingBox {
 public:
  double minX;
  double minY;
  double minZ;
  double maxX;
  double maxY;
  double maxZ;
  
  bool doesIntersect(ray3& t);
  void printMe() {
    printf("Bounding box: ");
    printf("  %f < x < %f\n",minX,maxX);
    printf("  %f < y < %f\n",minY,maxY);
    printf("  %f < z < %f\n",minZ,maxZ);
  }
};

class intersectionList {
 public:
  int intersectionCount;
  double dist[MAXDIST];
};

class clipPlane {
 public:
  double normalX;
  double normalY;
  double normalZ;
  double distance;
};

class baseObject {
 public:
  double centerX;
  double centerY;
  double centerZ;
  double density;
  clipPlane clips[MAXCLIP];
  int clipCount;

  baseObject() {centerX = centerY = centerZ = density = 0.0f; clipCount = 0;}
  virtual boundingBox getBoundingBox() = 0;
  virtual intersectionList getIntersectionList(ray3&) = 0;
};

class sphereObject : public baseObject {
public:
  double radius;
  
  sphereObject() {radius = 0.0f;}
  virtual boundingBox getBoundingBox();
  virtual intersectionList getIntersectionList(ray3&);
};

#endif

