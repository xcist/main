// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#include "Cube.h"
#include <math.h>

Cube::Cube() : BaseObject() {
}

Cube::~Cube() {
}

void Cube::WriteVRML(std::ostream& c) {
  xform.WriteVRMLCode(c);
  c << "Shape { appearance USE App1 geometry Box { size 1.0 1.0 1.0 } }\n";
  xform.CloseVRMLCode(c);  
}

void Cube::PrintMe(std::ostream& c) {
  c << "Object Cube \n";
  BaseObject::PrintMe(c);
}

double Cube::CalculateBoundingRadius() {
  return sqrt(3.0/4.0);
}

IntersectionSet Cube::ObjectVoxelIntersect(Vec3 src, Vec3 det) {
	return IntersectionSet();
}

void Cube::ObjectPathIntersect(Vec3 src, Vec3 delta, IntersectionSet &T) {
  double times[6];
  char planes[6];
  double coefAX = delta.x;
  double coefAY = delta.y;
  double coefAZ = delta.z;
  double coefBXMinus = src.x + 1.0/2.0;
  double coefBXPlus = src.x - 1.0/2.0;
  double coefBYMinus = src.y + 1.0/2.0;
  double coefBYPlus = src.y - 1.0/2.0;
  double coefBZMinus = src.z + 1.0/2.0;
  double coefBZPlus = src.z - 1.0/2.0;
  // Determine if the ray is parallel to one set of the planes.
  int noIntersectionFlag = 0;
  if (coefAX == 0.0)  // The ray is parallel to the plane.
    {
      if (((-1.0/2.0) <= src.x) && (src.x <= (1.0/2.0)))
	{
	  times[0] = -INFTY;
	  times[1] = INFTY;
	  planes[0] = 'x';
	  planes[1] = 'x';
	}
      else  // The ray is outside the 2 planes.
	{
	  noIntersectionFlag = 1;
	}
    }
  else  // The ray intersects the 2 planes.
    {
      times[0] = -coefBXMinus/coefAX;
      times[1] = -coefBXPlus/coefAX;
      planes[0] = 'x';
      planes[1] = 'x';
    }
  if (!noIntersectionFlag)
    {
      if (coefAY == 0.0)  // The ray is parallel to the planes.
	{
	  if (((-1.0/2.0) <= src.y) && (src.y <= (1.0/2.0)))
	    {
	      times[2] = -2.0 * INFTY;
	      times[3] = 2.0 * INFTY;
	      planes[2] = 'y';
	      planes[3] = 'y';
	    }
	  else  // The ray is outside the planes.
	    {
	      noIntersectionFlag = 1;
	    }
	}
      else  // The ray intersects the 2 planes.
	{
	  times[2] = -coefBYMinus/coefAY;
	  times[3] = -coefBYPlus/coefAY;
	  planes[2] = 'y';
	  planes[3] = 'y';
	}
    }

  if (!noIntersectionFlag)
    {
      if (coefAZ == 0.0)  // The ray is parallel to te 2 planes.
	{
	  if (((-1.0/2.0) <= src.z) && (src.z <= (1.0/2.0)))
	    {
	      times[4] = -3.0 * INFTY;
	      times[5] = 3.0 * INFTY;
	      planes[4] = 'z';
	      planes[5] = 'z';
	    }
	  else   // The ray is outside the 2 planes.
	    {
	      noIntersectionFlag = 1;
	    }
	}
      else  // The ray intersects the 2 planes.
	{
	  times[4] = -coefBZMinus/coefAZ;
	  times[5] = -coefBZPlus/coefAZ;
	  planes[4] = 'z';
	  planes[5] = 'z';
	}
    }
  // Order the times, if necessary.
  if (!noIntersectionFlag)
    {
      for (int i=0;i<5;i++)
	for (int j=0;j<5;j++)
	  {
	    if (times[j] > times[j+1])
	      {
		double tempTime;
		char tempPlane;
		
		tempTime = times[j];
		tempPlane = planes[j];
		times[j] = times[j+1];
		planes[j] = planes[j+1];
		times[j+1] = tempTime;
		planes[j+1] = tempPlane;
	      }
	  }
      // Check for valid combinations.  Basically, the first 3 letters of the "planes" string
      // cannot contain 2 of the same letters.
      if ((planes[0] != planes[1]) && (planes[1] != planes[2]) && (planes[0] != planes[2]))
	{
	  T.push_back(Intersection(my_ID,times[2],true,density,priority,material));
	  T.push_back(Intersection(my_ID,times[3],false,density,priority,material));
	}
    }
}



