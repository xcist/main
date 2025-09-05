// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#include "Cylinder.h"
#include <cmath>

Cylinder::Cylinder() : BaseObject() {
}

Cylinder::~Cylinder() {
}

void Cylinder::WriteVRML(std::ostream& c) {
  xform.WriteVRMLCode(c);
  c << "Transform { rotation 1 0 0 1.570779632679 children [ ";
  c << "Shape { appearance USE App1 geometry Cylinder { height 1.0 radius 1.0 } } ] }\n";
  xform.CloseVRMLCode(c);
}

void Cylinder::PrintMe(std::ostream& c) {
  c << "Object Cylinder\n";
  BaseObject::PrintMe(c);
}

IntersectionSet Cylinder::ObjectVoxelIntersect(Vec3 src, Vec3 det) {
	return IntersectionSet();
}

void Cylinder::ObjectPathIntersect(Vec3 src, Vec3 delta, IntersectionSet& T) {
  int cylinderIntersectFlag, planesIntersectFlag, insidePlanesFlag, insideCylinderFlag;
  double intersectTimes[4];
  char intersectSequence[4];
  double coefA = delta.x*delta.x + delta.y*delta.y;
  double coefB = 2.0*(delta.x*src.x + delta.y*src.y);
  double coefC = src.x*src.x+src.y*src.y-1.0;
  double discrim = coefB*coefB - 4.0*coefA*coefC;

  // Find the intersection times of the ray with the cylinder and the planes.  
  // First, look at the cylinder.

  if (discrim < 0.0)  // The ray doesn't intersect the cylinder.
    {
      insideCylinderFlag = 0;
      cylinderIntersectFlag = 0;
    }
  else if (discrim == 0.0) // The ray could: 1. Just touch the cylinder.
                           //                2. Be parallel to the axis of the cylinder.
    {
      if (coefA != 0.0)  // The ray just touched the surface of the cylinder.
	{
	  intersectTimes[0] = intersectTimes[1] = -coefB/(2.0*coefA);
	  intersectSequence[0] = intersectSequence[1] = 'c';
	  insideCylinderFlag = 0;
	  cylinderIntersectFlag = 1;
	}
      else               //  The ray is parallel to the axis of the cylinder.
	//  See if the ray is inside or outside of the cylinder.
	{
	  if (coefC < 0.0) // Ray is inside the cylinder.
	    {
	      intersectTimes[0] = -(intersectTimes[1] = INFTY);
	      intersectSequence[0] = intersectSequence[1] = 'c';
	      insideCylinderFlag = 1;
	      cylinderIntersectFlag = 0;
	    }
	  else             // Ray is outside of the cylinder
	    {
	      insideCylinderFlag = 0;
	      cylinderIntersectFlag = 0;
	    }
	}
    }
  else     // discrim > 0.0.
    {
      intersectTimes[0] = (-coefB - sqrt(discrim))/(2.0*coefA);
      intersectTimes[1] = (-coefB + sqrt(discrim))/(2.0*coefA);
      intersectSequence[0] = intersectSequence[1] = 'c';
      insideCylinderFlag = 0;
      cylinderIntersectFlag = 1;
    }

  // Consider the planes.

  double lt3 = delta.z;
  double ct3 = src.z;

  coefA = lt3;
  coefB = ct3;

  if (coefA == 0.0)  // The ray is parallel to the planes.
    {
      if ((coefB < (1.0/2.0)) && (coefB > (-1.0/2.0)))  // The ray is between the planes.
	{
	  intersectTimes[2] = -(intersectTimes[3] = INFTY);
	  intersectSequence[2] = intersectSequence[3] = 'p';
	  insidePlanesFlag = 1;
	  planesIntersectFlag = 0;
	}
      else           // The ray is outside of the planes.
	{
	  insidePlanesFlag = 0;
	  planesIntersectFlag = 0;
	}
    }
  else              // The ray intersects the planes.
    {
      intersectTimes[2] = (1.0/2.0 - coefB)/coefA;
      intersectTimes[3] = ((-1.0/2.0) - coefB)/coefA;
      intersectSequence[2] = intersectSequence[3] = 'p';
      insidePlanesFlag = 0;
      planesIntersectFlag = 1;
    }

  // Sort the times.

  if ((insideCylinderFlag || cylinderIntersectFlag) &&
      (insidePlanesFlag || planesIntersectFlag))
    {
      int ctr1, ctr2;
      double tempTime;
      char tempChar;

      for (ctr1=0;ctr1<3;ctr1++)
	for (ctr2=0;ctr2<3;ctr2++)
	  if (intersectTimes[ctr2] > intersectTimes[ctr2+1])
	    {
	      tempTime = intersectTimes[ctr2+1];
	      tempChar = intersectSequence[ctr2+1];
	      intersectTimes[ctr2+1] = intersectTimes[ctr2];
	      intersectSequence[ctr2+1] = intersectSequence[ctr2];
	      intersectTimes[ctr2] = tempTime;
	      intersectSequence[ctr2] = tempChar;
	    }

      // Two sequences are not an appropriate intersection sequence:
      //            ppcc and ccpp.
      // Otherwise, we are good to go.  Take the two middle times.

      if ((strncmp(intersectSequence,"ppcc",4) == 0) || 
	  (strncmp(intersectSequence,"ccpp",4) == 0))
	{
	}
      else                // Copy the intersection times. {
	{
	  T.push_back(Intersection(my_ID,intersectTimes[1],true,density,priority,material));
	  T.push_back(Intersection(my_ID,intersectTimes[2],false,density,priority,material));
	}
    }
  else                  // There was no intersection.
    {
    }
}

double Cylinder::CalculateBoundingRadius() {
  return sqrt(2.0);
}


