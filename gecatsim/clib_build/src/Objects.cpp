// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "Objects.hpp"

boundingBox sphereObject::getBoundingBox() {
  boundingBox r;

  r.minX = centerX - radius;
  r.maxX = centerX + radius;
  r.minY = centerY - radius;
  r.maxY = centerY + radius;
  r.minZ = centerZ - radius;
  r.maxZ = centerZ + radius;

  return r;
}


intersectionList sphereObject::getIntersectionList(ray3&) {
	return intersectionList();
}


