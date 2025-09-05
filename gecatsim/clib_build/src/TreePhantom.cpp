// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
#include "TreePhantom.h"

std::ostream& operator<<(std::ostream& c, LinearPhantom& p) {
  for (unsigned i=0;i<p.size();i++)
    p[i]->PrintMe(c);
  return c;
}

static int tabLevel = 0;

void emitTabs(std::ostream& c) {
  for (int i=0;i<tabLevel;i++) 
    c << "   ";
}

std::ostream& operator<<(std::ostream& c, TreePhantom* p) {
  if (p == NULL) return c;
  // Print the head
  emitTabs(c);
  c << p->obj->BoundCenter() << " radius = " << p->obj->BoundRadius() << "\n";
  // Increase the tab level
  tabLevel++;
  // Walk the children
  c << p->down;
  // Decrease the tab level
  tabLevel--;
  // Walk the neighbor
  c << p->right;
  return c;
}


TreePhantom::TreePhantom() {
  right = NULL;
  down = NULL;
  obj = NULL;
}

TreePhantom::~TreePhantom() {
  if (right)
    delete right;
  if (down)
    delete down;
  delete obj;
}

// The object does not intersect the tube if the distance from the centroid to the tube
// is greater than the tube radius + object bounding radius
void TreePathTubeTest(TreePhantom* tree, float radius, float *src, float *dir, float SDD, LinearPhantom &P) {
  //Inline the boundary test...
  bool retval;
  float x1, y1, z1;
  float mindist, ip;
  float brad;
  brad = tree->obj->boundRadius + radius;
  x1 = src[0] - tree->obj->bx;
  y1 = src[1] - tree->obj->by;
  z1 = src[2] - tree->obj->bz;
  mindist = x1*x1 + y1*y1 + z1*z1;
  ip = -(dir[0]*x1+dir[1]*y1+dir[2]*z1);
  mindist -= ip*ip;
  retval = ((mindist <= (brad*brad)) && (ip >= 0) && (ip <= SDD));
  if (retval) {
    P.push_back(tree->obj);
    if (tree->down)
      TreePathTubeTest(tree->down,radius,src,dir,SDD,P);
  }
  if (tree->right)
    TreePathTubeTest(tree->right,radius,src,dir,SDD,P);
}

void TreePathIntersect(TreePhantom* tree, float *src, float *det, float *dir, float SDD, IntersectionSet& T) {
  //Inline the boundary test...
  bool retval;
  float x1, y1, z1;
  float mindist, ip;
  float brad;
  brad = tree->obj->boundRadius;
  x1 = src[0] - tree->obj->bx;
  y1 = src[1] - tree->obj->by;
  z1 = src[2] - tree->obj->bz;
  mindist = x1*x1 + y1*y1 + z1*z1;
  ip = -(dir[0]*x1+dir[1]*y1+dir[2]*z1);
  mindist -= ip*ip;
  retval = ((mindist <= (brad*brad)) && (ip >= 0) && (ip <= SDD));
  if (retval) {
    tree->obj->PathIntersect(src,dir,T);
    if (tree->down)
      TreePathIntersect(tree->down,src,det,dir,SDD,T);
  }
  if (tree->right)
    TreePathIntersect(tree->right,src,det,dir,SDD,T);
}

TreePhantom* TreePhantom::BuildTreePhantomFromLinear(LinearPhantom &p) {
  TreePhantom* completed;
  TreePhantom* working;

  // First, convert the array cList to a linked list
  TreePhantom* hptr, *head, *node;
  working = NULL;
  completed = NULL;
  for (unsigned i=0;i<p.size();i++) {
    node = new TreePhantom;
    node->right = working;
    node->down = NULL;
    node->obj = p[i];
    working = node;
  }
  // Loop until the working queue is empty
  while (working != NULL) {
    // Start with the first object
    head = working;
    // Start searching through until we reach the end or find a container
    hptr = working->right;
    bool foundContainer = false;
    while ((hptr != NULL) && !foundContainer) {
      foundContainer = hptr->obj->Contains(head->obj);
      if (!foundContainer) 
	hptr = hptr->right;
    }
    // Either way, we can advance the working queue
    working = working->right;
    // Did we find a container?  If not, move head to completed list
    if (!foundContainer) {
      head->right = completed;
      completed = head;
    } else {
      // We found a container, so add us to its list of children
      head->right = hptr->down;
      hptr->down = head;
    }
  }
  return completed;  
}


