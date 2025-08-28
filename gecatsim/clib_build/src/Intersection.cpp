// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#include "Intersection.h"
#include <vector>
#include <algorithm>
#include "MatVec.h"
#include <math.h>

Intersection::Intersection() {
}

Intersection::~Intersection() {
}

Intersection::Intersection(int id, double a_time, bool a_enterSense, 
			   double a_rho, int a_priority, int a_material) {
  time = a_time;
  enterSense = a_enterSense;
  rho = a_rho;
  priority = a_priority;
  material = a_material;
  ID = id;
}

void Intersection::PrintMe(std::ostream& c) const {
  c << "Intersection time " << time << "\r\n";
  c << "Entry time (vs. exit time) " << enterSense << "\r\n";
  c << "priority = " << priority << ", ";
  c << "rho = " << rho << ", ";
  c << "matl = " << material << "\r\n";
  c << "object ID = " << ID << "\r\n";
}

std::ostream& operator<<(std::ostream& c, Intersection& p) {
  p.PrintMe(c);
  return c;
}

std::ostream& operator<<(std::ostream& c, IntersectionSet& p) {
  IntersectionSet::iterator i;
  for (i=p.begin();i!=p.end();i++) {
    i->PrintMe(c);
    c << "\r\n";
  }
  return c;
}

bool Icmp(const Intersection& a, const Intersection& b) {
  return (a.time < b.time) || ((a.time == b.time) && a.enterSense);
}


//typedef std::vector<Intersection> hitStackType;

typedef IntersectionSet hitStackType;

void IntersectionSet::GetHitList(double *ptable) {

  if (empty()) return;
  std::sort(begin(),end(),Icmp);

  hitStackType hitStack;

  //
  // Calculate pathlength
  //
  double previousposition=0.;
  int currentpriority = 0;
  double currentrho = 0;
  int currentmaterial = 0;
  for (IntersectionSet::iterator nr=begin();nr!=end();nr++) {
    //
    // get new data
    //
    double position=nr->time;
    int priority=nr->priority;
    bool sense=nr->enterSense;
    if (hitStack.empty()) {
      currentpriority=0;
      currentrho=0;
      currentmaterial=0;
    } else {
      currentpriority=hitStack.back().priority;
      currentrho=hitStack.back().rho;
      currentmaterial=hitStack.back().material;
      ptable[currentmaterial-1] += (position-previousposition)*currentrho;
    }
    previousposition=position;
    if (priority > currentpriority) {
      // PUSH
      hitStack.push_back(*nr);
    } else if (priority == currentpriority) {
      if (hitStack.empty())
	throw "Error: stack underflow in GetHitList.\r\n";
      hitStack.pop_back();
    }
    else {
      if (sense) {
	hitStackType::iterator cp = hitStack.begin();
	while ((cp != hitStack.end()) && cp->priority < priority)
	  cp++;
	// INSERT IN STACK
	hitStack.insert(cp,*nr);
      } else {
	// GET OUT STACK
	hitStackType::iterator cp = hitStack.begin();
	while ((cp != hitStack.end()) && cp->priority != priority)
	  cp++;
	hitStack.erase(cp);
      }
    }
  }
  if (!hitStack.empty()) {
    printf("ERROR\r\n");
    std::cout << (*this);
    printf("hitstack...\r\n");
    std::cout << hitStack;
    throw "Error: stack not empty on exit of GetHitList.\r\n";
  }
}

void IntersectionSet::RenderIntersections(int N, float dX, float* IDmap, int targID) {
  if (empty()) return;

  std::sort(begin(),end(),Icmp);

  hitStackType hitStack;

  //
  // Calculate pathlength
  //
  double previousposition=0.;
  int currentpriority = 0;
  double currentrho = 0;
  int currentmaterial = 0;
  int currentID = 0;
  for (IntersectionSet::iterator nr=begin();nr!=end();nr++) {
    //
    // get new data
    //
    double position=nr->time;
    int priority=nr->priority;
    bool sense=nr->enterSense;
    if (hitStack.empty()) {
      currentpriority=0;
      currentrho=0;
      currentmaterial=0;
      currentID = 0;
    } else {
      currentpriority=hitStack.back().priority;
      currentrho=hitStack.back().rho;
      currentmaterial=hitStack.back().material;
      currentID = hitStack.back().ID;
      if (currentmaterial == targID) {
	float startndx, stopndx;
	startndx = previousposition/dX;
	stopndx = position/dX;
	int firstwhole, lastwhole, firstfrac, lastfrac;
	firstwhole = (int) ceil(startndx+0.5);
	lastwhole = (int) floor(stopndx-0.5);
	firstfrac = (int) ceil(startndx-0.5);
	lastfrac = (int) floor(stopndx+0.5);
	IDmap[firstfrac] += currentrho*(firstfrac+0.5-startndx);
	IDmap[lastfrac] += currentrho*(stopndx-lastfrac+0.5);
	if (firstfrac == lastfrac)
	  IDmap[lastfrac] -= currentrho*(1.0);
	for (int p=firstwhole;p<=lastwhole;p++)
	  IDmap[p] += currentrho*1.0;
      }
    }
    previousposition=position;
    if (priority > currentpriority) {
      // PUSH
      hitStack.push_back(*nr);
    } else if (priority == currentpriority) {
      if (hitStack.empty())
	throw "Error: stack underflow in RenderIntersections.\n";
      hitStack.pop_back();
    }
    else {
      if (sense) {
	hitStackType::iterator cp = hitStack.begin();
	while ((cp != hitStack.end()) && cp->priority < priority)
	  cp++;
	// INSERT IN STACK
	hitStack.insert(cp,*nr);
      } else {
	// GET OUT STACK
	hitStackType::iterator cp = hitStack.begin();
	while ((cp != hitStack.end()) && cp->priority != priority)
	  cp++;
	hitStack.erase(cp);
      }
    }
  }
}


