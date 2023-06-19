// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __CrossSection_hpp__
#define __CrossSection_hpp__

#include <map>
#include <string>
#include <iostream>

using namespace std;

class CrossSection {
  typedef map < int, map < double, double > > CrossSectionData_t;
private:
  CrossSectionData_t TheData;
  string             name;
  unsigned int       MaxZ;

public:
  CrossSection();
  bool  load(string FileBase);
  bool  print(int Z, ostream& os);
  double GetValue(int Z, double Energy);
};

#endif

