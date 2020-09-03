// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "CrossSection.hpp"
#include <float.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

//--------------------
CrossSection::CrossSection(){
  this->MaxZ = 100;
}

//--------------------
bool  CrossSection::print(int Z, ostream& os){
  
  bool RetVal;
  if(this->TheData.count(Z)){
    
    map<double, double> ThisCrossSection = this->TheData[Z];
    
    for (map<double, double>::iterator iter = ThisCrossSection.begin(); 
	 iter != ThisCrossSection.end(); iter++){
      os << iter->first << " " << iter->second << endl;
    }
    RetVal = true;
  }
  else
    RetVal = false;

  return RetVal;
  
}
//--------------------
double CrossSection::GetValue(int Z, double Energy){
  
  double RetVal;
  map<double, double> ThisCrossSection = this->TheData[Z];
  map<double, double>::iterator LesserIter = ThisCrossSection.begin();
  map<double, double>::iterator GreaterIter = --(ThisCrossSection.end());
  double LesserStep = FLT_MAX, GreaterStep = FLT_MAX;
  for (map<double, double>::iterator iter = ThisCrossSection.begin(); 
       iter != ThisCrossSection.end(); iter++){
    double ThisEnergy = iter->first;
    if((ThisEnergy <= Energy) && (Energy - ThisEnergy) < LesserStep){
      LesserStep = Energy - ThisEnergy;
      LesserIter = iter;
    }
    if((ThisEnergy >= Energy) && (ThisEnergy - Energy) < GreaterStep){
      GreaterStep = ThisEnergy - Energy;
      GreaterIter = iter;
    }
  }
  // linearly interpolate
  if(LesserIter == GreaterIter)
    RetVal = LesserIter->second;
  else if(Energy > GreaterIter->first)
    RetVal = GreaterIter->second;
  else if(Energy < LesserIter->first)
    RetVal = LesserIter->second;
  else{
    double E2 = GreaterIter->first;
    double sigma2 = GreaterIter->second;
    double E1 = LesserIter->first;
    double sigma1 = LesserIter->second;
    double logSigma = (log(sigma1)*log(E2/Energy) + log(sigma2)*log(Energy/E1))/log(E2/E1);
    RetVal = exp(logSigma);
  }
  return RetVal;
}
//--------------------
bool CrossSection::load(string FileBase){

  double Energy, CS;
  for (unsigned int Z = 1; Z <= this->MaxZ; Z++){
    map < double, double > ThisEntry;
    stringstream FileName;
    FileName << FileBase << "-" << Z << ".dat";
    ifstream IF(FileName.str().c_str());
/*    while(IF >> Energy >> CS){
      if(Energy > 0)
	ThisEntry[Energy] = CS;
    } */
    string curLine;
    while(!IF.eof())
    {
        getline(IF, curLine);
		if (curLine[0] == '%')
		{
			continue;
		}
        std::stringstream ss(curLine, std::stringstream::in);
        ss >> Energy >> CS;
        if(Energy > 0)
        {
            ThisEntry[Energy] = CS;
        }
    }
    if( ThisEntry.size() == 0)
      cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl 
	   << "failed to load data for " << FileName.str() << "\n"
	   << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
    this->TheData[Z] = ThisEntry;
  }
  return true;
}


