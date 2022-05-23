// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "CrossSection.hpp"
#include <float.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

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
  stringstream CSDataPath;
  CSDataPath << FileBase << "-all.dat";
  //test if file already exists
  struct stat file_exists;
  if (0==stat(CSDataPath.str().c_str(), &file_exists)) {
      //cout << "Cross section database " << CSDataPath.str() << " already exists, will directly read." << endl;
      ifstream ifs;
      ifs.open(CSDataPath.str().c_str(), ios::binary);
      ifs.seekg(0, ifs.end);
      int db_size = ifs.tellg();
      ifs.seekg(0, ifs.beg);

      char* buffer = new char[db_size];
      ifs.read(buffer, db_size);
      char* this_pt = buffer;

      int Z, num;
      double eng, xsec;
      while(this_pt < buffer+db_size) {
          map<double, double> thisEntry;
          Z = *(int*) this_pt;
          this_pt += 4;
          num = *(int*) this_pt; 
          this_pt += 4;
          for (int i=0; i<num; ++i) {
              eng = *(double*) this_pt;
              this_pt += 8;
              xsec = *(double*) this_pt;
              this_pt += 8;
              thisEntry[eng] = xsec;
          }
          this->TheData[Z] = thisEntry;
      }
      ifs.close();
      delete buffer;
  }
  else {
    cout << "Generating cross section database: " << CSDataPath.str() << endl;
    // cout << "You should only see this notice once." << endl;
    ofstream ofs;
    ofs.open(CSDataPath.str().c_str(), ios::binary);
    ofs.seekp(0, ios::beg);
    //#pragma omp parallel for 
    for (unsigned int Z = 1; Z <= this->MaxZ; Z++){
      map < double, double > ThisEntry;
      stringstream FileName;
      FileName << FileBase << "-" << Z << ".dat";
      ifstream IF(FileName.str().c_str());
      /*      while(IF >> Energy >> CS){
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
      if( ThisEntry.size() == 0) {
        //ofs.close();
        cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl 
         << "failed to load data for " << FileName.str() << "\n"
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl ;
      }
      this->TheData[Z] = ThisEntry;

      //now write this->theData to file
      ofs.write((char*)&Z, sizeof(Z));
      int map_size = ThisEntry.size();
      ofs.write((char*)&map_size, sizeof(map_size));
      if (map_size>0) {
          map<double, double>::iterator it;
          for (it = ThisEntry.begin(); it != ThisEntry.end(); ++it) {
              ofs.write((char*)&(it->first), sizeof(double));
              ofs.write((char*)&(it->second), sizeof(double));
          }
      }
    }
    ofs.close();
  }
  return true;
}

//bool  CrossSection::write(ostream& f){
//  map < double, double >::iterator it2;
//  map <int, < double, double >>::iterator it1;
//  for (it1=this->TheData.begin(); it1!=this->TheData.end(); it1++) {
//    for (it2=it1->second.begin(); it2!=it2->second.end(); it2++) {
//        f.write((char*)&it2, sizeof(it2));
//    }
//  }
//  return true;
//}
