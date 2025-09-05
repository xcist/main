// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#include "BaseObject.h"
#include "TreePhantom.h"
#include "MatVec.h"
#include "CrossSection.hpp"
#include "CrossSectionHandler.hpp"
#include "Detector.hpp"
#include <iostream>
#include <string>
#include <string.h>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <functional>
#include <cmath>
#include <fstream>
#include <sstream>
extern "C" {
#include "spline.h"
}

#include "label.h"

#ifdef WIN32
#include <sys/time.h>
#include <sys/timeb.h>
#endif

FILE *phantomin;
extern LinearPhantom parsePhantomDefFile(double scaleFactor, MaterialTable &mtab);
extern void printMaterialTable(FILE *fp);

LinearPhantom phant;
TreePhantom* tree;
MaterialTable mtab;
float* srcpos = NULL;
float* detpos = NULL;
float* mu_table = NULL;
float* spectrum = NULL;
float spectrum_total;
double* spectrumPDF = NULL;
int spectrum_length;
int srccount;
int detcount;
int Ecount;
int macrodetcount;
int clustercount;
int maxclustersize;
int** clustermap;
int usePairProduction=0;

bool crossSectionsInitialized = false;
CrossSection ComptonCrossSection;
CrossSection RayleighCrossSection;
CrossSection PairProductionCrossSection;
CrossSection PhotoElectricCrossSection;
CrossSection ComptonScatterFunction;
CrossSection RayleighFormFactor;
CrossSectionHandler Compton;
CrossSectionHandler Rayleigh;
CrossSectionHandler PhotoE;
CrossSectionHandler PairP;
DiscreteTable Scatter;
DiscreteTable Form;

Phantom *DiscretePhantom;
GenericDetector *ScatterDet;
Vec sourcePDF;
Vec sourcePDFWeighted;
Vec sourceEnergies;
float sourceMaxE;
float source_fanangle;
float source_coneangle;
int numSourceSamples;
Vec3 sourcePosition;
Vec3 sourceDirection;
int sourcePDFPoints;
float *bowtie;
float *bowtie_al_tmp, *bowtie_gr_tmp, *bowtie_cu_tmp, *bowtie_ti_tmp;
float *bowtie_al_att, *bowtie_gr_att, *bowtie_cu_att, *bowtie_ti_att;
int bowtierows;
int nrdetrows, nrdetcols;
float detcolcenter, detrowcenter, detcolsize, detrowsize;
float sdd;
int isDetFlat;		// added for flat detector test

#define max(a,b) ((a) > (b) ? (a) : (b))
#define Report(x) cout << x; cout.flush();

#ifdef WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

extern "C" {

  EXPORT  void InitializeCrossSectionDB(char *dirpath, int pairProdFlag) {
    if ((crossSectionsInitialized) && (usePairProduction == pairProdFlag)) return;
    usePairProduction = pairProdFlag;
    // Report("Loading GEANT4/EDLP97 cross section database...");
    // Report("compton..");
    ComptonCrossSection.load(string(dirpath) + "edlp/comp/ce-cs");
    // Report("rayleigh..");
    RayleighCrossSection.load(string(dirpath) + "edlp/rayl/re-cs");
    // Report("photo..");
    if (usePairProduction)
      {
        PairProductionCrossSection.load(string(dirpath) + "edlp/pair/pp-cs");
        // Report("pair..");
      }    
    PhotoElectricCrossSection.load(string(dirpath) + "edlp/phot/pe-cs");
    // Report("scatter..");
    ComptonScatterFunction.load(string(dirpath) + "edlp/comp/ce-sf");
    // Report("form..");
    RayleighFormFactor.load(string(dirpath) + "edlp/rayl/re-ff");
    // Report("done\r\n");
    crossSectionsInitialized = true;
  }

  EXPORT  void GetCrossSectionMAC(int numZ, int *Zdata, float *massfrac, int numE, float *E, float *data) {
    if (!crossSectionsInitialized) {
      cerr << "Cross Sections DB not initialized!\r\n";
      exit(1);
    }
    for (int i=0;i<numE;i++) {
      float muTotalMAC = 0;
      for (int z=0;z<numZ;z++) {
	float muCompton = ComptonCrossSection.GetValue(Zdata[z],E[i]/1.e3);
	float muPhoto = PhotoElectricCrossSection.GetValue(Zdata[z],E[i]/1.e3);
	float muRayleigh = RayleighCrossSection.GetValue(Zdata[z],E[i]/1.e3);
        float muTotal;
        if (usePairProduction)
          {
            float muPair = PairProductionCrossSection.GetValue(Zdata[z],E[i]/1.e3);
            muTotal = muCompton + muPhoto + muRayleigh + muPair;
          }
        else
          muTotal = muCompton + muPhoto + muRayleigh;

	muTotalMAC += massfrac[z]*muTotal/GetAtomicMass(Zdata[z])*0.6022;
      }
      data[i] = muTotalMAC;
    }
  }

  // return the cross section of four individual physic processes: Rayleigh, Compton, Photoelectric, and pair production 
  EXPORT  void GetCrossSectionByProcessMAC(int numZ, int *Zdata, float *massfrac, int numE, float *E, float *data) {
    if (!crossSectionsInitialized) {
      cerr << "Cross Sections DB not initialized!\r\n";
      exit(1);
    }
    for (int i=0;i<numE;i++) {
      float muRayleighMAC = 0;
      float muComptonMAC = 0;
      float muPhotoMAC = 0;
      float muPairMAC = 0;
      for (int z=0;z<numZ;z++) {
        float muCompton = ComptonCrossSection.GetValue(Zdata[z],E[i]/1.e3);
        float muPhoto = PhotoElectricCrossSection.GetValue(Zdata[z],E[i]/1.e3);
        float muRayleigh = RayleighCrossSection.GetValue(Zdata[z],E[i]/1.e3);
        muRayleighMAC += massfrac[z]*muRayleigh/GetAtomicMass(Zdata[z])*0.6022;
        muComptonMAC += massfrac[z]*muCompton/GetAtomicMass(Zdata[z])*0.6022;
        muPhotoMAC += massfrac[z]*muPhoto/GetAtomicMass(Zdata[z])*0.6022;

        if (usePairProduction) {
          float muPair = PairProductionCrossSection.GetValue(Zdata[z],E[i]/1.e3);
          muPairMAC += massfrac[z]*muPair/GetAtomicMass(Zdata[z])*0.6022;
        }
      }
      data[numE*0+i] = muRayleighMAC;
      data[numE*1+i] = muComptonMAC;
      data[numE*2+i] = muPhotoMAC;
      data[numE*3+i] = muPairMAC;
    }
  }

  EXPORT  void SetSourcePositions(int length, float *data) {
    if (srcpos) delete srcpos;
    srcpos = new float[length*4];
    memcpy(srcpos,data,length*4*sizeof(float));
    srccount = length;
  }
  
  EXPORT  void SetDetectorPositions(int length, float *data) {
    if (detpos) delete detpos;
    detpos = new float[length*6];
    memcpy(detpos,data,length*6*sizeof(float));
    detcount = length;
    // Count the number of clusters
    clustercount = 0;
    macrodetcount = 0;
    int i, j;
    for (i=0;i<length;i++) {
      clustercount = max(clustercount,(int)(data[6*i+5]));
      macrodetcount = max(macrodetcount,(int)(data[6*i+4]));
    }
    macrodetcount++;
    clustercount++;
    // Find the size of the largest cluster
    int *clustersizes = new int[clustercount];
    for (i=0;i<clustercount;i++)
      clustersizes[i] = 0;
    for (j=0;j<length;j++)
      clustersizes[(int)(data[6*j+5])]++;
    maxclustersize = 0;
    for (i=0;i<clustercount;i++)
      maxclustersize = max(maxclustersize,clustersizes[i]);
    // Allocate the array for the cluster map
    clustermap = new int*[clustercount];
    for (i=0;i<clustercount;i++) {
      clustermap[i] = new int[maxclustersize+1];
      memset(clustermap[i],0,(maxclustersize+1)*sizeof(int));
    }
    // build the cluster map
    for (j=0;j<length;j++) {
      int cl = (int) (data[6*j+5]);
      clustermap[cl][++clustermap[cl][0]] = j;
      if ((clustermap[cl][0] > clustersizes[cl]) || (cl >= clustercount))
	printf("Uh oh!\r\n");
    } 
  }

  EXPORT void SetSpectrum(int energies, int length, float *data) {
    if (spectrum) delete spectrum;
    if (spectrumPDF) delete spectrumPDF;
    spectrum = new float[length];
    spectrumPDF = new double[length];		// Changed float to double (bug with precision during the integral)
    memcpy(spectrum,data,length*sizeof(float)); // Necessary to support bigger detectors
    spectrum_total = 0;
    for (int i=0;i<length;i++)
      spectrum_total += spectrum[i];
    spectrumPDF[0] = spectrum[0];
    for (int i=1;i<length;i++)
      spectrumPDF[i] = spectrumPDF[i-1] + spectrum[i];
    spectrum_length = length;
    Ecount = energies;
  }

  EXPORT  void SetMuTable(int length, float *data) {
    if (mu_table) delete mu_table;
    mu_table = new float[length];
    memcpy(mu_table,data,length*sizeof(float));
  }

  EXPORT int GetMaterialCount() {
    return mtab.size();
  }
  
  EXPORT  void GetMaterialName(int index, char* matname) {
    strcpy(matname,mtab[index].c_str());
  }
  
  EXPORT void ParsePhantom(double scalefact, char *filename) {
    FILE *fp = fopen(filename,"r");
    if (!fp) {
      std::cerr << "ERROR: unable to open " << filename << " for parsing...\r\n";
      return;
    }
    phantomin = fp;
    phant = parsePhantomDefFile(scalefact,mtab);
    printf("    Found a total of %d objects\r\n", (int)phant.size());
    tree = TreePhantom::BuildTreePhantomFromLinear(phant);
    
    fclose(fp);
  }

  EXPORT void TranslatePhantom_FORBILD_to_tmp(double scalefact, char *filename_in, char *filename_out) {

    FILE *fp = fopen(filename_in,"r");

    if (!fp) {
      std::cerr << "ERROR: unable to open " << filename_in << " for parsing...\r\n";
      return;
    }
    phantomin = fp;
	int result1 = ftell(fp);
	int result2 = ftell(phantomin);

    phant = parsePhantomDefFile(scalefact,mtab);
    printf("    Found a total of %d objects\r\n", (int)phant.size());
    fclose(fp);

    fp = fopen(filename_out,"w");
    fprintf(fp,"#################### MATERIAL TABLE START ####################\n");
    printMaterialTable(fp);
    fprintf(fp,"#################### MATERIAL TABLE END ######################\n");
    std::vector<BaseObject*>::iterator myIter;
    BaseObject* myObjectPtr;
    
    for(myIter=phant.begin();myIter!=phant.end();myIter++)
      {
	myObjectPtr = *myIter;
	ostringstream stream1;
	const char *ch;
	string string1;
	myObjectPtr->PrintMe(stream1);
	string1 = stream1.str();
	//cout << string1;
	ch=string1.data();
	fprintf(fp,"%s\n",ch);
      }
    fclose(fp);
  }

  EXPORT void SetupDiscretePhantom(int pixels, int slices, float voxelsize, int numZ, int *zlist) {
    if (DiscretePhantom) delete DiscretePhantom ;
    DiscretePhantom = new Phantom;
    DiscretePhantom->Initialize(pixels,slices,voxelsize,numZ,zlist);
  }

  EXPORT void RegisterPhantomDensityMap(int numZ, int pixels, int slices, float *rho) {
    DiscretePhantom->Load(numZ,rho);
  }

  // Compute the centroid and radius of an N x 4 array
  float ComputeCentroidRadius(float *data, int count, float *centroid) {
    int i;
    // Centroid calculation
    centroid[0] = 0; centroid[1] = 0; centroid[2] = 0;
    for (i=0;i<count;i++) {
      centroid[0] += data[4*i];
      centroid[1] += data[4*i+1];
      centroid[2] += data[4*i+2];
    }
    centroid[0] /= count;
    centroid[1] /= count;
    centroid[2] /= count;
    float radius = 0;
    // For each point, calculate the L2 distance to the centroid, record max distance
    for (i=0;i<count;i++) {
      float dist;
      dist = sqrt(pow((double)(data[4*i]-centroid[0]),2.0) + 
		  pow((double)(data[4*i+1]-centroid[1]),2.0) +
		  pow((double)(data[4*i+2]-centroid[2]),2.0));
      radius = max(radius,dist);
    }
    return radius;
  }

  EXPORT int ConstructScatterDetector(int colcount, int rowcount, 
					    float colsize, float rowsize,
					    float SDD, float SID, float colfill,
					    float rowfill, float colcenter,
					    float rowcenter, float gridheight,
					    int coldecimation, int rowdecimation,
					    int maxE, int collimator_type, int detectorFlat) {

    colcenter--;  // convert to zero-based from one-based
    rowcenter--;  // convert to zero-based from one-based

    if (ScatterDet) delete ScatterDet; 

    switch (collimator_type)
      {
      case 1:
	
	if (detectorFlat == 0){
	  ScatterDet = new FocallyAlignedXCollimatedDetector(SDD, SDD,
							    colsize, rowsize,
							    rowcount, colcount,
							    rowcenter, colcenter,
							    gridheight, coldecimation, 
							    rowdecimation, maxE);
	  Report("\n\nDetector is Curved...\n");
	  break;
	}
	else{
	  ScatterDet = new XAlignedZCollimatedDetectorFlat(SDD, SID,
							    colsize, rowsize,
							    rowcount, colcount,
							    rowcenter, colcenter,
							    gridheight, coldecimation, 
							    rowdecimation, maxE);
	Report("\n\nDetector is Flat...\n");  
	break;
	}
      
      case 2:
	
	if (detectorFlat == 0){
	ScatterDet = new FocallyAlignedXandZCollimatedDetector(SDD, SDD,
							       colsize, rowsize,
							       rowcount, colcount,
							       rowcenter, colcenter,
							       gridheight, coldecimation, 
							       rowdecimation, maxE);
	Report("\n\nDetector is Curved...\n");
	break;
	}
	else{
	  ScatterDet = new XAlignedZCollimatedDetectorFlat(SDD, SID,
							    colsize, rowsize,
							    rowcount, colcount,
							    rowcenter, colcenter,
							    gridheight, coldecimation, 
							    rowdecimation, maxE);
	Report("\n\nDetector is Flat...\n");
	break;
	}
      default:
	cerr << "Unknown collimator type specified in cfg.collimator_type!\r\n";
	exit(1);
      }
    
    nrdetcols = colcount;
    nrdetrows = rowcount;
    detrowcenter = rowcenter;
    detcolcenter = colcenter;
    sdd = SDD;
    detcolsize = colsize;
    detrowsize = rowsize;
    isDetFlat = detectorFlat;
    
    return ScatterDet->GetCellCount();
  }
					     
  EXPORT void SetSourceEvec(int npoints, float *Evec) {
    sourceEnergies = VecAllocate(npoints);
    sourceMaxE = 0;
    for (int i=0;i<npoints;i++) {
      sourceEnergies[i] = Evec[i];
      sourceMaxE = max(sourceMaxE,Evec[i]);
    }
    Report("Setting up CatSim cross section tables...");
    Report("compton..");
    Compton.InitializeHandler(DiscretePhantom,ComptonCrossSection,1,sourceMaxE*1.1E3,1E2);
    Report("rayleigh..");
    Rayleigh.InitializeHandler(DiscretePhantom,RayleighCrossSection,1,sourceMaxE*1.1E3,1E2);
    Report("photo..");
    PhotoE.InitializeHandler(DiscretePhantom,PhotoElectricCrossSection,1,sourceMaxE*1.1E3,1E2);
    Report("scatter..");
    float wlMax = sourceMaxE/12.43*1.1e8;
    Scatter.InitializeTable(DiscretePhantom,ComptonScatterFunction,0,wlMax,1e5);
    Report("form..");
    Form.InitializeTable(DiscretePhantom,RayleighFormFactor,0,wlMax,1e5);
    Report("done\r\n");
  }

  Vec LinearVec(int count, int delta) {
    Vec t = VecAllocate(count);
    for (int i=0;i<count;i++)
      t[i] = i*delta;
    return t;
  }

#ifdef WIN32
__declspec(dllexport)
#endif
  void SplineUpsample(float *input_signal, int effCols, int numCols,
		      int colDecimationFactor, int effRows, int numRows,
		      int rowDecimationFactor, float *output_signal) {
    Mat fullRows = MatrixAllocate(numRows,effCols);
    Vec rowBuf = VecAllocate(effRows);
    Vec fullRowBuf = VecAllocate(numRows);
    float *col_vector_old = LinearVec(effCols,colDecimationFactor);
    float *col_vector_new = LinearVec(numCols,1);
    float *row_vector_old = LinearVec(effRows,rowDecimationFactor);
    float *row_vector_new = LinearVec(numRows,1);
    if (rowDecimationFactor != 1) {
      for (int i=0;i<effCols;i++) {
	for (int j=0;j<effRows;j++) 
	  rowBuf[j] = input_signal[i+j*effCols];
	spline_interpolate(row_vector_old,rowBuf,effRows,row_vector_new,fullRowBuf,numRows);
	for (int j=0;j<numRows;j++)
	  fullRows[j][i] = fullRowBuf[j];
      }
    } else {
      for (int i=0;i<effCols;i++) 
	for (int j=0;j<numRows;j++)
	  fullRows[j][i] = input_signal[i+j*effCols];
    }
    if (colDecimationFactor != 1) {
      for (int i=0;i<numRows;i++) {
	spline_interpolate(col_vector_old,fullRows[i],effCols,col_vector_new,output_signal+i*numCols,numCols);
      }
    } else {
      for (int i=0;i<numRows;i++) 
	for (int j=0;j<numCols;j++)
	  output_signal[i*numCols+j] = fullRows[i][j];
    }
    VecFree(row_vector_old);
    VecFree(row_vector_new);
    VecFree(col_vector_old);
    VecFree(col_vector_new);
    MatrixFree(fullRows);
    VecFree(rowBuf);
    VecFree(fullRowBuf);
  }

}


float* readraw(char *fname,int len) {
  FILE *fp;
  float *d;
  fp = fopen(fname,"rb");
  if (!fp) {
    std::cerr << "Unable to open file " << fname << " for reading\r\n";
    exit(1);
  }
  d = new float[len];
  fread(d,sizeof(float),len,fp);
  fclose(fp);
  return d;
}

void writeraw(float *d, int len, char *oname) {
  FILE *fp;
  fp = fopen(oname,"wb");
  if (!fp) {
    std::cerr << "Unable to open file " << oname << " for writing\r\n";
    exit(1);
  }
  fwrite(d,sizeof(float),len,fp);
  fclose(fp);
}

#ifdef STANDALONE2
int main(int argc, char *argv[]) {
  float tmp;
  float src_xform[16];
  float det_xform[16];
  float *output, *src_samples, *det_samples, *bow_samples;
  
  if (argc < 2) {
    std::cerr << "Need phantom filename!\n";
    exit(1);
  }
  ////////////////ParsePhantom(10,argv[1]);
  tmp = 1;
  SetSpectrum(1,&tmp);
  tmp = 0.0205899997057020;
  SetMuTable(1,&tmp);
  src_samples = readraw("srcdat.dat",25*4);
  SetSourcePositions(25,src_samples);
  det_samples = readraw("detdat.dat",1420800*6);
  SetDetectorPositions(1420800,det_samples);
  bow_samples = readraw("bowdat.dat",56832);
  SetBowtieVector(56832,bow_samples);
  memset(src_xform,0,16*sizeof(float));
  memset(det_xform,0,16*sizeof(float));
  src_xform[0] = 1; src_xform[5] = 1; src_xform[10] = 1; src_xform[13] = 541; src_xform[15] = 1;
  det_xform[0] = 1; det_xform[5] = 1; det_xform[10] = 1; det_xform[13] = -408; det_xform[15] = 1;
  output = new float[64*888];
  GetIntensities(src_xform,det_xform,1,1,64*888,output);
  writeraw(output,64*888,"testout.dat");
  return 0;
}
#endif

#ifdef STANDALONE
int main(int argc, char *argv[]) {
  std::cout << "Test program\n";
  if (argc < 2) {
    std::cerr << "Need phantom filename!\n";
    exit(1);
  }
  /////////////////ParsePhantom(10,argv[1]);
  int i, j;
  float srcpos[3];
  float dirvec[3];
  float *detpos;
  detpos = new float[3];  
  double matlbuf[10];
  // Rotation
  // [cos -sin; sin cos]
  double *pdata = new double[888];
  FILE *fp;
  fp = fopen("pout.dat","wb");
  for (i=0;i<984;i++) {
    printf("view %d...\r",i);
    fflush(stdout);
    float beta = i/984.0*2*M_PI;
    srcpos[0] = -540*sin(beta);
    srcpos[1] = 540*cos(beta);
    srcpos[2] = 65;
    for (j=0;j<888;j++) {
      float x, y, z;
      float alpha = (j - 444.0)*1.0/949.0;
      x = 949*sin(alpha);
      y = 540-949*cos(alpha);
      detpos[0] = cos(beta)*x - sin(beta)*y;
      detpos[1] = sin(beta)*x + cos(beta)*y;
      detpos[2] = 65;
      dirvec[0] = detpos[0] - srcpos[0];
      dirvec[1] = detpos[1] - srcpos[1];
      dirvec[2] = detpos[2] - srcpos[2];
      float nrm;
      nrm = sqrt(dirvec[0]*dirvec[0]+dirvec[1]*dirvec[1]+dirvec[2]*dirvec[2]);
      dirvec[0]/=nrm;
      dirvec[1]/=nrm;
      dirvec[2]/=nrm;
      IntersectionSet T;
      TreePathIntersect(tree,srcpos,detpos,dirvec,nrm,T);
      matlbuf[0] = 0;
      T.GetHitList(matlbuf);
      pdata[j] = matlbuf[0];
    }
    fwrite(pdata, sizeof(double), 888, fp);
    memset(pdata,0,sizeof(double)*888);
  }
  fclose(fp);
  return 0;
}
#endif


