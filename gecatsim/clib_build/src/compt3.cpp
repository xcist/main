// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

//This version improves on the straight monte-carlo by doing the following:
//  Launch a photon
//  Track the photon through its life cycle.
//  When it is absorbed, back up and do a deterministic casting...
#include <iostream>
#include <math.h>
#include <vector>

#ifdef WITH_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x) {}
#endif

bool EnableAccel;

// Implements forced detection and MC walking.

// mec2 is the rest mass energy of an electron (.511 KeV)
#define mec2 511
#define maxE 75
#define TMAX(a,b) (((a) < (b)) ? (b) : (a))

float KleinNishinaCrossSectionCos(float costheta, float E) {
  float alpha = E/mec2;
  float c = alpha*(1-costheta);
  return (1/((1+c)*(1+c)))*(1+costheta*costheta + c*c/(1+c));
}

float KleinNishinaCrossSection(float theta, float E) {
  return KleinNishinaCrossSectionCos(cos(theta), E);
}

float KNNormalization[maxE];

void BuildKNNormalizationTable() {
  for (int Energy=0;Energy<maxE;Energy++) {
    float accum = 0;
    for (int k = 0; k < 1000; k++) 
      accum += KleinNishinaCrossSection(2*M_PI/1000.0*k,Energy);
    KNNormalization[Energy] = 1.0/(accum*2*M_PI/1000.0);
  }
}

float GetKleinNishinaProbability(float costheta, int Energy) {
  return KleinNishinaCrossSectionCos(costheta,Energy)*KNNormalization[Energy];
}

int AirPhotons = 0;
int PrimaryPhotons = 0;
int PhotoE = 0;
int InteractionCount = 0;
int comptCount = 0;

float MuAbsorption[maxE] = {4.08E+03,6.16E+02,1.92E+02,8.21E+01,4.20E+01,2.42E+01,1.51E+01,1.01E+01,7.03E+00,5.10E+00,3.82E+00,2.95E+00,2.32E+00,1.87E+00,1.54E+00,1.29E+00,1.09E+00,9.39E-01,8.18E-01,7.21E-01,6.43E-01,5.79E-01,5.26E-01,4.82E-01,4.45E-01,4.14E-01,3.88E-01,3.65E-01,3.45E-01,3.29E-01,3.14E-01,3.01E-01,2.90E-01,2.80E-01,2.71E-01,2.63E-01,2.56E-01,2.50E-01,2.45E-01,2.40E-01,2.35E-01,2.31E-01,2.27E-01,2.23E-01,2.20E-01,2.17E-01,2.15E-01,2.12E-01,2.10E-01,2.08E-01,2.06E-01,2.04E-01,2.02E-01,2.00E-01,1.99E-01,1.97E-01,1.96E-01,1.94E-01,1.93E-01,1.92E-01,1.91E-01,1.90E-01,1.89E-01,1.88E-01,1.87E-01,1.86E-01,1.85E-01,1.84E-01,1.83E-01,1.82E-01,1.82E-01,1.81E-01,1.80E-01,1.79E-01,1.79E-01};

float MuPhotoElectric[maxE] = {4.08E+03,6.16E+02,1.92E+02,8.20E+01,4.19E+01,2.41E+01,1.50E+01,9.92E+00,6.88E+00,4.94E+00,3.66E+00,2.78E+00,2.16E+00,1.71E+00,1.37E+00,1.11E+00,9.17E-01,7.64E-01,6.42E-01,5.44E-01,4.65E-01,4.00E-01,3.46E-01,3.01E-01,2.64E-01,2.32E-01,2.06E-01,1.83E-01,1.63E-01,1.46E-01,1.31E-01,1.18E-01,1.07E-01,9.68E-02,8.80E-02,8.03E-02,7.34E-02,6.72E-02,6.17E-02,5.68E-02,5.24E-02,4.84E-02,4.48E-02,4.15E-02,3.86E-02,3.59E-02,3.34E-02,3.12E-02,2.91E-02,2.72E-02,2.55E-02,2.39E-02,2.25E-02,2.11E-02,1.99E-02,1.87E-02,1.77E-02,1.67E-02,1.58E-02,1.49E-02,1.41E-02,1.34E-02,1.27E-02,1.21E-02,1.15E-02,1.09E-02,1.04E-02,9.87E-03,9.41E-03,8.97E-03,8.56E-03,8.17E-03,7.81E-03,7.46E-03,7.14E-03};

float MuCompton[maxE] = {1.32E-02,4.18E-02,7.07E-02,9.43E-02,1.12E-01,1.26E-01,1.36E-01,1.44E-01,1.50E-01,1.55E-01,1.59E-01,1.62E-01,1.65E-01,1.68E-01,1.70E-01,1.72E-01,1.73E-01,1.75E-01,1.76E-01,1.77E-01,1.78E-01,1.79E-01,1.80E-01,1.81E-01,1.81E-01,1.82E-01,1.82E-01,1.82E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.83E-01,1.82E-01,1.82E-01,1.82E-01,1.82E-01,1.81E-01,1.81E-01,1.81E-01,1.81E-01,1.80E-01,1.80E-01,1.80E-01,1.79E-01,1.79E-01,1.79E-01,1.78E-01,1.78E-01,1.78E-01,1.77E-01,1.77E-01,1.77E-01,1.76E-01,1.76E-01,1.76E-01,1.75E-01,1.75E-01,1.75E-01,1.74E-01,1.74E-01,1.73E-01,1.73E-01,1.73E-01,1.72E-01,1.72E-01,1.72E-01};

// Simple program to compute the scatter off a water cylinder...
float ShimDistance;

// Return a sample from the Klein Nishina cross-section...  
// Input energy is in KeV
float SampleKleinNishina(float Ein) {
  // Convert the input energy from KeV to rest mass energies.
  float E = Ein/511;
  float R;
  bool accept = false;
  while (!accept) {
    // Generate 3 random numbers
    float v1 = drand48();
    float v2 = drand48();
    float v3 = drand48();
    if (v1 > (1+2*E)/(9+2*E)) {
      R = (1+2*E)/(1+2*v2*E);
      float t = (1-R)/E+1;
      accept = !(v3 > 0.5*(t*t+1/R));
    } else {
      R = (1+2*v2*E);
      accept = !(v3 > 4/R*(1-1/R));
    }
  }
  float Eout = Ein/R;
  return Eout;
}

// Given the input and output energies in KeV, calculate the change
// in polar angle given by
//  cos theta = 1-(1/Eout - 1/Ein)*511
float GetScatterPolarAngleCosine(float Ein, float Eout) {
  float costheta = 1-(1/Eout - 1/Ein)*511;
  return costheta;
}

float GetScatterRadialAngle() {
  return (float)(2*M_PI*drand48());
}

float GetMuScatter(float E) {
  return MuCompton[(int) E];
}

float GetMuPhotoElectric(float E) {
  return MuPhotoElectric[(int) E];
}

float GetMuAbsorption(float E) {
  return MuAbsorption[(int) E];
}

float GetMeanFreePathLength(float Ein) {
  return 10.0/MuAbsorption[(int) Ein];
}

class Vec3 {
public:
  float x;
  float y;
  float z;
  void Normalize();
  Vec3() : x(0), y(0), z(0) {};
  Vec3(float a, float b, float c) : x(a), y(b), z(c) {};
  float Dot(Vec3 b) {return (x*b.x+y*b.y+z*b.z);}
  void Print(std::ostream& o);
};

std::ostream& operator<<(std::ostream& o, Vec3& v) {
  v.Print(o);
  return o;
}

void Vec3::Print(std::ostream& o) {
  o << "<" << x << "," << y << "," << z << ">";
}

void Vec3::Normalize() {
  float t = x*x + y*y + z*z;
  t = 1/sqrt(t);
  x *= t;
  y *= t;
  z *= t;  
}

Vec3 phantom;
Vec3 *Detector;
int DetectorCount;
float *DetectorEnergies;
float *MCDetector;

class Photon {
  Vec3 pos;
  Vec3 dir;
  float Energy;
  bool Escaped;
public:
  Photon();
  bool HasEscaped();
  void SetPosition(Vec3 ipos);
  void SetDirection(Vec3 idir);
  Vec3 GetPosition() {return pos;};
  Vec3 GetDirection() {return dir;};
  void SetEnergy(float E);
  float GetEnergy() {return Energy;};
  void Advance(float dist);
  void Intersect(Vec3 ellipse, float &tmin, float &tmax);
  void Step(Vec3 ellipse);
  void Setup(Vec3 ellipse);
  void Interact();
  void DoDeterministicCompton();
  void DoCompton();
  void DoPhotoElectric();
  void RotateDirectionVector(float theta, float phi);
  bool IsAlive();
  void HitDetector();
  void Print(std::ostream& o);
};

std::vector<Photon> Track;

std::ostream& operator<<(std::ostream& o, Photon& v) {
  v.Print(o);
  return o;
}

bool Photon::HasEscaped() {
  return Escaped;
}

void Photon::Print(std::ostream& o) {
  o << "Pos: " << pos << " Dir: " << dir << " Energy " << Energy << "\n";
}

Photon::Photon() {
  Escaped = false;
}

void Photon::SetPosition(Vec3 ipos) {
  pos = ipos;
}
void Photon::SetDirection(Vec3 idir) {
  dir = idir;
}

void Photon::SetEnergy(float E) {
  Energy = E;
}

bool Photon::IsAlive() {
  return (Energy>0);
}

void Photon::DoPhotoElectric() {
  // Kill the photon.
  Energy = 0;
  Escaped = true;
}

bool GetIntersections(Vec3 ellipse, Vec3 pos, Vec3 dir, float &tmin, float &tmax) {
  // Transform pos and dir by the axis scaling parameters
  Vec3 spos;
  Vec3 sdir;
  spos.x = pos.x/ellipse.x;
  spos.y = pos.y/ellipse.y;
  spos.z = pos.z/ellipse.z;
  sdir.x = dir.x/ellipse.x;
  sdir.y = dir.y/ellipse.y;
  sdir.z = dir.z/ellipse.z;
  // Compute the coeffs of the quadratic
  float a, b, c;
  a = sdir.x*sdir.x+sdir.y*sdir.y+sdir.z*sdir.z;
  b = 2*(sdir.x*spos.x+sdir.y*spos.y+sdir.z*spos.z);
  c = spos.x*spos.x+spos.y*spos.y+spos.z*spos.z - 1;
  // Compute the discriminant
  float disc = b*b - 4*a*c;
  bool Escaped = false;
  if (disc >= 0) {
    float sqdisc = sqrt(disc);
    tmin = (-b-sqdisc)/(2*a);
    tmax = (-b+sqdisc)/(2*a);
    if (tmax <= 0)
      Escaped = true;
  } else {
    Escaped = true;
  }
  return Escaped;
}

// Assumes stop is outside the phantom.
float GetTransmissionProbability(Vec3 start, Vec3 stop, int Energy) {
  // Calculate the intersection length.
  Vec3 dir = stop;
  dir.x -= start.x; dir.y -= start.y; dir.z -= start.z;
  dir.Normalize();
  float tmin, tmax;
  if (GetIntersections(phantom,start,dir,tmin,tmax)) 
    return 1;
  tmin = TMAX(0,tmin);
  float pathLength = tmax - tmin;
  return exp(-MuPhotoElectric[Energy]/10.0*pathLength);
}

void Photon::DoDeterministicCompton() {
  comptCount++;
  for (int i=0;i<DetectorCount;i++) {
    // Calculate the direction vector 
    Vec3 todet = Detector[i];
    todet.x -= pos.x; todet.y -= pos.y; todet.z -= pos.z;
    // Get the distance to the detector
    float R = sqrt(todet.Dot(todet));
    todet.Normalize();
    // The dot product of the todet and our current direction is
      // the cosine of the scatter angle
    float costheta = dir.Dot(todet);
    // Get the Klein-Nishina cross section for this scattering angle
    float pKN = GetKleinNishinaProbability(costheta,Energy);
    // Modify this by computing the intersection probability with 
    // the detector (assume a 1 mm detector characteristic length)
    pKN *= 1/(1+2*M_PI*R*sqrt(1-costheta*costheta));
    // Calculate the drop in energy for this scattering angle
    // costheta = 1 - (1/Eout - 1/Ein)*511
    // (1/Eout - 1/Ein)*511 = 1 - costheta
    // 1/Eout = 1/511*(1-costheta) + 1/Ein
    //  Eout = 1/(1/511*(1-costheta) + 1/Ein)
    int EnergyOut = (int)(1.0/(1/511.0*(1.0 - costheta) + 1.0/Energy));
    // Calculate the probability of transmission from us to the detector
    float pTransmit = GetTransmissionProbability(Detector[i],pos,EnergyOut);
    // Accumulate a fractional photon into the detector
    DetectorEnergies[i*maxE + EnergyOut] += pKN*pTransmit;
    DEBUG(std::cout << Energy << " -> " << EnergyOut << " p = " << pKN << "x" << pTransmit << "\n");
  }
}

void Photon::DoCompton() {
  // Sample the Klein-Nishina distribution to get our new energy
  float Enew = SampleKleinNishina(Energy);
  // Get the scatter polar angle
  float costheta = GetScatterPolarAngleCosine(Energy, Enew);
  float phi = GetScatterRadialAngle();
  // Change our direction based on the angles
  RotateDirectionVector(costheta,phi);
  // Adjust our energy
  Energy = Enew;
}

void Photon::RotateDirectionVector(float costheta, float phi) {
  // Now, compute a new direction based on these two angles...  This is
  // somewhat tricky.  For a given vector [a,b,c], there are three
  // vectors that are orthogonal to it:
  // [c,0,-a]
  // [b,-a,0]
  // [0,c,-b]
  // For numerical stability reasons, we use the pair of vectors that
  // includes the maximal value.
  Vec3 bvec;
  Vec3 cvec;
  if ((fabs(dir.x) >= fabs(dir.y)) && (fabs(dir.x) >= fabs(dir.z))) {
    // Use the following as basis vectors:
    // [c,0,-a]
    // [b,-a,0]
    bvec.x = dir.z;
    bvec.y = 0;
    bvec.z = -dir.x;
    cvec.x = dir.y;
    cvec.y = -dir.x;
    cvec.z = 0;
  } else if ((fabs(dir.y) >= fabs(dir.x)) && (fabs(dir.y) >= fabs(dir.z))) {
    // Use the following as basis vectors:
    // [b,-a,0]
    // [0,c,-b]
    bvec.x = dir.y;
    bvec.y = -dir.x;
    bvec.z = 0;
    cvec.x = 0;
    cvec.y = dir.z;
    cvec.z = -dir.y;
  } else {
    // Use the following as basis vectors:
    // [c,0,-a]
    // [0,c,-b]
    bvec.x = dir.z;
    bvec.y = 0;
    bvec.z = -dir.x;
    cvec.x = 0;
    cvec.y = dir.z;
    cvec.z = -dir.y;
  }
  // Normalize the two vectors
  bvec.Normalize();
  cvec.Normalize();
  // Compute a rotated vector in the plane spanned by bvec and cvec
  float cosang = cos(phi);
  float sinang = sin(phi);
  Vec3 dvec;
  dvec.x = cosang*bvec.x + sinang*cvec.x;
  dvec.y = cosang*bvec.y + sinang*cvec.y;
  dvec.z = cosang*bvec.z + sinang*cvec.z;
  // Add a scaled version of the current vector (the scaling should be 
  // cosine(theta))
  dvec.x += dir.x*costheta;
  dvec.y += dir.y*costheta;
  dvec.z += dir.z*costheta;
  // Normalize dvec
  dvec.Normalize();
  // This is the new direction
  dir.x = dvec.x;
  dir.y = dvec.y;
  dir.z = dvec.z;
}



void Photon::Interact() {
  // An interaction has taken place... We are only tracking Incoherent scattering
  // and photoelectric absorption.  So we can compute the relative probabilities
  // as p(scatter) = mu_s/(mu_s+mu_p)
  // and p(photo) = 1-p(scatter).
  float mu_s = GetMuScatter(Energy);
  float mu_p = GetMuPhotoElectric(Energy);
  if (drand48() < (mu_s/(mu_s+mu_p)))
    DoCompton();
  else
    DoPhotoElectric();
}


void Photon::Advance(float dist) {
  pos.x += dir.x*dist;
  pos.y += dir.y*dist;
  pos.z += dir.z*dist;
}


// The object is an elliptical phantom that is axis aligned.
// This routines calculates the intersection times (distance along
// the given ray) to the boundaries of the ellipse.  Returns false
// if the ray never intersects the ellipse.
//  ||Dx||^2 < 1
//  ||D(x0+t*d)||^2 = ||D*x0||^2+t^2*||D*d||^2+t<D*x0,D*d>
void Photon::Intersect(Vec3 e, float &tmin, float &tmax) {
  // Transform pos and dir by the axis scaling parameters
  Vec3 spos;
  Vec3 sdir;
  spos.x = pos.x*e.x;
  spos.y = pos.y*e.y;
  spos.z = pos.z*e.z;
  sdir.x = dir.x*e.x;
  sdir.y = dir.y*e.y;
  sdir.z = dir.z*e.z;
  // Compute the coeffs of the quadratic
  float a, b, c;
  a = sdir.x*sdir.x+sdir.y*sdir.y+sdir.z*sdir.z;
  b = 2*(sdir.x*spos.x+sdir.y*spos.y+sdir.z*spos.z);
  c = spos.x*spos.x+spos.y*spos.y+spos.z*spos.z - 1;
  // Compute the discriminant
  float disc = b*b - 4*a*c;
  if (disc >= 0) {
    float sqdisc = sqrt(disc);
    tmin = (-b-sqdisc)/(2*a);
    tmax = (-b+sqdisc)/(2*a);
    if (tmax <= 0)
      Escaped = true;
  } else {
    Escaped = true;
  }
}

// Step the photon - it has an energy, a position and direction
// Here's how the routine works....  We start by 
// calculating the distance in material as seen by the current photon.
void Photon::Step(Vec3 ellipse) {
  // Calculate the intersection times
  float tmin, tmax;
  Intersect(ellipse,tmin,tmax);
  if (Escaped) return;
  float meanFreePathCompton = 10.0/GetMuScatter(Energy);
  float meanFreePathPhoto = 10.0/GetMuPhotoElectric(Energy);
  float actualPathCompton = -meanFreePathCompton*log(1.0 - drand48());
  float actualPathPhoto = -meanFreePathPhoto*log(1.0 - drand48());
  float actualPath;
  bool isCompton;
  if (actualPathCompton < actualPathPhoto) {
    isCompton = true;
    actualPath = actualPathCompton;
  } else {
    isCompton = false;
    actualPath = actualPathPhoto;
  }
  Track.push_back(*this);
  // If the sample says to leave the 
  if (actualPath > tmax) {
    Advance(tmax+ShimDistance); 
    Escaped = true;
    DEBUG(std::cout << "Escaped at distance " << tmax << "\n");
  } else {
    Advance(actualPath);
    InteractionCount++;
    if (isCompton) {
      DoCompton();
      DEBUG(std::cout << "Compton at distance " << actualPath << " energy is now " << Energy << "\n");
    } else {
      DoPhotoElectric();
      DEBUG(std::cout << "PhotoE at distance " << actualPath << "\n");
    }
  }
}

void SetupTables() {
  ShimDistance = 1;
  BuildKNNormalizationTable();
}

void SetupGeometry() {
  phantom.x = 1.0/(350/2);
  phantom.y = 1.0/(350/2);
  phantom.z = 1.0/(10000);
}

// Advances a photon to the surface of the ellipse - returns false
// if the photon never intersects the surface..
void Photon::Setup(Vec3 e) {
  float tmin, tmax;
  Intersect(e,tmin, tmax);
  if (!Escaped)
  Advance(tmin+ShimDistance);
}

// Generate a photon.
Photon GeneratePhoton() {
  Photon p;
  // Set the energy - we use 80 kEV for now
  p.SetEnergy(74);
  // Set the position - emit from the focal position
#if 0
  p.SetPosition(Vec3(0,540,0));
  // Set the direction - randomly chosen in the fan
  // Get a random half-radian fan angle in radians
  float fanang = (2*drand48()-1)*asin(500/2.0/540.0);
  p.SetDirection(Vec3(sin(fanang),-cos(fanang),0));
#endif
  p.SetPosition(Vec3(-540,0,0));
  float fanang = (2*drand48()-1)*asin(500/2.0/540.0);
  p.SetDirection(Vec3(cos(fanang),-sin(fanang),0));
  return p;
}

void Photon::HitDetector() {
  if (Energy == 74) return;
  if (!Energy) return;
  float t = (-410 - pos.y)/dir.y;
  if (t>0) {
    float x = pos.x+t*dir.x;
    float z = pos.z+t*dir.z;
    if ((z>-10) && (z<10))
      if ((x >= -444) && (x < 443))
	MCDetector[(int)((x+444)*maxE + Energy)]++;
  }
}

void RunEvent() {
  // Start by generating a photon
  Photon p = GeneratePhoton();
  // Step the photon into the object
  p.Setup(phantom);
  if (p.HasEscaped()) {
    p.HitDetector();
    AirPhotons++;
    return;
  }
  DEBUG(std::cout << "start photon \n");
  while (!p.HasEscaped())
    p.Step(phantom);
  if (p.GetEnergy() == 74)
    PrimaryPhotons++;
  else if (p.GetEnergy() == 0)
    PhotoE++;
  // Photon escaped... Calculate intersection with detector.
  p.HitDetector();
  // Examine the track if it has at least two data points then a Compton scattering
  // occured.
  if (EnableAccel && (Track.size() >= 2)) {
    Photon q;
    // Position of the last event is the last scatter location
    q.SetPosition(Track[Track.size()].GetPosition());
    // Energy and direction of the next to last scatter location gives incoming
    q.SetEnergy(Track[Track.size()-1].GetEnergy());
    q.SetDirection(Track[Track.size()-1].GetDirection());
    // Do a deterministic compton...
    q.DoDeterministicCompton();
  }
  // it scattered 
  // Print the track
  //   std::cout << "Track\n";
  //   for (int i=0;i<Track.size();i++) {
  //     std::cout << Track[i] << "\n";
  //   }
  Track.clear();
}

void SetupDetectors() {
  DetectorCount = 888;
  Detector = new Vec3[888];
  for (int i=0;i<888;i++) {
    Detector[i].x = i - 444.0;
    Detector[i].y = -410.0;
    Detector[i].z = 0.0;
  }
  DetectorEnergies = new float[888*maxE];
  MCDetector = new float[888*maxE];
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Requires a count of events and MC/DET flag\n";
    exit(1);
  }
  EnableAccel = atoi(argv[2]);
  SetupGeometry();
  SetupDetectors();
  SetupTables();
  int ecount = atoi(argv[1]);
  for (int i=0;i<ecount;i++)
    RunEvent();
  std::cout << "Air Photons = " << AirPhotons << "\n";
  std::cout << "Primary Photons = " << PrimaryPhotons << "\n";
  std::cout << "PhotoElectric = " << PhotoE << "\n";
  std::cout << "Number of Interactions = " << InteractionCount << "\n";
  FILE *fp;
  if (EnableAccel) {
    for (int i=0;i<maxE*888;i++) {
      DetectorEnergies[i] /= comptCount;
    }
    fp = fopen("outdata.dat","wb");
    fwrite(DetectorEnergies,sizeof(float),maxE*888,fp);
    fclose(fp);
  }
  fp = fopen("mcdata.dat","wb");
  fwrite(MCDetector,sizeof(float),maxE*888,fp);
  fclose(fp);

  return 0;
}


