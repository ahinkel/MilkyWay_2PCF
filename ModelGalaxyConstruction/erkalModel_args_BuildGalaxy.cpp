//This is a program that will build a simulated stellar distribution
//via the Metropolis algorithm.
//
// to compile and run:
// g++ thisfilename
// ./a.out seedNumber fileName
// seedNumber is for seed used in random number generator
// fileName is name of file to write galaxy to
// file will eventually be used in a bash script, histogrammed, then deleted.
//Austin Hinkel
//24 - Jun - 2021
#include <iostream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <sstream> //aph
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "MersenneTwister.h"
using namespace std;

//consts for functionality
const float conv = 3.141592/180.0;
const float bigG = 1; //units:
const float fourPiG = 4*3.141592*bigG;

//Distribution Function Parameters (see text+table 2 Erkal et al. 2019):
//disk:
const float diskM = 0.68; //10^11 Msols
const float aa = 3.00; //kpc scale length
const float bb = 0.28; //kpc scale height
//halo: (PARAMETERS UPDATED TO MAXIMUM LIKLIHOOD VALUES ON 2021-Jul-14)
const float haloM = 8.2225; //10^11 Msols
const float c = 15.3; //haloConcentration
const float q = 1.272; // oblate/prolate factor
const float h_rs = 15.013; //halo radial scale length (kpc)
const float n_l = 94.759; //direction of n vector
const float n_b = 29.145; //direction of n vector
const float nx = cos(conv * n_l) * cos(conv * n_b); //halo unit vector for orientation, x component
const float ny = sin(conv * n_l) * cos(conv * n_b); //halo unit vector for orientation, y component
const float nz = sin(conv * n_b); //halo unit vector for orientation, z component
const float cStuff = log(1+c) - c/(1+c);
const float qStuff = (1/(q*q) - 1);
//bulge (negligible near sun, ignore for now, but is Hernquist profile):
const float b_rs = 0.500; //bulge scale length in kpc
const float bM = 0; //10^11 Msols

//Random Walk Step Sizes:
const float stepSize = 0.1; //kpc

//Mock Galaxy Parameters:
const int nPoints = 170000;     //Number of mock stars to generate
const int nIter = 20000;        //Number of iterations for random walk


//Functions:

//////////////////////////////
// compute R from x,y
//////////////////////////////
float computeR(float xx, float yy)
{
  return sqrt(xx*xx + yy*yy);
}


//////////////////////////////
// compute phi from x,y
//////////////////////////////
float computePhi(float xx, float yy)
{
  float phi = -999.9;
  if(yy >= 0){
    //atan2 returns a number on 0 to pi, for y>0, just convert to deg.
    phi = atan2(yy, xx)/conv; //y then x for args.
  } else if(yy < 0) {
    //atan2 returns  a number on  -pi to 0 for this region. Add 360.
    phi = 360.0 + atan2(yy, xx)/conv; //y then x for args.
  }
  return phi;
}


//////////////////////////////
//       erkalNFW()
// computes E-NFW halo distribution
//////////////////////////////
float erkalNFW(float xx, float yy, float zz)
{
  float dot = nx*xx + ny*yy + nz*zz;
  float rTilde = sqrt(xx*xx + yy*yy + zz*zz + qStuff*dot*dot);
  float fact1 = (bigG * haloM / cStuff);
  float fact2 = 1/fourPiG;
  float dphidr = fact1 * (log(1+rTilde/h_rs)/(rTilde*rTilde) - 1.0/(rTilde*(rTilde + h_rs)));
  float d2phidr2 = fact1*((3*rTilde+2*h_rs)/(pow(rTilde*(h_rs + rTilde), 2)) - (2*log(1+rTilde/h_rs))/pow(rTilde, 3));
  float drdx = (qStuff*nx*dot + xx)/rTilde;
  float drdy = (qStuff*ny*dot + yy)/rTilde;
  float drdz = (qStuff*nz*dot + zz)/rTilde;
  float d2rdx2 = (1+qStuff*nx*nx)/rTilde - (qStuff*nx*dot + xx)*(qStuff*nx*dot + xx)/(rTilde*rTilde*rTilde);
  float d2rdy2 = (1+qStuff*ny*ny)/rTilde - (qStuff*ny*dot + yy)*(qStuff*ny*dot + yy)/(rTilde*rTilde*rTilde);
  float d2rdz2 = (1+qStuff*nz*nz)/rTilde - (qStuff*nz*dot + zz)*(qStuff*nz*dot + zz)/(rTilde*rTilde*rTilde);
  //cout << fact2 * (dphidr*(d2rdx2 + d2rdy2 + d2rdz2) + (drdx*drdx + drdy*drdy + drdz*drdz)*d2phidr2) << endl;
  return fact2 * (dphidr*(d2rdx2 + d2rdy2 + d2rdz2) + (drdx*drdx + drdy*drdy + drdz*drdz)*d2phidr2);
}


//////////////////////////////
//       MN_disk()        
// computes MN disk distribution
//////////////////////////////
float MN_disk(float xx, float yy, float zz)
{
  float RR = computeR(xx, yy);
  float f1 = bb*bb*diskM/(4*3.141592);
  float tmp = aa + sqrt(zz*zz + bb*bb);
  float num = aa*RR*RR + (aa + 3*sqrt(zz*zz + bb*bb))*tmp*tmp;
  float denom1 = pow((RR*RR + tmp*tmp), 2.5);
  float denom2 = pow((zz*zz + bb*bb), 1.5); 
  //cout << f1 * num / (denom1 * denom2) << endl;
  return f1 * num / (denom1 * denom2);
}


//////////////////////////////
//      distributionFunction()
// computes distribution function
//////////////////////////////
float distributionFunction(float xIn, float yIn, float zIn)
{
  float haloTerm = erkalNFW(xIn, yIn, zIn);
  float diskTerm = MN_disk(xIn, yIn, zIn);
  float tmpR = computeR(xIn, yIn);
  float tmpPhi = computePhi(xIn, yIn);
  if((tmpR < 7.6) or (tmpR > 8.4) or (tmpPhi < 176) or (tmpPhi > 184) or (zIn < 0) or (abs(zIn) > 3.0) or abs(zIn) < 1.2){
    return 0;
  }
  else{
    return haloTerm + diskTerm;
  }
}

//////////////////////////////
//      acceptCriterion()
// decides if randomly walked point is accepted or not
//////////////////////////////
bool acceptCriterion(double trialX, double trialY, double trialZ, double xx, double yy, double zz, float random_)
{
  bool accept = 0;
  float ratio = distributionFunction(trialX, trialY, trialZ) / distributionFunction(xx, yy, zz);
  // b/lmc/smc cuts not implemented

  //check acceptance:
  if(ratio >= 1){
    accept = 1;
  }
  else {
    float randomDraw = random_;
    if(randomDraw < ratio){
      accept = 1;
    }
  }
  return accept;
}

//////////////////////////////
//       makeTrialPoint()
// steps one coordinate into a new, random coordinate.
//////////////////////////////
float makeTrialPoint(float xx, float step, float random_)
{
  //rand gives number b/t 0 and 1
  //cout << random_ << endl;
  return xx + step * (2.0 * random_ - 1); 
}


//////////////////////////////
//          main()
//////////////////////////////
int main(int argC, char** argv)
{

  char* SEED = argv[1];
  string fileName = argv[2];
  MTRand ran1 (*SEED);

  //initialize data
  vector<float> xCoord;
  vector<float> yCoord;
  vector<float> zCoord;
  vector<float> Rcoord;
  vector<float> Phicoord;

  //loop for random walk
  for(int i=0; i < nPoints; i++){
    //make new star, start at
    float ztemp = 1.4;
    float ytemp = 0.0;
    float xtemp = -8.0;  //****WARNING - place start point in range
    //Randomly walk the new star
    for(int j=0; j < nIter; j++) {
      bool acceptPoint = 0;
      float xNew = makeTrialPoint(xtemp, stepSize, ran1.rand());
      float yNew = makeTrialPoint(ytemp, stepSize, ran1.rand());
      float zNew = makeTrialPoint(ztemp, stepSize, ran1.rand());
      acceptPoint = acceptCriterion(xNew, yNew, zNew, xtemp, ytemp, ztemp, ran1.rand());
      if(acceptPoint == 1){
	ztemp = zNew;
        ytemp = yNew;
        xtemp = xNew;
      }
    }
    xCoord.push_back(xtemp);
    yCoord.push_back(ytemp);
    zCoord.push_back(ztemp);
  }

  //CONVERT TO CYL?



  //data out - write to file.
  ofstream fout;
  fout.open(fileName);
  for(int k=0; k < nPoints; k++){
    fout << xCoord[k] << " " << yCoord[k] << " " << zCoord[k] << endl;
  }
  fout.close();

  return 0;
}
