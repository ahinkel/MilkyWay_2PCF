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
//21 - Aug - 2020
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
//#include <time.h>
using namespace std;


//Distribution Function Parameters:
const float Rscale = 2.7;  //Radial scale length in kpc
const float hzThin = 0.20;  //vertical scale height, thin disk, in kpc
const float hzThick = 0.70; //vertical scale height, thick disk, in kpc
const float frac = 0.1;    //unitless, fraction of stars in thick disk?
//const float zSun = 0.015;  //z height of sun above plane is ~15pc. in kpc
//const double Norm = 1.0;    //normalization (# of stars) doesn't matter

//Random Walk Step Sizes:
const double Rstep = 0.40; //kpc
const float zStep = 0.25; //kpc


//Mock Galaxy Parameters:
const int nPoints = 100000;     //Number of mock stars to generate
const int nIter = 10000;        //Number of iterations for random walk


//Functions:

//////////////////////////////
//        sech()
// computes hyperbolic secant
//////////////////////////////
float sech(float input)
{
  return 2.0 / (exp(input) + exp(-input));
}

//////////////////////////////
//      distributionFunction()
// computes distribution function
//////////////////////////////
float distributionFunction(float RR, float zz)
{
  float thinTerm = sech((zz)/(2*hzThin));
  float thickTerm = sech((zz)/(2*hzThick));
  //can multiply below by Norm but it doesn't matter -- we take ratio.
  if((RR < 8.0) or (RR > 8.2) or (zz < 0) or (abs(zz) > 1.985) or abs(zz) < 0.215){
    return 0;
  }
  else{
    return exp(-RR/Rscale) * (thinTerm*thinTerm + frac * thickTerm*thickTerm);
  }
}

//////////////////////////////
//      acceptCriterion()
// decides if randomly walked point is accepted or not
//////////////////////////////
bool acceptCriterion(double trialR, double trialZ, double RR, double zz, float random_)
{
  bool accept = 0;
  float ratio = distributionFunction(trialR, trialZ) / distributionFunction(RR, zz);
  //CUTS here:
  // 7 < R < 9 kpc built into Distribution Function.
  // |z| < 3 kpc built into Distribution Function.
  // 174 < phi < 186 built into phi uniform generation.
  // |z| > 0.2 kpc built into Distribution Function.
  // b cuts not implemented

  //make pre-cuts:
  bool cutsSatisfied = true;
  //if(z < 0){
  //i.e. keep only north
  //  cutsSatisfied = false;
  //}

  //check acceptance:
  if((ratio >= 1) and (cutsSatisfied = true)){
    accept = 1;
  }
  else {
    float randomDraw = random_;
    if((randomDraw < ratio) and (cutsSatisfied = true)) {
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
  vector<float> zCoord;
  vector<float> Rcoord;
  vector<float> Phicoord;

  //loop for random walk
  for(int i=0; i < nPoints; i++){
    //make new star, start at sun's pos
    float ztemp = 0.0;
    float Rtemp = 8.1;  //****WARNING - place start point in R range
    //Randomly walk the new star
    for(int j=0; j < nIter; j++) {
      bool acceptPoint = 0;
      float Rnew = makeTrialPoint(Rtemp, Rstep, ran1.rand());
      float zNew = makeTrialPoint(ztemp, zStep, ran1.rand());
      acceptPoint = acceptCriterion(Rnew, zNew, Rtemp, ztemp, ran1.rand());
      if(acceptPoint == 1){
	ztemp = zNew;
        Rtemp = Rnew;
      }
    }

    //Assume axisymmetry, so phi uniformly distributed
    float Phitemp = 179.5 + 0.5 * (2*ran1.rand() - 1); //phi 179-180

    Rcoord.push_back(Rtemp);
    Phicoord.push_back(Phitemp);
    zCoord.push_back(ztemp);
  }


  //data out - write to file.
  ofstream fout;
  fout.open(fileName);
  //fout.open("MockGalseed45511771_10000strs_3000stps.txt");
  for(int k=0; k < nPoints; k++){
    fout << Rcoord[k] << " " << Phicoord[k] << " " << zCoord[k] << endl;
  }
  fout.close();

  return 0;
}
