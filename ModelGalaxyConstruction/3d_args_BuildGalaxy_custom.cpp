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
const float Zscale = 0.28; //scale height


//Random Walk Step Sizes:
const double Rstep = 0.01; //kpc
const float zStep = 0.01; //kpc


//Mock Galaxy Parameters:
const int nPoints = 300000;     //Number of mock stars to generate
const int nIter = 75000;        //Number of iterations for random walk


//Functions:

//////////////////////////////
//      verticalWave()
// creates model of wave in z
//////////////////////////////
float verticalWave(float zz)
{
  //No wave:
  /*
  float mu = 0.0; //UPDATE!
  float expScale = 0.4647;
  return exp(-abs(zz + mu)/expScale);
  */

  //sech^2
  float shift = 0.02;   //kpc (zSun)
  float scale = 0.28;   //kpc scale height
  float cosh1 = 0.5*(exp((zz+shift)/(2*scale)) + exp(-(zz+shift)/(2*scale)));
  return 1.0/(cosh1*cosh1);

  //sech^2 times wave here:
  //float shift = 0.00;   //kpc (zSun) 
  //float scale = 0.28;   //kpc scale height
  //float cosh1 = 0.5*(exp((zz+shift)/(2*scale)) + exp(-(zz+shift)/(2*scale)));
  //return (1.0 + 0.2*sin(8.0 * zz))/(cosh1*cosh1);
  //return (1.0 + 0.2*cos(8.0 * zz))/(cosh1*cosh1);


  //Gaussian sum wave:
  //return (exp(-abs(zz)/Zscale) + 0.04 + 0.55*(exp(-zz*zz/1.0) - 0.8*exp(-zz*zz/2.0) + 0.2*exp(-zz*zz/4.0)));

  //sine wave:
  //return exp(-abs(zz)/Zscale)*(1.0 + 0.2*sin(8.0 * zz));
}

//////////////////////////////
//      distributionFunction()
// computes distribution function
//////////////////////////////
float distributionFunction(float RR, float zz)
{
  //can multiply below by Norm but it doesn't matter -- we take ratio.
  if((RR < 7.6) or (RR > 8.2) or (zz > 0) or (abs(zz) < 0.2) or (abs(zz) > 2.0)){
    return 0;                    //UPDATE BELOW START POINT!!!
  }
  else{
    float vertTerm = verticalWave(zz);
    return exp(-RR/Rscale) * vertTerm;
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

  //check acceptance:
  if((ratio >= 1)){
    accept = 1;
  }
  else {
    float randomDraw = random_;
    if((randomDraw < ratio)) {
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
    float ztemp = -0.6;
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
    float Phitemp = 180.0 + 4.0 * (2*ran1.rand() - 1); //phi

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
