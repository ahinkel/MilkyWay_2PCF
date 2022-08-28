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
//https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
//http://www.foo.be/docs-free/Numerical_Recipe_In_C/c5-8.pdf
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

//input z limits of chebyshev fit!!!!!!!!!!!!!!!!!!! custom
const float zmin = 0.5;
const float zmax = 0.8;

//don't forget to change DF too

//start stars here:
const float Rstart = 8.0;
const float zStart = 0.7;

//Distribution Function Parameters:
const float Rscale = 2.7;  //Radial scale length in kpc
const float hzThin = 0.20;  //vertical scale height, thin disk, in kpc
const float hzThick = 0.70; //vertical scale height, thick disk, in kpc
const float frac = 0.1;    //unitless, fraction of stars in thick disk?
//const float zSun = 0.015;  //z height of sun above plane is ~15pc. in kpc
//const double Norm = 1.0;    //normalization (# of stars) doesn't matter

//Random Walk Step Sizes:
const double Rstep = 0.40; //kpc
const float zStep = 0.3; //kpc

//Mock Galaxy Parameters:
const int nPoints = 330000;     //Number of mock stars to generate
const int nIter = 10000;        //Number of iterations for random walk


//Functions:

//////////////////////////////          
//        readFile()                    
// reads input file and fills Chebyshev coefficients
//////////////////////////////
void readFile(string fileName, vector<double> &ChebyshevCoefficients){
  ifstream fin;
  fin.open(fileName);
  if(!fin.is_open()) {
    cout << "ERROR: file not found" << endl;
  }
  else {
    double currData;
    fin >> currData;
    while (!fin.eof()) {
      ChebyshevCoefficients.push_back(currData);
      fin >> currData;
    }
  }
  fin.close();
  return;
}


//////////////////////////////
//       chebyshevEval()
//   evaluates n(z) approx. 
//////////////////////////////
float chebyshevEval(vector<double> &coefficients, float zValue){
  float result;
  int m = coefficients.size(); //ex: size 3 means 0,1,2
  float y = (zValue-0.5*(zmax+zmin))/(0.5*(zmax-zmin));
  float d = 0.0;
  float dd = 0.0;
  float sv;
  //cout << "mvalue          : " << m << endl;
  for(int i = m-1; i>=1; i--){
    sv = d;
    d = 2*y*d - dd + coefficients.at(i);
    dd = sv;
    //cout << "ivalue          :  " << i << endl;
    //cout << "coefficient used:  " << coefficients.at(i) << endl;
    //cout << "d is            :  " << d << endl;
  }
  //NOTE!!!!  Watch out for how C_0*T_0 is defined!!!
  result = y*d - dd + coefficients.at(0);
  return result;
}


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
float distributionFunction(float RR, float zz, vector<double> &coefs)
{
  //can multiply below by Norm but it doesn't matter -- we take ratio.
  if((RR < 7.6) or (RR > 8.4) or (zz < 0) or (abs(zz) < 0.5 or abs(zz) > 0.8)){
    return 0;
  }
  else{
    float nofz = chebyshevEval(coefs, zz);
    return exp(-RR/Rscale) * nofz;
  }
}

//////////////////////////////
//      acceptCriterion()
// decides if randomly walked point is accepted or not
//////////////////////////////
bool acceptCriterion(double trialR, double trialZ, double RR, double zz, float random_, vector<double> &coefs)
{
  bool accept = 0;
  float ratio = distributionFunction(trialR, trialZ, coefs) / distributionFunction(RR, zz, coefs);
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
  string fileName = argv[2];      //file to write galaxy data to
  string chebyshevFile = argv[3]; //file with Cheb. coefs. stored in it 0-n
  MTRand ran1 (*SEED);

  //read in chebyshev file:
  vector<double> chebyshevCoefs;
  readFile(chebyshevFile, chebyshevCoefs);

  //unit testing:
  //check coef read in:
  //for(int n=0; n<chebyshevCoefs.size(); n++){
  //  cout << chebyshevCoefs.at(n) << endl;
  //}
  //evaluate at x=0.3 to test
  //float testtt = 1.5; 
  //cout << "The approximation of n at this point is: " 
  //     << chebyshevEval(chebyshevCoefs, testtt) << endl;

  
  //initialize data
  vector<float> zCoord;
  vector<float> Rcoord;
  vector<float> Phicoord;

  //loop for random walk
  for(int i=0; i < nPoints; i++){
    //make new star, start at sun's pos
    float ztemp = zStart;  // *****WARNING!!!!! - place start point in z range
    float Rtemp = Rstart;  // ****WARNING - place start point in R range
    //Randomly walk the new star
    for(int j=0; j < nIter; j++) {
      bool acceptPoint = 0;
      float Rnew = makeTrialPoint(Rtemp, Rstep, ran1.rand());
      float zNew = makeTrialPoint(ztemp, zStep, ran1.rand());
      acceptPoint = acceptCriterion(Rnew, zNew, Rtemp, ztemp, ran1.rand(), chebyshevCoefs);
      if(acceptPoint == 1){
	ztemp = zNew;
        Rtemp = Rnew;
      }
    }

    //Assume axisymmetry, so phi uniformly distributed
    float Phitemp = 180.0 + 2.0 * (2*ran1.rand() - 1); //phi 178-182

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
