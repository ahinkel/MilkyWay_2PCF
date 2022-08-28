//This program computes the 3D 2 point correlation function in phi separation only
//
//Austin Hinkel
//25 - Sept - 2020 
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
using namespace std;

const float conv = 3.141592/180.0;
const int numBins = 200; //number of bins in the histogram
const float binMin = 0.00; //in kpc, lowest length scale to probe
const float binMax = 0.40; //in kpc, largest length scale to probe
const float binSize = (float) (binMax-binMin)/numBins;

/////////////////////////////////////////
//         testModulus()
// Determines how far away two stars are
/////////////////////////////////////////
double testModulus(float q1, float q2){
  double modulus;
  modulus = abs(q1-q2); 
  return modulus;
}

/////////////////////////////////////////
//         readFile()
// read in file for galaxy or mock galaxy
// and appends data to existing vectors
/////////////////////////////////////////
void readFile(string fileName, vector<float> &xin, vector<float> &yin, vector<float> &zin){
  ifstream fin;
  fin.open(fileName);
  if(!fin.is_open()) {
    cout << "ERROR: file not found" << endl;
  }
  else {
    float currData;
    fin >> currData;
    while (!fin.eof()) {
      xin.push_back(currData);
      fin >> currData;
      yin.push_back(currData);
      fin >> currData;
      zin.push_back(currData);
      fin >> currData;
    }
  }
  fin.close();
  return;
}

/////////////////////////////////////////
//            makeHist()
// takes coordinate data (real or mock) 
// and counts distances to add to histogram
/////////////////////////////////////////
void makeHist(vector<unsigned long long int> &histogram, vector<float> Y1){
  int whichBin = numBins;
  for(int i=0; i<Y1.size(); i++){
    //avoid double counting:
    for(int j=0; j<i; j++){
      //distance divided by bin -> floor yields bin number:  
      double distance = testModulus(Y1.at(i), Y1.at(j));
      if(distance > 0.0001){
        whichBin = floor(distance/binSize);
      } else whichBin = 0;
      //above if/else statement avoids issue when distance is super small
      if(whichBin < numBins){
        //histogram[whichBin] += (unsigned long long int) 1;
        histogram.at(whichBin) = histogram.at(whichBin) + (unsigned long long int) 1;
      }
    }
  }
  return;
}

/////////////////////////////////////////
//            makeCrossHist()
// counts distances between RR and DD
// and adds them to histogram vector
/////////////////////////////////////////
void makeCrossHist(vector<unsigned long long int> &crossHist, vector<float> Y1, vector<float> Y2){
  int whichBin = numBins;
  for(int i=0; i < Y1.size(); i++){
    for(int j=0; j < Y2.size(); j++){  
      double distance = testModulus(Y1.at(i), Y2.at(j));
      if(distance > 0.0001){
        whichBin = floor(distance/binSize);
      } else whichBin = 0;
      if(whichBin < numBins){
        crossHist.at(whichBin) = crossHist.at(whichBin) + (unsigned long long int) 1;
      }
    }
  }
  return;
}


//-------------------------------MAIN():-----------------
int main(int argC, char** argv){
  //arg stuff:
  string file1 = argv[1];
  string file2 = argv[2];
  string coordUsed = argv[3]; //accepts x, y, z

  //Real galaxy data:
  vector<float> Xreal;
  vector<float> Yreal; //IN DEGREES!
  vector<float> Zreal;
  readFile(file1, Xreal, Yreal, Zreal);


  //Mock galaxy data:
  vector<float> Xmock;
  vector<float> Ymock;
  vector<float> Zmock;
  readFile(file2, Xmock, Ymock, Zmock);



  //Create histograms:
  vector<unsigned long long int> histRR; //random data histogram
  histRR.resize(numBins, (unsigned long long int) 0);
  vector<unsigned long long int> histDD; //real data histogram
  histDD.resize(numBins, (unsigned long long int) 0);
  vector<unsigned long long int> histDR; //cross-correlation histogram
  histDR.resize(numBins, (unsigned long long int) 0);

  vector<float> qReal;
  vector<float> qMock;
  if(coordUsed == "x") {
    qReal = Xreal;
    qMock = Xmock;
  } else if(coordUsed == "y"){
    qReal = Yreal;
    qMock = Ymock;
  } else if(coordUsed == "z"){
    qReal = Zreal;
    qMock = Zmock;
  }

  //Fill histograms:
  makeHist(histRR, qMock);
  makeHist(histDD, qReal);
  makeCrossHist(histDR, qReal, qMock);

  //cout << "Write data: " << endl;
  //write histogram data:
  for(int n=0; n < histDR.size(); n++){
    cout << histRR[n] << " " << histDD[n] << " " << histDR[n] << endl;
  }

  return 0;
}
