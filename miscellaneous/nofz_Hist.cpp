//This program computes the n(z) distribution for a wedge of Gaia data
//
//Austin Hinkel
//09 - Feb - 2021 
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


//UPDATE BIN SIZE AND LIMITS!!!!!!
const float conv = 3.141592/180.0;
const int numBins = 100; //number of bins in the n(z) histogram
const float binMin = 1.20; //in kpc, lowest z to probe
const float binMax = 3.00; //in kpc, largest z to probe
const float binSize = (float) (binMax-binMin)/numBins;


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
void makeHist(vector<unsigned long long int> &histogram, vector<float> zz){
  int whichBin = numBins;
  for(int i=0; i<zz.size(); i++){
    //distance divided by bin -> floor yields bin number:  
    if(abs(zz[i]) > 0.0001 and abs(zz[i])-binMin >= 0){
      whichBin = floor((abs(zz[i])-binMin)/binSize);
    } else whichBin = 0;
    //above if/else statement avoids issue when distance is super small
    if(whichBin < numBins and abs(zz[i]) - binMin >= 0){
      //histogram[whichBin] += (unsigned long long int) 1;
      histogram.at(whichBin) = histogram.at(whichBin) + (unsigned long long int) 1;
    }
  }
  return;
}



//-------------------------------MAIN():-----------------
int main(){
  //Real galaxy data:
  vector<float> Rsouth;
  vector<float> Phisouth; //IN DEGREES!
  vector<float> Zsouth;
  readFile("AAS_prodrun/GAIAwedges/GAIAdata_R7684_z1230N_phi176184.txt", Rsouth, Phisouth, Zsouth);
  //cout << "file1 read successfully" << endl;

  vector<float> Rnorth;
  vector<float> Phinorth; //IN DEGREES!
  vector<float> Znorth;
  readFile("AAS_prodrun/GAIAwedges/GAIAdata_R7684_z1230N_phi176184.txt", Rnorth, Phinorth, Znorth);


  //Create histograms:
  vector<unsigned long long int> histSouth; //random data histogram
  histSouth.resize(numBins, (unsigned long long int) 0);
  vector<unsigned long long int> histNorth; //random data histogram
  histNorth.resize(numBins, (unsigned long long int) 0);

  //cout << sizeof(histRR) << endl;
  //cout << sizeof(Rreal) << endl;

  //Fill histograms:
  //cout << "Making histograms..." << endl;
  makeHist(histSouth, Zsouth);
  makeHist(histNorth, Znorth);

  //cout << "Write data: " << endl;
  //write histogram data:
  for(int n=0; n < histNorth.size(); n++){
    cout << histSouth[n] << " " << histNorth[n] << endl;
  }

  return 0;
}
