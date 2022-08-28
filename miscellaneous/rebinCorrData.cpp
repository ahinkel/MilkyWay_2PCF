//This is a program that will build a simulated stellar distribution
//via the Metropolis algorithm.
//
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
using namespace std;

const int rebinFactor = 2; //3 means 1/3 of the bins. be sure to chose a number that goes evenly into the numBins

//Functions:

//////////////////////////////
//        readFile()
// reads input file and fills vectors
//////////////////////////////
void readFile(string fileName, vector<long long int> &rin, vector<long long int> &phiin, vector<long long int> &zin){
  ifstream fin;
  fin.open(fileName);
  if(!fin.is_open()) {
    cout << "ERROR: file not found" << endl;
  }
  else {
    double currData;
    fin >> currData;
    while (!fin.eof()) {
      rin.push_back(currData);
      fin >> currData;
      phiin.push_back(currData);
      fin >> currData;
      zin.push_back(currData);
      fin >> currData;
    }
  }
  fin.close();
  return;
}



//////////////////////////////
//          main()
//////////////////////////////
int main()
{
  //initialize data to read in
  vector<long long int> RR;
  vector<long long int> DD;
  vector<long long int> DR;


  //read in data from mock galaxy file:
  readFile("9Sep2021_400pc_zmin.txt", RR, DD, DR);
  //readFile("mockgalaxy_seed91_stars125k_steps10k_R8183_z0205_phi179180.txt", Rcoord, Phicoord, Zcoord);

  //cout << Rcoord.size() << endl;

  //get numPoints in existing histogram
  int numPoints = RR.size();

  //create newly rebinned vector
  vector<long long> rebinnedRR;
  rebinnedRR.resize(numPoints/rebinFactor);
  vector<long long> rebinnedDD;
  rebinnedDD.resize(numPoints/rebinFactor);
  vector<long long> rebinnedDR;
  rebinnedDR.resize(numPoints/rebinFactor);

  //loop through stars
  for(int i=0; i<RR.size(); i++){
    int rebinBin = i/rebinFactor; //integer division!
    rebinnedRR.at(rebinBin) += RR.at(i);
    rebinnedDD.at(rebinBin) += DD.at(i);
    rebinnedDR.at(rebinBin) += DR.at(i);
  }



  //store data that satisfies cuts into post cut galaxy file (use convention .gal?)
  ofstream fout;
  fout.open("9Sep2021_400pc_zmin_50bins.txt");
  for(int k=0; k < rebinnedRR.size(); k++){
    fout << rebinnedRR.at(k) << " " << rebinnedDD.at(k) << " " << rebinnedDR.at(k) << endl;
  }
  fout.close();
    

  return 0;
}
