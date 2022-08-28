// reflects azimuth about 180deg line
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

//Functions:

//////////////////////////////
//        readFile()
// reads input file and fills vectors
//////////////////////////////
void readFile(string fileName, vector<double> &rin, vector<double> &phiin, vector<double> &zin){
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
  vector<double> Rcoord;
  vector<double> Phicoord;
  vector<double> Zcoord;

  //read in data from mock galaxy file:
  readFile("AAS_prodrun/GAIAwedges/GAIAdata_R7684_z1230S_phi180184.txt", Rcoord, Phicoord, Zcoord);

  //cout << Rcoord.size() << endl;

  //loop through stars
  for(int i=0; i<Phicoord.size(); i++){
    Phicoord.at(i) = 360.0 - Phicoord.at(i);
  }

  //write Phi reflected data to new file, filename convention start with reflected
  ofstream fout;
  fout.open("AAS_prodrun/GAIAwedges/GAIAdata_R7684_z1230S_phi180184refl176180.txt");
  for(int k=0; k < Rcoord.size(); k++){
    fout << Rcoord.at(k) << " " << Phicoord.at(k) << " " << Zcoord.at(k) << endl;
  }
  fout.close();
    

  return 0;
}
