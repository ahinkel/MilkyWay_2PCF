//Austin Hinkel
//20 - Jan - 2021
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
  readFile("oct7_47.txt", Rcoord, Phicoord, Zcoord);
  //readFile("erkalXYZ_R7684_z0220S_phi178182_seed352.txt", Rcoord, Phicoord, Zcoord);

  //loop through stars
  for(int i=0; i<Zcoord.size(); i++){
    Zcoord.at(i) = 0.0 - Zcoord.at(i); //reflect stars from N to S or S to N.
  }

  //write Z-reflect data to new file, filename convention start with reflected
  ofstream fout;
  fout.open("oct7_47refl.txt");
  //fout.open("erkalXYZ_R7684_z0220SreflN_phi178182_seed352.txt");
  for(int k=0; k < Rcoord.size(); k++){
    fout << Rcoord.at(k) << " " << Phicoord.at(k) << " " << Zcoord.at(k) << endl;
  }
  fout.close();
    

  return 0;
}
