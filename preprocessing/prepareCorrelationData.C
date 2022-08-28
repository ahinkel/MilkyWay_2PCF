////////////////////////////////////////////////////////////////////////////
// Austin Hinkel --- Original code authored by: Susan Gardner/Wolfgang Korsch
// Modified c. 25-November-2020
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// RUNNING: (outdated)
// make -f Makefile.prepareCorrelationData
// ./prepareCorrelationData
// (optional) can pipe above command into a file, e.g. add to above command >>data.txt
////////////////////////////////////////////////////////////////////////////////////////

#include "prepareCorrelationData.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TMinuit.h"

//bool debug(0);

int main (int argc,char* argv[])
{
  ///////////////////////////////////////////////////////////
  // Initialization Stage:
  ///////////////////////////////////////////////////////////
  TApplication theApp("App", &argc, argv);
  TFile *f = new TFile("FINAL_GAIA_out.root","READ");
  TNtuple* nt = (TNtuple*)(f -> Get("ntuple"));
  //int nentries = (int)nt->GetEntries();
  //cout << nentries << endl;

  //Gaia data info d:d_err:l:b:pmra:pmra_err:pmdec:pmdec_err:gmag:bpmag:rpmag:bpmrp:bpmg:gmrp:R:phi:prlx:prlx_err

  // Open a file, save the ntuple and close the file
  //TFile in_file("../data/conductivity_experiment.root");
  //TNtuple* my_tuple;in_file.GetObject("cond_data",my_tuple);
  float d,d_err,l,b,pmra,pmra_err,pmdec,pmdec_err,gmag,bpmag,rpmag,bpmrp,bpmg,gmrp,R,phi,prlx,prlx_err; float* row_content;
  //cout << "Potential\tCurrent\tTemperature\tPressure\n"; 
  // for (int irow=0;irow<nt->GetEntries();++irow){
  for (int irow=0; irow<nt->GetEntries(); ++irow){
    nt->GetEntry(irow);
    row_content = nt->GetArgs();
    d = row_content[0]; 
    //d_err = row_content[1];
    l = row_content[2];
    b = row_content[3];
    //pmra = row_content[4];
    //pmra_err = row_content[5];
    //pmdec = row_content[6];
    //pmdec_err = row_content[7];
    gmag = row_content[8];
    //bpmag = row_content[9];
    //rpmag = row_content[10];
    bpmrp = row_content[11];
    //bpmg = row_content[12];
    //gmrp = row_content[13];
    R = row_content[14]; //This is R
    phi = row_content[15]; //This is phi
    prlx = row_content[16];
    //prlx_err = row_content[17];

    //b cuts already implemented via query
    if(gmag < 18 && gmag > 14 && d*sin(b*3.141592/180.0) > 0 && bpmrp < 2.5 && bpmrp > 0.5 && prlx > 0 && R > 7.6 && R < 8.4 && phi > 180 && phi < 184 && abs(d*sin(b*3.141592/180.0))>0.20 && abs(d*sin(b*3.141592/180.0)) < 0.3 && !(l>271&&l<287&&abs(b)<39) && !(l>299&&l<307&&abs(b)>41&&abs(b)<48) && !(l<89&&l>73&&abs(b)<39) && !(l<61&&l>53&&abs(b)<48&&abs(b)>41)){
      cout << R << " " << phi << " " << d*sin(b*3.141592/180.0) << endl;
      /*
      cout << d << " "  
	   << d_err << " "  
	   << l << " "  
	   << b << " "  
	   << pmra << " "  
	   << pmra_err << " "  
	   << pmdec << " "  
	   << pmdec_err << " "
  	   << gmag << " "
  	   << bpmag << " "
  	   << rpmag << " "
  	   << bpmrp << " "
  	   << bpmg << " "
  	   << gmrp << " "
  	   << R << " "
  	   << phi << " "
           << prlx << endl;
      */
    }
  }



  //for(size_t i=0; i<vecR.size() ;i++){ 
  //  cout << vecR[i] << endl;
  // }

  f->Close();
  //theApp.Run();
  return 0;
}

