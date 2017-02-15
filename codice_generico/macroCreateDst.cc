#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TMath.h"

#include "functions.hh"

using namespace std;

void macroCreateDst(string inFileNumber) {

  string inFileName = "out_root/out_"+inFileNumber+"_raw.root";
  string outFileName = "out_dst/"+inFileNumber+".dat";
  
  TFile* inFile=new TFile(inFileName.c_str());
    
  TDirectoryFile* iVfas3Dir = (TDirectoryFile*)inFile->Get("Ivfas_data3");
  TDirectoryFile* iVfas4Dir = (TDirectoryFile*)inFile->Get("Ivfas_data4");
 
  TTree* iVfas3Tree = (TTree*)iVfas3Dir->Get("Ivfas_data3_tree");
  TTree* iVfas4Tree = (TTree*)iVfas4Dir->Get("Ivfas_data4_tree");
   
  iVfas3Tree->AddFriend(iVfas4Tree);
   
  Int_t nEv3,nEv4;
  Int_t nClu3,nClu4;
  vector<Double_t> *xPos3=0;
  vector<Double_t> *xPos4=0;
  vector<Int_t> *nStrip3=0;
  vector<Int_t> *nStrip4=0;
 
  iVfas3Tree->SetBranchAddress("Ivfas_data3_nEv",&nEv3);
  iVfas3Tree->SetBranchAddress("Ivfas_data4_nEv",&nEv4);

  iVfas3Tree->SetBranchAddress("Ivfas_data3_nClu",&nClu3);
  iVfas3Tree->SetBranchAddress("Ivfas_data4_nClu",&nClu4);

  iVfas3Tree->SetBranchAddress("Ivfas_data3_xPos",&xPos3);
  iVfas3Tree->SetBranchAddress("Ivfas_data4_xPos",&xPos4);

  iVfas3Tree->SetBranchAddress("Ivfas_data3_nStrip",&nStrip3);
  iVfas3Tree->SetBranchAddress("Ivfas_data4_nStrip",&nStrip4);

  ofstream outFile;
  outFile.open(outFileName.c_str());

  if (outFile.is_open()) {
      
    for(Int_t i=0;i<iVfas3Tree->GetEntries();i++) {
     
      iVfas3Tree->GetEntry(i);
   
      if(nClu3==1 && nClu4==1 &&
         xPos3->at(0)>-1 && xPos4->at(0)>-1) {

        outFile.precision(5);
        
        outFile << scientific << xPos3->at(0) << "\t" << xPos4->at(0) << "\t" << nStrip3->at(0) << "\t"
                << nStrip4->at(0) << "\t" << 0.0 << "\t" << 0.1 <<"\t"<< nEv3 << endl;
      }
    }
  }
  outFile.close();
  
}
