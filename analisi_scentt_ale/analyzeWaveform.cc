#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TMath.h" 
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TFile.h>
#include "TROOT.h"
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TApplication.h>
#include <stdlib.h>
#include <vector>
#include "functions.hh"
#include "WaveForm3.hh"
#include "SiliconDetectorRT.hh"
#include "fitfunc.c"

using namespace std;

void analyzeWaveform(Int_t runNumber, Bool_t debug=0, Int_t maxEntries=-1) {

  const Int_t sibcRawStrips = 768;

  const Int_t nWord730 = 4152;
  const Int_t nSample730 = 512;
  const Int_t nChannels730 = 8;
  const Int_t nModules730 = 2;
  const Int_t baseline1730 = 16384-2400;
  const Int_t sampleToTime730 = 2;

  const Int_t nWord720 = 2104;
  const Int_t nSample720 = 256;
  const Int_t nChannels720 = 8;
  const Int_t nModules720 = 7;
  const Int_t baseline1720 = 4092-150;
  const Int_t sampleToTime720 = 4;
  
  const Int_t nWord730_16 = 8304;
  const Int_t nSample730_16 = 512;
  const Int_t nChannels730_16 = 16;
  const Int_t baseline1730_16 = 16384-2400;
  
  const Int_t headerTot = 0;
  const Int_t headerChan = 5;
  const Int_t trailerChan = 2;

  const Int_t nAvrg730 = 25;
  const Int_t thUp730 = 200;
  const Int_t thDown730 = -100;
  const Int_t nAvrg720 = 13;
  const Int_t thUp720 = 100;
  const Int_t thDown720 = -60;
  
  
  string inFileName="run0"+IntegerToString(runNumber)+".root";
  string pedeFileName;
  
  if(runNumber==60078 || runNumber==60079 || runNumber==60080 ||
     runNumber==60081 || runNumber==60082 || runNumber==60083) {
    pedeFileName="root_digi/run040251_pede";
  }
  else if(runNumber>=60026 && runNumber<=60034) {
    pedeFileName="root_digi/run040230_pede";
  }
  else {
    pedeFileName="root_digi/run011005_pede";
  }


  //SILICON DETECTOR INITIALIZATION//
  SiliconDetectorRT* sibc1 = new SiliconDetectorRT("sibc1","Ivfas_data3",384,0.0242,3,0);
  sibc1->ReadAndSetDeadStrips("strip_files/silicon_strip_status.dat");
  sibc1->ComputeAndSetPedestal(pedeFileName.c_str());

  SiliconDetectorRT* sibc2 = new SiliconDetectorRT("sibc2","Ivfas_data3",384,0.0242,3,384);
  sibc2->ReadAndSetDeadStrips("strip_files/silicon_strip_status.dat");
  sibc2->ComputeAndSetPedestal(pedeFileName.c_str());

  SiliconDetectorRT* sibc3 = new SiliconDetectorRT("sibc3","Ivfas_data4",384,0.0242,3,0);
  sibc3->ReadAndSetDeadStrips("strip_files/silicon_strip_status.dat");
  sibc3->ComputeAndSetPedestal(pedeFileName.c_str());

  SiliconDetectorRT* sibc4 = new SiliconDetectorRT("sibc4","Ivfas_data4",384,0.0242,3,384);
  sibc4->ReadAndSetDeadStrips("strip_files/silicon_strip_status.dat");
  sibc4->ComputeAndSetPedestal(pedeFileName.c_str());

  //DATA-IN DECLARATION//
  TFile* inFile= new TFile(inFileName.c_str());
  TTree* dataTree1 = (TTree*)inFile->Get("h1");
  TTree* dataTree2 = (TTree*)inFile->Get("h2");

  Int_t ivfasData3[sibcRawStrips];
  Int_t ivfasData4[sibcRawStrips];
  UShort_t iData730[nModules730][nWord730];
  UShort_t iData720[nModules720+3][nWord720];
  UShort_t iData730_16[nWord730_16];
  
  dataTree1->SetBranchAddress("Ivfas_data3", ivfasData3);
  dataTree1->SetBranchAddress("Ivfas_data4", ivfasData4);
  dataTree2->SetBranchAddress("Idigi_730", iData730);
  dataTree2->SetBranchAddress("Idigi_720", iData720);
  dataTree2->SetBranchAddress("Idigi_730_16", iData730_16);

  //DATA-OUT DECLARATION//
  string outFileName="out_root_test/run0"+IntegerToString(runNumber)+"_out.root";
  Float_t xPos[4];
  Float_t caloEmPhMax[56];
  Float_t caloEmIntegral[56];
  Float_t caloHadPhMax[18];
  Float_t caloHadIntegral[18];
  Float_t cherenkovPhMax[2];
  Float_t cherenkovIntegral[2];
  Float_t muCatchPhMax[4];
  Float_t muCatchIntegral[4];
  Float_t triggerIntegral,triggerPhMax;
  
  TFile* outRootFile = new TFile(outFileName.c_str(),"RECREATE");
  TTree* outTree = new TTree("outTree","outTree");
  outTree->Branch("xPos",xPos,"xPos[4]/F");
  outTree->Branch("caloEmIntegral",caloEmIntegral,"caloEmIntegral[56]/F");
  outTree->Branch("caloEmPhMax",caloEmPhMax,"caloEmPhMax[56]/F");
  outTree->Branch("caloHadIntegral",caloHadIntegral,"caloHadIntegral[18]/F");
  outTree->Branch("caloHadPhMax",caloHadPhMax,"caloHadPhMax[18]/F");
  outTree->Branch("cherenkovIntegral",cherenkovIntegral,"cherenkovIntegral[2]/F");
  outTree->Branch("cherenkovPhMax",cherenkovPhMax,"cherenkovPhMax[2]/F");
  outTree->Branch("muCatchIntegral",muCatchIntegral,"muCatchIntegral[4]/F");
  outTree->Branch("muCatchPhMax",muCatchPhMax,"muCatchPhMax[4]/F");
  outTree->Branch("triggerIntegral",&triggerIntegral,"triggerIntegral/F");
  outTree->Branch("triggerPhMax",&triggerPhMax,"triggerPhMax/F");

  //DEBUG AND VARIABLE DECLARATION//
  Int_t nEntries = dataTree1->GetEntries();
  if(maxEntries!=-1)
    nEntries=maxEntries;

  vector<WaveForm> waveformVec;
    
  Int_t eventNumber=0;

  vector<TCanvas*> canvasVec;
  TCanvas* cSibcPede1;
  TCanvas* cSibcPede2;
  TCanvas* cSibcPede3;
  TCanvas* cSibcPede4;

  if(debug) {
    for(Int_t j=0;j<nModules720+nModules730+nChannels730_16/8.;j++) {
      string name="cDebug"+IntegerToString(j+1);
      canvasVec.push_back(new TCanvas(name.c_str(),name.c_str(), 0, 0, 1600, 800));
      canvasVec.at(j)->Divide(4,4);
    }
  }

  //LOOP ON ALL ENTRIES//
  for(Int_t i = 0; i < nEntries; i++) {

    if(debug) {
      for(Int_t j=0;j<nModules720+nModules730+nChannels730_16/8.;j++) {
        canvasVec.at(j)->Clear("D");
        canvasVec.at(j)->Update();
      }
    }
    
    waveformVec.clear();
    
    dataTree1->GetEntry(i);
    dataTree2->GetEntry(i);

    eventNumber++;
    if(eventNumber%10000==0)
      cout<<"--> Analyzing event "<<eventNumber<<endl;

    //WAVEFORM ANALYSIS//
    for(Int_t n=0; n<nModules720; n++) {//7x8 channels of EM calorimeter
      for(Int_t m=0; m<nChannels720; m++) {
        
        waveformVec.push_back(WaveForm((UShort_t *)iData720,nSample720,nAvrg720,
                                       nWord720,thUp720,thDown720,sampleToTime720,
                                       m,n,baseline1720,headerTot,headerChan,trailerChan,debug));
                 
        if(debug) {
          Int_t vecIndex =n*nChannels720+m;
          Int_t padIndex = Int_t(m/4.)*4+m+1;
          Int_t canvasIndex = n;
          canvasVec.at(canvasIndex)->cd(padIndex); waveformVec.at(vecIndex).DrawWave(canvasVec.at(canvasIndex));
          canvasVec.at(canvasIndex)->cd(padIndex+4); waveformVec.at(vecIndex).DrawDelta(canvasVec.at(canvasIndex));
          canvasVec.at(canvasIndex)->Update();
          cout<<"NUMBER OF PULSES FOUND CHANNEL "<<vecIndex<<" = "<<waveformVec.at(vecIndex).GetNumberOfPulses()<<endl;
        }
      }
    }

    for(Int_t m=0; m<nChannels730_16; m++) {

      if(m<12) {//12 channels hadronic
        waveformVec.push_back(WaveForm((UShort_t *)iData730_16,nSample730_16,nAvrg730,
                                       nWord730_16,thUp730,thDown730,sampleToTime730,
                                       m,0,baseline1730_16,headerTot,headerChan,trailerChan,debug));
      }
      else {//4 channels muon catcher
        waveformVec.push_back(WaveForm((UShort_t *)iData730_16,nSample730_16,nAvrg730,
                                       nWord730_16,150,-100,sampleToTime730,
                                       m,0,baseline1730_16,headerTot,headerChan,trailerChan,debug));
      }
      
      if(debug) {
        Int_t vecIndex =nModules720*nChannels720+m;
        Int_t padIndex = Int_t(m/4.)*4+m+1-Int_t(m/8.)*16;
        Int_t canvasIndex = nModules720+Int_t(m/8.);
        canvasVec.at(canvasIndex)->cd(padIndex); waveformVec.at(vecIndex).DrawWave(canvasVec.at(canvasIndex));
        canvasVec.at(canvasIndex)->cd(padIndex+4); waveformVec.at(vecIndex).DrawDelta(canvasVec.at(canvasIndex));
        canvasVec.at(canvasIndex)->Update();
        cout<<"NUMBER OF PULSES FOUND CHANNEL "<<vecIndex<<" = "<<waveformVec.at(vecIndex).GetNumberOfPulses()<<endl;
      }
    }
    
    for(Int_t n=0; n<nModules730; n++) {
      for(Int_t m=0; m<nChannels730; m++) {

        if(n==0 && m<2) {//cherenkov detectors
          waveformVec.push_back(WaveForm((UShort_t *)iData730,nSample730,nAvrg730,
                                         nWord730,thUp730,thDown730,sampleToTime730,
                                         m,n,baseline1730,headerTot,headerChan,trailerChan,debug));
        }
        else if(n==0 && m==2) {//trigger scinti
          waveformVec.push_back(WaveForm((UShort_t *)iData730,nSample730,nAvrg730,
                                         nWord730,5000,3000,sampleToTime730,
                                         m,n,baseline1730,headerTot,headerChan,trailerChan,debug));
        }
        else if(n==1 && m<6) {//6 channels hadronic
          waveformVec.push_back(WaveForm((UShort_t *)iData730,nSample730,nAvrg730,
                                         nWord730,thUp730,thDown730,sampleToTime730,
                                         m,n,baseline1730,headerTot,headerChan,trailerChan,debug));
        }
        else {//everything else
          waveformVec.push_back(WaveForm((UShort_t *)iData730,nSample730,nAvrg730,
                                         nWord730,thUp730,thDown730,sampleToTime730,
                                         m,n,baseline1730,headerTot,headerChan,trailerChan,debug,kFALSE));
        }
        
        if(debug) {
          Int_t vecIndex = nModules720*nChannels720+nChannels730_16+n*nChannels730+m;
          Int_t padIndex = Int_t(m/4.)*4+m+1;
          Int_t canvasIndex = nModules720+2+n;
          canvasVec.at(canvasIndex)->cd(padIndex); waveformVec.at(vecIndex).DrawWave(canvasVec.at(canvasIndex));
          canvasVec.at(canvasIndex)->cd(padIndex+4); waveformVec.at(vecIndex).DrawDelta(canvasVec.at(canvasIndex));
          canvasVec.at(canvasIndex)->Update();
          cout<<"NUMBER OF PULSES FOUND CHANNEL "<<vecIndex<<" = "<<waveformVec.at(vecIndex).GetNumberOfPulses()<<endl;
        }
      }
    }

    //SILICON DETECTOR ANALYSIS//
    sibc1->Analyze(ivfasData3, 10., 5.);
    sibc2->Analyze(ivfasData3, 10., 5.);
    sibc3->Analyze(ivfasData4, 10., 5.);
    sibc4->Analyze(ivfasData4, 10., 5.);
    
    if(sibc1->GetNumberOfClusters()==1 && sibc2->GetNumberOfClusters()==1 &&
       sibc3->GetNumberOfClusters()==1 && sibc4->GetNumberOfClusters()==1) {
    
      xPos[0] = sibc1->GetXpos().at(0);
      xPos[1] = sibc2->GetXpos().at(0);
      xPos[2] = sibc3->GetXpos().at(0);
      xPos[3] = sibc4->GetXpos().at(0);

      for(Int_t j=0;j<56;j++) {
        caloEmIntegral[j]=waveformVec.at(j).GetPulseIntegral(0);
        caloEmPhMax[j]=waveformVec.at(j).GetPulseMax(0);
      }

      for(Int_t j=0;j<12;j++) {
        caloHadIntegral[j+6]=waveformVec.at(56+j).GetPulseIntegral(0);
        caloHadPhMax[j+6]=waveformVec.at(56+j).GetPulseMax(0);
      }

      for(Int_t j=0;j<6;j++) {
        caloHadIntegral[j]=waveformVec.at(80+j).GetPulseIntegral(0);
        caloHadPhMax[j]=waveformVec.at(80+j).GetPulseMax(0);
      }
        
      for(Int_t j=0;j<4;j++) {
        muCatchIntegral[j]=waveformVec.at(68+j).GetTotalIntegral();
        muCatchPhMax[j]=waveformVec.at(68+j).GetPulseMax(0);
      }
      for(Int_t j=0;j<2;j++) {
        cherenkovIntegral[j]=waveformVec.at(72+j).GetTotalIntegral();
        cherenkovPhMax[j]=waveformVec.at(72+j).GetPulseMax(0);
      }

      triggerIntegral=waveformVec.at(74).GetTotalIntegral();
      triggerPhMax=waveformVec.at(74).GetPulseMax(0);

      outTree->Fill();
    }
    
    if(debug) 
      getchar();
    
  }//END ENTRIES LOOP

  if(debug) {
    cSibcPede1 = new TCanvas("cSibcPede1", "cSibcPede1", 0, 0, 800, 600);
    sibc1->DrawPedestal(cSibcPede1);
    cSibcPede2 = new TCanvas("cSibcPede2", "cSibcPede2", 800, 0, 800, 600);
    sibc2->DrawPedestal(cSibcPede2);
    cSibcPede3 = new TCanvas("cSibcPede3", "cSibcPede3", 0, 600, 800, 600);
    sibc3->DrawPedestal(cSibcPede3);
    cSibcPede4 = new TCanvas("cSibcPede4", "cSibcPede4", 800, 600, 800, 600);
    sibc4->DrawPedestal(cSibcPede4);
  }

  //OUTPUT FILE WRITE//
  outRootFile->cd();
  outTree->Write();
  outRootFile->Close();

  return;
  
}
  
