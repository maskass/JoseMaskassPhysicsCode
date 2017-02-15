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
#include "TProfile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH1I.h"

static const Int_t maxRawDataStrips = 1000;//max number of array entries for Ivfas raw data

class SiliconDetectorRT {
  
private:
  std::string detectorName;
  std::string detectorBranchName;
  Int_t numOfStrips;
  Int_t numOfChips;
  Int_t stripOffset;
  std::vector<Double_t> stripStartPositions,stripStopPositions;
  std::vector<Double_t> pedestal,pedestalRmsCM;
  std::vector<Int_t> deadStrips;
  TProfile pedeProf,pedeProfCM,pedeProfRMS,pedeProfRMSCM;
  Double_t pullMaxStrip;
  Int_t numOfClusters;
  std::vector<Double_t> xPos;
  std::vector<Int_t> nStripClu;
    

public:
  SiliconDetectorRT(std::string,std::string,Int_t,Double_t,Int_t,Int_t);
  void ReadAndSetDeadStrips(std::string);
  void ComputeAndSetPedestal(std::string);
  std::vector<Double_t> ComputeCommonMode(Int_t* );
  void DrawPedestal(TCanvas*);
  void SetPedestalProfiles();
  void Analyze(Int_t*,Double_t,Double_t);
  Double_t ComputePullMaxStrip(std::vector<Double_t>);
  std::vector<Int_t> GetStripsOverThreshold(std::vector<Double_t>,Double_t);
  Int_t CountCluAndSetCluRanges(std::vector<Int_t>,std::vector<Int_t>&,std::vector<Int_t>&);
  void ComputeMaxStripsPerCluster(std::vector<Double_t>,Int_t,std::vector<Int_t>,std::vector<Int_t>,
				  std::vector<Int_t>&,std::vector<Int_t>&,std::vector<Int_t>&);
  std::vector<Double_t> ComputeClusterPositions(std::vector<Double_t>,Int_t,Double_t,std::vector<Int_t>,
                                                std::vector<Int_t>,std::vector<Double_t>&, std::vector<Int_t>&);
  std::vector<Double_t> GetXpos() {return xPos;};
  Int_t GetNumberOfClusters() {return numOfClusters;};
  Double_t GetPullMaxStrip() {return pullMaxStrip;};
};

SiliconDetectorRT::SiliconDetectorRT(std::string name,std::string branchName,
                                     Int_t strip, Double_t pitch, Int_t chip, Int_t offset) {
  detectorName=name;
  detectorBranchName=branchName;
  numOfStrips=strip;
  numOfChips=chip;
  stripOffset=offset;
  for(Int_t i=0;i<numOfStrips;i++) {
    stripStartPositions.push_back(i*pitch);
    stripStopPositions.push_back((i+1)*pitch);
  }
  SetPedestalProfiles();
  pullMaxStrip=0;
  numOfClusters=0;
}

void SiliconDetectorRT::ReadAndSetDeadStrips(std::string deadFileName) {
   Int_t lines=FileLineCounter(deadFileName);
   if(lines!=numOfStrips) {
    std::cout<<"WRONG NUMBER OF LINES IN THE DEAD STRIP FILE!!"<<std::endl;
    return;
  }
  else {
    ifstream deadStripFile(deadFileName.c_str());
    Double_t temp1,temp2;
    
    if (deadStripFile.is_open()) {
      for(Int_t i=0;i<numOfStrips;i++) {
	deadStripFile>>temp1>>temp2;
	deadStrips.push_back(Int_t(temp2));
      }
      deadStripFile.close();
    }
    else 
      std::cout<<"CAN'T OPEN DEAD STRIP FILE!!"<<std::endl;
  }
  return;
}

void SiliconDetectorRT::ComputeAndSetPedestal(std::string pedeFileName) {

  pedeFileName+=".root";
  TFile pedeFile(pedeFileName.c_str());
  TTree* inTree = (TTree*)pedeFile.Get("h1");

  Int_t rawData[maxRawDataStrips];
  inTree->SetBranchAddress(detectorBranchName.c_str(),rawData);
  for(Int_t i=0;i<inTree->GetEntries();i++) {
    inTree->GetEntry(i);
    for(Int_t j=0;j<numOfStrips;j++)
      pedeProf.Fill(Double_t(j),Double_t(rawData[j+stripOffset]));
  }
  for(Int_t i=0;i<numOfStrips;i++) 
    pedestal.push_back(pedeProf.GetBinContent(i+1));
  std::vector<Double_t> rawDataCMSub;
  for(Int_t i=0;i<inTree->GetEntries();i++) {
    inTree->GetEntry(i);
    rawDataCMSub=ComputeCommonMode(rawData);
    for(Int_t j=0;j<numOfStrips;j++) 
      pedeProfCM.Fill(Double_t(j),rawDataCMSub.at(j));
  }
  for(Int_t i=0;i<numOfStrips;i++) {
    pedeProfRMS.Fill(Double_t(i),pedeProf.GetBinError(i+1));
    pedeProfRMSCM.Fill(Double_t(i),pedeProfCM.GetBinError(i+1));
    pedestalRmsCM.push_back(pedeProfCM.GetBinError(i+1));
  }
  pedeFile.Close();
}

std::vector<Double_t> SiliconDetectorRT::ComputeCommonMode(Int_t* rawData) {
  
  Double_t sum;
  Int_t j,nstr;
  Int_t stripsPerChip=Int_t(numOfStrips/Double_t(numOfChips));
  std::vector<Double_t> cmChip,rawDataCMSub;
  for(Int_t ichip=0;ichip<numOfChips;ichip++) {
    sum=0.;
    nstr=0;
    for(Int_t istrip=0;istrip<stripsPerChip;istrip++) {
      j=ichip*stripsPerChip+istrip;
      if( (rawData[j+stripOffset]<4000) && (deadStrips.at(j)!=1)) {
	sum+=rawData[j+stripOffset]-pedestal.at(j);
	nstr++;
      }
    }
    if(nstr>0)
      cmChip.push_back(sum/Double_t(nstr));
    else 
      cmChip.push_back(0.);
  }

  for(Int_t ichip=0;ichip<numOfChips;ichip++) {
    for(Int_t istrip=0;istrip<stripsPerChip;istrip++){
      j=ichip*stripsPerChip+istrip;
      if( (rawData[j+stripOffset]<4000) && (deadStrips.at(j)!=1)) {
	rawDataCMSub.push_back(rawData[j+stripOffset]-cmChip.at(ichip));
      }
      else
	rawDataCMSub.push_back(0.);
    }
  }
  return rawDataCMSub;
}

void SiliconDetectorRT::DrawPedestal(TCanvas* canvas) {
  canvas->SetName((detectorName+"_pede").c_str());
  canvas->SetTitle(("Pedestal plots "+detectorName).c_str());
  
  canvas->Divide(1,2);
  canvas->cd(1);
  pedeProf.GetYaxis()->SetRangeUser(1500.,2500.);
  pedeProf.Draw();
  canvas->cd(2);
  pedeProfRMS.GetYaxis()->SetRangeUser(0.,30.);
  pedeProfRMS.Draw();
  pedeProfRMSCM.SetLineColor(kRed);
  pedeProfRMSCM.Draw("SAME");
}

void SiliconDetectorRT::SetPedestalProfiles() {

  pedeProf.BuildOptions(0.,4000.,"s");
  pedeProf.SetBins(numOfStrips,-0.5,numOfStrips-0.5);
  pedeProfCM.BuildOptions(0.,4000.,"s");
  pedeProfCM.SetBins(numOfStrips,-0.5,numOfStrips-0.5);  
  pedeProfRMS.BuildOptions(0.,4000.,"s");
  pedeProfRMS.SetBins(numOfStrips,-0.5,numOfStrips-0.5);
  pedeProfRMSCM.BuildOptions(0.,4000.,"s");
  pedeProfRMSCM.SetBins(numOfStrips,-0.5,numOfStrips-0.5);
}

void SiliconDetectorRT::Analyze(Int_t* rawData,Double_t stripPullCut,
                                Double_t stripPullLatCut) {
  xPos.clear();
  nStripClu.clear();
  
  std::vector<Double_t> rawDataCMSub,rawDataPedeCMSub;
  rawDataCMSub=ComputeCommonMode(rawData);

  for(UInt_t j=0;j<rawDataCMSub.size();j++)
    rawDataPedeCMSub.push_back(rawDataCMSub.at(j)-pedestal.at(j));
    
  pullMaxStrip = ComputePullMaxStrip(rawDataPedeCMSub);
  std::vector<Int_t> stripsOverThreshold;
  stripsOverThreshold=GetStripsOverThreshold(rawDataPedeCMSub,stripPullCut);
    
  std::vector<Int_t> clusterStripStart,clusterStripStop;
  numOfClusters=CountCluAndSetCluRanges(stripsOverThreshold,clusterStripStart,clusterStripStop);
    
  std::vector<Int_t> stripMaxClu, stripMaxCluLeft,stripMaxCluRight;
    
  if(numOfClusters>0) {
    ComputeMaxStripsPerCluster(rawDataPedeCMSub,numOfClusters,clusterStripStart,clusterStripStop,
                               stripMaxClu,stripMaxCluLeft,stripMaxCluRight);
  }
       
  std::vector<Double_t> clusterPH;
  if(numOfClusters>0) {
    xPos=ComputeClusterPositions(rawDataPedeCMSub,numOfClusters,stripPullLatCut,clusterStripStart,
                                 clusterStripStop,clusterPH,nStripClu);
  }
  else {
    xPos.push_back(-99);
    nStripClu.push_back(-99);
  }

}

Double_t SiliconDetectorRT::ComputePullMaxStrip(std::vector<Double_t> rawDataPedeCMSub) {
  
  Int_t iMax=std::distance(rawDataPedeCMSub.begin(),std::max_element(rawDataPedeCMSub.begin(),
								     rawDataPedeCMSub.end()));//.end() non e' incluso
  return rawDataPedeCMSub.at(iMax)/pedestalRmsCM.at(iMax);
}

std::vector<Int_t> SiliconDetectorRT::GetStripsOverThreshold(std::vector<Double_t> rawDataPedeCMSub,
                                                             Double_t stripPullCut) {
  std::vector<Int_t> stripsOverThreshold;
  for(Int_t i=0;i<numOfStrips;i++) {
    if(rawDataPedeCMSub.at(i)>pedestalRmsCM.at(i)*stripPullCut) 
      stripsOverThreshold.push_back(i);
  }
  stripsOverThreshold.push_back(-99);//aggiungo questo in modo che l'ultima componente del vettore sia presente, ma non utile. 
                                     //Mi serve dopo per non andare OutOfRange
  return stripsOverThreshold;
}

Int_t SiliconDetectorRT::CountCluAndSetCluRanges(std::vector<Int_t> stripsOverThreshold,
                                                 std::vector<Int_t>&clusterStripStart,
                                                 std::vector<Int_t>&clusterStripStop) {
  
  Int_t numOfClusters=0;
  Int_t spread;
  Int_t temp=0;
  
  if(stripsOverThreshold.size()>=2) {//una componente buona piu' il -99 del metodo GetStripsOverThreshold
    
    for(UInt_t k=0;k<stripsOverThreshold.size()-1;k++) {//conto i cluster. Metto -1 perche' altrimenti la k mi va OutOfRange sul vettore
      spread=stripsOverThreshold.at(k+1)-stripsOverThreshold.at(k);
      
      if(temp==0) {//setto start e stop alla strip corrente
	clusterStripStart.push_back(stripsOverThreshold.at(k));
	clusterStripStop.push_back(stripsOverThreshold.at(k));
      }
      if(TMath::Abs(spread)>=3) {//se spread è >=3 allora incremento di 1 i clusters
	numOfClusters++;
	temp=0;
      }
      else {//altrimenti setto lo stop alla strip (k+1)esima del vettore e setto temp=1 in modo che non resetti lo start al prox evento
	clusterStripStop.at(numOfClusters)=stripsOverThreshold.at(k+1);
	temp=1;
      }
    }
  }
  else {//se il vettore delle strip sopra soglia ha meno di 2 elementi, nn ho cluster. numOfClusters e' gia' zero
    clusterStripStart.push_back(-99);
    clusterStripStop.push_back(-99);
  }
  
  return numOfClusters;
}

void SiliconDetectorRT::ComputeMaxStripsPerCluster(std::vector<Double_t> rawDataPedeCMSub,Int_t numOfClusters,
                                                   std::vector<Int_t> clusterStripStart,std::vector<Int_t> clusterStripStop,
                                                   std::vector<Int_t>& stripMaxClu,
                                                   std::vector<Int_t>& stripMaxCluLeft,
                                                   std::vector<Int_t>& stripMaxCluRight) {
  
  for(Int_t k=0;k<numOfClusters;k++) {
    
    stripMaxClu.push_back(std::distance(rawDataPedeCMSub.begin(),std::max_element(rawDataPedeCMSub.begin()+clusterStripStart.at(k),
										  rawDataPedeCMSub.begin()+clusterStripStop.at(k)+1)));
    stripMaxCluLeft.push_back(TMath::Max(0,stripMaxClu.at(k)-1));
    stripMaxCluRight.push_back(TMath::Min(numOfStrips-1,stripMaxClu.at(k)+1));
  }
}

std::vector<Double_t> SiliconDetectorRT::ComputeClusterPositions(std::vector<Double_t> rawDataPedeCMSub,Int_t numOfClusters,Double_t stripPullLatCut,
                                                                 std::vector<Int_t> clusterStripStart,std::vector<Int_t> clusterStripStop,
                                                                 std::vector<Double_t>& clusterPH, std::vector<Int_t>& nStripClu) {
  std::vector<Double_t> xPos;
  
  for(Int_t k=0;k<numOfClusters;k++) {
    
    Int_t istart=TMath::Max(0,clusterStripStart.at(k)-2);
    Int_t istop=TMath::Min(numOfStrips-1,clusterStripStop.at(k)+2);
    
    Double_t clusterAdc=0;
    Int_t clusterStrip=0;
    Double_t wSum=0;
    
    for(Int_t ik=istart;ik<=istop;ik++) {//calcolo la posizione del cluster
      
      if(rawDataPedeCMSub.at(ik)>pedestalRmsCM.at(ik)*stripPullLatCut) {
	clusterAdc+=rawDataPedeCMSub.at(ik);
	clusterStrip++;
	wSum+=rawDataPedeCMSub.at(ik)*(stripStartPositions.at(ik)+
				       (stripStopPositions.at(ik)-stripStartPositions.at(ik))/2.);
      }
    }
    if(clusterAdc!=0.) {
      xPos.push_back(wSum/clusterAdc);
      clusterPH.push_back(clusterAdc);
      nStripClu.push_back(clusterStrip);
    }
    else {
      xPos.push_back(-99);
      clusterPH.push_back(0);
      nStripClu.push_back(0);
    }
  }
  
  return xPos;
}
