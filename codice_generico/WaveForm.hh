#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"

class WaveForm {
  
private:
  Int_t numOfPulses;
  std::vector<Int_t> gateStart;
  std::vector<Int_t> gateStop;
  std::vector<Float_t> baseline;
  std::vector<Float_t> integral;
  std::vector<Float_t> max;
  Float_t totalIntegral;
  std::vector<Float_t> deltaVec;
  Float_t waveRMS;
  Int_t avgWindPoints;
  Int_t thresUp;
  Int_t thresDown;
  Int_t sampleToTime;
  
public:
  WaveForm(Int_t*,Int_t,Int_t,Int_t,Int_t,Int_t);
  void ComputeDeltaVector(Int_t*,Int_t);
  void DrawDelta(TCanvas*);
  Int_t ComputePulsesAndGates(std::vector<Float_t>);
  Int_t ComputeBaseline(Int_t*,Int_t);
  Float_t ComputePulseIntegral(Int_t*,Int_t);
  Float_t ComputePulseMax(Int_t*,Int_t);
  Float_t FindRMS(Int_t*,Int_t);
  Float_t GetTotalIntegral() {return totalIntegral;};
  Float_t GetPulseIntegral(Int_t i) {return integral.at(i);};
  Float_t GetPulseMax(Int_t i) {return max.at(i);};
  Int_t GetNumberOfPulses() {return numOfPulses;};
  Float_t GetRMS() {return waveRMS;};
  Int_t GetGateStart(Int_t i) {return gateStart.at(i);};
  Int_t GetGateStop(Int_t i) {return gateStop.at(i);};
  Float_t GetBaseline(Int_t i) {return baseline.at(i);};  
  
};

WaveForm::WaveForm(Int_t* array,Int_t size,Int_t nPoints,
                   Int_t thUp,Int_t thDown,Int_t nTime) {

  avgWindPoints=nPoints;
  thresUp=thUp;
  thresDown=thDown;
  sampleToTime=nTime;
  totalIntegral=0;
  
  ComputeDeltaVector(array,size);
  numOfPulses=ComputePulsesAndGates(deltaVec);

  if(numOfPulses==0) {
    baseline.push_back(0);
    integral.push_back(0);
    max.push_back(0);
  }
  else {
    for(Int_t i=1; i <= numOfPulses; i++) {
      baseline.push_back(ComputeBaseline(array,i));
      integral.push_back(ComputePulseIntegral(array,i));
      totalIntegral+=integral.at(i-1);
      max.push_back(ComputePulseMax(array,i));
    }
  }
  waveRMS=FindRMS(array,size);
}

void WaveForm::ComputeDeltaVector(Int_t* array,Int_t size) {

  Float_t fwdWind,backWind;
  Float_t delta;
  
  for (Int_t i = avgWindPoints; i < size-avgWindPoints; i++) {
    fwdWind=0; backWind=0;
    
    for (Int_t k = 1; k < avgWindPoints; k++) {
      fwdWind+=array[i+k];
      backWind+=array[i-k];
    }
    
    delta=fwdWind-backWind;
    deltaVec.push_back(delta);
  }
}

void WaveForm::DrawDelta(TCanvas* canvas) {

  TGraph* deltaGraph = new TGraph(2*avgWindPoints+deltaVec.size());
  TLine* upThLine = new TLine(0,thresUp,(deltaVec.size()+2*avgWindPoints)*sampleToTime,thresUp);
  TLine* downThLine = new TLine(0,thresDown,(deltaVec.size()+2*avgWindPoints)*sampleToTime,thresDown);
  upThLine->SetLineColor(kRed);
  downThLine->SetLineColor(kRed);
  
  for (Int_t i = 0; i <avgWindPoints; i++)
    deltaGraph->SetPoint(i,i*sampleToTime,0);
  for (Int_t i = 0; i < deltaVec.size(); i++)
    deltaGraph->SetPoint(i+avgWindPoints,(i+avgWindPoints)*sampleToTime,deltaVec.at(i));
  for (Int_t i = 0; i <avgWindPoints; i++)
    deltaGraph->SetPoint(i+avgWindPoints+deltaVec.size(),(i+avgWindPoints+deltaVec.size())*sampleToTime,0);
    
  canvas->GetSelectedPad();
  deltaGraph->Draw("AWL");
  if(numOfPulses>0) {
    upThLine->Draw();
    downThLine->Draw();
  }
}

Int_t WaveForm::ComputePulsesAndGates(std::vector<Float_t> deltaVec) {
  
  Bool_t up1=kFALSE;
  Bool_t up2=kFALSE;
  Bool_t low1=kFALSE;
  Bool_t low2=kFALSE;

  Int_t peakCount=0;

  for (Int_t i=0; i < deltaVec.size(); i++) {

    if(deltaVec.at(i)>thresUp &&
       up1==kFALSE && up2==kFALSE &&
       low1==kFALSE && low2==kFALSE) {
      gateStart.push_back(i+avgWindPoints);
      up1=kTRUE;
    }
    else if(deltaVec.at(i)<thresUp &&
            up1==kTRUE && up2==kFALSE &&
            low1==kFALSE && low2==kFALSE) {
      up2=kTRUE;
    }
    else if(deltaVec.at(i)>thresUp &&
            up1==kTRUE && up2==kTRUE &&
            low1==kFALSE && low2==kFALSE) {
      gateStart.at(peakCount)=i+avgWindPoints;
      up2=kFALSE;
    }
    else if(deltaVec.at(i)<thresDown &&
            up1==kTRUE && up2==kTRUE &&
            low1==kFALSE && low2==kFALSE) {
      low1=kTRUE;
    }
    else if(deltaVec.at(i)>thresDown &&
            up1==kTRUE && up2==kTRUE &&
            low1==kTRUE && low2==kFALSE) {
      low2=kTRUE;
    }
    else if(deltaVec.at(i)>thresDown/2 &&
            up1==kTRUE && up2==kTRUE &&
            low1==kTRUE && low2==kTRUE) {
      up1=kFALSE;
      up2=kFALSE;
      low1=kFALSE;
      low2=kFALSE;
      gateStop.push_back(i+avgWindPoints);
      peakCount++;
    }
  }

  return peakCount;
}

Int_t WaveForm::ComputeBaseline(Int_t* array,Int_t peakNumber) {

  Float_t base=0;
  Int_t start,stop;
  
  if(peakNumber==1)
    start=0;
  else
    start=gateStop.at(peakNumber-2);

  stop=gateStart.at(peakNumber-1);
  
  for (Int_t i = start; i < stop; i++) {
    base+=array[i];
  }
  base/=(stop-start);
   
  return base;
 }

Float_t WaveForm::ComputePulseIntegral(Int_t* array, Int_t peakNumber) {

  Float_t sum=0;
 
  for (Int_t i = gateStart.at(peakNumber-1); i < gateStop.at(peakNumber-1); i++) {
    sum+=(array[i]-baseline.at(peakNumber-1));
  }
 
  return sum;
}

Float_t WaveForm::ComputePulseMax(Int_t* array, Int_t peakNumber) {

  Float_t max=0;
 
  for (Int_t i = gateStart.at(peakNumber-1); i < gateStop.at(peakNumber-1); i++) {
    if(array[i]-baseline.at(peakNumber-1)>max)
      max=array[i]-baseline.at(peakNumber-1);
  }
 
  return max;
}

Float_t WaveForm::FindRMS(Int_t* array,Int_t size) {

  Float_t mean=0;
  Float_t RMS;
  Int_t sumsquared=0;
  
  for (Int_t i = 0; i < size; i++)
    mean+=array[i];

  mean=mean/Float_t(size);
  
  for (Int_t i = 0; i < size; i++)
    sumsquared += (array[i]-mean)*(array[i]-mean);
  
  RMS = sqrt((Float_t(1)/Float_t(size))*(Float_t(sumsquared)));

  return RMS;
}
