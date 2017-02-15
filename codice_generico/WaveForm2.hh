#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"

class WaveForm {
  
private:
  Int_t numOfPulses;
  Int_t nSample;
  Float_t totalIntegral;
  Float_t waveRMS;
  Int_t avgWindPoints;
  Int_t thresUp;
  Int_t thresDown;
  Int_t sampleToTime;
  TGraph waveGraph;
  TGraph deltaGraph;
  TLine upThLine;
  TLine downThLine;
  
  std::vector<Int_t> waveVector;
  std::vector<Int_t> gateStart;
  std::vector<Int_t> gateStop;
  std::vector<Float_t> baseline;
  std::vector<Float_t> integral;
  std::vector<Float_t> max;
  std::vector<Float_t> deltaVec;
  std::vector<TLine*> lineStart;
  std::vector<TLine*> lineStop;
  
 public:
  WaveForm(Int_t*,Int_t,Int_t,Int_t,Int_t,Int_t,
           Int_t,Int_t,Int_t,Int_t,Int_t,Bool_t);
  ~WaveForm();
  void ComputeDeltaVector(std::vector<Int_t>);
  Int_t ComputePulsesAndGates(std::vector<Float_t>);
  Int_t ComputeBaseline(std::vector<Int_t>,Int_t);
  Float_t ComputePulseIntegral(std::vector<Int_t>,Int_t);
  Float_t ComputePulseMax(std::vector<Int_t>,Int_t);
  Float_t FindRMS(std::vector<Int_t>);
  void DrawDelta(TCanvas*);
  void DrawWave(TCanvas*);
  void CreateWaveGraph();
  void CreateDeltaGraph();
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
                   Int_t thUp,Int_t thDown,Int_t nTime,
                   Int_t channel,Int_t digiBaseline,Int_t headerTot,
                   Int_t headerChan,Int_t trailerChan, Bool_t debug) {

  nSample=size;
  avgWindPoints=nPoints;
  thresUp=thUp;
  thresDown=thDown;
  sampleToTime=nTime;
  totalIntegral=0;

  //create waveVector as subset of general digitizer array using header/trailer information
  for(Int_t i=0; i < nSample; i++) {
    Int_t index=headerTot+(channel+1)*headerChan + channel*(nSample+trailerChan) + i;
    waveVector.push_back(digiBaseline-array[index]);
  }

  ComputeDeltaVector(waveVector);
  numOfPulses=ComputePulsesAndGates(deltaVec);
  
  if(numOfPulses==0) {
    baseline.push_back(0);
    integral.push_back(0);
    max.push_back(0);
  }
  else {
    for(Int_t i=1; i <= numOfPulses; i++) {
      baseline.push_back(ComputeBaseline(waveVector,i));
      integral.push_back(ComputePulseIntegral(waveVector,i));
      totalIntegral+=integral.at(i-1);
      max.push_back(ComputePulseMax(waveVector,i));
    }
  }
  waveRMS=FindRMS(waveVector);

  if(debug) {
    CreateWaveGraph();
    CreateDeltaGraph();
  }
}

WaveForm::~WaveForm() {

  for(Int_t i=0; i<lineStart.size(); i++)
    delete lineStart.at(i);

  for(Int_t i=0; i<lineStop.size(); i++)
    delete lineStop.at(i);
  
}

void WaveForm::CreateWaveGraph() {
  
  waveGraph.Set(nSample);
  for(Int_t i=0; i < nSample; i++) 
    waveGraph.SetPoint(i,i*sampleToTime,waveVector.at(i));
}

void WaveForm::CreateDeltaGraph() {

  deltaGraph.Set(2*avgWindPoints+deltaVec.size());

  for (Int_t i = 0; i <avgWindPoints; i++)
    deltaGraph.SetPoint(i,i*sampleToTime,0);
  for (Int_t i = 0; i < deltaVec.size(); i++)
    deltaGraph.SetPoint(i+avgWindPoints,(i+avgWindPoints)*sampleToTime,deltaVec.at(i));
  for (Int_t i = 0; i <avgWindPoints; i++)
    deltaGraph.SetPoint(i+avgWindPoints+deltaVec.size(),(i+avgWindPoints+deltaVec.size())*sampleToTime,0);
}

void WaveForm::ComputeDeltaVector(std::vector<Int_t> array) {

  Float_t fwdWind,backWind;
  Float_t delta;
  
  for (Int_t i = avgWindPoints; i < nSample-avgWindPoints; i++) {
    fwdWind=0; backWind=0;
    
    for (Int_t k = 1; k < avgWindPoints; k++) {
      fwdWind+=array.at(i+k);
      backWind+=array.at(i-k);
    }
    
    delta=fwdWind-backWind;
    deltaVec.push_back(delta);
  }
}

void WaveForm::DrawDelta(TCanvas* canvas) {

  canvas->GetSelectedPad();
  deltaGraph.Draw("AWL");

  upThLine.SetX1(0);  upThLine.SetY1(thresUp); upThLine.SetY2(thresUp);
  upThLine.SetX2((deltaVec.size()+2*avgWindPoints)*sampleToTime); 
  downThLine.SetX1(0); downThLine.SetY1(thresDown); downThLine.SetY2(thresDown);
  downThLine.SetX2((deltaVec.size()+2*avgWindPoints)*sampleToTime); 
  
  upThLine.SetLineColor(kRed);
  downThLine.SetLineColor(kRed);
    
  if(numOfPulses>0) {
    upThLine.Draw("SAME");
    downThLine.Draw("SAME");
  }
}

void WaveForm::DrawWave(TCanvas* canvas) {

  canvas->GetSelectedPad();

  waveGraph.GetXaxis()->SetTitle("Sampling Time (ns)");
  waveGraph.GetYaxis()->SetTitle("ADC");
  waveGraph.Draw("AWL");
  
  for(Int_t nPuls=0; nPuls<numOfPulses; nPuls++) {
    lineStart.push_back(new TLine(gateStart.at(nPuls)*sampleToTime,baseline.at(nPuls),
                              gateStart.at(nPuls)*sampleToTime,*max_element(waveVector.begin(),waveVector.end())));
      
    lineStop.push_back(new TLine(gateStop.at(nPuls)*sampleToTime,baseline.at(nPuls),
                             gateStop.at(nPuls)*sampleToTime,*max_element(waveVector.begin(),waveVector.end())));
    lineStart.at(nPuls)->SetLineWidth(2);
    lineStart.at(nPuls)->SetLineColor(kGreen);
    lineStop.at(nPuls)->SetLineWidth(2);
    lineStop.at(nPuls)->SetLineColor(kRed);
    
    lineStart.at(nPuls)->Draw("SAME");
    lineStop.at(nPuls)->Draw("SAME");
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

Int_t WaveForm::ComputeBaseline(std::vector<Int_t> array,Int_t peakNumber) {

  Float_t base=0;
  Int_t start,stop;
  
  if(peakNumber==1)
    start=0;
  else
    start=gateStop.at(peakNumber-2);

  stop=gateStart.at(peakNumber-1);
  
  for (Int_t i = start; i < stop; i++) {
    base+=array.at(i);
  }
  base/=(stop-start);
   
  return base;
 }

Float_t WaveForm::ComputePulseIntegral(std::vector<Int_t> array, Int_t peakNumber) {

  Float_t sum=0;
 
  for (Int_t i = gateStart.at(peakNumber-1); i < gateStop.at(peakNumber-1); i++) {
    sum+=(array.at(i)-baseline.at(peakNumber-1));
  }
  return sum;
}

Float_t WaveForm::ComputePulseMax(std::vector<Int_t> array, Int_t peakNumber) {

  Float_t max=0;
 
  for (Int_t i = gateStart.at(peakNumber-1); i < gateStop.at(peakNumber-1); i++) {
    if(array.at(i)-baseline.at(peakNumber-1)>max)
      max=array.at(i)-baseline.at(peakNumber-1);
  }
  return max;
}

Float_t WaveForm::FindRMS(std::vector<Int_t> array) {

  Float_t mean=0;
  Float_t RMS;
  Int_t sumsquared=0;
  
  for (Int_t i = 0; i < nSample; i++)
    mean+=array.at(i);

  mean=mean/Float_t(nSample);
  
  for (Int_t i = 0; i < nSample; i++)
    sumsquared += (array.at(i)-mean)*(array.at(i)-mean);
  
  RMS = sqrt((Float_t(1)/Float_t(nSample))*(Float_t(sumsquared)));

  return RMS;
}
