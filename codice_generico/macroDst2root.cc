void macroDst2root(string nFile,string nRun) {

  const Int_t nSili=6;
  const Int_t nMaroc=64;
  const Int_t nDeva=6;
  const Int_t nXinfo=3;
  const Int_t nInfoPlus=2;
  
  string inFileName = "dst_ssh/"+nFile+"_ascii_";
  inFileName += nRun+".dat";
  cout<<"--->Converting file "<<inFileName<<endl;
  string outFileName = "root/"+nFile+"_";
  outFileName += nRun+".root";

  ifstream myfile;

  myfile.open(inFileName.c_str()); 
  
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE"); 
  TTree* outTree = new TTree("outTree","outTree");
  gROOT->cd();
  
  Float_t xPos[nSili];
  Int_t nStrip[nSili];
  Int_t nClu[nSili];
  Float_t maroc[nMaroc];
  Float_t deva[nDeva];
  Float_t xInfo[nXinfo];
  Int_t infoPlus[nInfoPlus];
  Int_t eventNumber;
  
  outTree->Branch("xPos",xPos,TString::Format("xPos[%i]/F",nSili));
  outTree->Branch("nStrip",nStrip,TString::Format("nStrip[%i]/I",nSili));
  outTree->Branch("nClu",nClu,TString::Format("nClu[%i]/I",nSili));
  outTree->Branch("maroc",maroc,TString::Format("maroc[%i]/F",nMaroc));
  outTree->Branch("deva",deva,TString::Format("deva[%i]/F",nDeva)); 
  outTree->Branch("xInfo",xInfo,TString::Format("xInfo[%i]/F",nXinfo));
  outTree->Branch("infoPlus",infoPlus,TString::Format("infoPlus[%i]/F",nInfoPlus));
  outTree->Branch("eventNumber",&eventNumber,"eventNumber/I");


  if (myfile.is_open()) {
    
    while (myfile>>xPos[0]) {
      
      for(Int_t i=1;i<nSili;i++)
        myfile>>xPos[i];
      
      for(Int_t i=0;i<nSili;i++)
	myfile>>nStrip[i];
      
      for(Int_t i=0;i<nSili;i++)
	myfile>>nClu[i];
      
      for(Int_t i=0;i<nMaroc;i++)
	myfile>>maroc[i];
      
      for(Int_t i=0;i<nDeva;i++)
	myfile>>deva[i];
      
      for(Int_t i=0;i<nXinfo;i++)
	myfile>>xInfo[i];
      
      for(Int_t i=0;i<nInfoPlus;i++)
	myfile>>infoPlus[i];
      
      myfile>>eventNumber;
            
      outTree->Fill();
    }
    myfile.close();
  }

  else {
    cout<<"Can't open file!!"<<endl;
  }
  
  outFile->cd();
  outTree->Write();
  
  outFile->Close();
    
  return;
}
