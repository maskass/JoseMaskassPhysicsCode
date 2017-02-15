//#include <string>
//#include <vector>

std::string IntegerToString(int iinput){
  std::stringstream converter;
  converter<<iinput;
  return converter.str();
}

std::string FloatToString(float iinput2){
  std::stringstream converter2;
  converter2<<iinput2;
  return converter2.str();
}

Int_t FileLineCounter(std::string filename) {
  ifstream inFile(filename.c_str()); 
  return count(std::istreambuf_iterator<char>(inFile), 
	       std::istreambuf_iterator<char>(), '\n');
}

// void ReadDeadStrip(std::string fileName,std::vector<Int_t> &vector) {

//   Int_t lines=FileLineCounter(fileName.c_str());
  
//   ifstream deadStripFile;

//   deadStripFile.open(fileName.c_str());

//   Float_t temp1,temp2;
  
//   if (deadStripFile.is_open()) {
            
//     for(Int_t i=0;i<lines;i++) {
      
//       deadStripFile>>temp1>>temp2;

//       vector.push_back(Int_t(temp2));
      
//     }
    
//     deadStripFile.close();
//   }
  
//   else 
//     std::cout<<"Can't open strip status file!!"<<std::endl;
  

//   return;
// }

// void BookHisto1d(int mod,vector<TH1F*> &histo,int numID,string title,
// 		 int nBins,float xMin,float xMax) {
  
//   int i;
//   string name1,name2;
  
//   for(i=0;i<mod;i++) {
    
//     name1="histo" + IntegerToString(numID+i+1);
//     name2=title + IntegerToString(i+1);
//     histo.push_back(new TH1F(name1.c_str(),name2.c_str(),nBins,xMin,xMax));
    
//   }
  
//   return;
// }



