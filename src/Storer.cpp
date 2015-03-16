#include "Storer.hpp"

using namespace std;
using namespace::Eigen;
using namespace boost::filesystem;

ostream& operator<<(ostream& output, const Storer& data){
  
  output<<"Files to store from:\n";
  for(const path& filep : data.getFilePaths()) output<<filep<<"\n";
  return output;  
  
}

Storer::Storer(const path& filePath):filePaths({filePath}){

}

Storer::Storer(const vector<path>& filePaths):filePaths(filePaths){

}

bool Storer::itemMatches(TObject* obj, const char* className){

  return string(obj->ClassName()).find(string(className)) != string::npos;
  
}

void Storer::pushAsHist(Data& data, TObject* obj) const{

  if(itemMatches(obj, "TH1D")) data.pushHist(*dynamic_cast<TH1D*>(obj)); //if "GetClassName" contains TH1, the function finds shouldn't return "string::npos"
  else if(itemMatches(obj, "TH1F")){
    
    TH1D htemp;
    dynamic_cast<TH1F*>(obj)->Copy(htemp);//copy the TH1F into a TH1D
    data.pushHist(htemp);
    
  }

}

void Storer::pushAsMatrix(Data& data, TObject* obj) const{//reads a key to retrieve a TH1 and store it in h

  TMatrixD mTemp;
  if(itemMatches(obj, "TMatrixT<double>")){
    
    mTemp.ResizeTo(*dynamic_cast<TMatrixD*>(obj));
    mTemp = *dynamic_cast<TMatrixD*>(obj);//ROOT doesn't support direct equality between matrices of different sizes, and 'Copy' is broken for matrices with different sizes
    
  }
  else if(itemMatches(obj, "TMatrixT<float>")){
    
    mTemp.ResizeTo(*dynamic_cast<TMatrixF*>(obj));
    dynamic_cast<TMatrixF*>(obj)->Copy(mTemp);//copy the TMatrixF into a TMatrixD
    
  }
  data.pushMatrix(Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(mTemp.GetMatrixArray(), mTemp.GetNrows(), mTemp.GetNcols()));

}

void Storer::pushFromCan(Data& data, TObject* readObject) const{//reads a key to retrieve a TH1 from a canvas and store it in h

  TCanvas* can = dynamic_cast<TCanvas*>(readObject);//if one of the 'items' is a canvas, we store the pointer to it in can
  TObject* obj; //then we can loop over the objects of this canvas
  
  TIter objit(can->GetListOfPrimitives()); //iterator over the contents of the canvas
  while((obj = objit())) pushAsHist(data, obj);
  
}

void Storer::pushPath(const path& filePath){
  
  filePaths.push_back(filePath);

}

void Storer::clear(){
  
  filePaths.clear();

}

void Storer::fill(Data& data) const{
  
  TKey* currentKey; //to point to the 'items' in a root file
  TObject* readObject;//to flush the key into a TObject
  TDirectory* CurrentDir = gDirectory->GetDirectory("");//shameful trick 
  
  for(const auto& p : filePaths){
  
    TFile file(p.string().c_str());
    CurrentDir->cd(); //shameful trick, second part
    
    TIter keyit(file.GetListOfKeys()); //iterator over the contents of filePaths
    while((currentKey = dynamic_cast<TKey*>(keyit()))){ //returns the pointer and increment it, when the pointer is not allocated it returns zero so the loop ends
      
      readObject = currentKey->ReadObj();//this is silly to cast a TObject to a TKey to eventually return a TObject with the ReadObj method, but heh, this is ROOT !
      if(itemMatches(readObject, "TH1")) pushAsHist(data, readObject);//reads a key to retrieve a TH1 and store it in h
      else if(itemMatches(readObject, "TMatrix")) pushAsMatrix(data, readObject);
      else if(itemMatches(readObject, "TCanvas")) pushFromCan(data, readObject);//reads a key to retrieve a TH1 from a canvas and store it in h
      
    }
    
  }
  
}

void Storer::setFilePaths(const vector<path>& filePaths){
  
  this->filePaths = filePaths;

}

const vector<path>& Storer::getFilePaths() const{

  return filePaths;
  
}