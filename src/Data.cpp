#include "Data.hpp"

ostream& operator<<(ostream& output, const Data& data){
  
  output<<"Histograms:\n";
  for(const TH1D& Hist : data.GetHistograms()) output<<Hist.GetName()<<"\n";
  if(! data.GetMatrices().empty()){
    
    output<<"And matrices:\n";
    for(const MatrixXd& m : data.GetMatrices()) output<<m<<"\n";
  }
  output<<"From files:\n";
  for(const path& filep : data.GetFilePaths()) output<<filep<<"\n";
  return output;  
  
}

Data join(const Data& d1, const Data& d2){

  Data d = d1;
  d.filepaths.insert(d.filepaths.end(),d2.filepaths.begin(),d2.filepaths.end());
  d.h.insert(d.h.end(),d2.h.begin(),d2.h.end());
  if(d.matrices.empty()) d.completeWithEmptyMatrices();
  d.matrices.insert(d.matrices.end(),d2.matrices.begin(),d2.matrices.end());
  return d;
  
}

Data::Data(){

}

Data::Data(const vector<TH1D>& h):h(h.begin(),h.end()){
  
}

Data::Data(const vector<TH1D>& h, const vector<TMatrixD>& matrices):h(h.begin(),h.end()),matrices(matrices.size()){
  
  FillMatricesFromRoot(matrices.begin(), matrices.end());
  
}

Data::Data(const path& search_path){
  
  PathGrabber(search_path);
  StoreData();
  
}

Data::Data(const path& search_path, const string& file_sorter){
  
  PathGrabber(search_path, file_sorter);
  StoreData();
  
}

Data::Data(const Data& other){
  
  filepaths = other.filepaths;
  h = other.h;
  matrices = other.matrices;
  
}

const Data& Data::operator=(const Data& other){
  
  filepaths = other.filepaths;
  h = other.h;
  matrices = other.matrices;
  return *this;

}

void Data::FillMatricesFromRoot(vector<TMatrixD>::const_iterator itStart, vector<TMatrixD>::const_iterator itEnd){

  for(auto it = itStart; it != itEnd; ++it) matrices[it-itStart] = Map<const MatrixXd>(it->GetMatrixArray(), it->GetNrows(), it->GetNcols());
  
}

void Data::completeWithEmptyMatrices(){

  matrices.insert(matrices.end(), h.size()-matrices.size(), MatrixXd(0,0));
  
}

bool Data::PathMatches(const directory_iterator& it, const string& file_sorter){
  
  return it->path().filename().string().find(file_sorter) != string::npos;
  
}

void Data::FillIfRoot(const directory_iterator& it){//fills 'found_paths' with the path pointed by 'it' if 'it' refers to a file
  
  if(is_directory(it->status()) == 0 && PathMatches(it, ".root")) filepaths.push_back(it->path()); ; //find returns the position of 'name' in the path, so if it returned a position indeed(namely not npos), we found something
  
}

void Data::FillIfRootAnd(const directory_iterator& it, const string& file_sorter){//fills 'found_paths' with the path pointed by 'it' if 'it' refers to a file
  
  if(is_directory(it->status()) == 0 && PathMatches(it, ".root") && PathMatches(it, file_sorter)) filepaths.push_back(it->path()); ; //find returns the position of 'name' in the path, so if it returned a position indeed(namely not npos), we found something
  
}

void Data::PathGrabber(const path& search_path){ //retrieves the path of all files in a directory

  directory_iterator end; //the default constructor creates an end iterator which cannot be reached unless nothing was found
  for(directory_iterator it(search_path); it!= end; ++it) FillIfRoot(it);
  
}

void Data::PathGrabber(const path& search_path, const string& file_sorter){ //retrieves the path of all the files in a directory that match the name file_sorter

  directory_iterator end; //the default constructor creates an end iterator which cannot be reached unless nothing was found
  for(directory_iterator it(search_path); it!= end; ++it) FillIfRootAnd(it, file_sorter);
  
}

bool Data::ItemMatches(TObject* obj, const char* className){

  return string(obj->ClassName()).find(string(className)) != string::npos;
  
}

void Data::StoreAsHist(TObject* obj){

  if(ItemMatches(obj, "TH1D")) h.push_back(*dynamic_cast<TH1D*>(obj)); //if "GetClassName" contains TH1, the function finds shouldn't return "string::npos"
  else if(ItemMatches(obj, "TH1F")){
    
    TH1D htemp;
    dynamic_cast<TH1F*>(obj)->Copy(htemp);//copy the TH1F into a TH1D
    h.push_back(htemp);
    
  }

}

void Data::StoreAsMatrix(TObject* obj){//reads a key to retrieve a TH1 and store it in h

  TMatrixD mTemp;
  if(ItemMatches(obj, "TMatrixT<double>")){
    
    mTemp.ResizeTo(*dynamic_cast<TMatrixD*>(obj));
    mTemp = *dynamic_cast<TMatrixD*>(obj);//ROOT doesn't support direct equality between matrices of different sizes, and 'Copy' is broken for matrices with different sizes
    
  }
  else if(ItemMatches(obj, "TMatrixT<float>")){
    
    mTemp.ResizeTo(*dynamic_cast<TMatrixF*>(obj));
    dynamic_cast<TMatrixF*>(obj)->Copy(mTemp);//copy the TMatrixF into a TMatrixD
    
  }
  matrices.push_back(Map<const MatrixXd>(mTemp.GetMatrixArray(), mTemp.GetNrows(), mTemp.GetNcols()));

}

void Data::StoreFromCan(TObject* read_object){//reads a key to retrieve a TH1 from a canvas and store it in h

  TCanvas* can = dynamic_cast<TCanvas*>(read_object);//if one of the 'items' is a canvas, we store the pointer to it in can
  TObject* obj; //then we can loop over the objects of this canvas
  
  TIter objit(can->GetListOfPrimitives()); //iterator over the contents of the canvas
  while((obj = objit())) StoreAsHist(obj);
  
}

void Data::StoreData(){ //store all TH1D Histograms of all files into the vector h
  
  TKey* current_key; //to point to the 'items' in a root file
  TObject* read_object;//to flush the key into a TObject
  TDirectory* CurrentDir = gDirectory->GetDirectory("");//shameful trick 
  
  for(const path& p : filepaths){
  
    TFile file(p.string().c_str());
    CurrentDir->cd(); //shameful trick, second part
    
    TIter keyit(file.GetListOfKeys()); //iterator over the contents of filepaths
    while((current_key = dynamic_cast<TKey*>(keyit()))){ //returns the pointer and increment it, when the pointer is not allocated it returns zero so the loop ends
      
      read_object = current_key->ReadObj();//this is silly to cast a TObject to a TKey to eventually return a TObject with the ReadObj method, but heh, this is ROOT !
      if(ItemMatches(read_object, "TH1")) StoreAsHist(read_object);//reads a key to retrieve a TH1 and store it in h
      else if(ItemMatches(read_object, "TMatrix")) StoreAsMatrix(read_object);
      else if(ItemMatches(read_object, "TCanvas")) StoreFromCan(read_object);//reads a key to retrieve a TH1 from a canvas and store it in h
      
    }
    
  }
  
}

const vector<path>& Data::GetFilePaths() const{

  return filepaths;
  
}

const vector<Hist>& Data::GetHistograms() const{
  
  return h;
  
}

const vector<TH1D> Data::GetTH1DCopies() const{
  
  return vector<TH1D>(h.begin(),h.end());
  
}

const vector<MatrixXd>& Data::GetMatrices() const{
  
  return matrices;
  
}

TMatrixD Data::GetRootMatrixCopy(unsigned i) const{
  
  auto itMat = matrices.begin()+i;
  if(itMat != matrices.end()) return TMatrixD (itMat->rows(), itMat->cols(), itMat->data());
  else return TMatrixD();
  
}

unsigned int Data::GetSize() const{

  return max({h.size(), matrices.size()});
  
}
