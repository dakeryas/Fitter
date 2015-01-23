# include "DataRebin.hpp"

using namespace std;
using namespace Eigen;

ostream& operator<<(ostream& output, const DataRebin& rebinner){

  output<<"rebinVector = \n";
  for(const double& binEdge : rebinner.rebin) output<<binEdge<<"\n";
  return output;
  
}

DataRebin::DataRebin():rebinArray(new double){
  
}

DataRebin::DataRebin(const Data& data):data(data){
  
  FillRebin();

}

DataRebin::DataRebin(const Data& data, const vector<double>& rebin):data(data),rebin(rebin){
  
  if(!AdmissibleRebin()) FillRebin();//if the rebin passed in argument is not admissible, rebin the data following the normal procedure
  else AllocAndInitArray();
  
}


DataRebin::DataRebin(const DataRebin& other){
  
  data = other.data;
  rebin = other.rebin;
  AllocAndInitArray();
  
}

DataRebin& DataRebin::operator=(const DataRebin& other){

  data = other.data;
  rebin = other.rebin;
  delete rebinArray;//release the memory
  AllocAndInitArray();//before allocating the new relevant size based on the updated rebin vector
  return *this;
  
}


DataRebin::~DataRebin(){

  delete rebinArray;
  
}

void DataRebin::AllocAndInitArray(){
  
  rebinArray = new double[rebin.size()];
  copy(rebin.begin(), rebin.end(), rebinArray);//also create an array copy
  
}

vector<double> DataRebin::GetCommonElements(const vector<double>& v1, const vector<double>& v2, const double epsilon){
  
  vector<double> common;
  
  for(unsigned k = 0; k<v1.size(); ++k)
    for(unsigned p = 0; p<v2.size(); ++p) if(abs(v1[k] - v2[p])<epsilon) common.push_back(v1[k]);
  
  return common;
  
}

vector<int> DataRebin::GetCommonIndices(const vector<double>& v1, const vector<double>& v2, const double epsilon) {

  vector<int> commonIndices;
  
  for(unsigned k = 0; k<v1.size(); ++k)
    for(unsigned p = 0; p<v2.size(); ++p) if(abs(v1[k] - v2[p])<epsilon) commonIndices.push_back(k);
  
  return commonIndices;
  
}


bool DataRebin::AdmissibleRebin() const{

  vector<double> intersection(rebin.size());
  
  for(const Hist& h : data.getHistograms())
    if(intersection.begin() == set_intersection(rebin.begin(), rebin.end(), h.getBins().begin(), h.getBins().end(), intersection.begin())) return false;//if there is no intersection with the proposed bining, return false

  return true;
  
}

void DataRebin::FillRebin(){
  
  rebin = data.getHistograms().front().getBins();
  for(const Hist& h : data.getHistograms()) rebin = GetCommonElements(rebin, h.getBins());
  AllocAndInitArray();
  
}

void DataRebin::ApplyRebin(){

  auto itHist = data.getHistogramStartIterator();
  auto itMat = data.getmatricesStartIterator();
  while(itHist != data.getHistograms().end()){

    if(itMat != data.getMatrices().end()){//Rebin the matrix first if it is valid

      vector<int> commonIndices = GetCommonIndices(itHist->getBins(), rebin);
      MatrixXd reducedMatrix(commonIndices.size()-1, commonIndices.size()-1);//create a temporary matrix
      for(unsigned i = 0; i<reducedMatrix.rows(); ++i)
	for(unsigned j = 0; j<reducedMatrix.cols(); ++j)
	  reducedMatrix(i,j) = itMat->block(commonIndices[i],commonIndices[j],commonIndices[i+1]-commonIndices[i]-1,commonIndices[j+1]-commonIndices[j]-1).sum();//we need to exclude one index (hence the '-1') to avoid overlap with the nth step
    
      *itMat = reducedMatrix;//replace the old matrix with the rebinned one
      ++itMat;
      
    }
    
    *itHist = *dynamic_cast<TH1D*>(itHist->Rebin(rebin.size()-1, itHist->GetName(), rebinArray));//rebin the Histogram
    itHist->Scale(1/itHist->Integral());//rescale to unit area
    ++itHist;
    
  }
  
}

void DataRebin::Normalise(Hist& h, const int ref){

  if(ref>0 && ref<h.GetNbinsX())
    for(int k = 0; k<h.GetNbinsX(); ++k) h.SetBinContent(k+1, h.GetBinContent(k+1)*h.GetBinWidth(ref)/h.GetBinWidth(k+1));
  
}

const Data& DataRebin::GetData() const{
  
  return data;
  
}

const vector<double>& DataRebin::GetRebin() const{
  
  return rebin;
  
}
  