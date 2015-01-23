# include "Rebinner.hpp"

using namespace std;
using namespace Eigen;

ostream& operator<<(ostream& output, const Rebinner& binsner){

  output<<"binsVector = \n";
  for(const double& binEdge : binsner.getRebin()) output<<binEdge<<"\n";
  return output;
  
}

Rebinner::Rebinner(const Data& data){
  
  buildFrom(data);

}

Rebinner::Rebinner(const vector<double>& bins):bins(bins){
  
}

void Rebinner::Normalise(Hist& h, const int ref){

  if(ref>0 && ref<h.GetNbinsX())
    for(int k = 0; k<h.GetNbinsX(); ++k) h.SetBinContent(k+1, h.GetBinContent(k+1)*h.GetBinWidth(ref)/h.GetBinWidth(k+1));
  
}

vector<double> Rebinner::getCommonElements(const vector<double>& v1, const vector<double>& v2, const double epsilon){
  
  vector<double> common;
  
  for(unsigned k = 0; k<v1.size(); ++k)
    for(unsigned p = 0; p<v2.size(); ++p) if(abs(v1[k] - v2[p])<epsilon) common.push_back(v1[k]);
  
  return common;
  
}

vector<int> Rebinner::getCommonIndices(const vector<double>& v1, const vector<double>& v2, const double epsilon) {

  vector<int> commonIndices;
  
  for(unsigned k = 0; k<v1.size(); ++k)
    for(unsigned p = 0; p<v2.size(); ++p) if(abs(v1[k] - v2[p])<epsilon) commonIndices.push_back(k);
  
  return commonIndices;
  
}

void Rebinner::buildFrom(const Data& data){
  
  bins = data.getHistograms().front().getBins();//overwrite the old bins
  for(const Hist& h : data.getHistograms()) bins = getCommonElements(bins, h.getBins());
  
}

void Rebinner::rebin(Data& data) const{

  double* binsArray = new double[bins.size()];
  copy(bins.begin(), bins.end(), binsArray);
  
  auto itHist = data.getHistogramStartIterator();
  auto itMat = data.getmatricesStartIterator();
  while(itHist != data.getHistograms().end()){

    if(itMat != data.getMatrices().end()){//Rebin the matrix first if it is valid

      vector<int> commonIndices = getCommonIndices(itHist->getBins(), bins);
      MatrixXd reducedMatrix(commonIndices.size()-1, commonIndices.size()-1);//create a temporary matrix
      for(unsigned i = 0; i<reducedMatrix.rows(); ++i)
	for(unsigned j = 0; j<reducedMatrix.cols(); ++j)
	  reducedMatrix(i,j) = itMat->block(commonIndices[i],commonIndices[j],commonIndices[i+1]-commonIndices[i]-1,commonIndices[j+1]-commonIndices[j]-1).sum();//we need to exclude one index (hence the '-1') to avoid overlap with the nth step
    
      *itMat = reducedMatrix;//replace the old matrix with the binsned one
      ++itMat;
      
    }
    
    *itHist = *dynamic_cast<TH1D*>(itHist->Rebin(bins.size()-1, itHist->GetName(), binsArray));//bins the Histogram
    itHist->Scale(1/itHist->Integral());//rescale to unit area
    ++itHist;
    
  }
  
}

bool Rebinner::admissibleRebinFor(Data& data) const{

  vector<double> intersection(bins.size());
  
  for(const Hist& h : data.getHistograms())
    if(intersection.begin() == set_intersection(bins.begin(), bins.end(), h.getBins().begin(), h.getBins().end(), intersection.begin())) return false;//if there is no intersection with the proposed bining, return false

  return true;
  
}

const vector<double>& Rebinner::getRebin() const{
  
  return bins;
  
}
  