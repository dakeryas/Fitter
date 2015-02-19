# include "Rebinner.hpp"

using namespace std;
using namespace Eigen;

ostream& operator<<(ostream& output, const Rebinner& edgener){

  output<<"Bin edges = \n";
  for(const double& binEdge : edgener.getRebin()) output<<binEdge<<"\n";
  return output;
  
}

Rebinner::Rebinner(const Data& data){
  
  buildFrom(data);

}

Rebinner::Rebinner(const vector<double>& edge):edge(edge){
  
}

void Rebinner::Normalise(Hist& h, const int ref){

  if(ref>0 && ref<h.GetNbinsX())
    for(int k = 0; k<h.GetNbinsX(); ++k) h.SetBinContent(k+1, h.GetBinContent(k+1)*h.GetBinWidth(ref)/h.GetBinWidth(k+1));
  
}

vector<double> Rebinner::getCommonElements(const vector<double>& v1, const vector<double>& v2, double epsilon){
  
  vector<double> common;
  set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(common), [&](double a, double b){if(abs(a-b)<epsilon) return false; else return a < b;});//if the numbers are too close to each other, decrete that none is strictly greater than the other

  return common;
  
}

vector<unsigned> Rebinner::getCommonIndices(const vector<double>& v1, const vector<double>& v2, double epsilon) {

  vector<unsigned> commonIndices;
  
  for(unsigned k = 0; k<v1.size(); ++k)
    for(unsigned p = 0; p<v2.size(); ++p) if(abs(v1[k] - v2[p])<epsilon) commonIndices.push_back(k);
  
  return commonIndices;
  
}

void Rebinner::buildFrom(const Data& data){
  
  edge = data.getHistograms().front().getEdge();//overwrite the old edge
  for(const Hist& h : data.getHistograms()) edge = getCommonElements(edge, h.getEdge());
  
}

void Rebinner::excludeBinsAbove(const double& newUpEdge){

  auto itGarbage = remove_if(edge.begin(), edge.end(), [&](const double& currentEdge){return currentEdge > newUpEdge;});
  edge.erase(itGarbage, edge.end());

}

void Rebinner::excludeBinsBelow(const double& newLowEdge){

  auto itGarbage = remove_if(edge.begin(), edge.end(), [&](const double& currentEdge){return currentEdge < newLowEdge;});//remove_if stores the garbage at the end of the vector and returns an iterator to the start of the garbage
  edge.erase(itGarbage, edge.end());//the garbage must then be removed with 'erase'
  
}

void Rebinner::squeezeBinning(unsigned factor){

  int counter = -1;
  auto itGarbage = remove_if(edge.begin(), edge.end(), [&](double){++counter; return counter % factor != 0;});//remove all the indices that are not multiples from factor
  edge.erase(itGarbage, edge.end());
  
}

void Rebinner::rebinMatrix(MatrixXd& matrix, const vector<double>& edgeForMatrix) const{
  
  vector<unsigned> commonIndices = getCommonIndices(edgeForMatrix, edge);

  MatrixXd reducedMatrix(commonIndices.size()-1, commonIndices.size()-1);//create a temporary matrix
  for(unsigned i = 0; i<reducedMatrix.rows(); ++i)
    for(unsigned j = 0; j<reducedMatrix.cols(); ++j)
      reducedMatrix(i,j) = matrix.block(commonIndices[i],commonIndices[j],commonIndices[i+1]-commonIndices[i],commonIndices[j+1]-commonIndices[j]).sum();
    
  matrix = reducedMatrix;//replace the old matrix with the rebinned one

}

void Rebinner::rebin(Data& data) const{
    
  auto itHist = data.getHistogramStartIterator();
  auto itMat = data.getmatricesStartIterator();
  while(itHist != data.getHistograms().end()){

=======

}

void Rebinner::rebin(Data& data) const{
    
  auto itHist = data.getHistogramStartIterator();
  auto itMat = data.getmatricesStartIterator();
  while(itHist != data.getHistograms().end()){

    if(itMat != data.getMatrices().end()){//rebin the matrices
      
      rebinMatrix(*itMat, itHist->getEdge());//rebin the matrix first before the edges of the hist change with its rebin
      *itHist = *dynamic_cast<TH1D*>(itHist->Rebin(edge.size()-1, itHist->GetName(), edge.data()));//rebin the Histogram
      itHist->setErrorsFrom(*itMat);//if there are covariance matrices available, they should be used for the rebin and the errors recomputed
      ++itMat;//we're now done with the matrix
      
    }
    else *itHist = *dynamic_cast<TH1D*>(itHist->Rebin(edge.size()-1, itHist->GetName(), edge.data()));//rebin the Histogram
    
    itHist->Scale(1/itHist->Integral());//rescale to unit area
    ++itHist;
    
  }
  
}

bool Rebinner::admissibleRebinFor(Data& data) const{

  vector<double> intersection(edge.size());
  
  for(const Hist& h : data.getHistograms())
    if(intersection.begin() == set_intersection(edge.begin(), edge.end(), h.getEdge().begin(), h.getEdge().end(), intersection.begin())) return false;//if there is no intersection with the proposed bining, return false

  return true;
  
}

const vector<double>& Rebinner::getRebin() const{
  
  return edge;
  
}

unsigned int Rebinner::getNumberOfBins() const{
  
  if(!edge.empty()) return edge.size() -1;
  else return 0;
  
}