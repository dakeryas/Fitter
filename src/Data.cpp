#include "Data.hpp"

using namespace std;
using namespace::Eigen;

enum Correlation {independent, correlated};
MatrixXd covariance(const MatrixXd& m1, const MatrixXd& m2, Correlation correlation = correlated){//shortcut from real maths: the variables X1 and X2 are represented by their matrices 'm1' and 'm2' to compute Cov(X1, X2)

  if(correlation == correlated){
    
    return m1.diagonal().array().sqrt().matrix()*m2.diagonal().array().sqrt().matrix().transpose();//take the square root of the variances to get the sigma's and return sigma1 * sigma2.transpose()
    
  }
  else return MatrixXd(m1.cols(), m1.cols()).setZero();
  
}

ostream& operator<<(ostream& output, const Data& data){
  
  output<<"Histograms:\n";
  for(const TH1D& Hist : data.getHistograms()) output<<Hist<<"\n";
  if(! data.getMatrices().empty()){
    
    output<<"And matrices:\n";
    for(const MatrixXd& m : data.getMatrices()) output<<m<<"\n";
    
  }
  return output;
  
}

Data operator+(Data d1, const Data& d2){
  
  d1 += d2;
  return d1;

}

Data operator*(Data data, double a){
  
  data *= a;
  return data;

}

Data operator*(double a, Data data){
  
  return data * a;

}

Data join(Data d1, const Data& d2){

  for(auto it = d2.getHistograms().begin(); it != d2.getHistograms().end(); ++it) d1.pushHist(*it);
  if(d1.getMatrices().empty() && !d2.getMatrices().empty()) d1.completeWithEmptyMatrices();
  else for(auto it = d2.getMatrices().begin(); it != d2.getMatrices().end(); ++it) d1.pushMatrix(*it);
  return d1;
  
}

Data::Data(){

}

Data::Data(const vector<TH1D>& histograms):histograms(histograms.begin(),histograms.end()){
  
}

Data::Data(const vector<TH1D>& histograms, const vector<TMatrixD>& matrices):histograms(histograms.begin(),histograms.end()),matrices(matrices.size()){
  
  fillMatricesFromRoot(matrices.begin(), matrices.end());
  
}

Data& Data::operator+=(const Data& other){
  
  for(pair<vector<MatrixXd>::iterator, vector<MatrixXd>::const_iterator> itPair(matrices.begin(), other.matrices.begin()); itPair.first != matrices.end() && itPair.second != other.matrices.end(); ++itPair.first, ++itPair.second)
    *itPair.first += *itPair.second + covariance(*itPair.first, *itPair.second) + covariance(*itPair.second, *itPair.first);
  
  for(pair<vector<Hist>::iterator, vector<Hist>::const_iterator> itPair(histograms.begin(), other.histograms.begin()); itPair.first != histograms.end() && itPair.second != other.histograms.end(); ++itPair.first, ++itPair.second)
    *itPair.first += *itPair.second;
  
  for(pair<vector<Hist>::iterator, vector<MatrixXd>::const_iterator> itPair(histograms.begin(), matrices.begin()); itPair.first != histograms.end() && itPair.second != matrices.end(); ++itPair.first, ++itPair.second)
    itPair.first->setErrorsFrom(*itPair.second);//if they are covariance matrices, once the matrices have been properly summed, set the errors on the histograms from them
  
  return *this;

}

Data& Data::operator*=(double a){
  
  for(auto it = histograms.begin(); it != histograms.end(); ++it) *it *= a;
  for(auto it = matrices.begin(); it != matrices.end(); ++it) *it *= pow(a,2);//for the covariance matrices be careful to use the square of a
  
  return *this;

}

void Data::completeWithEmptyMatrices(){

  matrices.insert(matrices.end(), histograms.size()-matrices.size(), MatrixXd(0,0));
  
}

void Data::pushHist(const Hist& hist){
  
  histograms.push_back(hist);

}

void Data::pushMatrix(const MatrixXd& matrix){
  
  matrices.push_back(matrix);

}

void Data::push(const Hist& hist, const MatrixXd& matrix){
  
  pushHist(hist);
  pushMatrix(matrix);

}

void Data::fillMatricesFromRoot(vector<TMatrixD>::const_iterator itStart, vector<TMatrixD>::const_iterator itEnd){
  
  for(auto it = itStart; it != itEnd; ++it) matrices[it-itStart] = Map<const Matrix<double, Dynamic, Dynamic, RowMajor>>(it->GetMatrixArray(), it->GetNrows(), it->GetNcols());
  
}

const vector<Hist>& Data::getHistograms() const{
  
  return histograms;
  
}

vector<Hist>::iterator Data::getHistogramStartIterator(){
  
  return histograms.begin();

}

const vector<TH1D> Data::getTH1DCopies() const{
  
  return vector<TH1D>(histograms.begin(),histograms.end());
  
}

const vector<MatrixXd>& Data::getMatrices() const{
  
  return matrices;
  
}

vector<MatrixXd>::iterator Data::getmatricesStartIterator(){

  return matrices.begin();
  
}

TMatrixD Data::getRootMatrixCopy(unsigned i) const{
  
  auto itMat = matrices.begin()+i;
  if(itMat != matrices.end()) return TMatrixD (itMat->rows(), itMat->cols(), itMat->data());
  else return TMatrixD();
  
}

unsigned int Data::getNumberOfBins() const{
  
  return (*min_element(histograms.begin(), histograms.end(), [](const Hist& h1, const Hist& h2){return h1.getNumberOfBins() < h2.getNumberOfBins();})).getNumberOfBins();

}

unsigned int Data::getSize() const{

  return max({histograms.size(), matrices.size()});
  
}

void Data::clear(){
  
  histograms.clear();
  matrices.clear();

}
