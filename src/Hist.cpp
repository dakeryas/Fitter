# include "Hist.hpp"

using namespace::std;

Hist operator+(Hist h1, const Hist& h2){
  
  h1 += h2;
  return h1;

}

Hist::Hist(){
  
}

Hist::Hist(const TH1D& h):TH1D(h){
  
  FillBins();
  
}

Hist::Hist(const Hist& other):TH1D(other){
  
  bins = other.bins;
  
}

Hist& Hist::operator+=(const Hist& other){
  
  if(isCompatibleWith(other)) Add(&other);
  return *this;

}

Hist& Hist::operator*=(double a){
  
  (*this) * a;
  return *this;

}

void Hist::FillBins(){//turn an array into a vector
   
  int n = GetNbinsX();
  double edge[n];
  GetLowEdge(edge);
  bins = vector<double>(edge, edge + n);
  bins.push_back(GetBinLowEdge(n+1));//the low edge of the overflow bin should be the upper edge of the previous bin
  
}

const vector<double>& Hist::GetBins() const{

  return bins;
  
}

bool Hist::isCompatibleWith(const Hist& other) const{
  
  return bins == other.bins;

}
