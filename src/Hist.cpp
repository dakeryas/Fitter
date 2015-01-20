# include "Hist.hpp"

Hist::Hist(){
  
}

Hist::Hist(const TH1D& h):TH1D(h){
  
  FillBins();
  
}

Hist::Hist(const Hist& other):TH1D(other){
  
  bins = other.bins;
  
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
