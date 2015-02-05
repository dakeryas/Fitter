# include "Hist.hpp"

using namespace::std;

ostream& operator<<(ostream& output, const Hist& h){
  
  output<<h.GetName()<<"\n";
  for(unsigned k = 0; k<h.getNumberOfBins(); ++k) output<<"["<<setw(6)<<internal<<h.GetXaxis()->GetBinLowEdge(k+1)<<", "<<setw(6)<<internal<<h.GetXaxis()->GetBinUpEdge(k+1)<<"]"<<setw(8)<<left<<" "
    <<"-->"<<setw(8)<<left<<" "<<setw(12)<<left<<h.GetBinContent(k+1)
    <<"+/-"<<setw(2)<<left<<" "<<setw(12)<<left<<h.GetBinError(k+1)<<"\n";
  return output;
  
}

Hist operator+(Hist h1, const Hist& h2){
  
  h1 += h2;
  return h1;

}

Hist::Hist(){
  
}

Hist::Hist(const TH1D& h):TH1D(h){
  
  fillEdge();
  
}

Hist::Hist(const Hist& other):TH1D(other){
  
  edge = other.edge;
  
}

Hist& Hist::operator+=(const Hist& other){
  
  if(isCompatibleWith(other)){
    
    Add(&other);
    SetName((GetName() + string("_+_") + other.GetName()).c_str());
    
  }
  return *this;

}

Hist& Hist::operator*=(double a){
  
  this->Scale(a);
  SetName((GetName() + string("_x_")+ to_string(a)).c_str());
  return *this;

}

void Hist::fillEdge(){//turn an array into a vector
   
  int n = GetNbinsX();
  double edgeArray[n];
  GetLowEdge(edgeArray);
  edge = vector<double>(edgeArray, edgeArray + n);
  edge.push_back(GetBinLowEdge(n+1));//the low edge of the overflow bin should be the upper edge of the previous bin
  
}

unsigned Hist::getNumberOfBins() const{
  
  return edge.size()-1;//we have one bin less than edges
  
}

const vector<double>& Hist::getEdge() const{

  return edge;
  
}

bool Hist::isCompatibleWith(const Hist& other) const{
  
  return edge == other.edge;

}
