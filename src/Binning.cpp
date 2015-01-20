#include "Binning.hpp"

std::ostream& operator<<(std::ostream& output, const Binning& binning){

  output<<"Steps per percent: "<<binning.getNumberOfStepsPerPercent()<<"\t Total number of steps: "<<binning.getNumberOfSteps()
    <<"\nMin percent: "<<binning.getMinPercent()<<"\t Max percent : "<<binning.getMaxPercent()<<"\n";
  return output;
  
}


Binning::Binning(unsigned int numberOfStepsPerPercent, double minPercent, double maxPercent):numberOfStepsPerPercent(numberOfStepsPerPercent),minPercent(minPercent),maxPercent(maxPercent), bins(getNumberOfSteps()){

  fillBins();
  
}

void Binning::fillBins(){
  
  for(auto it = bins.begin(); it != bins.end(); ++it) *it = static_cast<double>((it-bins.begin())+minPercent*numberOfStepsPerPercent)/100/numberOfStepsPerPercent;
  
}


const unsigned& Binning::getNumberOfStepsPerPercent() const{
  
  return numberOfStepsPerPercent;

}

const double& Binning::getMinPercent() const{
  
  return minPercent;

}

const double& Binning::getMaxPercent() const{
  
  return maxPercent;

}

double Binning::getRange() const{
  
  return maxPercent - minPercent;

}

unsigned Binning::getFlooredRange() const{
  
  return floor(maxPercent - minPercent);

}

unsigned Binning::getNumberOfSteps() const{
  
  return getFlooredRange()*numberOfStepsPerPercent + 1;

}

double Binning::getValue(unsigned int i) const{
  
  try{
    
    return bins.at(i);
    
  }
  catch(int e){
    
    return -1;
    
  }

}

const double* Binning::getDataBins() const{
  
  return bins.data();

}

void Binning::setBinning(unsigned int numberOfStepsPerPercent, double minPercent, double maxPercent){
  
  this->numberOfStepsPerPercent = numberOfStepsPerPercent;
  this->minPercent = minPercent;
  this->maxPercent = maxPercent;
  fillBins();

}











