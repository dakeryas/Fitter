#include "TimeEstimator.hpp"

std::ostream& operator<<(std::ostream& output, const TimeEstimator& timeEstimator){

  output<<"Multiplying the amount of data by: "<<timeEstimator.getDataFactor()<<"\n Increased the amout of data "<<timeEstimator.getDataIncrease()<<" times";
  return output;
  
}


TimeEstimator::TimeEstimator(double dataFactor):dataFactor(dataFactor),dataIncrease(0){
  
}

const double& TimeEstimator::getDataFactor() const{
  
  return dataFactor;

}

const double& TimeEstimator::getDataIncrease() const{
  
  return dataIncrease;

}

double TimeEstimator::getRelativeTime(Minimizer& min, Chi& chiSquared, unsigned int nSigma){
  
  min.Update(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));
  min.Process();
  
  double dataIncrease = 0;

  while(min.getSol().front()-nSigma*min.getErrors().front()<0){//test zero-ness to nSigma

    chiSquared.SetDataErr(chiSquared.getDataErr()/std::sqrt(dataFactor));//double the data if the errors are to large to exclude a zero fraction of the first element
    min.Update(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));
    min.Process();
    ++dataIncrease;
    
  }
  
  return std::pow(dataFactor, dataIncrease);

}

void TimeEstimator::setDataFactor(double dataFactor){
  
  this->dataFactor = dataFactor;
  dataIncrease = 0;

}



