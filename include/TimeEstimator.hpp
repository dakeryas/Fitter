#ifndef TIME_ESTIMATOR_H
#define TIME_ESTIMATOR_H

#include <cmath>
#include "Chi.hpp"
#include "Minimizer.hpp"

class TimeEstimator{//class to estimate the amount of additional data necessary to have non-zero minimum of chiSquared to nSigma

  double dataFactor;
  double dataIncrease;
  
public:  
  TimeEstimator(double dataFactor = 1.01);
  const double& getDataFactor() const;
  const double& getDataIncrease() const;
  double getRelativeTime(Minimizer& min, Chi& chiSquared, unsigned nSigma);//estimates the amount of additional data necessary to have non-zero minimum of chiSquared with a nSigma significance
  void setDataFactor(double dataFactor);
  
};

std::ostream& operator<<(std::ostream& output, const TimeEstimator& timeEstimator);

#endif


