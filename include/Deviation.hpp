#ifndef DEVIATION_H
#define DEVIATION_H

#include <cmath>

class Deviation{

public:
  template <class T, class K>
  static T getForFraction(const T& mean1, const T& mean2, const K& covarianceMatrix);
  template <class T, class K>
  static T getForRatio(const T& mean1, const T& mean2, const K& covarianceMatrix);
  
};

template <class T, class K>
T Deviation::getForFraction(const T& mean1, const T& mean2, const K& covarianceMatrix){

  return sqrt(covarianceMatrix(0, 0)/pow(mean1, 2) + covarianceMatrix(1, 1)/pow(mean2, 2) - 2*covarianceMatrix(0, 1)/mean1/mean2) * std::abs(mean1*mean2)/pow(mean1 + mean2, 2);
  
}

template <class T, class K>
T Deviation::getForRatio(const T& mean1, const T& mean2, const K& covarianceMatrix){

  return sqrt(covarianceMatrix(0, 0)/pow(mean1, 2) + covarianceMatrix(1, 1)/pow(mean2, 2) - 2*covarianceMatrix(0, 1)/mean1/mean2) * std::abs(mean1/mean2);
  
}


#endif
