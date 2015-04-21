#ifndef DEVIATION_H
#define DEVIATION_H

#include <cmath>

template <class T, class K>
class Deviation{

public:
  static T getForFraction(const T& mean1, const T& mean2, const K& covarianceMatrix);
  static T getForRatio(const T& mean1, const T& mean2, const K& covarianceMatrix);
  
};

template <class T, class K>
T Deviation<T, K>::getForFraction(const T& mean1, const T& mean2, const K& covarianceMatrix){

  return sqrt(covarianceMatrix(0, 0)/pow(mean1, 2) + covarianceMatrix(1, 1)/pow(mean2, 2) - 2*covarianceMatrix(0, 1)/mean1/mean2) * mean1*mean2/pow(mean1 + mean2, 2);
  
}

template <class T, class K>
T Deviation<T, K>::getForRatio(const T& mean1, const T& mean2, const K& covarianceMatrix){

  return sqrt(covarianceMatrix(0, 0)/pow(mean1, 2) + covarianceMatrix(1, 1)/pow(mean2, 2) - 2*covarianceMatrix(0, 1)/mean1/mean2) * mean1/mean2;
  
}


#endif
