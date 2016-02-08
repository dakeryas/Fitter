#ifndef MIN_H
#define MIN_H

#include <vector>
#include <iomanip>
#include <Eigen/Core>
#include "Math/Functor.h"
#include <Minuit2/Minuit2Minimizer.h>
#include <TError.h> //to access gErrorIgnoreLevel to prevent info messages printing in muli-threading

class Minimiser{
 
  ROOT::Math::Functor f;//functor to Minimize
  ROOT::Minuit2::Minuit2Minimizer minuit;//underlying Minuit2 Root Minimizer
  std::vector<double> step;
  std::vector<double> variable;
  std::vector<double> sol;//found solution
  std::vector<double> err;//errors on the solution
  Eigen::MatrixXd covariance;//covariance matrix
  double minVal;//value of f at 'sol'
  void setDefaultValues();
  void setMaths();//set f as the function and set the minuit2 variables to variable

public:
  Minimiser();
  Minimiser(ROOT::Math::Functor f);
  void Process();
  void Update(ROOT::Math::Functor f);
  void setInitialValues(std::vector<double> variable);
  const ROOT::Math::Functor& getFunctor() const;
  const std::vector<double>& getSol() const;
  const std::vector<double>& getErrors() const;
  Eigen::MatrixXd getCovariance() const;
  Eigen::MatrixXd getCorrelation() const;
  const double& getMinVal() const;

};

std::ostream& operator<<(std::ostream& output, const Minimiser& min);  

#endif