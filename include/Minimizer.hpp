#ifndef MIN_H
#define MIN_H

#include <vector>
#include "Math/Functor.h"
#include <Minuit2/Minuit2Minimizer.h>

class Minimizer{
 
  ROOT::Math::Functor f;//functor to Minimize
  ROOT::Minuit2::Minuit2Minimizer minuit;//underlying Minuit2 Root Minimizer
  std::vector<double> step;
  std::vector<double> variable;
  std::vector<double> sol;//found solution
  std::vector<double> err;//errors on the solution
  double minVal;//value of f at 'sol'
  void setDefaultValues();
  void setMaths();//set f as the function and set the minuit2 variables to variable

public:
  Minimizer();
  Minimizer(ROOT::Math::Functor f);
  void Process();
  void Update(ROOT::Math::Functor f);
  void setInitialValues(const std::vector<double>& variable);
  const ROOT::Math::Functor& getFunctor() const;
  const std::vector<double>& getSol() const;
  const std::vector<double>& getErrors() const;
  const double& getMinVal() const;

};

std::ostream& operator<<(std::ostream& output, const Minimizer& min);  

#endif