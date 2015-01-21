#ifndef MIN_H
#define MIN_H

#include <vector>
#include "Math/Functor.h"
#include <Minuit2/Minuit2Minimizer.h>

using namespace std;

class Minimizer{
 
  ROOT::Math::Functor f;//functor to Minimize
  ROOT::Minuit2::Minuit2Minimizer minuit;//underlying Minuit2 Root Minimizer
  vector<double> step;
  vector<double> variable;
  vector<double> sol;//found solution
  vector<double> err;//errors on the solution
  double minVal;//value of f at 'sol'
  void SetDefaultValues();
  void SetMaths();//set f as the function and set the minuit2 variables to variable

public:
  Minimizer();
  Minimizer(ROOT::Math::Functor f);
  void Process();
  void Update(ROOT::Math::Functor f);
  const ROOT::Math::Functor& getFunctor() const;
  const vector<double>& getSol() const;
  const vector<double>& getErrors() const;
  const double& getMinVal() const;

};

ostream& operator<<(ostream& output, const Minimizer& min);  

#endif