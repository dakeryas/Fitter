#include "Minimizer.hpp"

ostream& operator<<(ostream& output, const Minimizer& min){

  output<<"Solution:\n";
  for(unsigned k = 0; k<min.getSol().size(); ++k) output<<min.getSol().at(k)<<" +/ "<<min.getErrors().at(k)<<"\n";
  output<<"Value at solution:\n"<<min.getMinVal();
  return output;
  
}

Minimizer::Minimizer():minVal(0){

}

Minimizer::Minimizer(ROOT::Math::Functor f):f(f),minuit(ROOT::Minuit2::kMigrad),step(f.NDim()),variable(f.NDim()),sol(f.NDim()),err(f.NDim()),minVal(0){
  
  SetDefaultValues();
  SetMaths();
  
}

void Minimizer::SetDefaultValues(){
  
  for(double& s : step) s = 1e-6;//set the default step size
  minuit.SetMaxFunctionCalls(1e9);
  minuit.SetTolerance(1e-6);
  
}

void Minimizer::SetMaths(){
  
  minuit.SetFunction(f);
  for(unsigned k = 0; k<f.NDim(); ++k) minuit.SetVariable(k,"x_"+to_string(k),variable[k], step[k]);//set all the variables according to the dimension of the Functor

}

void Minimizer::Process(){
  
  minuit.Minimize();
  sol = vector<double>(minuit.X(), minuit.X()+f.NDim());
  err = vector<double>(minuit.Errors(), minuit.Errors()+f.NDim());
  minVal = minuit.MinValue();
  
}

void Minimizer::Update(ROOT::Math::Functor f){
  
  this->f = f;
  SetMaths();

}


const vector<double>& Minimizer::getSol() const{

  return sol;
  
}

const vector<double>& Minimizer::getErrors() const{

  return err;
  
}

const double& Minimizer::getMinVal() const{

  return minVal;
  
}

