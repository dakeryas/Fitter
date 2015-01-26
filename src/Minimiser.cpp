#include "Minimiser.hpp"

using namespace std;

ostream& operator<<(ostream& output, const Minimiser& min){

  output<<"Solution:\n";
  for(unsigned k = 0; k<min.getSol().size(); ++k) output<<min.getSol().at(k)<<" +/ "<<min.getErrors().at(k)<<"\n";
  output<<"Value at solution:\n"<<min.getMinVal();
  return output;
  
}

Minimiser::Minimiser():minVal(0){

}

Minimiser::Minimiser(ROOT::Math::Functor f):f(f),minuit(ROOT::Minuit2::kMigrad),step(f.NDim()),variable(f.NDim()),sol(f.NDim()),err(f.NDim()),minVal(0){
  
  setDefaultValues();
  setMaths();
  
}

void Minimiser::setDefaultValues(){
  
  for(double& s : step) s = 1e-6;//set the default step size
  minuit.SetMaxFunctionCalls(1e9);
  minuit.SetTolerance(1e-6);
  
}

void Minimiser::setMaths(){
  
  minuit.SetFunction(f);
  for(unsigned k = 0; k<f.NDim(); ++k) minuit.SetVariable(k,"x_"+to_string(k),variable[k], step[k]);//set all the variables according to the dimension of the Functor

}

void Minimiser::Process(){

  minuit.Minimize();
  sol = vector<double>(minuit.X(), minuit.X()+f.NDim());
  err = vector<double>(minuit.Errors(), minuit.Errors()+f.NDim());
  minVal = minuit.MinValue();
  
}

void Minimiser::Update(ROOT::Math::Functor f){
  
  this->f = f;
  setMaths();

}

void Minimiser::setInitialValues(const vector<double>& variable){
  
  if(this->variable.size() == variable.size()){
    
    this->variable = variable;
    setMaths();
    
  }

}

const ROOT::Math::Functor& Minimiser::getFunctor() const{
  
  return f;

}

const vector<double>& Minimiser::getSol() const{

  return sol;
  
}

const vector<double>& Minimiser::getErrors() const{

  return err;
  
}

const double& Minimiser::getMinVal() const{

  return minVal;
  
}

