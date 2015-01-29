#include "Minimiser.hpp"

using namespace std;
using namespace Eigen;

ostream& operator<<(ostream& output, const Minimiser& min){

  output<<"Solution:\n";
  for(unsigned k = 0; k<min.getSol().size(); ++k) output<<setw(8)<<left<<min.getSol().at(k)<<" +/ "<<setw(8)<<left<<min.getErrors().at(k)<<"\n";
  output<<"#######################\n"
    <<"Covariance matrix:\n"<<setfill(' ')<<min.getCovariance()<<"\n"
    <<"Correlation matrix:\n"<<min.getCorrelation()<<"\n"
    <<"#######################\n"
    <<"Value at solution:\n"<<min.getMinVal();
  return output;
  
}

Minimiser::Minimiser():minVal(0){

}

Minimiser::Minimiser(ROOT::Math::Functor f):f(f),minuit(ROOT::Minuit2::kMigrad),step(f.NDim()),variable(f.NDim()),sol(f.NDim()),err(f.NDim()),covariance(f.NDim(),f.NDim()),minVal(0){
  
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
  double covarianceArray[sol.size()];
  minuit.GetCovMatrix(covarianceArray);
  covariance = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(covarianceArray, sol.size(), sol.size());//the covariance matrix is diagonal anyway, so we don't need to read it as Row Major from the ROOT array
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

MatrixXd Minimiser::getCovariance() const{
  
  return covariance;

}

MatrixXd Minimiser::getCorrelation() const{
  
  MatrixXd inverseOfErrors = covariance.diagonal().array().sqrt().inverse().matrix().asDiagonal();
  return  inverseOfErrors * covariance * inverseOfErrors;

}

const double& Minimiser::getMinVal() const{

  return minVal;
  
}

