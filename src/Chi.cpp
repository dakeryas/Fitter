#include "Chi.hpp"

Chi::Chi(){

}

Chi::Chi(const VectorXd& dataToFit, const MatrixXd& simulations):dataToFit(dataToFit),dataErrors(dataToFit.array().sqrt()),simulations(simulations){
  
}

Chi::Chi(const VectorXd& dataToFit, const VectorXd& dataErrors, const MatrixXd& simulations):dataToFit(dataToFit),dataErrors(dataErrors),simulations(simulations){
  
}

Chi::Chi(const VectorXd& dataToFit, const VectorXd& dataErrors, const MatrixXd& simulations, const vector< MatrixXd >& covariances):dataToFit(dataToFit),dataErrors(dataErrors),simulations(simulations),covariances(covariances){

}

Chi::Chi(vector<TH1D>::const_iterator itDataToFit, vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations):dataToFit(Map<const VectorXd>(itDataToFit->GetArray()+1,itDataToFit->GetNbinsX())),dataErrors(itDataToFit->GetNbinsX()),simulations(getMaxNbins(itStartSimulations,itEndSimulations), itEndSimulations-itStartSimulations){

  FillErrors(*itDataToFit);
  FillSim(itStartSimulations, itEndSimulations);
  
}

Chi::Chi(const TH1D& dataToFit, const vector<TH1D>& simHist):dataToFit(Map<const VectorXd>(dataToFit.GetArray()+1, dataToFit.GetNbinsX())),dataErrors(dataToFit.GetNbinsX()),simulations(getMaxNbins(simHist.begin(),simHist.end()), simHist.size()){

  FillErrors(dataToFit);
  FillSim(simHist.begin(), simHist.end());
  
}

Chi::Chi(const Data& dataToFit, const Data& simulations):Chi(dataToFit.getHistograms().front(), simulations.getTH1DCopies()){

  covariances = simulations.getMatrices();
  
}

double Chi::operator()(const double* args) const{
 
  const VectorXd fractions = Map<const VectorXd>(args, simulations.cols());// there are as many fractions as there are different spectra, hence as there are different columns in 'simulations'
  MatrixXd Var = dataErrors.array().pow(2).matrix().asDiagonal();//create a diagonal matrix from the squares of the errors on the data
  for(unsigned k = 0; k<covariances.size(); ++k) Var+=covariances[k]*pow(fractions(k),2);//add the covariances for the simulations
  return (dataToFit-simulations*fractions).transpose()*Var.inverse()*(dataToFit-simulations*fractions);

}

void Chi::SetData(const VectorXd& dataToFit){
  
  this->dataToFit = dataToFit;

}

void Chi::SetDataErr(const VectorXd& dataErrors){
  
  this->dataErrors = dataErrors;

}

const VectorXd& Chi::getDataErr() const{

  return dataErrors;
  
}

const MatrixXd& Chi::getSimulations() const{

  return simulations;
  
}

int Chi::getMaxNbins(vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations){

  int maxNBinsX = itStartSimulations->GetNbinsX();
  for (auto it = itStartSimulations; it != itEndSimulations; ++it) if(it->GetNbinsX()>maxNBinsX) maxNBinsX=it->GetNbinsX();
  return maxNBinsX;
  
}

unsigned Chi::getNumberOfFreeParameters() const{
  
  return simulations.cols();
  
}

void Chi::FillErrors(const TH1D& dataToFit){
  
  for(int i = 0; i<dataErrors.size(); ++i) dataErrors(i) = dataToFit.GetBinError(i+1);
  
}

void Chi::FillSim(vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations){
  
  for(auto it = itStartSimulations; it != itEndSimulations; ++it) simulations.col(it-itStartSimulations) = Map<const VectorXd>(it->GetArray()+1, it->GetNbinsX());
  
}

