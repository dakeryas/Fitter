#ifndef CHI_H
#define CHI_H

#include <Eigen/LU>
#include "Data.hpp"

using namespace::Eigen;

class Chi{
  VectorXd dataToFit;
  VectorXd dataErrors;
  MatrixXd simulations;//each column stands for a spectrum (each to be weighted by the corresponding fraction)
  vector<MatrixXd> covariances;//covariances associated to the simulations, the sizes must match
  static int GetMaxNbins(vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations);
  void FillErrors(const TH1D& dataToFit);//fill dataErrors from the TH1D:GetBinError method
  void FillSim(vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations);//fills the simulations matrix from the simulations Histogram

public:
  Chi();
  Chi(const VectorXd& dataToFit, const MatrixXd& simulations);//if no errors are given, build the dataErrors from the coefficient-wise sqrt of dataToFit
  Chi(const VectorXd& dataToFit, const VectorXd& dataErrors, const MatrixXd& simulations);
  Chi(const VectorXd& dataToFit, const VectorXd& dataErrors, const MatrixXd& simulations, const vector<MatrixXd>& covariances);
  Chi(vector<TH1D>::const_iterator itDataToFit, vector<TH1D>::const_iterator itStartSimulations, vector<TH1D>::const_iterator itEndSimulations);
  Chi(const TH1D& dataToFit, const vector<TH1D>& simulations);
  Chi(const Data& dataToFit, const Data& simulations);//retrive the Histograms from the Data class to use the first constructor
  double operator ()(const double* args) const;
  void SetData(const VectorXd& dataToFit);
  void SetDataErr(const VectorXd& dataErrors);
  const VectorXd& GetDataErr() const;
  const MatrixXd& GetSimulations() const;

};

#endif
