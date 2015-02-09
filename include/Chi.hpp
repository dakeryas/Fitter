#ifndef CHI_H
#define CHI_H

#include <Eigen/LU>
#include "Data.hpp"

class Chi{
  Eigen::VectorXd dataToFit;
  Eigen::VectorXd dataErrors;
  Eigen::MatrixXd simulations;//each column stands for a spectrum (each to be weighted by the corresponding fraction)
  std::vector<Eigen::MatrixXd> covariances;//covariances associated to the simulations, the sizes must match
  static int getMaxNbins(std::vector<TH1D>::const_iterator itStartSimulations, std::vector<TH1D>::const_iterator itEndSimulations);
  void fillErrors(const TH1D& dataToFit);//fill dataErrors from the TH1D:getBinError method
  void fillSim(std::vector<TH1D>::const_iterator itStartSimulations, std::vector<TH1D>::const_iterator itEndSimulations);//fills the simulations matrix from the simulations Histogram

public:
  Chi();
  Chi(const Eigen::VectorXd& dataToFit, const Eigen::MatrixXd& simulations);//if no errors are given, build the dataErrors from the coefficient-wise sqrt of dataToFit
  Chi(const Eigen::VectorXd& dataToFit, const Eigen::VectorXd& dataErrors, const Eigen::MatrixXd& simulations);
  Chi(const Eigen::VectorXd& dataToFit, const Eigen::VectorXd& dataErrors, const Eigen::MatrixXd& simulations, const std::vector<Eigen::MatrixXd>& covariances);
  Chi(std::vector<TH1D>::const_iterator itDataToFit, std::vector<TH1D>::const_iterator itStartSimulations, std::vector<TH1D>::const_iterator itEndSimulations);
  Chi(const TH1D& dataToFit, const std::vector<TH1D>& simulations);
  Chi(const Data& dataToFit, const Data& simulations);//retrive the Histograms from the Data class to use the first constructor
  double operator ()(const double* args) const;
  void setData(const Eigen::VectorXd& dataToFit);//keep in mind that it only sets the data itself, without its errors (which keep their former values; they can be set with the setDataErr method)
  void setDataErr(const Eigen::VectorXd& dataErrors);
  const Eigen::VectorXd& getDataToFit() const;
  const Eigen::VectorXd& getDataErr() const;
  const Eigen::MatrixXd& getSimulations() const;
  unsigned getNumberOfFreeParameters() const;
  Eigen::VectorXd getRelativeDispersion(const std::vector<double>& args) const;//return the vector (data - best fit model)/-statistical error only)
  
};

#endif
