#ifndef DATA_H
#define DATA_H

#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <TH1.h>
#include <TMatrixF.h>
#include <TMatrixD.h>
#include "Hist.hpp"

class Data{

  std::vector<Hist> histograms;
  std::vector<Eigen::MatrixXd> matrices;

public:
  Data();
  Data(const std::vector<TH1D>& histograms);
  Data(const std::vector<TH1D>& histograms, const std::vector<TMatrixD>& matrices);
  Data& operator+=(const Data& other);//add the bin contents and the matrices coefficient wise
  Data& operator*=(double a);//multiplies the histograms and the matrices by 'a'
  void completeWithEmptyMatrices();//if matrices is not the same size as h, fill up the matrices std::vector with empty matrices
  void pushHist(const Hist& hist);//pushes hist into 'histograms'
  void pushMatrix(const Eigen::MatrixXd& matrix);//pushes 'matrix' into matrices
  void push(const Hist& hist, const Eigen::MatrixXd& matrix);//pushes using pushHist and pushMatrix
  void fillMatricesFromRoot(std::vector<TMatrixD>::const_iterator itStart, std::vector<TMatrixD>::const_iterator itEnd);//convert TMatrices to fill Eigen matrices
  const std::vector<Hist>& getHistograms() const;
  std::vector<Hist>::iterator getHistogramStartIterator();//awfully ugly solution
  const std::vector<TH1D> getTH1DCopies() const;
  const std::vector<Eigen::MatrixXd>& getMatrices() const;
  std::vector<Eigen::MatrixXd>::iterator getmatricesStartIterator();//awfully ugly solution
  TMatrixD getRootMatrixCopy(unsigned i) const;
  unsigned getNumberOfBins() const;//returns the minimum of the Hist::getNumberOfBins
  unsigned getSize() const;//returns the largest size of the std::vectors in Data
  void clear();//resizes all vectors to zero

};

Data join(Data d1, const Data& d2);//utilitary function to join Data objects and return the joined data
std::ostream& operator<<(std::ostream& output, const Data& data);
Data operator+(Data d1, const Data& d2);
Data operator*(Data data, double a);
Data operator*(double a, Data data);

#endif