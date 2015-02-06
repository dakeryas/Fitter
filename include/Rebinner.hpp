#ifndef HREBIN_H
#define HREBIN_H

#include <cmath>
#include <algorithm>
#include "Data.hpp"

class Rebinner{ 
  
  std::vector<double> edge;//to store the bin edges result
  static void Normalise(Hist& h, const int ref);//normalise the bin contents of h according to the BinWidth of bin number 'normalise_bin_ref'
  static std::vector<double> getCommonElements(const std::vector<double>& v1, const std::vector<double>& v2, double epsilon = 1e-3);
  static std::vector<unsigned> getCommonIndices(const std::vector<double>& v1, const std::vector<double>& v2, double epsilon = 1e-3);//find the indices in v1 whose values match that of v2

public:
  Rebinner(const Data& data);//calls buildFrom(data)
  Rebinner(const std::vector<double>& edge);
  void buildFrom(const Data& data);//compute the new 'rebin' compatible with 'data' and the current rebin
  void excludeBinsAbove(const double& newUpEdge);//removes bins above newUpEdge in 'edge'
  void excludeBinsBelow(const double& newLowEdge);//removes bins below newLowEdge in 'edge'
  void squeezeBinning(unsigned factor);
  void rebin(Data& data) const;//actually rebin 'data'
  bool admissibleRebinFor(Data& data) const;//check whether 'edge' is admissible or not for 'data'
  const std::vector<double>& getRebin() const;
  unsigned getNumberOfBins() const;

};

std::ostream& operator<<(std::ostream& output, const Rebinner& rebinner);

#endif