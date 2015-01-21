#ifndef EXCLUSION_H
#define EXCLUSION_H

#include <thread>
#include <TGraph.h>
#include <TLatex.h>
#include "Binning.hpp"
#include "TimeEstimator.hpp"

class Exclusion{
  
  unsigned nSigma;//desired significance
  Binning heFraction;//Binning to explore for the 8He fraction
  std::vector<double> time;//relative time increase to exclude the corresponding fraction
  TGraph graph;//graph to store the required time vs heBins
  
public:  
  Exclusion(unsigned nSigma, const Binning& heFraction);
  void buildExclusionGraph(const Data& dataToFit, const Data& simulations);
  void makeUpGraph(unsigned colourNumber = 4);//sets the graph to an exclusion plot style of colour = colourNumber (defaults to blue)
  void setSignificance(unsigned nSigma);//sets nSigma
  void setBinning(const Binning& heFraction);
  const Binning& getHeFraction() const;
  const std::vector<double>& getTime() const;
  const TGraph& getExclusionGraph() const;
  
};

std::ostream& operator<<(std::ostream& output, const Exclusion& exclusion);

#endif
