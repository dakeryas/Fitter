#ifndef EXCLUSION_H
#define EXCLUSION_H

#include <TGraph.h>
#include <TLatex.h>
#include "Binning.hpp"
#include "TimeEstimator.hpp"

class Exclusion{
  
  unsigned nSigma;//desired significance
  Binning heFraction;//Binning to explore for the 8He fraction
  TGraph graph;//graph to store the required time vs heBins
  
public:  
  Exclusion(unsigned nSigma, const Binning& heFraction);
  void buildExclusionGraph(const Data& dataToFit, const Data& simulations);
  void makeUpGraph(unsigned colourNumber = 4);//sets the graph to an exclusion plot style of colour = colourNumber (defaults to blue)
  const TGraph& getExclusionGraph() const;
  
};

#endif
