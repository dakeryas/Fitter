#ifndef BINNING_H
#define BINNING_H

#include <iostream>
#include <vector>
#include <cmath>

class Binning{
  
  unsigned numberOfStepsPerPercent;//number of steps per percent (not equal to the number of steps between min and max)
  double minPercent;//starting percent
  double maxPercent;//ending percent
  std::vector<double> bins;
  void fillBins();//fill bins from numberOfStepsPerPercent, minPercent, and maxPercent
  
public:  
  Binning(unsigned numberOfStepsPerPercent, double minPercent, double maxPercent);
  const unsigned& getNumberOfStepsPerPercent() const ;
  const double& getMinPercent() const;
  const double& getMaxPercent() const;
  double getRange() const;
  unsigned getFlooredRange() const; //uses the floor function
  unsigned getNumberOfSteps() const;
  double getValue(unsigned i) const;//returns -1 if the value is not in the range
  const double* getDataBins() const;
  void setBinning(unsigned numberOfStepsPerPercent, double minPercent, double maxPercent);
  
};

std::ostream& operator <<(std::ostream& output, const Binning& binning);

#endif


