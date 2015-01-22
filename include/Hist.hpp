#ifndef HIST_H
#define HIST_H

#include <vector>
#include <iostream>
#include <TH1.h>

class Hist: public TH1D{ 

  std::vector<double> bins;
  void FillBins();

public:
  Hist();
  Hist(const TH1D& h);
  Hist(const Hist& other);
  const std::vector<double>& GetBins() const;

};

#endif