#ifndef HIST_H
#define HIST_H

#include <vector>
#include <iostream>
#include <TH1.h>

using namespace::std;

class Hist: public TH1D{ 

  vector<double> bins;
  void FillBins();

public:
  Hist();
  Hist(const TH1D& h);
  Hist(const Hist& other);
  const vector<double>& GetBins() const;

};

#endif