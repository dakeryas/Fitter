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
  Hist& operator+=(const Hist& other);
  Hist& operator*=(double a);
  unsigned getNumberOfBins() const;
  const std::vector<double>& getBins() const;
  bool isCompatibleWith(const Hist& other) const; //checks if the binnings match 
  
};

Hist operator+(Hist h1, const Hist& h2);

#endif