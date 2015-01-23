#ifndef HIST_H
#define HIST_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <TH1.h>

class Hist: public TH1D{ 

  std::vector<double> edge;
  void fillEdge();

public:
  Hist();
  Hist(const TH1D& h);
  Hist(const Hist& other);
  Hist& operator+=(const Hist& other);
  Hist& operator*=(double a);
  unsigned getNumberOfBins() const;
  const std::vector<double>& getEdge() const;
  bool isCompatibleWith(const Hist& other) const; //checks if the binnings match 
  
};

Hist operator+(Hist h1, const Hist& h2);
std::ostream& operator<<(std::ostream& output, const Hist& hist);

#endif