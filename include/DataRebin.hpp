#ifndef HREBIN_H
#define HREBIN_H

#include <cmath>
#include <algorithm>
#include "Data.hpp"

class DataRebin{ 

friend std::ostream& operator<<(std::ostream& output, const DataRebin& rebinner);
  
  Data data;//local copy of the data to rebin
  std::vector<double> rebin;//to store the rebin result
  double* rebinArray;//rebin array in a C-style version
  static std::vector<double> GetCommonElements(const std::vector<double>& v1, const std::vector<double>& v2, const double epsilon = 1e-3);
  static vector<int> GetCommonIndices(const std::vector<double>& v1, const std::vector<double>& v2, const double epsilon = 1e-3);//find the indices in v1 whose values match that of v2
  static void Normalise(Hist& h, const int ref);//normalise the bin contents of h according to the BinWidth of bin number 'normalise_bin_ref'
  void AllocAndInitArray();//allocate the rebinArray from the size of the vector rebin and initilase rebinArray from the contents of rebin
  void FillRebin();//compute the rebin vector(and fill the array)
  bool AdmissibleRebin() const;//check whether the rebin is admissible or not

public:
  DataRebin();
  DataRebin(const Data& data);
  DataRebin(const Data& data, const std::vector<double>& rebin);//pass the rebin vector and check if it is admissible, otherwise, overlook it
  DataRebin(const DataRebin& other);
  DataRebin& operator = (const DataRebin& other);
  ~DataRebin();
  void ApplyRebin();//actually rebin the Histograms
  const Data& GetData() const;
  const std::vector<double>& GetRebin() const;

};

#endif