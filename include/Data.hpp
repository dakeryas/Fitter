#ifndef DATA_H
#define DATA_H

#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <boost/filesystem.hpp>
#include <TH1.h>
#include <TMatrixF.h>
#include <TMatrixD.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TObject.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TIterator.h>
#include "Hist.hpp"

class DataRebin;
using namespace::std;
using namespace::Eigen;
using namespace boost::filesystem;

class Data{

friend class DataRebin;
friend Data join(const Data& d1, const Data& d2);//utilitary function to join Data objects and return the joined data

  vector<path> filepaths;
  vector<Hist> h;
  vector<MatrixXd> matrices;
  void FillMatricesFromRoot(vector<TMatrixD>::const_iterator itStart, vector<TMatrixD>::const_iterator itEnd);//convert TMatrices to fill Eigen matrices
  void completeWithEmptyMatrices();//if matrices is not the same size as h, fill up the matrices vector with empty matrices
  void PathGrabber(const path& search_path);
  void PathGrabber(const path& search_path, const string& file_sorter);
  static bool PathMatches(const directory_iterator& it, const string& file_sorter);
  void FillIfRoot(const directory_iterator& it);
  void FillIfRootAnd(const directory_iterator& it, const string& file_sorter);
  static bool ItemMatches(TObject* obj, const char* className);
  void StoreAsHist(TObject* obj);
  void StoreAsMatrix(TObject* obj);
  void StoreFromCan(TObject* read_object);
  void StoreData();

public:
  Data();
  Data(const vector<TH1D>& h);
  Data(const vector<TH1D>& h, const vector<TMatrixD>& matrices);
  Data(const path& search_path);
  Data(const path& search_path, const string& file_sorter);//fills any root file whose name matches the string sorter
  Data(const Data& other);
  const Data& operator=(const Data& other);
  const vector<path>& GetFilePaths() const;
  const vector<Hist>& GetHistograms() const;
  const vector<TH1D> GetTH1DCopies() const;
  const vector<MatrixXd>& GetMatrices() const;
  TMatrixD GetRootMatrixCopy(unsigned i) const;
  unsigned GetSize() const;//returns the largest size of the vectors in Data

};

ostream& operator<<(ostream& output, const Data& data);

#endif