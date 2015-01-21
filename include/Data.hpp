#ifndef DATA_H
#define DATA_H

#include <vector>
#include <ostream>
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

class Data{

friend class DataRebin;
friend Data join(const Data& d1, const Data& d2);//utilitary function to join Data objects and return the joined data

  std::vector<boost::filesystem::path> filepaths;
  std::vector<Hist> h;
  std::vector<Eigen::MatrixXd> matrices;
  void FillMatricesFromRoot(std::vector<TMatrixD>::const_iterator itStart, std::vector<TMatrixD>::const_iterator itEnd);//convert TMatrices to fill Eigen matrices
  void completeWithEmptyMatrices();//if matrices is not the same size as h, fill up the matrices std::vector with empty matrices
  void PathGrabber(const boost::filesystem::path& search_path);
  void PathGrabber(const boost::filesystem::path& search_path, const string& file_sorter);
  static bool PathMatches(const boost::filesystem::directory_iterator& it, const string& file_sorter);
  void FillIfRoot(const boost::filesystem::directory_iterator& it);
  void FillIfRootAnd(const boost::filesystem::directory_iterator& it, const string& file_sorter);
  static bool ItemMatches(TObject* obj, const char* className);
  void StoreAsHist(TObject* obj);
  void StoreAsMatrix(TObject* obj);
  void StoreFromCan(TObject* read_object);
  void StoreData();

public:
  Data();
  Data(const std::vector<TH1D>& h);
  Data(const std::vector<TH1D>& h, const std::vector<TMatrixD>& matrices);
  Data(const boost::filesystem::path& search_path);
  Data(const boost::filesystem::path& search_path, const string& file_sorter);//fills any root file whose name matches the string sorter
  Data(const Data& other);
  const Data& operator=(const Data& other);
  const std::vector<boost::filesystem::path>& getFilePaths() const;
  const std::vector<Hist>& getHistograms() const;
  const std::vector<TH1D> getTH1DCopies() const;
  const std::vector<Eigen::MatrixXd>& getMatrices() const;
  TMatrixD getRootMatrixCopy(unsigned i) const;
  unsigned getSize() const;//returns the largest size of the std::vectors in Data

};

std::ostream& operator<<(std::ostream& output, const Data& data);

#endif