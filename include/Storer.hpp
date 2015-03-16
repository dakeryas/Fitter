#ifndef STORER_H
#define STORER_H

#include <ostream>
#include <boost/filesystem.hpp>
#include <TCanvas.h>
#include <TFile.h>
#include <TObject.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TIterator.h>
#include "Data.hpp"

class Storer{

  std::vector<boost::filesystem::path> filePaths;//paths of the ROOT files where to extract histograms and matrices
  static bool itemMatches(TObject* obj, const char* className);
  void pushAsHist(Data& data, TObject* obj) const;//pushes the TObject as a Hist into data
  void pushAsMatrix(Data& data, TObject* obj) const;//pushes the TObject as a Matrix into data
  void pushFromCan(Data& data, TObject* read_object) const;//pushes the contents of the TObject into data

public:
  Storer(const boost::filesystem::path& filePath);//in case you only need to store from one path
  Storer(const std::vector<boost::filesystem::path>& filePaths);
  void setFilePaths(const std::vector<boost::filesystem::path>& filePaths);
  const std::vector<boost::filesystem::path>& getFilePaths() const;
  void pushPath(const boost::filesystem::path& filePath);
  void clear();//empties filePaths
  void fill(Data& data) const;//store what can be found in filePaths into data

};

std::ostream& operator<<(std::ostream& output, const Storer& data);

#endif