#ifndef PATHGRABBER_H
#define PATHGRABBER_H

#include <ostream>
#include <boost/filesystem.hpp>

class PathGrabber{

  std::vector<boost::filesystem::path> filePaths;//paths of the ROOT files where to extract histograms and matrices
  static bool pathMatches(const boost::filesystem::directory_iterator& it, const std::string& fileSorter);
  void pushIfRoot(const boost::filesystem::directory_iterator& it);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file
  void pushIfRootAnd(const boost::filesystem::directory_iterator& it, const std::string& fileSorter);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file

public:
  void pushPathsFrom(const boost::filesystem::path& searchPath);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file
  void pushPathsFrom(const boost::filesystem::path& searchPath, const std::string& fileSorter);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file
  void clear();//resizes the filePaths to zero
  const std::vector<boost::filesystem::path>& getFilePaths() const;
  
};

std::ostream& operator<<(std::ostream& output, const PathGrabber& pathGrabber);

#endif