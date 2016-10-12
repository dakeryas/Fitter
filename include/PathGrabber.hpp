#ifndef PATHGRABBER_H
#define PATHGRABBER_H

#include <iostream>
#include <boost/filesystem.hpp>

class PathGrabber{

  std::vector<boost::filesystem::path> filePaths;//paths of the ROOT files where to extract histograms and matrices
  std::string fileSorter;
  std::string extension;//extension the matched paths should have
  bool pathMatches(boost::filesystem::directory_iterator it);//check if the path matches 'fileSorter'

public:
  PathGrabber() = default;
  PathGrabber(const std::string& extension);//the extension the files should have
  const std::string& getMatchingExtension() const;
  void setMatchingExtension(const std::string& extension);
  void pushPathsFrom(const boost::filesystem::path& searchPath);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file
  void pushPathsFrom(const boost::filesystem::path& searchPath, const std::string& fileSorter);//pushes back into 'filePaths' the path pointed to by 'it' if that path references a ROOT file
  void clearPaths();//resizes the filePaths to zero
  const std::vector<boost::filesystem::path>& getFilePaths() const;
  bool foundPaths() const;
  
};

std::ostream& operator<<(std::ostream& output, const PathGrabber& pathGrabber);

#endif