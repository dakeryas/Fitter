#include "PathGrabber.hpp"

using namespace boost::filesystem;

std::ostream& operator<<(std::ostream& output, const PathGrabber& pathGrabber){
  
  output<<"Grabbed the files:\n";
  for(const path& filep : pathGrabber.getFilePaths()) output<<filep<<"\n";
  return output;  
  
}

PathGrabber::PathGrabber(const std::string& extension):extension(extension){

}

const std::string& PathGrabber::getMatchingExtension() const{
  
  return extension;

}

void PathGrabber::setMatchingExtension(const std::string& extension){
  
  this->extension = extension;

}

bool PathGrabber::pathMatches(directory_iterator it){

  return it->path().filename().string().find(fileSorter) != std::string::npos;
  
}

void PathGrabber::pushPathsFrom(const path& searchPath){ //retrieves the path of all files in a directory matching 'extension'

  pushPathsFrom(searchPath, "");
  
}

void PathGrabber::pushPathsFrom(const path& searchPath, const std::string& fileSorter){ //retrieves the path of all the files in a directory that match the name fileSorter

  this->fileSorter = fileSorter;
  
  directory_iterator end; //the default constructor creates an end iterator which cannot be reached unless nothing was found
  for(directory_iterator it(searchPath); it!= end; ++it){

    if(is_regular_file(it->status()) && (it->path().extension() == extension) && pathMatches(it)){
      
      filePaths.push_back(it->path());
      
    }
    
  }
  
  std::sort(filePaths.begin(), filePaths.end());//sort the paths
  
}

void PathGrabber::clearPaths(){
  
  filePaths.clear();

}

const std::vector<path>& PathGrabber::getFilePaths() const{
  
  return filePaths;

}

bool PathGrabber::foundPaths() const{
  
  return !filePaths.empty();

}

