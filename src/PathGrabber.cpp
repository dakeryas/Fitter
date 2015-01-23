#include "PathGrabber.hpp"

using namespace std;
using namespace boost::filesystem;

ostream& operator<<(ostream& output, const PathGrabber& pathGrabber){
  
  output<<"Grabbed the files:\n";
  for(const path& filep : pathGrabber.getFilePaths()) output<<filep<<"\n";
  return output;  
  
}

bool PathGrabber::pathMatches(const directory_iterator& it, const string& fileSorter){
  
  return it->path().filename().string().find(fileSorter) != string::npos;
  
}

void PathGrabber::pushIfRoot(const directory_iterator& it){//fills 'found_paths' with the path pointed to by 'it' if 'it' refers to a file
  
  if(is_directory(it->status()) == 0 && pathMatches(it, ".root")) filePaths.push_back(it->path()); ; //find returns the position of 'name' in the path, so if it returned a position indeed(namely not npos), we found something
  
}

void PathGrabber::pushIfRootAnd(const directory_iterator& it, const string& fileSorter){//fills 'found_paths' with the path pointed to by 'it' if 'it' refers to a file
  
  if(is_directory(it->status()) == 0 && pathMatches(it, ".root") && pathMatches(it, fileSorter)) filePaths.push_back(it->path()); ; //find returns the position of 'name' in the path, so if it returned a position indeed(namely not npos), we found something
  
}

void PathGrabber::pushPathsFrom(const path& searchPath){ //retrieves the path of all files in a directory

  directory_iterator end; //the default constructor creates an end iterator which cannot be reached unless nothing was found
  for(directory_iterator it(searchPath); it!= end; ++it) pushIfRoot(it);
  
}

void PathGrabber::pushPathsFrom(const path& searchPath, const string& fileSorter){ //retrieves the path of all the files in a directory that match the name fileSorter

  directory_iterator end; //the default constructor creates an end iterator which cannot be reached unless nothing was found
  for(directory_iterator it(searchPath); it!= end; ++it) pushIfRootAnd(it, fileSorter);
  
}

void PathGrabber::clear(){
  
  filePaths.clear();

}

const vector<path>& PathGrabber::getFilePaths() const{
  
  return filePaths;

}

