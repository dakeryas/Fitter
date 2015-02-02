#include "PathGrabber.hpp"
#include "Storer.hpp"
#include "Rebinner.hpp"
#include "Exclusion.hpp"

void saveExclusion(const Data& dataToFit, const Data& simulations, const Binning& heFraction, const char* outname){
 
  TCanvas can("can");
  can.SetGrid();
  can.SetLogx();
  
  Exclusion exclusionTwo(2, heFraction);//exclusion object to two sigma for the 8He Fractions to test
  Exclusion exclusionThree(3, heFraction);
  Exclusion exclusionFour(4, heFraction);
  
  exclusionTwo.buildExclusionGraph(dataToFit, simulations);
  exclusionThree.buildExclusionGraph(dataToFit, simulations);
  exclusionFour.buildExclusionGraph(dataToFit, simulations);
  
  exclusionTwo.makeUpGraph(4);
  exclusionThree.makeUpGraph(2);
  exclusionFour.makeUpGraph(416+2);
  
  TGraph exclusionGraphTwo = exclusionTwo.getExclusionGraph();
  TGraph exclusionGraphThree = exclusionThree.getExclusionGraph();
  TGraph exclusionGraphFour = exclusionFour.getExclusionGraph();

  exclusionGraphTwo.Draw("AL");//draw with axis and line to enable exclusions...
  exclusionGraphThree.Draw("same");
  exclusionGraphFour.Draw("same");
  TLatex textTwo(3, 0.03, "2 #sigma");//draw for 2 times as much data and a 0.025 He fraction
  textTwo.SetTextColor(4);
  textTwo.Draw();
  TLatex textThree(3, 0.067, "3 #sigma");
  textThree.SetTextColor(2);
  textThree.Draw();
  TLatex textFour(3, 0.105, "4 #sigma");
  textFour.SetTextColor(416+2);
  textFour.Draw();
  
  TFile file(outname, "recreate");
  can.Write();
  
}

void fitFirstToRest(const Data& dataToFit, const Data& simulations){

  Chi chiSquared(dataToFit, simulations); //Data first and a vector of simulations secondly (with the Data removed first)

  Minimiser min(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));
  min.setInitialValues({0, dataToFit.getHistograms().front().Integral()});//have it start close to the actual solution
  min.Process();
  
  std::cout<<min<<"\n";
  std::cout<<"NDF = "<<dataToFit.getNumberOfBins()-2<<"\n";
  const double heFraction = min.getSol().front()/min.getSol().back();
  const double heFractionErr = heFraction * (min.getErrors().front()/min.getSol().front() + min.getErrors().back()/min.getSol().back());//use relative errors, i.e. df = f (ds1/s1 + ds2/s2)
  std::cout<<"He fraction = "<<heFraction<<" +/- "<<heFractionErr<<std::endl;
  
}

void Fitter(const boost::filesystem::path& directory, const std::string& dataSorterGd, const std::string dataSorterH, const std::string& simuSorterGd, const std::string& simuSorterH){

  Data measuresGd, measuresH, simuGd, simuH;
  
  PathGrabber pathGrabber;
  pathGrabber.pushPathsFrom(directory, dataSorterGd);
  Storer storer(pathGrabber.getFilePaths());
  storer.fill(measuresGd);
  
  pathGrabber.clear();
  pathGrabber.pushPathsFrom(directory, dataSorterH);
  storer.setFilePaths(pathGrabber.getFilePaths());
  storer.fill(measuresH);
  
  pathGrabber.clear();
  pathGrabber.pushPathsFrom(directory, simuSorterGd);//retrieve the paths that the contain the ROOT simulation files
  storer.setFilePaths(pathGrabber.getFilePaths());
  storer.fill(simuGd);//fill the ROOT objects from the paths into 'simu'
  
  pathGrabber.clear();
  pathGrabber.pushPathsFrom(directory, simuSorterH);
  storer.setFilePaths(pathGrabber.getFilePaths());
  storer.fill(simuH);
  
  Data measures = .5*measuresGd + .5*measuresH;
  double fractionGd = 4.176542e-01;
  Data simu = fractionGd * simuGd + (1-fractionGd) * simuH;
  
  Rebinner rebinner(join(measures, simu));//join measures and simulations to get the right rebin
  rebinner.rebin(measures);
  rebinner.rebin(measuresGd);
  rebinner.rebin(measuresH);
  rebinner.rebin(simu);
  rebinner.rebin(simuGd);
  rebinner.rebin(simuH);
  
  fitFirstToRest(measuresGd, simuGd);
  fitFirstToRest(measuresH, simuH);
  fitFirstToRest(measures, simu);
  Binning heFracBinning = {5, 2, 20};//steps per percent, starting percent, ending percent
  saveExclusion(measures, simu, heFracBinning, "helium_exclusion.root");//number of steps per percent first, then min frac to test, then last frac to test
  
}

int main (int argc, char* argv[]){
  
  if (argc == 5 && boost::filesystem::is_directory(boost::filesystem::path("./ToFit"))) Fitter(boost::filesystem::path("./ToFit"), argv[1], argv[2], argv[3], argv[4]);
  else if (argc == 6 && boost::filesystem::is_directory(boost::filesystem::path(argv[1]))) Fitter(boost::filesystem::path(argv[1]), argv[2], argv[3], argv[4], argv[5]);
  else std::cout<<"Error: you must provide a valid target directory, a data file to fit, and simulation files for Gd and Hydrogen"<<std::endl;
  return 0;
  
}
