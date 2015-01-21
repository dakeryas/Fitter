#include "DataRebin.hpp"
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

  Minimizer min(ROOT::Math::Functor(chiSquared, simulations.getHistograms().size()));
  min.Process();
  cout<<min<<"\n";
  cout<<"NDF = "<<dataToFit.getHistograms().front().GetNbinsX()-2<<endl;
  
}

void Fitter(const boost::filesystem::path& directory, const string& data_sorter, const string& simu_sorter){

  Data measures = Data(directory, data_sorter);//retrieve the Data Histogram
  Data simu = Data(directory, simu_sorter);//retrieve the simulations Histograms
  
  DataRebin rebinner(join(measures, simu));//join measures and simulations to get the right rebin
  DataRebin measuresRebinner(measures, rebinner.GetRebin());//define a rebinner with the newly found Binning
  measuresRebinner.ApplyRebin();//apply this Binning to the measures
  DataRebin simuRebinner(simu, rebinner.GetRebin());
  simuRebinner.ApplyRebin();//apply this Binning to the simulations

  fitFirstToRest(measuresRebinner.GetData(), simuRebinner.GetData());
  Binning heFracBinning = {5, 2.8, 20};//steps per percent, starting percent, ending percent
  saveExclusion(measuresRebinner.GetData(), simuRebinner.GetData(), heFracBinning, "helium_exclusion.root");//number of steps per percent first, then min frac to test, then last frac to test
  
}

int main (int argc, char* argv[]){
  
  if (argc == 3 && boost::filesystem::is_directory(boost::filesystem::path("./ToFit"))) Fitter(boost::filesystem::path("./ToFit"), argv[1], argv[2]);
  else if (argc == 4 && boost::filesystem::is_directory(boost::filesystem::path(argv[1]))) Fitter(boost::filesystem::path(argv[1]), argv[2], argv[3]);
  else cout<<"Error: you must provide a valid target directory, a data file to fit, and simulation files"<<endl;
  return 0;
  
}
