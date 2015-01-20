#include <TGraph.h>
#include <TLatex.h>
#include "DataRebin.hpp"
#include "Chi.hpp"
#include "Minimizer.hpp"

struct binning{
  
  unsigned steps;//number of steps per unit (not equal to the number of steps between min and max)
  double min;
  double max;
  
};

TGraph exclusion(const Data& dataToFit, const Data& simulations, const unsigned nSigma, const binning& heBins){//modify h to fill the first component with fake data taken from the simulations

  VectorXd fractions(simulations.GetSize());//we have one simulation less the total number of histograms
  Chi chiSquared(dataToFit, simulations);//intialisation with the real values
  VectorXd dataError = chiSquared.GetDataErr();//save the data error
  Minimizer min(ROOT::Math::Functor(chiSquared, simulations.GetSize()));
  
  unsigned rangeFrac = floor(heBins.max - heBins.min);//rangeFrac will actually stand for the fraction range
  rangeFrac *= heBins.steps;//we as many points as we have steps per percent times the total number of percents
  rangeFrac += 1;
  double HeFraction[rangeFrac];
  for(unsigned k = 0; k<rangeFrac; ++k) HeFraction[k] =  static_cast<double>(k+heBins.min*heBins.steps)/100/heBins.steps; //divide each percent into four points
  double time[rangeFrac];
  double dataFactor = 1.01;
  
  for(unsigned k = 0; k<rangeFrac; ++k){
    
    fractions(0) = HeFraction[k];
    fractions(1) = 1 - fractions(0);
    
    chiSquared.SetData(chiSquared.GetSimulations()*fractions);
    chiSquared.SetDataErr(dataError);//reset the data error for every new fraction to test
    min.Update(ROOT::Math::Functor(chiSquared,simulations.GetSize()));
    min.Process();
    
    double dataIncrease = 0;

    while(min.GetSol().front()-nSigma*min.GetErrors().front()<0){//test zero-ness to nSigma

      chiSquared.SetDataErr(chiSquared.GetDataErr()/sqrt(dataFactor));//double the data if the errors are to large to exclude a zero fraction of the first element
      min.Update(ROOT::Math::Functor(chiSquared,simulations.GetSize()));
      min.Process();
      ++dataIncrease;
      
    }
    
    time[k] = pow(dataFactor, dataIncrease);
    cout<<"For "<<100*HeFraction[k]<<"% of Helium, we must have "<<time[k]<<" time(s) as much data to exclude a 0% fraction of Helium with a "<<nSigma<<"-sigma confidence"<<endl;

  }
  
  return TGraph(rangeFrac, time, HeFraction);
  
}

void makeUpExclusion(TGraph& graph, const unsigned colourNumber){
  
  graph.SetTitle("Helium exclusion zones");
  graph.GetXaxis()->SetTitle("Data (times current amount)");
  graph.GetXaxis()->SetTitleOffset(1.25);
  graph.GetYaxis()->SetTitle("Actual Helium fraction");
  graph.GetYaxis()->SetTitleOffset(1.25);
  graph.SetFillColor(colourNumber);
  graph.SetFillStyle(3004);
  graph.SetLineColor(colourNumber);
  graph.SetLineWidth(-602);//write an number greater than 99 otherwise the exclusion area is not drawn ! //add a minus sign to reverse the direction of the dashes !
  
}

void saveExclusion(const Data& dataToFit, const Data& simulations, const binning& heBins, const char* outname){
 
  TCanvas can("can");
  can.SetGrid();
  can.SetLogx();
  
  TGraph exclusionGraphTwo = exclusion(dataToFit, simulations, 2, heBins);
  TGraph exclusionGraphThree = exclusion(dataToFit, simulations, 3, heBins);
  TGraph exclusionGraphFour = exclusion(dataToFit, simulations, 4, heBins);
  makeUpExclusion(exclusionGraphTwo, 4);
  makeUpExclusion(exclusionGraphThree, 2);
  makeUpExclusion(exclusionGraphFour, 416+2);

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

  Minimizer min(ROOT::Math::Functor(chiSquared, simulations.GetHistograms().size()));
  min.Process();
  cout<<min<<endl;
  cout<<"NDF = "<<dataToFit.GetHistograms().front().GetNbinsX()-2<<endl;
  
}

void Fitter(const path& directory, const string& data_sorter, const string& simu_sorter){

  Data measures = Data(directory, data_sorter);//retrieve the Data Histogram
  Data simu = Data(directory, simu_sorter);//retrieve the simulations Histograms
  
  DataRebin rebinner(join(measures, simu));//join measures and simulations to get the right rebin
  DataRebin measuresRebinner(measures, rebinner.GetRebin());//define a rebinner with the newly found binning
  measuresRebinner.ApplyRebin();//apply this binning to the measures
  DataRebin simuRebinner(simu, rebinner.GetRebin());
  simuRebinner.ApplyRebin();//apply this binning to the simulations

  fitFirstToRest(measuresRebinner.GetData(), simuRebinner.GetData());
  binning heFracBinning = {5, 3, 20};//steps per percent, starting percent, ending percent
  saveExclusion(measuresRebinner.GetData(), simuRebinner.GetData(), heFracBinning, "helium_exclusion.root");//number of steps per percent first, then min frac to test, then last frac to test
  
}

int main (int argc, char* argv[]){
  
  if (argc == 3 && is_directory(path("./ToFit"))) Fitter(path("./ToFit"), argv[1], argv[2]);
  else if (argc == 4 && is_directory(path(argv[1]))) Fitter(path(argv[1]), argv[2], argv[3]);
  else cout<<"Error: you must provide a valid target directory, a Datasfile to fit, and simulation files"<<endl;
  return 0;
  
}
