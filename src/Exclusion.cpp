#include "Exclusion.hpp"

template <class T> void launchThreads(vector<T>& workers){

  unsigned nThreads = thread::hardware_concurrency();//get the number of working cores
  if(nThreads == 0) nThreads = 1;
  vector<thread> threads(nThreads);
  
  auto itWk = workers.begin();
  while(itWk != workers.end()){
    
    auto itTh = threads.begin(); 
    while(itTh != threads.end() && itWk != workers.end()){

      *itTh = thread(*itWk);//start each thread
      ++itTh;
      ++itWk;
      
    }
    
    for(thread& th : threads) if(th.joinable()) th.join();//join them all to the current thread
    
  }
  
}

struct relativeKthTimeEstimator{
  
  double& kthTime;
  ROOT::Math::Functor f;
  Chi chiSquared;
  unsigned nSigma;
  
  relativeKthTimeEstimator(double& kthTime, const ROOT::Math::Functor& f, const Chi& chiSquared, unsigned nSigma):kthTime(kthTime),f(f),chiSquared(chiSquared),nSigma(nSigma){
    
  };
  void operator()(){
    
    Minimizer min(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));
    TimeEstimator timeEstimator;
    kthTime = timeEstimator.getRelativeTime(min, chiSquared, nSigma);
    
  };
  
};

Exclusion::Exclusion(unsigned int nSigma, const Binning& heFraction):nSigma(nSigma),heFraction(heFraction){
  
}

void Exclusion::buildExclusionGraph(const Data& dataToFit, const Data& simulations){//modify h to fill the first component with fake data taken from the simulations

  VectorXd fractions(simulations.getSize());//we have one simulation less the total number of histograms
  Chi chiSquared(dataToFit, simulations);//intialisation with the real values
  VectorXd dataError = chiSquared.getDataErr();//save the data error
  Minimizer min(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));

  double time[heFraction.getNumberOfSteps()];

  vector<relativeKthTimeEstimator> workers;
  
  for(unsigned k = 0; k<heFraction.getNumberOfSteps(); ++k){
    
    fractions(0) = heFraction.getValue(k);
    fractions(1) = 1 - fractions(0);
    
    chiSquared.SetData(chiSquared.getSimulations()*fractions);
    chiSquared.SetDataErr(dataError);//reset the data error for every new fraction to test
    
    workers.push_back(relativeKthTimeEstimator(time[k],min.getFunctor(), chiSquared, nSigma));

  }
  
  launchThreads(workers);
  
  graph = TGraph(heFraction.getNumberOfSteps(), time, heFraction.getDataBins());
  
}

void Exclusion::makeUpGraph(unsigned colourNumber){
  
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

const TGraph& Exclusion::getExclusionGraph() const{
  
  return graph;

}
