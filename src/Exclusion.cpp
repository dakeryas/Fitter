#include "Exclusion.hpp"

using namespace std;

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
  Chi chiSquared;
  vector<double> initialValues;
  unsigned nSigma;
  
  relativeKthTimeEstimator(double& kthTime, const Chi& chiSquared, const vector<double>& initialValues, unsigned nSigma):kthTime(kthTime),chiSquared(chiSquared),initialValues(initialValues),nSigma(nSigma){
    
  };
  void operator()(){
    
    Minimizer min(ROOT::Math::Functor(chiSquared, chiSquared.getNumberOfFreeParameters()));
    TimeEstimator timeEstimator;
    kthTime = timeEstimator.getRelativeTime(min, chiSquared, initialValues, nSigma);
    
  };
  
};

ostream& operator<<(ostream& output, const Exclusion& exclusion){
  
  for(unsigned k = 0; k<exclusion.getHeFraction().getNumberOfSteps(); ++k){
    
    output<<"To exclude the "<<exclusion.getHeFraction().getValue(k)<<" fraction to "<<exclusion.getSignificance()<<"-sigma we need "<<exclusion.getTime().at(k)<<" time(s) as much data as we currently have.\n";
    
  }
  
  return output;

}

Exclusion::Exclusion(unsigned int nSigma, const Binning& heFraction):nSigma(nSigma),heFraction(heFraction),time(heFraction.getNumberOfSteps()){
  
}

void Exclusion::buildExclusionGraph(const Data& dataToFit, const Data& simulations){//modify h to fill the first component with fake data taken from the simulations

  VectorXd fractions(simulations.getSize());//we have one simulation less the total number of histograms
  Chi chiSquared(dataToFit, simulations);//intialisation with the real values
  VectorXd dataError = chiSquared.getDataErr();//save the data error

  vector<relativeKthTimeEstimator> workers;
  
  for(unsigned k = 0; k<heFraction.getNumberOfSteps(); ++k){
    
    fractions(0) = heFraction.getValue(k);
    fractions(1) = 1 - fractions(0);
    
    chiSquared.SetData(chiSquared.getSimulations()*fractions);
    chiSquared.SetDataErr(dataError);//reset the data error for every new fraction to test
    
    workers.push_back(relativeKthTimeEstimator(time.at(k), chiSquared, {fractions(0),fractions(1)}, nSigma)); //the initial values for the chiSquared are the fractions you put in 

  }
  
  launchThreads(workers);
  cout<<*this<<endl;
  
  graph = TGraph(heFraction.getNumberOfSteps(), time.data(), heFraction.getDataBins());
  
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

void Exclusion::setSignificance(unsigned int nSigma){
  
  this->nSigma = nSigma;

}

void Exclusion::setBinning(const Binning& heFraction){
  
  this->heFraction = heFraction;
  time.resize(heFraction.getNumberOfSteps(), 0);

}

const unsigned int& Exclusion::getSignificance() const{
  
  return nSigma;

}

const Binning& Exclusion::getHeFraction() const{
  
  return heFraction;

}

const vector<double>& Exclusion::getTime() const{
  
  return time;

}

const TGraph& Exclusion::getExclusionGraph() const{
  
  return graph;

}
