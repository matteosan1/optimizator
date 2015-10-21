#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TMath.h"

#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <string>
#include <vector>
#include <array>
#include <iostream>

#include "parser.h"

static int binning=5000;

class Variable {
public:
  Variable(std::string n) { name = n; };
  ~Variable() {};

  void setVar(std::string t, float l, float h, std::string dir) {
    type = t;
    low = l;
    high = h;
    direction = dir;
  };
  void setLow(float l) { low = l; };
  void setHigh(float h) { high = h; };
  
  void setLimit(float cut) {
    if (getDir() == "<")
      setHigh(cut);
    
    if (getDir() == ">")
      setLow(cut);
  }
  
  std::string getCutString() {
    std::ostringstream sel;
    if (getDir() == "<")
      sel << getName() << "<"  << getHigh(); 

    if (getDir() == ">")
      sel << getName() << ">" << getLow(); 

    return sel.str();
  }
  
  std::string getName() { return name; };
  std::string getType() { return type; };
  double getLow() { return double(low); };
  double getHigh() { return double(high); };
  std::string getDir() { return direction; };

private:
  std::string name;
  std::string type;
  float low;
  float high;
  std::string direction;
};

bool is_file_exist(const char *fileName) {
  std::ifstream infile(fileName);
  return infile.good();
}

Double_t* ComputeIntegral(TH1F* h, Bool_t direction=true) {
  Int_t nbins  = h->GetNbinsX();

  Double_t* fIntegral = new Double_t[nbins + 2];
  Int_t ibin = 0; 
  fIntegral[ibin] = 0;

  if (direction) { // cut <
    for (Int_t binx=1; binx <= nbins; ++binx) {
      ++ibin;
      Double_t y = h->GetBinContent(binx);
      fIntegral[ibin] = fIntegral[ibin - 1] + y;
    }
    //   - Normalize integral to 1
    if (fIntegral[nbins] == 0) {
      Error("ComputeIntegral", "Integral = zero"); 
      return 0;
    }

    for (Int_t bin=1; bin <= nbins; ++bin)  
      fIntegral[bin] /= fIntegral[nbins];
    fIntegral[nbins+1] = h->GetEntries();
  } else {
    ibin = nbins;
    for (Int_t binx=nbins; binx >= 1; --binx) {
      --ibin;
      Double_t y = h->GetBinContent(binx);
      fIntegral[ibin] = fIntegral[ibin + 1] + y;
      //std::cout << "INTEGRAL:" << h->GetBinLowEdge(binx) << " " << y << " " << fIntegral[ibin] << std::endl;
    }

    //   - Normalize integral to 1
    if (fIntegral[0] == 0) {
      Error("ComputeIntegral", "Integral = zero"); 
      return 0;
    }

    for (Int_t bin=nbins; bin >=0; --bin)  {
      fIntegral[bin] /= fIntegral[0];
      //std::cout << "INTEGRAL:" << bin << " " << h->GetBinLowEdge(bin) << " " << fIntegral[bin] << std::endl;
    }
    fIntegral[0] = h->GetEntries();
  }
  
  return fIntegral;
}

Long64_t BinarySearch(Long64_t n, const double *array, double value) {

  for (int i=0; i<n; i++) { 
    if (array[i] < value) 
      return i;
  }
}

Int_t GetQuantiles(TH1F* h, Int_t nprobSum, Double_t *q, const Double_t *probSum, bool direction) {

  Int_t i, ibin;
  Double_t *prob = (Double_t*)probSum;
  Int_t nq = nprobSum;
  const Int_t nbins = h->GetXaxis()->GetNbins();

  if (direction) { // cut <
    double* fIntegral = ComputeIntegral(h);//h->ComputeIntegral();
    //double* fIntegral = h->GetIntegral();
    //fIntegral[nbins] = signalEntries;

    for (i = 0; i < nq; i++) {
      ibin = TMath::BinarySearch(nbins, fIntegral, prob[i]);
      q[i] = h->GetBinLowEdge(ibin+1);
      const Double_t dint = fIntegral[ibin+1] - fIntegral[ibin];
      if (dint > 0)  {
	q[i] += h->GetBinWidth(ibin+1)*(prob[i]-fIntegral[ibin])/dint;
      }
    }
  } else {
    //h->ComputeIntegral();
    double* fIntegral = ComputeIntegral(h, false);//h->GetIntegral();
    //double* fIntegral = h->GetIntegral();

    for (i = 0; i < nq; i++) {
      ibin = BinarySearch(nbins, fIntegral, prob[i]);
      q[i] = h->GetBinLowEdge(ibin-1);
      const Double_t dint = fIntegral[ibin] - fIntegral[ibin-1];
      if (dint > 0)  {
	q[i] += h->GetBinWidth(ibin-1)*(prob[i]-fIntegral[ibin])/dint;
      }
    }
  }
  
  if (!probSum) delete [] prob;
  return nq;

}

void getVariables(std::vector<Variable>& variables, std::map<std::string, std::string> values) {

  std::string buf = getString(values, "vars");
  std::vector<std::string> stringValues;
  boost::split(stringValues, buf, boost::is_any_of(","));
  
  for (unsigned int i=0; i<stringValues.size(); i+=5) {
    std::vector<std::string> strs;
    std::vector<float> limits;
    getVar(i, stringValues, limits, strs);
  
    Variable v(strs[0]);
    v.setVar(strs[1], limits[0], limits[1], strs[2]);
    variables.push_back(v);
  }
}
	   
std::pair<float, float> makeHistos(RooDataSet* sig, RooDataSet* bkg, Variable v, float target, std::string weightVar) {

  TH1F* hsig = new TH1F("hsig", "", binning, v.getLow(), v.getHigh());
  TH1F* hbkg = new TH1F("hbkg", "", binning, v.getLow(), v.getHigh());
  
  for (unsigned int i=0; i<sig->numEntries(); i++) {
    const RooArgSet* set = sig->get(i);
    RooRealVar* var = (RooRealVar*)set->find(v.getName().c_str());
    RooRealVar* w = (RooRealVar*)set->find(weightVar.c_str());
    if (w != 0)
      hsig->Fill(var->getVal(), w->getVal());
    else
      hsig->Fill(var->getVal());
  }

  for (unsigned int i=0; i<bkg->numEntries(); i++) {
    const RooArgSet* set = bkg->get(i);
    RooRealVar* var = (RooRealVar*)set->find(v.getName().c_str());
    RooRealVar* w = (RooRealVar*)set->find(weightVar.c_str());
    if (w != 0)
      hbkg->Fill(var->getVal(), w->getVal()); 
    else     
      hbkg->Fill(var->getVal());
  }
  
  //TFile* out = new TFile("out.root", "recreate");
  //hsig->Write();
  //hbkg->Write();
  //out->Close();
    
  double q[1];
  double p[1];
  p[0]= target;
  if (v.getDir() == "<") 
    GetQuantiles(hsig, 1, q, p, 1);
  else if (v.getDir() == ">")
    GetQuantiles(hsig, 1, q, p, 0);

  hsig->Divide(hbkg);
  //
  int bin = hsig->FindBin(q[0]);
  float sob = hsig->GetBinContent(bin);
  
  delete hsig;
  delete hbkg;
  return std::pair<float, float>(sob, q[0]);  
}

int main(int argc, char** argv) {
  
  // read the configuration file and book histograms 
  std::map<std::string, std::string> values = parseConfigFile("config.cfg");

  // get list of variables to use + weight 
  std::vector<Variable> variables;
  getVariables(variables, values);
  
  std::string weightVar("");
  std::vector<RooRealVar> rooVars;
  for (unsigned int v=0; v<variables.size(); v++) {
    rooVars.push_back(RooRealVar(variables[v].getName().c_str(), "", double(variables[v].getLow()), double(variables[v].getHigh())));
    if (variables[v].getType() == "w")
      weightVar = variables[v].getName();
  }
  
  // SETUP TREES AND VARIABLES
  RooDataSet* sigDataSet, *bkgDataSet;
  std::string preselCut;
  std::string outputFile = getString(values, "outputFile");
  binning = getUint(values, "binning");

  if (!is_file_exist(outputFile.c_str())) { 
    std::cout << "Creating dataset for optimization..." << std::endl;
    std::string fileName, treeName;
  
    outputFile = getString(values, "outputFile");
    fileName = getString(values, "signalFile");
    treeName = getString(values, "signalTree");
    preselCut = getString(values, "signalPreselCut");
    
    TFile* f = new TFile(fileName.c_str());
    TTree* t = (TTree*)f->Get(treeName.c_str());
    RooArgSet set;
    for (unsigned int i=0; i<rooVars.size(); i++) 
      set.add(rooVars[i]);
    
    sigDataSet = new RooDataSet("sigdata", "", t, set, preselCut.c_str(), weightVar.c_str());
    
    fileName = getString(values,  "bkgFile");
    treeName = getString(values,  "bkgTree");
    preselCut = getString(values, "bkgPreselCut");

    // SETUP TREES AND VARIABLES
    f = new TFile(fileName.c_str());
    t = (TTree*)f->Get(treeName.c_str());

    bkgDataSet = new RooDataSet("bkgdata", "", t, set, preselCut.c_str(), weightVar.c_str());
    TFile* out = new TFile(outputFile.c_str(), "recreate");
    RooWorkspace* w = new RooWorkspace("id_opt");
    w->import(*bkgDataSet);
    w->import(*sigDataSet);
    w->Write();
    out->Close();
    
    std::cout << "Dataset for optimization created. Please run again the same command." << std::endl;
    return 0;

  } else {
    TFile* out = TFile::Open(outputFile.c_str());
    RooWorkspace* w = (RooWorkspace*)out->Get("id_opt");
    bkgDataSet = (RooDataSet*)w->data("bkgdata");
    sigDataSet = (RooDataSet*)w->data("sigdata");
    out->Close();
  }

  // ADD PRELIMINARY STEP TO MOVE ALL VARS TO 100% EFFICIENCY POINT
  for (unsigned int v=0; v<variables.size(); v++) {
    if (variables[v].getType() == "v") {
      std::pair<float, float> result = makeHistos(sigDataSet, bkgDataSet, variables[v], 1., weightVar);
      variables[v].setLimit(result.second);
    }
  } 

  // MAIN LOOP OVER EFF TARGETS
  std::ofstream myfile;
  myfile.open("results.txt", ios::out);
  
  float oldTarget = 1.;
  for (float t=1.; t>.50; t-=.005) {
    std::cout << "ITERATION: " << t << std::endl;
    float realTarget = t/oldTarget;
    
    // LOOP OVER EACH VARIABLE
    float highestSOB = -99999;
    std::string maxVar = "";
    float maxCut = 99999.;
    for (unsigned int v=0; v<variables.size(); v++) {
      if (variables[v].getType() == "v") {
	std::pair<float, float> result = makeHistos(sigDataSet, bkgDataSet, variables[v], realTarget, weightVar); 
	std::cout << variables[v].getName() << " " << result.first << " " << result.second << std::endl;
	if (result.first > highestSOB) {
	  highestSOB = result.first;
	  maxVar = variables[v].getName();
	  maxCut = result.second;
	}
      }
    }
    
    myfile << "+++++++++++++++++++++++++++++\n";
    std::ostringstream new_selection; 
    if (maxVar == "") {
      std::cerr << "WARNING: no improvement found." << std::endl;
      return 0;
    } else {
      for (unsigned int v=0; v<variables.size(); v++) {
	if (variables[v].getName() == maxVar) {
	  variables[v].setLimit(maxCut);
	  myfile << "Chosen: " << variables[v].getName() << " Eff: " << t << " SoB: " << highestSOB << "\n";	  
	}    

	if (variables[v].getType() == "v") {
	  if (v > 0) 
	    new_selection << " && ";
	  new_selection << variables[v].getCutString();
	}
      }

      myfile << new_selection.str() << "\n";
    }
      
    // PRUNING SIGNAL AND BACKGROUND 
    RooDataSet* newSig = (RooDataSet*)sigDataSet->reduce(new_selection.str().c_str());
    RooDataSet* newBkg = (RooDataSet*)bkgDataSet->reduce(new_selection.str().c_str());
    
    sigDataSet = (RooDataSet*)newSig->Clone();
    bkgDataSet = (RooDataSet*)newBkg->Clone();

    delete newSig;
    delete newBkg;

    oldTarget = t;
    // FILL POSSIBLE CONTROL PLOTS
    //g1 = ROOT.TGraph(len(x), x, y)
    //
    //g1.Draw("AP")    
  }

  myfile.close();
}
