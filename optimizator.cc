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
#include <math.h>

#include "parser.h"

static int binning=5000;

double round_to_digits(double value, int digits) {
  if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
    return 0.0;

  double factor = pow(10.0, digits - ceil(log10(fabs(value))));
  return round(value * factor) / factor;   
}

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
  
  bool checkTollerance(float cut) {
    float diff = 0;
    if (getDir() == "<" or getDir() == "><")
      diff = fabs(cut - getHigh());
    
    if (getDir() == ">")
      diff = fabs(cut - getLow());

    if (diff/cut < 1e-5) 
      return false;
    else
      return true;
  }

  void setLimit(float cut) {
    if (getDir() == "<")
      setHigh(cut);
    
    if (getDir() == ">")
      setLow(cut);

    if (getDir() == "><") {
      setHigh(cut);
      setLow(-cut);
    }
  }
  
  std::string getCutString() {
    std::ostringstream sel;
    if (getDir() == "<")
      sel << getName() << "<"  << getHigh(); 

    if (getDir() == ">")
      sel << getName() << ">" << getLow(); 
    
    if (getDir() == "><")
      sel << "abs(" << getName() << ")<" << getHigh(); 

    return sel.str();
  }

  std::string getRoundedCutString() {
    std::ostringstream sel;
    if (getDir() == "<")
      sel << getName() << "<"  << round_to_digits(getHigh(), 4); 

    if (getDir() == ">")
      sel << getName() << ">" << round_to_digits(getLow(), 4); 
    
    if (getDir() == "><")
      sel << "abs(" << getName() << ")<" << round_to_digits(getHigh(), 4); 

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

void ComputeIntegral(TH1F* h, double* fIntegral, Bool_t direction=true) {
  
  Int_t nbins  = h->GetNbinsX();

  //Double_t* fIntegral = new Double_t[nbins + 2];
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
      return;
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
      return;
    }

    for (Int_t bin=nbins; bin >=0; --bin)  {
      fIntegral[bin] /= fIntegral[0];
      //std::cout << "INTEGRAL:" << bin << " " << h->GetBinLowEdge(bin) << " " << fIntegral[bin] << std::endl;
    }
    fIntegral[0] = h->GetEntries();
  }
  
  //return fIntegral;
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
    Double_t fIntegral[h->GetNbinsX()+ 2];
    ComputeIntegral(h, fIntegral);//h->ComputeIntegral();
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
    Double_t fIntegral[h->GetNbinsX()+ 2];
    ComputeIntegral(h, fIntegral, false);//h->GetIntegral();
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
  
  std::cout << "INPUT VARIABLES:" << std::endl;
  for (unsigned int i=0; i<stringValues.size(); i+=5) {    
    std::vector<std::string> strs;
    std::vector<float> limits;
    getVar(i, stringValues, limits, strs);
  
    Variable v(strs[0]);
    v.setVar(strs[1], limits[0], limits[1], strs[2]);
    variables.push_back(v);

    if (v.getType() == "v")
      std::cout << v.getName() << std::endl;
  }
}
	   
std::pair<float, float> makeHistos(RooDataSet* sig, RooDataSet* bkg, Variable v, float target, std::string weightVar) {

  TH1F* hsig = new TH1F("hsig", "", binning, v.getLow(), v.getHigh());
  TH1F* hbkg = new TH1F("hbkg", "", binning, v.getLow(), v.getHigh());
  
  for (unsigned int i=0; i<sig->numEntries(); i++) {
    const RooArgSet* set = sig->get(i);
    RooRealVar* var = (RooRealVar*)set->find(v.getName().c_str());
    RooRealVar* w = (RooRealVar*)set->find(weightVar.c_str());

    float val = var->getVal();
    if (v.getDir() == "><")
      val = fabs(val);

    if (w != 0)
      hsig->Fill(val, w->getVal());
    else
      hsig->Fill(val);
  }

  for (unsigned int i=0; i<bkg->numEntries(); i++) {
    const RooArgSet* set = bkg->get(i);
    RooRealVar* var = (RooRealVar*)set->find(v.getName().c_str());
    RooRealVar* w = (RooRealVar*)set->find(weightVar.c_str());

    float val = var->getVal();
    if (v.getDir() == "><")
      val = fabs(val);

    if (w != 0)
      hbkg->Fill(val, w->getVal()); 
    else     
      hbkg->Fill(val);
  }
    
  //TFile* out = new TFile("out.root", "recreate");
  //hsig->Write();
  //hbkg->Write();
  //out->Close();

  int bin;
  float sob;
  double q[1];
  double p[1];
  p[0]= target;
  if (v.getDir() == "<" or v.getDir() == "><") {
    GetQuantiles(hsig, 1, q, p, 1);
    bin = hsig->FindBin(q[0]);
    sob = hsig->Integral(0, bin)/hbkg->Integral(0, bin);
  }
  else if (v.getDir() == ">") {
    GetQuantiles(hsig, 1, q, p, 0);
    bin = hsig->FindBin(q[0]);
    sob = hsig->Integral(bin, binning)/hbkg->Integral(bin, binning);
  }
  
  delete hsig;
  delete hbkg;
  return std::pair<float, float>(sob, q[0]);  
}

int main(int argc, char** argv) {
  
  if (argc != 2) {
    std::cerr << "You need to specify the configuration file." << std::endl;
    return 1;
  }

  // read the configuration file and book histograms 
  std::map<std::string, std::string> values = parseConfigFile(argv[1]);

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
  std::vector<float> efficiencyRange = getVfloat(values, "efficiencyRange");
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
    delete out;
  }

  //// ADD PRELIMINARY STEP TO MOVE ALL VARS TO 100% EFFICIENCY POINT
  //for (unsigned int v=0; v<variables.size(); v++) {
  //  if (variables[v].getType() == "v") {
  //    std::pair<float, float> result = makeHistos(sigDataSet, bkgDataSet, variables[v], 1., weightVar);
  //    variables[v].setLimit(result.second);
  //  }
  //} 

  // MAIN LOOP OVER EFF TARGETS
  std::ofstream myfile;
  myfile.open("results.txt", ios::out);
  
  float oldTarget = efficiencyRange[0];
  std::cout << efficiencyRange[0] << " " << efficiencyRange[1] << " " << efficiencyRange[2] << std::endl;
  for (float t=efficiencyRange[0]; t>efficiencyRange[1]; t-=efficiencyRange[2]) {
    
    // ADD PRELIMINARY STEP TO MOVE ALL VARS TO 100% EFFICIENCY POINT
    for (unsigned int v=0; v<variables.size(); v++) {
      if (variables[v].getType() == "v") {
	std::pair<float, float> result = makeHistos(sigDataSet, bkgDataSet, variables[v], 1., weightVar);
	variables[v].setLimit(result.second);
      }
    }

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
	  if (variables[v].checkTollerance(result.second)) {	    
	    highestSOB = result.first;
	    maxVar = variables[v].getName();
	    maxCut = result.second;
	  }
	}
      }
    }
    
    myfile << "+++++++++++++++++++++++++++++\n";
    std::ostringstream new_selection, rounded_selection; 
    if (maxVar == "") {
      std::cerr << "WARNING: no improvement found." << std::endl;
      //return 0;
    } else {
      for (unsigned int v=0; v<variables.size(); v++) {
	if (variables[v].getName() == maxVar) {
	  variables[v].setLimit(maxCut);
	  std::cout << "Chosen: " << variables[v].getName() << " Eff: " << t << " SoB: " << highestSOB << std::endl;;	  
	  myfile << "Chosen: " << variables[v].getName() << " Eff: " << t << " SoB: " << highestSOB << "\n";	  
	}    

	if (variables[v].getType() == "v") {
	  if (v > 0) {
	    new_selection << " && ";
	    rounded_selection << " && ";
	  }
	  new_selection << variables[v].getCutString();
	  rounded_selection << variables[v].getRoundedCutString();
	}
      }

      myfile << rounded_selection.str() << "\n";
    
      
      // PRUNING SIGNAL AND BACKGROUND 
      RooDataSet* newSig = (RooDataSet*)sigDataSet->reduce(new_selection.str().c_str());
      RooDataSet* newBkg = (RooDataSet*)bkgDataSet->reduce(new_selection.str().c_str());
      
      delete sigDataSet;
      delete bkgDataSet;
      
      sigDataSet = newSig;
      bkgDataSet = newBkg;

      myfile.flush();
    }
    
    oldTarget = t;
    // FILL POSSIBLE CONTROL PLOTS
    //g1 = ROOT.TGraph(len(x), x, y)
    //
    //g1.Draw("AP")    
  }

  delete sigDataSet;
  delete bkgDataSet;

  myfile.close();
}
