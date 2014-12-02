#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

Float_t deltaPhi(Float_t p1, Float_t p2) {
  float dp = p1-p2;
  if (dp > TMath::Pi()) 
    dp -= (2*TMath::Pi());

  return dp;
}

Float_t deltaR(float e1, float e2, float p1, float p2) {
  float dp = fabs(p1-p2);
  if (dp > TMath::Pi()) 
    dp -= (2*TMath::Pi());

  return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

float getBestPair(int nPairs, float* sumpts) {
  float bestSumpt = 0;
  float bestPair = -1;
  for (int p=0; p<nPairs; p++) {
    if (sumpts[p] > bestSumpt) {
      bestSumpt = sumpts[p];
      bestPair = p;
    }
  }
   
  return bestPair;
}
 
void reduceOptree(const char* outputFile="test.root") {

  TFile* out_file = new TFile(outputFile, "recreate");
  TTree* out_tree = new TTree("new_opttree", "new_opttree");

  Int_t npf, gpn, itype;
  Float_t epf[10], etpf[10], etapf[10], phipf[10], sieiepf[10];
  Float_t ecalpf[10], dphipf[10], detapf[10], hoepf[10], hcalpf[10];
  Float_t tkisopf[10], eoppf[10], gppt[10], gpeta[10], gpphi[10], chi2pf[10];
  Float_t et, eta, phi, sieie, ecal, dphi, deta, hoe, hcal, tkiso, eop, weight, chi2, e;

  out_tree->Branch("e "   , &e    , "e/F");
  out_tree->Branch("et"   , &et   , "et/F");
  out_tree->Branch("eta"  , &eta  , "eta/F");
  out_tree->Branch("phi"  , &phi  , "phi/F");
  out_tree->Branch("sieie", &sieie, "sieie/F");
  out_tree->Branch("ecal" , &ecal , "ecal/F");
  out_tree->Branch("dphi" , &dphi , "dphi/F");
  out_tree->Branch("deta" , &deta , "deta/F");
  out_tree->Branch("hoe"  , &hoe  , "hoe/F");
  out_tree->Branch("hcal" , &hcal , "hcal/F");
  out_tree->Branch("tkiso", &tkiso, "tkiso/F");
  out_tree->Branch("eop"  , &eop  , "eop/F");
  out_tree->Branch("chi2" , &chi2 , "chi2/F");
  out_tree->Branch("itype"  , &itype  , "itype/I");
  out_tree->Branch("weight" , &weight , "weight/F");
  
  const int type[5] = {-1, 1, 2, 3, 4};
  const float weights[5] = {1., .5421, .0737, .0100, .00168};
  const char* names[5] = {"wev.root", "qcd30.root", "qcd50.root", "qcd80.root", "qcd120.root"};
    
  for (int nfiles=0; nfiles<5; nfiles++) {
  
    TFile* file = TFile::Open(names[nfiles]);
    TTree* inputTree = (TTree*)file->Get("tree");

    inputTree->SetBranchStatus("*", 0);
    inputTree->SetBranchStatus("npf"    , 1);
    inputTree->SetBranchStatus("epf"    , 1);
    inputTree->SetBranchStatus("etpf"   , 1);
    inputTree->SetBranchStatus("etapf"  , 1);
    inputTree->SetBranchStatus("phipf"  , 1);
    inputTree->SetBranchStatus("sieiepf", 1);
    inputTree->SetBranchStatus("ecalpf" , 1);
    inputTree->SetBranchStatus("dphipf" , 1);
    inputTree->SetBranchStatus("detapf" , 1);
    inputTree->SetBranchStatus("hoepf"  , 1);
    inputTree->SetBranchStatus("hcalpf" , 1);
    inputTree->SetBranchStatus("tkisopf", 1);
    inputTree->SetBranchStatus("eoppf"  , 1);
    inputTree->SetBranchStatus("chi2pf" , 1);
    inputTree->SetBranchStatus("gpn"    , 1);
    inputTree->SetBranchStatus("gppt"   , 1);
    inputTree->SetBranchStatus("gpeta"  , 1);
    inputTree->SetBranchStatus("gpphi"  , 1);

    inputTree->SetBranchAddress("npf"    , &npf);
    inputTree->SetBranchAddress("epf"    , &epf);
    inputTree->SetBranchAddress("etpf"   , &etpf);
    inputTree->SetBranchAddress("etapf"  , &etapf);
    inputTree->SetBranchAddress("phipf"  , &phipf);
    inputTree->SetBranchAddress("sieiepf", &sieiepf);
    inputTree->SetBranchAddress("ecalpf" , &ecalpf);
    inputTree->SetBranchAddress("dphipf" , &dphipf);
    inputTree->SetBranchAddress("detapf" , &detapf);
    inputTree->SetBranchAddress("hoepf"  , &hoepf);
    inputTree->SetBranchAddress("hcalpf" , &hcalpf);
    inputTree->SetBranchAddress("tkisopf", &tkisopf);
    inputTree->SetBranchAddress("eoppf"  , &eoppf);
    inputTree->SetBranchAddress("chi2pf" , &chi2pf);
    inputTree->SetBranchAddress("gpn"    , &gpn);
    inputTree->SetBranchAddress("gppt"   , &gppt);
    inputTree->SetBranchAddress("gpeta"  , &gpeta);
    inputTree->SetBranchAddress("gpphi"  , &gpphi);
    
    Int_t entries = inputTree->GetEntries();
 
    for (int z=0; z<entries; z++) {
      inputTree->GetEntry(z);
   
      itype = type[nfiles];
      weight = weights[nfiles];

      for (int i=0; i<npf; i++) {
	
	if (etpf[i] < 20. or fabs(etapf[i]>2.1))
	  continue;

	//if (itype < 0) {
	//  // it is signal so match candidates to MC truth
	//  bool isMCSignal = false;
	//  float drmax = 0.2;
	//  for (int j=0; j<gpn; j++) {
	//    float dr = deltaR(etapf[i], gpeta[j], phipf[i], gpphi[j]);
	//      if (dr < drmax) {
	//	isMCSignal = true;
	//	break;
	//      }
	//  }
	//
	//  if (!isMCSignal)
	//    continue;
	//}

	e = epf[i];
	et = etpf[i];
	eta = etapf[i];
	phi = phipf[i];
	sieie = sieiepf[i];
	ecal = ecalpf[i]/etpf[i];
	dphi = dphipf[i];
	deta = detapf[i];
	hoe = hoepf[i];//epf[i];
	hcal = hcalpf[i]/etpf[i];
	tkiso = tkisopf[i]/etpf[i];
	eop = eoppf[i];
	chi2 = chi2pf[i];

	out_tree->Fill();
      }
    }
  }
  
  //print "Write otuput"
  out_file->cd();
  out_tree->Write();
  out_file->Close();
  //file->Close();
}
