#include "Variables.hh"
#include "VarCut.hh"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TEntryList.h"

#include <iostream>

int checkWP(const char* prefix="", int cat=0) {

  char a[100];
  VarCut cut[10];
  sprintf(a, "%s_Category%d.root", prefix, cat);
  TFile* cutFile = new TFile(a);
  for (int i=0; i<20; i++) {
    sprintf(a, "cuts_iter%d", i);
    cut[i] = *((VarCut*)cutFile->Get(a));
    cut[i].getCut()->Print();
  }

  return 0;
  //float weight, id2;
  //Int_t itype, cat, njets20;
  TFile* tuple = new TFile("test.root");
  TTree* tree = (TTree*)tuple->Get("new_opttree");
  tree->SetBranchStatus("*", 1);
  //tree->SetBranchAddress("weight", &weight);
  //tree->SetBranchAddress("itype", &itype);
  //tree->SetBranchAddress("cat", &cat);
  //tree->SetBranchAddress("id2", &id2);
  //tree->SetBranchAddress("njets20", &njets20);

  for (int i=0; i<20; i++) {
    TTree* newTree = tree->CopyTree(*(cut[i].getCut()), "", 10000000, 0);
    newTree->SetBranchStatus("*", 1);
    //newTree->SetBranchAddress("weight", &weight);
    //newTree->SetBranchAddress("itype", &itype);
    //newTree->SetBranchAddress("cat", &cat);
    //newTree->SetBranchAddress("id2", &id2);
    //tree->SetBranchAddress("njets20", &njets20);

    sprintf(a, "out_cat%d_iter%d.root", cat, i);
    TFile* out = new TFile(a, "recreate");
    newTree->Write();
    out->Close();
  }

  return 0;
}

