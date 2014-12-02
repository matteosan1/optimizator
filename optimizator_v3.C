#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCuts.h"
//#endif

#include "Variables.hh"
#include "VarCut.hh"

void writeWorkingPoints(const TMVA::Factory *factory, TString cutsOutFileNameBase, int category, int iter) {

  TString cutsFileName = "weights";
  char a[100];
  sprintf(a, "Category%d", category);
  cutsFileName += "/";
  cutsFileName += cutsOutFileNameBase+"_"+TString(a);
  cutsFileName += ".root";
  
  // Loop over four working points
  std::cout << "The working points being saved:\n" << std::endl;
  
  TFile *cutsFile;
  if (iter == 0)
    cutsFile = new TFile(cutsFileName, "recreate");
  else
    cutsFile = new TFile(cutsFileName, "update");

  if( !cutsFile )
    assert(0);
  VarCut *cuts = new VarCut();
    
  sprintf(a, "Category%d_iter%d", category, iter);
  const TMVA::MethodCuts *method = dynamic_cast<TMVA::MethodCuts*> (factory->GetMethod(a));
  if( method == 0 )
    assert(0);
  
  std::vector <double> cutLo;
  std::vector <double> cutHi;

  double signalEff = 0.95-iter*0.05;
  method->GetCuts(signalEff, cutLo, cutHi);
  
  for(uint ivar=0; ivar<cutHi.size(); ivar++)
    cuts->setCutValue(Vars::variables[ivar]->name, cutLo.at(ivar), cutHi.at(ivar));

  //std::cout << "   working point %s\n", .9 << std::endl;
  cuts->print();
  sprintf(a, "cuts_iter%d", iter);
  cuts->Write(a);
  cutsFile->Close();
}

TString getMethodOptions(TString cutMaxFileName, int category, int iter) {

  TString methodOptions = "!H:!V:EffMethod=EffSel:FitMethod=GA";

  TString cutsFileName = "weights";
  char a[100];
  sprintf(a, "Category%d", category);
  cutsFileName += "/";
  cutsFileName += cutMaxFileName+"_"+TString(a);
  cutsFileName += ".root";
  
  VarCut *cuts;
  if (iter == 0)
    cuts = new VarCut();
  else {
    TFile* cutsFile = new TFile(cutsFileName);
    if( !cutsFile )
      assert(0);

    sprintf (a, "cuts_iter%d", iter-1);
    cuts = (VarCut*)cutsFile->Get(a);
    if( !cuts )
      assert(0);
  }

  // Make sure the user defined cut limits array is consistent with the optimization
  // variables set
  //bool checkPassed = true;
  //if( Vars::nVariables != VarLims::nVarLimits )
  //  checkPassed = false;
  //for(int i=0; i<Vars::nVariables; i++){
  //  if( Vars::variables[i]->name != userDefinedCutLimits[i]->name )
  //    checkPassed = false;
  //}
  //if( !checkPassed ){
  //  printf("ERROR: the list of optimization variables is not consistent with the list\n");
  //  printf("       of the variables with user defined cut limits.\n");
  //  assert(0);
  //}

  for(int i=0; i<Vars::nVariables; i++){
    if (not Vars::variables[i]->isUpper)
      methodOptions += TString::Format(":VarProp[%d]=FMax",i);
    else
      methodOptions += TString::Format(":VarProp[%d]=FMin",i);
  }
 
  //Add all cut ranges:
  for(int i=0; i<Vars::nVariables; i++) {
    //std::cout << Vars::variables[i]->name << " " << cutMax->getCutValue(Vars::variables[i]->name) <<  " " << userDefinedCutLimits[i]->isUpper << " " << userDefinedCutLimits[i]->max << std::endl;
    float max = cuts->getLimits(Vars::variables[i]->name, true);
    float min = cuts->getLimits(Vars::variables[i]->name, false);
    methodOptions += TString::Format(":CutRangeMax[%d]=%.6f", i, max);
    methodOptions += TString::Format(":CutRangeMin[%d]=%.6f", i, min);
  }
  
  printf("\nMethod configuration: method options are\n");
  printf("%s\n", methodOptions.Data());
  printf("\n");

  return methodOptions;
}

void optimizator_v3(const char* outputFileName = "TMVA_cic", int category = 0) {

  TFile* input = TFile::Open("../test.root");
  TTree* inputTree = (TTree*)input->Get("new_opttree");

  int nIterations = 20;
  for (int z=0; z<nIterations; z++) {

    TMVA::Tools::Instance();
    char a[100];
    sprintf(a, "%s.root", outputFileName); 
    TFile* outputFile = TFile::Open(a, "RECREATE");
    outputFile->cd();
    TMVA::Factory *factory = new TMVA::Factory(outputFileName, outputFile, "!V:!Silent");
    
    for (int i=0; i<Vars::nVariables; i++) 
      factory->AddVariable(Vars::variables[i]->nameTmva, Vars::variables[i]->type);
    
    for (int i=0; i<Vars::nSpectatorVariables; i++) 
      factory->AddSpectator(Vars::spectatorVariables[i]->name, Vars::spectatorVariables[i]->type);
    
    factory->SetInputTrees(inputTree, "itype < 0", "itype > 0");  
    factory->SetWeightExpression("weight");
    
    if (category == 0)
      sprintf(a, "sieie<0.02 && dphi < 1.0 && deta < 1.0 && eop < 1. && ecal<1.0 && hcal < 1.0 && tkiso < 1.0 && abs(eta)<1.479 && chi2<10.");
    else 
      sprintf(a, "sieie<0.04 && dphi < 1.0 && deta < 1.0 && eop < 1. && ecal<1.0 && hcal < 1.0 && tkiso < 1.0 && abs(eta)>1.479 && chi2<10.");

    factory->PrepareTrainingAndTestTree(a, "nTrain_Signal=0:nTrain_Background=0::nTest_Signal=0:nTest_Background=0:SplitMode=Random:!V");

    sprintf(a, "Category%d_iter%d", category, z);    
    TString methodOptions = getMethodOptions(outputFileName, category, z);

    factory->BookMethod(TMVA::Types::kCuts, a, methodOptions);
    
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    writeWorkingPoints(factory,  outputFileName, category, z);    
    delete factory;
  }
}
