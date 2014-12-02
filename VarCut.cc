#include "VarCut.hh"

const int UNDEFCUT = -999;

VarCut::VarCut()
{
  // All cuts are initialized to an unlikely value
  for(int i=0; i<Vars::nVariables; i++) {
    cutMax_[i] = Vars::variables[i]->upLimit;
    cutMin_[i] = Vars::variables[i]->lowLimit;
  }
};

// Construct TCut object for all cuts joined with &&
TCut *VarCut::getCut(){
  
  TCut *cut = 0;
  
  for(int i=0; i<Vars::nVariables; i++){
    if(cutMin_[i] == UNDEFCUT || cutMax_[i] == UNDEFCUT){
      printf("VarCut:: not all cuts are set! Die!\n");
      assert(0);
    }
  }
  
  cut = new TCut("");
  for(int i=0; i<Vars::nVariables; i++) {
    if (Vars::variables[i]->isUpper)
      (*cut) += TString::Format(" %s < %f ", Vars::variables[i]->nameTmva.Data(), cutMax_[i]);
    else
      (*cut) += TString::Format(" %s > %f ", Vars::variables[i]->nameTmva.Data(), cutMin_[i]);
  }
  
  return cut;
}

void VarCut::setCutValue(TString varName, float valMin, float valMax){

  int index = getVariableIndex(varName);

  if( index != -1 ){
    cutMin_[index] = valMin;
    cutMax_[index] = valMax;
  }else{
    printf("VarCut::setCutValue: requested variable is not known!!!\n");
  }

  return;
}

void VarCut::setCutValueTmvaName(TString varNameTmva, float valMin, float valMax){

  int index = getVariableIndexTmvaName(varNameTmva);

  if( index != -1 ){
    cutMin_[index] = valMin;
    cutMax_[index] = valMax;
  } else {
    printf("VarCut::setCutValue: requested variable is not known!!!\n");
  }

  return;
}


float VarCut::getLimits(TString variable, bool isMax) {

  float cutVal = UNDEFCUT;
  int index = getVariableIndex(variable);

  if( index != -1 ){
    if (isMax)
      cutVal = cutMax_[index];
    else
      cutVal = cutMin_[index];
  } else {
    printf("VarCut::getCutValue: requested variable is not known!!!\n");
  }

  return cutVal;
}

float VarCut::getCutValue(TString variable){

  float cutVal = UNDEFCUT;
  int index = getVariableIndex(variable);

  if( index != -1 ){
    if (Vars::variables[index]->isUpper)
      cutVal = cutMax_[index];
    else
      cutVal = cutMin_[index];
  }else{
    printf("VarCut::getCutValue: requested variable is not known!!!\n");
  }

  return cutVal;
}

int VarCut::getVariableIndex(TString variable){

  int index = -1;
  for(int i=0; i<Vars::nVariables; i++){
    if( variable == Vars::variables[i]->name ){
      index = i;
      break;
    } 
  }

  return index;
}

int VarCut::getVariableIndexTmvaName(TString variableTmva){

  int index = -1;
  for(int i=0; i<Vars::nVariables; i++){
    if( variableTmva == Vars::variables[i]->nameTmva ){
      index = i;
      break;
    } 
  }

  return index;
}

bool VarCut::isSymmetric(TString variable){

  bool result = false;
  int index = getVariableIndex(variable);

  if( index != -1 ){
    result = Vars::variables[index]->symmetric;
  }else{
    printf("VarCut::isSymmetric: requestd variable is not known!!!\n");
  }

  return result;
}

// Print all cut values to stdout
void VarCut::print(){

  printf("VarCut::print: Cut values are\n");
  for(int i=0; i<Vars::nVariables; i++) {
    if (Vars::variables[i]->isUpper)
      printf("  %30s < %.3f\n", Vars::variables[i]->nameTmva.Data(), cutMax_[i]);
    else
      printf("  %30s > %.3f\n", Vars::variables[i]->nameTmva.Data(), cutMin_[i]);
  }

}


  
