#ifndef CUTVARIABLES_H
#define CUTVARIABLES_H

#include "TString.h"

namespace Vars {

  const int nVariables = 9;
  const int nSpectatorVariables = 3;
  
  struct Variables{
    TString name;
    TString nameTmva; // for addition of abs
    char type;        // "F" for float, "I" for int
    bool symmetric;   // cuts symmetric or one-sided 
    bool isUpper;     // upper bound cut (e.g. iso)
    float lowLimit;
    float upLimit;
    
    // Constructor
    Variables(TString nameIn, TString nameTmvaIn, char typeIn, bool symIn, bool isUp, float llim, float ulim):
      name(nameIn), nameTmva(nameTmvaIn), type(typeIn), symmetric(symIn), isUpper(isUp), lowLimit(llim), upLimit(ulim) 
    {};
  };
  
  Variables *variables [nVariables] = {
    new Variables("sieie", "sieie", 'F', false, true, 0., .04),
    new Variables("hoe",   "hoe+0.01*e",   'F', false, true, -1.e30, 100.),
    new Variables("ecal",  "ecal",  'F', false, true, -1.e30, 1.0),
    new Variables("hcal",  "hcal",  'F', false, true, -1.e30, 1.0),
    new Variables("eop",   "eop",   'F', false, true, 0, 1.0),
    new Variables("dphi",  "dphi",  'F', false, true, 0, 1.0),
    new Variables("deta",  "deta",  'F', false, true, 0, 1.0),
    new Variables("tkiso", "tkiso", 'F', false, true, -1.e30, 1.0),
    new Variables("chi2",  "chi2" , 'F', false, true, -1.e30, 10.0),
  };

  Variables *spectatorVariables [nSpectatorVariables] = {
    new Variables("itype",  "itype",  'I', false, false, 0, 0),
    new Variables("weight", "weight", 'F', false, false, 0, 0),
    new Variables("eta",    "eta",    'F', false, false, 0, 0)
  };
  
};

#endif
