// Anaysis for ascii file from Iwatsu a3400 moudle,
// which create /home/iwatsu/bin/tof.txt for auto counting of TOF.
// created: Nov 14, 2001 KK
// last modified: 2014-03-24 Kenji Mishima
// last modified: 2014-12-18 Kenji Mishima for Nikiglass

#ifndef MAKENIKI_H
#define MAKENIKI_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <iostream>
//#include <iostream.h>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
using namespace std;   //to use <iostream> not <iostream.h>
  
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TCollection.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"    
#include "TH3.h"    
#include "TF1.h" 
#include "TF2.h" 
#include "TGraph.h"   
#include "TGraphErrors.h"   
#include "TMinuit.h"   
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"                                                         
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TCut.h"
#include "TGaxis.h"

#define PI 3.14159265
#define sq(x)  ((x)*(x))

#define KPROOT 1 //1:Add kp event in root file 0:not

const Int_t board_ch = 0;
const Int_t monitor_ch = 6;
const Int_t kp_ch = 0 ;
const Int_t tp_ch = 1;
const Int_t kp_LLD = 512;
const Int_t kp_HLD = 8190;
const Int_t rpmt_LLD = 400;
const Int_t rpmt_HLD = 7400;
const Int_t monitor_thr = 150;
const Int_t rpmt_ch1= 2;
const Int_t rpmt_ch2= 3;
const Int_t rpmt_ch3= 4;
const Int_t rpmt_ch4= 5;
const Int_t timing_ch = 7;

Int_t nOption = 0;

Int_t *MakeNikiCh(char* DataFileName);
Int_t *MakeNikiRPMT(char* DataFileName);
Int_t *DrawNikiRPMT(char* DataFileName);

struct sCntl {
  long long int sett; //setting time
  Int_t iSumFlag; //setting time
  Int_t iPulse[4];
  Int_t iFlag[4];
};

struct sCntl ClearCntl(struct sCntl sss)
{
  sss.sett=0;
  sss.iSumFlag=0;
  for(Int_t i=0;i<4;i++) {
    sss.iPulse[i]=0;
    sss.iFlag[i]=0;
  } 
  return sss;
}

#endif
