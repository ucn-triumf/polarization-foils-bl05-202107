// File path: 
// polarization-foils-bl05-202107/codes/simple/Draw2D.C
//
// Run this code under polarization-foils-bl05-202107/ as 
// > root -l codes/simple/Draw2D.C

#include <iostream>
#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TEventList.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>

#include "../../tools/ichikawa/RPMT.h"
#include "../../tools/ichikawa/NikiControllerX.C"

TString path_R = "results/"; // path to the results directory 
TString path_D = "data/210713_SiFe/"; // path to the data directory 
TString rootfile  = "20210714193654_list.root"; // name of the target root file

void Draw2D() {
  TString rootfile_num    = path_D + rootfile;
  
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TFile *rootfile  = TFile::Open(rootfile_num.Data());
  TTree *tree = rootfile->Get<TTree>(tree_name);
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas();
  TH2D *h_xy = new TH2D("h_xy","XY from downstr.",640,0,128,640,0,128);
  tree->Draw("y*128:(1-x)*128>>h_xy",cut_str*cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  c_xy->Update();
  
  const Double_t sum  = h_xy->GetSumOfWeights();
  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;
  

}
