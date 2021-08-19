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
TString rootfile  = "20210714000204_list.root"; // name of the target root file
void InitColor(){
  
  //set default color
  gROOT->GetColor(2)->SetRGB(220./255.,  50./255.,  47./255.); // Red
  gROOT->GetColor(3)->SetRGB(135./255., 194./255.,  63./255.); // Green
  gROOT->GetColor(4)->SetRGB( 38./255., 139./255., 210./255.); // Blue
  gROOT->GetColor(5)->SetRGB(250./255., 202./255.,  18./255.); // Yellow
  gROOT->GetColor(6)->SetRGB(236./255.,   0./255., 140./255.); // Magenta
  gROOT->GetColor(7)->SetRGB(135./255., 206./255., 250./255.); // Cyan
  gROOT->GetColor(8)->SetRGB(102./255., 205./255., 170./255.); // Lightgreen
  return;
}

TTree* GetTree(TString filestr){

  TString ROOTstr = filestr(0,19);//filestr　ファイルを番号によって読むものを変える
  ROOTstr += ".root";
  TString path = "data/210713_SiFe/";
  TString ROOTstr_path = path+ROOTstr;
  TFile *file = TFile::Open(ROOTstr_path.Data());

  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}


void FeMirror_1() {
  TString rootfile_num    = path_D + rootfile;
  
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TFile *rootfile  = TFile::Open(rootfile_num.Data());
  TTree *tree = rootfile->Get<TTree>(tree_name);
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas();
  TH2D *h_xy = new TH2D("h_xy","XY from downstr.",640,0,128,640,0,128);
  tree->Draw("y*128:(1-x)*128>>h_xy",
	     cut_rpmt*TCut(Form("%f",1./time)), "colz");
  c_xy->Update();
  
  const Double_t sum  = h_xy->GetSumOfWeights();
  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;
  




#if 1
  TFile *outfile = TFile::Open("FeMirrorhist.root","RECREATE");
  for(Int_t i=0; i<num; i++){
    hx[i]->Write();
    hlambda[i]->Write();
    hratio[i]->Write();
    hq[i]->Write();
    hxylambda[i]->Write();
  }
  outfile->Close();
#endif

  return 0 ;
}