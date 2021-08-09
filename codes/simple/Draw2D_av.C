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
TString rootfile  = "20210717054725_list.root"; // name of the target root file


// Double_t x_cut_low = 60;
// Double_t x_cut_up =  68;
Double_t x_cut_low = 40;
Double_t x_cut_up =  55;
Double_t y_cut_low = 50;
Double_t y_cut_up = 80;
Double_t range=128;

TCut cut_xy =Form("x*%f>%f && x*%f<%f && y*%f>%f && y*%f<%f && f==4", range, x_cut_low,range,x_cut_up,range,y_cut_low, range,y_cut_up);


void Draw2D_av() {
  TString rootfile_num    = path_D + rootfile;
  
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TFile *rootfile  = TFile::Open(rootfile_num.Data());
  TTree *tree = rootfile->Get<TTree>(tree_name);
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas("c_xy","c_xy", 800, 1000);
  c_xy->Divide(1,2);
  TH2D *h_xy = new TH2D("h_xy","XY from upstream",640,0,128,640,0,128);
  TH1D *h_x = new TH1D("h_x","X histogram",200,x_cut_low,x_cut_up);
  
  c_xy->cd(1);
  tree->Draw(Form("y*%f:x*%f>>h_xy", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  // TLine *line_xy1= new TLine(x_cut_low, y_cut_low, x_cut_low, y_cut_up);
  // TLine *line_xy2= new TLine(x_cut_low, y_cut_up, x_cut_up, y_cut_up);
  // TLine *line_xy3= new TLine(x_cut_up, y_cut_low, x_cut_up, y_cut_up);
  // TLine *line_xy4= new TLine(x_cut_up, y_cut_up, x_cut_low, y_cut_up);
  // line_xy1->SetLineColor(kRed);
  // line_xy1->Draw();
  // line_xy2->SetLineColor(kRed);
  // line_xy2->Draw();
  // line_xy3->SetLineColor(kRed);
  // line_xy3->Draw();
  // line_xy4->SetLineColor(kRed);
  // line_xy4->Draw();

  c_xy->cd(2);
  tree->Draw(Form("x*%f>>h_x", range), cut_xy * cut_rpmt* TCut(Form("%f",1./time)) , "EH"); //
  // h_x->Fit("gaus");
  Double_t x_peak = h_x->GetMean();
  Double_t count_max = h_x->GetMaximum();
  TLine *line = new TLine(x_peak, 0, x_peak, count_max);
  line->SetLineColor(kRed);
  line->Draw();

  c_xy->Update();

  
  const Double_t sum  = h_xy->GetSumOfWeights();
  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;
  std::cout << "average x: " << x_peak << " mm" << std::endl;
  

}
