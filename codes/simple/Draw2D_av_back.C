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

/*
namestr[0]="20210714193654_list.root";
  namestr[1]="20210716232122_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  namestr[2]="20210717004515_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr[3]="20210717022140_list.root"; 
*/

TString path_R = "results/"; // path to the results directory 
TString path_D = "data/210713_SiFe/"; // path to the data directory 
//TString rootfile  = "20210714185238_list.root"; // name of the target root file
//TString rootfile  = "_list.root"; // name of the target root file
TString rootfile  = "20210716232122_list.root"; // name of the target root file
//up 20210715085144_list
//down 20210715083141_list 20210715073913_list 20210715073913_list 20210715085714_list
/*
  Double_t xbegin=40.;
  Double_t xcenter=55.;
  Double_t xend=71.;
  Double_t ybegin=65.;
  Double_t yend=82.;
*/
///*
//Double_t x_cut_low = 40; // for transmission wave 
//Double_t x_cut_up =  55; 
Double_t x_cut_low = 40; // for transmission wave 
Double_t x_cut_up =  55; // for transmission wave
Double_t x_cut_low3 = 40; // for transmission wave 
Double_t x_cut_up3 =  55; // for transmission wave
Double_t y_cut_low3 = 85;
Double_t y_cut_up3 = 102;
Double_t y_cut_low = 65;
Double_t y_cut_up = 82;
//*/

//BG range (x 20,35 y 65,82)  (40,55 85,102) (40,55 43,60)

/*
Double_t x_cut_low = 0; // for transmission wave 
Double_t x_cut_up =  120; // for transmission wave
Double_t y_cut_low = 0;
Double_t y_cut_up = 120;
*/
Double_t range=128;

TCut cut_xy =Form("x*%f>%f && x*%f<%f && y*%f>%f && y*%f<%f && f==4", range, x_cut_low,range,x_cut_up,range,y_cut_low, range,y_cut_up);


const Int_t nBinXY = 640;
const Double_t startX = 0;
const Double_t endX = 128;
const Double_t startY = 0;
const Double_t endY = 128;
const Int_t nBinCut = 200;


void Draw2D_av_back() {
  TString rootfile_num    = path_D + rootfile;
  
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TFile *rootfile  = TFile::Open(rootfile_num.Data());
  TTree *tree = rootfile->Get<TTree>(tree_name);
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas("c_xy","c_xy", 800, 1000);
  c_xy->Divide(1,3);
  TH2D *h_xy = new TH2D("h_xy","XY from upstream",nBinXY,startX,endX,nBinXY,startY,endY);
  TH1D *h_x = new TH1D("h_x","X histogram",nBinCut,x_cut_low,x_cut_up);
  TH1D *h_y = new TH1D("h_y","Y histogram",nBinCut,y_cut_low,y_cut_up);
  TH2D *h_xy2 = new TH2D("h_xy2","XY from upstream2",nBinXY,x_cut_low,x_cut_up,nBinXY,y_cut_low,y_cut_up);
  TH2D *h_xy3 = new TH2D("h_xy3","XY from upstream3",nBinXY,x_cut_low3,x_cut_up3,nBinXY,y_cut_low3,y_cut_up3);

  c_xy->cd(1);
  tree->Draw(Form("y*%f:x*%f>>h_xy2", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  tree->Draw(Form("y*%f:x*%f>>h_xy3", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)

  tree->Draw(Form("y*%f:x*%f>>h_xy", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  /*
  TLine *line_xy1= new TLine(x_cut_low, y_cut_low, x_cut_low, y_cut_up);
  TLine *line_xy2= new TLine(x_cut_low, y_cut_up, x_cut_up, y_cut_up);
  TLine *line_xy3= new TLine(x_cut_up, y_cut_low, x_cut_up, y_cut_up);
  TLine *line_xy4= new TLine(x_cut_up, y_cut_low, x_cut_low, y_cut_low);
  */
  TLine *line_xy1= new TLine(x_cut_low3, y_cut_low3, x_cut_low3, y_cut_up3);
  TLine *line_xy2= new TLine(x_cut_low3, y_cut_up3, x_cut_up3, y_cut_up3);
  TLine *line_xy3= new TLine(x_cut_up3, y_cut_low3, x_cut_up3, y_cut_up3);
  TLine *line_xy4= new TLine(x_cut_up3, y_cut_low3, x_cut_low3, y_cut_low3);
  line_xy1->SetLineColor(kGreen);
  line_xy1->Draw();
  line_xy2->SetLineColor(kGreen);
  line_xy2->Draw();
  line_xy3->SetLineColor(kGreen);
  line_xy3->Draw();
  line_xy4->SetLineColor(kGreen);
  line_xy4->Draw();

  c_xy->cd(2);
  tree->Draw(Form("x*%f>>h_x", range), cut_xy * cut_rpmt* TCut(Form("%f",1./time)) , "EH"); //
  // h_x->Fit("gaus");
  Double_t x_peak = h_x->GetMean();
  Double_t count_max_x = h_x->GetMaximum();
  TLine *line_x = new TLine(x_peak, 0, x_peak, count_max_x*1.2);
  line_x->SetLineColor(kRed);
  line_x->Draw();

  c_xy->cd(3);
  tree->Draw(Form("y*%f>>h_y", range), cut_xy * cut_rpmt* TCut(Form("%f",1./time)) , "EH"); //
  // h_x->Fit("gaus");
  Double_t y_peak = h_y->GetMean();
  Double_t count_max_y = h_y->GetMaximum();
  TLine *line_y = new TLine(y_peak, 0, y_peak, count_max_y*1.2);
  line_y->SetLineColor(kRed);
  line_y->Draw();
  c_xy->Update();
  
  const Double_t sum  = h_xy->GetSumOfWeights();
  const Double_t sum2  = h_xy2->GetSumOfWeights();
  const Double_t sum3  = h_xy3->GetSumOfWeights();
  const Double_t sum32=sum3/2.;

  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;
  std::cout << "average x: " << x_peak << " mm" << std::endl;
  
  std::cout << "average y: " << y_peak << " mm" << std::endl;

  std::cout << "Green_area sum: "  << sum2     << std::endl;
  std::cout << "Green_area time: " << time     << " s"   << std::endl;
  std::cout << "Green_area rate: " << sum2/time << " cps" << std::endl;

  std::cout << "det_area sum: "  << sum3     << std::endl;
  std::cout << "det_area time: " << time     << " s"   << std::endl;
  std::cout << "det_area rate: " << sum32/time << " cps" << std::endl;

}
