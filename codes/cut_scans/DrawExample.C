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

TString scan_path = "data/gatenetlog/";
TString data_path = "data/210713_SiFe/";

void DrawExample() {
  // TString rootfile_num    = "20210716235619";
  TString rootfile_num    = "20210717002421";
  TString rootfile_path = data_path + rootfile_num;
  TString scanfile_name   = "scan20210713_flipper_agilent_scan_rough_1";
  TString gatenetlog_name =
    Form(scan_path + "gatenetlog_%s.txt", scanfile_name.Data());
  Int_t iscan = 5;

  const Int_t kp_lag       = 5;
  const Int_t kp_kill_head = 3;
  const Int_t kp_kill_tail = 2;
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TObjArray *arr_elist
    = NikiControllerX::GetArrayElist(rootfile_path, gatenetlog_name,
				     kp_lag, kp_kill_head, kp_kill_tail,
				     tree_name, cut_str);
  TFile *rootfile  = TFile::Open(Form("%s_list.root",  rootfile_path.Data()));
  TTree *tree = rootfile->Get<TTree>(tree_name);
  TEventList *elist = dynamic_cast<TEventList*>(arr_elist->At(iscan));
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas();
  TH2D *h_xy = new TH2D("h_xy","XY from downstr.",640,0,128,640,0,128);
  tree->Draw("y*128:(1-x)*128>>h_xy",
	     cut_rpmt*TCut(Form("%f",1./time)), "colz");
  c_xy->Update();
  TCanvas *c_tof = new TCanvas();
  TH1D *h_tof = new TH1D("h_tof","TOF",400,0.,40.);
  tree->Draw("tof*1e-3>>h_tof",
	     cut_rpmt*TCut(Form("%f",1./time)), "colz");
  c_tof->Update();

  const Double_t sum  = h_xy->GetSumOfWeights();
  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;

}
