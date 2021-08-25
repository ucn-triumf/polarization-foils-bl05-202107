#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TEventList.h>
#include <TH1.h>
#include <TH2.h>
#include "RPMT.h"
#include "TCanvas.h"

// void DrawOneScan(TString rootfile_num ="20210612181251", Int_t iscan = 44*10+10) {
// void DrawOneScan(TString rootfile_num ="20210613082746", Int_t iscan = 24+44*19) {
void DrawOneScan(TString rootfile_num ="20210714000204", Int_t iscan = 5) {
  auto *rootfile  = TFile::Open(Form("%s_list.root",  rootfile_num.Data()));
  auto *elistfile = TFile::Open(Form("%s_elist.root", rootfile_num.Data()));
  auto *tree      = rootfile->Get<TTree>("T");
  auto *arr_elist = elistfile->Get<TObjArray>("arr_elist");

  TEventList *elist = dynamic_cast<TEventList*>(arr_elist->At(iscan));
  tree->SetEventList(elist);
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
  // tree->Draw("tof*1e-3>>h_tof",
	    //  TCut(Form("%f",1./time)), "colz");
c_tof->Update();
  Double_t tree_size = tree->GetEntries();
  std::cout << tree_size << std::endl;
  const Double_t sum  = h_xy->GetSumOfWeights();
  std::cout << "sum: "  << sum      << std::endl;
  std::cout << "time: " << time     << " s"   << std::endl;
  std::cout << "rate: " << sum/time << " cps" << std::endl;
  TTree *T =tree->CopyTree("");
  T->SaveAs(Form("../../data_scans/%s_list_%02d.root",rootfile_num.Data(), iscan));
}
