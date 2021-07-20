
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TEventList.h>
#include <TH1.h>
#include <TH2.h>

void DrawScan1D(const TString rootfile_num) {
  auto *rootfile  = TFile::Open(Form("%s_list.root",  rootfile_num.Data()));
  auto *elistfile = TFile::Open(Form("%s_elist.root", rootfile_num.Data()));
  auto *tree      = rootfile->Get<TTree>("T");
  auto *arr_elist = elistfile->Get<TObjArray>("arr_elist");

  const Int_t    num_step = arr_elist->GetEntries();

  // 1D
  TH1D *h1d = new TH1D("h1d", ";step",
  		       num_step, 0, num_step);
  for (Int_t i = 0; i < num_step; i++) {
    TEventList *elist = dynamic_cast<TEventList*>(arr_elist->At(i));
    tree->SetEventList(elist);
    TH1D *htof = new TH1D("htof","",100,0,40);
    tree->Draw("tof*1e-3>>htof","","goff");
    htof->Scale(1./(tree->GetMaximum("kp")/25.));
    h1d->SetBinContent(i+1, htof->GetSumOfWeights());
    htof->Delete();
  }
  TCanvas *c = new TCanvas();
  c->SetGrid();
  h1d->Draw("eh");
}
void DrawScan2D(const Int_t step_nx, const Int_t step_ny, const TString rootfile_num) {
  auto *rootfile  = TFile::Open(Form("%s_list.root",  rootfile_num.Data()));
  auto *elistfile = TFile::Open(Form("%s_elist.root", rootfile_num.Data()));
  auto *tree      = rootfile->Get<TTree>("T");
  auto *arr_elist = elistfile->Get<TObjArray>("arr_elist");
  std::cout << "arr_elist " << arr_elist->GetEntries() << " entries" << std::endl;
  const Int_t    num_step = arr_elist->GetEntries();
  TH2D *h2d = new TH2D("h2d", ";x step;y step",
		       step_nx, 0, step_nx,
		       step_ny, 0, step_ny);
  for (Int_t ix = 0; ix < step_nx; ix++) {
    for (Int_t iy = 0; iy < step_ny; iy++) {
      Int_t i = iy*step_ny + ix;
      TEventList *elist = dynamic_cast<TEventList*>(arr_elist->At(i));
      tree->SetEventList(elist);
      TH2D *hxy = new TH2D("hxy","",128,0.,128.,128,0.,128.);
      tree->Draw("y*128.:(1-x)*128.>>hxy","f==4","goff");
      hxy->Scale(1./(tree->GetMaximum("kp")/25.));
      h2d->SetBinContent(ix+1, iy+1, hxy->GetSumOfWeights());
      hxy->Delete();
    }
  }
  TCanvas *c = new TCanvas();
  c->SetGrid();
  h2d->Draw("colz");
}
void DrawScan(const Int_t step_nx, const Int_t step_ny, const TString rootfile_num ="20210612173802") {
  DrawScan2D(step_nx, step_ny, rootfile_num);
}
void DrawScan(const TString rootfile_num ="20210714000204") {
  DrawScan1D(rootfile_num);
}
