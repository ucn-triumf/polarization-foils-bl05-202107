#ifndef RPMTCALIBRATOR_C_
#define RPMTCALIBRATOR_C_

#include <iostream>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TROOT.h>

class RPMTCalibrator {
public:
  RPMTCalibrator(const TString rootfile_name_,
                 const TString gr_xcalib_name_,    // xcalib vs. xy
                 const TString gr_ycalib_name_,    // ycalib vs. xy
                 const TString gr_inv_wrel_name_,  // 1/wrel vs. xy
                 const Double_t source_to_rpmt_m_, // for wabs
                 const Double_t thre_overlap_ms_,  // for wabs
                 const Double_t eff_rpmt_nm_inv_   // for wabs
                 );
  ~RPMTCalibrator();
  // void Calibrate(TTree *tree); //
  TTree* GetFriendTree(TTree *tree_orig,
                       const TString tree_name,
                       const TString tree_title);

private:
  void Init();
  const Double_t conv_factor; // lambda_nm = tof*1e-6*conv_factor/source_to_rpmt_m
  // 395.6034 pdg2021
  const Double_t thre_inv_wrel;
  TString rootfile_name;
  TString gr_xcalib_name;
  TString gr_ycalib_name;
  TString gr_inv_wrel_name;
  TFile    *rootfile;
  TGraph2D *gr_xcalib;
  TGraph2D *gr_ycalib;
  TGraph2D *gr_inv_wrel;
  Double_t source_to_rpmt_m;
  Double_t thre_overlap_ms;
  Double_t eff_rpmt_nm_inv;
};

RPMTCalibrator::RPMTCalibrator(const TString rootfile_name_,
                               const TString gr_xcalib_name_,    // xcalib vs. xy
                               const TString gr_ycalib_name_,    // ycalib vs. xy
                               const TString gr_inv_wrel_name_,  // 1/wrel vs. xy
                               const Double_t source_to_rpmt_m_, // for wabs
                               const Double_t thre_overlap_ms_,  // for wabs
                               const Double_t eff_rpmt_nm_inv_   // for wabs
                               ) : thre_inv_wrel(0.6), conv_factor(395.6) {
  rootfile_name    = rootfile_name_;
  gr_xcalib_name   = gr_xcalib_name_;
  gr_ycalib_name   = gr_ycalib_name_;
  gr_inv_wrel_name = gr_inv_wrel_name_;
  source_to_rpmt_m = source_to_rpmt_m_;
  thre_overlap_ms  = thre_overlap_ms_;
  eff_rpmt_nm_inv  = eff_rpmt_nm_inv_;
  Init();
}

RPMTCalibrator::~RPMTCalibrator() {
  rootfile->Close();
}

void RPMTCalibrator::Init() {
  if (gROOT->GetVersion()[0] == '5') {
    std::cerr
      << "[RPMTCalibrator] !!! Interpolation is too slow in ROOT5 !!!"
      << std::endl;
  }
  rootfile = TFile::Open(rootfile_name);
  if (!rootfile->IsOpen()) {
    std::cerr << "[RPMTCalibrator] cannot open " << rootfile_name << std::endl;
    exit(1);
  }
  gr_xcalib     = (TGraph2D*)rootfile->Get(gr_xcalib_name);
  if (gr_xcalib == NULL) {
    std::cerr << "[RPMTCalibrator] cannot get " << gr_xcalib_name << std::endl;
    exit(1);
  }
  gr_ycalib     = (TGraph2D*)rootfile->Get(gr_ycalib_name);
  if (gr_ycalib == NULL) {
    std::cerr << "[RPMTCalibrator] cannot get " << gr_ycalib_name << std::endl;
    exit(1);
  }
  gr_inv_wrel = (TGraph2D*)rootfile->Get(gr_inv_wrel_name);
  if (gr_inv_wrel == NULL) {
    std::cerr << "[RPMTCalibrator] cannot get " << gr_inv_wrel_name << std::endl;
    exit(1);
  }
  std::cout << "[RPMTCalibrator] Initialize end" << std::endl;
}
#if 0
void RPMTCalibrator::Calibrate(TTree *tree) {
  Double_t xcalib, ycalib, inv_wcalib, wcalib, x, y;
  Int_t    f;
  TBranch *xcalib_branch = tree->Branch("xcalib", &xcalib, "xcalib/D");
  TBranch *ycalib_branch = tree->Branch("ycalib", &ycalib, "ycalib/D");
  TBranch *wcalib_branch = tree->Branch("wcalib", &wcalib, "wcalib/D");
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("f", &f);
  std::cout << "[RPMTCalibrator] Calibrating tree" << std::endl;
  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    if ((i+1)*100%tree->GetEntries() == 0 || i==0) {
      // if (true) {
      std::cout << "\r" << "                                " << std::flush;
      std::cout << "\r" << (i+1) << "/" << tree->GetEntries()
                << " "  << (i+1)*100/tree->GetEntries() << "%" << std::flush;
    }
    tree->GetEvent(i);
    xcalib = gr_xcalib->Interpolate(x, y);
    ycalib = gr_ycalib->Interpolate(x, y);
    inv_wcalib = gr_inv_wrel->Interpolate(x, y);
    if (inv_wcalib < thre_inv_wcalib) {
      wcalib = 0.0;
    } else {
      wcalib = 1./inv_wcalib;
    }
    xcalib_branch->Fill();
    ycalib_branch->Fill();
    wcalib_branch->Fill();
  }
  std::cout << std::endl;
}
#endif

TTree* RPMTCalibrator::GetFriendTree(TTree *tree_orig,
                                     const TString tree_name,
                                     const TString tree_title) {
  TTree *tree_friend = new TTree(tree_name, tree_title);
  // Double_t xcalib, ycalib, inv_wcalib, wcalib, x, y;
  Double_t xcalib, ycalib, wcalib, inv_wrel, wrel, wabs, tcalib, x, y, tof;
  // int f; // I don't know why but cause seg fault in tree->Draw()
  tree_orig->SetBranchAddress("x",   &x);
  tree_orig->SetBranchAddress("y",   &y);
  tree_orig->SetBranchAddress("tof", &tof);
  // tree_orig->SetBranchAddress("f", &f);
  tree_friend->Branch("xcalib", &xcalib, "xcalib/D");
  tree_friend->Branch("ycalib", &ycalib, "ycalib/D");
  tree_friend->Branch("wcalib", &wcalib, "wcalib/D");
  tree_friend->Branch("wrel",   &wrel,   "wrel/D");
  tree_friend->Branch("wabs",   &wabs,   "wabs/D");
  tree_friend->Branch("tcalib", &tcalib, "tcalib/D");
  std::cout << "[RPMTCalibrator] Filling friend tree" << std::endl;
  for (Long64_t i = 0; i < tree_orig->GetEntries(); i++) {
    if ((i+1)*100%tree_orig->GetEntries() == 0 || i==0) {
      std::cout << "\r" << "                                " << std::flush;
      std::cout << "\r" << (i+1) << "/" << tree_orig->GetEntries()
                << " "  << (i+1)*100/tree_orig->GetEntries() << "%" << std::flush;
    }
    tree_orig->GetEvent(i);
    xcalib = gr_xcalib->Interpolate(x, y);
    ycalib = gr_ycalib->Interpolate(x, y);
    inv_wrel = gr_inv_wrel->Interpolate(x, y);
    if (inv_wrel < thre_inv_wrel) wrel = 0.0;
    else                          wrel = 1./inv_wrel;
    // Double_t lambda_nm;
    // if (tof*1e-3 < thre_overlap_ms) {
    //   lambda_nm = (tof+40.e3)*1e-6 * conv_factor / source_to_rpmt_m;
    // } else {
    //   lambda_nm = tof*1e-6 * conv_factor / source_to_rpmt_m;
    // }
    if (tof*1e-3 < thre_overlap_ms) tcalib = tof*1e-3 + 40.;
    else                            tcalib = tof*1e-3;
    Double_t lambda_nm = (tcalib*1e-3) * conv_factor / source_to_rpmt_m;
    Double_t eff       = 1. - exp(-eff_rpmt_nm_inv * lambda_nm);
    wabs = 1./eff;
    wcalib = wabs * wrel;
    tree_friend->Fill();
  }
  std::cout << std::endl;
  tree_orig->AddFriend(tree_friend);
  std::cout << "[RPMTCalibrator] add friend" << std::endl;
  return tree_friend;
}
#endif // RPMTCALIBRATOR_C_
