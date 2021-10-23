#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TPad.h>
#include "Common.C"
#include "RPMTCalibrator.C"

void MakeRPMTCalibrate(TString graphfile_name,
                       TString treefile_name,
                       const Double_t source_to_rpmt_m,
                       const Double_t thre_overlap_ms = 10.0,
                       const Double_t eff_rpmt_nm_inv = 2.528
                       // eff=1-exp(-2.528 λ[nm]),see RPMTEfficiency.C .png
                       ) {
  Common::InitROOT();
  TString outfile_name = treefile_name;
  outfile_name.ReplaceAll("_list.root", "_calib.root");
  const TString gr_xcalib_name     = "gr_xstage_xyrpmt";
  const TString gr_ycalib_name     = "gr_ystage_xyrpmt";
  const TString gr_inv_wcalib_name = "gr_integ_xyrpmt";
  TFile *treefile = TFile::Open(treefile_name);
  if (!treefile->IsOpen()) {
    std::cerr << "cannot open " << treefile_name << std::endl;
    exit(1);
  }
  TTree *tree = (TTree*)treefile->Get("T");
  RPMTCalibrator *rpmtcalib
    = new RPMTCalibrator(graphfile_name, gr_xcalib_name, gr_ycalib_name, gr_inv_wcalib_name,
                         source_to_rpmt_m, thre_overlap_ms, eff_rpmt_nm_inv);
  // A new file must be opened before make a new friend tree.
  TFile *outfile = TFile::Open(outfile_name, "recreate");
  std::cout << "Created " << outfile_name << std::endl;
  // This tree has xcalib, ycalib, wrel, wabs, wcalib, and tcalib branches.
  TTree *tree_fr = rpmtcalib->GetFriendTree(tree, "Tcalib", "calibration tree");
  tree_fr->Write();
  std::cout << "Wrote " << outfile_name << std::endl;
}
