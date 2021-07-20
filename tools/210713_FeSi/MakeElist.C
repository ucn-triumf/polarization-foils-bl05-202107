#include <iostream>
#include "NikiControllerX.C"
#include "TLine.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TEventList.h"
#include "TDirectory.h"

// void MakeElist(const TString rootfile_num ="20210613213139",
// 	       const TString gatenetlog_name =
// 	       "../gatenetlog/gatenetlog_scan20210612_chop0_1.txt",
void MakeElist(const TString rootfile_num ="20210714000204",
	       const TString gatenetlog_name =
	       "../gatenetlog/gatenetlog_scan20210713_x_m2_scan_1.txt",
	       const Int_t kp_lag       = 5,
	       const Int_t kp_kill_head = 3,
	       const Int_t kp_kill_tail = 2,
	       const TString tree_name  = "T",
	       const TString cut_str    = "f==4",
	       Int_t step_first = 0,
	       Int_t step_last  = 0
	       // Int_t step_first = 0,
	       // Int_t step_last  = 20
	       ) {
  const TString rootfile_name = Form("%s_list.root",  rootfile_num.Data());
  const TString outfile_name  = Form("%s_elist.root", rootfile_num.Data());

  TTree *tree_whole = NikiControllerX::GetWholeTree(rootfile_name, tree_name);
  std::cout << "GetWholeTree done" << std::endl;
  std::vector<Int_t> kp_head_v, kp_tail_v;
  NikiControllerX::GetAdjustedKPVectors
    (gatenetlog_name, kp_kill_head, kp_kill_tail, kp_lag, kp_head_v, kp_tail_v);
  if (step_first == 0 && step_last == 0) {
    step_last = kp_head_v.size() - 1;
    std::cout << "step_last is set to be " << step_last << std::endl;
  }
  TFile     *outfile   = TFile::Open(outfile_name, "recreate");
  TObjArray *arr_elist = new TObjArray();
  for (Int_t istep = step_first; istep <= step_last; istep++) {
    std::cout << "\r" << "step: " << istep << std::flush;
    TCut cut_base(cut_str);
    TCut cut_step = Form("kp>=%d && kp<=%d", kp_head_v.at(istep), kp_tail_v.at(istep));
    TEventList *elist = new TEventList(Form("elist_%d", istep));
    tree_whole->Draw(Form(">>elist_%d", istep), cut_base&&cut_step);
    arr_elist->Add(elist);
  }
  arr_elist->Write("arr_elist", 1);
  std::cout << "Saved " << outfile_name << std::endl;
  outfile->Close();
}
