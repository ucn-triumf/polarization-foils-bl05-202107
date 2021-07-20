#ifndef NIKICONTROLLERX_C_
#define NIKICONTROLLERX_C_

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <TObjArray.h>
#include <vector>
#include <fstream>

namespace NikiControllerX {
Bool_t use_set_max_entry_loop = false;
// Bool_t use_set_max_entry_loop = true;
Int_t  max_entry_loop         = 1E4;

void GetGatenetLog(TString gatenetlog_name,
                   std::vector<Int_t> &kp_head_v,
                   std::vector<Int_t> &kp_tail_v) {
  // read start and end kp of gatenetlog_filename
  kp_head_v.clear();
  kp_tail_v.clear();
  ifstream gatenetlog;
  TString  str_start = "start";
  TString  str_end   = "end";
  TString  read_str;
  Int_t    read_val;
  gatenetlog.open(gatenetlog_name.Data());
  if (gatenetlog.fail()) {
    std::cerr << "cannot open " << gatenetlog_name << std::endl;
    exit(1);
  }
  std::cout << "opened " << gatenetlog_name << std::endl;
  while(gatenetlog >> read_str >> read_val) {
    if      (read_str == str_start) kp_head_v.push_back(read_val);
    else if (read_str == str_end)   kp_tail_v.push_back(read_val);
    else {
      std::cerr << "read_str=" << read_str << " is not valid" << std::endl;
      exit(1);
    }
  }
  if (kp_head_v.size() != kp_tail_v.size()) {
    std::cerr << "kp_head_v(" << kp_head_v.size()
              << ") and kp_tail_v(" << kp_tail_v.size()
              << ") have different size()" << std::endl;
    exit(1);
  }
  std::cout << kp_head_v.size() << " kp windows are loaded." << std::endl;
}

void AdjustKPVectors(Int_t kp_kill_head, Int_t kp_kill_tail, Int_t kp_lag,
                     std::vector<Int_t> &kp_head_v, std::vector<Int_t> &kp_tail_v) {
  // Subtract all kp by the first start kp.
  // kp_lag is the difference between the first kp of
  // Niki and gatenet.
  // And then kill head and tail extra kp just in case.
  Int_t the_first_kp = kp_head_v.at(0);
  for (UInt_t i = 0; i < kp_head_v.size(); i++) {
    kp_head_v.at(i) =  kp_head_v.at(i) - the_first_kp;
    kp_head_v.at(i) += kp_kill_head;
  }
  for (UInt_t i = 0; i < kp_tail_v.size(); i++) {
    kp_tail_v.at(i) =  kp_tail_v.at(i) - the_first_kp;
    kp_tail_v.at(i) -= kp_kill_tail;
  }
  if (kp_lag != 0) {
    std::cout << "kp is shifted due to kp_lag=" << kp_lag << std::endl;
    for (UInt_t i = 0; i < kp_head_v.size(); i++) {
      kp_head_v.at(i) += kp_lag;
      kp_tail_v.at(i) += kp_lag;
    }
  }

  std::cout << the_first_kp << " is the first kp and set to be zero." << std::endl;
  std::cout << kp_kill_head << " kp (head) and "
            << kp_kill_tail << " kp (tail) are killed." << std::endl;
}

void AdjustKPVectors(Int_t kp_kill_head, Int_t kp_kill_tail,
                     std::vector<Int_t> &kp_head_v, std::vector<Int_t> &kp_tail_v) {
  AdjustKPVectors(kp_kill_head, kp_kill_tail, 0, kp_head_v, kp_tail_v);
}

void GetAdjustedKPVectors(TString gatenetlog_name, Int_t kp_kill_head, Int_t kp_kill_tail,
                          Int_t kp_lag,
                          std::vector<Int_t> &kp_head_v, std::vector<Int_t> &kp_tail_v) {
  GetGatenetLog(gatenetlog_name, kp_head_v, kp_tail_v);
  AdjustKPVectors(kp_kill_head, kp_kill_tail, kp_lag, kp_head_v, kp_tail_v);
}
void GetAdjustedKPVectors(TString gatenetlog_name, Int_t kp_kill_head, Int_t kp_kill_tail,
                          std::vector<Int_t> &kp_head_v, std::vector<Int_t> &kp_tail_v) {
  GetAdjustedKPVectors(gatenetlog_name, kp_kill_head, kp_kill_tail, 0,
                       kp_head_v, kp_tail_v);
}
TTree *GetTree(TString rootfile_name, TString tree_name = "T") {
  TFile *file = TFile::Open(rootfile_name);
  if ( file->IsOpen() ) {
    std::cout << rootfile_name << " opened successfully" << std::endl;
  } else {
    std::cerr << rootfile_name << " is not opened" << std::endl;
    exit(1);
  }
  TTree* tup=(TTree*)file->Get(tree_name);
  if(use_set_max_entry_loop==1) {
    std::cout << "use max_entry_loop=" << max_entry_loop << std::endl;
    tup->SetMaxEntryLoop(max_entry_loop);
  }
  std::cout << "tree in " << rootfile_name << " is returned." << std::endl;
  return tup;
}
TTree *GetWholeTree(TString rootfile_name, TString tree_name = "T") {
  return GetTree(rootfile_name, tree_name);
}
TObjArray* GetArrayTCut(TString gatenetlog_name,
                        Int_t kp_kill_head, Int_t kp_kill_tail,
                        Int_t kp_lag = 0) {
  std::vector<Int_t> kp_head_v, kp_tail_v;
  GetAdjustedKPVectors(gatenetlog_name, kp_kill_head, kp_kill_tail,
                       kp_lag, kp_head_v, kp_tail_v);
  TObjArray* arr_cut = new TObjArray();
  for (size_t i = 0; i < kp_head_v.size(); i++) {
    TCut *cut = new TCut(Form("kp>=%d && kp<=%d",
                              kp_head_v.at(i), kp_tail_v.at(i)));
    arr_cut->Add(cut);
  }
  return arr_cut;
}
}
#endif // NIKICONTROLLERX_C_
