#include <iostream> 
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

using namespace std;

int loadCSV() {
    TTree* T = new TTree("T", "T");
    T->ReadFile("results/30nm_mT_P.csv", "H/D:P/D:Er/D", " ");
    Int_t num=T->GetEntries();
    TGraphErrors *gr = new TGraphErrors(num, T->, pol_at_qcut,0,error_pol_at_qcut);
 return 0;   // ifstream infile("results/30nm_mT_P.csv")
}