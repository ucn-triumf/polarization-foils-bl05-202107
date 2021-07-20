#include "MakeUCN.h"

Int_t *DrawNikiRPMT(char* DataFileName){

  TString RootFileName = Form("%s.rpmt.root",DataFileName);
  TString FigFileName2D = Form("%s.rpmt.2D.png",DataFileName);

  TFile *file = TFile::Open(RootFileName);
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TNtuple* tup=(TNtuple*)file->Get("T");
  Int_t kp = tup->GetMaximum("kp");
  
  TCanvas *c1 = new TCanvas("c1","",1800,600); 
  c1->Divide(3,1); 
  Int_t nbin = 512;
  Double_t range = 128.;
  TH2F *h2D = new TH2F("h2D","XY view from Upstream;X [mm]; Y[mm]",nbin,0.,range,nbin,0.,range);
  TH1F *hTof = new TH1F("hTof","TOF;TOF [ms];counts [1/0.1ms]",400,0,40);
  TH2F *hTof2D = new TH2F("hTof2D","TOF;TOF [ms];X [1/0.1ms]",400,0,40,nbin,0,range);
  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     rpmt_LLD,rpmt_LLD,rpmt_LLD,rpmt_LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  c1->cd(1);
  tup->Draw(Form("y*%f:x*%f>>h2D",range,range),cut_rpmt_basic,"colz");

  h2D->GetYaxis()->SetTitleOffset(1.2);
  c1->cd(2);
  tup->Draw("tof*1e-3>>hTof",cut_rpmt_basic,"eh");
  hTof->GetYaxis()->SetTitleOffset(1.2);
  hTof->SetMinimum(0.);
  c1->cd(3);
  tup->Draw(Form("x*%f:tof*1e-3>>hTof2D",range),cut_rpmt_basic,"colz");
  c1->SaveAs(FigFileName2D);

#if 1
  FILE *tof2D;
  tof2D = fopen("/home/nikiglass/bin/tof.txt", "w");
  if (!tof2D) {
    cerr << "/home/iwatsu/bin/tof.txt could not be opened"<<endl;
    return 0;
  }

  fprintf(tof2D,"%d\n",kp); // number of T0 counts
  fprintf(tof2D,"100\n"); // bin width is 100 usec

  int n = hTof->GetNbinsX();
  float *tmp = hTof->GetArray();
  for(int i=0; i<n; i++){
    fprintf(tof2D,"%lf\n",tmp[i]); 
  }
  fclose(tof2D);
#endif

  return 0 ;
}

