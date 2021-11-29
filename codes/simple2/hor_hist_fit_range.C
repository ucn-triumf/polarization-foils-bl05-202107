//積分する時に面倒なため、ビン幅の中央値にプロットするのをやめる

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"

double AC_power_freqency=3.e5;//Hz
double duty=0.2;
double peak_width=1.*duty/AC_power_freqency;

double func(double *xx,double *par){
  double x1=xx[0];
  double peak_start0=par[0];
  double peak_start1=par[1];
  double peak_hight0=par[2];
  double peak_hight1=par[3];

  double peak_g1=par[4];
  double peak_g2=par[5];

  double peak_end0=peak_start0+peak_width;
  double peak_end1=peak_start1+peak_width;
  double FF;

  if(x1<peak_start0){
    FF=0.;
  }

  if(x1>peak_start0){
    if(x1<peak_end0){
      FF=peak_hight0*x1+peak_g1;
    }
    if(x1>peak_end0){
      if(x1<peak_start1){
        FF=0.;
      }
    }
    if(x1>peak_start1){//if(x1>peak_end0){  
      if(x1<peak_end1){
        FF=peak_hight1+x1+peak_g2;
      }
    }
    if(x1>peak_end1){//if(x1>peak_end0){  
      
      FF=0.;
      
    }
    
  
  }

  return FF;
}

TF1 *f0 = new TF1("f0","func",-1.,1.,6);//TF1("",定義した関数の名前,f0の定義範囲min,max.,パラメータ数)


double func1(double *xx,double *par){
  double x1=xx[0];
  double peak_start0=par[0];
  double peak_start1=par[1];
  double peak_hight0=par[2];
  double peak_hight1=par[3];

  double peak_g1=par[4];
  double peak_g2=par[5];

  double peak_end0=peak_start0+peak_width;
  double peak_end1=peak_start1+peak_width;
  double FF;

  if(x1<peak_start0){
    FF=0.;
  }

  if(x1>peak_start0){
    if(x1<peak_end0){
      FF=peak_hight0*x1+peak_g1;
    }
    if(x1>peak_end0){
      if(x1<peak_start1){
        FF=0.;
      }
    }
    if(x1>peak_start1){//if(x1>peak_end0){  
      if(x1<peak_end1){
        FF=peak_hight1+x1+peak_g2;
      }
    }
    if(x1>peak_end1){//if(x1>peak_end0){  
      
      FF=0.;
      
    }
    
  
  }

  return FF;
}

TF1 *f1 = new TF1("f1","func1",-1.,1.,6);//TF1("",定義した関数の名前,f0の定義範囲min,max.,パラメータ数)

double func2(double *xx,double *par){
  double x1=xx[0];
  double peak_start0=par[0];
  double peak_start1=par[1];
  //double peak_hight0=par[2];
  //double peak_hight1=par[3];

  //double peak_g1=par[4];
  //double peak_g2=par[5];

  double peak_end0=peak_start0+peak_width;
  double peak_end1=peak_start1+peak_width;
  double FF;
  FF=peak_start0*x1+peak_start1;

  return FF;
}
TF1 *f2 = new TF1("f2","func2",-1.,1.,2);
TF1 *f3 = new TF1("f3","func2",-1.,1.,2);

void hor_hist_fit_range(){

  TGraph *gr1 = new TGraph("hor1gsc10212100.csv","%lf,%lf");//%sは文字列
  TGraph *gr2 = new TGraph("horgsc10212100.csv","%lf,%lf");//%sは文字列
  //TH1D *h1 = new TH1D ("","",1000,0.,100.e-6);
  TCanvas* c = new TCanvas("c","c",800,600);
  
  TH1F *h1 = new TH1F("hist1","",1000,-0.25075e-6,0.24925e-6);
  TH1F *h2 = new TH1F("hist2","",1000,-0.25075e-6,0.24925e-6);

  int nn=1001;
  double E1[nn],E2[nn];
  for(int i=0; i < 1001; i++) {

    double x1,y1;
    double x2,y2;
    double y1E[nn],y2E[nn];
    
    gr1->GetPoint(i, x1, y1);
    gr2->GetPoint(i, x2, y2);

    //cout<<x1<<"_"<<y1<<endl;
    h1->Fill(x1,y1); // ?
    h2->Fill(x2,y2); // ? 

    double yy1=h1->GetBinContent(i);
    //cout<<"_"<<yy1<<endl;

    y1E[i]=sqrt(y1);
    y2E[i]=sqrt(y2);

    //cout<<y1E[i]<<"_"<<y2E[i]<<endl;
    
    h1->SetError(y1E);
    h2->SetError(y2E);
  } 
  
  /*
  for(int i=0; i < 1000; i++){
    E1[i]=h1->GetBinError(i);
    E2[i]=h2->GetBinError(i);
    cout<<"error1_"<<E1<<endl;
    cout<<"error2_"<<E2<<endl;
  }
  */

  //gr1->Draw("AP");
  h1->Draw("");
  h2->Draw("same");
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.5);
 
  gPad->SetGridx();
  gPad->SetGridy();

  //TCanvas* c1 = new TCanvas("c","c",800,600);
  //gr2->Draw("AP");
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(0.5);
  //gPad->SetGridx();
  //gPad->SetGridy();
  

  gr1->GetXaxis()->SetTitle("sec[s]");//元のプログラムが単位を上の行のように間違えている？
  //mg->GetYaxis()->SetTitle("Neutron surface flux [n/cm^{2}/(10%EnergyWidth)/s/#muA]");
  gr1->GetYaxis()->SetTitle("count/bin");

  f0->SetParameter(0.,1000.);//m //第１引数が変数の番号、第２引数がその値
  //f0->SetParLimits(0.,1.e-6,2.e-6);//m //第１引数が変数の番号、第２引数がその値
  //f0->FixParameter(0.,0.);
  f0->SetParameter(1.,3.e-6);//m //第１引数が変数の番号、第２引数がその値
  //f0->FixParameter(1.,1.5e-6);
  f0->SetParameter(2.,1.);//m //第１引数が変数の番号、第２引数がその値
  //f0->FixParameter(2,1.);
  f0->SetParameter(3.,8.);//m //第１引数が変数の番号、第２引数がその値
  f0->SetNpx(10000);
  //f0->FixParameter(3,6.);
  //f1->SetParameter(0.,30.);//nm 

  //gr1->Fit("f0","","",1.3e-6,5.e-6);

  h1->GetXaxis()->SetRangeUser(-0.2e-6,0.);
  h2->GetXaxis()->SetRangeUser(-0.2e-6,0.);
  
  
  //+10%
  h1->Fit("f2","+","",-0.168e-6,-0.1595e-6);
  h2->Fit("f3","+","",-0.081e-6,-0.075e-6);
  //0%
  //h1->Fit("f2","+","",-0.168e-6,-0.161e-6);
  //h2->Fit("f3","+","",-0.081e-6,-0.07605e-6);
  //-10%
  //h1->Fit("f2","+","",-0.168e-6,-0.162e-6);
  //h2->Fit("f3","+","",-0.081e-6,-0.077e-6);

  //gr1->Fit("f3","+","",3.465e-6,3.515e-6);

  //gr1->SetStats(1); //表示
  gStyle->SetOptFit(0); //fit情報(だけ)を表示

  //gr1->Fit("f2","","",3.465e-6,3.515e-6);
  //*
  Double_t p0a = f2->GetParameter(0);
  Double_t p1a = f2->GetParameter(1);
  Double_t E0a = f2->GetParError(0);
  Double_t E1a = f2->GetParError(1);
  
  Double_t p0b = f3->GetParameter(0);
  Double_t p1b = f3->GetParameter(1);
  Double_t E0b = f3->GetParError(0);
  Double_t E1b = f3->GetParError(1);

  Double_t peakto=-(p1b/p0b)+(p1a/p0a);

  Double_t peaktoEa=sqrt(pow(E1a/p0a,2)+pow(p1a*E0a/pow(p0a,2),2));
  Double_t peaktoEb=sqrt(pow(E1b/p0b,2)+pow(p1b*E0b/pow(p0b,2),2));

  Double_t peaktoE=sqrt(pow(peaktoEa,2)+pow(peaktoEb,2));
    //*/
  //double p_to_p=p1-p0;

  Double_t time_lag=2.e-6-peakto;
  cout<<"peak_to_peak(プロットの範囲で)_"<<peakto<<"_pm_"<<peaktoE<<"[s]"<<endl;
  cout<<"E0a_"<<E0a<<"E0b_"<<E0b<<"[s]"<<endl;
  cout<<"E1a_"<<E1a<<"E1b_"<<E1b<<"[s]"<<endl;
  cout<<"peaktoEa_"<<peaktoEa<<"peaktoEb_"<<peaktoEb<<"[s]"<<endl;
  cout<<"peak_to_peak(プロット範囲外で)_"<<time_lag<<"_pm_"<<peaktoE<<"[s]"<<endl;

  
  //  cout<<"fit probability : "<< gr1->GetFunction("f0")->GetProb() <<endl;

    return;
}
