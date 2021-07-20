// Anaysis for ascii file from Iwatsu a3100 moudle,
// which create /home/iwatsu/bin/tof.txt for auto counting of TOF.
// created: Nov 14, 2001 KK
// last modified: 2014-03-24 Kenji Mishima
// last modified: 2014-12-18 Kenji Mishima


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <iostream>
//#include <iostream.h>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
using namespace std;   //to use <iostream> not <iostream.h>
  
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TClassTable.h"
#include "TCollection.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"    
#include "TH3.h"    
#include "TF1.h" 
#include "TF2.h" 
#include "TGraph.h"   
#include "TGraphErrors.h"   
#include "TMinuit.h"   
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"                                                         
#include "TPaveLabel.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TCut.h"

#define PI 3.14159265
#define sq(x)  ((x)*(x))

#define KPROOT 1 //1:Add kp event in root file 0:not

const Int_t board_ch = 0;
const Int_t monitor_ch = 6;
const Int_t kp_ch = 0 ;
const Int_t tp_ch = 1;
const Int_t kp_LLD = 512;
const Int_t kp_HLD = 8190;
const Int_t rpmt_LLD = 400;
const Int_t rpmt_HLD = 7400;
const Int_t monitor_thr = 150;
const Int_t rpmt_ch1= 2;
const Int_t rpmt_ch2= 3;
const Int_t rpmt_ch3= 4;
const Int_t rpmt_ch4= 5;

extern void InitGui();
 
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
Int_t *MakeIwatsuCh(char* DataFileName);
Int_t *MakeIwatsuRPMT(char* DataFileName);
Int_t *DrawIwatsuRPMT(char* DataFileName);

struct sCntl {
  long long int sett; //setting time
  Int_t iSumFlag; //setting time
  Int_t iPulse[4];
  Int_t iFlag[4];
};

struct sCntl ClearCntl(struct sCntl sss)
{
  sss.sett=0;
  sss.iSumFlag=0;
  for(Int_t i=0;i<4;i++) {
    sss.iPulse[i]=0;
    sss.iFlag[i]=0;
  } 
  return sss;
}

TROOT api("proton_spectrum","proton spectrum display",initfuncs);

int main(int argc, char *argv[])
{

   gStyle->SetPalette(1);

  // ******************************************************************************
  // timer for information about running time of code
  TStopwatch timer;
  timer.Start();
  
  if(argc!=2){
    std::cerr<<"Usage: ./MakeIwatsu [DataFileName]" <<std::endl;
    return -1;
  }

  //  std::ostringstream DataFileName;
  std::ostringstream os;
  os <<argv[1];
  char DataFileName[128];

  sprintf(DataFileName,os.str().c_str());  
  std::cout <<"Making ROOT file from "<<DataFileName<<std::endl;

  MakeIwatsuRPMT(DataFileName); //for RPMT 
  DrawIwatsuRPMT(DataFileName); //for RPMT 

  //  MakeIwatsuCh(DataFileName); //for beam monitor

  // ******************************************************************************
  // timer for information about running time of code
  cout<<"*********************************"<<endl;
  cout<<"ROOT - Time at the end of job:"<<endl;
  cout<<"CpuTime = "<<timer.CpuTime()<<" seconds"<<endl;
  cout<<"RealTime = "<<timer.RealTime()<<" seconds"<<endl;
  cout<<"*********************************"<<endl<<flush;
  timer.Stop();

  return 0;

 }

//for single channel data
Int_t *MakeIwatsuCh(char* DataFileName)
{
  Int_t nlines=0;
  Int_t kp=0,tp=0;
  Int_t board,ch,e,dummy; 
  long long int t=0,tzero=0;    // use long long t t!!
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp:ch:board"); //"f" is coincidence number of event: use f==4

  FILE *fp;
  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }

  // data read start
  while(fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy)!=EOF){
    
    if (nlines%1000000==0) cout <<" nlines =  "<<nlines<<endl;
    //for noise reduction
    if(1){
      //for kicker pulse (beam timing)
      //if data is kp trigger, then reset tzero to 0.
      if(ch==kp_ch && e > kp_LLD && e<kp_HLD){
	tzero=t;
	kp++;
#if KPROOT
	//add kp in root
	six->Fill(e,t,t-tzero,kp,tp,ch,board);
#endif
      }
      //for timing pulse(25Hz)
      else if(ch==tp_ch && e > kp_LLD && e<kp_HLD){
	//	tzero=t;
	tp++;
#if KPROOT
	//add tp in root
	six->Fill(e,t,t-tzero,kp,tp,ch,board);
#endif
      }
      else{
	six->Fill(e,t,t-tzero,kp,tp,ch,board);
      }
    }
    nlines++;
  }
  
  //for last event
  if(ch==kp_ch){
    kp++;
#if KPROOT
    six->Fill(e,t,t-tzero,kp,tp,ch,board);
#endif
  }
  else if(ch==tp_ch){
    tp++;
#if KPROOT
    six->Fill(e,t,t-tzero,kp,tp,ch,board);
#endif
  }
  else six->Fill(e,t,t-tzero,kp,tp,ch,board);
  fclose(fp);
  cout <<"KP= "<<kp<<endl;
  
  //create tof.txt file for Taketani automatic scan system.
  FILE *tofp;
  tofp = fopen("/home/iwatsu/bin/tof.txt", "w");
  if (!tofp) {
    cerr << "/home/iwatsu/bin/tof.txt could not be opened"<<endl;
    return 0;
  }

  fprintf(tofp,"%d\n",kp); // number of T0 counts
  fprintf(tofp,"100\n"); // bin width is 100 usec

  TH1F *htof = new TH1F("htof","TOF",400,0,40);
  six->Draw("tof*1e-3>>htof",Form("e > %d && ch == %d", monitor_thr, monitor_ch),"off");
  int n = htof->GetNbinsX();
  float *tmp = htof->GetArray();
  for(int i=0; i<n; i++){
    fprintf(tofp,"%lf\n",tmp[i]); 
  }
  fclose(tofp);
 
  char RootFileName[64];
  sprintf(RootFileName,"%s.kp.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  six->Write();
  cout << RootFileName << " has been created" << endl;
  file.Close();

  return 0;
}


//for RPMT and other data
Int_t *MakeIwatsuRPMT(char* DataFileName)
{
  Int_t i=0,nlines=0;
  Int_t x,ch,e,dummy;  
  Int_t kp=0,tp=0;
  Int_t t0flag=0,t1flag=0;
  long long int t=0;            // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  long long int tzero_six=0;       // use long long t t!!
  Double_t dt=1.1; //coincidence time for RPMT channels

  struct sCntl sRpmt;
  //  sRpmt=ClearCntl(sRpmt);
  sRpmt.sett=0;
  sRpmt.iSumFlag=0;
  for(i=0;i<4;i++) {
    sRpmt.iPulse[i]=0;
    sRpmt.iFlag[i]=0;
  } 

  //RPMT root file
  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:kp:tp:x:y"); //"f" is coincidence number of event: use f==4
  //other root file
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp"); 
  
  FILE *fp;
  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }
  // data read start
  while(fscanf(fp, "%d,%d,%d,%lld,%d",&x,&ch,&e,&t,&dummy)!=EOF){
  
    if (nlines%100000==0){
      cout <<" nlines =  "<<nlines<<", kp = "<<kp<<endl;
    }

    //This counts the number of kicker pulse.
    if(ch==kp_ch && e > kp_LLD && e<kp_HLD){
      kp++;
      tzero=t; 
      t0flag=1;
    }
    //This counts the number of timing pulse.
    else if(ch==tp_ch && e > kp_LLD && e<kp_HLD){
      tp++;
      t1flag=1;
    }

    //for other ch events
    six->Fill(e,t,tzero_six,kp,tp);

    //ignore events before first kp.
    if(kp==0 && t0flag==0){
    }
    //if time is later than sett+dt, then Fill data in tuple
    else if(t-sRpmt.sett>dt || t0flag==1 || t1flag==1){

      Double_t RPMT_X, RPMT_Y;
      for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
      if((sRpmt.iPulse[0]+sRpmt.iPulse[1])!=0)  RPMT_X = Double_t(sRpmt.iPulse[0])/Double_t(sRpmt.iPulse[0]+sRpmt.iPulse[1]);
      else  RPMT_X = -1.;
      if((sRpmt.iPulse[2]+sRpmt.iPulse[3])!=0)  RPMT_Y = Double_t(sRpmt.iPulse[2])/Double_t(sRpmt.iPulse[2]+sRpmt.iPulse[3]);
      else  RPMT_Y = -1.;
      tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],
		sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,kp,tp,RPMT_X,RPMT_Y);
      
      t0flag=0;
      t1flag=0;

      // initialize structure
      //  sRpmt = ClearCntl(sRpmt); //<- this routine takes time very long! do not use.
      sRpmt.sett=0;
      sRpmt.iSumFlag=0;
      for(i=0;i<4;i++) {
	sRpmt.iPulse[i]=0;
	sRpmt.iFlag[i]=0;
      }
	
      //reseting sett
      sRpmt.sett=t;
    }
      
    //fill new data
    Int_t rpmtch = -1;
    if(ch==rpmt_ch1) rpmtch=0;
    else if(ch==rpmt_ch2) rpmtch=1;
    else if(ch==rpmt_ch3) rpmtch=2;
    else if(ch==rpmt_ch4) rpmtch=3;
    if(rpmtch >= 0){
      sRpmt.iFlag[rpmtch]++;
      sRpmt.iPulse[rpmtch]=e;
    }

    //next line
    nlines++;   
  }

  fclose(fp);

  //for last is not taken into account
  
#define CREATE_ROOT 1
  char RootFileName[64];
#if CREATE_ROOT
  sprintf(RootFileName,"%s.rpmt.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  tup->Write();
  six->Write();
  cout << RootFileName << " has been created" << endl;
  file.Close();
#endif
  
  std::cout << "kp = " << kp << std::endl;
  
  return 0;
}

Int_t *DrawIwatsuRPMT(char* DataFileName){

  TString RootFileName = Form("%s.rpmt.root",DataFileName);
  TString FigFileNameTOF = Form("%s.rpmt.TOF.png",DataFileName);
  TString FigFileName2D = Form("%s.rpmt.2D.png",DataFileName);

  cout << RootFileName.Data() << endl;
  TFile file(RootFileName.Data());
  if ( file.IsOpen() ) printf("ROOT file opened successfully\n");
  TNtuple* tup=(TNtuple*)file.Get("T");

  TCanvas *c1 = new TCanvas("c1","",700,700); 
  TH2F *h2D = new TH2F("h2D","X-Y",128,0.,128.,128,0.,128.);
  TH1F *hTof = new TH1F("hTof2D","X vs TOF",400,0,40);

  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     rpmt_LLD,rpmt_LLD,rpmt_LLD,rpmt_LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  tup->Draw("x*128.:y*128.>>h2D",cut_rpmt_basic,"colz");
  c1->SaveAs(FigFileName2D);

  //  tup->Draw("tof*1e-3>>htof2D","a>400 && b>400 && c>400 && d>400 && f==4&&a<7400&&b<7400&&c<7400&&d<7400 && c/(c+d)*128<80","off");
  //  c1->SaveAs(FigFileNameTOF);


#if 0
  FILE *tof2D;
  tof2D = fopen("/home/nikiglass/bin/tof.txt", "w");
  if (!tof2D) {
    cerr << "/home/iwatsu/bin/tof.txt could not be opened"<<endl;
    return 0;
  }

  fprintf(tof2D,"%d\n",kp); // number of T0 counts
  fprintf(tof2D,"100\n"); // bin width is 100 usec

  int n = htof2D->GetNbinsX();
  float *tmp = htof2D->GetArray();
  for(int i=0; i<n; i++){
    fprintf(tof2D,"%lf\n",tmp[i]); 
  }
  fclose(tof2D);
#endif
   
  return 0 ;
}



