// raman_spectra.cxx
//
// program reads and displays ascii FT-Raman spectra
// of the 2 column format: wavenumber   value
//
// it then fits the line intensities and fits the
// temperature and ortho/para concentration to the
// measured intensities
//
// created: Nov 14, 2001 KK
//
// last modified: 12-13-01 KK
// last modified: 02-10-02 KK


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

#define PI 3.14159265
#define sq(x)  ((x)*(x))

extern void InitGui();
 
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
Int_t *MakeIwatsuUCN(char* DataFileName);
Int_t *MakeIwatsuRPMTucn(char* DataFileName);
Int_t *MakeIwatsuRPMTkp(char* DataFileName);

struct sCntl {
  long long int sett; //setting time
  Int_t iSumFlag; //setting time
  Int_t iPulse[4];
  Int_t iFlag[4];
};
/*
struct sCntl ClearCntl(struct sCntl sss);

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
*/


TROOT api("proton_spectrum","proton spectrum display",initfuncs);


//int main(int argc, char **argv)
int main(int argc, char *argv[])
{

  //  TApplication theApp("App",&argc,argv); //??? why was the filename which exists ignored???
  //  gROOT->Reset();
  //  gStyle->SetPalette(1);
  //  gStyle->SetScreenFactor(2);
  // ******************************************************************************
  // timer for information about running time of code
  TStopwatch timer;
  timer.Start();
  
  if(argc!=2){
    std::cerr<<"Usage: ./MakeIwatsuUCN [DataFileName]" <<std::endl;
    return -1;
  }

  //  std::ostringstream DataFileName;
  // DataFileName << argv[1];
  //  std::cout <<"Making ROOT file from "<<DataFileName.str().c_str()<<std::endl;
  std::ostringstream os;
  os <<argv[1];
  char DataFileName[128];

  sprintf(DataFileName,os.str().c_str());
  std::cout <<"Making ROOT file from "<<DataFileName<<std::endl;
  MakeIwatsuUCN(DataFileName);
  //  MakeIwatsuRPMTucn(DataFileName);
  MakeIwatsuRPMTkp(DataFileName);
  
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

Int_t *MakeIwatsuUCN(char* DataFileName)
{
  Int_t i=0,k=0,nlines=0;
  Int_t kp=0,tp=0;
  Int_t x,ch,e,dummy; 
  long long int t=0,tzero=0;    // use long long t t!!
  
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp"); //"f" is coincidence number of event: use f==4
  
  //fscan is used because of token ","
  FILE *fp;
  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }

  // data read start
  while(fscanf(fp, "%d,%d,%d,%lld,%d",&x,&ch,&e,&t,&dummy)!=EOF){
    
    if (nlines%1000000==0) cout <<" nlines =  "<<nlines<<endl;
    //for noise reduction
    if(!(nlines<20 && t>1e8)){
      //if data is Trigger(ch4), then change Tzero.
      if(ch==4 && e > 512 && e<8190){
	tzero=t;
	kp++;
      }
      else if(ch==5 && e > 512 && e<8190){
	tp++;
      }
      else if(ch==6){
	six->Fill(e,t,t-tzero,kp,tp);
      }
    }
    nlines++;
  }
  
  //for last event
  if(ch==4) kp++;
  else if(ch==5) tp++;
  else if(ch==6) six->Fill(e,t,t-tzero,kp,tp);
  fclose(fp);
  cout <<"KP= "<<kp<<endl;
  
  char RootFileName[64];
  sprintf(RootFileName,"%s.kp.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  six->Write();
  cout << RootFileName << " has been created" << endl;
  file.Close();

  //  gROOT->LoadMacro("~/data/ReadUCN.C");
  //  Int_t readFlag =  ReadUCN(RootFileName);

  return 0;
}

Int_t *MakeIwatsuRPMTucn(char* DataFileName)
{
  Int_t i=0,k=0,nlines=0;
  Int_t x,ch,e,dummy;  
  Int_t kp=0,tp=0;
  Int_t t0flag=0,t1flag=0;
  long long int t=0;            // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  //  Long_t t;               // use Long_t t!!
  Double_t dt=1.1; //coincidence time for RPMT channels

  struct sCntl sRpmt;
  //  sRpmt=ClearCntl(sRpmt);
  sRpmt.sett=0;
  sRpmt.iSumFlag=0;
  for(i=0;i<4;i++) {
    sRpmt.iPulse[i]=0;
    sRpmt.iFlag[i]=0;
  } 

  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:kp:tp"); //"f" is coincidence number of event: use f==4
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp"); //"f" is coincidence number of event: use f==4
  
  //fscan is used because of token ","
  FILE *fp;
  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }
  // data read start
  while(fscanf(fp, "%d,%d,%d,%lld,%d",&x,&ch,&e,&t,&dummy)!=EOF){
    //   if(nlines > 2.e5 && nlines < 3.e5 )    cout <<x<<" "<<ch<<" "<<e<<" "<<t<<" "<<dummy<<endl;
  
    if (nlines%1000000==0){
      cout <<" nlines =  "<<nlines<<endl;
    }
    //This counts the number of kicker pulse.
    if(ch==4 && e > 512) {
      t0flag=1;
    }
    else if(ch==5 && e > 512) {
      t1flag=1;
    }
    //////   if (kp>60000) break;
    //for noise reduction
    if(kp==0 && t0flag==0){
    }
    else{
      //if time is later than sett+dt, then Fill data in tuple
      //      if(t-sRpmt.sett>dt || t0flag==1){
      if(t-sRpmt.sett>dt || t0flag==1 || t1flag==1){
	for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
	tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,kp,tp);
	
	if(t0flag==1) {
	  kp++;
	  tzero=t; 
	  t0flag=0;
	}

	if(t1flag==1) {
	  tp++;
	  t1flag=0;
	}
	
	// initialize structure
	//  sRpmt = ClearCntl(sRpmt); //<- this routine takes time very long! do not use.
	sRpmt.sett=0;
	sRpmt.iSumFlag=0;
	for(i=0;i<4;i++) {
	  sRpmt.iPulse[i]=0;
	  sRpmt.iFlag[i]=0;
	}
	
	//if data is Trigger(ch4), then change Tzero.
	//  if(ch==4 && e > 512) tzero=t; 
	
	//for 3He events
	if(ch==6) six->Fill(e,sRpmt.sett,sRpmt.sett-tzero,kp);
	
	//reseting sett
	sRpmt.sett=t;
      }
      
      //fill new data
      if(ch>=0 && ch<4){
	sRpmt.iFlag[ch]++;
	sRpmt.iPulse[ch]=e;
      }
      
    }//if (kp!=0)
    nlines++;   
    
  }
  
   //for last event
  for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
  //reject the last event (noise)
   tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],sRpmt.sett,sRpmt.iSumFlag,kp);
   
   fclose(fp);
   
   //  cout <<"nlines="<<nlines<<endl;
   
#define CREATE_ROOT 1
   char RootFileName[64];
#if CREATE_ROOT
   sprintf(RootFileName,"%s.rpmt.root",DataFileName);
   TFile file(RootFileName,"RECREATE");
   tup->Write();
   cout << RootFileName << " has been created" << endl;
   file.Close();
#endif

   std::cout << "kp = " << kp << std::endl;
   
   //   Int_t readFlag =  ReadIwatsuRPMTkp(RootFileName);
   //   Int_t readFlag =  ReadIwatsuRPMTpol(RootFileName);
   
   return 0;
}


//new
Int_t *MakeIwatsuRPMTkp(char* DataFileName)
{
  Int_t i=0,k=0,nlines=0;
  Int_t x,ch,e,dummy;  
  Int_t kp=0,tp=0;
  Int_t t0flag=0,t1flag=0;
  long long int t=0;            // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  //  Long_t t;               // use Long_t t!!
  Double_t dt=1.1; //coincidence time for RPMT channels

  struct sCntl sRpmt;
  //  sRpmt=ClearCntl(sRpmt);
  sRpmt.sett=0;
  sRpmt.iSumFlag=0;
  for(i=0;i<4;i++) {
    sRpmt.iPulse[i]=0;
    sRpmt.iFlag[i]=0;
  } 

  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:kp:tp"); //"f" is coincidence number of event: use f==4
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp"); //"f" is coincidence number of event: use f==4
  
  //fscan is used because of token ","
  FILE *fp;
  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }
  // data read start
  while(fscanf(fp, "%d,%d,%d,%lld,%d",&x,&ch,&e,&t,&dummy)!=EOF){
    //   if(nlines > 2.e5 && nlines < 3.e5 )    cout <<x<<" "<<ch<<" "<<e<<" "<<t<<" "<<dummy<<endl;
  
    if (nlines%1000000==0){
      cout <<" nlines =  "<<nlines<<endl;
    }
    //This counts the number of kicker pulse.
    if(ch==5 && e > 512) {
      t0flag=1;
    }
    else if(ch==4 && e > 512) {
      t1flag=1;
    }
    //////   if (kp>60000) break;
    //for noise reduction
    if(kp==0 && t0flag==0){
    }
    else{
      //if time is later than sett+dt, then Fill data in tuple
      //      if(t-sRpmt.sett>dt || t0flag==1){
      if(t-sRpmt.sett>dt || t0flag==1 || t1flag==1){
	for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
	tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,kp,tp);
	
	if(t0flag==1) {
	  kp++;
	  tzero=t; 
	  t0flag=0;
	}

	if(t1flag==1) {
	  tp++;
	  t1flag=0;
	}
	
	// initialize structure
	//  sRpmt = ClearCntl(sRpmt); //<- this routine takes time very long! do not use.
	sRpmt.sett=0;
	sRpmt.iSumFlag=0;
	for(i=0;i<4;i++) {
	  sRpmt.iPulse[i]=0;
	  sRpmt.iFlag[i]=0;
	}
	
	//if data is Trigger(ch4), then change Tzero.
	//  if(ch==4 && e > 512) tzero=t; 
	
	//for 3He events
	if(ch==6) six->Fill(e,sRpmt.sett,sRpmt.sett-tzero,kp);
	
	//reseting sett
	sRpmt.sett=t;
      }
      
      //fill new data
      if(ch>=0 && ch<4){
	sRpmt.iFlag[ch]++;
	sRpmt.iPulse[ch]=e;
      }
      
    }//if (kp!=0)
    nlines++;   
    
  }
  
   //for last event
  for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
  //reject the last event (noise)
   tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],sRpmt.sett,sRpmt.iSumFlag,kp);
   
   fclose(fp);
   
   //  cout <<"nlines="<<nlines<<endl;
   
#define CREATE_ROOT 1
   char RootFileName[64];
#if CREATE_ROOT
   sprintf(RootFileName,"%s.rpmt.root",DataFileName);
   TFile file(RootFileName,"RECREATE");
   tup->Write();
   cout << RootFileName << " has been created" << endl;
   file.Close();
#endif

   std::cout << "kp = " << kp << std::endl;
   
   
   return 0;
}
