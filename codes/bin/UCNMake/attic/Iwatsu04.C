#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "TMath.h"



//TNtuple *MakeIwatsuTuple(char* DataFileName)
TH1D *MakeIwatsuTuple(char* DataFileName)
{

  Int_t i=0,k=0,nlines=0;
  Int_t x,ch,e,dummy;
  long long int t;
  //  Float_t e,wt,t;

 TNtuple *tup = new  TNtuple("T","T","x:ch:e:t"); 
 // ifstream FD1(DataFileName,ios::in);
 FILE *fp;
 fp = fopen(DataFileName, "r");

 if (!fp) {
   cerr << "File could not be opened \n"<<DataFileName<<endl;
   return 0;
   }

 while(fscanf(fp, "%d,%d,%d,%lld,%d",&x,&ch,&e,&t,&dummy)!=EOF){

     tup->Fill(x,ch,e,t);
     // cout <<x<<" "<<ch<<" "<<e<<" "<<t<<" "<<dummy<<endl;
    nlines ++;
 }

 fclose(fp);

  //  cout <<"nlines="<<nlines<<endl;

 TH1D *h1 = new TH1D(DataFileName,DataFileName,2048,0,8192);
 //  tup->Draw("t","");
 tup->Draw("e>>h1","");

 // TF1 *f1 = new TF1("f1","exp(0)+pol0(1)",0,2048);


 gStyle->SetPalette(1);
#define CREATE_ROOT 1
 char RootFileName[64];
#if CREATE_ROOT
    sprintf(RootFileName,"%s.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  tup->Write();
  cout << RootFileName << " has been created" << endl;
  
  char ClassFileName[64];
  sprintf(ClassFileName,"%s.C",DataFileName);
  T->MakeClass(DataFileName);

#endif

  //  return tup;
  return h1;
}

struct sCntl {
  long long int sett; //setting time
  Int_t iSumFlag; //setting time
  Int_t iPulse[4];
  Int_t iFlag[4];
};
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


//new
TNtuple *MakeIwatsuRPMTkp(char* DataFileName)
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
  
    if (nlines%100000==0){
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
   
   Int_t readFlag =  ReadIwatsuRPMTkp(RootFileName);
   //   Int_t readFlag =  ReadIwatsuRPMTpol(RootFileName);
   
   return six;
}

TNtuple *MakeIwatsuRPMTucn(char* DataFileName)
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
  
    if (nlines%100000==0){
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
   
   Int_t readFlag =  ReadIwatsuRPMTkp(RootFileName);
   //   Int_t readFlag =  ReadIwatsuRPMTpol(RootFileName);
   
   return six;
}


Int_t *ReadIwatsuRPMTkp(char* RootFileName)
{
TCanvas *c1 = new TCanvas("c1","",700,700); 
 TFile *f = TFile::Open(RootFileName);
 //  TNtuple *tup = (TNtuple *)gROOT->FindObject("T");

 //  gPad->SetLogz(1);
  gPad->SetRightMargin(0.15);
  Int_t iii = 2;
  // TH2F *h2 = new TH2F("h2","h2",128,0.,1.,128,0.,1.);
 TH2F *h2 = new TH2F("h2","h2",128*iii,0.,128.,128*iii,0.,128.);
 //  T->Draw("a/(a+b):c/(c+d)>>h2","a>120 && b>60 && c>120 && d>60 && f==4","colz");
 //  T->Draw("a/(a+b):c/(c+d)>>h2","a>60 && b>30 && c>60 && d>30 && f==4","colz");
 //  T->Draw("a/(a+b):c/(c+d)>>h2","a>200 && b>200 && c>200 && d>200 && f==4","colz");
 /* 
    T->Draw("a/(a+b)*128:c/(c+d)*128>>h2",
             "a>200 && b>200 && c>200 && d>200 && f==4 && a/(a+b)>0.3 && a/(a+b)<0.7 && c/(c+d)>0.4 && c/(c+d)<0.7",
             "colz");
 */
 //     T->Draw("a/(a+b)*128:c/(c+d)*128>>h2","a>200 && b>200 && c>200 && d>200 && f==4","colz");
    T->Draw("a/(a+b)*128:c/(c+d)*128>>h2","a>400 && b>400 && c>400 && d>400 && f==4&&a<7400&&b<7400&&c<7400&&d<7400","colz");
 //    T->Draw("a/(a+b)*128:c/(c+d)*128>>h2","a>800 && b>800 && c>800 && d>800 && f==4","colz");
    std::cout << "cut cut cut" << std::endl;
 char GifFileName[64];
    sprintf(GifFileName,"%s.gif",RootFileName);
 c1->SaveAs(GifFileName);

#if 0
TCanvas *c2 = new TCanvas("c2","",700,700); 
 TH1D *ha = new TH1D("ha","ha",2048,0,8192);
 TH1D *hb = new TH1D("hb","hb",2048,0,8192);
 TH1D *hc = new TH1D("hc","hc",2048,0,8192);
 TH1D *hd = new TH1D("hd","hd",2048,0,8192);

   c2->Divide(2,2);
   c2->cd(1);
   gPad->SetLogy(1);
   T->Draw("a>>ha");
   c2->cd(2);
   gPad->SetLogy(1);
   T->Draw("b>>hb");
   c2->cd(3);
   gPad->SetLogy(1);
   T->Draw("c>>hc");
   c2->cd(4);
   gPad->SetLogy(1);
   T->Draw("d>>hd");
  c2->SaveAs("1D.gif");

#endif
  /*
TCanvas *c3 = new TCanvas("c3","",700,700); 
// T->Draw("f","a>80 && b>80 && c>80 && d>80");
// h2->ProjectionY("_py",32,96)->Draw();
 h2->ProjectionY("_py")->Draw();

TCanvas *c4 = new TCanvas("c4","",700,700); 
// h2->ProjectionX("_px",32,96)->Draw();
 h2->ProjectionX("_px")->Draw();
  */
  /*
TCanvas *c5 = new TCanvas("c5","",700,700); 
 T->Draw("a/(a+b)*128:c/(c+d)*128>>h3","a>200 && b>200 && c>200 && d>200 && f==4 && a/(a+b)*128<50","colz");
  */

 TCanvas *c6 = new TCanvas("c6","",700,700); 
 TH1D *htof = new TH1D("htofPMT","htofPMT",40,0,1.0);
 htof -> GetXaxis() -> SetTitle("WaveLength [nm]");
 htof -> GetYaxis() -> SetTitle("Counts/0.025nm");
 T->Draw("tof*395.6*1e-6/16.6>>htofPMT","a>400 && b>400 && c>400 && d>400 && f==4 && a<7400&&b<7400&&c<7400&&d<7400");

 TF1 *fERF = new TF1("fERF","[0]+[1]*TMath::Erf((x-[3])/(2.*[2]*[2]))");


 TCanvas *c7 = new TCanvas("c7","",700,700); 
 TH2F *h3 = new TH2F("h3","h3",256,0.,128.,100,0.,10.);
 T->Draw("tof*3956.*1e-6/16.6:a/(a+b)*128>>h3","a>400 && b>400 && c>400 && d>400 && f==4&&a<7400&&b<7400&&c<7400&&d<7400","colz");
 h3 -> GetYaxis() -> SetTitle("WaveLength [#AA]");
 h3 -> GetXaxis() -> SetTitle("y[mm]");
 h3 -> SetStats(0);

  char PHFileName[64];
#if 0
  sprintf(PHFileName,"%s.tof.root",RootFileName);
  TFile file(PHFileName,"RECREATE");
  htof->Write();
  cout << PHFileName << " has been created" << endl;
#endif

  return 0;
}


//for ch4 kicker
Int_t *MakeIwatsuKPucn(char* DataFileName)
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

    if (nlines%100000==0) cout <<" nlines =  "<<nlines<<endl;
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
    //     cout <<"t= "<<t<<" tzero="<<tzero<<" tof="<<t-tzero<<endl;
    }
  nlines++;
  }  

  //for last event
  if(ch==4) kp++;
  else if(ch==5) tp++;
  else if(ch==6) six->Fill(e,t,t-tzero,kp,tp);
  fclose(fp);
  cout <<"KP= "<<kp<<endl;
  //  cout <<"nlines="<<nlines<<endl;
  // TF1 *f1 = new TF1("f1","exp(0)+pol0(1)",0,2048);
  
#if 1
  char RootFileName[64];
  sprintf(RootFileName,"%s.kp.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  six->Write();
  cout << RootFileName << " has been created" << endl;
  file.Close();
  /*  
  char ClassFileName[64];
  sprintf(ClassFileName,"%s.C",DataFileName);
  Six->MakeClass(DataFileName);
  */
#endif
  
  // Int_t readFlag =  ReadIwatsuKP(RootFileName);
  //    cout <<"total count = "<< readFlag<<endl;
  gROOT->LoadMacro("~/data/ReadUCN.C");
  Int_t readFlag =  ReadUCN(RootFileName);
  return 0;
}


//for ch5 kicker
Int_t *MakeIwatsuKP1(char* DataFileName)
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

    if (nlines%100000==0) cout <<" nlines =  "<<nlines<<endl;
    //for noise reduction
    if(!(nlines<20 && t>1e8)){
      //if data is Trigger(ch4), then change Tzero.
    if(ch==4 && e > 512 && e<8190){
      tp++;
    }
    else if(ch==5 && e > 512 && e<8190){
      tzero=t;
      kp++;
    }
    else if(ch==6){
      six->Fill(e,t,t-tzero,kp,tp);
    }
    //     cout <<"t= "<<t<<" tzero="<<tzero<<" tof="<<t-tzero<<endl;
    }
  nlines++;
  }  

  //for last event
  if(ch==4)tp++;
  else if(ch==5)  kp++;
  else if(ch==6) six->Fill(e,t,t-tzero,kp,tp);
  fclose(fp);
  cout <<"KP= "<<kp<<endl;
  //  cout <<"nlines="<<nlines<<endl;
  // TF1 *f1 = new TF1("f1","exp(0)+pol0(1)",0,2048);
  
#if 1
  char RootFileName[64];
  sprintf(RootFileName,"%s.kp.root",DataFileName);
  TFile file(RootFileName,"RECREATE");
  six->Write();
  cout << RootFileName << " has been created" << endl;
  file.Close();
  /*  
  char ClassFileName[64];
  sprintf(ClassFileName,"%s.C",DataFileName);
  Six->MakeClass(DataFileName);
  */
#endif
  
 Int_t readFlag =  ReadIwatsuKP(RootFileName);
  //    cout <<"total count = "<< readFlag<<endl;
  
  return 0;
}


Int_t *ReadIwatsuKP(char* RootFileName)
{
  //    gStyle->SetPalette(1);
  TFile *f = TFile::Open(RootFileName);
  //  TNtuple *tup = (TNtuple *)gROOT->FindObject("T");
  
  //  TCanvas *c1 = new TCanvas("c1","",700,700); 
  TH2F *h2 = new TH2F("h2","h2",40,0.,40.e3,64,0.,8192.);
  //  TH2F *h2 = new TH2F("h2","h2",100,0.,40.e3,256,0.,512);
  // Six->Draw("e:tof>>h2","","colz");
  //  gPad->SetLogz(1);
  
  TCanvas *c2 = new TCanvas("c2","",700,700); 
  TH1D *ha = new TH1D("ha","ha",1024,0,8192);
  Six->Draw("e>>ha","");

  TCanvas *c4 = new TCanvas("c4","",700,700); 
  Six->Draw("kp>>kicker");

  TCanvas *c5 = new TCanvas("c5","",700,700); 
  Six->Draw("t>>time");

  TCanvas *c6 = new TCanvas("c6","",700,700); 
  TH1F *htof = new TH1F("htof","htof",40,0,1.0);
  // Six->Draw("tof>>htof","e>530");
  Six->Draw("tof*395.6*1e-6/16.4>>htof","e>150");

  char GifFileName[64];
  sprintf(GifFileName,"%s.gif",RootFileName);
  c6->SaveAs(GifFileName);
  
  // Int_t TotalCount=0;
  // TotalCount =  htof->Integral();

#if 0
  char PHFileName[64];
  sprintf(PHFileName,"%s.ph.root",RootFileName);
  TFile file(PHFileName,"RECREATE");
  ha->Write();
  h2->Write();
  htof->Write();
  cout << PHFileName << " has been created" << endl;
  /*  
      char ClassFileName[64];
      sprintf(ClassFileName,"%s.C",DataFileName);
      Six->MakeClass(DataFileName);
  */
#endif

  return 0;

}

Int_t *Go100RPMT(void){
  int i=0;
  char DataFileName[64];
  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f"); 
  for(i=0;i<10;i++){
    sprintf(DataFileName,"%d_list_000.dat",i);
    cout <<DataFileName<<endl;
    tup=MakeIwatsuRPMTtof(DataFileName);
    delete tup;
  }
  return 0;
} 

