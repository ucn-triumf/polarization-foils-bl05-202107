#include "MakeNiki.h"
#include "DrawNikiRPMT.cxx"

extern void InitGui();
Bool_t drawflag = 1;
 
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT api("nikiglass","nikiglass display",initfuncs);

int main(int argc, char *argv[])
{

  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1001111);//to delete number pallete
  gStyle->SetOptFit(1110);
  TGaxis::SetMaxDigits(3);

  // ******************************************************************************
  // timer for information about running time of code
  TStopwatch timer;
  timer.Start();
  
  std::ostringstream os, osopt;
  char DataFileName[128];

  if(argc==1 || argc>3){
    std::cerr<<"Usage: ./MakeNiki [DataFileName] [Option Make Both:0(default), Only RPMT:1, Only single:2 (including UCN case)]" <<std::endl;
    return -1;
  } else {
    os <<argv[1];
    sprintf(DataFileName,os.str().c_str());
    std::cout <<"Making ROOT file from "<<DataFileName<<std::endl;
    if(argc==2){
      nOption = 0;
    }else if(argc==3){
      osopt <<argv[2];
      nOption = atoi(osopt.str().c_str());
      std::cout <<"Option " <<nOption <<" was chosen."<<std::endl;
    } 
  }

  if(nOption==0)   std::cout <<"Create RPMT and Single channel trees."<<std::endl;
  if(nOption==1)   std::cout <<"Create Only a RPMT tree."<<std::endl;
  if(nOption==2)   std::cout <<"Create Only a Single channel tree."<<std::endl;
  if(nOption>2) { 
    std::cerr <<"No such a option."<<std::endl;
    exit(1);
  }

  MakeNiki(DataFileName); //for RPMT 
  if(nOption<2 && drawflag)  DrawNikiRPMT(DataFileName); //for RPMT 

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



//for RPMT, single and UCN data
Int_t *MakeNiki(char* DataFileName)
{
  Int_t i=0,nlines=0;
  Int_t board,ch,e,dummy;  
  Int_t kp=0,tp=0,up=0, mp=0;
  Int_t t0flag=0,t1flag=0,MRflag=0;
  Double_t tof,ucntof;
  Double_t tof_rpmt;
  long long int t=0;            // use long long t t!!
  long long int t_rpmt=0;       // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  long long int tzeroucn=0;     // use long long t t!!
  Int_t a=0, b=0, c=0, d=0, f=0;
  Double_t x=0, y=0;

  TString filestr = gSystem->BaseName(DataFileName);
  Int_t numstr =  filestr.Last('_') ;
  TString dummystr = filestr(numstr+1,3);
  Int_t filenum = dummystr.Atof();

  struct sCntl sRpmt;
  sRpmt=ClearCntl(sRpmt);

  //RPMT root file
  TTree *tup = new TTree("T","");
  tup->Branch("a",     &a,     "a/I");
  tup->Branch("b",     &b,     "b/I");
  tup->Branch("c",     &c,     "c/I");
  tup->Branch("d",     &d,     "d/I");
  tup->Branch("f",     &f,     "f/I");
  tup->Branch("t",     &t_rpmt,"t/L");
  tup->Branch("tof",   &tof_rpmt,"tof/D");
  tup->Branch("kp",    &kp,    "kp/I");
  tup->Branch("tp",    &tp,    "tp/I");
  tup->Branch("up",    &up,    "up/I");
  tup->Branch("mp",    &mp,    "mp/I");
  tup->Branch("MRflag",&MRflag,"MRflag/I");
  tup->Branch("x",     &x,      "x/D");
  tup->Branch("y",     &y,      "y/D");
  //  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:kp:tp:up:mp:MRflag:x:y"); //"f" is coincidence number of event: use f==4

  //other root file
  TTree *six = new  TTree("Six","");
  six->Branch("e",     &e,     "e/I");
  six->Branch("t",     &t,     "t/L");
  six->Branch("tof",   &tof,   "tof/D");
  six->Branch("ucntof",&ucntof,"ucntof/D");
  six->Branch("kp",    &kp,    "kp/I");
  six->Branch("tp",    &tp,    "tp/I");
  six->Branch("up",    &up,    "up/I");
  six->Branch("mp",    &mp,    "mp/I");
  six->Branch("MRflag",&MRflag,"MRflag/I");
  six->Branch("ch",    &ch,    "ch/I");
  six->Branch("board", &board, "board/I");

  FILE *fp;

  fp = fopen(DataFileName, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName<<endl;
    return 0;
  }

  // data read start
  while(1){
  
    Int_t retvalue = fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy) ;
    //    printf("%d,%d,%d,%d,%d,%lld,%d\n",nlines,retvalue,board,ch,e,t,dummy);

    //open new _xxx.dat file
    if(retvalue==EOF){
      filenum++;
      fclose(fp);
      TString newfilestr  = filestr(0,numstr+1) + Form("%03d.dat",filenum);
      fp = fopen(newfilestr.Data(), "r");
      if(fp){
	cout  << newfilestr.Data() <<" was opened continuously"<<endl;
      } else{
	cout  << newfilestr.Data() <<" was not opened. No more files."<<endl;
	break;
      }     
    }
    if (nlines%1000000==0){
      cout <<" nlines =  "<<nlines<<", kp = "<<kp<<endl;
    }
    
    //This counts the number of timing pulse.
    if(ch==tp_ch && e > kp_LLD && e<kp_HLD){
      tp++;
      t1flag=1;
    }
    //This counts the number of kicker pulse.
    else if(ch==kp_ch && e > kp_LLD && e<kp_HLD){
      kp++;
      tzero=t; 
      t0flag=1;
      //      cout << tof*timebin<<" "<<repitition+dt*timebin<<endl;
      if(MRflag == 2) {
	MRflag = 0;
      } else if(tof*timebin > repitition+dt*timebin) {
	MRflag = 2; //first kp after MR
	mp++;
      }
    }
    //This counts the number of ucn phase pulse.
    else if(ch==ucn_ch && e > kp_LLD && e<kp_HLD){
      up++;
      tzeroucn=t; 
    }

    //for other ch events
    tof = (Double_t)(t-tzero);
    ucntof = (Double_t)(t-tzeroucn);
    if(tof*timebin>repitition) MRflag = 1; //for tof tail flag

    if(nOption==0 || nOption==2)six->Fill();

    //ignore events before first kp.
    if(kp==0 && t0flag==0){
    }

    //if time is later than sett+dt, then Fill data in tuple
    else if(t-sRpmt.sett>dt || t0flag==1 || t1flag==1){
      //      cout <<sRpmt.sett<<endl;
      Double_t RPMT_X, RPMT_Y;
      for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
      if((sRpmt.iPulse[0]+sRpmt.iPulse[1])!=0)  RPMT_X = Double_t(sRpmt.iPulse[0])/Double_t(sRpmt.iPulse[0]+sRpmt.iPulse[1]);
      else  RPMT_X = -1.;
      if((sRpmt.iPulse[2]+sRpmt.iPulse[3])!=0)  RPMT_Y = Double_t(sRpmt.iPulse[2])/Double_t(sRpmt.iPulse[2]+sRpmt.iPulse[3]);
      else  RPMT_Y = -1.;

      a=sRpmt.iPulse[0]; b=sRpmt.iPulse[1]; c=sRpmt.iPulse[2]; d=sRpmt.iPulse[3]; f=sRpmt.iSumFlag;
      t_rpmt = sRpmt.sett; tof_rpmt=sRpmt.sett-tzero; x=RPMT_X; y=RPMT_Y; 
      //      cout << t_rpmt <<endl;
      //fill when flag is
      if(nOption==0 || nOption==1)tup->Fill();
	//	tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],
	//		  sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,kp,tp,RPMT_X,RPMT_Y);
      t0flag=0;
      t1flag=0;

      // initialize structure
      sRpmt = ClearCntl(sRpmt);
	
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
  
  //for last is not taken into account
  
#define CREATE_ROOT 1
#if CREATE_ROOT
  TString ROOTstr = filestr(0,19);
  ROOTstr += ".root";
  TFile file(ROOTstr.Data(),"RECREATE");
  if(nOption==0 || nOption==1)tup->Write();
  if(nOption==0 || nOption==2)six->Write();
  cout << ROOTstr.Data() << " has been created" << endl;
  file.Close();
#endif
  
  std::cout << "kp = " << kp << std::endl;


  return 0;

}
