#include "MakeNiki.h"
#include "DrawNikiRPMT.cxx"

extern void InitGui();
 
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
    std::cerr<<"Usage: ./MakeNiki [DataFileName] [Option Make Both:0(default), Only RPMT:1, Only single:2, Only ucn:3]" <<std::endl;
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

  if(nOption==0)   std::cout <<"Create RPMT and Single channel ntuples."<<std::endl;
  if(nOption==1)   std::cout <<"Create Only a RPMT ntuples."<<std::endl;
  if(nOption==2)   std::cout <<"Create Only a Single channel ntuple."<<std::endl;
  //  if(nOption==3)   std::cout <<"for UCN, Create Only a Single channel ntuple."<<std::endl;
  if(nOption>2) { 
    std::cerr <<"No such a option."<<std::endl;
    exit(1);
  }

  MakeNiki(DataFileName); //for RPMT 
  if(nOption<2)  DrawNikiRPMT(DataFileName); //for RPMT 

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
  const char* DataFileName2;

  Int_t i=0,nlines=0;
  Int_t board,ch,e,dummy;  
  Int_t kp=0,tp=0,up=0;
  Int_t sp=0,lp=0;
  Int_t t0flag=0,t1flag=0;
  Double_t tof,ucntof;
  Double_t lasertof, shuttertof; 

  //  Int_t dun_ch=5,bm_ch=4;  //local parameter
  Int_t sp_ch=3,lp_ch=6;  //local parameter

  long long int t=0;            // use long long t t!!
  long long int tzero=0;        // use long long t t!!
  long long int tzeroucn=0;     // use long long t t!!
  long long int tzerolaser=0; // use long long t t!!
  long long int tzeroshutter=0; // use long long t t!!
  Double_t dt=1.1; //coincidence time for RPMT channels

  TString filestr = DataFileName;
  Int_t numstr =  filestr.Last('_') ;
  TString dummystr = filestr(numstr+1,3);
  Int_t filenum = dummystr.Atof();

  struct sCntl sRpmt;
  sRpmt=ClearCntl(sRpmt);

  //RPMT root file
  TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:kp:tp:x:y"); //"f" is coincidence number of event: use f==4
  //other root file
  TTree *six = new  TTree("Six","");
  six->Branch("e",    &e,     "e/I");
  six->Branch("t",    &t,     "t/L");
  six->Branch("tof",  &tof,   "tof/D");
  six->Branch("ucntof",&ucntof,  "ucntof/D");
  six->Branch("lasertof",  &lasertof,    "lasertof/D");
  six->Branch("shuttertof",  &shuttertof,    "shuttertof/D");
  six->Branch("kp",   &kp,    "kp/I");
  six->Branch("tp",   &tp,    "tp/I");
  six->Branch("up",   &up,    "up/I");
  six->Branch("ch",   &ch,    "ch/I");
  six->Branch("sp",   &sp,   "sp/I");
  six->Branch("lp",   &lp,   "lp/I");
  six->Branch("board",&board, "board/I");

  Int_t findex = 0;
  for(findex = 0; findex < 7; findex++){

  if(findex == 0) DataFileName2 = "20170413191844_list_000.dat";
  if(findex == 1) DataFileName2 = "20170414122624_list_000.dat";
  if(findex == 2) DataFileName2 = "20170414135219_list_000.dat";
  if(findex == 3) DataFileName2 = "20170414145814_list_000.dat";
  if(findex == 4) DataFileName2 = "20170414150945_list_000.dat";
  if(findex == 5) DataFileName2 = "20170414160939_list_000.dat";
  if(findex == 6) DataFileName2 = "20170414164059_list_000.dat";

  FILE *fp;

  fp = fopen(DataFileName2, "r");
  if (!fp) {
    cerr << "File could not be opened \n"<<DataFileName2<<endl;
    return 0;
  }

  // data read start
  while(1){
  
    Int_t retvalue = fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy) ;
    //    printf("%d,%d,%d,%d,%d,%lld,%d\n",nlines,retvalue,board,ch,e,t,dummy);

    //open new _xxx.dat file
    if(retvalue==EOF){
      fclose(fp);
 	break;
 //     }     
    }
    if (nlines%1000000==0){
      cout <<" nlines =  "<<nlines<<", kp = "<<kp<< " " << findex << endl;
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
    //This counts the number of ucn phase pulse.
    else if(ch==ucn_ch && e > kp_LLD && e<kp_HLD){
      up++;
      tzeroucn=t; 
    }
    else if(ch==lp_ch && e > kp_LLD && e<kp_HLD){
      lp++;
      tzerolaser=t;
    }
    else if(ch==sp_ch && e > kp_LLD && e<kp_HLD){
      sp++;
      tzeroshutter=t;
    }

    //for other ch events
    tof = (Double_t)(t-tzero);
    ucntof = (Double_t)(t-tzeroucn);
    lasertof = (Double_t)(t-tzerolaser);
    shuttertof = (Double_t)(t-tzeroshutter);

    if(nOption==0 || nOption==2)six->Fill();

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
      //fill when flag is
      if(nOption==0 || nOption==1){
	tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],
		  sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,kp,tp,RPMT_X,RPMT_Y);
      }
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

    if(findex == 0 && tp * 0.04 > 3600. * 4.){
       cout<< "timecut" <<endl;
      break;
    }
  }

  }
  
  //for last is not taken into account
  
#define CREATE_ROOT 1
#if CREATE_ROOT
  TString ROOTstr = filestr(0,17);
  ROOTstr += "_sum.root";
  TFile file(ROOTstr.Data(),"RECREATE");
  if(nOption==0 || nOption==1)tup->Write();
  if(nOption==0 || nOption==2)six->Write();
  cout << ROOTstr.Data() << " has been created" << endl;
  file.Close();
#endif
  
  std::cout << "kp = " << kp << std::endl;


  return 0;

}
