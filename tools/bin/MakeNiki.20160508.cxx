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
    std::cerr<<"Usage: ./MakeNiki [DataFileName] [Option Make Both:0(default), Only RPMT:1, Only single:2]" <<std::endl;
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
  if(nOption==3)   std::cout <<"for UCN, Create Only a Single channel ntuple."<<std::endl;

  MakeNikiRPMT(DataFileName); //for RPMT 
  if(nOption<2)  DrawNikiRPMT(DataFileName); //for RPMT 

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



//for RPMT and other data
Int_t *MakeNikiRPMT(char* DataFileName)
{
  Int_t i=0,nlines=0;
  Int_t board,ch,e,dummy;  
  Int_t kp=0,tp=0;
  Int_t t0flag=0,t1flag=0;
  long long int t=0;            // use long long t t!!
  long long int tzero=0;        // use long long t t!!
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
  TNtuple *six = new  TNtuple("Six","Six","e:t:tof:kp:tp:ch:board"); //"f" is coincidence number of event: use f==4
  
  FILE *fp;
  //  FILE fp[100];

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
    six->Fill(e,t,t-tzero,kp,tp,ch,board);

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
  char RootFileName[64];
#if CREATE_ROOT
  TString ROOTstr = filestr(0,17);
  ROOTstr += ".root";
  TFile file(ROOTstr.Data(),"RECREATE");
  if(nOption!=2)tup->Write();
  if(nOption!=1)six->Write();
  cout << ROOTstr.Data() << " has been created" << endl;
  file.Close();
#endif
  
  std::cout << "kp = " << kp << std::endl;


  return 0;

}
