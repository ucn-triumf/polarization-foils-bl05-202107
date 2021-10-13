#include "../bin/MakeNiki.h"
//#include "/home/nop/data/Tools/GetMRcut.C"
#include <TH3.h>
TString path_R = "results/";

// Double_t Distance = 20.000; //tentative
// Double_t Conversion = 395.6;
// Double_t dist_det   = 337.; //sample to detector [mm]
// Double_t xdirect    = 64.71;
// Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
// Bool_t useThinout = 0; //thinning out the event <1e4.

Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 1439.; //sample to detector [mm]
//Double_t xdirect    = 92.4685;//92.4685 92.4736
Double_t xdirect    =92.4736;
//Double_t xdirect    = 93.4728;
Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.

void InitColor(){
  //set default color
  gROOT->GetColor(2)->SetRGB(220./255.,  50./255.,  47./255.); // Red
  gROOT->GetColor(3)->SetRGB(135./255., 194./255.,  63./255.); // Green
  gROOT->GetColor(4)->SetRGB( 38./255., 139./255., 210./255.); // Blue
  gROOT->GetColor(5)->SetRGB(250./255., 202./255.,  18./255.); // Yellow
  gROOT->GetColor(6)->SetRGB(236./255.,   0./255., 140./255.); // Magenta
  gROOT->GetColor(7)->SetRGB(135./255., 206./255., 250./255.); // Cyan
  gROOT->GetColor(8)->SetRGB(102./255., 205./255., 170./255.); // Lightgreen
  return;
}

TTree* GetTree(TString filestr){

  TString ROOTstr = filestr(0,19);//filestr　ファイルを番号によって読むものを変える
  ROOTstr += ".root";
  TString path = "data/210713_SiFe/";
  TString ROOTstr_path = path+ROOTstr;
  TFile *file = TFile::Open(ROOTstr_path.Data());

  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}

double func_R0(double *qqq,double *par){
  double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double qc=par[0];
  double ww=2.5E-03;//par[1];
  double R0=par[1];
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;

  double mm=5.2;
  double mm2=par[2];
  double R1;

  double Rup;
    if(q1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(q1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }

    double Rdown;
    if(q1<qc){
       Rdown=downrate*R0;  
    }
    
    else{
      if(q1>=qc){
        //double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        double down_R=downrate*R0/pow((1.+mm2*(q1-qc)),4);
        Rdown=down_R;//+down_R;
      }
    }
    R1=Rup+Rdown;

    return R1;

  }
TF1 *f0 = new TF1("",func_R0,0.01,0.5,3);

Int_t M1_pol_check(){

  InitColor();
  TH1::SetDefaultSumw2();

  const Int_t num = 7;
  Int_t kp[num];
  TTree* tup[num];
  TH1F* hx[num];
  TH1F* hx2[num];
  TH1F* hlambda[num];
  TH1F* hratio[num];
  TH1F* hq[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TH1F* hlambda2[num];
  TH1F* hratio2[num];
  TH1F* hq2[num];
  TH1F* hq3[num];
  TH1F* hq4[num];
  TH1F* hq02[num];
  TH3F* hxylambda2[num];
  TString namestr[num];

  TH1F* hqplus[num];


  namestr[0]="20210713202138_list.root"; 
  //namestr[0]="20210713212050_list.root"; //M1 reflect (direct) 1hour
  //namestr[0]="20210713203803_list.root";
  namestr[1]="20210713225533_list.root";

  //namestr[2]="20210713215303_list.root";
  namestr[2]="20210713224531_list.root"; 

  namestr[3]="20210713230326_list.root"; //Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1 mT from -8 mT  with AFP 760 mV
  namestr[4]="20210713230941_list.root"; //Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT  with AFP 760 mV
  namestr[5]="20210713231325_list.root"; 
  
  
  namestr[6]="20210715084052_list.root";

  // Hysteresis1&2
  /*
  namestr[1]="20210715072653_list.root"; //1.97A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 8 mT
  namestr[2]="20210715075452_list.root"; //0A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT
  namestr[3]="20210715081447_list.root"; //0A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = -8 mT -> 0 mT
  namestr[4]="20210715084835_list.root"; //0.15A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1 mT from -8 mT
  namestr[5]="20210715085349_list.root"; //0.264A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT
  namestr[6]="20210715082606_list.root"; //0.378A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT  -> 2 mT
  namestr[7]="20210715083711_list.root"; //0.6A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 2 mT  -> 3 mT
  */ 

  TString degstr[num];
  degstr[0]="Direct";
  // degstr[1]="M2 reflect(0.49 deg.)";
  // degstr[2]="M2 reflect(1.00 deg.)";
  // degstr[3]="M2 reflect(0.49 deg.) with AFP 760 mV";
  // degstr[4]="M2 reflect(1.00 deg.) with AFP 760 mV";

  // degstr[1]="Fe 30 nm, #theta = 0.69 deg. x =  0.0 mm";
  // degstr[2]="Fe 30 nm, #theta = 0.69 deg. x = -0.1 mm";
  // degstr[2]="Fe 30 nm, #theta = 0.69 deg. with AFP";
  // degstr[4]="Fe 30 nm, #theta = 0.69 deg. x = -0.3 mm";
  // degstr[5]="Fe 30 nm, #theta = 0.69 deg. x = +0.1 mm";
  // degstr[6]="Fe 30 nm, #theta = 0.69 deg. x = +0.0 mm";
  // degstr[7]="Fe 30 nm, #theta = 0.69 deg. x = -0.1 mm";

  // Hysteresis1
  // degstr[1]="Fe 30 nm, #theta = 0.69 deg., B = 8 mT";
  // degstr[2]="Fe 30 nm, #theta = 0.69 deg., B = 8 mT with AFP";
  // degstr[3]="Fe 30 nm, #theta = 0.69 deg., B = 0 mT" ;
  // degstr[4]="Fe 30 nm, #theta = 0.69 deg., B = 0 mT with AFP" ;
  // degstr[5]="Fe 30 nm, #theta = 0.69 deg., B = 0 mT from -8 mT";
  // degstr[6]="Fe 30 nm, #theta = 0.69 deg., B = 0 mT from -8 mT with AFP";
  // degstr[7]="Fe 30 nm, #theta = 0.69 deg., B = 2 mT from  0 mT";
  // degstr[8]="Fe 30 nm, #theta = 0.69 deg., B = 2 mT from  0 mT with AFP";
  // degstr[9]="Fe 30 nm, #theta = 0.69 deg., B = 3 mT from  2 mT";
  // degstr[10]="Fe 30 nm, #theta = 0.69 deg., B = 3 mT from  2 mT with AFP";

  // Hysteresis2
  // degstr[1]="Fe 30 nm, #theta = 0.69 deg., B = 8 mT";
  // degstr[2]="Fe 30 nm, #theta = 0.69 deg., B = 8 mT with AFP";
  // degstr[3]="Fe 30 nm, #theta = 0.69 deg., B = 1 mT from -8 mT";
  // degstr[4]="Fe 30 nm, #theta = 0.69 deg., B = 1 mT from -8 mT with AFP";
  // degstr[5]="Fe 30 nm, #theta = 0.69 deg., B = 1.5 mT from -8 mT";
  // degstr[6]="Fe 30 nm, #theta = 0.69 deg., B = 1.5 mT from -8 mT with AFP";

  // Hysteresis1&2
  /*
  degstr[1]="B = 8 mT";
  degstr[2]="B = 0 mT from 8 mT" ;
  degstr[3]="B = 0 mT from -8 mT";
  degstr[4]="B = 1 mT from -8mT";
  degstr[5]="B = 1.5 mT from -8mT";
  degstr[6]="B = 2 mT from -8mT";
  degstr[7]="B = 3 mT from -8mT";
  */
/*
  degstr[1]="AFP OFF B = 8.01 mT";
  degstr[2]="B = 0.322 mT from 8 mT" ;
  degstr[3]="B = 0.322 mT from -8 mT";
  degstr[4]="B = 0.908 mT from -8mT";
  degstr[5]="B = 1.35 mT from -8mT";
  degstr[6]="B = 1.80 mT from -8mT";
  degstr[7]="B = 2.66 mT from -8mT";
*/
  degstr[1]=namestr[1];
  degstr[2]=namestr[2];
  degstr[3]=namestr[3];
  degstr[4]=namestr[4];
  degstr[5]=namestr[5];
  degstr[6]=namestr[6];
  Double_t angle[num];
  /*

  */
 /*
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(47.07 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(47.13 - xdirect)/dist_det; //rad
  angle[3] = TMath::Abs(47.1 - xdirect)/dist_det; //rad
  angle[4] = TMath::Abs(47.19 - xdirect)/dist_det; //rad
  angle[5] = TMath::Abs(47.24 - xdirect)/dist_det; //rad
  angle[6] = TMath::Abs(47.21 - xdirect)/dist_det; //rad
  angle[7] = TMath::Abs(47.11 - xdirect)/dist_det; //rad
*/
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  
  //angle[2] = TMath::Abs(68.1399 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(68.0634 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(68.4522 - xdirect)/dist_det; //rad
  angle[3] = TMath::Abs(68.0634 - xdirect)/dist_det; //rad
  angle[4] = TMath::Abs(68.0634 - xdirect)/dist_det; //rad
  angle[5] = TMath::Abs(68.0634 - xdirect)/dist_det; //rad
  angle[6] = TMath::Abs(68.0634 - xdirect)/dist_det; //rad

  Double_t angledeg[num];
  Double_t angledeg2[num];
  
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  angledeg[3]=angle[3]*180./TMath::Pi()/2.;
  angledeg[4]=angle[4]*180./TMath::Pi()/2.;
  angledeg[5]=angle[5]*180./TMath::Pi()/2.;
  angledeg[6]=angle[6]*180./TMath::Pi()/2.;
  angledeg[7]=angle[7]*180./TMath::Pi()/2.;

  

  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  TLegend* leg = new TLegend(0.73, 0.30, 1.00, 0.70,"");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Int_t nbin_q  = 300;
  Double_t q_max  = 1.0;//0.6
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;
  //  Double_t HLD  = 7400.;

  //  Double_t xbegin=54.;
 
  Double_t xbegin=63.;
  Double_t xcenter=73.;
  Double_t xbegin2=86.5;
  Double_t xcenter2=96.5;
  Double_t xcenter1=86.5;
  Double_t xend=96.5;
  Double_t ybegin=55.;
  Double_t yend=85.;
  
  //  Double_t ybegin=70.;
  //  Double_t yend=77.;

  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     LLD,LLD,LLD,LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  TCut cut_x = Form("x*%f>20 && x*%f<100",range,range);
  TCut cut_y = Form("y*%f>%f && y*%f<%f",range,ybegin,range,yend);
  TCut cut_dir = Form("x*%f>%f && x*%f<%f",range,xcenter1,range,xend);
  TCut cut_ref = Form("x*%f>%f && x*%f<%f",range,xbegin,range,xcenter);
  TCut cut_ref2 = Form("x*%f>%f && x*%f<%f",range,xbegin2,range,xcenter2);
  //  TCut cut_tof = Form("tof>1.0e3 && tof<39.9e3");
  TCut cut_tof = "";
  TCut MRcut = "MRflag>0";
  TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  if(useMRfirst) thecut = thecut && MRcut;
  TCut thecut0;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();

  // for(Int_t i=0; i<2; i++){
  for(Int_t i=0; i<num; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    tup[i] = GetTree(namestr[i]);
    // tup[i]->SetAlias("toffo","(tof>9.e3)*(tof)+(tof<9.e3)*(tof+40.e3)");
    tup[i]->SetAlias("toffo","(tof>9.3e3)*(tof)+(tof<9.3e3)*(tof+40.e3)"); // edited based on suggestion by KM on the August 3rd

    //    tup[i]->SetAlias("toffo","tof");
    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
    else kp[i] = tup[i]->GetMaximum("kp");
    cout << kp[i]<<endl;
    hx[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    //hx2[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);

    hx2[i] = new TH1F(Form("hx2%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    //hx2[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    hlambda2[i] = new TH1F(Form("hlambda2%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq02[i] = new TH1F(Form("hq02%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq2[i] = new TH1F(Form("hq2%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq3[i] = new TH1F(Form("hq3%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq4[i] = new TH1F(Form("hq4%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hxylambda2[i] = new TH3F(Form("hxylambda2%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);
    tup[i]->Draw(Form("x*%f>>hx2%d",range,i), thecut,"goff");
    
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda2%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda2%d",lambda_coeff,i), thecut && cut_ref2,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq02%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq2%d",twopirad,lambda_coeff,i), thecut && cut_ref2,"goff");
    tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda2%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    //tup[i]->Draw(Form("x*%f>>hx2%d",range,i), thecut2,"goff");
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq4%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    
    
    if(i!=6){
      leg->AddEntry(hx[i],degstr[i],"l");
    }
    

    hx[i]->Scale(25./kp[i]);
    hlambda[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq4[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    hxylambda[i]->Scale(25./kp[i]);

    hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    hratio[i]->Divide(hlambda[0]);
    hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    
    hq[i]->GetYaxis()->SetTitle("Reflectivity");

    hx2[i]->Scale(25./kp[i]);
    hlambda2[i]->Scale(25./kp[i]);
    hq2[i]->Scale(25./kp[i]);
    hq02[i]->Scale(25./kp[0]);
    hxylambda2[i]->Scale(25./kp[i]);

    hratio2[i]=(TH1F*)hlambda2[i]->Clone(Form("hratio%d",i));
    hratio2[i]->Divide(hlambda2[0]);
    hratio2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq2[i]->Divide(hq02[i]);
    hq[i]->Divide(hq0[i]);
    hq4[i]->Divide(hq0[i]);
    hq3[i]->Add(hq2[i], hq[i],1., 1.);
    //hq4[i]=hq[i];
    //hq4[i]->Divide(hq3[i]);
    



    Double_t E1[nbin_q];//={0.001};
    if(i==2){
      double AA[nbin_q];
      double AAE[nbin_q];
      double BB[nbin_q];
      double BBE[nbin_q];
      for(Int_t i1=0; i1<nbin_q; i1++){
        AA[i1]=hq4[i]->GetBinContent(i1);
        AAE[i1]=hq4[i]->GetBinError(i1);
        BB[i1]=hq2[i]->GetBinContent(i1);
        BBE[i1]=hq2[i]->GetBinError(i1);
        
        //E1[i1]=sqrt(pow(AA[i1]*BBE[i1],2)/pow((AA[i1]+BB[i1]),4)+pow(BB[i1]*AAE[i1],2)/pow((AA[i1]+BB[i1]),4));
        
        E1[i1]=sqrt(pow(AA[i1]*BBE[i1],2)+pow(BB[i1]*AAE[i1],2))/pow((AA[i1]+BB[i1]),2);
        
        //cout<<"AA_"<<AA[i1]<<"_BB_"<<BB[i1]<<endl;
        //cout<<"AA_"<<AAE[i1]<<endl;
        cout<<"AAE_"<<AAE[i1]<<"_BBE_"<<BBE[i1]<<"_E1_"<<E1[i1]<<endl;
        //cout<<"AA_"<<AA[i1]+BB[i1]<<endl;
        //cout<<"E1_"<<E1[i1]<<endl;

        hq[i]->SetError(E1);
      }
    }

    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq2[i]->Add(hq2[i], hq[i],1., 1.);

    hq[i]->Divide(hq3[i]);
    hq2[i]->Divide(hq3[i]);

    if(i==9){
      hx[i]->SetLineColor(i+2);
      hlambda[i]->SetLineColor(i+2);
      hratio[i]->SetLineColor(i+2);
      hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);

      hx2[i]->SetLineColor(i+2);
      hlambda2[i]->SetLineColor(i+2);
      hratio2[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hq02[i]->SetLineColor(i+2);
    } else {
      hx[i]->SetLineColor(i+1);
      hlambda[i]->SetLineColor(i+1);
      hratio[i]->SetLineColor(i+1);
      hq[i]->SetLineColor(i+1);
      hq0[i]->SetLineColor(i+1);
      hx2[i]->SetLineColor(i+1);
      hlambda2[i]->SetLineColor(i+1);
      hratio2[i]->SetLineColor(i+1);
      hq2[i]->SetLineColor(i+1);
      hq02[i]->SetLineColor(i+1);
    }
    c1->cd(1);
    //if(i==0)hx[i]->Draw("eh");
    //else hx[i]->Draw("ehsames");

    
    //TLine *l2 = new TLine (xend,1e-3,xend, 1e3);
    
    if(i==0)hx[i]->Draw("eh");
    
    else if(i!=6)hx[i]->Draw("ehsames");
    

    TLine *l1 = new TLine (xbegin,1e-3,xbegin, 1e3);
    TLine *l2 = new TLine (xcenter,1e-3,xcenter, 1e3);
    l1->SetLineColor(6);
    l2->SetLineColor(6);
    l1->Draw("ehsames");
    l2->Draw("ehsames");
    
    TLine *l3 = new TLine (xend,1e-3,xend, 1e3);
    TLine *l4 = new TLine (xcenter1,1e-3,xcenter1, 1e3);
    l3->SetLineColor(1);
    l4->SetLineColor(1);
    l3->Draw("ehsames");
    l4->Draw("ehsames");

    //if(i==7)hx[i]->Draw("ehsames");
    leg->Draw();

    c1->cd(2);
    //if(i==0)hlambda[i]->Draw("eh");
    //else hlambda[i]->Draw("ehsames");
    //if(i==0)hlambda[i]->Draw("eh");
    //if(i==2)hlambda[i]->Draw("ehsames");
    if(i==0)hlambda[i]->Draw("eh");
    else if(i!=6)hlambda[i]->Draw("ehsames");
    //if(i==7)hlambda[i]->Draw("ehsames");
    leg->Draw();

    c1->cd(3);
    //if(i==2)hq[i]->Draw("eh");
    f0->SetParLimits(0,0.11,0.15);
    f0->FixParameter(1,1.);
    //f0->SetParLimits(1,0.9,1.0);

    f0->SetParLimits(2,1.,10.);
    //hq[i]->Fit("f0","R","10000",0.173,0.5);
    f0->SetNpx(10000);
    //f0->SetParLimits(1,1.8,1.999);

    if(i==2) hq[i]->Fit(f0,"+","",0.1,0.47);
    

    //if(i==2) hratio[i]->Draw("eh");
    //if(i==7) hratio[i]->Draw("ehsames");
    //leg->Draw();

    c1->cd(4);
    if(i==2)hq4[i]->Draw("eh");
    
    
    hratio[i]->GetYaxis()->SetRangeUser(0.,2.);
    //hq[i]->GetYaxis()->SetRangeUser(0.1,0.9);
    //hq[i]->GetYaxis()->SetRangeUser(0.,2.);
    //hq[i]->GetXaxis()->SetRangeUser(0.1,0.9);

    hq[i]->GetYaxis()->SetRangeUser(0.,1.2);
    hq[i]->GetXaxis()->SetRangeUser(0.1,0.6);

    hq2[i]->GetYaxis()->SetRangeUser(0.,1.2);  
    hq2[i]->GetXaxis()->SetRangeUser(0.1,0.6);

    

  }

  hratio[1]->GetYaxis()->SetRangeUser(0.,2.);
  hq[1]->GetYaxis()->SetRangeUser(0.,1.);

  


  c1->cd(1); gPad->SetGrid();
  c1->cd(2); gPad->SetGrid(); gPad->SetLogy();
  c1->cd(3); gPad->SetGrid();
  c1->cd(4); gPad->SetGrid();

  c1->SaveAs(path_R+"q_R_I_on.png");
  c1->SaveAs(path_R+"q_R_I_on.root");

#if 1
  TFile *outfile = TFile::Open(path_R+"FeMirrorhist.root","RECREATE");
  for(Int_t i=0; i<num; i++){
    hx[i]->Write();
    hlambda[i]->Write();
    hratio[i]->Write();
    hq[i]->Write();
    hxylambda[i]->Write();
  }
  outfile->Close();
#endif


  return 0 ;
}
