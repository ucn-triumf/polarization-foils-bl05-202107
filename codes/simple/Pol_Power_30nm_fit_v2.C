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
Double_t dist_det   = 666.; //sample to detector [mm]
Double_t xdirect    = 63.29;
Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.




//8/29 add
const double c=299792458; //[m/s]
const double c_nm=299792458.e9; //[nm/s]
const double m1=939.5654133; //[MeV/c^2]
const double m2=939.5654133e15; //[neV/c^2]
//Double_t m=939.5654133e9/(pow(c,2)); //[meV/(m/s)^2]
double m_nc2=939.5654133e15/(pow(c,2)); //[neV/(m/s)^2]
double m_nc2nm=m2/(pow(c_nm,2)); //[neV/(nm/s)^2]


const double h1=4.135667696e-15; //[eV.s]
const double h=4.135667696e-12; //[meV.s]
const double hbar=4.135667696e-6/(2.*TMath::Pi());//neV.s

const double V_Fe=209.0602;//neV
const double V_Si=54.0078;//neV
// const double mu_n=60.3;//[neV T^-1] 
const double mu_n=-60.3;//[neV T^-1] 

const double in_T=2.;//T(tesla)
const double V_Fe_p=V_Fe+mu_n*in_T;//neV
const double V_Fe_m=V_Fe-mu_n*in_T;//neV
  /*double k1b=sqrt(2*m_nc2*(E1-(V_Fe-mu_n*in_T)))/hbar;//m^-1
  double k2=sqrt(2*m_nc2*(E1-(V_Si))/hbar;//m^-1
  double alpha1=k1b=sqrt(2*m_nc2*((V_Fe-mu_n*in_T)-E1))/hbar;
  double alpha2=k1b=sqrt(2*m_nc2*((V_Fe+mu_n*in_T)-E1))/hbar;
*/
  //E>V1(+V)
double func_R0(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double R1;
  double V_div_p=(V_Fe+mu_n*in_mag);

    if(E1>V_div_p){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        R1=1.;
      }
    }

    return R1;

  }
  TF1 *f0 = new TF1("",func_R0,0.05,0.5,2);
  
  //E>V1(-V)
  double func_R1(double *qqq,double *par){
    double q1=qqq[0];
    double E1=pow(hbar*q1,2)/8./m_nc2nm;
    double d1=par[0];
    double in_mag=par[1];
    double R1;
    double V_div_m=(V_Fe-mu_n*in_mag);

    //if(E1>V_Fe_m){
    if(E1>V_div_m){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe-mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      //if(E1>V_Si){
      if(E1>V_Si){  
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe-mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        R1=1.;
      }
    }

    return R1;

  }
  TF1 *f1 = new TF1("",func_R1,0.05,0.5,2);


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

//8/29 add BG

//TString path_R = "results/"; // path to the results directory 
TString path_D = "data/210713_SiFe/"; // path to the data directory 
//TString rootfile  = "20210714185238_list.root"; // name of the target root file
//TString rootfile  = "_list.root"; // name of the target root file
//TString rootfile  = "20210715072653_list.root"; // name of the target root file

/*
  Double_t xbegin=40.;
  Double_t xcenter=55.;
  Double_t xend=71.;
  Double_t ybegin=65.;
  Double_t yend=82.;
*/
///*
Double_t x_cut_low = 40; // for transmission wave 
Double_t x_cut_up =  55; // for transmission wave
Double_t x_cut_low3 = 25; // for transmission wave 
Double_t x_cut_up3 =  55; // for transmission wave
Double_t y_cut_low3 = 85;
Double_t y_cut_up3 = 102;
Double_t y_cut_low = 65;
Double_t y_cut_up = 82;
//*/
/*
Double_t x_cut_low = 0; // for transmission wave 
Double_t x_cut_up =  120; // for transmission wave
Double_t y_cut_low = 0;
Double_t y_cut_up = 120;
*/
Double_t range=128;

TCut cut_xy =Form("x*%f>%f && x*%f<%f && y*%f>%f && y*%f<%f && f==4", range, x_cut_low,range,x_cut_up,range,y_cut_low, range,y_cut_up);


const Int_t nBinXY = 640;
const Double_t startX = 0;
const Double_t endX = 128;
const Double_t startY = 0;
const Double_t endY = 128;
const Int_t nBinCut = 200;


Int_t Pol_Power_30nm_fit_v2(){

  InitColor();
  TH1::SetDefaultSumw2();

  const Int_t num = 7;
  Int_t kp[num];
  Int_t kp2[num];
  TTree* tup[num];
  TTree* tup2[num];
  TH1F* hx[num];
  TH1F* hlambda[num];
  TH1F* hratio[num];
  TH1F* hq[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TH1F* hpolratio[num];
  TH1F* hpolratio2[num];
  TH1F* hpolratio3[num];
  TH1F* hq2[num];
  TH1F* hq02[num];

  TString namestr[num];
  TString namestr2[num];
  //off
  namestr[0]="20210714193654_list.root";
  namestr[1]="20210715072653_list.root"; //1.97A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 8 mT
  //namestr[2]="20210715075452_list.root"; //0A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT
  namestr[2]="20210715081447_list.root"; //0A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = -8 mT -> 0 mT
  namestr[3]="20210715084835_list.root"; //0.15A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1 mT from -8 mT
  namestr[4]="20210715085349_list.root"; //0.264A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT
  namestr[5]="20210715082606_list.root"; //0.378A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT  -> 2 mT
  namestr[6]="20210715083711_list.root"; //0.6A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 2 mT  -> 3 mT

  //on
  namestr2[0]="20210714193654_list.root"; //M1 reflect (direct) 1hour
  namestr2[1]="20210715073913_list.root";
  namestr2[2]="20210715082018_list.root";
  namestr2[3]="20210715085144_list.root"; //Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1 mT from -8 mT  with AFP 760 mV
  namestr2[4]="20210715085714_list.root"; //Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT  with AFP 760 mV
  namestr2[5]="20210715083141_list.root"; 
  namestr2[6]="20210715084052_list.root";
 
  TString degstr[num];
  TString degstr2[num];
  //off
  degstr[0]="Direct(M1 reflect)";
  degstr[1]="B = 8.01 mT";
  degstr[2]="B = 0.322 mT";
  degstr[3]="B = 0.908 mT";
  degstr[4]="B = 1.35 mT";
  degstr[5]="B = 1.80 mT";
  degstr[6]="B = 2.66 mT";

  //on
  degstr2[0]="Direct(M1 reflect)";
  degstr2[1]="B = 8.01 mT";
  degstr2[2]="B = 0.322 mT";
  degstr2[3]="B = 0.908 mT";
  degstr2[4]="B = 1.35 mT";
  degstr2[5]="B = 1.80 mT";
  degstr2[6]="B = 2.66 mT";
  

  Double_t angle[num];
  Double_t angle2[num];
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(47.1868 - xdirect)/dist_det; //rad 47.1868
  //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(47.04 - xdirect)/dist_det; //rad
  angle[3] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  angle[4] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  angle[5] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  angle[6] = TMath::Abs(47.2 - xdirect)/dist_det; //rad

  angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle2[1] = TMath::Abs(47.07 - xdirect)/dist_det; //rad
  angle2[2] = TMath::Abs(47.1 - xdirect)/dist_det; //rad
  angle2[3] = TMath::Abs(47.19 - xdirect)/dist_det; //rad
  angle2[4] = TMath::Abs(47.24 - xdirect)/dist_det; //rad
  angle2[5] = TMath::Abs(47.21 - xdirect)/dist_det; //rad
  angle2[6] = TMath::Abs(47.11 - xdirect)/dist_det; //rad

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
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 30 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 30 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 30 nm (OFF-ON)/(OFF+ON)");
  // TLegend* leg = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm OFF");
  // TLegend* leg2 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm ON");
  // TLegend* leg3 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.8, 0.5, 1.0, 1,"");
  TLegend* leg2 = new TLegend(0.8, 0.5, 1.0, 1,"");
  TLegend* leg3 = new TLegend(0.8, 0.5, 1.0, 1,"");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 60;//300
  Double_t q_min  = 0.05;//0.6 
  Double_t q_max  = 0.50;//0.6
  //Double_t q_max  = 1.0;//0.6
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;
  //  Double_t HLD  = 7400.;

  //  Double_t xbegin=54.;
  Double_t xbegin=40.;
  Double_t xcenter=55.;
  Double_t xend=71.;
  Double_t ybegin=65.;
  Double_t yend=82.;
  //  Double_t ybegin=70.;
  //  Double_t yend=77.;



  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     LLD,LLD,LLD,LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  TCut cut_x = Form("x*%f>20 && x*%f<100",range,range);
  TCut cut_y = Form("y*%f>%f && y*%f<%f",range,ybegin,range,yend);
  TCut cut_dir = Form("x*%f>%f && x*%f<%f",range,xcenter,range,xend);
  TCut cut_ref = Form("x*%f>%f && x*%f<%f",range,xbegin,range,xcenter);
  //  TCut cut_tof = Form("tof>1.0e3 && tof<39.9e3");
  TCut cut_tof = "";
  TCut MRcut = "MRflag>0";
  TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  if(useMRfirst) thecut = thecut && MRcut;
  TCut thecut0;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,2);
  c1->cd(1);


  


  // for(Int_t i=0; i<2; i++){


  for(Int_t i=0; i<num; i++){

    // 8/29 add BG
    /*
  TString rootfile[num];
  TString rootfile_num[num];
  rootfile[i]= namestr[i];
  TString rootfile_num    = path_D + rootfile;
  
  const TString tree_name  = "T";
  const TString cut_str    = "f==4";

  TFile *rootfile  = TFile::Open(rootfile_num.Data());
  TTree *tree = rootfile->Get<TTree>(tree_name);
  const Double_t time = (tree->GetMaximum("kp")-tree->GetMinimum("kp"))/25.;
  TCanvas *c_xy = new TCanvas("c_xy","c_xy", 800, 1000);
  //c_xy->Divide(1,3);
  TH2D *h_xy = new TH2D("h_xy","XY from upstream",nBinXY,startX,endX,nBinXY,startY,endY);
  TH1D *h_x = new TH1D("h_x","X histogram",nBinCut,x_cut_low,x_cut_up);
  TH1D *h_y = new TH1D("h_y","Y histogram",nBinCut,y_cut_low,y_cut_up);
  TH2D *h_xy2 = new TH2D("h_xy2","XY from upstream2",nBinXY,x_cut_low,x_cut_up,nBinXY,y_cut_low,y_cut_up);
  TH2D *h_xy3 = new TH2D("h_xy3","XY from upstream3",nBinXY,x_cut_low3,x_cut_up3,nBinXY,y_cut_low3,y_cut_up3);

  //c_xy->cd(1);
  tree->Draw(Form("y*%f:x*%f>>h_xy2", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  tree->Draw(Form("y*%f:x*%f>>h_xy3", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  tree->Draw(Form("y*%f:x*%f>>h_xy", range, range), cut_rpmt*TCut(Form("%f",1./time)), "colz");  // change to count rate (cps)
  
  const Double_t sum  = h_xy->GetSumOfWeights();
  const Double_t sum2  = h_xy2->GetSumOfWeights();
  const Double_t sum3  = h_xy3->GetSumOfWeights();
  const Double_t  sum32=sum3/2.;
*/
 

    thecut.Print();
    if(i==0) thecut0=thecut;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t twopirad2 = 2*TMath::Pi()*angle2[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    tup[i] = GetTree(namestr[i]);
    // tup[i]->SetAlias("toffo","(tof>9.e3)*(tof)+(tof<9.e3)*(tof+40.e3)");
    tup[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    
    tup2[i] = GetTree(namestr2[i]);
    tup2[i]->SetAlias("toffo2","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // edited based on suggestion by KM on the August 3rd
    

    //tup[i]->SetAlias("toffo","tof");
    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
    else kp[i] = tup[i]->GetMaximum("kp");
    cout << kp[i]<<endl;

    //hx[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    //hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    // hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);
  
    //tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    //if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    //else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    // tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    
    //hpolratio[i]->Rebin(10);
    //hpolratio2[i]->Rebin(10);

    //hpolratio[i]->Divide(hpolratio2[i]);


    leg->AddEntry(hq[i],degstr[i],"l");
    leg2->AddEntry(hq[i],degstr[i],"l");
    if(i!=0)leg3->AddEntry(hq[i],degstr[i],"l");
    //leg3->AddEntry(hq[i],degstr[i],"l");
    //leg->AddEntry(hq2[i],degstr2[i],"l2");

    //hx[i]->Scale(25./kp[i]);
    //hlambda[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    //hxylambda[i]->Scale(25./kp[i]);

    //BG 8/29add
    double BG = -0.22/nbin_q ;
    hq[i]->Add(hq[i],BG);
    hq0[i]->Add(hq0[i],BG);
    
    
    
    
    //hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    //hratio[i]->Divide(hlambda[0]);
    //hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    

    hq[i]->Divide(hq0[i]);
    hq[i]->GetYaxis()->SetTitle("Reflectivity");
    hq[i]->SetTitle("Reflectivity (SF OFF)");
    //hq2[i]->Divide(hq02[i]);

    if(useMRfirst) kp2[i] = tup2[i]->GetMaximum("mp");
    else kp2[i] = tup2[i]->GetMaximum("kp");
    //cout << kp2[i]<<endl;
    hq02[i] = new TH1F(Form("hq02%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    hq2[i] = new TH1F(Form("hq2%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    
    tup2[0]->Draw(Form("%f/(toffo2*%f)>>hq02%d",twopirad2,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup2[i]->Draw(Form("%f/(toffo2*%f)>>hq2%d",twopirad2,lambda_coeff,i), thecut && cut_ref,"goff");
    
    hq2[i]->Scale(25./kp2[i]);
    hq02[i]->Scale(25./kp2[0]);
    //BG 8/29add
    //double BG = -0.22/nbin_q ;
    hq2[i]->Add(hq2[i],BG);
    hq02[i]->Add(hq02[i],BG);
    

    hq2[i]->Divide(hq02[i]);

    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq2[i]->SetTitle("Reflectivity (SF ON)");
    //hq02[i]->GetYaxis()->SetTitle("Reflectivity02");

    hpolratio[i]=(TH1F*)hq[i]->Clone(Form("hpolratio%d",i));
    hpolratio[i]->Add(hq[i], hq2[i],1., -1.); // hq: OFF, hq2: ON, calculate Non - Noff
    // hpolratio[i]->Add(hq2[i],-1.);//各binの値=x*hist1の値+y*hist2の値
    hpolratio2[i]=(TH1F*)hq[i]->Clone(Form("hpolratio2%d",i));
    hpolratio2[i]->Add(hq[i], hq2[i],1., 1.); // hq: OFF, hq2: ON, calculate Non + Noff
    hpolratio[i]->Divide(hpolratio2[i]);
    hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{off}-R_{on})/(R_{off}+R_{on})");
    hpolratio[i]->SetTitle("Polarization power");

    // double pol_array[num];
    // //hq[i]->FindBin(nbb[i]);
    // double qq=0.25;
    // double nbin_qq=qq*nbin_q/(q_max-q_min);
    // entry[i]= hpolratio[i]->GetBinContent(nbin_qq); 
    // cout<<"[i]_"<<i<<"entry_"<<entry[i]<<endl;

    // Double_t pol_q[7] = {0.,8.01,0.322,0.908,1.35,1.80,2.66};
    // Double_t q_cuts[5] = {0.25, 0.3, 0.35};

    

    // g1->SetMarkerStyle(22);
    // g1->SetMarkerColor(2);
    // g1->SetMarkerSize(1);
    // g1->GetXaxis()->SetTitle("B [mT]");
    // g1->GetYaxis()->SetTitle("Polarization power");



    if(i==9){
      //hx[i]->SetLineColor(i+2);
      //hlambda[i]->SetLineColor(i+2);
      //hratio[i]->SetLineColor(i+2);
      hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hpolratio[i]->SetLineColor(i+2);
      hpolratio2[i]->SetLineColor(i+2);
    } else {
      //hx[i]->SetLineColor(i+1);
      //hlambda[i]->SetLineColor(i+1);
      //hratio[i]->SetLineColor(i+1);
      hq[i]->SetLineColor(i+1);
      hq0[i]->SetLineColor(i+1);
      hq2[i]->SetLineColor(i+1);
      hpolratio[i]->SetLineColor(i+1);
      hpolratio2[i]->SetLineColor(i+1);
    }

    //hpolratio[i]->Rebin(10);
    //hpolratio[i]->Scale(10);
    //hpolratio3[i]->hpolratio[i]/10.;
/*
    c1->cd(1);
    if(i==1)hpolratio2[i]->Draw("eh");
    else hpolratio2[i]->Draw("ehsames");
    leg->Draw();
 */  
    double qm=0.5;
    double ymax=2.;//gPad->GetUymax();
    double ymin=0.;//gPad->GetUymin();
    double xmax1=0.2;
    double xmin1=0.4;
    c1->cd(1);
    if(i!=0){
      

      TBox* b = new TBox(0.287,0,qm,2); 
      b->SetFillColor( 7 ); 
      b->SetFillStyle(3004); 
      b->Draw(); 
      TBox* b1 = new TBox(0.15,0,0.173,2); 
      b1->SetFillColor( kOrange ); 
      b1->SetFillStyle(3004); 
      b1->Draw("same");
      leg->Draw();
      
      hq[i]->SetStats(0);  
      if(i==1)hq[i]->Draw("eh");
      else hq[i]->Draw("ehsames");
      // if(i==1)hq[i]->Draw("ah");
      // else hq[i]->Draw("ahsames");


      ///8/29 add
      f0->SetParameter(0.,30.e-9);//m //第１引数が変数の番号、第２引数がその値
      //f0->SetParameter(0.,30.);//nm 
      f0->SetParameter(1,2.);
      //hq[i]->Fit("f0","R","10000",0.173,0.5);
      f0->SetNpx(1000);
      //f0->SetParLimits(1,1.8,1.999);
      //
      f0->Draw("sames");

      f0->SetParameter(0.,30.e-9);//m //第１引数が変数の番号、第２引数がその値
      //f1->SetParameter(0.,30.);//nm 
      f1->SetParameter(1,2.);
      //f1->SetParLimits(1,1.8,1.999);
      f1->SetNpx(1000);
      f1->Draw("sames");

      if(i==1){
        hq[1]->Fit("f0","","",0.252,0.5);
      }
      if(i==2){
        hq[2]->Fit("f1","","",0.252,0.5);
      }
      

    
    }

    c1->cd(2);
    if(i!=0){
     

      TBox* b2 = new TBox(0.287,0,qm,2.); 
      b2->SetFillColor( 7  ); 
      b2->SetFillStyle(3004); 
      b2->Draw();
      TBox* b3 = new TBox(0.15,0,0.173,2.); 
      b3->SetFillColor( kOrange); 
      b3->SetFillStyle(3004); 
      b3->Draw("same"); 
      leg2->Draw();
      hq2[i]->SetStats(0);
     
     if(i==1)hq2[i]->Draw("eh");
     else hq2[i]->Draw("ehsames");
    // if(i==1)hq2[i]->Draw("ah");
    // else hq2[i]->Draw("ahsames");
    

    
    
    }
    

    c1->cd(3);
    if(i!=0){
      TBox* b4 = new TBox(0.15,-1.2,0.173,1.2); 
      b4->SetFillColor( kOrange ); 
      b4->SetFillStyle(3004); 
      b4->Draw();
      TBox* b5 = new TBox(0.287,-1.2,qm,1.2); 
      b5->SetFillColor(7); 
      b5->SetFillStyle(3004); 
      b5->Draw("same");
      hpolratio[i]->SetStats(0);
      if(i==1)hpolratio[i]->Draw("eh");
      else hpolratio[i]->Draw("ehsames");

    leg3->Draw();
    }

    c1->cd(4);
    if(i!=0){
    // g1->Draw("AP");
    //leg->Draw();

    
    }
/*
    c1->cd(3);
    if(i==1)hpolrhqtio[i]->Draw("eh");
    else hpolratio[i]->Draw("ehsames");
    leg->Draw();
*/
    
    hpolratio[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hpolratio[i]->GetYaxis()->SetRangeUser(-1.2,1.2);
    hpolratio2[i]->GetYaxis()->SetRangeUser(-1.2,1.2);
    hq[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq[i]->GetYaxis()->SetRangeUser(1.e-3,2.);
    hq2[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq2[i]->GetYaxis()->SetRangeUser(1.e-3,2.);
    // hq[i]->SaveAs(path_R + Form("hq_off_%d.root", i));
    // hq2[i]->SaveAs(path_R + Form("hq_on_%d.root", i));
    


    cout<<"in_"<<angledeg[i]<<"_deg"<<endl;

  }


  const Int_t num_pol = 3;
  TGraphErrors* gr[num_pol];
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t q_cuts[num_pol] = {0.25, 0.3, 0.35};
  
  Double_t B[num-1] = {8.01,0.322,0.908,1.35,1.80,2.66};
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t pol_at_qcut[num-1];
  Double_t error_pol_at_qcut[num-1];

  Double_t ibin_pol[num_pol];

  ofstream ofs(path_R+"30nm_mT_P.csv");  
  c1->cd(4); 
  TLegend *leg4 = new TLegend(0.8, 0.8, 0.95, 1, "");
  // leg4->SetFillStyle(0);
  for (Int_t i=0; i<num_pol; i++){
    ibin_pol[i]= Int_t((q_cuts[i]-q_min)*nbin_q/(q_max-q_min));
    for (Int_t j=0; j<num-1; j++){
      pol_at_qcut[j] = hpolratio[j+1]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut[j] = hpolratio[j+1]->GetBinError(ibin_pol[i]);
      // cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
      
      ofs << B[j] << " "<< pol_at_qcut[j] << " "<< error_pol_at_qcut[j] << endl;
    }
    gr[i]= new TGraphErrors(num-1, B, pol_at_qcut,0,error_pol_at_qcut);
    gr[i]->GetXaxis()->SetRangeUser(0.2, 9);
    gr[i]->GetYaxis()->SetRangeUser(-1.2, 1.2);
    if (i==0) gr[i]->Draw("AP");
    else gr[i]->Draw("P");
    gr[i]->SetMarkerColor(i+1);
    gr[i]->SetLineColor(i+1);
    gr[i]->SetMarkerStyle(i+3);
    gr[i]->SetMarkerSize(1);
    gr[i]->GetXaxis()->SetTitle("B (mT)");
    gr[i]->GetYaxis()->SetTitle("Polarization power");
    gr[i]->SetTitle("");
    leg4->AddEntry(gr[i],Form("q=%.3f nm^{-1}",q_min + (q_max-q_min)*ibin_pol[i]/nbin_q),"p");

    /*
    ofstream ofs(path_R+"30nm_mT_P.csv");  // ファイルパスを指定する
    for(Int_t i=0; i<num-1; i++){
      ofs << B[i] << " "<< pol_at_qcut[i] << " "<< error_pol_at_qcut[i] << endl;
      
    
    return 0;
    }*/
  }
  leg4->Draw();
    
    


  // double qq=0.25;
  // double nbin_qq=qq*nbin_q/(q_max-q_min);
  // entry[i]= hpolratio[i]->GetBinContent(nbin_qq); 
  // cout<<"[i]_"<<i<<"entry_"<<entry[i]<<endl;

  

  c1->cd(1); gPad->SetGrid();//gPad->SetLogy();
  c1->cd(2); gPad->SetGrid();//gPad->SetLogy();
  
  c1->cd(3); gPad->SetGrid();//gPad->SetLogy();
  c1->cd(4); gPad->SetGrid();//gPad->SetLogy();
    

  c1->SaveAs(path_R+"pol_30nm.png");
  c1->SaveAs(path_R+"pol_30nm.root");

#if 1
  TFile *outfile = TFile::Open(path_R+"pol.root","RECREATE");
  for(Int_t i=0; i<num; i++){
    //hx[i]->Write();
    //hlambda[i]->Write();
    //hratio[i]->Write();
    hq[i]->Write();
    //hxylambda[i]->Write();
    hq2[i]->Write();
  }
  outfile->Close();
#endif


  return 0 ;
}
