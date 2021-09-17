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
const double mu_n=60.3;//[neV T^-1] 
const double mu_n2=-60.3;//[neV T^-1] 

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

    if(E1>V_Fe_p){
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
  TF1 *f0 = new TF1("",func_R0,0.01,0.5,2);
  
  //E>V1(-V)
  double func_R1(double *qqq,double *par){
    double q1=qqq[0];
    double E1=pow(hbar*q1,2)/8./m_nc2nm;
    double d1=par[0];
    double in_mag=par[1];
    double R1;

    if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n2*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n2*in_mag)-E1))/hbar;//nm^-1
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
  TF1 *f1 = new TF1("",func_R1,0.01,0.5,2);

  //E>V1(-V)
  double func_R2(double *qqq,double *par){
    double q1=qqq[0];
    double E1=pow(hbar*q1,2)/8./m_nc2nm;
    double d1=par[0];
    double in_mag=par[1];
    double R1;

    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      //double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k0-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k0)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k0)*sin(k1a*d1),2);
      R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        //double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*k0+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k0)*cosh(alpha*d1),2)+pow((k0*k0-alpha*alpha)*sinh(alpha*d1),2);
        R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }

    return R1;

  }
  TF1 *f2 = new TF1("",func_R2,0.01,0.5,2);


Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 666.; //sample to detector [mm]
Double_t xdirect    = 63.29;
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

//9/8
TTree* GetTree1(TString ROOTstr_path1){

  TFile *file = TFile::Open(ROOTstr_path1.Data());
  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}
////////

const string scan_id = "90nm_scan_fine_3";
const string run_id="20210717002421";
const Int_t num1 = 5; // this should be the half of the number of the files obtained by the scan 


Int_t P90nm_BG6(){
  // load CSV file with magnetic field
  vector<Int_t> vec_index;
  vector<Double_t> vec_I, vec_H; 
  // vec_index.clear(); vec_I.clear(); vec_H.clear();
  ifstream fcsv(Form("data_scans/magnetic_%s.csv", run_id.c_str())); 
  if (!fcsv.is_open())
  {
      exit(EXIT_FAILURE);
  }
  string str;
  getline(fcsv, str);
  cout << str << endl; // print out the first row
  while (getline(fcsv, str)){
      Int_t csv_i;     
      Double_t current;
      Double_t magfield;
      sscanf(str.c_str(), "%d,%lf,%lf", &csv_i, &current, &magfield);
      // str >> index >> ",'" >> current >>",'">> magfield;
      vec_index.push_back(csv_i); 
      vec_I.push_back(current);
      vec_H.push_back(magfield);
      // cout << " " << index
    //  << " " << current
      //   << " " << magfield << endl;
      // }
  // i_csv++;
  }

  InitColor();
  TH1::SetDefaultSumw2();

  const Int_t num = 4;


  Int_t kp11[num1];
  Int_t kp211[num1];
  TTree* tup11[num1];
  TTree* tup211[num1];
  //TTree* tup02[num1];

  TH1F* hx11[num1];
  TH1F* hlambda11[num1];
  TH1F* hratio11[num1];
  TH1F* hq11[num+5];
  TH1F* hq011[num+5];
  TH3F* hxylambda11[num1];


  TTree* tup0;
  Int_t kp0;


  Double_t angle11[num+5];
  //Double_t angle211[num1];
  for(Int_t i=0; i< num; i++){
    angle11[i]=TMath::Abs(47.1868 - xdirect)/dist_det;
    //angle211[i]=TMath::Abs(47.1868 - xdirect)/dist_det;
  }  



  TString namestr[num];
  TString namestr2[num];
  
  TString degstr[num];
  TString degstr2[num];
  //off
  degstr[0]="Direct(M1 reflect)";
  degstr[1]="B = 8.13 mT";
  degstr[2]="B = 0.322 mT";
  //degstr[3]="B = 0.908 mT";
  degstr[3]="B = 1.35 mT";
  //degstr[5]="B = 1.80 mT";
  //degstr[6]="B = 2.66 mT";

  //on
  degstr2[0]="Direct(M1 reflect)";
  degstr2[1]="B = 8.13 mT";
  degstr2[2]="B = 0.322 mT";
  //degstr2[3]="B = 0.908 mT";
  degstr2[3]="B = 1.35 mT";
  //degstr2[5]="B = 1.80 mT";
  //degstr2[6]="B = 2.66 mT";
  

  Double_t angle[num];
  Double_t angle2[num];
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(47.4493 - xdirect)/dist_det; //rad 47.1868(30nm)
  //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(47.3999 - xdirect)/dist_det; //rad 47.04
  angle[3] = TMath::Abs(47.4013 - xdirect)/dist_det; //rad 47.2
  //angle[4] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[5] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[6] = TMath::Abs(47.2 - xdirect)/dist_det; //rad

  angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle2[1] = TMath::Abs(47.409 - xdirect)/dist_det; //rad 47.07
  angle2[2] = TMath::Abs(47.479- xdirect)/dist_det; //rad 47.1 
  angle2[3] = TMath::Abs(47.4845 - xdirect)/dist_det; //rad 47.19
  //angle2[4] = TMath::Abs(47.24 - xdirect)/dist_det; //rad
  //angle2[5] = TMath::Abs(47.21 - xdirect)/dist_det; //rad
  //angle2[6] = TMath::Abs(47.11 - xdirect)/dist_det; //rad
  
  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm OFF");
  TLegend* leg2 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm ON");
  TLegend* leg3 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 60;//300 60
  Double_t q_min  = 0.1;//0.6
  Double_t q_max  = 0.50;//0.6
  //Double_t q_max  = 1.0;//0.6
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;

  Double_t xbegin1=40.;
  Double_t xcenter1=55.;

  Double_t xcenter=55.;
  Double_t xend=71.;
  //Double_t xcenter=55.;
  //Double_t xend=71.;
  
  //Double_t ybegin1=43.;
  //Double_t yend1=60.;
  Double_t ybegin1=85.;
  Double_t yend1=102.;
  //Double_t ybegin1=55.;
  //Double_t yend1=71.;

  
  Double_t ybegin=55.;
  Double_t yend=71.;

  //  Double_t ybegin=70.;
  //  Double_t yend=77.;

  TString degstr11[num1];
  TString degstr211[num1];
  TString degstr_p11[num1];

  for(Int_t i=0; i< num+5; i++){
    degstr11[i]=Form("SF:ON,  B=%lf mT", vec_H[i]);
    //degstr211[i]=Form("SF:OFF, B=%lf mT", vec_H[i]);
    //degstr_p11[i]=Form("B=%lf mT", vec_H[i]);

  } 

  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     LLD,LLD,LLD,LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  TCut cut_x = Form("x*%f>20 && x*%f<100",range,range);
  TCut cut_y = Form("y*%f>%f && y*%f<%f",range,ybegin,range,yend);
  TCut cut_y1 = Form("y*%f>%f && y*%f<%f",range,ybegin1,range,yend1);
  
  TCut cut_dir = Form("x*%f>%f && x*%f<%f",range,xcenter,range,xend);
  TCut cut_ref = Form("x*%f>%f && x*%f<%f",range,xbegin1,range,xcenter);
  TCut cut_ref1 = Form("x*%f>%f && x*%f<%f",range,xbegin1,range,xcenter1);
  //  TCut cut_tof = Form("tof>1.0e3 && tof<39.9e3");
  TCut cut_tof = "";
  TCut MRcut = "MRflag>0";
  //TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  TCut thecut1 = cut_rpmt_basic && cut_x && cut_y1 && cut_tof;
  TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  if(useMRfirst) thecut = thecut && MRcut;
  TCut thecut0;
  TCut thecut01;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,2);
  c1->cd(1);

  TString namestr11[num+5];
  //TString namestr211[num+5];
  
  TString namestr_ref= "data/210713_SiFe/20210714193654_list.root"; // file path of the direct data
  
  for(Int_t i=0; i<num+5; i++){
    int iscan=i;
    namestr11[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  /*
  for(Int_t i=0; i<num1; i++){
    int iscan=i*2+1;
    namestr211[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  */

  tup0 = GetTree1(namestr_ref); // direct data 
    tup0->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); 
    if(useMRfirst) kp0 = tup0->GetMaximum("mp");
      else kp0= (tup0->GetMaximum("kp") - tup0->GetMinimum("kp"));
  cout << "direct data # of kp: "<< kp0 <<endl;


  // for(Int_t i=0; i<2; i++){

  for(Int_t i=0; i<num+5; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;
    if(i==0) thecut01=thecut1;

    Double_t twopirad = 2*TMath::Pi()*angle[0];
    Double_t twopirad2 = 2*TMath::Pi()*angle2[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    Double_t twopirad11 = 2*TMath::Pi()*angle11[i];
 
    // Direct data

    hq011[i] = new TH1F(Form("hq011%d",i),Form("%s;q [nm^{-1}];count/bin/25kp","Direct"),nbin_q,q_min,q_max);
    tup0->Draw(Form("%f/(toffo*%f)>>hq011%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    hq011[i]->Scale(25./kp0);

     // Data with SF:ON
    tup11[i] = GetTree1(namestr11[i]);
    tup11[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    if(useMRfirst) kp11[i] = tup11[i]->GetMaximum("mp");
      else kp11[i] = (tup11[i]->GetMaximum("kp") - tup11[i]->GetMinimum("kp"));
    cout << "SF:ON data # of kp: "<< kp11[i] <<endl;

    hq11[i] = new TH1F(Form("hq11%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr11[i].Data()),nbin_q,q_min,q_max);
    // tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup11[i]->Draw(Form("%f/(toffo*%f)>>hq11%d",twopirad11,lambda_coeff,i), thecut && cut_ref,"goff");
    hq11[i]->Scale(25./kp11[i]);
    hq11[i]->Divide(hq011[i]); //Divide by the direct data
    hq11[i]->GetYaxis()->SetTitle("Reflectivity");
    hq11[i]->SetTitle("Reflectivity (SF ON)");
    //leg->AddEntry(hq1[i],degstr[i],"l");

  }
  const Int_t num_pol = 3;
  TGraphErrors* gr[num_pol];
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  //Double_t q_cuts[num_pol] = {0.25, 0.3, 0.35};
  
  //Double_t B[num-1] = {8.01,0.322,1.35};
  Double_t B11[num+5] = {1.96615,2.0101,2.05405,2.098,2.14195,2.1859,2.22985,2.2738,2.31775};

  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  //Double_t pol_at_qcut[num-1];
  //Double_t error_pol_at_qcut[num-1];

  Double_t ibin_pol[num_pol];

  //Int_t num11=10;

  Double_t q_cuts11[num_pol] = {0.25, 0.3, 0.35};
  
  //Double_t B11[num-1] = {8.01,0.322,1.35};
  //Double_t B11[num11-1]={1.96615,2.0101,2.05405,2.098,2.14195,2.1859,2.22985,2.2738,2.31775};
  
  // ={1,2,3,4,5,6,7,8,9};
  //B11[num11-1]={1.96615,2.0101,2.05405,2.098,2.14195,2.1859,2.22985,2.2738,2.31775};
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t pol_at_qcut11[num+5];
  Double_t error_pol_at_qcut11[num+5];

  //Double_t ibin_pol11[num_pol];

  c1->cd(4); 

  ofstream ofs(path_R+Form("B_Rdata_all_OFF_%s.csv", scan_id.c_str()));  // ファイルパスを指定する


  for (Int_t i=0; i<num_pol; i++){
    ibin_pol[i]= Int_t((q_cuts11[i]-q_min)*nbin_q/(q_max-q_min));
    Double_t q_cut_i = q_min + (q_max-q_min)*ibin_pol[i]/nbin_q;
    /*
    for (Int_t j=0; j<num-1; j++){
      pol_at_qcut[j] = hq[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut[j] = hq[j]->GetBinError(ibin_pol[i]);
      ofs << q_cut_i << "," << B[j] << ","<< pol_at_qcut[j] << ","<< error_pol_at_qcut[j] << endl;
      // cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;

    }
    */
    for (Int_t j=0; j<num+5; j++){
      pol_at_qcut11[j] = hq11[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut11[j] = hq11[j]->GetBinError(ibin_pol[i]);
      ofs << q_cut_i << "," << B11[j] << ","<< pol_at_qcut11[j] << ","<< error_pol_at_qcut11[j] << endl;
    }
  

    /*
    ofstream ofs1(path_R+Form("B_RdataOFF_scan_%s.csv", scan_id.c_str()));  // ファイルパスを指定する
  
    //for (Int_t i=0; i<num_pol; i++){
    ibin_pol11[i]= Int_t((q_cuts11[i]-q_min)*nbin_q/(q_max-q_min));
    Double_t q_cut_i11 = q_min + (q_max-q_min)*ibin_pol11[i]/nbin_q;
    for (Int_t j=0; j<num11; j++){
      pol_at_qcut11[j] = hq11[j]->GetBinContent(ibin_pol11[i]);
      error_pol_at_qcut11[j] = hq11[j]->GetBinError(ibin_pol11[i]);
      ofs1 << q_cut_i11 << "," << B11[j] << ","<< pol_at_qcut11[j] << ","<< error_pol_at_qcut11[j] << endl;
      // cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }
    */
  }



  return 0 ;
}
