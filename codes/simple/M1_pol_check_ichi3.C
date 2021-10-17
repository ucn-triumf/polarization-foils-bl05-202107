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


// definition of shared parameter
// background function
int iparB[3] = { 1,      // exp amplitude in B histo
                 2,
                 3    // exp common parameter
};
 
// signal + background function
int iparSB[3] = { 1, //thickness
                  2, //B(T)
                  3  //V_Fe
         
};
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f11,
                ROOT::Math::IMultiGenFunction & f22) :
      fChi2_1(&f11), fChi2_2(&f22) {}
 
   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[3];
      for (int i = 0; i < 3; ++i) p1[i] = par[iparB[i] ];
 
      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparSB[i] ];
 
      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }
 
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};
 


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

const double V_Fe1=209.0602;//neV
const double V_Si=54.0078;//neV
const double mu_n=60.3;//[neV T^-1] 
const double mu_n2=-60.3;//[neV T^-1] 

const double in_T=2.;//T(tesla)
const double V_Fe_p1=V_Fe1+mu_n*in_T;//neV
const double V_Fe_m1=V_Fe1-mu_n*in_T;//neV
  /*double k1b=sqrt(2*m_nc2*(E1-(V_Fe-mu_n*in_T)))/hbar;//m^-1
  double k2=sqrt(2*m_nc2*(E1-(V_Si))/hbar;//m^-1
  double alpha1=k1b=sqrt(2*m_nc2*((V_Fe-mu_n*in_T)-E1))/hbar;
  double alpha2=k1b=sqrt(2*m_nc2*((V_Fe+mu_n*in_T)-E1))/hbar;
*/
  //E>V1(+V)
double qc=0.13;//par[2];
double ww=7.16389E-02;//par[3];
double R0=1.;//par[4];

//double mm=par[3];
double alpha=-8.03288E-02;//par[5];

double uprate=0.5;//par[3];
double downrate=0.5;//par[4];

double mm=5.;
double qc_up=0.22;
double mm2=6.24502;

double func_R0(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double V_Fe=par[2];
  double V_Fe_p=V_Fe+mu_n*in_mag;
  double V_Fe_m=V_Fe-mu_n*in_mag;
  double R1;
  //double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  //double qc=1.36239E-01;
  double qc=1.30030e-01;
  double mm2=5.88427e+00 ;

  double ww=2.5E-03;//par[1];
  double R0=1.;
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;

  double mm=5.2;
  // 5.88427e+00 
  //double mm2=5.83252E+00;
  //double R1;


  double funcdown;
  if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n2*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n2*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        funcdown=1.;
      }
    }
    

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
    
    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);


      R1=(R1_nume/R1_deno)*(1.-Rdown)+funcdown*Rdown;
      //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        //R1=(R1_nume/R1_deno)*(funcP-0.5);
        R1=(R1_nume/R1_deno)*(1.-Rdown)+funcdown*Rdown;
        //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      }
      else{
        //R1=1.*(1.5-funcP);
        //R1=1*(funcP-0.5);


        R1=1.*(1.-Rdown)+funcdown*Rdown;
        //R1=1.*(1.5-funcP)+funcdown*(funcP-0.5);
        
      
      }
    }

    return R1;

  }
  TF1 *f0 = new TF1("",func_R0,0.01,0.5,3);
  
double func_R00(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double V_Fe=par[2];
  double R1;
  double V_Fe_p=V_Fe+mu_n*in_mag;
  double V_Fe_m=V_Fe-mu_n*in_mag;
  //double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  //double qc=1.36239E-01;
  double qc=1.30030e-01;
  double mm2=5.88427e+00 ;

  double ww=2.5E-03;//par[1];
  double R0=1.;
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;

  double mm=5.2;
  //double mm2=5.83252E+00;
  
  //double R1;

  double funcdown;
  if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n2*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n2*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        funcdown=1.;
      }
    }
    

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

    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);


      //R1=(R1_nume/R1_deno)*(funcP-0.5)+funcdown*(1.5-funcP);
      R1=(R1_nume/R1_deno)*Rdown+funcdown*(1.-Rdown);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        //R1=(R1_nume/R1_deno)*(funcP-0.5);
        //R1=(R1_nume/R1_deno)*(funcP-0.5)+funcdown*(1.5-funcP);
        R1=(R1_nume/R1_deno)*Rdown+funcdown*(1.-Rdown);
      }
      else{
        //R1=1.*(1.5-funcP);
        //R1=1*(funcP-0.5);


        //R1=1.*(funcP-0.5)+funcdown*(1.5-funcP);
        R1=1.*Rdown+funcdown*(1.-Rdown);
        
      
      }
    }
    //R1=Rdown;

    return R1;

  }
  TF1 *f1 = new TF1("",func_R00,0.01,0.5,3);

  //E>V1(-V)
  /*
  double func_R1(double *qqq,double *par){
    double q1=qqq[0];
    double E1=pow(hbar*q1,2)/8./m_nc2nm;
    double d1=par[0];
    double in_mag=par[1];
    double R1;

    double funcP;
    if(q1<qc){
      funcP=R0;  
    }
    
    else{
      if(q1>qc){
        double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        funcP=up_R+down_R;
      }
    }

    if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n2*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      R1=(R1_nume/R1_deno)*(funcP-0.5);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n2*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        R1=(R1_nume/R1_deno)*(funcP-0.5);
        //R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        R1=1.*(funcP-0.5);
      }
    }

    return R1;

  }
  //TF1 *f1 = new TF1("",func_R1,0.01,0.5,2);

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
      //R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      R1=(R1_nume/R1_deno);//  *(funcP-0.5);
      
    }
    else{
      
        double k0=sqrt(2.*m_nc2*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        //double k2=sqrt(2.*m_nc2*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*k0+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k0)*cosh(alpha*d1),2)+pow((k0*k0-alpha*alpha)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        R1=1.;//  *(funcP-0.5);
    }

    return R1;

  }
  TF1 *f2 = new TF1("",func_R2,0.01,0.5,2);
*/

Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 666.; //sample to detector [mm]
//Double_t xdirect    = 63.29;//62.9442
Double_t xdirect    = 60.4483;//62.9442;
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


Int_t M1_pol_check_ichi3(){
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

  TH1F* hq0BG[num];
  TH1F* hqBG[num];
  TH1F* hq02BG[num];
  TH1F* hq2BG[num];

  Int_t kp11[num1];
  Int_t kp211[num1];
  TTree* tup11[num1];
  TTree* tup211[num1];
  //TTree* tup02[num1];

  TH1F* hx11[num1];
  TH1F* hlambda11[num1];
  TH1F* hratio11[num1];
  TH1F* hq11[num1];
  TH1F* hq011[num1];
  TH3F* hxylambda11[num1];

  TH1F* hpolratio11[num1];
  TH1F* hpolratio211[num1];
  TH1F* hpolratio311[num1];
  TH1F* hq211[num1];
  TH1F* hq0211[num1];

  TH1F* hq0BG11[num1];
  TH1F* hqBG11[num1];
  TH1F* hq02BG11[num1];
  TH1F* hq2BG11[num1];

  TTree* tup0;
  Int_t kp0;


  Double_t angle11[num1];
  Double_t angle211[num1];
  for(Int_t i=0; i< num; i++){
    angle11[i]=TMath::Abs(45.0572 - xdirect)/dist_det;
    angle211[i]=TMath::Abs(45.0572 - xdirect)/dist_det;
  }  

  Double_t angledeg11[num];
  Double_t angledeg211[num];
  
  angledeg11[0]=angle11[0]*180./TMath::Pi()/2.;
  angledeg11[1]=angle11[1]*180./TMath::Pi()/2.;
  angledeg11[2]=angle11[2]*180./TMath::Pi()/2.;
  


  TString namestr[num];
  TString namestr2[num];
  //off
  namestr[0]="20210714193654_list.root";
  namestr[1]="20210716232122_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  //namestr[2]="20210717004515_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr[2]="20210716233530_list.root"; //AFP ON
  
  namestr[3]="20210717022140_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  //namestr[4]="20210715085349_list.root"; //0.264A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT
  //namestr[5]="20210715082606_list.root"; //0.378A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT  -> 2 mT
  //namestr[6]="20210715083711_list.root"; //0.6A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 2 mT  -> 3 mT

  //on
  namestr2[0]="20210714193654_list.root";
  namestr2[1]="20210716233530_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  //namestr2[2]="20210717005023_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr2[2]="20210716233530_list.root"; //AFP ON
  namestr2[3]="20210717022920_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  TString degstr[num];
  TString degstr2[num];
  //off
  degstr[0]="Direct(M1 reflect)";
  degstr[1]="B = 8.13 mT, SF OFF";
  degstr[2]="B = 8.13 mT, SF ON";
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
  //angle[1] = TMath::Abs(47.4493 - xdirect)/dist_det; //rad 47.1868(30nm)47.2685
  angle[1] = TMath::Abs(45.1074 - xdirect)/dist_det; 
  
  //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(45.0572 - xdirect)/dist_det; //rad 47.04
  angle[3] = TMath::Abs(47.4013 - xdirect)/dist_det; //rad 47.2
  //angle[4] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[5] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[6] = TMath::Abs(47.2 - xdirect)/dist_det; //rad

  angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle2[1] = TMath::Abs(45.1074 - xdirect)/dist_det; //rad 47.07
  angle2[2] = TMath::Abs(45.0572- xdirect)/dist_det; //rad 47.1 
  angle2[3] = TMath::Abs(47.4845 - xdirect)/dist_det; //rad 47.19
  //angle2[4] = TMath::Abs(47.24 - xdirect)/dist_det; //rad
  //angle2[5] = TMath::Abs(47.21 - xdirect)/dist_det; //rad
  //angle2[6] = TMath::Abs(47.11 - xdirect)/dist_det; //rad

  Double_t angledeg[num1];
  Double_t angledeg2[num1];
  
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  angledeg[3]=angle[3]*180./TMath::Pi()/2.;


  
  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.8, 0.20, 1.0, 0.40,"");
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
  //  Double_t HLD  = 7400.;

  //  Double_t xbegin=54.;
  /*
  Double_t xbegin=19.;
  Double_t xcenter=29.;
  Double_t xend=40.;
  */
  /*
  Double_t xbegin1=19.;
  Double_t xcenter1=29.;
  Double_t xend1=40.;
  */
  
  //Double_t xbegin1=20.;
  //Double_t xcenter1=35.;
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

  for(Int_t i=0; i< num1; i++){
    degstr11[i]=Form("SF:ON,  B=%lf mT", vec_H[i]);
    degstr211[i]=Form("SF:OFF, B=%lf mT", vec_H[i]);
    degstr_p11[i]=Form("B=%lf mT", vec_H[i]);

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

  TString namestr11[num1];
  TString namestr211[num1];
  
  TString namestr_ref= "data/210713_SiFe/20210714193654_list.root"; // file path of the direct data
  //0N
  for(Int_t i=0; i<num1; i++){
    int iscan=i*2;
    namestr11[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  //OFF
  for(Int_t i=0; i<num1; i++){
    int iscan=i*2+1;
    namestr211[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  
  tup0 = GetTree1(namestr_ref); // direct data 
    tup0->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); 
    if(useMRfirst) kp0 = tup0->GetMaximum("mp");
      else kp0= (tup0->GetMaximum("kp") - tup0->GetMinimum("kp"));
  cout << "direct data # of kp: "<< kp0 <<endl;


  // for(Int_t i=0; i<2; i++){

  for(Int_t i=0; i<num; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;
    if(i==0) thecut01=thecut1;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t twopirad2 = 2*TMath::Pi()*angle2[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    Double_t twopirad11 = 2*TMath::Pi()*angle11[i];
    Double_t twopirad211 = 2*TMath::Pi()*angle211[i];

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
    //hq11[i]->SetTitle("Reflectivity (SF ON)");
    //leg->AddEntry(hq1[i],degstr[i],"l");


    // Direct data

    hq0211[i] = new TH1F(Form("hq0211%d",i),Form("%s;q [nm^{-1}];count/bin/25kp","Direct"),nbin_q,q_min,q_max);
    tup0->Draw(Form("%f/(toffo*%f)>>hq0211%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    hq0211[i]->Scale(25./kp0);

     // Data with SF:ON
    tup211[i] = GetTree1(namestr211[i]);
    tup211[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    if(useMRfirst) kp211[i] = tup211[i]->GetMaximum("mp");
      else kp211[i] = (tup211[i]->GetMaximum("kp") - tup211[i]->GetMinimum("kp"));
    cout << "SF:ON data # of kp: "<< kp211[i] <<endl;

    hq211[i] = new TH1F(Form("hq211%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr211[i].Data()),nbin_q,q_min,q_max);
    // tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup211[i]->Draw(Form("%f/(toffo*%f)>>hq211%d",twopirad211,lambda_coeff,i), thecut && cut_ref,"goff");
    hq211[i]->Scale(25./kp211[i]);
    hq211[i]->Divide(hq0211[i]); //Divide by the direct data
    hq211[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq211[i]->SetTitle("Reflectivity (SF ON)");
    //leg->AddEntry(hq1[i],degstr[i],"l");



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
    hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq0BG[i] = new TH1F(Form("hq0BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hqBG[i] = new TH1F(Form("hqBG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    
    // hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);
  
    //tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut1 && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0BG%d",twopirad,lambda_coeff,i), thecut01 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hqBG%d",twopirad,lambda_coeff,i), thecut1 && cut_ref1,"goff");
    
    // tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    
    //hpolratio[i]->Rebin(10);
    //hpolratio2[i]->Rebin(10);

    //hpolratio[i]->Divide(hpolratio2[i]);


    if(i!=0){

      if(i!=3){
        leg->AddEntry(hq[i],degstr[i],"l");
      }
      
    }
    //leg->AddEntry(hq[i],degstr[i],"l");}
    leg2->AddEntry(hq[i],degstr[i],"l");
    if(i!=0)leg3->AddEntry(hq[i],degstr[i],"l");
    //leg->AddEntry(hq2[i],degstr2[i],"l2");

    //hx[i]->Scale(25./kp[i]);

    //9/6add
    hq[i]->Add(hq[i], hqBG[i],1., -1.);
    hq0[i]->Add(hq0[i], hq0BG[i],1., -1.);

    hlambda[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    hqBG[i]->Scale(25./kp[i]);
    hq0BG[i]->Scale(25./kp[0]);

    
    //hxylambda[i]->Scale(25./kp[i]);
    
    
    
    //hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    //hratio[i]->Divide(hlambda[0]);
    //hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    

    hq[i]->Divide(hq0[i]);
    hq[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq[i]->SetTitle("Reflectivity (SF OFF)");
    //hq2[i]->Divide(hq02[i]);

    if(useMRfirst) kp2[i] = tup2[i]->GetMaximum("mp");
    else kp2[i] = tup2[i]->GetMaximum("kp");
    //cout << kp2[i]<<endl;
    hq02[i] = new TH1F(Form("hq02%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    hq2[i] = new TH1F(Form("hq2%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    hq02BG[i] = new TH1F(Form("hq02BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq2BG[i] = new TH1F(Form("hq2BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    
    tup2[0]->Draw(Form("%f/(toffo2*%f)>>hq02%d",twopirad2,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup2[i]->Draw(Form("%f/(toffo2*%f)>>hq2%d",twopirad2,lambda_coeff,i), thecut && cut_ref,"goff");
    tup2[0]->Draw(Form("%f/(toffo*%f)>>hq02BG%d",twopirad2,lambda_coeff,i), thecut01 && cut_dir,"goff");
    tup2[i]->Draw(Form("%f/(toffo*%f)>>hq2BG%d",twopirad2,lambda_coeff,i), thecut1 && cut_ref1,"goff");

    hq2[i]->Add(hq2[i], hq2BG[i],1., -1.);
    hq02[i]->Add(hq02[i], hq02BG[i],1., -1.);
    hq2[i]->Scale(25./kp2[i]);
    hq02[i]->Scale(25./kp2[0]);
    hq2BG[i]->Scale(25./kp2[i]);
    hq02BG[i]->Scale(25./kp2[0]);

    hq2[i]->Divide(hq02[i]);

    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq[i]->SetTitle("");
    hq2[i]->SetTitle("");
    //hq02[i]->GetYaxis()->SetTitle("Reflectivity02");

    hpolratio[i]=(TH1F*)hq[i]->Clone(Form("hpolratio%d",i));
    hpolratio[i]->Add(hq[i], hq2[i],1., -1.); // hq: OFF, hq2: ON, calculate Non - Noff
    // hpolratio[i]->Add(hq2[i],-1.);//各binの値=x*hist1の値+y*hist2の値
    hpolratio2[i]=(TH1F*)hq[i]->Clone(Form("hpolratio2%d",i));
    hpolratio2[i]->Add(hq[i], hq2[i],1., 1.); // hq: OFF, hq2: ON, calculate Non + Noff
    hpolratio[i]->Divide(hpolratio2[i]);
    //hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{on}-R_{off})/(R_{on}+R_{off})");
    hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{off}-R_{on})/(R_{off}+R_{on})");    
    hpolratio[i]->SetTitle("Polarization power");

    if(i==9){
      //hx[i]->SetLineColor(i+2);
      hlambda[i]->SetLineColor(i+2);
      //hratio[i]->SetLineColor(i+2);
      hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hqBG[i]->SetLineColor(i+2);
      hq0BG[i]->SetLineColor(i+2);
      hq02BG[i]->SetLineColor(i+2);
      hq2BG[i]->SetLineColor(i+2);

      hpolratio[i]->SetLineColor(i+2);
      hpolratio2[i]->SetLineColor(i+2);
    } else {
      //hx[i]->SetLineColor(i+1);
      hlambda[i]->SetLineColor(i+1);
      //hratio[i]->SetLineColor(i+1);
      hq[i]->SetLineColor(i+1);
      hq0[i]->SetLineColor(i+1);
      hq2[i]->SetLineColor(i+1);
      hpolratio[i]->SetLineColor(i+1);
      hpolratio2[i]->SetLineColor(i+1);
      hqBG[i]->SetLineColor(i+1);
      hq0BG[i]->SetLineColor(i+1);
      hq02BG[i]->SetLineColor(i+1);
      hq2BG[i]->SetLineColor(i+1);
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
    c1->cd(1);
    leg->Draw();
    if(i!=0){
      /*
    if(i==1){
      hq[i]->Draw("eh");
      //hq11[i]->Draw("ehsame");
    }
    if(i==2){
      hq[i]->Draw("ehsame");
      //hq11[i]->Draw("ehsame");
    }*/
    

    /*if(i==0){
      //hq[i]->Draw("ehsame");
      hq11[i]->Draw("ehsame");
    }*/
    
    //else hq[i]->Draw("ehsames");
   
    // 10/17
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    //fitter.FitFCN(6,globalChi2,0,dataB.Size()+dataSB.Size(),true);
    //TH1D * hB = new TH1D("hB","histo B",100,0,100);
    //TH1D * hSB = new TH1D("hSB","histo S+B",100, 0,100);
  
    //TF1 * fB = new TF1("fB","expo",0,100);
    //fB->SetParameters(1,-0.05);
    //hB->FillRandom("fB");
  
    //TF1 * fS = new TF1("fS","gaus",0,100);
    //fS->SetParameters(1,30,5);
  
    //hSB->FillRandom("fB",2000);
    //hSB->FillRandom("fS",1000);
  
    // perform now global fit
  
    //TF1 * fSB = new TF1("fSB","expo + gaus(2)",0,100);
  
    ROOT::Math::WrappedMultiTF1 wfB(*f0,3);//??
    ROOT::Math::WrappedMultiTF1 wfSB(*f1,3);
  
    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange rangeB;
    // set the data range
    rangeB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataB(opt,rangeB);
    ROOT::Fit::FillData(dataB, hq[1]);
  
    ROOT::Fit::DataRange rangeSB;
    rangeSB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    ROOT::Fit::FillData(dataSB, hq[2]);
  
    ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
    ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);
  
    GlobalChi2 globalChi2(chi2_B, chi2_SB);
  
    ROOT::Fit::Fitter fitter;
  
    const int Npar = 3;
    double par0[Npar] = { 90.e-9,2.,209.0602};
  
    // create before the parameter settings in order to fix or set range on them
    fitter.Config().SetParamsSettings(3,par0);
    // fix 5-th parameter
    //fitter.Config().ParSettings(4).Fix();

    // set limits on the third and 4-th parameter
    //fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
    //fitter.Config().ParSettings(3).SetLimits(0,10000);
    //fitter.Config().ParSettings(3).SetStepSize(5);
  
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
  
    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(3,globalChi2,0,dataB.Size()+dataSB.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
  
    //TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",10,10,700,700);
    
    //c1->Divide(1,2);
    //c1->cd(1);
    //gStyle->SetOptFit(1111);
  
    f0->SetFitResult( result, iparB);
    f0->SetRange(rangeB().first, rangeB().second);
    //f0->SetLineColor(kBlue);
    if(i==1)hq[i]->GetListOfFunctions()->Add(f0);
    if(i==1)hq[i]->Draw();


  
    //c1->cd(2);
    f1->SetFitResult( result, iparSB);
    f1->SetRange(rangeSB().first, rangeSB().second);
    //fSB->SetLineColor(kRed);
    if(i==2)hq[i]->GetListOfFunctions()->Add(f1);
    if(i==2)hq[i]->Draw();
    result.Print(std::cout);

    /*
  
    f0->SetFitResult( result, iparB);
    
    f0->SetParLimits(0,80.e-9,100.e-9);
    
    //f0->SetParameter(0.,30.);//nm 
    //f0->SetParameter(1,2.);
    //f0->SetParLimits(0,.e-9,100.e-9);
    f0->SetParLimits(1,1.8,2.3);
    //f0->SetParLimits(2,190,211);
    //f0->FixParameter(0.,94.37e-9);//nm 


    //f0->FixParameter(1.,2.);//nm 
    f0->FixParameter(2.,209.0602);//nm 

    //f0->SetParLimits(2,190,211);
    //hq[i]->Fit("f0","R","10000",0.173,0.5);
    f0->SetNpx(10000);
    //f0->SetParLimits(1,1.8,1.999);
    //f0->Draw("sames");

    //f1->SetParameter(0.,90.e-9);//m //第１引数が変数の番号、第２引数がその値
    //f1->SetParameter(0.,30.);//nm 
    //f1->SetParameter(1,2.);

    f1->SetFitResult( result, iparSB);
    f1->SetParLimits(0,80.e-9,100.e-9);
    //f1->FixParameter(0.,94.37e-9);
    f1->SetParLimits(1,1.5,2.2);
    //f1->SetParLimits(2,150,220);
    //f1->FixParameter(1.,2.);//nm 
    //
    f1->FixParameter(2.,209.0602);//nm 

    //f0->SetParameter(0.,30.);//nm 
    //f0->SetParameter(1,2.);
    
    //f1->SetParLimits(1,1.8,1.999);
    f1->SetNpx(10000);
    f1->SetLineColor(4);
    //f1->Draw("sames");

    //f2->SetParameter(0.,90.e-9);//m //第１引数が変数の番号、第２引数がその値
    //f1->SetParameter(0.,30.);//nm 
    //f2->SetParameter(1,2.);
    
    //f2->SetNpx(10000);
    //f2->Draw("sames");
    */
    
    f0->SetNpx(10000);
    //f0->SetLineColor();
    f1->SetNpx(10000);
    f1->SetLineColor(4);
    /*
    if(i==1){
      //hq[1]->Fit("f0","","",0.252,0.5);
      hq[i]->Fit(f0,"+","",0.1,0.47);
    }
    if(i==2){
      //hq[2]->Fit("f1","","",0.252,0.5);
      hq[i]->Fit(f1,"+","",0.1,0.47);
    }
    */
    gStyle->SetOptFit(1111);
    
  

  }

    c1->cd(2);
    leg->Draw();
    if(i!=0){
    if(i==1){
      hq[i]->Draw("eh");
      //hq11[i]->Draw("ehsame");
    }
    if(i==2){
      hq[i]->Draw("ehsame");
      //hq11[i]->Draw("ehsame");
    }
    
    if(i==1){
      //hq[1]->Fit("f0","","",0.252,0.5);
      hq[i]->Fit(f0,"+","",0.1,0.5);
    }
    if(i==2){
      //hq[2]->Fit("f1","","",0.252,0.5);
      hq[i]->Fit(f1,"+","",0.1,0.5);
    }
    gStyle->SetOptFit(1111);
    
  

  }

    /*if(i==0)hlambda[i]->Draw("eh");
    else hlambda[i]->Draw("ehsames");
    leg->Draw();*/

    
/*
    TBox* b = new TBox(0.287,0,0.5,2.); 
    b->SetFillColor( 7 ); 
    b->SetFillStyle(3004); 
    b->Draw(); 
    TBox* b1 = new TBox(0.15,0,0.173,2.); 
    b1->SetFillColor( kOrange ); 
    b1->SetFillStyle(3004); 
    b1->Draw("same");
    leg->Draw();
    */

//9/6
/*
    if(i==1)hq2[i]->Draw("eh");
    else hq2[i]->Draw("ehsames");
    // if(i==1)hq2[i]->Draw("ah");
    // else hq2[i]->Draw("ahsames");
  */  

    double ymax=1.1;//gPad->GetUymax();
    double ymin=0.;//gPad->GetUymin();
    double xmax1=0.2;
    double xmin1=0.4;


    c1->cd(3);
    /*if(i==1){
      hq2[i]->Draw("eh");
    }*/
    if(i==1){
      hq2[i]->Draw("ehsame");
      hq211[i]->Draw("ehsame");
    }
    /*if(i==0){
      //hq[i]->Draw("ehsame");
      hq11[i]->Draw("ehsame");
    }*/
    
    else hq2[i]->Draw("ehsames");
    /*
    if(i!=0){
    if(i==1)hpolratio[i]->Draw("eh");
    else hpolratio[i]->Draw("ehsames");
    */
/*
    if(i==0){
      hq0BG[i]->Draw("eh");
      //hq02BG[i]->Draw("ehsame");
    }
    else hqBG[i]->Draw("ehsames");
    //hq2BG[i]->Draw("ehsames");
    leg->Draw();
   */ 

    TBox* b4 = new TBox(0.15,-1,0.173,1.); 
    b4->SetFillColor( kOrange ); 
    b4->SetFillStyle(3004); 
    //b4->Draw();
    TBox* b5 = new TBox(0.287,-1,0.5,1.); 
    b5->SetFillColor(7); 
    b5->SetFillStyle(3004); 
    //b5->Draw("same");
    leg3->Draw();

    //}
    
/*
    c1->cd(3);
    if(i==1)hpolrhqtio[i]->Draw("eh");
    else hpolratio[i]->Draw("ehsames");
    leg->Draw();
*/
    
    hpolratio[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hpolratio[i]->GetYaxis()->SetRangeUser(-1.,1.);
    hpolratio2[i]->GetYaxis()->SetRangeUser(-1.,1.);
    hq[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    hq2[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq2[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    // hq[i]->SaveAs(path_R + Form("hq_off_%d.root", i));
    // hq2[i]->SaveAs(path_R + Form("hq_on_%d.root", i));
    cout<<"in_"<<angledeg[i]<<"_deg"<<endl;
  }

  const Int_t num_pol = 3;
  TGraphErrors* gr[num_pol];
  TGraphErrors* gr11[num_pol];
  TGraphErrors* gr211[num_pol];
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t q_cuts[num_pol] = {0.25, 0.3, 0.225};
  
  Double_t B[num-1] = {8.01,0.322,1.35};
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t pol_at_qcut[num-1];
  Double_t error_pol_at_qcut[num-1];

  Double_t ibin_pol[num_pol];

  Int_t num11=10;

  Double_t q_cuts11[num_pol] = {0.25, 0.3, 0.35};
  
  //Double_t B11[num-1] = {8.01,0.322,1.35};
  Double_t B211[num]={2.0101,2.098,2.1859,2.2738};
  Double_t B11[num]={1.96615,2.05405,2.14195,2.22985};//,2.31775};

  // ={1,2,3,4,5,6,7,8,9};
  //B11[num11-1]={1.96615,2.0101,2.05405,2.098,2.14195,2.1859,2.22985,2.2738,2.31775};
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t pol_at_qcut11[num1];
  Double_t error_pol_at_qcut11[num1];

  Double_t pol_at_qcut211[num1];
  Double_t error_pol_at_qcut211[num1];

  Double_t ibin_pol11[num_pol];

  c1->cd(4); 

  ofstream ofs(path_R+Form("B_RdataOFF_q225_%s.csv", scan_id.c_str()));  // ファイルパスを指定する


  for (Int_t i=0; i<num_pol; i++){
    ibin_pol[i]= Int_t((q_cuts[i]-q_min)*nbin_q/(q_max-q_min));
    Double_t q_cut_i = q_min + (q_max-q_min)*ibin_pol[i]/nbin_q;
    for (Int_t j=0; j<num; j++){
      pol_at_qcut[j] = hq11[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut[j] = hq11[j]->GetBinError(ibin_pol[i]);
      ofs << q_cut_i << "," << B11[j] << ","<< pol_at_qcut[j] << ","<< error_pol_at_qcut[j] << endl;
      cout<<"B11_"<<B11[j]<<endl;

      // cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }

    for (Int_t j=0; j<num; j++){
      pol_at_qcut211[j] = hq211[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut211[j] = hq211[j]->GetBinError(ibin_pol[i]);
      ofs << q_cut_i << "," << B211[j] << ","<< pol_at_qcut211[j] << ","<< error_pol_at_qcut211[j] << endl;
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



  TLegend *leg4 = new TLegend(0.8, 0.8, 0.95, 1, "");
  // leg4->SetFillStyle(0);
  for (Int_t i=0; i<num_pol; i++){
    ibin_pol[i]= Int_t((q_cuts[i]-q_min)*nbin_q/(q_max-q_min));//bin number
    for (Int_t j=0; j<num; j++){
      pol_at_qcut11[j] = hq11[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut11[j] = hq11[j]->GetBinError(ibin_pol[i]);

      cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }
    for (Int_t j=0; j<num; j++){
      pol_at_qcut211[j] = hq211[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut211[j] = hq211[j]->GetBinError(ibin_pol[i]);

      cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }

    for (Int_t j=1; j<num; j++){
      pol_at_qcut[j] = hq[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut[j] = hq[j]->GetBinError(ibin_pol[i]);

      cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }

    gr11[i]= new TGraphErrors(num, B11, pol_at_qcut11,0,error_pol_at_qcut11);
    gr11[i]->GetXaxis()->SetRangeUser(0.2, 9);
    gr11[i]->GetYaxis()->SetRangeUser(0., 1.5);
    if (i==0) gr11[i]->Draw("AP");
    else gr11[i]->Draw("P");
    
    gr211[i]= new TGraphErrors(num, B211, pol_at_qcut211,0,error_pol_at_qcut211);
    gr211[i]->GetXaxis()->SetRangeUser(0.2, 9);
    gr211[i]->GetYaxis()->SetRangeUser(0., 1.5);
    if (i==0) gr211[i]->Draw("APsame");
    else gr211[i]->Draw("Psame");

    gr[i]= new TGraphErrors(num, B, pol_at_qcut,0,error_pol_at_qcut);
    gr[i]->GetXaxis()->SetRangeUser(0.2, 9);
    gr[i]->GetYaxis()->SetRangeUser(0., 1.5);
    if (i==0) gr[i]->Draw("APsame");
    else gr[i]->Draw("Psame");
    
    

    gr[i]->SetMarkerColor(i+1);
    gr[i]->SetLineColor(i+1);
    gr[i]->SetMarkerStyle(i+3);
    gr[i]->SetMarkerSize(1);
    gr[i]->GetXaxis()->SetTitle("B (mT)");
    gr[i]->GetYaxis()->SetTitle("Polarization power");
    gr[i]->SetTitle("");

    gr211[i]->SetMarkerColor(i+1);
    gr211[i]->SetLineColor(i+1);
    gr211[i]->SetMarkerStyle(i+3);
    gr211[i]->SetMarkerSize(1);
    gr211[i]->GetXaxis()->SetTitle("B (mT)");
    gr211[i]->GetYaxis()->SetTitle("Polarization power");
    gr211[i]->SetTitle("");
    leg4->AddEntry(gr[i],Form("q=%.3f nm^{-1}",q_min + (q_max-q_min)*ibin_pol[i]/nbin_q),"p");

  }
  leg4->Draw();

  c1->cd(1); gPad->SetGrid();//gPad->SetLogy();
  c1->cd(2); gPad->SetGrid();gPad->SetLogy();
  
  c1->cd(3); gPad->SetGrid();//gPad->SetLogy();
  // c1->cd(4); gPad->SetGrid();//gPad->SetLogy();

  c1->SaveAs(path_R+"pol_90nm.png");
  c1->SaveAs(path_R+"pol_90nm.root");

#if 1
  TFile *outfile = TFile::Open(path_R+"pol_90nm_BG.root","RECREATE");
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
