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
int iparB[3] = { 0,      // exp amplitude in B histo
                 1,
                 2    // exp common parameter
};
 
// signal + background function
int iparSB[3] = { 0,
                  1, // exp amplitude in S+B histo
                  2 // exp common parameter
                  
};
 
// Create the GlobalCHi2 structure
 
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
//double qc=0.13;//par[2];
//double ww=7.16389E-02;//par[3];
//double R0=1.;//par[4];

//double mm=par[3];
//double alpha=-8.03288E-02;//par[5];

double uprate=0.5;//par[3];
double downrate=0.5;//par[4];

double mm=5.;
double qc_up=0.22;
//double mm2=6.24502;

//double qc=1.88418e-1;
double qc=0.154971;

//double qc=1.58418e-1;
double mm2=6.72226;
//double mm2=16.72226;

double q_coff_samp=2*TMath::Pi()* ((60.4483-45.1074)/666.);
double q_coff=2*TMath::Pi()* ((89.3295-65.4023)/1439.);//((89.3295-60.4483)/1439.);
double q_coff_true=2*TMath::Pi()* ((89.3295-60.4583)/1439.);

double Spl1(double x1) {

  //double x2=q_coff_samp/x1;
  //double x=q_coff/q_coff_true*x2;
  //double lam=4*TMath::Pi()*((60.4483-45.1074)/666.)/x;

  //double x2=q_coff_samp/x1;
  //double x=q_coff/q_coff_true*x2;
  
  //double q1=qqq[0];
  //double lam1=qqq[0];
  Double_t dist_det   = 1439.; 
  Double_t xdirect    = 89.3295;
  //double q1=4*TMath::Pi()*(xdirect1-45.1074)/(2*dist_det1*lam1);
  
  double lam=q_coff_samp/x1;
  double q22=q_coff_true/lam;
  //double q11=0.8*q_coff/q_coff_true*q22;
  double q1=q22;

  //4*TMath::Pi()*(xdirect-65.4023)/(2*dist_det*lam);
  //double q1=4*TMath::Pi()*(xdirect-60.4583)/(2*dist_det*lam1);
  //60.4583

  //(xdirect-65.4023)/(2*dist_det*lam1);//upの角度を用いた
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  
  //double lmc=par[0];//critical value (lambda )
  
  //double qc=4*TMath::Pi()*(xdirect-65.4023)/(2*dist_det*lmc);//par[0];
  double ww=3.E-03;//par[1];
  
  double qc=1.35980e-01;//par[0];
  double R0=9.9e-01;//par[1];
  double mm2= 6.56329e+00;//par[2];
  double alpha=0.2;//par[3];
  
  //double mm=par[3];
  //double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;
  double q13=0.13*0.8284697312;
  double q145=0.145*0.8284697312;

  double mm=5.;
  
  double R1;

  double Rup;
    if(q1<qc_up){
      if(q1<q13){
       Rup=uprate*R0;
      }
      if(q1>=q13){
       if(q1<q145){
        Rup=0.99*uprate;
       }
       else if(q1>=q145){
        //
        Rup=uprate*R0;
        //Rdown=0.99*downrate;
       }
      }
       //Rup=uprate*R0;  
    }
    
    else{
      if(q1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }

    double Rdown;
    if(q1<q13){
       Rdown=downrate*R0;
    }
    if(q1>=q13){
       if(q1<q145){
        Rdown=0.99*downrate;
       }
       else if(q1<qc){
        //
        Rdown=downrate*R0;
        //Rdown=0.99*downrate;
        
       }
       else{
        if(q1>=qc){
          //double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
          Rdown=downrate*R0/pow((1.+mm2*(q1-qc)),4);
          //Rdown=downrate*pow((q1-mm2*pow((q1*q1-qc*qc),0.5))/(q1+mm2*pow((q1*q1-qc*qc),0.5)),2);
          //double down_R=downrate*pow((q1-R0*pow((q1*q1-qc*qc),0.5))/(q1+R0*pow((q1*q1-qc*qc),0.5)),2);
          //Rdown=down_R;//+down_R;
        }
      }
    }
    R1=Rup+Rdown;
    //R1=Rdown;//+Rdown;
    //R1=Rdown;
    return R1;

    
}



double func_R0(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double V_Fe=par[2];
  double V_Fe_p=V_Fe+mu_n*in_mag;
  double V_Fe_m=V_Fe-mu_n*in_mag;
  double R1,R11;
  //double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  //double qc=1.36239E-01;
  //double qc=1.30030e-01*sin(0.0231098)/sin(0.0159505);
  //double qc=1.24590e-01*sin(0.0231098/2.)/sin(0.0159505/2.);
  //double mm2=5.88427e+00 ;
  //double mm2=6.25737;
  //double qc=1.70382e-01;//
  
  

  double ww=3.E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.2;//par[3];
  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];

  //model parameter
  double qc_up=0.217;
  double mm=5.;


  // 5.88427e+00 
  //double mm2=5.83252E+00;
  //double R1;


  double funcdown;
  if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe-mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe-mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

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

    /*
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
    */
    double Rdown;
    if(q1>0.13){
      if(q1<0.145){
        Rdown=(R0-Rup)/R0;//1.-Rup;
      }
      else{
        //RR=Spl1(x1);
        Rdown=(Spl1(q1)-Rup)/R0;
      }
    }
    else{
      //RR=Spl1(x1);
      Rdown=(Spl1(q1)-Rup)/R0;
    }
    //Rdown=Spl1(q1)-Rup;
    
    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);

      R11=(R1_nume/R1_deno);//*(1.-Rdown)+funcdown*Rdown;
      //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        //R1=(R1_nume/R1_deno)*(funcP-0.5);
        R11=(R1_nume/R1_deno);//*(1.-Rdown)+funcdown*Rdown;
        //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      }
      else{
        //R1=1.*(1.5-funcP);
        //R1=1*(funcP-0.5);
        R11=1.;//*(1.-Rdown)+funcdown*Rdown;
        //R1=1.*(1.5-funcP)+funcdown*(funcP-0.5);   
      }
    }

    R1=R11*(1.-Rdown)+funcdown*Rdown;

    return R1;

  }
  TF1 *f0 = new TF1("",func_R0,0.01,0.5,3);

double Spl2(double *x, double *par){
  double x1=x[0];
  double RR;
  //double R0=0.99;
  double ww=3.E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.2;//par[3];
  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];

  //model parameter
  double qc_up=0.217;
  double mm=5.;

  double Rup;
    if(x1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(x1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((x1-mm*qc_up)/ww))*(1.-alpha*(x1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }
  
  if(x1>0.13){
      if(x1<0.145){
        RR=(R0-Rup-Rup)/R0;
        //RR=(R0-Rup)/R0;;
      }
      else{
        //RR=(Spl1(x1)-Rup)/R0;
        RR=(Spl1(x1)-Rup-Rup)/R0;
      }
  }
  else{
    //RR=(Spl1(x1)-Rup)/R0;
    RR=(Spl1(x1)-Rup-Rup)/R0;
  }
  //RR=Spl1(x1);
  return RR;
}
double Spl3(double *x, double *par){
  double x1=x[0];
  double RR;
  //double R0=0.99;
  double ww=3.E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.2;//par[3];
  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];

  //model parameter
  double qc_up=0.217;
  double mm=5.;

  double Rup;
    if(x1<qc_up){
       Rup=uprate*R0;
    }
    
    else{
      if(x1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((x1-mm*qc_up)/ww))*(1.-alpha*(x1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }
  
  if(x1>0.13){
      if(x1<0.145){
        //RR=(R0-Rup-Rup)/R0;
        RR=(Rup-R0+Rup)/R0;
        //RR=(Rup-Spl1(x1)-Rup)/R0;
        //RR=(R0-Rup)/R0;
      }
      else{
        //RR=(Spl1(x1)-Rup-Rup)/R0;
        RR=(Rup-(Spl1(x1)-Rup))/R0;
        //RR=(Spl1(x1)-Rup)/R0;
      }
  }
  else{
    //RR=(Spl1(x1)-Rup)/R0;
    RR=(Rup-(Spl1(x1)-Rup))/R0;
  }
  //RR=Spl1(x1);
  //double RR1=1.-RR;
  double RR1=RR;
  return RR1;
}

TF1 *sp2 = new TF1("",Spl2,0.1,0.5,0);
TF1 *sp3 = new TF1("",Spl3,0.1,0.5,0);
  
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
  //double qc=1.30030e-01;
  
  //double qc=1.24590e-01*sin(0.0231098/2.)/sin(0.0159505/2.);
  //double mm2=5.88427e+00 ;
  //double mm2=6.24502;
  

  double ww=2.5E-03;//par[1];
  double R0=0.99;
  
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
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe-mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe-mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

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

/*
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
    */
   double Rdown;
    if(q1>0.13){
      if(q1<0.145){
        Rdown=(R0-Rup)/R0;///1.-Rup;
      }
      else{
        //RR=Spl1(x1);
        Rdown=(Spl1(q1)-Rup)/R0;
      }
    }
    else{
      //RR=Spl1(x1);
      Rdown=(Spl1(q1)-Rup)/R0;
    }

    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);


      //R1=(R1_nume/R1_deno)*(funcP-0.5)+funcdown*(1.5-funcP);
      R1=(R1_nume/R1_deno)*Rdown+funcdown*(1.-Rdown);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

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
  //if(useThinout==1)tup->SetMaxEntryLoop(10000);
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


Int_t global_proc_nospline(){

  
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
  TH1F* hlm0[num];
  TH1F* hlm[num];
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

  for(int i=0; i<3; i++){
    cout<<angle[i]<<"_[rad]_"<<angledeg[i]<<"_[deg]_"<<endl;
  }

  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.6, 0.60, 1.0, 0.80,"");
  TLegend* leg2 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm ON");
  TLegend* leg3 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm");

  TLegend* leg4 = new TLegend(0.6, 0.40, 1.0, 0.60,"");


  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Int_t nbin_q  = 60;//300 60
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
    
    hlm0[i] = new TH1F(Form("hlm0%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hlm[i] = new TH1F(Form("hlm%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    
    
    
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq0BG[i] = new TH1F(Form("hq0BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hqBG[i] = new TH1F(Form("hqBG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    
    // hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);
  
    //tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0BG%d",twopirad,lambda_coeff,i), thecut01 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hqBG%d",twopirad,lambda_coeff,i), thecut1 && cut_ref1,"goff");
    
    tup[0]->Draw(Form("toffo*%f>>hlm0%d",lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("toffo*%f>>hlm%d",lambda_coeff,i), thecut && cut_ref,"goff");
    


    // tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    
    //hpolratio[i]->Rebin(10);
    //hpolratio2[i]->Rebin(10);

    //hpolratio[i]->Divide(hpolratio2[i]);


    if(i!=0){

      if(i!=3){
        leg->AddEntry(hq[i],degstr[i],"l");
      }
      
    }

   
    if(i==1) leg4->AddEntry(hratio[i],"upspin","l");
    if(i==2) leg4->AddEntry(hratio[i],"downspin","l");


    //leg->AddEntry(hq[i],degstr[i],"l");}
    leg2->AddEntry(hq[i],degstr[i],"l");
    if(i!=0)leg3->AddEntry(hq[i],degstr[i],"l");
    //leg->AddEntry(hq2[i],degstr2[i],"l2");

    //hx[i]->Scale(25./kp[i]);

    //9/6add
    hq[i]->Add(hq[i], hqBG[i],1., -1.);
    hq0[i]->Add(hq0[i], hq0BG[i],1., -1.);

    hlambda[i]->Scale(25./kp[i]);
    hlm0[i]->Scale(25./kp[i]);
    hlm[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    hqBG[i]->Scale(25./kp[i]);
    hq0BG[i]->Scale(25./kp[0]);

    for(int i11=0; i11<nbin; i11++){
      double ff[nbin];
      double ff1[nbin];
      if(i==1){
        ff[i11]=hlm[i]->GetBinContent(i11);
        ff1[i11]=hlm0[i]->GetBinContent(i11);
      }
      //cout<<"hlm_"<<ff[i11]<<"_hlm0_"<<ff1[i11]<<endl;
    }

    
    


    
    //hxylambda[i]->Scale(25./kp[i]);
    
    
    
    hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    hratio[i]->Divide(hlambda[0]);
    hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    
    hlm[i]->Divide(hlm0[i]);
    hlm[i]->GetYaxis()->SetTitle("Reflectivity");
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
      hlm0[i]->SetLineColor(i+2);
      hlm[i]->SetLineColor(i+2);
      hratio[i]->SetLineColor(i+2);
      //hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hqBG[i]->SetLineColor(i+2);
      hq0BG[i]->SetLineColor(i+2);
      hq02BG[i]->SetLineColor(i+2);
      hq2BG[i]->SetLineColor(i+2);

      hpolratio[i]->SetLineColor(i+2);
      hpolratio2[i]->SetLineColor(i+2);
    } 
    else {
      //hx[i]->SetLineColor(i+1);
      hlambda[i]->SetLineColor(i+1);
      hlm0[i]->SetLineColor(i+1);
      hlm[i]->SetLineColor(i+1);
      hratio[i]->SetLineColor(i+1);
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
    TH1F* hq11;
    if(i==1)hq11=(TH1F*)hq[i]->Clone(Form("hq%d",i));
    TH1F* hq22;
    if(i==2)hq22=(TH1F*)hq[i]->Clone(Form("hq%d",i));

    ROOT::Math::WrappedMultiTF1 wfB(*f0,1);
    ROOT::Math::WrappedMultiTF1 wfSB(*f1,1);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange rangeB;
    // set the data range
    //rangeB.SetRange(10,90);
    rangeB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataB(opt,rangeB);
    if(i==2)ROOT::Fit::FillData(dataB, hq[i-1]);

    ROOT::Fit::DataRange rangeSB;
    //rangeSB.SetRange(10,50);
    rangeSB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    if(i==2)ROOT::Fit::FillData(dataSB, hq[i]);

    ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
    ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);

    GlobalChi2 globalChi2(chi2_B, chi2_SB);

    ROOT::Fit::Fitter fitter;

    const int Npar = 3;
    double par0[Npar] = {89.066,2., 209.0602};
  
    // create before the parameter settings in order to fix or set range on them
    fitter.Config().SetParamsSettings(3,par0);
    // fix 5-th parameter
    //fitter.Config().ParSettings(2).Fix();
    // set limits on the third and 4-th parameter
    //fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
    //fitter.Config().ParSettings(3).SetLimits(0,10000);

    fitter.Config().ParSettings(0).Fix();
    //fitter.Config().ParSettings(0).SetLimits(88.,100.);
    fitter.Config().ParSettings(1).SetLimits(1.8,2.2);
    fitter.Config().ParSettings(2).SetLimits(150.,220.);
   
    
    //fitter.Config().ParSettings(3).SetStepSize(5);
  
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
  
    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(3,globalChi2,0,dataB.Size()+dataSB.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

    c1->cd(1);
    gStyle->SetOptFit(1111);

    f0->SetNpx(10000);
    f0->SetFitResult( result, iparB);
    f0->SetRange(rangeB().first, rangeB().second);
    f0->SetLineColor(kRed);
    if(i==2)hq[i-1]->GetListOfFunctions()->Add(f0);
    //if(i==1)f0->Draw("same");
    if(i==2)hq[i-1]->Draw("eh");

    f1->SetNpx(10000);
    f1->SetFitResult( result, iparSB);
    f1->SetRange(rangeSB().first, rangeSB().second);
    f1->SetLineColor(kBlue);
    if(i==2)hq[i]->SetLineColor(4);
    if(i==2)hq[i]->GetListOfFunctions()->Add(f1);
    if(i==2)hq[i]->Draw("ehsame");

    
    
    
    //leg->Draw();

    //if(i!=0){
    
    
      if(i==1){
        hq[i]->SetLineColor(2);
        hq[i]->Draw("eh");

        //hq11[i]->Draw("ehsame");
      }
      if(i==2){
        hq[i]->SetLineColor(4);
        hq[i]->Draw("ehsame");
        //hq11[i]->Draw("ehsame");
      }
      
      
      //else hq[i]->Draw("ehsames");
      
      //f0->SetParLimits(0,80.e-9,100.e-9);
      
      //f0->SetParameter(0.,30.);//nm 
      //f0->SetParameter(1,2.);
      //f0->SetParLimits(0,.e-9,100.e-9);
      
      //f0->SetParLimits(1,1.8,2.3);
      //f1->SetParLimits(1,1.5,2.2);
      f0->FixParameter(1.,2.);//nm 
      f1->FixParameter(1.,2.);//nm 
      

      //f0->SetParLimits(2,190,211);
      f0->FixParameter(0.,94.37);//nm 


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

      f1->SetParLimits(0,80.,100.);
      //f1->FixParameter(0.,94.37e-9);
      
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

      double qq;
      double qq1;
      double EE1;
      qq=sqrt(8*m_nc2nm*(209.062+60.3*1.98))/hbar;
      qq1=sqrt(8*m_nc2nm*(60.3*0.007))/hbar;
      EE1=pow(hbar*0.25,2)/(8*m_nc2nm);
      cout<<"qqqqq"<<qq<<endl;
      cout<<"qqqqq1"<<qq1<<endl;
      cout<<"EE1"<<EE1<<endl;


      if(i==1){
        //hq[1]->Fit("f0","","",0.252,0.5);
        //hq[i]->Fit(f0,"+","",0.1,0.47);
      }
      if(i==2){
        //hq[2]->Fit("f1","","",0.252,0.5);
        //hq[i]->Fit(f1,"+","",0.1,0.47);
      }
      //gStyle->SetOptFit(0101);
      
    

    //}
    
    
    c1->cd(2);
    //gStyle->SetOptFit(1111);
    //gPad->SetLogy();
    //leg->Draw();
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
      //hq[i]->Fit(f0,"","",0.1,0.47);
    }
    if(i==2){
      //hq[2]->Fit("f1","","",0.252,0.5);
      //hq[i]->Fit(f1,"+","",0.1,0.47);
    }
    //gStyle->SetOptFit(0101);
    
    hq[i]->SetStats(0); //非表示
  

  }
    

    

   
    


    double ymax=1.1;//gPad->GetUymax();
    double ymin=0.;//gPad->GetUymin();
    double xmax1=0.2;
    double xmin1=0.4;


    c1->cd(3);
    hratio[i]->SetStats(0); //非表示
    if(i!=0){
      if(i==1){
        //hlm[i]->Draw("eh");
        
        hratio[i]->Draw("eh");
        hratio[i]->SetLineColor(4);
        //hq11[i]->Draw("ehsame");
      }
      if(i==2){
        hratio[i]->SetLineColor(2);
        hratio[i]->Draw("ehsame");
        //hlm[i]->Draw("ehsame");
        //hq11[i]->Draw("ehsame");
      }
    }
  
    leg->Draw();

    
    c1->cd(4);
    
    leg4->Draw();
    c1->DrawFrame(0.1,0.,0.5,1.1);
    sp2->SetNpx(10000);
    sp2->SetLineColor(4);
    sp2->Draw("");
    sp3->SetNpx(10000);
    //sp3->SetLineColor(2);
    sp3->Draw("same");

    //}
    
    hratio[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hratio[i]->GetYaxis()->SetRangeUser(0.,1.);
    hpolratio2[i]->GetYaxis()->SetRangeUser(-1.,1.);
    hq[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    hq2[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq2[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    // hq[i]->SaveAs(path_R + Form("hq_off_%d.root", i));
    // hq2[i]->SaveAs(path_R + Form("hq_on_%d.root", i));
    //cout<<"in_"<<angledeg[i]<<"_deg"<<endl;
  }

  
#if 1
  TFile *outfile = TFile::Open(path_R+"pol_check_ichi2_90nm.root","RECREATE");
  //TFile *outfile = TFile::Open(path_R+"pol_check_ichi2_90nm.csv","RECREATE");
  for(Int_t i=0; i<num; i++){
    //hx[i]->Write();
    //hlambda[i]->Write();
    //hratio[i]->Write();
    if(i==1)hq[i]->Write();
    if(i==2)hq[i]->Write();
    //hxylambda[i]->Write();
    //hq2[i]->Write();
  }
  outfile->Close();
#endif
for(int i=0; i<3;i++){
    ofstream ofs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));  // ファイルパスを指定する
    
    for(int i11=0; i11<nbin_q; i11++){
      double xc[nbin_q];
      double xc1[nbin_q];
      double ff[nbin_q];
      double ffE[nbin_q];
      double ff1[nbin_q];
      double ff1E[nbin_q];
      if(i==1){
        ff[i11]=hq[i]->GetBinContent(i11);
        ffE[i11]=hq[i]->GetBinError(i11);
        xc[i11]=hq[i]->GetXaxis()->GetBinCenter(i11); // 
      }
      if(i==2){
        ff1[i11]=hq[i]->GetBinContent(i11);
        ff1E[i11]=hq[i]->GetBinError(i11);
        xc1[i11]=hq[i]->GetXaxis()->GetBinCenter(i11); // 
      }
      //cout<<"hlm_"<<ff[i11]<<"_hlm0_"<<ff1[i11]<<endl;
      //ofs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
      
      ofs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
    }
}


  return 0 ;
}

