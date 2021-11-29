#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 666.; //sample to detector [mm]
//Double_t xdirect    = 63.29;//62.9442
Double_t xdirect    = 60.4483;//62.9442;
Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.
   

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

 
 
// definition of shared parameter
// background function
int iparB[2] = { 0,      // exp amplitude in B histo
                 2    // exp common parameter
};
 
// signal + background function
int iparSB[5] = { 1, // exp amplitude in S+B histo
                  2, // exp common parameter
                  3, // gaussian amplitude
                  4, // gaussian mean
                  5  // gaussian sigma
};
 
// Create the GlobalCHi2 structure
 
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}
 
   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[iparB[i] ];
 
      double p2[5];
      for (int i = 0; i < 5; ++i) p2[i] = par[iparSB[i] ];
 
      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }
 
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};


void global_hist_add() {
   int num=3;
   TTree* tup[num];
   TTree* tup2[num];
   
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

   TTree* tup0;
   Int_t kp0;

   TString namestr[num];
   TString namestr2[num];
   //off
   namestr[0]="20210714193654_list.root";
   namestr[1]="20210716232122_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
   //namestr[2]="20210717004515_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
   namestr[2]="20210716233530_list.root"; //AFP ON

   //namestr[3]="20210717022140_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 

   //on
   namestr2[0]="20210714193654_list.root";
   namestr2[1]="20210716233530_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
   //namestr2[2]="20210717005023_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
   namestr2[2]="20210716233530_list.root"; //AFP ON
   //namestr2[3]="20210717022920_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 

   TString degstr[num];
   TString degstr2[num];
   //off
   degstr[0]="Direct(M1 reflect)";
   degstr[1]="B = 8.13 mT, SF OFF";
   degstr[2]="B = 8.13 mT, SF ON";
   //degstr[3]="B = 0.908 mT";
   //degstr[3]="B = 1.35 mT";
   //degstr[5]="B = 1.80 mT";
   //degstr[6]="B = 2.66 mT";

   //on
   degstr2[0]="Direct(M1 reflect)";
   degstr2[1]="B = 8.13 mT";
   degstr2[2]="B = 0.322 mT";
   //degstr2[3]="B = 0.908 mT";
   ////degstr2[3]="B = 1.35 mT";
   //degstr2[5]="B = 1.80 mT";
   //degstr2[6]="B = 2.66 mT";


   Double_t angle[num];
   Double_t angle2[num];
   angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
   //angle[1] = TMath::Abs(47.4493 - xdirect)/dist_det; //rad 47.1868(30nm)47.2685
   angle[1] = TMath::Abs(45.1074 - xdirect)/dist_det; 

   //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
   angle[2] = TMath::Abs(45.0572 - xdirect)/dist_det; //rad 47.04
   //angle[3] = TMath::Abs(47.4013 - xdirect)/dist_det; //rad 47.2

   angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
   angle2[1] = TMath::Abs(45.1074 - xdirect)/dist_det; //rad 47.07
   angle2[2] = TMath::Abs(45.0572- xdirect)/dist_det; //rad 47.1 
   //angle2[3] = TMath::Abs(47.4845 - xdirect)/dist_det; //rad 47.19

   Double_t angledeg[num];
   Double_t angledeg2[num];

   angledeg[0]=angle[0]*180./TMath::Pi()/2.;
   angledeg[1]=angle[1]*180./TMath::Pi()/2.;
   angledeg[2]=angle[2]*180./TMath::Pi()/2.;
   //angledeg[3]=angle[3]*180./TMath::Pi()/2.;




   TString path_R = "results/";
   
   TH1F* hq[num];
   TH1F* hq0[num];
   Int_t kp[num];
   //Int_t kp0;

   //Int_t nbin_q  = 60;//300 60
   //Double_t q_min  = 0.1;//0.6
   //Double_t q_max  = 0.50;//0.6
   //const string scan_id = "90nm_scan_fine_3";
   ifstream ifs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));
   double xc[nbin_q];
   double xc1[nbin_q];
   double ff[nbin_q];
   double ffE[nbin_q];
   double ff1[nbin_q];
   double ff1E[nbin_q];
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
   TString namestr_ref= "data/210713_SiFe/20210714193654_list.root"; // file path of the direct data


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

   for(Int_t i=1; i<num; i++){

      tup0 = GetTree1(namestr_ref); // direct data 
      tup0->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); 
      if(useMRfirst) kp0 = tup0->GetMaximum("mp");
      else kp0= (tup0->GetMaximum("kp") - tup0->GetMinimum("kp"));
      cout << "direct data # of kp: "<< kp0 <<endl;

      hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
      hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
      tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
      tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
      hq[i]->Scale(25./kp[i]);
      hq0[i]->Scale(25./kp[0]);
      hq[i]->Divide(hq0[i]);
      hq[i]->GetYaxis()->SetTitle("Reflectivity");



      //ofstream ofs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));  

      TGraphErrors *grD2O = new TGraphErrors("trend.dat","%lg %*lg %*lg %lg %lg");
   
      TH1D * hB = new TH1D("hB","histo B",100,0,100);
      TH1D * hSB = new TH1D("hSB","histo S+B",100, 0,100);

      
      // // Using TFile
      // TTree* T = new TTree("T", "T");
      // T->ReadFile("data_scans/magnetic_20210716210153.csv", "index/D:I/D:B/D");
      // Int_t num=T->GetEntries();
      // cout << num << endl;

      // Using vector
      vector<Int_t> vec_index;
      vector<Double_t> vec_I, vec_H; 
      vec_index.clear(); vec_I.clear(); vec_I.clear();
      //ifstream fcsv("data_scans/magnetic_20210716210153.csv"); 
      int i_csv = 0;
      Int_t index;     
      Double_t current, magfield;
      
      if (!ifs.is_open())
      {
         exit(EXIT_FAILURE);
      }
      
      string str;
      getline(ifs, str);
      cout << str << endl; // print out the first row
      while (getline(ifs, str))
      {
      // while(!incsv.eof()) {
      // if (i_csv>0){
         
         // cout << str << endl;
      //     double r = 0., z = 0., v =0.;
         sscanf(str.c_str(), "%d,%lf,%lf", &index, &current, &magfield);
         // str >> index >> ",'" >> current >>",'">> magfield;
         vec_index.push_back(index); 
         vec_I.push_back(current);
         vec_H.push_back(magfield);
         // cout << " " << index
         //  << " " << current
         //   << " " << magfield << endl;
         // }
      // i_csv++;
         }
      cout << "vec_"<<vec_index[0] << endl;
      cout << "vec_"<<vec_I[0] << endl;

   
      TF1 * fB = new TF1("fB","expo",0,100);
      fB->SetParameters(1,-0.05);
      hB->FillRandom("fB");
   
      TF1 * fS = new TF1("fS","gaus",0,100);
      fS->SetParameters(1,30,5);
   
      hSB->FillRandom("fB",2000);
      hSB->FillRandom("fS",1000);
   
      // perform now global fit
   
      TF1 * fSB = new TF1("fSB","expo + gaus(2)",0,100);
   
      ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
      ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
   
      ROOT::Fit::DataOptions opt;
      ROOT::Fit::DataRange rangeB;
      // set the data range
      rangeB.SetRange(10,90);
      ROOT::Fit::BinData dataB(opt,rangeB);
      ROOT::Fit::FillData(dataB, hB);
   
      ROOT::Fit::DataRange rangeSB;
      rangeSB.SetRange(10,50);
      ROOT::Fit::BinData dataSB(opt,rangeSB);
      ROOT::Fit::FillData(dataSB, hSB);
   
      ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
      ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);
   
      GlobalChi2 globalChi2(chi2_B, chi2_SB);
   
      ROOT::Fit::Fitter fitter;
   
      const int Npar = 6;
      double par0[Npar] = { 5,5,-0.1,100, 30,10};
   
      // create before the parameter settings in order to fix or set range on them
      fitter.Config().SetParamsSettings(6,par0);
      // fix 5-th parameter
      fitter.Config().ParSettings(4).Fix();
      // set limits on the third and 4-th parameter
      fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
      fitter.Config().ParSettings(3).SetLimits(0,10000);
      fitter.Config().ParSettings(3).SetStepSize(5);
   
      fitter.Config().MinimizerOptions().SetPrintLevel(0);
      fitter.Config().SetMinimizer("Minuit2","Migrad");
   
      // fit FCN function directly
      // (specify optionally data size and flag to indicate that is a chi2 fit)
      fitter.FitFCN(6,globalChi2,0,dataB.Size()+dataSB.Size(),true);
      ROOT::Fit::FitResult result = fitter.Result();
      result.Print(std::cout);
   
      TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
                                 10,10,700,700);
      c1->Divide(1,2);
      c1->cd(1);
      gStyle->SetOptFit(1111);
   
      fB->SetFitResult( result, iparB);
      fB->SetRange(rangeB().first, rangeB().second);
      fB->SetLineColor(kBlue);
      hB->GetListOfFunctions()->Add(fB);
      hB->Draw();

      if(i!=0){
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
   }
 
   c1->cd(2);
   fSB->SetFitResult( result, iparSB);
   fSB->SetRange(rangeSB().first, rangeSB().second);
   fSB->SetLineColor(kRed);
   hSB->GetListOfFunctions()->Add(fSB);
   hSB->Draw();
 
 
}