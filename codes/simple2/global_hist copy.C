#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
 
 
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



void global_hist() {

   TString path_R = "results/";
   int num=2;
   TH1F* hq[num];
   Int_t nbin_q  = 60;//300 60
   Double_t q_min  = 0.1;//0.6
   Double_t q_max  = 0.50;//0.6
   const string scan_id = "90nm_scan_fine_3";
   ifstream ifs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));
   double xc[nbin_q];
   double xc1[nbin_q];
   double ff[nbin_q];
   double ffE[nbin_q];
   double ff1[nbin_q];
   double ff1E[nbin_q];

   //ofstream ofs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));  
   
   /*
   for(int i=1; i<3;i++){
      hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",i),nbin_q,q_min,q_max);

      for(int i11=0; i11<nbin_q; i11++){
      double xc[nbin_q];
      double xc1[nbin_q];
      double ff[nbin_q];
      double ffE[nbin_q];
      double ff1[nbin_q];
      double ff1E[nbin_q];
      
      //cout<<"hlm_"<<ff[i11]<<"_hlm0_"<<ff1[i11]<<endl;
      //ofs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
      
      //ifs << xc[i11]<<ff[i11]<<ffE[i11]<< ff1[i11]<< ff1E[i11] << endl;
      ifs << xc<<ff<<ffE<< ff1<< ff1E << endl;
      //ifs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
      }
   }
   */
   

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
 
   c1->cd(2);
   fSB->SetFitResult( result, iparSB);
   fSB->SetRange(rangeSB().first, rangeSB().second);
   fSB->SetLineColor(kRed);
   hSB->GetListOfFunctions()->Add(fSB);
   hSB->Draw();
 
 
}