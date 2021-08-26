#include "../bin/MakeNiki.h"
#include <iostream> 
#include <fstream>
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TAxis.h"

using namespace std;

using namespace std;

Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 666.; //sample to detector [mm]
Double_t xdirect    = 63.29;
Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.

TString path_R = "results/";
const string scan_id = "30nm_scan_1";
const string run_id="20210716210153";
const Int_t num = 5; // this should be the half of the number of the files obtained by the scan 


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

TTree* GetTree(TString ROOTstr_path){

  TFile *file = TFile::Open(ROOTstr_path.Data());
  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}


Int_t Pol_Power_30nm_scan_1(){

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

  Int_t kp[num];
  Int_t kp2[num];
  Int_t kp0;
  TTree* tup[num];
  TTree* tup2[num];
  TTree* tup0;
  TH1F* hx[num];
  TH1F* hlambda[num];
  TH1F* hratio[num];
  TH1F* hq1[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TH1F* hpolratio[num];
  TH1F* hpolratio2[num];
  TH1F* hpolratio3[num];
  TH1F* hq2[num];
  TH1F* hq02[num];

  TString namestr[num];
  TString namestr2[num];
  
  TString namestr_ref= "data/210713_SiFe/20210714193654_list.root"; // file path of the direct data
  //0N
  for(Int_t i=0; i<num; i++){
    int iscan=i*2;
    namestr[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  //OFF
  for(Int_t i=0; i<num; i++){
    int iscan=i*2+1;
    namestr2[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }


  TString degstr[num];
  TString degstr2[num];
  TString degstr_p[num];

  for(Int_t i=0; i< num; i++){
    degstr[i]=Form("SF:ON,  B=%lf mT", vec_H[i]);
    degstr2[i]=Form("SF:OFF, B=%lf mT", vec_H[i]);
    degstr_p[i]=Form("B=%lf mT", vec_H[i]);

  }  
  
  Double_t angle[num];
  Double_t angle2[num];
  for(Int_t i=0; i< num; i++){
    angle[i]=TMath::Abs(47.1868 - xdirect)/dist_det;
    angle2[i]=TMath::Abs(47.1868 - xdirect)/dist_det;
  }  
  
  Double_t angledeg[num];
  Double_t angledeg2[num];
  
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  
  TLegend* leg = new TLegend(0.8, 0.5, 1.0, 1,"");
  TLegend* leg2 = new TLegend(0.8, 0.5, 1.0, 1,"");
  TLegend* leg3 = new TLegend(0.8, 0.5, 1.0, 1,"");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 60;//300
  Double_t q_min  = 0.15;//0.6 
  Double_t q_max  = 0.50;//0.6
  //Double_t q_max  = 1.0;//0.6


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
  
 

  tup0 = GetTree(namestr_ref); // direct data 
    tup0->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); 
    if(useMRfirst) kp0 = tup0->GetMaximum("mp");
      else kp0= (tup0->GetMaximum("kp") - tup0->GetMinimum("kp"));
  cout << "direct data # of kp: "<< kp0 <<endl;



  for(Int_t i=0; i<num; i++){
    thecut.Print();
    // if(i==0) thecut0=thecut;
    thecut0=thecut;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t twopirad2 = 2*TMath::Pi()*angle2[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

  // Direct data

    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp","Direct"),nbin_q,q_min,q_max);
    tup0->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    hq0[i]->Scale(25./kp0);
    
    // Data with SF:ON
    tup[i] = GetTree(namestr[i]);
    tup[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
      else kp[i] = (tup[i]->GetMaximum("kp") - tup[i]->GetMinimum("kp"));
    cout << "SF:ON data # of kp: "<< kp[i] <<endl;

    hq1[i] = new TH1F(Form("hq1%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    // tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq1%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    hq1[i]->Scale(25./kp[i]);
    hq1[i]->Divide(hq0[i]); //Divide by the direct data
    hq1[i]->GetYaxis()->SetTitle("Reflectivity");
    hq1[i]->SetTitle("Reflectivity (SF ON)");
    leg->AddEntry(hq1[i],degstr[i],"l");


  // Data with SF:OFF
    tup2[i] = GetTree(namestr2[i]);
    tup2[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // edited based on suggestion by KM on the August 3rd
    
    if(useMRfirst) kp2[i] = tup2[i]->GetMaximum("mp");
    else kp2[i] = (tup2[i]->GetMaximum("kp") - tup2[i]->GetMinimum("kp"));
    cout << "SF:OFF data # of kp: "<< kp2[i] <<endl;

    hq2[i] = new TH1F(Form("hq2%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    tup2[i]->Draw(Form("%f/(toffo*%f)>>hq2%d",twopirad2,lambda_coeff,i), thecut && cut_ref,"goff");
    hq2[i]->Scale(25./kp2[i]); 
    hq2[i]->Divide(hq0[i]); //Divide by the direct data
    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq2[i]->SetTitle("Reflectivity (SF OFF)");
    leg2->AddEntry(hq2[i],degstr2[i],"l");
    

// Calculate Polarization Power
  hpolratio[i]=(TH1F*)hq1[i]->Clone(Form("hpolratio%d",i));
  // hpolratio[i]->Add(hq[i], hq2[i],1., -1.); // hq: OFF, hq2: ON, calculate Non - Noff
  hpolratio[i]->Add(hq1[i], hq2[i],-1., 1.); // hq: ON, hq2: OFF, calculate Noff - Non
  hpolratio2[i]=(TH1F*)hq1[i]->Clone(Form("hpolratio2%d",i));
  // hpolratio2[i]->Add(hq[i], hq2[i],1., 1.); // hq: OFF, hq2: ON, calculate Non + Noff
  hpolratio2[i]->Add(hq1[i], hq2[i],1., 1.); // hq: OFF, hq2: ON, calculate Non + Noff
  hpolratio[i]->Divide(hpolratio2[i]);
  hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{off}-R_{on})/(R_{off}+R_{on})");
  hpolratio[i]->SetTitle("Polarization power");
  leg3->AddEntry(hpolratio[i], degstr_p[i],"l");


  if(i==9){
      //hx[i]->SetLineColor(i+2);
      //hlambda[i]->SetLineColor(i+2);
      //hratio[i]->SetLineColor(i+2);
      hq1[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hpolratio[i]->SetLineColor(i+2);
      hpolratio2[i]->SetLineColor(i+2);
    } else {
      //hx[i]->SetLineColor(i+1);
      //hlambda[i]->SetLineColor(i+1);
      //hratio[i]->SetLineColor(i+1);
      hq1[i]->SetLineColor(i+1);
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
      // TBox* b = new TBox(0.287,0,qm,2); 
      // b->SetFillColor( 7 ); 
      // b->SetFillStyle(3004); 
      // b->Draw(); 
      // TBox* b1 = new TBox(0.15,0,0.173,2); 
      // b1->SetFillColor( kOrange ); 
      // b1->SetFillStyle(3004); 
      // b1->Draw("same");
      // leg->Draw();
      hq1[i]->SetStats(0);  
      if(i==0)hq1[i]->Draw("eh");
      else hq1[i]->Draw("ehsames");
      leg->Draw();
      hq1[i]->SetStats(0);
    }

    c1->cd(2);
    if(i!=0){
     

      // TBox* b2 = new TBox(0.287,0,qm,2.); 
      // b2->SetFillColor( 7  ); 
      // b2->SetFillStyle(3004); 
      // b2->Draw();
      // TBox* b3 = new TBox(0.15,0,0.173,2.); 
      // b3->SetFillColor( kOrange); 
      // b3->SetFillStyle(3004); 
      // b3->Draw("same"); 
      
     if(i==0)hq2[i]->Draw("eh");
     else hq2[i]->Draw("ehsames");
    // if(i==1)hq2[i]->Draw("ah");
    // else hq2[i]->Draw("ahsames");
    leg2->Draw();
    hq2[i]->SetStats(0);
    
    
    }
    

    c1->cd(3);
    if(i!=0){
      // TBox* b4 = new TBox(0.15,-1.2,0.173,1.2); 
      // b4->SetFillColor( kOrange ); 
      // b4->SetFillStyle(3004); 
      // b4->Draw();
      // TBox* b5 = new TBox(0.287,-1.2,qm,1.2); 
      // b5->SetFillColor(7); 
      // b5->SetFillStyle(3004); 
      // b5->Draw("same");
      hpolratio[i]->SetStats(0);
      if(i==0)hpolratio[i]->Draw("eh");
      else hpolratio[i]->Draw("ehsames");

      leg3->Draw();
    }

    c1->cd(4);

    hpolratio[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hpolratio[i]->GetYaxis()->SetRangeUser(-1.2,1.2);
    hpolratio2[i]->GetYaxis()->SetRangeUser(-1.2,1.2);
    hq1[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq1[i]->GetYaxis()->SetRangeUser(1.e-3,2.);
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
  
  Double_t B[num];
  for(int i=0;i<num;i++){
    B[i] = vec_H[i];
  }
  Double_t pol_at_qcut[num];
  Double_t error_pol_at_qcut[num];

  Double_t ibin_pol[num_pol];

  c1->cd(4); 
  TLegend *leg4 = new TLegend(0.8, 0.8, 0.95, 1, "");
  // leg4->SetFillStyle(0);
  ofstream ofs(path_R+Form("poldata_%s.csv", scan_id.c_str()));  // ファイルパスを指定する

  for (Int_t i=0; i<num_pol; i++){
    ibin_pol[i]= Int_t((q_cuts[i]-q_min)*nbin_q/(q_max-q_min));
    Double_t q_cut_i = q_min + (q_max-q_min)*ibin_pol[i]/nbin_q;
    for (Int_t j=0; j<num; j++){
      pol_at_qcut[j] = hpolratio[j]->GetBinContent(ibin_pol[i]);
      error_pol_at_qcut[j] = hpolratio[j]->GetBinError(ibin_pol[i]);
      ofs << q_cut_i << "," << B[j] << ","<< pol_at_qcut[j] << ","<< error_pol_at_qcut[j] << endl;
      // cout <<   pol_at_qcut[j] << endl;
      // count << Form("At q=%.2f nm^{-1}, B=%.3f mT: Polarization power=", q_cuts[i], B[i])<<  pol_at_qcut[i] << endl;
    }
    gr[i]= new TGraphErrors(num, B, pol_at_qcut,0,error_pol_at_qcut);
    if (i==0) gr[i]->Draw("AP");
    else gr[i]->Draw("P");
    gr[i]->SetMarkerColor(i+1);
    gr[i]->SetLineColor(i+1);
    gr[i]->SetMarkerStyle(i+3);
    gr[i]->SetMarkerSize(1);
    gr[i]->GetXaxis()->SetRange(0,9);
    gr[i]->GetXaxis()->SetTitle("B (mT)");
    gr[i]->GetYaxis()->SetRangeUser(-1.2, 1.2);
    gr[i]->GetXaxis()->SetRange(0,9);
    gr[i]->GetXaxis()->SetRangeUser(0,9);
    gr[i]->GetYaxis()->SetTitle("Polarization power");
    gr[i]->SetTitle("");
    
    // leg4->AddEntry(gr[i],Form("q=%.3f nm^{-1}",q_cut_i));
    leg4->AddEntry(gr[i],Form("q=%.3f nm^{-1}",q_min + (q_max-q_min)*ibin_pol[i]/nbin_q),"p");
    // ofstream ofs(path_R+Form("30nm_long1_%.3f.csv",q_cut_i));  // ファイルパスを指定する
    // for(Int_t i=0; i<num-1; i++){
    //   ofs << B[i] << ","<< pol_at_qcut[i] << ","<< error_pol_at_qcut[i] << endl;
    // }
    // for(Int_t i=0; i<num; i++){
    //   ofs << q_cut_i << "," << B[i] << ","<< pol_at_qcut[i] << ","<< error_pol_at_qcut[i] << endl;
    // }
    
  }


  leg4->Draw();


    


  // double qq=0.25;
  // double nbin_qq=qq*nbin_q/(q_max-q_min);
  // entry[i]= hpolratio[i]->GetBinContent(nbin_qq); 
  // cout<<"[i]_"<<i<<"entry_"<<entry[i]<<endl;

  

  c1->cd(1); gPad->SetGrid();gPad->SetLogy();
  c1->cd(2); gPad->SetGrid();gPad->SetLogy();
  
  c1->cd(3); gPad->SetGrid();//gPad->SetLogy();
  c1->cd(4); gPad->SetGrid();//gPad->SetLogy();
    

  c1->SaveAs(path_R+Form("pol_%s.png", scan_id.c_str()));
  

#if 1
  TFile *outfile = TFile::Open(path_R+Form("pol_%s.root", scan_id.c_str()),"RECREATE");
  for(Int_t i=0; i<num; i++){
    //hx[i]->Write();
    //hlambda[i]->Write();
    //hratio[i]->Write();
    hq1[i]->Write();
    //hxylambda[i]->Write();
    hq2[i]->Write();
  }
  outfile->Close();
#endif


  return 0 ;
}
