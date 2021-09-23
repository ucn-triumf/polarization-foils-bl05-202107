#include "../bin/MakeNiki.h"
//#include "/home/nop/data/Tools/GetMRcut.C"
#include <TH3.h>
TString path_R = "results/";

Double_t Distance = 20.000; //tentative
Double_t Conversion = 395.6;
Double_t dist_det   = 337.; //sample to detector [mm]
Double_t xdirect    = 64.71;
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

Int_t M1_R_192944(){

  TH1::SetDefaultSumw2();

  const Int_t num = 8;
  Int_t kp[num];
  TTree* tup[num];
  TH1F* hx[num];
  TH1F* hlambda[num];
  TH1F* hratio[num];
  TH1F* hq[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TString namestr[num];
  /*
  //  namestr[0]="20210714184125_list.root"; //M1 reflect (direct)
  namestr[0]="20210714193654_list.root"; //M1 reflect (direct) 1hour
  namestr[1]="20210714204714_list.root"; //M2 reflect theta = 0.49 deg.
  //  namestr[2]="20210714191741_list.root"; //M2 reflect theta = 1.00 deg.
  namestr[2]="20210714204714_list.root"; //M2 reflect theta = 0.49 deg. with AFP 100 mV
  namestr[3]="20210714211642_list.root"; //M2 reflect theta = 0.49 deg. with AFP 500 mV
  namestr[4]="20210714211037_list.root"; //M2 reflect theta = 0.49 deg. with AFP 500 mV
  namestr[5]="20210714205602_list.root"; //M2 reflect theta = 0.49 deg. with AFP 760 mV
  namestr[6]="20210714210221_list.root"; //M2 reflect theta = 0.49 deg. with AFP 1000 mV
  namestr[7]="20210714214337_list.root"; //M2 reflect theta = 0.49 deg. with AFP 760 mV magnet 0A
*/
  namestr[0]="20210713192944_list.root"; 
  namestr[1]="20210713202138_list.root"; //M1 reflect (direct) 1hour
  namestr[2]="20210714204714_list.root"; //100mV -8.01mT
  //namestr[2]="20210714205602_list.root"; //760mV
  //namestr[2]="20210713220438_list.root"; //760mV
  
  namestr[3]="20210714210221_list.root"; //1000mV
  namestr[4]="20210714211037_list.root"; //500mV
  namestr[5]="20210714211642_list.root"; //300mV
  namestr[6]="20210714214337_list.root"; //760mV I_LV=0A
  //namestr[7]="20210714215803_list.root"; //760mV 276.2deg
  namestr[7]="20210714191741_list.root"; //760mV 276.2deg
  

  TString degstr[num];
  degstr[0]="Direct(M1 trans)";
  degstr[1]="M1 reflect(0.454 deg.)";
  //  degstr[2]="M2 reflect(1.00 deg.)";
  degstr[2]="M2 reflect(0.48 deg.) with AFP 760 mV";
  degstr[3]="M2 reflect(0.49 deg.) with AFP 300 mV";
  degstr[4]="M2 reflect(0.49 deg.) with AFP 500 mV";
  degstr[5]="M2 reflect(0.49 deg.) with AFP 760 mV";
  degstr[6]="M2 reflect(0.49 deg.) with AFP 1000 mV";
  degstr[7]="M2 reflect(0.97 deg.) with AFP 760 mV ";

  Double_t angle[num];
  /*
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(58.1188 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(51. - xdirect)/dist_det; //rad
  */
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(57.8928 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(58.1188 - xdirect)/dist_det; //rad
  angle[3] = TMath::Abs(58.1635 - xdirect)/dist_det; //rad
  angle[4] = TMath::Abs(58.0323 - xdirect)/dist_det; //rad
  angle[5] = TMath::Abs(57.9249 - xdirect)/dist_det; //rad
  angle[6] = TMath::Abs(58.1081 - xdirect)/dist_det; //rad
  angle[7] = TMath::Abs(51.6895- xdirect)/dist_det;

  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 300;
  Double_t q_max  = 0.6;
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;
  //  Double_t HLD  = 7400.;

  //  Double_t xbegin=54.;
  /*Double_t xbegin=48.;
  Double_t xcenter=61.;
  Double_t xend=71.;
  Double_t ybegin=65.;
  Double_t yend=82.;*/
  Double_t xbegin=60.;
  Double_t xcenter=85.;
  Double_t xend=95.;
  Double_t ybegin=65.;
  Double_t yend=82.;

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
  gPad->SetLogy();

  // for(Int_t i=0; i<2; i++){
  for(Int_t i=0; i<num; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    tup[i] = GetTree(namestr[i]);
    // tup[i]->SetAlias("toffo","(tof>9.e3)*(tof)+(tof<9.e3)*(tof+40.e3)");
    tup[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // edited based on suggestion by KM on the August 3rd

    //    tup[i]->SetAlias("toffo","tof");
    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
    else kp[i] = tup[i]->GetMaximum("kp");
    cout << kp[i]<<endl;
    hx[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);

    tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    if(i==2)leg->AddEntry(hx[i],degstr[i],"l");
    //if(i==7)leg->AddEntry(hx[i],degstr[i],"l");

    hx[i]->Scale(25./kp[i]);
    hlambda[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    hxylambda[i]->Scale(25./kp[i]);

    hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    hratio[i]->Divide(hlambda[0]);
    hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    hq[i]->Divide(hq0[i]);
    hq[i]->GetYaxis()->SetTitle("Reflectivity");

    hratio[i]->GetYaxis()->SetRangeUser(0.,2.);
    hratio[i]->GetXaxis()->SetTitle("wave length [nm]");
    hq[i]->GetYaxis()->SetRangeUser(0.,2.);
    hq[i]->GetXaxis()->SetTitle("q [nm^-1]");

    if(i==3){
      hx[i]->SetLineColor(i+2);
      hlambda[i]->SetLineColor(i+2);
      hratio[i]->SetLineColor(i+2);
      hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
    } else {
      hx[i]->SetLineColor(i+1);
      hlambda[i]->SetLineColor(i+1);
      hratio[i]->SetLineColor(i+1);
      hq[i]->SetLineColor(i+1);
      hq0[i]->SetLineColor(i+1);
    }

    c1->cd(1);
    //if(i==0)hx[i]->Draw("eh");
    //else hx[i]->Draw("ehsames");
    if(i==0)hx[i]->Draw("eh");
    if(i==2)hx[i]->Draw("ehsames");
    //if(i==7)hx[i]->Draw("ehsames");
    leg->Draw();

    c1->cd(2);
    //if(i==0)hlambda[i]->Draw("eh");
    //else hlambda[i]->Draw("ehsames");
    if(i==0)hlambda[i]->Draw("eh");
    if(i==2)hlambda[i]->Draw("ehsames");
    //if(i==7)hlambda[i]->Draw("ehsames");
    leg->Draw();

    c1->cd(3);
    if(i==2) hratio[i]->Draw("eh");
    //if(i==7) hratio[i]->Draw("ehsames");
    leg->Draw();

    c1->cd(4);
    if(i==2) hq[i]->Draw("eh");
    //if(i==7) hq[i]->Draw("ehsames");
    leg->Draw();
  }

  
  //  hq[1]->GetYaxis()->SetRangeUser(0.,2.);
  c1->cd(2); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
  c1->cd(3); gPad->SetGridx(); gPad->SetGridy();
  c1->cd(4); gPad->SetGridx(); gPad->SetGridy();

  c1->SaveAs(Form(path_R+"AFP.png"));
  c1->SaveAs(Form(path_R+"AFP.root"));

#if 1
  TFile *outfile = TFile::Open("AFPhist.root","RECREATE");
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
