#include "/home/nikiglass/bin/MakeNiki.h"
#include "/home/nikiglass/data/Tools/GetMRcut.C"
#include <TH3.h>

Double_t Distance = 17.880;
Double_t Conversion = 395.6;
Double_t dist_det   = 203.; //sample to detector [mm]
Double_t xdirect    = 51.25;
Bool_t useMRfirst = 1; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.

TTree* GetTree(TString filestr){

  TString ROOTstr = filestr(0,19);
  ROOTstr += ".root";

  TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}

Int_t *UCNMirror(){

  TH1::SetDefaultSumw2();
  
  const Int_t num = 11;
  Int_t kp[num];
  TTree* tup[num];
  TH1F* hx[num];
  TH1F* hlambda[num];
  TH1F* hratio[num];
  TH1F* hq[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TString namestr[num];
  //    namestr[0]="20200410121622_list.root"; //direct short
  //    namestr[0]="20200412172130_list.root"; //direct long
  namestr[0]="20200414191516_list.root"; //direct det. pos. correct
  namestr[1]="20200409195419_list.root"; //Crystal1L
  namestr[2]="20200415192018_list.root"; //Crystal1T
  namestr[3]="20200414104521_list.root"; //Crystal2L
  namestr[4]="20200418184223_list.root"; //Crystal2T
  namestr[5]="20200413174905_list.root"; //Nichizo 
  namestr[6]="20200410184756_list.root"; //TRIUMF1
  namestr[7]="20200417125126_list.root"; //TRIUMF2
  namestr[8]="20200416140927_list.root"; //TRIUMF3
  namestr[9]="20200411162553_list.root"; //TRIUMF1 +0.4deg
  namestr[10]="20200419013408_list.root"; //Background

  TString degstr[num];
  degstr[0]="Direct";
  degstr[1]="Plate Crystal1L: 20mrad";
  degstr[2]="Plate Crystal1T: 20mrad";
  degstr[3]="Plate Crystal2L: 20mrad";
  degstr[4]="Plate Crystal2T: 20mrad";
  degstr[5]="Cylinder Nichizo: 20mrad";
  degstr[6]="Cylinder TRIUMF_S1: 20mrad";
  degstr[7]="Cylinder TRIUMF_S3: 20mrad";
  degstr[8]="Cylinder TRIUMF_S4: 20mrad";
  degstr[9]="Cylinder TRIUMF_S1: 27mrad";
  degstr[10]="Background";

  Double_t angle[num];
  angle[0] = (55.4 - xdirect)/dist_det; //rad
  angle[1] = (55.4 - xdirect)/dist_det; //rad
  angle[2] = (55.4 - xdirect)/dist_det; //rad
  angle[3] = (55.6 - xdirect)/dist_det; //rad 
  angle[4] = (55.4 - xdirect)/dist_det; //rad
  angle[5] = (55.4 - xdirect)/dist_det; //rad
  angle[6] = (55.4 - xdirect)/dist_det; //rad
  angle[7] = (55.4 - xdirect)/dist_det; //rad
  angle[8] = (55.4 - xdirect)/dist_det; //rad
  angle[9] = (57.5 - xdirect)/dist_det; //rad 
  //  angle[9] = (55.4 - xdirect)/dist_det + 0.4*TMath::DegToRad(); //rad
  angle[10] = (55.4 - xdirect)/dist_det; //rad

  TLegend* leg = new TLegend(0.5, 0.75, 0.75, 0.98,"");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 300;
  Double_t lambda_max  = 3.;
  Double_t nbin_q  = 300;
  Double_t q_max  = 0.6;
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;
  //  Double_t HLD  = 7400.;
  
  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     LLD,LLD,LLD,LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  TCut cut_x = Form("x*%f>20 && x*%f<100",range,range);
  TCut cut_y = Form("y*%f>62 && y*%f<67",range,range);
  TCut cut_dir = Form("x*%f<53.5 && x*%f>40.",range,range);
  TCut cut_ref = Form("x*%f>53.5 && x*%f<65.",range,range);
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

  for(Int_t i=0; i<num; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    tup[i] = GetTree(namestr[i]);

    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
    else kp[i] = tup[i]->GetMaximum("kp");
    cout << kp[i]<<endl;
    hx[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,0.,q_max);
    hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);

    tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    tup[i]->Draw(Form("tof*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(tof*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_ref,"goff");
    tup[i]->Draw(Form("%f/(tof*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[i]->Draw(Form("tof*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    leg->AddEntry(hx[i],degstr[i],"l");
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

    if(i==9){
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
    if(i==0)hx[i]->Draw("eh");
    else hx[i]->Draw("ehsames");
    leg->Draw();
    c1->cd(2);
    if(i==0)hlambda[i]->Draw("eh");
    else hlambda[i]->Draw("ehsames");    
    leg->Draw();

    c1->cd(3);
    if(i==1)hratio[i]->Draw("eh");
    else hratio[i]->Draw("ehsames");    
    leg->Draw();

    c1->cd(4);
    if(i==1)hq[i]->Draw("eh");
    else hq[i]->Draw("ehsames");    
    leg->Draw();
  }

  hratio[1]->GetYaxis()->SetRangeUser(0.,2.);
  //  hq[1]->GetYaxis()->SetRangeUser(0.,2.);
  c1->cd(2); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
  c1->cd(3); gPad->SetGridx(); gPad->SetGridy();
  c1->cd(4); gPad->SetGridx(); gPad->SetGridy();

  c1->SaveAs(Form("UCNMirror.png"));
  c1->SaveAs(Form("UCNMirror.root"));

#if 1
  TFile *outfile = TFile::Open("UCNMirrorHists.root","RECREATE");
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

