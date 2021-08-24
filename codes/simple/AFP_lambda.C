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

Int_t AFP_lambda(){

  InitColor();
  TH1::SetDefaultSumw2();

    const Int_t num = 8;
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
  namestr[1]="20210714185238_list.root"; 
  namestr[2]="20210714185238_list.root"; 
  namestr[3]="20210714185238_list.root"; 
  namestr[4]="20210714185238_list.root";
  namestr[5]="20210714185238_list.root"; 
  namestr[6]="20210714185238_list.root"; 
  namestr[7]="20210714191741_list.root"; //off 276.2 deg

  //on
  namestr2[0]="20210714193654_list.root"; //M1 reflect (direct) 1hour
  namestr2[1]="20210714204714_list.root"; //100mV -8.01mT
  namestr2[2]="20210714205602_list.root"; //760mV
  namestr2[3]="20210714210221_list.root"; //1000mV
  namestr2[4]="20210714211037_list.root"; //500mV
  namestr2[5]="20210714211642_list.root"; //300mV
  namestr2[6]="20210714214337_list.root"; //760mV I_LV=0A
  namestr2[7]="20210714215803_list.root"; //760mV 276.2deg

  Double_t angle[num];
  Double_t angle2[num];
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle[1] = TMath::Abs(57.563 - xdirect)/dist_det; //rad 47.1868
  //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(57.563 - xdirect)/dist_det; //rad
  angle[3] = TMath::Abs(57.563 - xdirect)/dist_det; //rad
  angle[4] = TMath::Abs(57.563 - xdirect)/dist_det; //rad
  angle[5] = TMath::Abs(57.563- xdirect)/dist_det; //rad
  angle[6] = TMath::Abs(57.563- xdirect)/dist_det; //rad
  angle[7] = TMath::Abs(51.634 - xdirect)/dist_det; 

    /*
namestr[0]="20210714193654_list.root";
  namestr[1]="20210714185238_list.root"; 
  namestr[2]="20210714185238_list.root"; 
  namestr[3]="20210714185238_list.root"; 
  namestr[4]="20210714185238_list.root";
  namestr[5]="20210714185238_list.root"; 
  namestr[6]="20210714185238_list.root"; 
  namestr[7]="20210714191741_list.root"; //off 276.2 deg
*/

  angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle2[1] = TMath::Abs(57.8928 - xdirect)/dist_det; //rad
  angle2[2] = TMath::Abs(58.1188 - xdirect)/dist_det; //rad
  angle2[3] = TMath::Abs(58.1635 - xdirect)/dist_det; //rad
  angle2[4] = TMath::Abs(58.0323 - xdirect)/dist_det; //rad
  angle2[5] = TMath::Abs(57.9249 - xdirect)/dist_det; //rad
  angle2[6] = TMath::Abs(58.1081 - xdirect)/dist_det; //rad
  angle2[7] = TMath::Abs(51.6895- xdirect)/dist_det;
  /*
namestr2[0]="20210714193654_list.root"; //M1 reflect (direct) 1hour
  namestr2[1]="20210714204714_list.root"; //100mV -8.01mT
  namestr2[2]="20210714205602_list.root"; //760mV
  namestr2[3]="20210714210221_list.root"; //1000mV
  namestr2[4]="20210714211037_list.root"; //500mV
  namestr2[5]="20210714211642_list.root"; //300mV
  namestr2[6]="20210714214337_list.root"; //760mV I_LV=0A
  namestr2[7]="20210714215803_list.root"; //760mV 276.2deg
  */
  double angledeg[num];
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  angledeg[3]=angle[3]*180./TMath::Pi()/2.;
  angledeg[4]=angle[4]*180./TMath::Pi()/2.;
  angledeg[5]=angle[5]*180./TMath::Pi()/2.;
  angledeg[6]=angle[6]*180./TMath::Pi()/2.;
  angledeg[7]=angle[7]*180./TMath::Pi()/2.;

  TString degstr[num];
  TString degstr2[num];
  //off
  degstr[0]="Direct(M1 reflect)";
  degstr[1]="theta_m2 275.7 deg";
  degstr[2]=Form("theta_m2 275.7 deg %f_deg",angledeg[2]);
  degstr[3]="theta_m2 275.7 deg";
  degstr[4]="theta_m2 275.7 deg";
  degstr[5]="theta_m2 275.7 deg";
  degstr[6]="theta_m2 275.7 deg";
  degstr[7]=Form("theta_m2 276.2 deg %f_deg",angledeg[7]);

  //on
  degstr2[0]="Direct(M1 reflect)";
  degstr2[1]="SF-RF 100mV";
  degstr2[2]="SF-RF 760mV";
  degstr2[3]="SF-RF 1000mV";
  degstr2[4]="SF-RF 500mV";
  degstr2[5]="SF-RF 300mV";
  degstr2[6]="I_LV 0A";
  degstr2[7]="theta_m2 276.2 deg";


  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 30 nm");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 300;
  Double_t q_max  = 1.0;//0.6
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

    hratio[i]->GetYaxis()->SetRangeUser(0.,2.);
    hq[i]->GetYaxis()->SetRangeUser(0.1,0.9);
    hq[i]->GetYaxis()->SetRangeUser(0.,1.);
    hq[i]->GetXaxis()->SetRangeUser(0.1,0.9);

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
