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

Int_t Pol_Power_90nm(){

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

  TString namestr[num];
  TString namestr2[num];
  //off
  namestr[0]="20210714193654_list.root";
  namestr[1]="20210716232122_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  namestr[2]="20210717004515_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr[3]="20210717022140_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  //namestr[4]="20210715085349_list.root"; //0.264A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT
  //namestr[5]="20210715082606_list.root"; //0.378A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT  -> 2 mT
  //namestr[6]="20210715083711_list.root"; //0.6A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 2 mT  -> 3 mT

  //on
  namestr2[0]="20210714193654_list.root";
  namestr2[1]="20210716233530_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  namestr2[2]="20210717005023_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr2[3]="20210717022920_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  
 
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

  Double_t angledeg[num];
  Double_t angledeg2[num];
  
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  angledeg[3]=angle[3]*180./TMath::Pi()/2.;
  
  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm OFF");
  TLegend* leg2 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm ON");
  TLegend* leg3 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 30 nm (ON-OFF)/(OFF+ON)");
  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Double_t nbin_q  = 60;//300
  Double_t q_min  = 0.15;//0.6
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
    leg3->AddEntry(hq[i],degstr[i],"l");
    //leg->AddEntry(hq2[i],degstr2[i],"l2");

    //hx[i]->Scale(25./kp[i]);
    //hlambda[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    //hxylambda[i]->Scale(25./kp[i]);
    
    
    
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

    hq2[i]->Divide(hq02[i]);

    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq2[i]->SetTitle("Reflectivity (SF ON)");
    //hq02[i]->GetYaxis()->SetTitle("Reflectivity02");

    hpolratio[i]=(TH1F*)hq[i]->Clone(Form("hpolratio%d",i));
    hpolratio[i]->Add(hq[i], hq2[i],-1, 1); // hq: OFF, hq2: ON, calculate Non - Noff
    // hpolratio[i]->Add(hq2[i],-1.);//各binの値=x*hist1の値+y*hist2の値
    hpolratio2[i]=(TH1F*)hq[i]->Clone(Form("hpolratio2%d",i));
    hpolratio2[i]->Add(hq[i], hq2[i],1, 1); // hq: OFF, hq2: ON, calculate Non + Noff
    hpolratio[i]->Divide(hpolratio2[i]);
    hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{on}-R_{off})/(R_{on}+R_{off})");
    hpolratio[i]->SetTitle("Polarization power");

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
    c1->cd(1);
    if(i!=0){
    if(i==1)hq[i]->Draw("eh");
    else hq[i]->Draw("ehsames");
    // if(i==1)hq[i]->Draw("ah");
    // else hq[i]->Draw("ahsames");
    leg->Draw();
    }

    c1->cd(2);
    if(i!=0){
    if(i==1)hq2[i]->Draw("eh");
    else hq2[i]->Draw("ehsames");
    // if(i==1)hq2[i]->Draw("ah");
    // else hq2[i]->Draw("ahsames");
    leg2->Draw();
    }

    c1->cd(3);
    if(i!=0){
    if(i==1)hpolratio[i]->Draw("eh");
    else hpolratio[i]->Draw("ehsames");
    leg3->Draw();
    }
    
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
  
  Double_t B[num-1] = {8.01,0.322,1.35};
  // Double_t q_cuts[num_pol] = {0.2, 0.25, 0.3, 0.35, 0.4};
  Double_t pol_at_qcut[num-1];
  Double_t error_pol_at_qcut[num-1];

  Double_t ibin_pol[num_pol];

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

  }
  leg4->Draw();

  c1->cd(1); gPad->SetGrid();gPad->SetLogy();
  c1->cd(2); gPad->SetGrid();gPad->SetLogy();
  
  c1->cd(3); gPad->SetGrid();//gPad->SetLogy();
  // c1->cd(4); gPad->SetGrid();//gPad->SetLogy();

  c1->SaveAs(path_R+"pol_90nm.png");
  c1->SaveAs(path_R+"pol_90nm.root");

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
