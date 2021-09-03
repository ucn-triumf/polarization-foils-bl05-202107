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

int all_thickness() {
    // // Using TFile
    // TTree* T = new TTree("T", "T");
    // T->ReadFile("data_scans/magnetic_20210716210153.csv", "index/D:I/D:B/D");
    // Int_t num=T->GetEntries();
    // cout << num << endl;

    // Using vector
    vector<Int_t> vec_index;
    vector<Double_t> vec_q,vec_I, vec_H,vec_H_E; 
    vec_index.clear(); vec_I.clear(); vec_I.clear();
    ifstream fcsv("results/poldata_30nm_scan_1.csv"); 
    int i_csv = 0;
    
    if (!fcsv.is_open())
    {
        exit(EXIT_FAILURE);
    }
    string str;
    //getline(fcsv, str);
    cout << str << endl; // print out the first row
    while (getline(fcsv, str))
    {
    // while(!incsv.eof()) {
    // if (i_csv>0){
       //Int_t index;     
       Double_t qq, current, magfield, magfield_error;
        // cout << str << endl;
    //     double r = 0., z = 0., v =0.;
        sscanf(str.c_str(), "%lf,%lf,%lf,%lf", &qq, &current, &magfield, &magfield_error);
        // str >> index >> ",'" >> current >>",'">> magfield;
        //vec_index.push_back(index); 
        vec_q.push_back(qq); 
        vec_I.push_back(current);
        vec_H.push_back(magfield);
        vec_H_E.push_back(magfield_error);
        // cout << " " << index
    	//  << " " << current
        //   << " " << magfield << endl;
        // }
    // i_csv++;
  }
  

  int num=15;
  double q1[num],I1[num],H1[num],HE1[num];
  for(int i=0;i<num;i++){
    cout << vec_q[i]<< vec_I[i]<< vec_H[i]  << vec_H_E[i]<< endl;
    //double *q1[i]=&vec_q[i];

      
  }
  for(int i=0;i<num;i++){
    if(vec_q[i]==0.249167){
        //*
        TGraphErrors *gr1= new TGraphErrors(num,vec_I,vec_H,0,vec_H_E);
        gr1->SetMarkerColor(1);
        gr1->SetMarkerColor(1);
        gr1->SetLineColor(1);
        gr1->SetMarkerStyle(3);
        gr1->SetMarkerSize(1);
        gr1->GetXaxis()->SetRange(0,9);
        gr1->GetXaxis()->SetTitle("B (mT)");
        gr1->GetYaxis()->SetRangeUser(-1.2, 1.2);
        gr1->GetXaxis()->SetRange(0,9);
        gr1->GetXaxis()->SetRangeUser(0,9);
        gr1->GetYaxis()->SetTitle("Polarization power");
        gr1->SetTitle("");
        gr1->Draw("AP");
        ///*/
        /*
        gr[i]=new TGraphErrors(num,vec_I[i],vec_H[i],0.,vec_H_E[i]);
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
        */

    }
    if(vec_q[i]==0.295833){
        
    }
    else{

    }
  }

 return 0;   // ifstream infile("results/30nm_mT_P.csv")
}