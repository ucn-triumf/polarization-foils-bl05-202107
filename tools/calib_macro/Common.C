#ifndef COMMON_C_
#define COMMON_C_
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TColor.h>
#include <TAttLine.h>
#include <TAttFill.h>
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>

namespace Common {
void InitROOT() {
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH1::StatOverflows(kTRUE);
  gStyle->SetOptStat("irMen"); // "ksiourmen"
  // gStyle->SetOptStat("iourMen"); // "ksiourmen"
  gStyle->SetOptFit(1110);   // "pcev"
  //set default color
  gROOT->GetColor(2)->SetRGB(220./255.,  50./255.,  47./255.); // Red
  gROOT->GetColor(3)->SetRGB(135./255., 194./255.,  63./255.); // Green
  gROOT->GetColor(4)->SetRGB( 38./255., 139./255., 210./255.); // Blue
  gROOT->GetColor(5)->SetRGB(250./255., 202./255.,  18./255.); // Yellow
  gROOT->GetColor(6)->SetRGB(236./255.,   0./255., 140./255.); // Magenta
  gROOT->GetColor(7)->SetRGB(135./255., 206./255., 250./255.); // Cyan
  gROOT->GetColor(8)->SetRGB(102./255., 205./255., 170./255.); // Lightgreen
  //
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.13);
}

void SetColors(TObjArray *obj_array, Int_t i_first=0) {
  Int_t num_object = obj_array->GetEntries();
  Int_t light = 127;
  Int_t satu  = 200;
  for (Int_t i = i_first; i < num_object; i++) {
    Int_t hue = int(255/num_object*(i-i_first));
    Int_t red, green, blue;
    TColor::HLS2RGB(hue, light, satu, red, green, blue);
    Int_t color = TColor::GetColor(red, green, blue);
    dynamic_cast<TAttLine*>(obj_array->At(i))
      ->SetLineColor(color);
    dynamic_cast<TAttFill*>(obj_array->At(i))
      ->SetFillColor(color);
    dynamic_cast<TAttMarker*>(obj_array->At(i))
      ->SetMarkerColor(color);
  }
}

TCanvas* GetSquarePlot (TString name, TString title,
                        int w, float l, float r, float b, float t) {
  // https://root-forum.cern.ch/t/draw-square-th2/29135/6
  // This function creates a canvas making sure the plot is square.
  // The first parameter w is the canvas width and l, r, b and t are the
  // left, right, bottom and top magins. The canvas height is computed with these
  // parameters in order to have a square plot.
  int h = ((1.-(l+r))*w)/(1.-(b+t));
  TCanvas *c = new TCanvas(name,title,w,h);
  c->SetLeftMargin(l),
    c->SetRightMargin(r),
    c->SetBottomMargin(b),
    c->SetTopMargin(t);
  // c->Draw();
  return c;
}

// from ROOT6 TCanvas
void SetRealAspectRatio(TCanvas *c, const Int_t axis=1) {
  c->Update();
  //Get how many pixels are occupied by the canvas
  Int_t npx = c->GetWw();
  Int_t npy = c->GetWh();

  //Get x-y coordinates at the edges of the canvas (extrapolating outside the axes, NOT at the edges of the histogram)
  Double_t x1 = c->GetX1();
  Double_t y1 = c->GetY1();
  Double_t x2 = c->GetX2();
  Double_t y2 = c->GetY2();

  //Get the length of extrapolated x and y axes
  Double_t xlength2 = x2 - x1;
  Double_t ylength2 = y2 - y1;
  Double_t ratio2   = xlength2/ylength2;

  //Now get the number of pixels including the canvas borders
  Int_t bnpx = c->GetWindowWidth();
  Int_t bnpy = c->GetWindowHeight();

  if (axis==1) {
    c->SetCanvasSize(TMath::Nint(npy*ratio2), npy);
    c->SetWindowSize((bnpx-npx)+TMath::Nint(npy*ratio2), bnpy);
  } else if (axis==2) {
    c->SetCanvasSize(npx, TMath::Nint(npx/ratio2));
    c->SetWindowSize(bnpx, (bnpy-npy)+TMath::Nint(npx/ratio2));
  } else {
    std::cerr << "axis value " << axis
              << " is neither 1 (resize along x-axis) nor 2 (resize along y-axis)."
              << std::endl;
    exit(1);
  }

  //Check now that resizing has worked
  c->Update();

  //Get how many pixels are occupied by the canvas
  npx = c->GetWw();
  npy = c->GetWh();

  //Get x-y coordinates at the edges of the canvas (extrapolating outside the axes,
  //NOT at the edges of the histogram)
  x1 = c->GetX1();
  y1 = c->GetY1();
  x2 = c->GetX2();
  y2 = c->GetY2();

  //Get the length of extrapolated x and y axes
  xlength2 = x2 - x1;
  ylength2 = y2 - y1;
  ratio2 = xlength2/ylength2;

  //Check accuracy +/-1 pixel due to rounding
  if (!(abs(TMath::Nint(npy*ratio2) - npx)<2)) {
    std::cerr << "Resizing failed." << std::endl;
    // exit(1);
  }
}
}

#endif // COMMON_C_
