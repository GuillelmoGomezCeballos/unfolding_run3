#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include "TLegend.h"
#include <iostream>
#include "CMS_lumi.C"
#include "TTree.h"

Float_t GetMaximumIncludingErrors(TH1D* h, bool doApplyBinWidth)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if(doApplyBinWidth) binHeight = binHeight/h->GetBinWidth(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}

void eraselabel(TPad *p,Double_t h){
  p->cd();
  TPad *pe = new TPad("pe","pe",0.02,0,p->GetLeftMargin()-0.007,h);
  pe->Draw();
  pe->SetFillColor(p->GetFillColor()); 
  pe->SetBorderMode(0);
}

void atributes(TH1D *histo, TString xtitle = "", TString ytitle = "Fraction", TString units = "", Int_t color = 1){

  histo->SetTitle("");
  //histo->SetMarkerStyle(20);
  //histo->SetMarkerSize(0.8);
  //histo->SetLineWidth(4);
  if(strcmp(units.Data(),"")==0){
    histo->GetXaxis()->SetTitle(xtitle.Data());
  } else {
    histo->GetXaxis()->SetTitle(Form("%s [%s]",xtitle.Data(),units.Data()));
  }
  histo->GetXaxis()->SetLabelFont  (   42);
  histo->GetXaxis()->SetLabelOffset(0.015);
  histo->GetXaxis()->SetLabelSize  (0.110);
  histo->GetXaxis()->SetNdivisions (  505);
  histo->GetXaxis()->SetTitleFont  (   42);
  histo->GetXaxis()->SetTitleOffset(  0.9);
  histo->GetXaxis()->SetTitleSize  (0.150);
  histo->GetXaxis()->SetTickLength (0.07 );

  histo->GetYaxis()->SetTitle(ytitle.Data());
  histo->GetYaxis()->SetLabelFont  (   42);
  histo->GetYaxis()->SetLabelOffset(0.015);
  histo->GetYaxis()->SetLabelSize  (0.160);
  histo->GetYaxis()->SetNdivisions (  505);
  histo->GetYaxis()->SetTitleFont  (   42);
  histo->GetYaxis()->SetTitleOffset(  0.7);
  histo->GetYaxis()->SetTitleSize  (0.110);
  histo->GetYaxis()->SetTickLength (0.03 );
  histo->SetFillColor(color);
  histo->SetMarkerColor(color);
  histo->SetFillStyle(3345);
  histo->SetMarkerSize(0.8);
  histo->SetMarkerStyle(kFullCircle);
}

void finalPlotVBSVVUnfolding(TString keyLabel0 = "MLL", bool isNormalized = false) {

  TString XTitle = "X";
  TString units = "GeV";
  bool isLogY = false;
  bool isLogX = false;

  TString theYTitle = "#sigma / GeV [pb]";
  if     (!isNormalized && keyLabel0.Contains("MJJ"))    theYTitle = "d#sigma/dm_{jj} [fb]";
  else if( isNormalized && keyLabel0.Contains("MJJ"))    theYTitle = "1/#sigma d#sigma/dm_{jj} [1/bin]";
  else if(!isNormalized && keyLabel0.Contains("MLL"))    theYTitle = "d#sigma/dm_{ll} [fb]";
  else if( isNormalized && keyLabel0.Contains("MLL"))    theYTitle = "1/#sigma d#sigma/dm_{ll} [1/bin]";
  else if(!isNormalized && keyLabel0.Contains("NJET"))   theYTitle = "d#sigma/dN_{j} [fb]";
  else if( isNormalized && keyLabel0.Contains("NJET"))   theYTitle = "1/#sigma d#sigma/dN_{j}";
  else if(!isNormalized && keyLabel0.Contains("DELTAETAJJ")) theYTitle = "d#sigma/d#Delta#eta_{jj} [fb]";
  else if( isNormalized && keyLabel0.Contains("DELTAETAJJ")) theYTitle = "1/#sigma d#sigma/d#Delta#eta_{jj} [1/bin]";
  else if(!isNormalized && keyLabel0.Contains("DELTAPHIJJ")) theYTitle = "d#sigma/d#Delta#phi_{jj} [fb]";
  else if( isNormalized && keyLabel0.Contains("DELTAPHIJJ")) theYTitle = "1/#sigma d#sigma/d#Delta#phi_{jj} [1/bin]";
  else {printf("PROBLEM!\n"); return;}

  if     (keyLabel0 == "EWKWZMJJ"       ) {XTitle = "m_{jj}";}
  else if(keyLabel0 == "EWKWWMJJ"       ) {XTitle = "m_{jj}";}
  else if(keyLabel0 == "EWKWWMLL"       ) {XTitle = "m_{ll}";}
  else if(keyLabel0 == "EWKWWNJET"      ) {XTitle = "Number of jets"; units = "";}
  else if(keyLabel0 == "EWKWWDELTAETAJJ") {XTitle = "#Delta#eta_{jj}"; units = "";}
  else if(keyLabel0 == "EWKWWDELTAPHIJJ") {XTitle = "#Delta#phi_{jj}"; units = "";}

  gInterpreter->ExecuteMacro("PaperStyle.C");
  gStyle->SetOptStat(0);

  bool isDebug = true;

  double scaleDFSF = 1.0;

  const int theFillColor1 = 12;
  const int theFillStyle1 = 3345;
  const int theFillColor2 = 27;
  const int theFillStyle2 = 3004;
  const int theFillColor3 = 46;
  const int theFillStyle3 = 3007;

  TString genFileName1 = "histogen_vbsvv_mapgraph_noewkcorr.root";
  TFile *_fileGenWW1 = TFile::Open(genFileName1.Data());
  TH1D* hPred1     = (TH1D*)_fileGenWW1->Get(Form("hD%s",keyLabel0.Data()));     hPred1    ->SetDirectory(0);
  TH1D* hPred1_PDF = (TH1D*)_fileGenWW1->Get(Form("hD%s_PDF",keyLabel0.Data())); hPred1_PDF->SetDirectory(0);
  TH1D* hPred1_QCD = (TH1D*)_fileGenWW1->Get(Form("hD%s_QCD",keyLabel0.Data())); hPred1_QCD->SetDirectory(0);
  TH1D* hPred1_PS  = (TH1D*)_fileGenWW1->Get(Form("hD%s_PS",keyLabel0.Data()));  hPred1_PS ->SetDirectory(0);
  _fileGenWW1->Close();

  TString genFileName2 = "histogen_vbsvv_mapgraph.root";
  TFile *_fileGenWW2 = TFile::Open(genFileName2.Data());
  TH1D* hPred2     = (TH1D*)_fileGenWW2->Get(Form("hD%s",keyLabel0.Data()));     hPred2    ->SetDirectory(0);
  TH1D* hPred2_PDF = (TH1D*)_fileGenWW2->Get(Form("hD%s_PDF",keyLabel0.Data())); hPred2_PDF->SetDirectory(0);
  TH1D* hPred2_QCD = (TH1D*)_fileGenWW2->Get(Form("hD%s_QCD",keyLabel0.Data())); hPred2_QCD->SetDirectory(0);
  TH1D* hPred2_PS  = (TH1D*)_fileGenWW2->Get(Form("hD%s_PS",keyLabel0.Data()));  hPred2_PS ->SetDirectory(0);
  _fileGenWW2->Close();

  TString genFileName3 = "histogen_vbsvv_sherpa_withWZ.root";
  TFile *_fileGenWW3 = TFile::Open(genFileName3.Data());
  TH1D* hPred3     = (TH1D*)_fileGenWW3->Get(Form("hD%s",keyLabel0.Data()));     hPred3    ->SetDirectory(0);
  TH1D* hPred3_PDF = (TH1D*)_fileGenWW3->Get(Form("hD%s_PDF",keyLabel0.Data())); hPred3_PDF->SetDirectory(0);
  TH1D* hPred3_QCD = (TH1D*)_fileGenWW3->Get(Form("hD%s_QCD",keyLabel0.Data())); hPred3_QCD->SetDirectory(0);
  TH1D* hPred3_PS  = (TH1D*)_fileGenWW3->Get(Form("hD%s_PS",keyLabel0.Data()));  hPred3_PS ->SetDirectory(0);
  _fileGenWW3->Close();

  TString plotName = Form("input_files/xs_%s_normalized%d.root", keyLabel0.Data(), isNormalized);
  TFile *_file0 = TFile::Open(plotName.Data());
  TH1D* hData = (TH1D*)_file0->Get(Form("hD%s",keyLabel0.Data())); hData->SetDirectory(0);
  _file0->Close();

  for(Int_t i=1;i<=hPred1->GetNbinsX();++i){

    hData->SetBinContent(i,hData->GetBinContent(i)*hPred2->GetBinContent(i));
    hData->SetBinError  (i,hData->GetBinError  (i)*hPred2->GetBinContent(i));

    // Begin Pred1
    double diff[5] = {hPred1->GetBinError(i)/hPred1->GetBinContent(i),
                      TMath::Abs(1-hPred1_PDF->GetBinContent(i)/hPred1->GetBinContent(i)),
                      TMath::Abs(1-hPred1_QCD->GetBinContent(i)/hPred1->GetBinContent(i)),
                      TMath::Abs(1-hPred1_PS ->GetBinContent(i)/hPred1->GetBinContent(i)),
                      0.0};

    if(isNormalized) {
      diff[1] = TMath::Abs(1-(hPred1_PDF->GetBinContent(i)/hPred1_PDF->GetSumOfWeights())/(hPred1->GetBinContent(i)/hPred1->GetSumOfWeights()));
      diff[2] = TMath::Abs(1-(hPred1_QCD->GetBinContent(i)/hPred1_QCD->GetSumOfWeights())/(hPred1->GetBinContent(i)/hPred1->GetSumOfWeights()));
      diff[3] = TMath::Abs(1-(hPred1_PS ->GetBinContent(i)/hPred1_PS ->GetSumOfWeights())/(hPred1->GetBinContent(i)/hPred1->GetSumOfWeights()));
    }

    hPred1->SetBinError(i,sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3]+diff[4]*diff[4])*hPred1->GetBinContent(i));
    if(isDebug) printf("hPredSyst1 (%2d) %5.2f %5.2f %5.2f %5.2f %5.2f -> %5.2f\n",i,100*diff[0],100*diff[1],100*diff[2],100*diff[3],100*diff[4],100*hPred1->GetBinError(i)/hPred1->GetBinContent(i));
    // End Pred1

    // Pred2
    diff[0] = hPred2->GetBinError(i)/hPred2->GetBinContent(i);
    diff[1] = TMath::Abs(1-hPred2_PDF->GetBinContent(i)/hPred2->GetBinContent(i));
    diff[2] = TMath::Abs(1-hPred2_QCD->GetBinContent(i)/hPred2->GetBinContent(i));
    diff[3] = TMath::Abs(1-hPred2_PS ->GetBinContent(i)/hPred2->GetBinContent(i));
    diff[4] = 0.0;

    if(isNormalized) {
      diff[1] = TMath::Abs(1-(hPred2_PDF->GetBinContent(i)/hPred2_PDF->GetSumOfWeights())/(hPred2->GetBinContent(i)/hPred2->GetSumOfWeights()));
      diff[2] = TMath::Abs(1-(hPred2_QCD->GetBinContent(i)/hPred2_QCD->GetSumOfWeights())/(hPred2->GetBinContent(i)/hPred2->GetSumOfWeights()));
      diff[3] = TMath::Abs(1-(hPred2_PS ->GetBinContent(i)/hPred2_PS ->GetSumOfWeights())/(hPred2->GetBinContent(i)/hPred2->GetSumOfWeights()));
    }

    hPred2->SetBinError(i,sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3]+diff[4]*diff[4])*hPred2->GetBinContent(i));
    if(isDebug) printf("hPredSyst2 (%2d) %5.2f %5.2f %5.2f %5.2f %5.2f -> %5.2f\n",i,100*diff[0],100*diff[1],100*diff[2],100*diff[3],100*diff[4],100*hPred2->GetBinError(i)/hPred2->GetBinContent(i));
    // End Pred2

    // Pred3
    diff[0] = hPred3->GetBinError(i)/hPred3->GetBinContent(i);
    diff[1] = TMath::Abs(1-hPred3_PDF->GetBinContent(i)/hPred3->GetBinContent(i));
    diff[2] = TMath::Abs(1-hPred3_QCD->GetBinContent(i)/hPred3->GetBinContent(i));
    diff[3] = TMath::Abs(1-hPred3_PS ->GetBinContent(i)/hPred3->GetBinContent(i));
    diff[4] = 0.0;

    if(isNormalized) {
      diff[1] = TMath::Abs(1-(hPred3_PDF->GetBinContent(i)/hPred3_PDF->GetSumOfWeights())/(hPred3->GetBinContent(i)/hPred3->GetSumOfWeights()));
      diff[2] = TMath::Abs(1-(hPred3_QCD->GetBinContent(i)/hPred3_QCD->GetSumOfWeights())/(hPred3->GetBinContent(i)/hPred3->GetSumOfWeights()));
      diff[3] = TMath::Abs(1-(hPred3_PS ->GetBinContent(i)/hPred3_PS ->GetSumOfWeights())/(hPred3->GetBinContent(i)/hPred3->GetSumOfWeights()));
    }

    hPred3->SetBinError(i,sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]+diff[3]*diff[3]+diff[4]*diff[4])*hPred3->GetBinContent(i));
    if(isDebug) printf("hPredSyst3 (%2d) %5.2f %5.2f %5.2f %5.2f %5.2f -> %5.2f\n",i,100*diff[0],100*diff[1],100*diff[2],100*diff[3],100*diff[4],100*hPred3->GetBinError(i)/hPred3->GetBinContent(i));
    // End Pred3
  }

  hData ->Scale(scaleDFSF);
  hPred1->Scale(scaleDFSF);
  hPred2->Scale(scaleDFSF);
  hPred3->Scale(scaleDFSF);

  //hData ->Scale(1,"width");
  //hPred1->Scale(1,"width");
  //hPred2->Scale(1,"width");
  //hPred3->Scale(1,"width");

  Int_t ww = 800;
  Int_t wh = 800;

  TCanvas *c1 = new TCanvas("c1", "c1", ww, wh);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.350, 1, 0.975);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.000, 1, 0.345);

  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.00);  // 0.02

  pad2->SetTopMargin   (0.05);  // 0.08
  pad2->SetBottomMargin(0.30);  // 0.35

  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  gStyle->SetOptStat(0);
  if(isLogY == true) pad1->SetLogy();
  if(isLogX == true) {pad1->SetLogx(); pad2->SetLogx();}

  // draw pad1
  if(strcmp(units.Data(),"")==0){
    hPred1->GetXaxis()->SetTitle(XTitle.Data());
    hPred1->GetXaxis()->SetLabelOffset(0.005);
    hPred1->GetXaxis()->SetTitleOffset(  0.9);
    hPred2->GetXaxis()->SetTitle(XTitle.Data());
    hPred2->GetXaxis()->SetLabelOffset(0.005);
    hPred2->GetXaxis()->SetTitleOffset(  0.9);
    hPred3->GetXaxis()->SetTitle(XTitle.Data());
    hPred3->GetXaxis()->SetLabelOffset(0.005);
    hPred3->GetXaxis()->SetTitleOffset(  0.9);
  } else {
    hPred1->GetXaxis()->SetTitle(Form("%s [%s]",XTitle.Data(),units.Data()));
    hPred1->GetXaxis()->SetLabelOffset(0.00);
    hPred1->GetXaxis()->SetTitleOffset(  1.1);
    hPred2->GetXaxis()->SetTitle(Form("%s [%s]",XTitle.Data(),units.Data()));
    hPred2->GetXaxis()->SetLabelOffset(0.00);
    hPred2->GetXaxis()->SetTitleOffset(  1.1);
    hPred3->GetXaxis()->SetTitle(Form("%s [%s]",XTitle.Data(),units.Data()));
    hPred3->GetXaxis()->SetLabelOffset(0.00);
    hPred3->GetXaxis()->SetTitleOffset(  1.1);
  }

  hPred1->GetYaxis()->SetTitle(theYTitle.Data());
  hPred1->GetYaxis()->SetLabelFont  (   42);
  hPred1->GetYaxis()->SetLabelOffset(0.015);
  hPred1->GetYaxis()->SetLabelSize  (0.050);
  hPred1->GetYaxis()->SetNdivisions (  505);
  hPred1->GetYaxis()->SetTitleFont  (   42);
  hPred1->GetYaxis()->SetTitleOffset(  1.0);
  hPred1->GetYaxis()->SetTitleSize  (0.080);
  hPred1->GetYaxis()->SetTickLength (0.03 );
  hPred1->GetXaxis()->SetLabelFont  (   42);
  hPred1->GetXaxis()->SetLabelSize  (0.040);
  hPred1->GetXaxis()->SetNdivisions (  505);
  hPred1->GetXaxis()->SetTitleFont  (   42);
  hPred1->GetXaxis()->SetTitleSize  (0.060);
  hPred1->GetXaxis()->SetTickLength (0.07 );
 
  hPred2->GetYaxis()->SetTitle(theYTitle.Data());
  hPred2->GetYaxis()->SetLabelFont  (   42);
  hPred2->GetYaxis()->SetLabelOffset(0.015);
  hPred2->GetYaxis()->SetLabelSize  (0.050);
  hPred2->GetYaxis()->SetNdivisions (  505);
  hPred2->GetYaxis()->SetTitleFont  (   42);
  hPred2->GetYaxis()->SetTitleOffset(  1.2);
  hPred2->GetYaxis()->SetTitleSize  (0.060);
  hPred2->GetYaxis()->SetTickLength (0.03 );
  hPred2->GetXaxis()->SetLabelFont  (   42);
  hPred2->GetXaxis()->SetLabelSize  (0.040);
  hPred2->GetXaxis()->SetNdivisions (  505);
  hPred2->GetXaxis()->SetTitleFont  (   42);
  hPred2->GetXaxis()->SetTitleSize  (0.060);
  hPred2->GetXaxis()->SetTickLength (0.07 );
 
  hPred3->GetYaxis()->SetTitle(theYTitle.Data());
  hPred3->GetYaxis()->SetLabelFont  (   42);
  hPred3->GetYaxis()->SetLabelOffset(0.015);
  hPred3->GetYaxis()->SetLabelSize  (0.050);
  hPred3->GetYaxis()->SetNdivisions (  505);
  hPred3->GetYaxis()->SetTitleFont  (   42);
  hPred3->GetYaxis()->SetTitleOffset(  1.2);
  hPred3->GetYaxis()->SetTitleSize  (0.060);
  hPred3->GetYaxis()->SetTickLength (0.03 );
  hPred3->GetXaxis()->SetLabelFont  (   42);
  hPred3->GetXaxis()->SetLabelSize  (0.040);
  hPred3->GetXaxis()->SetNdivisions (  505);
  hPred3->GetXaxis()->SetTitleFont  (   42);
  hPred3->GetXaxis()->SetTitleSize  (0.060);
  hPred3->GetXaxis()->SetTickLength (0.07 );
 
  hData->SetFillColor(12);
  hData->SetFillStyle(3003);
  hData->GetYaxis()->SetTitleFont(42);
  hData->GetYaxis()->SetLabelFont(42);
  hData->GetXaxis()->SetTitleFont(42);
  hData->GetXaxis()->SetLabelFont(42);
  hData->GetYaxis()->SetTitleSize(0.055);
  hData->GetXaxis()->SetTitleSize(0.055);
  hData->GetYaxis()->SetLabelSize(0.039);
  hData->GetXaxis()->SetLabelSize(0.039);
  hData->GetYaxis()->SetTitleOffset(1.30);
  hData->GetXaxis()->SetTitleOffset(1.20);
  hData->GetYaxis()->SetLabelOffset(0.015);
  hData->GetXaxis()->SetLabelOffset(0.015);

  hPred1->SetLineColor(kRed+2);
  hPred1->SetMarkerStyle(3);
  hPred1->SetMarkerColor(kRed+2);
  hPred1->SetLineWidth(3);

  hPred2->SetLineColor(kBlue);
  hPred2->SetLineStyle(2);
  hPred2->SetMarkerStyle(5);
  hPred2->SetMarkerColor(kBlue);
  hPred2->SetLineWidth(3);

  hPred3->SetLineColor(kMagenta);
  hPred3->SetLineStyle(3);
  hPred3->SetMarkerStyle(6);
  hPred3->SetMarkerColor(kMagenta);
  hPred3->SetLineWidth(3);

  TAxis *xa = hData->GetXaxis();
  hPred1->SetTitle("");
  hPred2->SetTitle("");
  hPred3->SetTitle("");
  hData ->SetTitle("");
  double normalization[4] = {1.0, 1.0, 1.0, 1.0};
  if(isNormalized) {normalization[0] = hData->GetSumOfWeights(); normalization[1] = hPred1->GetSumOfWeights(); normalization[2] = hPred2->GetSumOfWeights(); normalization[3] = hPred3->GetSumOfWeights();};
  hData ->Scale(1./normalization[0]);
  hPred1->Scale(1./normalization[1]);
  hPred2->Scale(1./normalization[2]);
  hPred3->Scale(1./normalization[3]);

  Float_t dataMax = GetMaximumIncludingErrors(hPred1,false);
  hPred1->SetMaximum(1.40 * dataMax);

  double theEdges[2] = {TMath::Min(hPred1->GetMinimum(),hPred2->GetMinimum()),
                        TMath::Max(hPred1->GetMaximum(),hPred2->GetMaximum())};
  if(isLogY == true) hPred1->GetYaxis()->SetRangeUser(theEdges[0]/10,theEdges[1]*100);
  else               hPred1->GetYaxis()->SetRangeUser(0.0,theEdges[1]*1.5);
  hPred1->Draw("hist");
  hPred2->Draw("hist,same");
  if(keyLabel0 != "EWKWZMJJ") hPred3->Draw("hist,same");
  hData->Draw("ep,same");

  gStyle->SetOptStat(0);
  TLegend* legend = new TLegend(0.20,0.65,0.80,0.85);
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   62);
  legend->SetTextSize  (0.040);
  legend->AddEntry(hData,  "Data", "pfl");
  legend->AddEntry(hPred1, "MADGRAPH5_aMC@NLO+Pythia8 without NLO corr.", "l");
  legend->AddEntry(hPred2, "MADGRAPH5_aMC@NLO+Pythia8 with NLO corr.", "l");
  if(keyLabel0 != "EWKWZMJJ") legend->AddEntry(hPred3, "SHERPA 2", "l");

  bool plotSystErrorBars = true;
  if(plotSystErrorBars == true) {
    TGraphAsymmErrors * gsyst1 = new TGraphAsymmErrors(hData);
    for (int i = 0; i < gsyst1->GetN(); ++i) {
      gsyst1->SetPointEYlow (i,hData->GetBinError(i+1));
      gsyst1->SetPointEYhigh(i,hData->GetBinError(i+1));
    }
    gsyst1->SetFillColor(12);
    gsyst1->SetFillStyle(3003);
    gsyst1->SetMarkerSize(0);
    gsyst1->SetLineWidth(0);
    gsyst1->SetLineColor(kWhite);
    gsyst1->Draw("E2same");

    TGraphAsymmErrors * gsyst2 = new TGraphAsymmErrors(hData);
    for (int i = 0; i < gsyst2->GetN(); ++i) {
      gsyst2->SetPointEYlow (i,hData->GetBinError(i+1));
      gsyst2->SetPointEYhigh(i,hData->GetBinError(i+1));
    }
    gsyst2->SetFillColor(12);
    gsyst2->SetFillStyle(3003);
    gsyst2->SetMarkerSize(0);
    gsyst2->SetLineWidth(0);
    gsyst2->SetLineColor(kWhite);
    //gsyst2->Draw("E2same");
    //legend->AddEntry(gsyst1, "Theoretical uncertainty", "f");
  }
  legend->Draw();
  // plotting again
  hPred1->Draw("hist,same");
  hPred2->Draw("hist,same");
  if(keyLabel0 != "EWKWZMJJ") hPred3->Draw("hist,same");
  hData->Draw("ep,same");

  CMS_lumi( pad1, 2027, 11 );

  pad2->cd();
  gStyle->SetOptStat(0);

  TH1D* hNum1 = (TH1D*) hPred1->Clone(); hNum1->Reset();
  TH1D* hNum2 = (TH1D*) hPred2->Clone(); hNum2->Reset();
  TH1D* hNum3 = (TH1D*) hPred3->Clone(); hNum3->Reset();
  TH1D* hDen  = (TH1D*) hData ->Clone(); hDen->Reset();

  hNum1->Add(hPred1);
  hNum2->Add(hPred2);
  hNum3->Add(hPred3);
  hDen ->Add(hData);

  double pull; 
  double pullerr;
  double pullinv; 
  double pullinverr;

  // hPred1 w.r.t. hData
  TH1D* hRatio1 = (TH1D*) hPred1->Clone(); hRatio1->Reset();
  TH1D* hBand   = (TH1D*) hData ->Clone(); hBand  ->Reset();
  for(int i=1; i<=hDen->GetNbinsX(); i++){
      pull = 1.0; pullerr = 0.0;
      pullinv = 1.0; pullinverr = 0.0;
      if(hNum1->GetBinContent(i) > 0 && hDen->GetBinContent(i) > 0){
        pull    = (hNum1->GetBinContent(i)/hDen ->GetBinContent(i));
        pullinv = (hDen ->GetBinContent(i)/hNum1->GetBinContent(i));
	pullerr    = pull*   hNum1->GetBinError(i)/hNum1->GetBinContent(i);
	pullinverr = pullinv*hNum1->GetBinError(i)/hNum1->GetBinContent(i);
      }
      else {
        printf("0 events in %d\n",i);
      }
      if(isDebug) printf("ratio(%2d): pred/data = %.3f +/- %.3f predUnc: %.3f\n",i,pull,pullerr,hNum1->GetBinError(i)/hNum1->GetBinContent(i));
      if(isDebug) printf("ratio(%2d): data/pred = %.3f +/- %.3f, sigma = %.3f fb\n",i,pullinv,pullinverr,hNum1->GetBinContent(i));
      hRatio1->SetBinContent(i,pull);
      hRatio1->SetBinError(i,pullerr);
      hBand->SetBinContent(i,1);
      hBand->SetBinError  (i,hDen->GetBinError(i)/hDen->GetBinContent(i)); 
  }
  units = units.ReplaceAll("BIN","");
  atributes(hRatio1,XTitle.Data(),"#frac{Theory}{Data}",units.Data(),kRed+2);

  hRatio1->Draw("e");
  hBand->SetFillColor(12);
  hBand->SetFillStyle(3003);
  hBand->SetMarkerSize(0);
  hBand->SetLineWidth(0);
  hBand->Draw("E2same");
  
  // hPred2 w.r.t. hData
  TH1D* hRatio2  = (TH1D*) hPred2->Clone(); hRatio2->Reset();
  for(int i=1; i<=hDen->GetNbinsX(); i++){
      pull = 1.0; pullerr = 0.0;
      pullinv = 1.0; pullinverr = 0.0;
      if(hNum2->GetBinContent(i) > 0 && hDen->GetBinContent(i) > 0){
        pull = (hNum2->GetBinContent(i)/hDen->GetBinContent(i));
	pullerr = pull*hNum2->GetBinError(i)/hNum2->GetBinContent(i);
      }
      else {
        printf("0 events in %d\n",i);
      }
      if(isDebug) printf("ratio(%2d): pred/pred2 = %.3f +/- %.3f predUnc: %.3f\n",i,pull,pullerr,hNum2->GetBinError(i)/hNum2->GetBinContent(i));
      hRatio2->SetBinContent(i,pull);
      hRatio2->SetBinError(i,pullerr);
  }
  
  // hPred3 w.r.t. hData
  TH1D* hRatio3  = (TH1D*) hPred3->Clone(); hRatio3->Reset();
  for(int i=1; i<=hDen->GetNbinsX(); i++){
      pull = 1.0; pullerr = 0.0;
      pullinv = 1.0; pullinverr = 0.0;
      if(hNum3->GetBinContent(i) > 0 && hDen->GetBinContent(i) > 0){
        pull = (hNum3->GetBinContent(i)/hDen->GetBinContent(i));
	pullerr = pull*hNum3->GetBinError(i)/hNum3->GetBinContent(i);
      }
      else {
        printf("0 events in %d\n",i);
      }
      if(isDebug) printf("ratio(%2d): pred/Pred3 = %.3f +/- %.3f predUnc: %.3f\n",i,pull,pullerr,hNum3->GetBinError(i)/hNum3->GetBinContent(i));
      hRatio3->SetBinContent(i,pull);
      hRatio3->SetBinError(i,pullerr);
  }

  hRatio2->GetYaxis()->SetLabelFont  (   42);
  hRatio2->GetYaxis()->SetLabelOffset(0.015);
  hRatio2->GetYaxis()->SetLabelSize  (0.050);
  hRatio2->GetYaxis()->SetNdivisions (  505);
  hRatio2->GetYaxis()->SetTitleFont  (   42);
  hRatio2->GetYaxis()->SetTitleOffset(  1.2);
  hRatio2->GetYaxis()->SetTitleSize  (0.060);
  hRatio2->GetYaxis()->SetTickLength (0.03 );
  hRatio2->GetXaxis()->SetLabelFont  (   42);
  hRatio2->GetXaxis()->SetLabelSize  (0.040);
  hRatio2->GetXaxis()->SetNdivisions (  505);
  hRatio2->GetXaxis()->SetTitleFont  (   42);
  hRatio2->GetXaxis()->SetTitleSize  (0.060);
  hRatio2->GetXaxis()->SetTickLength (0.07 );

  hRatio3->GetYaxis()->SetLabelFont  (   42);
  hRatio3->GetYaxis()->SetLabelOffset(0.015);
  hRatio3->GetYaxis()->SetLabelSize  (0.050);
  hRatio3->GetYaxis()->SetNdivisions (  505);
  hRatio3->GetYaxis()->SetTitleFont  (   42);
  hRatio3->GetYaxis()->SetTitleOffset(  1.2);
  hRatio3->GetYaxis()->SetTitleSize  (0.060);
  hRatio3->GetYaxis()->SetTickLength (0.03 );
  hRatio3->GetXaxis()->SetLabelFont  (   42);
  hRatio3->GetXaxis()->SetLabelSize  (0.040);
  hRatio3->GetXaxis()->SetNdivisions (  505);
  hRatio3->GetXaxis()->SetTitleFont  (   42);
  hRatio3->GetXaxis()->SetTitleSize  (0.060);
  hRatio3->GetXaxis()->SetTickLength (0.07 );
 
  hRatio2->SetLineWidth(3);
  hRatio2->SetLineStyle(2);
  hRatio2->SetLineColor(kBlue);
  hRatio2->Draw("same,hist");
 
  hRatio3->SetLineWidth(3);
  hRatio3->SetLineStyle(3);
  hRatio3->SetLineColor(kMagenta);
  if(keyLabel0 != "EWKWZMJJ") hRatio3->Draw("same,hist");
  
  TLegend* leg = new TLegend(0.20,0.70,0.30,0.85);                                                    
  leg ->SetFillStyle(0);
  leg ->SetFillColor(kWhite);
  leg ->SetBorderSize(0);
  leg->SetTextSize(0.05);                                                                         
  leg->AddEntry(hRatio1,"Theoretical prediction 1","f");
  leg->AddEntry(hRatio2,"Theoretical prediction 2","l");
  leg->AddEntry(hRatio3,"Theoretical prediction 3","l");
  leg->AddEntry(hBand,"Experimental data","pe");
  //leg->Draw();
  // plotting again
  hRatio1->Draw("e,same");
  hRatio2->Draw("same,hist,e");

  // Draw a line throgh y=0
  double theLines[2] = {1.0, 0.5};
  TLine* baseline = new TLine(hRatio1->GetXaxis()->GetXmin(), theLines[0],
                              hRatio1->GetXaxis()->GetXmax(), theLines[0]);
  baseline->SetLineStyle(kDashed);
  baseline->Draw();
  // Set the y-axis range symmetric around y=0
  Double_t dy = TMath::Max(TMath::Abs(hRatio1->GetMaximum()),
                           TMath::Abs(hRatio1->GetMinimum())) + theLines[1];
  // Double_t dy = TMath::Max(TMath::Abs(TMath::Abs(hRatio1->GetMaximum())-1),TMath::Abs(TMath::Abs(hRatio1->GetMinimum()))-1);
  hRatio1->GetYaxis()->SetRangeUser(0.201,1.999);
  hRatio1->GetYaxis()->CenterTitle();
  eraselabel(pad1,hData->GetXaxis()->GetLabelSize());

  printf("Total yields: %f - %f -  %f - %f\n", hData->GetSumOfWeights(),hPred1->GetSumOfWeights(),hPred2->GetSumOfWeights(),hPred3->GetSumOfWeights());

  char CommandToExec[300];
  sprintf(CommandToExec,"mkdir -p plotsvbsvv");
  gSystem->Exec(CommandToExec);  

  TH1D* unfold = (TH1D*) hData->Clone("unfold");

  TString outputName = Form("unf_WW%s_normalized%d",keyLabel0.Data(),isNormalized);
  if(strcmp(outputName.Data(),"") != 0){
    TString myOutputFile;
    myOutputFile = Form("plotsvbsvv/%s.png",outputName.Data());
    c1->SaveAs(myOutputFile.Data());
    myOutputFile = Form("plotsvbsvv/%s.pdf",outputName.Data());
    c1->SaveAs(myOutputFile.Data());
    //myOutputFile = Form("plotsvbsvv/%s.root",outputName.Data());
    //TFile *outRoot = TFile::Open(myOutputFile.Data(),"recreate");
    //outRoot->cd();
    //unfold->Write();
    //outRoot->Close();
  }
}
