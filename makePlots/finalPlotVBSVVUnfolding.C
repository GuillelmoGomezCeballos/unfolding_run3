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
  const int theFillStyle2 = 3005;
  const int theFillColor3 = 46;
  const int theFillStyle3 = 3007;

  TString genFileName1 = "histogen_vbsvv_mapgraph.root";
  TFile *_fileGenWW1 = TFile::Open(genFileName1.Data());
  TH1D* hPred1     = (TH1D*)_fileGenWW1->Get(Form("hD%s",keyLabel0.Data()));     hPred1    ->SetDirectory(0);
  TH1D* hPred1_PDF = (TH1D*)_fileGenWW1->Get(Form("hD%s_PDF",keyLabel0.Data())); hPred1_PDF->SetDirectory(0);
  TH1D* hPred1_QCD = (TH1D*)_fileGenWW1->Get(Form("hD%s_QCD",keyLabel0.Data())); hPred1_QCD->SetDirectory(0);
  TH1D* hPred1_PS  = (TH1D*)_fileGenWW1->Get(Form("hD%s_PS",keyLabel0.Data()));  hPred1_PS ->SetDirectory(0);
  _fileGenWW1->Close();

  TString genFileName2 = "histogen_vbsvv_sherpa.root";
  TFile *_fileGenWW2 = TFile::Open(genFileName2.Data());
  TH1D* hPred2     = (TH1D*)_fileGenWW2->Get(Form("hD%s",keyLabel0.Data()));     hPred2    ->SetDirectory(0);
  TH1D* hPred2_PDF = (TH1D*)_fileGenWW2->Get(Form("hD%s_PDF",keyLabel0.Data())); hPred2_PDF->SetDirectory(0);
  TH1D* hPred2_QCD = (TH1D*)_fileGenWW2->Get(Form("hD%s_QCD",keyLabel0.Data())); hPred2_QCD->SetDirectory(0);
  TH1D* hPred2_PS  = (TH1D*)_fileGenWW2->Get(Form("hD%s_PS",keyLabel0.Data()));  hPred2_PS ->SetDirectory(0);
  _fileGenWW2->Close();

  TString plotName = Form("input_files/xs_%s_normalized%d.root", keyLabel0.Data(), isNormalized);
  TFile *_file0 = TFile::Open(plotName.Data());
  TH1D* hData = (TH1D*)_file0->Get(Form("hD%s",keyLabel0.Data())); hData->SetDirectory(0);
  _file0->Close();

  for(Int_t i=1;i<=hPred1->GetNbinsX();++i){
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

    hData->SetBinContent(i,hData->GetBinContent(i)*hPred1->GetBinContent(i));
    hData->SetBinError  (i,hData->GetBinError  (i)*hPred1->GetBinContent(i));

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
  }

  hData ->Scale(scaleDFSF);
  hPred1->Scale(scaleDFSF);
  hPred2->Scale(scaleDFSF);

  //hData ->Scale(1,"width");
  //hPred1->Scale(1,"width");
  //hPred2->Scale(1,"width");

  Int_t ww = 800;
  Int_t wh = 800;

  TCanvas *c1 = new TCanvas("c1", "c1", ww, wh);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.400, 1, 0.975);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.200, 1, 0.400);
  TPad* pad3 = new TPad("pad3", "pad3", 0, 0.000, 1, 0.200);

  pad1->SetTopMargin   (0.08);
  pad1->SetBottomMargin(0.00);

  pad2->SetTopMargin   (0.00);
  pad2->SetBottomMargin(0.00);

  pad3->SetTopMargin   (0.00);
  pad3->SetBottomMargin(0.30);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();

  pad1->cd();
  gStyle->SetOptStat(0);
  if(isLogY == true) pad1->SetLogy();
  if(isLogX == true) {pad1->SetLogx();}

  // draw pad1
  if(strcmp(units.Data(),"")==0){
    hPred1->GetXaxis()->SetTitle(XTitle.Data());
    hPred1->GetXaxis()->SetLabelOffset(0.005);
    hPred1->GetXaxis()->SetTitleOffset(  0.9);
  } 
  else {
    hPred1->GetXaxis()->SetTitle(Form("%s [%s]",XTitle.Data(),units.Data()));
    hPred1->GetXaxis()->SetLabelOffset(0.00);
    hPred1->GetXaxis()->SetTitleOffset(  1.1);
  }

  hPred1->GetYaxis()->SetTitle(theYTitle.Data());
  hPred1->GetYaxis()->SetLabelFont  (   42);
  hPred1->GetYaxis()->SetLabelOffset(0.015);
  hPred1->GetYaxis()->SetLabelSize  (0.060);
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
  hPred2->GetYaxis()->SetLabelSize  (0.060);
  hPred2->GetYaxis()->SetNdivisions (  505);
  hPred2->GetYaxis()->SetTitleFont  (   42);
  hPred2->GetYaxis()->SetTitleOffset(  1.0);
  hPred2->GetYaxis()->SetTitleSize  (0.080);
  hPred2->GetYaxis()->SetTickLength (0.03 );
  hPred2->GetXaxis()->SetLabelFont  (   42);
  hPred2->GetXaxis()->SetLabelSize  (0.040);
  hPred2->GetXaxis()->SetNdivisions (  505);
  hPred2->GetXaxis()->SetTitleFont  (   42);
  hPred2->GetXaxis()->SetTitleSize  (0.060);
  hPred2->GetXaxis()->SetTickLength (0.07 );
 
  hData->SetMarkerSize(0.8);
  hData->SetMarkerStyle(kFullCircle);
  hData->SetLineColor  (kBlack);

  hPred1->SetLineColor(kBlack);
  hPred1->SetMarkerStyle(3);
  hPred1->SetMarkerColor(kBlack);

  hPred2->SetLineColor(kBlue);
  hPred2->SetMarkerStyle(2);
  hPred2->SetMarkerColor(kBlue);

  TAxis *xa = hData->GetXaxis(); TAxis *xb = hPred1->GetXaxis();
  hPred1->SetTitle("");
  hPred2->SetTitle("");
  hData ->SetTitle("");
  double normalization[3] = {1.0, 1.0, 1.0};
  if(isNormalized) {
    normalization[0] = hPred1->GetSumOfWeights();
    normalization[1] = hPred2->GetSumOfWeights();
    normalization[2] = hData ->GetSumOfWeights();
  }
  hPred1->Scale(1./normalization[0]);
  hPred2->Scale(1./normalization[1]);
  hData ->Scale(1./normalization[2]);

  for(Int_t i=1;i<=hPred1->GetNbinsX();++i){
    if(isDebug) printf("hData (%2d) %7.3f %7.3f\n",i,hData->GetBinContent(i),hData->GetBinError(i));
  }

  if(isLogY == true) hPred1->GetYaxis()->SetRangeUser(hPred1->GetMinimum()/10,hPred1->GetMaximum()*100);
  else               hPred1->GetYaxis()->SetRangeUser(0.0,hPred1->GetMaximum()*1.4);
  hPred1->Draw("hist,x0");
  hPred2->Draw("hist,x0,same");
  hData->Draw("epx0,same");

  bool plotSystErrorBars = true;
  TGraphAsymmErrors * gsyst1 = new TGraphAsymmErrors(hPred1);
  if(plotSystErrorBars == true) {
    for (int i = 0; i < gsyst1->GetN(); ++i) {
      double systBck = 0;
      gsyst1->SetPointEYlow (i, sqrt(hPred1->GetBinError(i+1)*hPred1->GetBinError(i+1)+hPred1->GetBinContent(i+1)*hPred1->GetBinContent(i+1)*systBck*systBck));
      gsyst1->SetPointEYhigh(i, sqrt(hPred1->GetBinError(i+1)*hPred1->GetBinError(i+1)+hPred1->GetBinContent(i+1)*hPred1->GetBinContent(i+1)*systBck*systBck));
    }
    gsyst1->SetFillColor(theFillColor1);
    gsyst1->SetFillStyle(theFillStyle1);
    gsyst1->SetMarkerSize(0);
    gsyst1->SetLineWidth(0);
    gsyst1->SetLineColor(kWhite);
    gsyst1->Draw("E2same");
    //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
    //setex1->Draw();
  }
  TGraphAsymmErrors * gsyst2 = new TGraphAsymmErrors(hPred2);
  if(plotSystErrorBars == true) {
    for (int i = 0; i < gsyst2->GetN(); ++i) {
      double systBck = 0;
      gsyst2->SetPointEYlow (i, sqrt(hPred2->GetBinError(i+1)*hPred2->GetBinError(i+1)+hPred2->GetBinContent(i+1)*hPred2->GetBinContent(i+1)*systBck*systBck));
      gsyst2->SetPointEYhigh(i, sqrt(hPred2->GetBinError(i+1)*hPred2->GetBinError(i+1)+hPred2->GetBinContent(i+1)*hPred2->GetBinContent(i+1)*systBck*systBck));
    }
    gsyst2->SetFillColor(theFillColor2);
    gsyst2->SetFillStyle(theFillStyle2);
    gsyst2->SetMarkerSize(0);
    gsyst2->SetLineWidth(0);
    gsyst2->SetLineColor(kWhite);
    gsyst2->Draw("E2same");
    //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
    //setex1->Draw();
  }
  hPred1->Draw("hist,x0,same");
  hPred2->Draw("hist,x0,same");
  hData ->Draw("epx0,same");

  TH1D* hColorDummy1 = (TH1D*) hPred1->Clone();
  hColorDummy1->SetFillColor(theFillColor1);
  hColorDummy1->SetFillStyle(theFillStyle1);

  TH1D* hColorDummy2 = (TH1D*) hPred2->Clone();
  hColorDummy2->SetFillColor(theFillColor2);
  hColorDummy2->SetFillStyle(theFillStyle2);

  gStyle->SetOptStat(0);
  TLegend* legend = new TLegend(0.58,0.70,0.78,0.90);
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   62);
  legend->SetTextSize  (0.040);
  legend->AddEntry(hData,  "Data", "ep");
  legend->AddEntry(hColorDummy1, "MADGRAPH+PYTHIA", "lf");
  legend->AddEntry(hColorDummy2, "SHERPA", "lf");
  legend->Draw();

  CMS_lumi( pad1, 2027, 11 );

  { // Pred1
  pad2->cd();
  gStyle->SetOptStat(0);

  TH1D* hNum = (TH1D*) hPred1->Clone(); hNum->Reset();
  TH1D* hDen = (TH1D*) hData ->Clone(); hDen->Reset();

  TH1D* hRatio = (TH1D*) hPred1->Clone(); hRatio->Reset();
  TH1D* hBand  = (TH1D*) hData ->Clone(); hBand ->Reset();

  hNum->Add(hPred1);
  hDen->Add(hData);

  double pull; 
  double pullerr;
  double pullinv; 
  double pullinverr;

  for(int i=1; i<=hNum->GetNbinsX(); i++){
      pull = 1.0; pullerr = 0.0;
      pullinv = 1.0; pullinverr = 0.0;
      if(hNum->GetBinContent(i) > 0 && hDen->GetBinContent(i) > 0){
        pull = (hNum->GetBinContent(i)/hDen->GetBinContent(i));
        pullerr = hDen->GetBinError(i)/hDen->GetBinContent(i);
        pullinv = (hDen->GetBinContent(i)/hNum->GetBinContent(i));
        pullinverr = pullinv*hDen->GetBinError(i)/hDen->GetBinContent(i);
      }
      else {
        printf("0 events in %d\n",i);
      }
      if(isDebug) printf("ratio(%2d): pred/data = %.3f +/- %.3f predUnc: %.3f\n",i,pull,pullerr,hNum->GetBinError(i)/hNum->GetBinContent(i));
      if(isDebug) printf("ratio(%2d): data/pred = %.3f +/- %.3f, sigma = %.3f pb\n",i,pullinv,pullinverr,hNum->GetBinContent(i));
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,pullerr);
      hBand->SetBinContent(i,1);
      hBand->SetBinError  (i,hNum->GetBinError(i)/hNum->GetBinContent(i)); 
  }
  units = units.ReplaceAll("BIN","");
  atributes(hRatio,XTitle.Data(),"#frac{MADGRAPH}{Data}",units.Data());

  hRatio->Draw("ex0");
  hBand->SetFillColor(theFillColor1);
  hBand->SetFillStyle(theFillStyle1);
  hBand->SetMarkerSize(0);
  hBand->SetLineWidth(0);
  hBand->Draw("E2same");

  if(1){
  TLegend* leg = new TLegend(0.20,0.80,0.30,0.90);                                                    
  leg->SetBorderSize(	 0);
  leg->SetFillColor (	 0);
  leg->SetTextAlign (	12);
  leg->SetTextFont  (	62);
  leg->SetTextSize  (0.065);
  leg->AddEntry(hBand,"Theoretical uncertainty","f");
  leg->AddEntry(hRatio,"Theoretical prediction / measurement","pe");
  leg->Draw();
  }

  // Draw a line throgh y=0
  double theLines[2] = {1.0, 0.5};
  TLine* baseline = new TLine(hRatio->GetXaxis()->GetXmin(), theLines[0],
                              hRatio->GetXaxis()->GetXmax(), theLines[0]);
  baseline->SetLineStyle(kDashed);
  baseline->Draw();
  hRatio->Draw("ex0,same");
  // Set the y-axis range symmetric around y=0
  Double_t dy = TMath::Max(TMath::Abs(hRatio->GetMaximum()),
                           TMath::Abs(hRatio->GetMinimum())) + theLines[1];
  // Double_t dy = TMath::Max(TMath::Abs(TMath::Abs(hRatio->GetMaximum())-1),TMath::Abs(TMath::Abs(hRatio->GetMinimum()))-1);
  double maxValue = 1.499;
  double minValue = 0.601;
  if     (keyLabel0.Contains("WZ")) {maxValue = 1.850; minValue = 0.150;}
  else if(keyLabel0.Contains("NJET")) {maxValue = 1.999; minValue = 0.601;}
  hRatio->GetYaxis()->SetRangeUser(minValue,maxValue);
  hRatio->GetYaxis()->CenterTitle();
  eraselabel(pad1,hData->GetXaxis()->GetLabelSize());
  eraselabel(pad2,hRatio->GetXaxis()->GetLabelSize());
  } // End Pred1

  { // Pred2
  pad3->cd();
  gStyle->SetOptStat(0);

  TH1D* hNum = (TH1D*) hPred2->Clone(); hNum->Reset();
  TH1D* hDen = (TH1D*) hData ->Clone(); hDen->Reset();

  TH1D* hRatio = (TH1D*) hPred2->Clone(); hRatio->Reset();
  TH1D* hBand  = (TH1D*) hData ->Clone(); hBand ->Reset();

  hNum->Add(hPred2);
  hDen->Add(hData);

  double pull; 
  double pullerr;
  double pullinv; 
  double pullinverr;

  for(int i=1; i<=hNum->GetNbinsX(); i++){
      pull = 1.0; pullerr = 0.0;
      pullinv = 1.0; pullinverr = 0.0;
      if(hNum->GetBinContent(i) > 0 && hDen->GetBinContent(i) > 0){
        pull = (hNum->GetBinContent(i)/hDen->GetBinContent(i));
        pullerr = hDen->GetBinError(i)/hDen->GetBinContent(i);
        pullinv = (hDen->GetBinContent(i)/hNum->GetBinContent(i));
        pullinverr = pullinv*hDen->GetBinError(i)/hDen->GetBinContent(i);
      }
      else {
        printf("0 events in %d\n",i);
      }
      if(isDebug) printf("ratio(%2d): pred/data = %.3f +/- %.3f predUnc: %.3f\n",i,pull,pullerr,hNum->GetBinError(i)/hNum->GetBinContent(i));
      if(isDebug) printf("ratio(%2d): data/pred = %.3f +/- %.3f, sigma = %.3f pb\n",i,pullinv,pullinverr,hNum->GetBinContent(i));
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,pullerr);
      hBand->SetBinContent(i,1);
      hBand->SetBinError  (i,hNum->GetBinError(i)/hNum->GetBinContent(i)); 
  }
  units = units.ReplaceAll("BIN","");
  atributes(hRatio,XTitle.Data(),"#frac{SHERPA}{Data}",units.Data());

  hRatio->Draw("ex0");
  hBand->SetFillColor(theFillColor1);
  hBand->SetFillStyle(theFillStyle1);
  hBand->SetMarkerSize(0);
  hBand->SetLineWidth(0);
  hBand->Draw("E2same");

  if(1){
  TLegend* leg = new TLegend(0.20,0.80,0.30,0.90);                                                    
  leg->SetBorderSize(	 0);
  leg->SetFillColor (	 0);
  leg->SetTextAlign (	12);
  leg->SetTextFont  (	62);
  leg->SetTextSize  (0.065);
  leg->AddEntry(hBand,"Theoretical uncertainty","f");
  leg->AddEntry(hRatio,"Theoretical prediction / measurement","pe");
  leg->Draw();
  }

  // Draw a line throgh y=0
  double theLines[2] = {1.0, 0.5};
  TLine* baseline = new TLine(hRatio->GetXaxis()->GetXmin(), theLines[0],
                              hRatio->GetXaxis()->GetXmax(), theLines[0]);
  baseline->SetLineStyle(kDashed);
  baseline->Draw();
  hRatio->Draw("ex0,same");
  // Set the y-axis range symmetric around y=0
  Double_t dy = TMath::Max(TMath::Abs(hRatio->GetMaximum()),
                           TMath::Abs(hRatio->GetMinimum())) + theLines[1];
  // Double_t dy = TMath::Max(TMath::Abs(TMath::Abs(hRatio->GetMaximum())-1),TMath::Abs(TMath::Abs(hRatio->GetMinimum()))-1);
  double maxValue = 1.499;
  double minValue = 0.601;
  if     (keyLabel0.Contains("WZ")) {maxValue = 1.850; minValue = 0.150;}
  else if(keyLabel0.Contains("NJET")) {maxValue = 1.999; minValue = 0.601;}
  hRatio->GetYaxis()->SetRangeUser(minValue,maxValue);
  hRatio->GetYaxis()->CenterTitle();
  eraselabel(pad1,hData->GetXaxis()->GetLabelSize());
  eraselabel(pad2,hRatio->GetXaxis()->GetLabelSize());
  } // End Pred2

  printf("Total yields: %f - %f - %f\n", hData->GetSumOfWeights(),hPred1->GetSumOfWeights(),hPred2->GetSumOfWeights());

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
