#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>                  // Math
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TFile.h>                  // File handle
#include <TH1D.h>                   // 1D histogram class
#endif

void makeVBSVVResult(TString type, bool isNormalized = false){

  const int nBin4 = 4; const int nBin2 = 2;
  const Float_t xbinsWZMJJ       [nBin4+1] = {500, 900, 1300, 1900, 2500};
  const Float_t xbinsWWMJJ       [nBin4+1] = {500, 900, 1300, 1900, 2500};
  const Float_t xbinsWWMLL       [nBin4+1] = {20, 80, 140, 220, 300};
  const Float_t xbinsWWNJET      [nBin2+1] = {1.5, 2.5, 3.5};
  const Float_t xbinsWWDELTAETAJJ[nBin4+1] = {2.5,3.6,4.5,5.5,8.0};
  const Float_t xbinsWWDELTAPHIJJ[nBin4+1] = {0.0, 1.8, 2.5, 2.9, TMath::Pi()};

  TString xsfname(Form("input_files/%s",type.Data()));
  TH1D* histoResult;
  if     (type == "EWKWZMJJ")        {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWZMJJ       );}
  else if(type == "EWKWWMJJ")        {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWWMJJ       );}
  else if(type == "EWKWWMLL")        {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWWMLL       );}
  else if(type == "EWKWWNJET")       {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin2, xbinsWWNJET      );}
  else if(type == "EWKWWDELTAETAJJ") {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWWDELTAETAJJ);}
  else if(type == "EWKWWDELTAPHIJJ") {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWWDELTAPHIJJ);}
  else if(type == "QCDWZMJJ")        {histoResult = new TH1D(Form("hD%s",type.Data()), Form("hD%s",type.Data()), nBin4, xbinsWZMJJ       );}
  else {printf("WRONG TYPE\n"); return;}

  printf("Making %s ...\n",type.Data());
  if(isNormalized) xsfname = xsfname + "_normalized.txt";
  else             xsfname = xsfname + ".txt";

  int count = 0;
  ifstream ifs;
  ifs.open(xsfname.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    Double_t r,rup,rdown;
    stringstream ss(line);
    ss >> r >> rup >> rdown;
    count++;
    histoResult->SetBinContent(count, r);
    histoResult->SetBinError  (count, (rup+rdown)/2.0);
  }
  ifs.close();

  if(histoResult->GetNbinsX() != count) {printf("DIFFERENT NUMBER OF BINS IN HISTOGRAM AND INPUT FILE: %d %d\n",histoResult->GetNbinsX(),count); return;}

  TString outNtuplename = Form("input_files/xs_%s_normalized%d.root",type.Data(),isNormalized);
  TFile *outtuple = TFile::Open(outNtuplename.Data(),"recreate");
  outtuple->cd();
  histoResult->Write();
  outtuple->Close();
}

void makeVBSVVNtuple(){
 makeVBSVVResult("EWKWZMJJ"        , false);
 makeVBSVVResult("EWKWWMJJ"        , false);
 makeVBSVVResult("EWKWWMLL"        , false);
 makeVBSVVResult("EWKWWNJET"       , false);
 makeVBSVVResult("EWKWWDELTAETAJJ" , false);
 makeVBSVVResult("EWKWWDELTAPHIJJ" , false);
}
