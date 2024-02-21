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

void makeWWResult(TString type, TString the0J = "", bool isNormalized = false){

  const int nBinWWMLL = 13;
  Float_t xbinsWWMLL[nBinWWMLL+1];
    xbinsWWMLL[ 0] =  55;      xbinsWWMLL[ 1] =  75;      xbinsWWMLL[ 2] =  85;      xbinsWWMLL[ 3] =   95;     xbinsWWMLL[ 4] = 110;      
    xbinsWWMLL[ 5] = 125;      xbinsWWMLL[ 6] = 140;      xbinsWWMLL[ 7] = 160;      xbinsWWMLL[ 8] =  185;     xbinsWWMLL[ 9] = 220;      
    xbinsWWMLL[10] = 280;      xbinsWWMLL[11] = 380;      xbinsWWMLL[12] = 600;      xbinsWWMLL[13] = 1500;

  const int nBinWWDPHILL = 9;
  Float_t xbinsWWDPHILL[nBinWWDPHILL+1];
    xbinsWWDPHILL[ 0] =   0*TMath::Pi()/180.;      xbinsWWDPHILL[ 1] =  20*TMath::Pi()/180.;      xbinsWWDPHILL[ 2] =  40*TMath::Pi()/180.;      xbinsWWDPHILL[ 3] =  60*TMath::Pi()/180.;      xbinsWWDPHILL[ 4] =  80*TMath::Pi()/180.;
    xbinsWWDPHILL[ 5] = 100*TMath::Pi()/180.;      xbinsWWDPHILL[ 6] = 120*TMath::Pi()/180.;      xbinsWWDPHILL[ 7] = 140*TMath::Pi()/180.;      xbinsWWDPHILL[ 8] = 160*TMath::Pi()/180.;      xbinsWWDPHILL[ 9] = 180*TMath::Pi()/180.;

  const int nBinWWPTL1 = 14;
  Float_t xbinsWWPTL1[nBinWWPTL1+1];
    xbinsWWPTL1[ 0] =  27;      xbinsWWPTL1[ 1] =  40;      xbinsWWPTL1[ 2] =  50;      xbinsWWPTL1[ 3] =  60;      xbinsWWPTL1[ 4] =  70;
    xbinsWWPTL1[ 5] =  80;      xbinsWWPTL1[ 6] =  90;      xbinsWWPTL1[ 7] = 100;      xbinsWWPTL1[ 8] = 110;      xbinsWWPTL1[ 9] = 130;      
    xbinsWWPTL1[10] = 150;      xbinsWWPTL1[11] = 175;      xbinsWWPTL1[12] = 220;      xbinsWWPTL1[13] = 300;      xbinsWWPTL1[14] = 400;

  const int nBinWWPTL2 = 8;
  Float_t xbinsWWPTL2[nBinWWPTL2+1];
    xbinsWWPTL2[ 0] =  25;      xbinsWWPTL2[ 1] =  30;      xbinsWWPTL2[ 2] =  35;      xbinsWWPTL2[ 3] =  40;      xbinsWWPTL2[ 4] =  45;      
    xbinsWWPTL2[ 5] =  50;      xbinsWWPTL2[ 6] =  75;      xbinsWWPTL2[ 7] = 100;      xbinsWWPTL2[ 8] = 150;

  const int nBinWWPTLL = 15;
  Float_t xbinsWWPTLL[nBinWWPTLL+1];
    xbinsWWPTLL[ 0] =  30;      xbinsWWPTLL[ 1] =  35;      xbinsWWPTLL[ 2] =  40;      xbinsWWPTLL[ 3] =  45;      xbinsWWPTLL[ 4] =  50;
    xbinsWWPTLL[ 5] =  55;      xbinsWWPTLL[ 6] =  60;      xbinsWWPTLL[ 7] =  65;      xbinsWWPTLL[ 8] =  70;      xbinsWWPTLL[ 9] =  75;
    xbinsWWPTLL[10] =  80;      xbinsWWPTLL[11] =  90;      xbinsWWPTLL[12] = 105;      xbinsWWPTLL[13] = 140;      xbinsWWPTLL[14] = 200;
    xbinsWWPTLL[15] = 300;

  const int nBinWWNJET = 3;
  Float_t xbinsWWNJET[nBinWWNJET+1];
    xbinsWWNJET[ 0] =-0.5;      xbinsWWNJET[ 1] = 0.5;      xbinsWWNJET[ 2] = 1.5;      xbinsWWNJET[ 3] = 2.5;

  const int nBinWWN0JET = 5;
  Float_t xbinsWWN0JET[nBinWWN0JET+1];
    xbinsWWN0JET[ 0] =-0.5;      xbinsWWN0JET[ 1] = 0.5;      xbinsWWN0JET[ 2] = 1.5;
    xbinsWWN0JET[ 3] = 2.5;      xbinsWWN0JET[ 4] = 3.5;      xbinsWWN0JET[ 5] = 4.5;

  TString xsfname("input_files/");
  TH1D* histoResult;
  if     (type == "MLL")    {histoResult = new TH1D(Form("hDWWMLL%s",the0J.Data()),    Form("hDWWMLL%s",the0J.Data()),    nBinWWMLL,    xbinsWWMLL   ); xsfname = xsfname + Form("WWMLL%s",the0J.Data());}
  else if(type == "DPHILL") {histoResult = new TH1D(Form("hDWWDPHILL%s",the0J.Data()), Form("hDWWDPHILL%s",the0J.Data()), nBinWWDPHILL, xbinsWWDPHILL); xsfname = xsfname + Form("WWDPHILL%s",the0J.Data());}
  else if(type == "PTL1")   {histoResult = new TH1D(Form("hDWWPTL1%s",the0J.Data()),   Form("hDWWPTL1%s",the0J.Data()),   nBinWWPTL1,   xbinsWWPTL1  ); xsfname = xsfname + Form("WWPTL1%s",the0J.Data());}
  else if(type == "PTL2")   {histoResult = new TH1D(Form("hDWWPTL2%s",the0J.Data()),   Form("hDWWPTL2%s",the0J.Data()),   nBinWWPTL2,   xbinsWWPTL2  ); xsfname = xsfname + Form("WWPTL2%s",the0J.Data());}
  else if(type == "PTLL")   {histoResult = new TH1D(Form("hDWWPTLL%s",the0J.Data()),   Form("hDWWPTLL%s",the0J.Data()),   nBinWWPTLL,   xbinsWWPTLL  ); xsfname = xsfname + Form("WWPTLL%s",the0J.Data());}
  else if(type == "NJET")   {histoResult = new TH1D(Form("hDWWNJET"),                  Form("hDWWNJET"),                  nBinWWNJET,   xbinsWWNJET  ); xsfname = xsfname + "WWNJET";}
  else if(type == "NJETS")  {histoResult = new TH1D(Form("hDWWNJETS"),                 Form("hDWWNJETS"),                 nBinWWNJET,   xbinsWWNJET  ); xsfname = xsfname + "WWNJETS";}
  else if(type == "N0JET")  {histoResult = new TH1D(Form("hDWWN0JET"),                 Form("hDWWN0JET"),                 nBinWWN0JET,  xbinsWWN0JET ); xsfname = xsfname + "WWN0JET";}
  else {printf("WRONG TYPE\n"); return;}

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

  TString outNtuplename = Form("input_files/xs_WW%s%s_normalized%d.root",type.Data(),the0J.Data(),isNormalized);
  TFile *outtuple = TFile::Open(outNtuplename.Data(),"recreate");
  outtuple->cd();
  histoResult->Write();
  outtuple->Close();
}

void makeWWNtuple(){
 //makeWWResult("MLL", ""    , false);
 //makeWWResult("MLL", "0JET", false);
 //makeWWResult("MLL", ""    , true);
 //makeWWResult("MLL", "0JET", true);

 //makeWWResult("DPHILL", ""    , false);
 //makeWWResult("DPHILL", "0JET", false);
 //makeWWResult("DPHILL", ""    , true);
 //makeWWResult("DPHILL", "0JET", true);

 //makeWWResult("PTL1", ""    , false);
 //makeWWResult("PTL1", "0JET", false);
 //makeWWResult("PTL1", ""    , true);
 //makeWWResult("PTL1", "0JET", true);

 //makeWWResult("PTL2", ""    , false);
 //makeWWResult("PTL2", "0JET", false);
 //makeWWResult("PTL2", ""    , true);
 //makeWWResult("PTL2", "0JET", true);

 //makeWWResult("PTLL", ""    , false);
 //makeWWResult("PTLL", "0JET", false);
 //makeWWResult("PTLL", ""    , true);
 //makeWWResult("PTLL", "0JET", true);

 //makeWWResult("NJET", ""    , false);
 //makeWWResult("NJET", ""    , true);
 makeWWResult("NJETS", ""   , false);
 makeWWResult("NJETS", ""   , true);
 //makeWWResult("N0JET", ""   , false);

}
