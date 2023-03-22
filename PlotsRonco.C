{
// ROOT includes
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TArrow.h>
#include <TLegend.h>
#include "rootStyle.h"
#include "TPave.h"
#include "TBox.h"
//standard header
#include "TObject.h"
#include <vector>
#include <string>
#include <cmath>
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom3.h"

//custom headers
#include "TEfficiency.h"


//TFile *pippo = new TFile("~/work/work_new/work_Gap/HistosGap2019/histoTime.root");

//TFile *pippo = new TFile("~/work/work_new/histogramsTime55Deg6T5Events.root");
//TFile *pippo2 = new TFile("~/work/work_new/histoTime.root");


TFile *miofileBefore = new TFile("before.root","read");
TH1D* hBefore;
TCanvas* c1Before = (TCanvas *)miofileBefore.Get("Eventsnumber");
hBefore = (TH1D )c1Before->GetPrimitive("h6t5__1");
TFile *miofileAfter = new TFile("after.root","read");
TH1D* hAfter;
TCanvas* c1After = (TCanvas *)miofileAfter.Get("Eventsnumber");
hAfter = (TH1D )c1After->GetPrimitive("h6t5__1");
//TH1F *gTriggerE     = (TH1F *) pippo->Get("h6T5");
//TH1F *gTriggerEKathy     = (TH1F *) pippo2->Get("hKathy");

//TFile *fileExpo = new TFile("~/work/work_new/work_Gap/HistosGap2019/exposure.root");
/*TGraph* gExpo;
TCanvas* cExpo = (TCanvas *)fileExpo.Get("c1");
gExpo0 = (TGraph* )cExpo->GetPrimitive("Graph0");
gExpo1 = (TGraph* )cExpo->GetPrimitive("Graph1");

TFile *fileEn = new TFile("~/work/work_new/work_Gap/HistosGap2019/histoenergy.root");
TH1D* hEn;
TCanvas*cEn = (TCanvas *)fileEn.Get("c1");
hEn = (TH1D* )cEn->GetPrimitive("hevent");
          
TFile *file1 = new TFile("~/work/work_new/work_Gap/HistosGap2019/h6T5Old.root");
TFile *file2 = new TFile("~/work/work_new/work_Gap/HistosGap2019/h6T5After.root");
TH1D *h6T5Before     = (TH1D *) file1->Get("h6t5__1");
TH1D *h6T5After     = (TH1D *) file2->Get("h6t5__1");
*/
  // Define the ratio plot
 //  h6T5Before->Divide(h6T5After);


	  gStyle->SetOptStat(0);
          TCanvas  *tProfCanvasZenith      = new TCanvas("TProf6T5FullEff","TProf6T5FullEff",900,700);
	  tProfCanvasZenith->Range(-8.82213e+15,-0.2114754,1.112959e+17,1.192623);
	  tProfCanvasZenith->SetFillColor(0);
	  tProfCanvasZenith->SetFillStyle(4000);
	  tProfCanvasZenith->SetBorderMode(0);
	  tProfCanvasZenith->SetBorderSize(2);
	  tProfCanvasZenith->SetLeftMargin(0.1510025);
	  tProfCanvasZenith->SetRightMargin(0.08834586);
	  tProfCanvasZenith->SetBottomMargin(0.150613);
	  tProfCanvasZenith->SetFrameFillStyle(4000);
	  tProfCanvasZenith->SetFrameBorderMode(0);
	  tProfCanvasZenith->SetFrameBorderMode(0);
          hAfter->GetXaxis()->SetTitle("Year ");
          hAfter->GetXaxis()->SetTitleSize(0.05);
          hAfter->GetYaxis()->SetTitleSize(0.05);
          hAfter->GetXaxis()->SetTickLength(0.02);
          hAfter->GetYaxis()->SetTickLength(0.02);
          hAfter->GetYaxis()->SetLabelOffset(0.008);
          hAfter->GetXaxis()->SetLabelOffset(0.008);
          hAfter->GetYaxis()->SetLabelSize(0.041);
          hAfter->GetXaxis()->SetNdivisions(10);
          hAfter->GetYaxis()->SetNdivisions(10);
          hAfter->GetXaxis()->SetLabelSize(0.041);
          hAfter->GetXaxis()->SetTitleOffset(1.3);
          hAfter->GetYaxis()->SetTitleOffset(1.2);
          hAfter->GetYaxis()->SetTitle("6T5 Events / day");
          hBefore->SetLineColor(kRed-3);
          hBefore->SetLineWidth(1);
          hAfter->SetLineColor(kBlue-3);
          hAfter->SetLineWidth(1);
          hBefore->GetXaxis()->SetNdivisions(10);
          hBefore->GetYaxis()->SetNdivisions(10);
	  tProfCanvasZenith->cd();
	  hBefore->SetFillColor(kRed-3);
	  hAfter->SetFillColor(kBlue-3);
          hAfter->Draw( "");
	  hBefore->Draw("SAME");
          TLegend *leg = new TLegend(0.7318296,0.2872154,0.9310777,0.408056,NULL,"brNDC");
          leg->SetLineColor(kWhite);
          leg->SetFillStyle(0);
          gStyle->SetStatStyle(0);
          gStyle->SetTitleStyle(0);
          gROOT->ForceStyle();
          gStyle->SetOptStat(0);
          leg->AddEntry(hBefore,"Previous events : 33530","f");
          leg->AddEntry(hAfter,"Recovered events: 133416","f");
          leg->Draw();
  
  /*      TCanvas  *tProfCanvasZenith      = new TCanvas("TProf6T5FullEff","TProf6T5FullEff",900,700);
	  tProfCanvasZenith->Range(-8.82213e+15,-0.2114754,1.112959e+17,1.192623);
	  tProfCanvasZenith->SetFillColor(0);
	  tProfCanvasZenith->SetFillStyle(4000);
	  tProfCanvasZenith->SetBorderMode(0);
	  tProfCanvasZenith->SetBorderSize(2);
	  tProfCanvasZenith->SetLeftMargin(0.1510025);
	  tProfCanvasZenith->SetRightMargin(0.08834586);
	  tProfCanvasZenith->SetBottomMargin(0.150613);
	  tProfCanvasZenith->SetFrameFillStyle(4000);
	  tProfCanvasZenith->SetFrameBorderMode(0);
	  tProfCanvasZenith->SetFrameBorderMode(0);
          hAfter->GetYaxis()->CenterTitle();
          hAfter->GetXaxis()->CenterTitle();
          hAfter->GetXaxis()->SetTitle("Year ");
          hAfter->GetXaxis()->SetTitleSize(0.05);
          hAfter->GetYaxis()->SetTitleSize(0.05);
          hAfter->GetXaxis()->SetRangeUser(1e+16,1e+17);
          hAfter->GetXaxis()->SetTickLength(0.02);
          hAfter->GetYaxis()->SetTickLength(0.02);
          hAfter->GetYaxis()->SetLabelOffset(0.008);
          hAfter->GetXaxis()->SetLabelOffset(0.008);
          hAfter->GetYaxis()->SetLabelSize(0.041);
          hAfter->GetXaxis()->SetNdivisions(10);
          hAfter->GetYaxis()->SetNdivisions(10);
          hAfter->GetXaxis()->SetLabelSize(0.041);
          hAfter->GetXaxis()->SetTitleOffset(1.3);
          hAfter->GetYaxis()->SetTitleOffset(1.2);
          hAfter->GetYaxis()->SetTitle("6T5 Events / day");
          hBefore->SetLineColor(kViolet-6);
          hBefore->SetLineWidth(1);
          hAfter->SetLineColor(kGreen+2);
          hAfter->SetLineWidth(1);
          gTriggerE->SetLineColor(kOrange-3);
          gTriggerE->SetLineWidth(1);
          hBefore->GetXaxis()->SetNdivisions(10);
          hBefore->GetYaxis()->SetNdivisions(10);
	  tProfCanvasZenith->cd();
	  hBefore->SetFillColor(kViolet-6);
	  //hBefore->SetFillStyle(3002);
	  hAfter->SetFillColor(kGreen+2);
	  hAfter->SetFillStyle(3018);
          hAfter->Draw( "");
          gTriggerE->Draw( "SAME");
	  hBefore->Draw("SAME");
//	  tProfCanvasZenith->Update();

          TLegend *leg = new TLegend(0.7318296,0.2872154,0.9310777,0.408056,NULL,"brNDC");
          leg->SetLineColor(kWhite);
          leg->SetFillStyle(0);
          gStyle->SetStatStyle(0);
          gStyle->SetTitleStyle(0);
          gROOT->ForceStyle();
          gStyle->SetOptStat(0);
          leg->AddEntry(hBefore," Old events ","f");
          leg->AddEntry(gTriggerE,"Recovered events method A","f");
          leg->AddEntry(hAfter,"Recovered events method B","f");
          leg->Draw();
*/

/*
          TCanvas  *tProfCanvas      = new TCanvas("prova","prova",900,700);
          gTriggerE->GetYaxis()->CenterTitle();
          gTriggerE->GetXaxis()->CenterTitle();
          gTriggerE->GetXaxis()->SetTitle("Year ");
          gTriggerE->GetXaxis()->SetTitleSize(0.05);
          gTriggerE->GetYaxis()->SetTitleSize(0.05);
          gTriggerE->GetXaxis()->SetRangeUser(1e+16,1e+17);
          gTriggerE->GetXaxis()->SetTickLength(0.02);
          gTriggerE->GetYaxis()->SetTickLength(0.02);
          gTriggerE->GetYaxis()->SetLabelOffset(0.008);
          gTriggerE->GetXaxis()->SetLabelOffset(0.008);
          gTriggerE->GetYaxis()->SetLabelSize(0.041);
          gTriggerE->GetXaxis()->SetLabelSize(0.041);
          gTriggerE->GetXaxis()->SetTitleOffset(1.3);
          gTriggerE->GetYaxis()->SetTitleOffset(1.2);
          gTriggerE->GetYaxis()->SetTitle("Events");
          gTriggerE->SetLineColor(kOrange-3);
          gTriggerE->SetMarkerColor(kOrange-3);
          gTriggerE->SetMarkerStyle(21);
          gTriggerE->SetMarkerSize(1.2);
          gTriggerE->SetLineWidth(1);
	  gTriggerE->Draw( ""); 
          gTriggerEKathy->SetLineColor(kBlue-3);
          gTriggerEKathy->SetMarkerColor(kBlue-3);
          gTriggerEKathy->SetMarkerStyle(21);
          gTriggerEKathy->SetMarkerSize(1.2);
          gTriggerEKathy->SetLineWidth(1);
	  gTriggerEKathy->Draw( "SAME"); 
  //      gTriggerERoncoBefore->Draw( "SAME"); 
          TLegend *leg = new TLegend(0.7318296,0.2872154,0.9310777,0.408056,NULL,"brNDC"); 
          leg->SetLineColor(kWhite); 
          leg->SetFillStyle(0); 
          gStyle->SetStatStyle(0); 
          gStyle->SetTitleStyle(0); 
          gROOT->ForceStyle(); 
          gStyle->SetOptStat(0); 
          leg->AddEntry(gTriggerE,"6T5 events 2^{nd} method","l"); 
          leg->Draw(); 

	  TCanvas  *tEn      = new TCanvas("Energy","Energy",900,700);
	  tEn->Range(-8.82213e+15,-0.2114754,1.112959e+17,1.192623);
	  tEn->SetFillColor(0);
	  tEn->SetFillStyle(4000);
	  tEn->SetBorderMode(0);
	  tEn->SetBorderSize(2);
	  tEn->SetLeftMargin(0.1510025);
	  tEn->SetRightMargin(0.08834586);
	  tEn->SetBottomMargin(0.150613);
	  tEn->SetFrameFillStyle(4000);
	  tEn->SetFrameBorderMode(0);
	  tEn->SetFrameBorderMode(0);
          tEn->SetLogy();
          hEn->GetYaxis()->CenterTitle();
          hEn->GetXaxis()->CenterTitle();
          hEn->GetXaxis()->SetTitle("log_{10} E / eV ");
          hEn->GetXaxis()->SetTitleSize(0.05);
          hEn->GetYaxis()->SetTitleSize(0.05);
          hEn->GetXaxis()->SetRangeUser(1e+16,1e+19);
          hEn->GetXaxis()->SetTickLength(0.02);
          hEn->GetXaxis()->SetNdivisions(10);
          hEn->GetYaxis()->SetNdivisions(10);
	  hEn->GetYaxis()->SetTickLength(0.02);
          hEn->GetYaxis()->SetLabelOffset(0.008);
          hEn->GetXaxis()->SetLabelOffset(0.008);
          hEn->GetYaxis()->SetLabelSize(0.041);
          hEn->GetXaxis()->SetLabelSize(0.041);
          hEn->GetXaxis()->SetTitleOffset(1.3);
          hEn->GetYaxis()->SetTitleOffset(1.2);
          hEn->GetYaxis()->SetTitle("Events");
          hEn->SetLineColor(kGreen+2);
          hEn->SetMarkerColor(kGreen+2);
          hEn->SetMarkerSize(1.2);
          hEn->SetLineWidth(4);
	  hEn->Draw( ""); 
	  


	  TCanvas  *tExpo      = new TCanvas("Esposizione","Esposizione",900,700);
	  tExpo->Range(-8.82213e+15,-0.2114754,1.112959e+17,1.192623);
	  tExpo->SetFillColor(0);
	  tExpo->SetFillStyle(4000);
	  tExpo->SetBorderMode(0);
	  tExpo->SetBorderSize(2);
	  tExpo->SetLeftMargin(0.1510025);
	  tExpo->SetRightMargin(0.08834586);
	  tExpo->SetBottomMargin(0.150613);
	  tExpo->SetFrameFillStyle(4000);
	  tExpo->SetFrameBorderMode(0);
	  tExpo->SetFrameBorderMode(0);
          gExpo0->GetYaxis()->CenterTitle();
          gExpo0->GetXaxis()->CenterTitle();
          gExpo0->GetXaxis()->SetTitle("Year ");
          gExpo0->GetXaxis()->SetTitleSize(0.05);
          gExpo0->GetYaxis()->SetTitleSize(0.05);
          gExpo0->GetXaxis()->SetRangeUser(1e+16,1e+17);
          gExpo0->GetXaxis()->SetTickLength(0.02);
          gExpo0->GetXaxis()->SetNdivisions(10);
          gExpo0->GetYaxis()->SetNdivisions(10);
	  gExpo0->GetYaxis()->SetTickLength(0.02);
          gExpo0->GetYaxis()->SetLabelOffset(0.008);
          gExpo0->GetXaxis()->SetLabelOffset(0.008);
          gExpo0->GetYaxis()->SetLabelSize(0.041);
          gExpo0->GetXaxis()->SetLabelSize(0.041);
          gExpo0->GetXaxis()->SetTitleOffset(1.3);
          gExpo0->GetYaxis()->SetTitleOffset(1.2);
          gExpo0->GetYaxis()->SetTitle("Exposure / km^{2} sr yr");
          gExpo1->SetLineColor(kBlack);
          gExpo1->SetMarkerColor(kBlack);
          gExpo0->SetLineColor(kGreen+2);
          gExpo0->SetMarkerColor(kGreen+2);
          gExpo0->SetMarkerSize(1.2);
          gExpo0->SetLineWidth(4);
	  gExpo0->Draw( ""); 
          gExpo1->SetMarkerSize(1.2);
          gExpo1->SetLineWidth(4);
	  gExpo1->Draw( "SAME"); 
          TLegend *leg = new TLegend(0.1603563,0.7451565,0.3596882,0.828614,NULL,"brNDC");
	  leg->SetLineColor(kWhite); 
          leg->SetFillStyle(0); 
          gStyle->SetStatStyle(0); 
          gStyle->SetTitleStyle(0); 
          gROOT->ForceStyle(); 
          gStyle->SetOptStat(0); 
          leg->AddEntry(gExpo0,"#theta/#circ < 45","l"); 
          leg->AddEntry(gExpo1,"#theta/#circ < 55","l"); 
          leg->Draw(); 
          
*/}
