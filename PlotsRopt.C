#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"
// ROOT includes
 #include <TGraph2DErrors.h>
 #include <TGraph2D.h>
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
 #include "TGraphErrors.h"
 #include "TH1.h"
 #include "TList.h"
 #include "TMath.h"
 #include <math.h>
 #include "TROOT.h"
 #include "TRandom3.h"
 #include <cmath>
 #include <iostream>
 using namespace std;
 void PlotsRopt() {
	gStyle->SetOptStat(0);
//	gStyle->SetOptTitle(kFALSE);

//
//	TFile *inputFileldf = new TFile("histogramsPerDiegoFixedBetaAllyears.root");//nuovo bineado beta fisso 2013-2020(dicembre)
	TFile *inputFileldf = new TFile("histogramsPerDiegoFreeBeta.root");//nuovo bineado beta fisso 2013-2020(dicembre)
        TFile *histoPerDiego2 = new TFile("histogramsPerDiego2.root");//per Ropt a 250
        TFile *histoPerDiegoFixedBeta = new TFile("histogramsPerDiegoFixedBetaAllyears.root");//per Ropt a 250 beta fisso
        TFile *histoFixedBetaFullEff_simProton = new TFile("histogramsFixedBetaFullEff_simProton.root");//per Ropt a 250 beta fisso full efficiency
        TFile *histoFixedBetaFullEff_simIron = new TFile("histogramsFixedBetaFullEff_simIron.root");//per Ropt a 250 beta fisso full efficiency
        TFile *histoFixedBetaFullEff_simHadrons = new TFile("histogramsFixedBetaFullEff_simHadrons.root");//per Ropt a 250 beta fisso full efficiency
        TFile *histoPerDiegoFixedBetaFullEff = new TFile("histogramsPerDiegoFixedBetaFullEff.root");//per Ropt a 250 beta fisso full efficiency
        TFile *histoPerDiegoFreeBetaFullEff = new TFile("histogramsPerDiegoFreeBetaFullEff.root");//per Ropt a 250 beta libero full Efficiency
	TFile *inputFilemappa = new TFile("mappaRopt_infill4332019.root");//mappa con array  e eventi con ricosgtruiti con Tagli in Ropt
	//TFile *roptFile = new TFile("histogramsFixedBetaFullEffNewBin.root");//nuovo binnato per grafici ThRopt2D e hprofRopt
	TFile *roptFile = new TFile("histogramsLDFRoptInfill.root");//infill 2013-2019
	TFile *roptFileLogSTheta250 = new TFile("histogramsLDFRopt250NewBin.root");//FreeBeta logSThetaRopt
	TFile *roptFileLogSTheta300 = new TFile("histogramsLDFRopt300NewBin.root");// logSThetaRopt
        TFile *histoLDF300 = new TFile("histogramsLDF300.root");//per Ropt a 300 beta free full eff
        TFile *histoLDF250 = new TFile("histogramsLDF250.root");//per Ropt a 250 beta free full eff
        TFile *histoLDF350 = new TFile("histogramsLDF350.root");//per Ropt a 250 beta free full eff
        TFile *histoLDF200 = new TFile("histogramsLDF200.root");//per Ropt a 200 beta free full eff
	TCanvas  *canvasMappa = new TCanvas("mappa433EventiRecRoptCut","mappa RoptCut",963,800);
	
	
//mappa Eventi rec Con tagli in Ropt	
	TGraph *mappaMarkers = (TGraph *) inputFilemappa->Get("graphmarkers");
	TH2D *ThMappaEvt = (TH2D *)inputFilemappa->Get("ThRopt3");//ThRopt 0, 1, 2
	gStyle->SetPalette(54);
	gStyle->SetNumberContours(255);
	//ThMappaEvt->SetTitle(" R_{opt} / m < 260 "); 
//	ThMappaEvt->SetTitle("260 < R_{opt} / m < 280 "); 
//	ThMappaEvt->SetTitle("280 < R_{opt} / m < 300 "); 
	ThMappaEvt->SetTitle(" R_{opt} / m > 300 "); 
	ThMappaEvt->SetTitleSize(0.031); 
	ThMappaEvt->GetXaxis()->SetTitleOffset(1.1);
	ThMappaEvt->GetYaxis()->SetTitleOffset(1.1);
	ThMappaEvt->GetXaxis()->CenterTitle();
	ThMappaEvt->GetZaxis()->CenterTitle();
	ThMappaEvt->GetYaxis()->CenterTitle();
	ThMappaEvt->GetYaxis()->SetTitle("y / m ");
	ThMappaEvt->GetXaxis()->SetTitle("x / m");
	ThMappaEvt->GetZaxis()->SetTitle("N events ");
	mappaMarkers->SetMarkerColor(kGray+2);
	mappaMarkers->SetLineColor(kGray+3);
	mappaMarkers->SetMarkerStyle(8);
	mappaMarkers->SetMarkerSize(2);
	ThMappaEvt->GetYaxis()->SetTitleSize(0.031);
	ThMappaEvt->GetXaxis()->SetTitleSize(0.031);
	ThMappaEvt->GetZaxis()->SetTitleSize(0.031);
	ThMappaEvt->GetYaxis()->SetLabelSize(0.021);
	ThMappaEvt->GetXaxis()->SetLabelSize(0.021);
	ThMappaEvt->GetZaxis()->SetLabelSize(0.021);
	ThMappaEvt->GetXaxis()->SetTickLength(0.02);
	ThMappaEvt->GetYaxis()->SetTickLength(0.02);
	ThMappaEvt->GetXaxis()->SetNdivisions(5);
	ThMappaEvt->GetYaxis()->SetNdivisions(5);
	canvasMappa->cd();
	ThMappaEvt->Draw("colz");
	mappaMarkers->Draw("SAME P");

         //usare root6 per il invertpalette
	 TCanvas  *canvas2D          = new TCanvas("canvas2D","Th2D",900,700);
	 TH2D *ldf2D = (TH2D *) roptFile->Get("ThRopt2D");
	 //TColor::InvertPalette(); 
	 TProfile *prof = (TProfile *) roptFile->Get("hprofRopt"); 
	 ldf2D->GetXaxis()->SetTitleOffset(1.3);
	 ldf2D->GetYaxis()->SetTitleOffset(1.3);
	 ldf2D->GetXaxis()->CenterTitle();
	 ldf2D->GetZaxis()->CenterTitle();
	 ldf2D->GetYaxis()->CenterTitle();
	 ldf2D->GetYaxis()->SetTitle("r_{opt} / m ");
	 ldf2D->GetXaxis()->SetTitle("r_{S_{max}} / m");
	 ldf2D->GetZaxis()->SetTitle("N events ");
	 ldf2D->GetYaxis()->SetRangeUser(0,500);
	 ldf2D->GetXaxis()->SetRangeUser(50,260);
	 ldf2D->SetMarkerColor(kGray+2);
	 ldf2D->SetLineColor(kGray+2);
	 ldf2D->SetMarkerStyle(8);
	 ldf2D->SetMarkerSize(1);
	 prof->SetMarkerColor(kBlack);
	 prof->SetLineColor(kBlack);
	 prof->SetLineWidth(2);
	 prof->SetMarkerSize(1.2);
	 prof->SetMarkerStyle(24);
	 ldf2D->GetYaxis()->SetLabelSize(0.031);
	 ldf2D->GetXaxis()->SetLabelSize(0.031);
	 ldf2D->GetXaxis()->SetTickLength(0.02);
	 ldf2D->GetYaxis()->SetTickLength(0.02);
	 ldf2D->GetXaxis()->SetNdivisions(10);
	 ldf2D->GetYaxis()->SetNdivisions(10);
	 ldf2D->SetLineWidth(1);
	 ldf2D->Draw("colz");
	 prof->Draw("SAME P");	 
           TLine *line1 = new TLine(50,300,260,300);
	   line1->SetLineColor(kBlue+3);
           line1->SetLineStyle(1);
           line1->SetLineWidth(3);
           line1->Draw("SAME");
           TLine *line6 = new TLine(50,270,260,270);
	   line6->SetLineColor(kRed-3);
           line6->SetLineStyle(7);
           line6->SetLineWidth(3);
           line6->Draw("SAME");
  
        
	   
	   /*
	   
	TCanvas  *canvasRoptMult = new TCanvas("RoptMult","cRoptMult",900,700); //
	std::string name = "hRopt";
	std::string interval [] = { "< 5","5 - 10"," > 10" };
        TLegend *leg = new TLegend(0.7126949,0.6538462,0.9977728,0.8786982,NULL,"brNDC");
	leg->SetBorderSize(0);
        leg->SetLineColor(kWhite);
        leg->SetFillStyle(0);
        gStyle->SetStatStyle(0);
        gStyle->SetTitleStyle(0);
	int counter = 1;
        gROOT->ForceStyle();
	double nmax = 0;
	double n_Entries = 0;
	for(int bin=0; bin<=2; bin++) {
	   counter = counter+2;
	   stringstream ss;
	   ss << bin;
	   string str = ss.str();
	   std::string nameTh = name + str;
	   stringstream label;
	   stringstream mean;
	   label << interval[bin] << " fold";
           cout << nameTh << " = name TH1F" << endl;
	   TH1F *hRopt  = (TH1F *) roptFile->Get(nameTh.c_str());
           hRopt->SetTitle("r_{opt} / m infill");
           hRopt->GetXaxis()->SetTitleOffset(1.30);
           hRopt->GetYaxis()->SetTitleOffset(1.6);
           hRopt->GetXaxis()->CenterTitle();
           hRopt->GetYaxis()->CenterTitle();
           hRopt->GetYaxis()->SetRangeUser(0,22000);
           //hRopt->GetYaxis()->SetTitle("Events ");
           hRopt->GetYaxis()->SetTitle("Nomalized entries ");
           hRopt->GetXaxis()->SetTitle("Optimal distance (m)");
           hRopt->GetYaxis()->SetLabelOffset(0.014);
           hRopt->GetXaxis()->SetLabelOffset(0.014);
           hRopt->GetYaxis()->SetLabelSize(0.031);
           hRopt->GetXaxis()->SetLabelSize(0.031);
           hRopt->GetXaxis()->SetTickLength(-0.010);
           hRopt->GetYaxis()->SetTickLength(-0.010);
           hRopt->GetXaxis()->SetNdivisions(10);         
           hRopt->GetYaxis()->SetNdivisions(10);
           hRopt->SetLineColor(kAzure-bin-counter);
	   Double_t norm = 1;
           hRopt->Scale(norm/hRopt->Integral());//normalize by the “integral” to show the frequency probability in each bin
	   cout << " bin+counter=" << bin+counter << endl;
	   hRopt->SetLineWidth(2.0);
           //mean value
	   cout << " mean = " << hRopt->GetMean() << " Entries " << hRopt->GetEntries()  <<endl;
           double meanVal = hRopt->GetMean();
	   cout << "mean Val " << meanVal << endl;
//	   mean << "r_{opt} = " << int(hRopt->GetMean());
	   mean << int(hRopt->GetMean()) << " m "; 
	   //Draw
	   if (bin>0) hRopt->Draw("SAME");
	   else hRopt->Draw("");
	   n_Entries = hRopt->GetBinContent(hRopt->GetMaximumBin());
	   if (n_Entries>nmax) nmax = n_Entries;
	   //line (dentro for perchè ogni volta creo una linea nuova)
           TLine *line = new TLine(meanVal,0,meanVal,22000);
	   line->SetLineColor(kAzure-bin-counter);
           line->SetLineStyle(7-counter+bin);
           line->SetLineWidth(3);
           line->Draw("SAME");
	   //text (dentro for perchè ogni volta creo un text nuovo)
           TLatex *   tex = new TLatex(332,2583.207, mean.str().c_str());
     	   tex->SetTextAlign(13);
	   tex->SetTextColor(kAzure-bin-counter);
 	   tex->SetTextFont(42);
 	   tex->SetTextSize(0.03);
 	   tex->SetLineWidth(2);
//        tex->Draw("SAME");
	   //legend (fuori del for perchè la creo una volta sola e qui la riempio)
           leg->AddEntry(hRopt, label.str().c_str(), "f");//interval[bin].c_str(),"f");
           leg->AddEntry(line,  mean.str().c_str(), "l");//interval[bin].c_str(),"f");
           leg->Draw("SAME");
  	}

	   cout <<"maxBin content" << nmax << endl;

	   cout << "---------------------------------------------------------------" << endl;
*/

	   /*
	    if ((thetaDeg_433>0 &&  thetaDeg_433<=30) && (log10(S250) >= 1 && log10(S250)<=1.5)) myVectorRoptLogSTheta[0]->Fill(Ropt);
	    if ((thetaDeg_433>30 &&  thetaDeg_433<=45) && (log10(S250) >= 1 && log10(S250)<=1.5)) myVectorRoptLogSTheta[1]->Fill(Ropt);
	    if ((thetaDeg_433>0 &&  thetaDeg_433<=30) && (log10(S250) > 1.5 && log10(S250)<=2)) myVectorRoptLogSTheta[2]->Fill(Ropt);
	    if ((thetaDeg_433>30 &&  thetaDeg_433<=45) && (log10(S250) > 1.5 && log10(S250)<=2)) myVectorRoptLogSTheta[3]->Fill(Ropt);
	    */


        TCanvas  *canvasRoptLogSTheta = new TCanvas("RoptLogSTheta","cRoptlogSTheta",900,700); //
	short colour[6] = {kBlue-4, kRed-4, kBlue-9, kRed-9};
	std::string name = "hRoptlogSTheta";
	std::string intervalTheta [] = { "#theta/#circ #leq 30 ","30 #leq #theta/#circ #leq 45","#theta/#circ #leq 30 ","30 #leq #theta/#circ #leq 45"};
	std::string intervalLogS [] = { "1 #leq log_{10}S_{250} #leq 1.5", "1 #leqlog_{10}S_{250} #leq 1.5","1.5 #leq log_{10}S_{250} #leq 2", "1.5 #leq log_{10}S_{250} #leq 2" };
        TLegend *leg2 = new TLegend(0.7126949,0.6538462,0.9977728,0.8786982,NULL,"brNDC");
	leg2->SetBorderSize(0);
        leg2->SetLineColor(kWhite);
        leg2->SetFillStyle(0);
        gStyle->SetStatStyle(0);
        gStyle->SetTitleStyle(0);
	int counter = 0;
        gROOT->ForceStyle();
	double nmax = 0;
	double n_Entries = 0;
	for(int bin=0; bin<=3; bin++) {
	   cout << "bin = ---- " << bin << endl;
	   stringstream label2;
	   label2  << intervalTheta[bin]  << "   " <<  intervalLogS[bin];
	   counter = counter+2;
	   stringstream ss;
	   ss << bin;
	   string str = ss.str();
	   std::string nameTh = name + str;
	   stringstream mean;
           cout << nameTh << " = name TH1F" << endl;
	   TH1F *hRoptLogSTheta  = (TH1F *) roptFileLogSTheta250->Get(nameTh.c_str());
           //hRoptLogSTheta->SetTitle("r_{opt} / m ");
           hRoptLogSTheta->GetXaxis()->SetTitleOffset(1.30);
           hRoptLogSTheta->GetYaxis()->SetTitleOffset(1.6);
           hRoptLogSTheta->GetXaxis()->CenterTitle();
           hRoptLogSTheta->GetYaxis()->CenterTitle();
           hRoptLogSTheta->GetYaxis()->SetRangeUser(0,8000);
           hRoptLogSTheta->GetXaxis()->SetRangeUser(0,450);
           hRoptLogSTheta->GetYaxis()->SetTitle("Events ");
           //hRoptLogSTheta->GetYaxis()->SetTitle("Nomalized entries ");
           hRoptLogSTheta->GetXaxis()->SetTitle("r_{opt} / m");
	   hRoptLogSTheta->GetXaxis()->SetTitleSize(0.05);
           hRoptLogSTheta->GetYaxis()->SetTitleSize(0.050);
	   hRoptLogSTheta->GetYaxis()->SetLabelOffset(0.014);
           hRoptLogSTheta->GetXaxis()->SetLabelOffset(0.014);
           hRoptLogSTheta->GetYaxis()->SetLabelSize(0.045);
           hRoptLogSTheta->GetXaxis()->SetLabelSize(0.045);
           hRoptLogSTheta->GetXaxis()->SetTickLength(-0.010);
           hRoptLogSTheta->GetYaxis()->SetTickLength(-0.010);
           hRoptLogSTheta->GetXaxis()->SetNdivisions(10);         
           hRoptLogSTheta->GetYaxis()->SetNdivisions(10);
           hRoptLogSTheta->SetLineColor(colour[bin]);
           hRoptLogSTheta->SetLineWidth(3);
	   //Double_t norm = 1;
           //hRoptLogSTheta->Scale(norm/hRoptLogSTheta->Integral());//normalize by the “integral” to show the frequency probability in each bin
	   cout << " bin+counter=" << bin+counter << endl;
	   hRoptLogSTheta->SetLineWidth(3.0);
           //mean value
	   cout  << intervalTheta[bin]  << "   " <<  intervalLogS[bin];
	   cout << " mean = " << hRoptLogSTheta->GetMean() << " Entries " << hRoptLogSTheta->GetEntries()  <<endl;
           double meanVal = hRoptLogSTheta->GetMean();
	   cout << "mean Val " << meanVal << endl;
//	   mean << "r_{opt} = " << int(hRoptLogSTheta->GetMean());
	   mean << int(hRoptLogSTheta->GetMean()) << " m "; 
	   //Draw
	   if (bin>0) hRoptLogSTheta->Draw("SAME");
	   else hRoptLogSTheta->Draw("");
	   n_Entries = hRoptLogSTheta->GetBinContent(hRoptLogSTheta->GetMaximumBin());
	   if (n_Entries>nmax) nmax = n_Entries;
	   //line (dentro for perchè ogni volta creo una linea nuova)
	   TLine *line2 = new TLine(meanVal,0,meanVal,8000);
	   line2->SetLineColor(colour[bin]);
           line2->SetLineStyle(7);
           line2->SetLineWidth(3);
           line2->Draw("SAME");
	   TLatex *   tex9 = new TLatex(376.662,1723.66,"Real data");
	   tex9->SetTextAlign(13);
           tex9->SetTextColor(kBlack);
           tex9->SetTextFont(42);
           tex9->SetTextSize(0.05);
           tex9->SetLineWidth(2);
 	   tex9->Draw();
	   //text (dentro for perchè ogni volta creo un text nuovo)
           TLatex *   tex2 = new TLatex(332,2583.207, mean.str().c_str());
     	   tex2->SetTextAlign(13);
	   tex2->SetTextColor(colour[bin]);
 	   tex2->SetTextFont(42);
 	   tex2->SetTextSize(0.05);
 	   tex2->SetLineWidth(2);
	   tex2->Draw("SAME");
	   //legend (fuori del for perchè la creo una volta sola e qui la riempio)
           leg2->AddEntry(hRoptLogSTheta, label2.str().c_str(), "f");//interval[bin].c_str(),"f");
           //leg2->AddEntry(line2,  mean.str().c_str(), "l");//interval[bin].c_str(),"f");
           leg2->Draw("SAME");
  	
	}




/*    
 
        TCanvas  *canvasRoptPerDiego          = new TCanvas("canvasDiego","cRoptDiego",900,700); //hRoptSatur0 solo saturati, 1 solo no saturati, 2 tutti
        TH1F *hRoptFixed = (TH1F *) histoPerDiegoFixedBetaFullEff->Get("hRoptSatur1");//non saturati beta fisso
        TH1F *hRoptSatur = (TH1F *) histoPerDiegoFixedBetaFullEff->Get("hRoptSatur0"); //solo saturati beta fisso
        TH1F *hRoptFree = (TH1F *) histoPerDiegoFreeBetaFullEff->Get("hRoptSatur1");//non saturati beta libero
        TH1F *hRoptSaturFree = (TH1F *) histoPerDiego2->Get("hRoptSatur0"); //solo saturati
        hRoptFixed->GetXaxis()->SetTitleOffset(1.30);
        hRoptFixed->GetYaxis()->SetTitleOffset(1.6);
        hRoptFixed->GetXaxis()->CenterTitle();
        hRoptFixed->GetYaxis()->CenterTitle();
        hRoptFixed->GetYaxis()->SetTitle("Events ");
        hRoptFixed->GetXaxis()->SetTitle("Optimal distance (m)");
        hRoptFixed->GetYaxis()->SetLabelOffset(0.014);
        hRoptFixed->GetXaxis()->SetLabelOffset(0.014);
        hRoptFixed->GetYaxis()->SetLabelSize(0.031);
        hRoptFixed->GetXaxis()->SetLabelSize(0.031);
        hRoptFixed->GetXaxis()->SetTickLength(-0.010);
        hRoptFixed->GetYaxis()->SetTickLength(-0.010);
        hRoptFixed->GetXaxis()->SetNdivisions(10);         
        hRoptFixed->GetYaxis()->SetNdivisions(10);
        hRoptFixed->SetLineColor(kBlue-4);
        cout << "Fixed mean = " << hRoptFixed->GetMean() << " Entries " << hRoptFixed->GetEntries() << "Free mean = " << hRoptFree->GetMean() << " Entries " << hRoptFree->GetEntries() <<endl;
        hRoptFixed->SetLineWidth(2.0);
        hRoptFree->SetLineWidth(2.0);
        hRoptFree->SetLineColor(kOrange-4);
        hRoptFixed->Draw("");
        hRoptFree->Draw("SAME");
	gStyle->SetStatStyle(0);
	gStyle->SetTitleStyle(0);
	gROOT->ForceStyle();
	TLine *line1 = new TLine(296.,0,296.,23500);
	line1->SetLineColor(kBlue-4);
	line1->SetLineStyle(7);
	line1->SetLineWidth(3);
        line1->Draw("SAME");
        TLine *line2 = new TLine(289.,0,287.,23500);
        line2->SetLineColor(kOrange-4);
        line2->SetLineStyle(9);
        line2->SetLineWidth(3);
        line2->Draw("SAME");
        TLegend *leg2 = new TLegend(0.72049,0.7884615,0.9387528,0.8727811,NULL,"brNDC");
        leg2->SetBorderSize(0);
        leg2->SetLineColor(kWhite);
        leg2->SetFillStyle(0);
        gStyle->SetStatStyle(0);
        gStyle->SetTitleStyle(0);
        gROOT->ForceStyle();
        leg2->AddEntry(hRoptFree,"Free #beta","f");
        leg2->AddEntry(hRoptFixed,"Fixed #beta","f");
        leg2->Draw();
	TLatex *   tex2 = new TLatex(292.6503,8583.207,"r_{opt} = 296 m");
	tex2->SetTextAlign(13);
	tex2->SetTextColor(kBlue-4);
	tex2->SetTextFont(42);
	tex2->SetTextSize(0.03);
	tex2->SetLineWidth(2);
	tex2->Draw();
	TLatex *   tex1 = new TLatex(289.650,8500,"r_{opt} = 289 m");
	tex1->SetTextColor(kOrange-4);
	tex1->SetTextAlign(13);
	tex1->SetTextFont(42);
	tex1->SetTextSize(0.03);
	tex1->SetLineWidth(2);
	tex1->Draw("SAME");
   */ 
    
          //Ropt con simulazioni di Corsika
    
          TCanvas  *canvasRoptSim          = new TCanvas("canvasSim","cRoptSim",900,700);
          TH1F *hRoptProt = (TH1F *) histoFixedBetaFullEff_simProton->Get("hRoptSatur1");//non saturati beta fisso
          TH1F *hRoptIron = (TH1F *) histoFixedBetaFullEff_simIron->Get("hRoptSatur1");//non saturati beta fisso
          TH1F *hRoptHadrons = (TH1F *) histoFixedBetaFullEff_simHadrons->Get("hRoptSatur1");//non saturati beta fisso
	  hRoptProt->GetXaxis()->SetTitleOffset(1.30);
          hRoptProt->GetYaxis()->SetTitleOffset(1.6);
          hRoptProt->GetXaxis()->CenterTitle();
          hRoptProt->GetYaxis()->SetRangeUser(0,5000);
          hRoptProt->GetXaxis()->SetRangeUser(0,450);
          hRoptProt->GetYaxis()->CenterTitle();
          hRoptProt->GetYaxis()->SetTitle("Events ");
          hRoptProt->GetXaxis()->SetTitle("r_{opt} / m");
 	  hRoptProt->GetYaxis()->SetLabelOffset(0.014);
          hRoptProt->GetXaxis()->SetLabelOffset(0.014);
 	  hRoptProt->GetYaxis()->SetTitleSize(0.05);
          hRoptProt->GetXaxis()->SetTitleSize(0.05);
 	  hRoptProt->GetYaxis()->SetLabelSize(0.045);
          hRoptProt->GetXaxis()->SetLabelSize(0.045);
          hRoptProt->GetXaxis()->SetTickLength(-0.010);
          hRoptProt->GetYaxis()->SetTickLength(-0.010);
	  hRoptProt->GetXaxis()->SetNdivisions(10);	  
	  hRoptProt->GetYaxis()->SetNdivisions(10);
  	  hRoptProt->SetLineWidth(3);
  	  hRoptIron->SetLineWidth(3);
  	  hRoptHadrons->SetLineWidth(3);
  	  hRoptProt->SetLineColor(kRed-7);
  	  hRoptIron->SetLineColor(kBlue-7);
  	  hRoptHadrons->SetLineColor(kTeal-5);
  	  cout << "Protons mean = " << hRoptProt->GetMean() << " Entries " << hRoptProt->GetEntries() << "Irons mean = " << hRoptIron->GetMean() << " Entries " << hRoptIron->GetEntries() <<  " hadrons mean = " << hRoptHadrons->GetMean() << " Entries " << hRoptHadrons->GetEntries() << endl;
	  hRoptProt->Draw("");
	  hRoptIron->Draw("SAME");
	  hRoptHadrons->Draw("SAME");
          gStyle->SetStatStyle(0);
          gStyle->SetTitleStyle(0);
          gROOT->ForceStyle();
          TLine *line1 = new TLine(307.,0,307.,5000);
          line1->SetLineColor(kTeal-5);
          line1->SetLineStyle(9);
          line1->SetLineWidth(3);
          line1->Draw("SAME");
          TLine *line2 = new TLine(307.,0,309.,5000);
          line2->SetLineColor(kRed-7);
          line2->SetLineStyle(7);
          line2->SetLineWidth(3);
          line2->Draw("SAME");
          TLine *line3 = new TLine(308.,0,308.,5000);
          line3->SetLineColor(kBlue-7);
          line3->SetLineStyle(9);
          line3->SetLineWidth(3);
          line3->Draw("SAME");
      TLegend *leg2 = new TLegend(0.72049,0.7884615,0.9387528,0.8727811,NULL,"brNDC");
	  leg2->SetBorderSize(0);
	  leg2->SetLineColor(kWhite);
          leg2->SetFillStyle(0);
          gStyle->SetStatStyle(0);
          gStyle->SetTitleStyle(0);
          gROOT->ForceStyle();
          leg2->AddEntry(hRoptProt,"p","f");
          leg2->AddEntry(hRoptIron,"Fe","f");
          leg2->AddEntry(hRoptHadrons,"p + Fe","f");
          leg2->Draw();
	  
	  //protons
	  TLatex *   tex7 = new TLatex(376.662,1723.66,"Simulations");
	  tex7->SetTextAlign(13);
          tex7->SetTextColor(kBlack);
          tex7->SetTextFont(42);
          tex7->SetTextSize(0.05);
          tex7->SetLineWidth(2);
 	  tex7->Draw();
	  TLatex *   tex2 = new TLatex(292.6503,2583.207,"r_{opt} = 308 m");
	  tex2->SetTextAlign(13);
          tex2->SetTextColor(kRed-7);
          tex2->SetTextFont(42);
          tex2->SetTextSize(0.05);
          tex2->SetLineWidth(2);
 	  tex2->Draw();
	  //iron
	  TLatex *   tex1 = new TLatex(292.6503,2583.207,"r_{opt} = 307 m");
          tex1->SetTextColor(kBlue-7);
	  tex1->SetTextAlign(13);
          tex1->SetTextFont(42);
          tex1->SetTextSize(0.05);
          tex1->SetLineWidth(2);
 	  tex1->Draw("SAME");
	  //hadrons
	  TLatex *   tex3 = new TLatex(292.6503,2583.207,"r_{opt} = 308 m");
          tex3->SetTextColor(kTeal-5);
	  tex3->SetTextAlign(13);
          tex3->SetTextFont(42);
          tex3->SetTextSize(0.05);
          tex3->SetLineWidth(2);
 	  tex3->Draw("SAME");
 
      
        /*
	
          TCanvas  *canvasRoptSatur          = new TCanvas("canvasRoptSatur","cRoptSatur",900,700);
	  
          TH1F *hRoptSatur = (TH1F *) file250->Get("hRoptSatur0"); //solo saturati
          TH1F *hRoptNoSatur = (TH1F *) file250->Get("hRoptSatur1");//non saturati 
          TH1F *hRoptAllSatur = (TH1F *) file250->Get("hRoptSatur2"); //tutti gli eventi
	  hRoptAllSatur->GetXaxis()->SetTitleOffset(1.3);
          hRoptAllSatur->GetYaxis()->SetTitleOffset(1.5);
          hRoptAllSatur->GetYaxis()->SetRangeUser(0,4000);
          hRoptAllSatur->GetXaxis()->SetRangeUser(150,550);
          hRoptAllSatur->GetXaxis()->CenterTitle();
          hRoptAllSatur->GetYaxis()->CenterTitle();
          hRoptAllSatur->GetYaxis()->SetTitle("normalized entries ");
          hRoptAllSatur->GetXaxis()->SetTitle("r_{opt} / m");
 	  hRoptAllSatur->GetYaxis()->SetLabelSize(0.031);
          hRoptAllSatur->GetXaxis()->SetLabelSize(0.031);
          hRoptAllSatur->GetXaxis()->SetTickLength(0.02);
          hRoptAllSatur->GetYaxis()->SetTickLength(0.02);
	  Double_t norm = 1;
          hRoptAllSatur->Scale(norm/hRoptAllSatur->Integral());//normalize by the “integral” to show the frequency probability in each bin
          hRoptSatur->Scale(norm/hRoptSatur->Integral());//normalizar por el numero de eventos...si normalizas por Integral(), cada bin te va a quedar con altura igual a la proporcion del total que represente ese bin
          Double_t norm = 1;
          hRoptAllSatur->Scale(norm/hRoptAllSatur->Integral("width"));//normalize by the “integral * bin width” to show the estimated probability density function
          hRoptNoSatur->Scale(norm/hRoptNoSatur->Integral("width"));//normalizar por el area del histograma
          hRoptSatur->Scale(norm/hRoptSatur->Integral("width"));
	  hRoptAllSatur->GetXaxis()->SetNdivisions(10);	  
	  hRoptAllSatur->GetYaxis()->SetNdivisions(10);
  	  hRoptAllSatur->SetLineColor(kOrange-4);
  	  hRoptAllSatur->SetLineWidth(3.0);
  	  hRoptNoSatur->SetLineWidth(3.0);
  	  hRoptSatur->SetLineWidth(3.0);
//  	  hRoptNoSatur->SetLineColor(kBlack);
  	  hRoptSatur->SetLineColor(kBlue-9);
	  hRoptAllSatur->Draw("");
//	  hRoptNoSatur->Draw("SAME");
	  hRoptSatur->Draw("SAME");
          TLine *line1 = new TLine(286.2,0,286.2,0.0269);
          line1->SetLineColor(kOrange-4);
          line1->SetLineStyle(7);
          line1->SetLineWidth(3);
          line1->Draw();
          TLine *line = new TLine(365.1,0,365.1,0.0269);
          line->SetLineColor(kBlue-9);
          line->SetLineStyle(7);
          line->SetLineWidth(3);
          line->Draw();

          TLatex *   tex = new TLatex(260.7925,-168.8438,"#color[2]{270}");
	  tex->SetTextAlign(13);
          tex->SetTextFont(42);
          tex->SetTextSize(0.03);
          tex->SetLineWidth(2);
//        tex->Draw();
          tex = new TLatex(0,0,"");

          TPaveText *pt = new TPaveText(0.4312607,0.9342405,0.5687393,0.995,"blNDC");
	  pt->SetName("title");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextFont(42);
	  text = pt->AddText("blabla");
	  pt->Draw();
	  canvasRoptSatur->Modified();
	  canvasRoptSatur->cd();
	  TLegend *leg = new TLegend(0.6536748,0.7292899,0.8797327,0.8786982,NULL,"brNDC");
	  leg->SetLineColor(kWhite);
          leg->SetFillStyle(0);
          gStyle->SetStatStyle(0);
          gStyle->SetTitleStyle(0);
          gROOT->ForceStyle();
          leg->AddEntry(hRoptAllSatur,"57342 events","f");
  //        leg->AddEntry(hRoptNoSatur,"57051 no saturated events","f");
          leg->AddEntry(hRoptSatur,"331 saturated events","f");
          leg->AddEntry(line1,"<r_{opt}> = 286.2 m","l");
          leg->AddEntry(line,"<r_{opt}> = 365.1 m","l");
          leg->Draw();

*/
    
 }
