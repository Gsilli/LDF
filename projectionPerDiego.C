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
 #include "TROOT.h"
 #include "TRandom3.h"
//
 void projectionPerDiego() {
	gStyle->SetOptStat(0);

//
//	TFile *inputFileldf = new TFile("histogramsPerDiego.root");
//	TFile *inputFileldf = new TFile("histogramsPerDiegoUltimate.root");//nuovo bineado beta libero
//	TFile *inputFileldf = new TFile("histogramsPerDiegoFixedBetaAllyears.root");//nuovo bineado beta fisso 2013-2020(dicembre)
	TFile *inputFileldf = new TFile("histogramsPerDiegoFreeBeta.root");//nuovo bineado beta fisso 2013-2020(dicembre)
	TFile *inputFilemappa = new TFile("mappaRopt_infill4332019.root");//mappa con array  e eventi con ricosgtruiti con Tagli in Ropt
//	TFile *inputFileldf = new TFile("histogramsPerDiego2.root");//nuovo bineado
//	TFile *inputFileldf = new TFile("histogramsPerDiego3.root");//per plot con colori
	TFile *betaErr = new TFile("histogramsLDFRopt300.root");//Ropt300
	TFile *betaErrTh2 = new TFile("histogramsLDF_freeBeta300BetaTh2Bin.root");//Ropt300
	//TFile *ldfFile300 = new TFile("histogramsLDFRopt300BinnatoAlto.root");//Ropt300
	TFile *ldfFile300 = new TFile("histogramsLDFRopt300.root");//Ropt300

	TCanvas *c1 = new TCanvas("c1","Beta vs logS300", 1800, 1500);
	TCanvas *c2 = new TCanvas("c2","Beta vs logS300_2", 1800, 1500);
	TCanvas  *canvasSigma   = new TCanvas("canvasSigma","Sigma",900,700);
	
	
        gStyle->SetPalette(54);
        gStyle->SetNumberContours(255);

	TH2D *ldf2D300          = (TH2D *) ldfFile300->Get("ldfpar2D250");
/*	TCanvas  *canvas2DDiego = new TCanvas("betaVsS300","Th2DPerDiego",900,700);
        TCanvas  *canvasProjx   = new TCanvas("canvasProjx","projx",900,700);
	ldf2D300->GetXaxis()->SetTitleOffset(1.3);
	ldf2D300->GetYaxis()->SetTitleOffset(1.3);
	ldf2D300->GetXaxis()->CenterTitle();
	ldf2D300->GetZaxis()->CenterTitle();
	ldf2D300->GetYaxis()->CenterTitle();
	ldf2D300->GetYaxis()->SetTitle("#beta ");
	ldf2D300->GetXaxis()->SetTitle("log_{10}(S_{300}) / VEM");
        ldf2D300->GetZaxis()->SetTitle("N events ");
	ldf2D300->GetYaxis()->SetRangeUser(-2.4, -1.9);
	ldf2D300->GetXaxis()->SetRangeUser(1,2.2);
	ldf2D300->SetMarkerColor(kGray+2);
	ldf2D300->SetLineColor(kGray+2);
	ldf2D300->SetMarkerStyle(8);
	ldf2D300->SetMarkerSize(1);
	ldf2D300->GetZaxis()->SetLabelSize(0.025);
	ldf2D300->GetYaxis()->SetLabelSize(0.025);
	ldf2D300->GetXaxis()->SetLabelSize(0.025);
	ldf2D300->GetXaxis()->SetTickLength(0.018);
	ldf2D300->GetYaxis()->SetTickLength(0.018);
	ldf2D300->GetXaxis()->SetNdivisions(10);
	ldf2D300->GetYaxis()->SetNdivisions(10);
	ldf2D300->SetLineWidth(1);
	ldf2D300->SetLineColor(kBlack);
	canvas2DDiego->cd();
	ldf2D300->Draw("colz");
	gStyle->SetStatStyle(0);
	gStyle->SetTitleStyle(0);
	gROOT->ForceStyle();
	TH1D *px = ldf2D300->ProjectionX("px", 0, 9); // where firstYbin = 0 and lastYbin = 9
        canvasProjx->cd();
        px->Draw();

*/
	// y: beta  x: logS300 
	int nxbins = ldf2D300->GetXaxis()->GetNbins();
	int nybins = ldf2D300->GetYaxis()->GetNbins();
	vector<double> vectorCenterLogS, vector_sigma, vector_sigmaErr;
cout << "nxbins ---------------------" << nxbins << endl;
	TF1 *gauss = new TF1("gaus","gaus",-3,-1);
	TF1 *f = new TF1("f","expo",1,2.2); //exp(constant + slope * x)
//	TF1 *f = new TF1("f","x[0]*exp([1]*x)",1,2);//Nico
	c1->Divide(3,2);
	c2->Divide(3,2);

	for (int xbin=1; xbin <= nxbins; xbin++) {
//		if(xbin>6) c2->cd(xbin-6);
//		else c1->cd(xbin);
		c1->cd(xbin);
		stringstream label;
		double logS300_1 = ldf2D300->GetXaxis()->GetBinLowEdge(xbin); 
		double logS300_2 = ldf2D300->GetXaxis()->GetBinUpEdge(xbin); 
		label << "[" << logS300_1 << ", " << logS300_2 << "]";
                double centerLog   = (logS300_1+logS300_2)/2;
                vectorCenterLogS.push_back(centerLog);

		cout << "----------------->>>>    log10(S300) bin = " << label.str() << endl;
		stringstream hname;
		hname << "_py" << xbin; 	
		TH1D *hpy = ldf2D300->ProjectionY(hname.str().c_str(), xbin, xbin);
	//	hpy->Print("all");
		hpy->GetXaxis()->SetTitleOffset(1.3);
		hpy->GetYaxis()->SetTitleOffset(1.3);
		hpy->GetXaxis()->CenterTitle();
		hpy->GetZaxis()->CenterTitle();
		hpy->GetYaxis()->CenterTitle();
		hpy->GetYaxis()->SetTitleSize(0.045);
		hpy->GetYaxis()->SetTitleSize(0.045);
		hpy->GetYaxis()->SetLabelSize(0.040);
		hpy->GetXaxis()->SetLabelSize(0.040);
		hpy->GetXaxis()->SetTickLength(0.02);
		hpy->GetYaxis()->SetTickLength(0.02);
		hpy->GetXaxis()->SetNdivisions(10);
		hpy->GetYaxis()->SetNdivisions(10);
		hpy->GetXaxis()->SetTitle("#beta");
		hpy->GetYaxis()->SetTitle("Entries");
		//hpy->SetMinimum(-2.5);
	        TRandom3 rng;
	        Double_t x,y;
	        rng.Rannor(x,y);
	  	hpy->SetMarkerStyle(kFullSquare);
	  	hpy->SetMarkerSize(1.5);
		hpy->SetMarkerColor(kTeal-5);
		hpy->SetLineColor(kTeal-5);	
		hpy->SetLineWidth(2);	
	        gauss->SetLineColor(kBlue-9);		
	        gauss->SetLineWidth(3);		
	        gauss->SetLineStyle(1);		
		hpy->Fit("gaus", "V", "E1", -3, -1);	
		hpy->Draw(" EP");
		if(xbin==1) {   
	          TLegend *leg = new TLegend(0.1544624,0.7469774,0.3532649,0.8727483,NULL,"brNDC");
		  hpy->SetMaximum(850);
	          leg->SetLineColor(kWhite);
		  leg->SetFillStyle(0);
		  gStyle->SetStatStyle(0);
		  gStyle->SetTitleStyle(0);
		  leg->SetTextSize(.045);
		  gROOT->ForceStyle();
		  leg->AddEntry(gaus,"fit","l");
		  leg->AddEntry(hpy,"data","p");
		  leg->Draw();
		}
		// Show logS300 bin
		if(xbin==1){ 
			TLatex *text=new TLatex();
			text->SetTextFont(42);
			text->SetTextSize(.045);	
			text->DrawLatexNDC(.45, .83, "log_{10}(S_{300})/VEM ");
			}
		TLatex *text1=new TLatex();
		text1->SetTextFont(42);
		text1->SetTextSize(.045);
		text1->DrawLatexNDC(.73, .83, label.str().c_str());


//		if(xbin>1){
			//f(x) = p0*exp(-0.5*((x-p1)/p2)^2); il sigma è il parametro 2
			vector_sigma.push_back(gauss->GetParameter(2));
			vector_sigmaErr.push_back(gauss->GetParError(2));
//			}
		}
	
	//plotteo del sigma vs LogS300 
	//x:logS300   y:sigma Gauss
	TGraphErrors *gr    = new TGraphErrors();
	TProfile *deltaBeta = (TProfile *) ldfFile300->Get("betaErr");
	canvasSigma->cd();
	cout << " vector_sigma.size() " << vector_sigma.size() << endl;
	for (int u=0; u<vector_sigma.size(); u++) {
		gr->SetPoint(u, vectorCenterLogS[u], vector_sigma[u]); 
		gr->SetPointError(u, 0, vector_sigmaErr[u]); 
		cout << " bin " << u << " " << vectorCenterLogS[u] << " sigma= " << vector_sigma[u] << endl; }
	gr->GetXaxis()->SetTitleOffset(1.3);
        gr->GetYaxis()->SetTitleOffset(1.2);
	gr->GetXaxis()->SetNdivisions(10);
	gr->GetYaxis()->SetNdivisions(10);
	gr->GetXaxis()->CenterTitle();
	gr->GetYaxis()->CenterTitle();
	gr->GetYaxis()->SetTitleSize(0.06);
	gr->GetXaxis()->SetTitleSize(0.050);
	gr->GetYaxis()->SetLabelSize(0.045);
	gr->GetXaxis()->SetLabelSize(0.045);
	gr->GetXaxis()->SetTickLength(0.02);
	gr->GetYaxis()->SetTitle("#sigma_{#beta}");
	gr->GetXaxis()->SetTitle("log_{10}(S_{300})/VEM");
	gr->GetYaxis()->SetTickLength(0.02);
	gr->SetLineColor(kTeal-5);
	gr->SetMarkerColor(kTeal-5);
	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(2.7);
	gr->SetLineWidth(3);
	f->SetLineColor(kBlue-9);		
	f->SetLineWidth(4);		
	f->SetLineStyle(1);		
	cout << "--------------Fit sigma-----------------" << endl;
	gr->Fit(f);	
	gr->Draw("AP");
	c2->cd();
	deltaBeta->Draw("P");
	f->Draw("SAME l");
	TLegend *leg = new TLegend(0.1544624,0.7469774,0.3532649,0.8727483,NULL,"brNDC");
        leg->SetBorderSize(0);
	leg->SetLineColor(kWhite);
	leg->SetFillStyle(0);
        gStyle->SetStatStyle(0);
	gStyle->SetTitleStyle(0);
	leg->SetTextSize(.045);
	gROOT->ForceStyle();
	leg->AddEntry(f,"fit","l");
	leg->Draw();


}
