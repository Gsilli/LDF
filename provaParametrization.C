#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"
#include "TGraph.h"
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
 #include "TObject.h"
 #include <vector>
 #include <string>
 #include <cmath>
 #include "TF1.h"
 #include "TH1.h"
 #include "TList.h"
 #include "TMath.h"
 #include "TROOT.h"
 #include <sstream>
 #include <fstream>




void provaParametrization() {
vector<double> par;
vector<double> parError;
ifstream  inFile;
vector<double> parS;
vector<double> parErrorS;
ifstream  inFileS;
vector<double> par1D;
vector<double> parError1D;
ifstream  inFile1D;
inFile.open("minuitPar.txt");// parametri provenmienti da minuit su profile2D
inFileS.open("fitParRoot.txt");//parametri da root su profile1D
inFile1D.open("minuitParPol2.txt"); //parametri provenienti da minuit2 usando profile1D
double value, valueError, valueS, valueErrorS, value1D, valueError1D;
if (!inFile.is_open())
   {
   cout << "Error opening the file..." << endl;
   }
else cout << "------- File opened -------" << endl;
//double par[6] = {-4.550, 1.160, 3.414, -1.946, -0.994, 0.736};

std::string line;
while (std::getline(inFile, line)) {
      std::istringstream iss(line);   
      iss >> value >> valueError;
      par.push_back(value); 
      parError.push_back(valueError); 
      }

cout << "parameters minuit su profile2D "<< endl;
for (int m=0; m<3; m++){ //par.size(); m++){
	cout << par[m] << " " << parError[m] << endl;
}
std::string lineS;
while (std::getline(inFileS, lineS)) {
      std::istringstream issS(lineS);   
      issS >> valueS >> valueErrorS;
      parS.push_back(valueS); 
      parErrorS.push_back(valueErrorS); 
      }

cout << "parameters root da profile 1D "<< endl;
for (int s=0; s<3; s++){ //par.size(); m++){
	cout << parS[s] << " "  << parErrorS[s] << endl;
}

std::string line1D;
while (std::getline(inFile1D, line1D)) {
      std::istringstream iss1D(line1D);   
      iss1D >> value1D >> valueError1D;
      par1D.push_back(value1D); 
      parError1D.push_back(valueError1D); 
      }

cout << " parameters minuit su profile1D " << endl;
for (int d=0; d<3; d++){ //par.size(); m++){
	cout << par1D[d]  << " "  << parError1D[d] << endl;

}

     gStyle->SetTitleX(0.8); //title X location  
     gStyle->SetErrorX(0);
     gStyle->SetOptTitle(0);
     gStyle->SetStatStyle(0);
     gStyle->SetOptFit(0);
     gStyle->SetOptStat(0);



	TFile *inputFile = new TFile("HistogramsLDF4.root");


	// x: sec(θ) y: logS250 z: beta -----> con projection x ottengo x: secTheta & y:beta
	hprof2d->Print();
	int nybins = hprof2d->GetYaxis()->GetNbins();
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << "---------- ProjectionX ----------nybins = " << nybins  << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	cout << "---------------------------------- " << endl;
	TCanvas *c1 = new TCanvas("c1", "profile1D", 1200, 300);
	TCanvas *c2 = new TCanvas("c2","beta vs secTheta", 1200, 300);
	c2->Divide(4,2);
        TProfile *Prof_S = (TProfile *) inputFile->Get("BetaS1");
        TF1  *f0_pol2 = new TF1("f0_pol2","[0] + [1]*x + [2]*x*x",1,1.4);
        TF1  *fS_pol2 = new TF1("fS_pol2","[0] + [1]*x + [2]*x*x",1,1.4);
        TF1  *f1D_pol2 = new TF1("f1D_pol2","[0] + [1]*x + [2]*x*x",1,1.4);
        TF1  *f2 = new TF1("f","pol2",1,1.4);
        vector<double> vectorCenterBin250;

	for (int ybin=1; ybin <=nybins; ybin++) {
                double logS250_1 = hprof2d->GetYaxis()->GetBinLowEdge(ybin);
                double logS250_2 = hprof2d->GetYaxis()->GetBinUpEdge(ybin);
		c2->cd(ybin);
         	TGraph* gbetaData6 = new TGraph();
	        TGraph* gbetaPar6  = new TGraph();
		stringstream label2;
		double logS250_1 = hprof2d->GetYaxis()->GetBinLowEdge(ybin); 
		double logS250_2 = hprof2d->GetYaxis()->GetBinUpEdge(ybin); 
                vectorCenterBin250.push_back((logS250_1 + logS250_2)/2);	
                double centerLogS = (logS250_1+logS250_2)/2;
		f0_pol2->SetParameter(0,par[0]);
                f0_pol2->SetParameter(1,par[1]);
                f0_pol2->SetParameter(2,par[2]);
		fS_pol2->SetParameter(0,parS[0]);
                fS_pol2->SetParameter(1,parS[1]);
                fS_pol2->SetParameter(2,parS[2]);
		f1D_pol2->SetParameter(0,par1D[0]);
                f1D_pol2->SetParameter(1,par1D[1]);
                f1D_pol2->SetParameter(2,par1D[2]);
		stringstream hname2;
		hname2 << "_px" << ybin; 	
		label2 << "log(S_{250}) = [" << logS250_1 << ", " << logS250_2 << "]";
	
		TH1D *hpx = hprof2d->ProjectionX(hname2.str().c_str(), ybin, ybin);
		//hpx->Print("all");
                hpx->GetXaxis()->SetTitleOffset(1.3);
                hpx->GetYaxis()->SetTitleOffset(1.5);
	        hpx->GetXaxis()->CenterTitle();
		hpx->GetYaxis()->CenterTitle();
		hpx->GetXaxis()->SetTitle("sec(#theta)");
		hpx->GetYaxis()->SetTitle("#beta");
                hpx->SetMaximum(-1.5);
		hpx->SetMinimum(-2.5);
	        hpx->SetMarkerStyle(1);
		hpx->GetYaxis()->SetLabelSize(0.031);
		hpx->GetXaxis()->SetLabelSize(0.031);
		hpx->GetXaxis()->SetTickLength(0.02);
		hpx->GetYaxis()->SetTickLength(0.02);
	        f0_pol2->SetLineColor(kRed-4);
	        fS_pol2->SetLineColor(kOrange-3);
	        f1D_pol2->SetLineStyle(2);
	        f1D_pol2->SetLineColor(kGreen-1);
	        f1D_pol2->SetLineWidth(3);
		f0_pol2->GetXaxis()->SetTitleOffset(1.3);
                f0_pol2->GetYaxis()->SetTitleOffset(1.5);
	        f0_pol2->GetXaxis()->CenterTitle();
		f0_pol2->GetYaxis()->CenterTitle();
		f0_pol2->GetXaxis()->SetNdivisions(10);
		f0_pol2->GetYaxis()->SetNdivisions(10);
		hpx->SetLineColor(kBlue+3);
		hpx->SetMarkerColor(kBlue+3);
		hpx->SetMarkerStyle(21);
		hpx->SetMarkerSize(1.);
		hpx->SetLineWidth(1);
		hpx->GetXaxis()->SetNdivisions(10);
		hpx->GetYaxis()->SetNdivisions(10);
	        f2->SetLineColor(kBlue-7); 
		hpx->Fit(f2,"R+");	
		hpx->Draw();
                f0_pol2->Draw("SAME");
                fS_pol2->Draw("SAME");
                f1D_pol2->Draw("SAME");
		hpx->GetXaxis()->SetTitle("sec(#theta)");
		hpx->GetYaxis()->SetTitle("#beta");

		// Show sec(θ) bin
		TLatex *text2=new TLatex();
		text2->SetTextFont(42);
		text2->SetTextSize(.035);
		text2->DrawLatexNDC(.57, .8, label2.str().c_str());
                leg = new TLegend(0.116335,0.6773374,0.472231,0.8704268,NULL,"brNDC");
		leg->SetLineColor(kWhite);
                leg->SetFillStyle(0);
                gStyle->SetStatStyle(0);
                gStyle->SetTitleStyle(0);
                gROOT->ForceStyle();
                leg->AddEntry(f2,"fit","l");
                leg->AddEntry(f0_pol2,"profile2D minuit par","l"); 
                leg->AddEntry(fS_pol2,"profile1D root fit par","l"); 
                leg->AddEntry(f1D_pol2,"profile1D minuit par","l"); 
                leg->Draw();


	        

	}

		c1->cd();
                Prof_S->GetXaxis()->SetTitleOffset(1.3);
                Prof_S->GetYaxis()->SetTitleOffset(1.5);
	        Prof_S->GetXaxis()->CenterTitle();
		Prof_S->GetYaxis()->CenterTitle();
		Prof_S->GetXaxis()->SetTitle("sec(#theta)");
		Prof_S->GetYaxis()->SetTitle("#beta");
                Prof_S->SetMaximum(-1.5);
		Prof_S->SetMinimum(-2.5);
	        Prof_S->SetMarkerStyle(1);
		Prof_S->GetYaxis()->SetLabelSize(0.031);
		Prof_S->GetXaxis()->SetLabelSize(0.031);
		Prof_S->GetXaxis()->SetTickLength(0.02);
		Prof_S->GetYaxis()->SetTickLength(0.02);
	        f0_pol2->SetLineColor(kRed-4);
	        fS_pol2->SetLineColor(kOrange-3);
	        f1D_pol2->SetLineStyle(2);
	        f1D_pol2->SetLineColor(kGreen-1);
	        f1D_pol2->SetLineWidth(3);
		f0_pol2->GetXaxis()->SetTitleOffset(1.3);
                f0_pol2->GetYaxis()->SetTitleOffset(1.5);
	        f0_pol2->GetXaxis()->CenterTitle();
		f0_pol2->GetYaxis()->CenterTitle();
		f0_pol2->GetXaxis()->SetNdivisions(10);
		f0_pol2->GetYaxis()->SetNdivisions(10);
		Prof_S->SetLineColor(kBlue+3);
		Prof_S->SetMarkerColor(kBlue+3);
		Prof_S->SetMarkerStyle(21);
		Prof_S->SetMarkerSize(1.);
		Prof_S->SetLineWidth(1);
		Prof_S->GetXaxis()->SetNdivisions(10);
		Prof_S->GetYaxis()->SetNdivisions(10);
	        f2->SetLineColor(kBlue-7); 
		Prof_S->Fit(f2,"R+");	
		Prof_S->Draw();
                f0_pol2->Draw("SAME");
                fS_pol2->Draw("SAME");
                f1D_pol2->Draw("SAME");
		Prof_S->GetXaxis()->SetTitle("sec(#theta)");
		Prof_S->GetYaxis()->SetTitle("#beta");

		// Show sec(θ) bin
                leg = new TLegend(0.116335,0.6773374,0.472231,0.8704268,NULL,"brNDC");
		leg->SetLineColor(kWhite);
                leg->SetFillStyle(0);
                gStyle->SetStatStyle(0);
                gStyle->SetTitleStyle(0);
                gROOT->ForceStyle();
                leg->AddEntry(f2,"fit","l");
                leg->AddEntry(f0_pol2,"profile2D minuit par","l"); 
                leg->AddEntry(fS_pol2,"profile1D root fit par","l"); 
                leg->AddEntry(f1D_pol2,"profile1D minuit par","l"); 
                leg->Draw();


}
