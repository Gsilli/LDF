#include <map>
#include "TMatrixTBase.h"
#include "TMatrixDBasefwd.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TMatrixTLazy.h"
#include "TMatrixTCramerInv.h"
#include "TDecompLU.h"
#include "TMatrixDEigen.h"
#include "TClass.h"
#include "TMath.h"
#include "Math/BinaryOpPolicy.h"
#include "Math/Expression.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"
#include "TGraph.h"
#include <string.h>

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
 #include <TH2.h>
 #include "TMath.h"
 #include "TROOT.h"
 #include "TProfile.h"
 #include "TFile.h"
 #include <cmath>
 #include "TMatrixDfwd.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TMatrixTLazy.h"
#include "TMatrixTCramerInv.h"
#include "TDecompLU.h"
#include "TMatrixDEigen.h"
#include "TClass.h"
#include "TMath.h"
	
gStyle->SetTitleX(0.8); //title X location  
gStyle->SetErrorX(0);
gStyle->SetOptTitle(0);
gStyle->SetStatStyle(0);
gStyle->SetOptFit(0);
gStyle->SetOptStat(0);

void fitBetaPerDiego (){
       /* vector<double> vecPar;
	vecPar.push_back(-1.64518);
	vecPar.push_back(-0.282937);
	vecPar.push_back(-1.73762);
	vecPar.push_back(0.560785);
	vecPar.push_back(1.20272);
	vecPar.push_back(-0.345354);
	vecPar.push_back(0.0118886);
	*/ 
        vector<double> vecPar;
	vecPar.push_back(-2.02998);
	vecPar.push_back(-0.0563938);
	vecPar.push_back(1.05906);
	vecPar.push_back(-0.150488);
	vecPar.push_back(3.08801);
	vecPar.push_back(-3.63978);
	vecPar.push_back(1.24104);

	int nybins = 8;
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << "------------------ ProvaFit BetavsSecTheta  -------------nybins = " << nybins  << endl;
	cout << "------------------------------------------------------- " << endl;
cout << " Uso i parametri :\n p0=" << vecPar[0] <<  "\n p1=" << vecPar[1] << "\n p2=" << vecPar[2] << "\n p3=" << vecPar[3] << "\n p4=" << vecPar[4] << "\n p5=" << vecPar[5] << "\n p6=" << vecPar[6]endl;
	cout << "---------------------------------------------------------------------------- " << endl;
	cout << "---------------------------------------------------------------------------- " << endl;
	cout << "---------------------------------------------------------------------------- " << endl;

	TCanvas *cProva = new TCanvas("fitBeta","betaVsSecTheta", 1200, 300);
	cProva->Divide(4,2);
        TF1  *fProva = new TF1("fProva","[0] + [1]*(x-1.2) + [2]*pow((x-1.2),2)",1,1.4);
        vector<double> logS300;
        logS300.push_back(0.6);
        logS300.push_back(0.8);
        logS300.push_back(1.0);
        logS300.push_back(1.2);
        logS300.push_back(1.4);
        logS300.push_back(1.6);
        logS300.push_back(1.8);
        logS300.push_back(2.0);
        logS300.push_back(1.2);
//	{0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2};


	for (int ybin=1; ybin<=8; ybin++) {

		cProva->cd(ybin);
		stringstream labelProva;
		double logS300_1 = logS300[ybin-1]; 
		double logS300_2 = logS300[ybin]; 
                double centerLogS = (logS300_1+logS300_2)/2;
		double A = vecPar[0] + vecPar[1]*((logS300_1+logS300_2)/2);
                double B = vecPar[2] + vecPar[3]*((logS300_1+logS300_2)/2);
        	double C = vecPar[4] + vecPar[5]*((logS300_1+logS300_2)/2)+ vecPar[6]*pow((logS300_1+logS300_2)/2,2);
		labelProva << "log(S_{300}) = [" << logS300_1 << ", " << logS300_2 << "]";
		cout << " intervallo " << "bin " << labelProva.str() << endl;
		cout << "A = p[0] + p[1]LogS " << A << endl;
		cout << "B = p[2] + p[3]LogS " << B << endl;
		cout << "C = p[4] + p[5]LogS + p[6]logS^2 " << C << endl;
		cout << " " << endl;
		cout << " " << endl;
		fProva->SetParameter(0,A);
                fProva->SetParameter(1,B);
                fProva->SetParameter(2,C);
		cout << " getParameter "<< fProva->GetParameter(0);
		cout << " eval f0Pol2 in 1.2 "<< fProva->Eval(1.2) << endl;
		double evalFunc = fProva->Eval(1.2);
                fProva->GetXaxis()->SetTitleOffset(1.3);
                fProva->GetYaxis()->SetTitleOffset(1.5);
	        fProva->GetXaxis()->CenterTitle();
		fProva->GetYaxis()->CenterTitle();
		fProva->GetXaxis()->SetTitle("sec(#theta)");
		fProva->GetYaxis()->SetTitle("#beta");
                fProva->SetMaximum(-1.5);
		fProva->SetMinimum(-2.5);
		fProva->GetYaxis()->SetLabelSize(0.031);
		fProva->GetXaxis()->SetLabelSize(0.031);
		fProva->GetXaxis()->SetTickLength(0.02);
		fProva->GetYaxis()->SetTickLength(0.02);
	        fProva->SetLineColor(kPink+5);
		fProva->SetMarkerStyle(21);
		fProva->SetMarkerSize(1.);
		fProva->SetLineWidth(2);
		fProva->GetXaxis()->SetNdivisions(10);
		fProva->GetYaxis()->SetNdivisions(10);
		fProva->DrawCopy();
		TLine *l=new TLine(1.2,-2.5,1.2,-1.5);
		l->SetLineWidth(2);
		l->SetLineColor(kOrange-3);
		l->SetLineStyle(7);
		l->Draw();
		TLine *l2=new TLine(1,evalFunc,1.4,evalFunc);
		l2->SetLineColor(kBlue-7);
		l2->SetLineWidth(2);
		l2->SetLineStyle(7);
	        l2->Draw();
		stringstream ss;
		ss <<  evalFunc << endl;
   		TLegend *leg = new TLegend(0.1152015,0.7449489,0.4121594,0.8707197,NULL,"brNDC");
		leg->SetLineColor(kWhite);
		leg->SetFillStyle(0);
		gStyle->SetStatStyle(0);
		gStyle->SetTitleStyle(0);
		gROOT->ForceStyle();
		leg->AddEntry(l2,ss.str().c_str(),"l");
		leg->Draw();

	}
}


