// Ajuste de datos con una función de likelihood usando el paquete Minuit2
// ROOT includes
#include <algorithm>
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include <sstream>
#include <fstream>
#include "TProfile2D.h"
#include <TStyle.h>
 #include <TLegend.h>
 #include "rootStyle.h"
#include "TPave.h"
#include "TBox.h"
#include "TFitResult.h"
#include "TMatrixD.h"
//standard header
#include "TObject.h"
#include <vector>
#include <string>
#include <cmath>
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <TLegend.h>
#include "rootStyle.h"
#include "TPave.h"
#include "TBox.h"
//standard header
#include "TObject.h"
#include <vector>
#include <string>
#include <cmath>
#include <TF1.h>
#include <TF2.h>
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
//double par[6] = {-4.550, 1.160, 3.414, -1.946, -0.994, 0.736};
double par[6] = {-2.02998,-0.0563938,1.05906,-0.150488,3.08801,-3.63978};
//double par[6] = {-1.6, -0.5, -1.1, 0.6, 0.9, -0.4};
double parErr[6] = {0.6, 0.4, 1.0, 0.7, 0.4, 0.3};
TFile *pippo  = new TFile("covarianceMatrixMinuitProfModifiedLDFFinder.root","RECREATE");

struct ldfInfos  {
     double *betaS      = new double[10];
     double *secThetaS  = new double[10];
     double *betaErrorS = new double[10];
     double *beta       = new double[80];
     double *secTheta   = new double[10];
     double *betaError  = new double[80];
     double *logS250    = new double[8];
     double *parMin     = new double[3]; //due parametri
     double *parMax     = new double[3]; //due parametri
     unsigned int        ndataX;
     unsigned int        ndataY;
     double *par        = new double[3];
     double *parErr     = new double[3];
     vector<double> binLowEntries;

} profileInfos;

/**
 * Get data from histograms
 *
 * @return {ldfInfos} data 
 */
ldfInfos getData() {
	ldfInfos p;
	TFile *a = new TFile("histogramsLDFRopt300ModifiedLDFFinder.root");//250m
	TFile *a100 = new TFile("histogramsLDFRopt300ModifiedLDFFinder.root");//250m
	TProfile *Prof_S = (TProfile *) a->Get("BetaS1");
	TH2D *ldf2D = (TH2D *) a100->Get("ldfpar2D");
	TProfile2D *Prof_2D = (TProfile2D *) a->Get("hprof2d");
	
        double nX = ldf2D->GetNbinsX();
        double nY = ldf2D->GetNbinsY();
	int counter = 0;
	for (unsigned int i=1; i<=nX; i++) {
	        for (unsigned int j=1; j<=nY; j++) {
				cout << " BinContent "<< ldf2D->GetBinContent(i,j) << " x=" << i << " y=" << j << " " <<(i-1)*nY+j<<endl;
		        /*if(ldf2D->GetBinContent(i,j)<=50){//levo commento se voglio tagli per entrate nei bin
			       	counter++;
				p.binLowEntries.push_back((i-1)*nY+j);
			}*/
		}
	}
	//cout << counter << " counter bin rifiutati per il fit " << endl;
	TF1  *fpol2 = new TF1("fpol2","pol2",1,1.4);
	//Prof_S->Fit(fpol2,"B+");
	cout << "chi2Fit root " << fpol2->GetChisquare() << endl;
	for (unsigned int j=0; j<3; j++) {
    		std::string parName = "par" + std::to_string(j);
    		p.par[j] = fpol2->GetParameter(j);
    		p.parMin[j] = p.par[j] -3*(fpol2->GetParError(j));
    		p.parMax[j] = p.par[j] +3*(fpol2->GetParError(j));
    	}
	p.ndataX = Prof_2D->GetNbinsX();
	p.ndataY = Prof_2D->GetNbinsY();
	for (unsigned int i=1; i<=p.ndataX; i++) {
		for (unsigned int j=1; j<=p.ndataY; j++) {
			p.beta[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinContent(i,j);
	    		p.betaError[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinError(i,j);
	    		p.logS250[j-1] = Prof_2D->GetYaxis()->GetBinCenter(j);//valori di logS sono solo 8 (8bin)
	//		cout << "profile--------------- " <<(i-1)*p.ndataY+j-1 << endl;
		}	      
	    	p.secTheta[i-1] = Prof_2D->GetXaxis()->GetBinCenter(i);//valori di secTheta sono solo 10 (10bin)
	}
	return p;		
}


// Función usada para ajustar los datos 
double fmodel(double *x, double *par)
{
    	double sectheta = x[0];
    	double lgSRef = x[1];
    	double a0 = par[0];
    	double a1 = par[1];
    	double b0 = par[2];
    	double b1 = par[3];
    	double c0 = par[4];
    	double c1 = par[5];

	double y = a0 + a1*lgSRef + sectheta*(b0 + b1*lgSRef)+ pow(sectheta,2)*(c0 + c1*lgSRef);

  	return y;
}

// Función de peso del fit
double fitchi2(const double *fitpar) {
	
	double par0 = fitpar[0];
	double par1 = fitpar[1];
	double par2 = fitpar[2];
	double par3 = fitpar[3];
	double par4 = fitpar[4];
	double par5 = fitpar[5];
	TF2 fitf("fitf",fmodel, 1, 1.4, 1, 2.2, 6); //per Ropt300 ho bin logSref[0.6,2.2] mentre per Ropt250 ho tra Sropt[1,2.6]
	//int tra 1 e 1.4 della secTheta i le passo due parametri

	fitf.SetParameter(0,par0);
	fitf.SetParameter(1,par1);	
	fitf.SetParameter(2,par2);	
	fitf.SetParameter(3,par3);	
	fitf.SetParameter(4,par4);	
	fitf.SetParameter(5,par5);	
	
	double chi2 = 0;
	for (unsigned int i=0; i<profileInfos.ndataX; i++) { //ascissa secTheta da 1 a 1.4 10 bin
		for (unsigned int j=0; j<profileInfos.ndataY; j++) { // ordinata da 1 a 2.6 logS250 8bin
// per cut nei bin con poche entrate---------------------		
//         	if (std::find(profileInfos.binLowEntries.begin(), profileInfos.binLowEntries.end(), (i*profileInfos.ndataY+j)) != profileInfos.binLowEntries.end()) continue; 
	
			if (profileInfos.betaError[i*profileInfos.ndataY+j]==0) continue;
			//cout <<"--------------" << i*profileInfos.ndataY+j << endl;
			double yfit = fitf.Eval(profileInfos.secTheta[i], profileInfos.logS250[j]);
	                double chi2SingleData = pow( (yfit-profileInfos.beta[i*profileInfos.ndataY+j]) / profileInfos.betaError[i*profileInfos.ndataY+j], 2); //i*p.ndataY+j valore di beta per ciascun secTheta e logS250
	//		double chi2SingleData = flike1w(&yfit, &ydata[i]);
			chi2 += chi2SingleData;
		}	
	}
        //cout << " chi2 " << chi2 << endl;
	return chi2;
}

void minuit2Tprofile() {
	std::fstream outFile;
	outFile.open("fitParMinProf6parModifiedLDFFinder.txt", outFile.out);

        profileInfos = getData();

	// Crear el minimizador
	ROOT::Minuit2::Minuit2Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer();  

	// Setear la función a minimizar
        ROOT::Math::Functor f(&fitchi2,2);  //  pasar función por referencia! 
	min->SetFunction(f);
 	min->SetMaxFunctionCalls(100000000); // used by Minuit and Minuit2
       	min->SetMaxIterations(100000000); // used by GSL*
       	min->SetTolerance(0.000001);

	min->SetVariable(0,"par0", par[0], 0.01);
        min->SetVariable(1,"par1", par[1], 0.01);
        min->SetVariable(2,"par2", par[2], 0.01);
        min->SetVariable(3,"par3", par[3], 0.01);
        min->SetVariable(4,"par4", par[4], 0.01);
        min->SetVariable(5,"par5", par[5], 0.01);

/*	min->SetVariableLimits(0,par[0]-3*parErr[0],par[0]+3*parErr[0]);
	min->SetVariableLimits(1,par[1]-3*parErr[1],par[1]+3*parErr[1]);
	min->SetVariableLimits(2,par[2]-3*parErr[2],par[2]+3*parErr[2]);
	min->SetVariableLimits(3,par[3]-3*parErr[3],par[3]+3*parErr[3]);
	min->SetVariableLimits(4,par[4]-3*parErr[4],par[4]+3*parErr[4]);
	min->SetVariableLimits(5,par[5]-3*parErr[5],par[5]+3*parErr[5]);
	
	min->SetVariableLimits(0,-5,-2);
	min->SetVariableLimits(1,0,3);
	min->SetVariableLimits(2,2,5);
	min->SetVariableLimits(3,-4,0);
	min->SetVariableLimits(4,-3, 0);
	min->SetVariableLimits(5,0,2);
	
	
	min->SetUpperLimitedVariable(0, "par0", profileInfos.parMax[0], .1, -3);    	//set par limits 
	min->SetLowerLimitedVariable(0, "par0", profileInfos.parMin[0], .1, -5);    	// setear los parámetros 
	min->SetUpperLimitedVariable(1, "par1", profileInfos.parMax[1], .1, 2);
	min->SetLowerLimitedVariable(1, "par1", profileInfos.parMin[1], .1, 0);
	min->SetUpperLimitedVariable(2, "par2", profileInfos.parMax[2], .1, 4);    	 
	min->SetLowerLimitedVariable(2, "par2", profileInfos.parMin[2], .1, 2);    	 
	min->SetUpperLimitedVariable(3, "par3", profileInfos.parMax[3], .1, 0);
	min->SetLowerLimitedVariable(3, "par3", profileInfos.parMin[3], .1, -3);
	min->SetUpperLimitedVariable(4, "par4", profileInfos.parMax[4], .1, 0);    	 
	min->SetLowerLimitedVariable(4, "par4", profileInfos.parMin[4], .1, -2);    	 
	min->SetUpperLimitedVariable(5, "par5", profileInfos.parMax[5], .1, 2);
	min->SetLowerLimitedVariable(5, "par5", profileInfos.parMin[5], .1, 0);
        *///Set the free variables to be minimized!

    
	min->SetPrintLevel(4);   // incrementar para debuggear

	// Ejecutar la minimización
 	min->Minimize();

	// Obtener los parámetros
	ROOT::Minuit2::MnUserParameterState parameters = min->State();
	double mu0hat     = parameters.Value(0);
	double alphahat   = parameters.Value(1);	
	double gammahat   = parameters.Value(2);	
	double deltahat   = parameters.Value(3);	
	double epsilonhat = parameters.Value(4);	
	double zetahat    = parameters.Value(5);	
	double mu0hatError     = parameters.Error(0);
	double alphahatError   = parameters.Error(1);	
	double gammahatError   = parameters.Error(2);	
	double deltahatError   = parameters.Error(3);	
	double epsilonhatError = parameters.Error(4);	
	double zetahatError    = parameters.Error(5);	
        outFile << mu0hat << " " << mu0hatError << "\n"<<  alphahat << " " << alphahatError <<"\n" << gammahat << " " << gammahatError<< "\n" << deltahat << " " << deltahatError<< "\n" << epsilonhat << " " << epsilonhatError<< "\n" << zetahat << " " << zetahatError << endl; 



	// Matriz de covariancia - Only supported from ROOT version 6
//         TMatrixD matrix(npar,npar,fitter->GetCovarianceMatrix());
	const ROOT::Minuit2::MnUserCovariance& cova = parameters.Covariance();
	TMatrixD covMatrix(cova.Nrow(),cova.Nrow());

	cout << "Covariance Matrix" << endl;
	for (int i = 0; i < cova.Nrow(); i++) {
		//cout << "par" << i;
		for ( int j = 0; j < cova.Nrow(); j++) {
			covMatrix(i,j) = cova.operator()(i,j);
			//cout << "    " << cova.operator()(i,j);
		}
		//cout << endl;
	}
        cout <<  "-----------  covariance matrix ----------- "<< endl;
	covMatrix.Print();
        TMatrixD transpCov(6,6);
        transpCov = covMatrix.T();
        cout <<  "-----------  covariance  trasposta matrix ----------- "<< endl;
        transpCov.Print();
        pippo->WriteObject(&covMatrix, "covMatr_minuitProf");
	outFile.close();

}

