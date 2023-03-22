// Ajuste de datos con una función de likelihood usando el paquete Minuit2
// ROOT includes
#include "TObject.h"
#include "TVirtualFitter.h" 
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
#include <TH2.h>
#include <TF2.h>
 #include "TH1.h"
 #include "TMath.h"
 #include "TROOT.h"
 #include "TProfile.h"
 #include "TFile.h"
#include "TFitResult.h"
#include "TMatrixD.h"

//double par[6] = {-4.550, 1.160, 3.414, -1.946, -0.994, 0.736};//messina
double par[6] = {-2.02998,-0.0563938,1.05906,-0.150488,3.08801,-3.63978};
//   TFile *a = new TFile("HistogramsLDFNewOffline_LDF2.root");//250m
   //TFile *a = new TFile("histogramsLDFRopt300.root");//300 ropt con ldffinder modificato billoir condition cambiata a else if 300 (prima stava in 250)
   TFile *a = new TFile("histogramsLDFRopt300ModifiedLDFFinder.root");//300m
//vecchia rec con vecchio offline (feb 2020)
//   TFile *a = new TFile("HistogramsLDF4.root");//250m
   TProfile *Prof_S = (TProfile *) a->Get("BetaS1");
   TProfile2D *Prof_2D = (TProfile2D *) a->Get("hprof2d");
   TProfile *Prof_S2 = (TProfile *) a->Get("BetaS2"); 
   TFile *pippo  = new TFile("covarianceMatrixRootProfModifieLDFFinder.root","RECREATE");

// Función usada para ajustar los datos - ley de potencias
double fmodel(double *x, double *par)
{

    double sectheta = x[0];
    double lgSRef = x[1];
    double par0 = par[0];
    double par1 = par[1];
    double par2 = par[2];
    double par3 = par[3];
    double par4 = par[4];
    double par5 = par[5];

  double y = par0 + par1*lgSRef + sectheta*(par2 + par3*lgSRef) + (par4 + par5*lgSRef)*pow(sectheta,2);
//	double y = par0 + par1 * sectheta + par2*pow(sectheta,2);
  return y;
}
double fmodel1(double *x, double *par)
{

    double sectheta = x[0];
    double lgSRef = 1.15; //centro intervallo fisso di beta_S2
    double par0 = par[0];
    double par1 = par[1];
    double par2 = par[2];
    double par3 = par[3];
    double par4 = par[4];
    double par5 = par[5];

  double y1 = par0 + par1*lgSRef + sectheta*(par2 + par3*lgSRef) + (par4 + par5*lgSRef)*pow(sectheta,2);
  return y1;
}



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
	unsigned int        ndataY;
	unsigned int        ndataX;
	double *par        = new double[3];
	double *parErr     = new double[3];
} profileInfos;


ldfInfos getData() {
std::fstream outFile2;
outFile2.open("fitParRootpol2.txt", outFile2.out);
   ldfInfos p;
   TF1 *fpol2 = new TF1("fpol2","pol2", 1, 1.4);
   Prof_S->Fit(fpol2);
   p.ndataY = Prof_2D->GetNbinsY();
   p.ndataX = Prof_2D->GetNbinsX();									                                                      
   cout << "chi2Fit root " << fpol2->GetChisquare() << endl;
   for (unsigned int j=0; j<3; j++) {
       std::string parName = "par" + std::to_string(j);
       p.par[j] = fpol2->GetParameter(j);
       p.parErr[j] = fpol2->GetParError(j);
       outFile2 << p.par[j] <<" " << p.parErr[j] << " " <<  parName << endl;
       p.parMin[j] = p.par[j] -3*(fpol2->GetParError(j));
       p.parMax[j] = p.par[j] +3*(fpol2->GetParError(j));
       cout << "..........." << parName <<" "<< p.par[j]<< endl;
       cout << parName << " min=" << p.parMin[j] << endl;
       cout << parName << " max=" << p.parMax[j] << endl;
       }   


       for (unsigned int i=1; i<=p.ndataX; i++) {
	       for (unsigned int j=1; j<=p.ndataY; j++) {
		       p.beta[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinContent(i,j);
		       p.betaError[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinError(i,j);
		       p.logS250[j-1] = Prof_2D->GetYaxis()->GetBinCenter(j);//valori di logS sono solo 8 (8bin)
		       cout << (i-1)*p.ndataY+j-1 << endl;
	  	}
		p.secTheta[i-1] = Prof_2D->GetXaxis()->GetBinCenter(i);//valori di secTheta sono solo 10 (10bin)
	}

   for (unsigned int i=1; i<=p.ndataX; i++) {
       p.betaS[i-1] = Prof_S->GetBinContent(i);
       p.betaErrorS[i-1] = Prof_S->GetBinError(i);
       p.secThetaS[i-1] = Prof_S->GetBinCenter(i);
       cout << "betaS---profile1D " << p.betaS[i-1] << "  betaErrorS " << p.betaErrorS[i-1] << " secThetaS "  <<p.secThetaS[i-1] << endl;
       }
  
       return p;

  }













void fit2() {
std::fstream outFile;
outFile.open("fitParRootProf6parModifieLDFFinder.txt", outFile.out);
profileInfos = getData();

   TF2 *f = new TF2("myfunc",fmodel, 1, 1.4, 1, 2.2, 6);//Ropt300 allora uso LogS[0.6,2.2] mentre per Ropt250 usavo LSref[1,2.6]
   TF1 *f1= new TF1("myfunc1",fmodel1, 1, 1.4, 6);
   f->SetParameter(0,par[0]);
   f->SetParameter(1,par[1]);
   f->SetParameter(2,par[2]);
   f->SetParameter(3,par[3]);
   f->SetParameter(4,par[4]);
   f->SetParameter(5,par[5]);
   f1->SetParameter(0,par[0]);
   f1->SetParameter(1,par[1]);
   f1->SetParameter(2,par[2]);
   f1->SetParameter(3,par[3]);
   f1->SetParameter(4,par[4]);
   f1->SetParameter(5,par[5]);
   //Prof_2D->Fit(f);
   Prof_2D->Fit(f,"W");

   Prof_S2->Fit(f1,"S"); 

   TFitResultPtr r = Prof_2D->Fit(f, "S");
   TFitResultPtr p = Prof_S2->Fit(f1, "S");
   TMatrixD corBetaS2 = p->GetCorrelationMatrix();
   TMatrixD cor = r->GetCorrelationMatrix();
   TMatrixD cov = r->GetCovarianceMatrix();
   cout <<  "-----------  covariance matrix ----------- "<< endl;
   cov.Print();
   TMatrixD transpCov(6,6);
   transpCov = cov.T();
   cout <<  "-----------  covariance  trasposta matrix ----------- "<< endl;
   transpCov.Print();
   cout <<  "-----------  correlation matrix ----------- "<< endl;
   cor.Print();
   cout <<  "-----------  correlation matrix betaS2----------- "<< endl;
   corBetaS2.Print();
   
   TVirtualFitter *fitter = TVirtualFitter::GetFitter();
   TMatrixD matrix(6,6,fitter->GetCovarianceMatrix());
   Double_t errorFirstPar = fitter->GetCovarianceMatrixElement(0,0);
   pippo->WriteObject(&cov, "covMatr_rootProf");

   for (int i=0; i<6; i++) {
       outFile << f->GetParameter(i) << " " <<  f->GetParError(i) <<  endl; 
       }
 
  double chi2 = 0;
  for (unsigned int i=0; i<profileInfos.ndataX; i++) { //ascissa secTheta da 1 a 1.4 10 bin
 	 for (unsigned int j=0; j<profileInfos.ndataY; j++) { // ordinata da 1 a 2.6 logS250 8bin
	 	if (profileInfos.betaError[i*profileInfos.ndataY+j]==0) continue;
	     	double yfit = f->Eval(profileInfos.secTheta[i], profileInfos.logS250[j]);
		double chi2SingleData = pow( (yfit-profileInfos.beta[i*profileInfos.ndataY+j]) / profileInfos.betaError[i*profileInfos.ndataY+j], 2);
                chi2 += chi2SingleData;
         }
  }
  cout << " chi2 a mano " << chi2 << endl;
  cout << f->GetChisquare()<<" chi quadrato root " << endl;
  cout << f->GetNDF()<<" ndf root " << endl;


  outFile.close();

}

