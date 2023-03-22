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

//double par[6] = {-4.550, 1.160, 3.414, -1.946, -0.994, 0.736};
double par[7] = {-2.02998,-0.0563938,1.05906,-0.150488,3.08801,-3.63978,1.24104};
//   TFile *a = new TFile("HistogramsLDFNewOffline_LDF2.root");//250m
   TFile *a = new TFile("histogramsLDFRopt300.root");//300m
//vecchia rec con vecchio offline (feb 2020)
//   TFile *a = new TFile("HistogramsLDF4.root");//250m
   TProfile *Prof_S = (TProfile *) a->Get("BetaS1");
   TProfile2D *Prof_2D = (TProfile2D *) a->Get("hprof2d");
   TProfile *Prof_S2 = (TProfile *) a->Get("BetaS2"); 
   TFile *pippo  = new TFile("covarianceMatrixRootProf.root","RECREATE");

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
    double par6 = par[6];
    double y =  par0 + par1*lgSRef + par2*(sectheta-1.2) + par3*lgSRef*(sectheta-1.2) + par4* (sectheta-1.2)*(sectheta-1.2) + par5*lgSRef*(sectheta-1.2)*(sectheta-1.2) + par6*lgSRef*lgSRef*(sectheta-1.2)*(sectheta-1.2); 
 // double y = par0 + par1*lgSRef + sectheta-1.2*(par2 + par3*lgSRef) + (par4 + par5*lgSRef + par6*pow(lgSRef,2))*pow(sectheta-1.2,2);
//	double y = par0 + par1 * sectheta + par2*pow(sectheta,2);
  return y;
}



struct ldfInfos  {
	double *betaS      = new double[10];
	double *secThetaS  = new double[10];
	double *betaErrorS = new double[10];
	double *beta       = new double[80];
	double *secTheta   = new double[10];
	double *betaError  = new double[80];
	double *logS250    = new double[8];
	unsigned int        ndataY;
	unsigned int        ndataX;
} profileInfos;


ldfInfos getData() {
   ldfInfos p;
   p.ndataY = Prof_2D->GetNbinsY();
   p.ndataX = Prof_2D->GetNbinsX();									                                                      

       for (unsigned int i=1; i<=p.ndataX; i++) {
	       for (unsigned int j=1; j<=p.ndataY; j++) {
		       p.beta[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinContent(i,j);
		       p.betaError[(i-1)*p.ndataY+j-1] = Prof_2D->GetBinError(i,j);
		       p.logS250[j-1] = Prof_2D->GetYaxis()->GetBinCenter(j);//valori di logS sono solo 8 (8bin)
		       cout << (i-1)*p.ndataY+j-1 << endl;
	  	}
		p.secTheta[i-1] = Prof_2D->GetXaxis()->GetBinCenter(i);//valori di secTheta sono solo 10 (10bin)
	}

       return p;

  }













void fit2New() {
std::fstream outFile;
outFile.open("fitParRootProf7par_W.txt", outFile.out);
profileInfos = getData();

   TF2 *f = new TF2("myfunc",fmodel, 1, 1.4, 0.6, 2.2, 7);//Ropt300 allora uso LogS[0.6,2.2] mentre per Ropt250 usavo LSref[1,2.6]
   f->SetParameter(0,par[0]);
   f->SetParameter(1,par[1]);
   f->SetParameter(2,par[2]);
   f->SetParameter(3,par[3]);
   f->SetParameter(4,par[4]);
   f->SetParameter(5,par[5]);
   f->SetParameter(6,par[6]);
   Prof_2D->Fit(f,"VW"); 


   TFitResultPtr r = Prof_2D->Fit(f, "S");
   TMatrixD cor = r->GetCorrelationMatrix();
   TMatrixD cov = r->GetCovarianceMatrix();
   cout <<  "-----------  covariance matrix ----------- "<< endl;
   cov.Print();
   TMatrixD transpCov(7,7);
   transpCov = cov.T();
   cout <<  "-----------  covariance  trasposta matrix ----------- "<< endl;
   transpCov.Print();
   cout <<  "-----------  correlation matrix ----------- "<< endl;
   cor.Print();
   
   TVirtualFitter *fitter = TVirtualFitter::GetFitter();
   TMatrixD matrix(7,7,fitter->GetCovarianceMatrix());
   Double_t errorFirstPar = fitter->GetCovarianceMatrixElement(0,0);
   pippo->WriteObject(&cov, "covMatr_rootProf");

   for (int i=0; i<7; i++) {
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

