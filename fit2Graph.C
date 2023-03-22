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
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"



double par[6] = {-4.550, 1.160, 3.414, -1.946, -0.994, 0.736};
//double par[6] = {-1.6, -0.5, -1.1, 0.6, 0.9, -0.4};
double parErr[6] = {0.6, 0.4, 1.0, 0.7, 0.4, 0.3};
TFile *pippo  = new TFile("covarianceMatrixRootGraph.root","RECREATE");

struct ldfInfos  {
                    double *betaS      = new double[32419];
                    double *secThetaS  = new double[32419];
                    double *betaErrorS = new double[32419];
                    double *beta       = new double[32419];
                    double *secTheta   = new double[32419];
                    double *betaError  = new double[32419];
                    double *logS250    = new double[32419];
                    double *parMin     = new double[3]; //due parametri
                    double *parMax     = new double[3]; //due parametri
                    double *par        = new double[3];
                    unsigned int ndata;
                     } profileInfos;


ldfInfos getData() {
	        ldfInfos p;
		std::fstream outFile2;
	        TFile *a = new TFile("HistogramsLDFNewOffline_LDF2.root");//250m
		TGraph2DErrors *gr = (TGraph2DErrors *) a->Get("graphErrLdf2D");
		int n = gr->GetN();
		p.ndata = n;
                cout << "N= " << n << endl;
                for (int j=0; j<=n; j++) {
                    p.logS250[j] = (gr->GetY())[j];
                    p.secTheta[j] = (gr->GetX())[j];
                    p.beta[j] = (gr->GetZ())[j];
                    p.betaError[j] = (gr->GetEZ())[j];
                    cout << p.logS250[j]  << " " << p.secTheta[j] << " " << p.beta[j] << " " << p.betaError[j] << endl;
                    }

		return p;
}

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
  return y;
}




void fit2Graph() {
std::fstream outFile;
outFile.open("fitParRootGraph.txt", outFile.out);
profileInfos = getData();

   TF2 *f = new TF2("myfunc",fmodel, 1, 1.4, 1, 2.6, 6);
   f->SetParameter(0,par[0]);
   f->SetParameter(1,par[1]);
   f->SetParameter(2,par[2]);
   f->SetParameter(3,par[3]);
   f->SetParameter(4,par[4]);
   f->SetParameter(5,par[5]);
   TFile *a = new TFile("HistogramsLDFNewOffline_LDF2.root");//250m
   TGraph2DErrors *gr = (TGraph2DErrors *) a->Get("graphErrLdf2D");
   gr->Fit(f);
   TFitResultPtr r = gr->Fit(f, "S");
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
   
   TVirtualFitter *fitter = TVirtualFitter::GetFitter();
   TMatrixD matrix(6,6,fitter->GetCovarianceMatrix());
   Double_t errorFirstPar = fitter->GetCovarianceMatrixElement(0,0);
   pippo->WriteObject(&cov, "covMatr_rootGraph");

   cout << "element 00 "<< fitter->GetCovarianceMatrixElement(1,1) << endl;
   for (int i=0; i<6; i++) {
       outFile << f->GetParameter(i) << " " <<  f->GetParError(i) <<  endl; 
       }
 
   double chi2 = 0;
   for (unsigned int i=0; i<profileInfos.ndata; i++) {
       if(profileInfos.beta[i]>=0)continue;
  //     cout << profileInfos.beta[i] << " profileInfos.beta " << i << endl;
       double yfit = f->Eval(profileInfos.secTheta[i], profileInfos.logS250[i]);
       double chi2SingleData = pow( (yfit-profileInfos.beta[i]) / profileInfos.betaError[i], 2);
       chi2 += chi2SingleData;
       }
  
  cout << " chi2 " << chi2 << " " << f->GetNDF() << " " << f->GetChisquare() <<  endl;











  outFile.close();

}

