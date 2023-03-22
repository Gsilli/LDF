#include "TProfile2D.h"
#include <TProfile.h>
#include <RecEvent.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <vector>
#include <list>
#include <TH1F.h>
#include <TH3.h>
#include "TGraph.h" 
#include "TEfficiency.h"
#include <string>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
using namespace std;


// Fixed size dimensions of array or collections stored in the TTree if any.

   std::vector<TGraphAsymmErrors*> createGraphAsymm (int n);
   vector <double> deltaR_v_433, deltaEn_v_433, deltaS250_v, enRec_v_433, signal_v, rSP_v, signal0_v, rSP0_v, signal1_v, rSP1_v, signal2_v, rSP2_v, signal3_v, rSP3_v;
   TProfile     hprof, hprof0, hprof1, hprof2, hprof3;
   TProfile2D hprof2d, hprof_betaSecTheta, hprof_betaLogS;
   TProfile     hprofRopt, hprofEn0, hprofEn1, hprofEn2, hprofEn3, profBetaError;
   TH1F *DeltaEn, *DeltaCore, *DeltaS250, *hRRel, *hTheta;
   TH3D *ldf3D;
   TH2D *evtDistrib,*ldf2D_S250, *ldf2D, *ldf2DRopt;
   TGraph2DErrors *ldfGraph2D;
   TGraphErrors *ldfPura, *ldfModel;
   std::vector<TGraphAsymmErrors*> myVectorGraphAsymm;
   std::vector<TH1F*> createTH1F (int n, double inf, double sup, int bin, std::string  name);
   std::vector<TProfile*> myVectorAlpha, myVectorProfiles, myVectorProfilesNew, myVectorBetaS, myVectorProfBeta_EnBin, myVectorProfBeta_EnBinNozenithCut;
   std::vector<TH1F*> myVectorHistos, myVectorRoptLogSTheta,  myVectorRopt, myVectorRopt_saturation, myVectorHistosMultEn, myVectorHistosMultZenith;
   vector <double> beta_vector, sec_vector, s250_y, s450_y, s250_tot, s450_tot, s250_bin0, s450_bin0, s250_bin1, s450_bin1, s250_bin2, s450_bin2, s250_bin3, s450_bin3, EnRec_750, EnRec_433, EnRec_FD, En_Err_y_6T5_750, En_Err_y_6T5_433, En_Err_y_750, En_Err_y_433, Total_distRel_x_6T5_750, Total_distRel_x_6T5_433, Total_distRel_x_750, Total_distRel_x_433, sRef33_tot, sRef35_tot, sRef33_bin0, sRef35_bin0, sRef33_bin1, sRef35_bin1,sRef33_bin2, sRef35_bin2,sRef33_bin3, sRef35_bin3;
   vector <double> zenithVector4, zenithVector0, zenithVector1, zenithVector2, zenithVector3, vector_logS250, vector_beta, vector_secTheta, vectorErr_logS250, vectorErr_beta, vectorErr_secTheta;
   double energyCalib(std::string showerSizeLabel, double showerSize, double theta);
vector <int> recEvents750, recEvents433;
   TFile *file;
   TProfile     pTriggerET4, pTriggerETot;
   double zenithBin(vector <double>);
   void book_histograms();
   void graphFunc(vector <double>, vector <double>, char*);
   void graphErrorFunc(vector <double> var_x, vector <double> var_y, vector <double> var_z, vector <double> err_y, vector <double> err_z, char *name);
   void write_histos(std::vector<TProfile*> prof_vector, std::vector<TProfile*> prof_vector2, std::vector<TH1F*> h_vector, std::vector<TH1F*> vec_Ropt, std::vector<TH1F*> vec_RoptSatur, std::vector<TProfile*> prof_vector3, std::vector<TProfile*> prof_vector4, std::vector<TH1F*> vectorMultEn, std::vector<TH1F*> vectorMultZenith, std::vector<TH1F*> vec_RoptLogSTheta, std::vector<TGraphAsymmErrors*> myVectorGraphAsymmTrigSt);
   std::vector<TProfile*> createProfiles(int n, double inf, double sup, int bin, std::string name);
   void FdInfoEye(const RecEvent *recEvent, const unsigned int eyeId, std::fstream &outFile, struct var *Eye_var);
   bool FdCuts(const RecEvent *recEvent, const unsigned int eyeId);
   void binomialErrors(TProfile prof,TGraphAsymmErrors*& g);
   TH1F  *hRopt, *hEnRec_750, *hEnRec_433, *hEn_Fd;
   double energyRec_433, energyRec_750, s450, s250;
   TGraph *graph_Sopt, *gEn_Err433, *gEn6T5_Err433, *gEn_Err750, *gEn6T5_Err750;
   struct var{
     double FdInfo_xMax;
     double FdInfo_eTot;
     double FdInfo_zenith;
     double FdInfo_azimuth;
     double FdInfo_xCore;
     double FdInfo_yCore;
     };

   struct var Eye_info;
   string timeStamp(SDEvent& theSdEvent);
   void  stationFunc(const std::vector<SdRecStation>& fStation, TProfile *hprof, double S250, double beta750);
  /* struct calib{
     };
 */struct var_SD{
     double S_opt750;
     double Theta750;
     double Energy750;
     double S_refAngle750;
     double beta750;
     double S_opt433;
     double Theta433;
     double Energy433;
     double S_refAngle433;
     double beta433;
     };

   struct varCher{
     double Energy433;
     double EnergyFD;
     };
   struct multVar{
     vector<double> energyVector;
     vector<double> zenithVector;
     };

   std::map<int,varCher> mymapCher;
   std::map<int,var_SD> mymap;
   std::map<char *,multVar> mapMult;
  bool pointInTriangle(double x1, double y1, double x2, double y2, double x3, double y3, double px, double py);
  void plotResidual(std::map<int,var_SD> mymap, char *name);
  struct limits{
	double max;
	double min;
        };
   std::vector<limits> AngularBinCalculator();
  
   string timeStamp(SDEvent& theSdEvent);

double computeBetaFromThetaAndSref(double *x, double *par);
double computeNKG(double dist, double sref, double beta);
double calculateMin(vector<double> en);
double calculateMax(vector<double> en);
vector<double> en_bin1, en_bin2, en_bin3, dist_vec;
double ShapeModel(const double cosTheta, const double showerSize);


std::vector<TGraphAsymmErrors*> myVectorGraphAsymmTrigSt;
std::vector<TProfile*> myVectorProfilesTrigSt;


//void plotResidual(vector <double>, vector <double>, char*);
