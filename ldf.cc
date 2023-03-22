// Calculate the time differences between the signals in the modules of the same counter

// using the leading edge of large signals

#include "TGraph2DErrors.h" 
#include "TProfile2D.h"
#include <iterator>
#include "TProfile.h"
#include "TFile.h"
#include <utl/ErrorLogger.h>
#include <utl/AugerException.h>
#include <utl/Point.h>
#include <utl/Vector.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AugerUnits.h>
#include <utl/Branch.h>
#include "ldf.h"
#include <fwk/CentralConfig.h>
#include "TGraph.h"
#include "TEfficiency.h"
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH3.h>
#include <evt/Event.h>
#include <evt/ShowerRecData.h>
#include <evt/ShowerSRecData.h>
#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <det/Detector.h>
#include <sdet/SDetector.h>
#include <sdet/Station.h>
#include <map>
#include <string.h>

#include "TBox.h"
#include <evt/Event.h>
#include <evt/ShowerRecData.h>
#include <evt/ShowerSRecData.h>

#include <sevt/SEvent.h>
#include <sevt/EventTrigger.h>
#include <sevt/Station.h>
#include <sevt/StationTriggerData.h>



#include <sstream>
#include <fstream>
#include <cassert>
#include <numeric>
#include <iomanip>

#include <RecEvent.h>
#include <DetectorGeometry.h>
#include <RecEventFile.h>
#include <math.h>
#include <TAttMarker.h>

#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TList.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include "TGaxis.h" 
#include <stdio.h>

#include "TH1.h"
#include "TTrace.h"
#include "rootStyle.h"
#include "utils.h"

// Maximum zenith angle of events selected 
//const double maxThetaDeg(45);


// Detector geometry defined as a global variable visible from everywhere in this code
DetectorGeometry* theGeometry;

using std::cout;
using std::endl;
unsigned long nEvents(0); 



/*========================================================

print some usage hints

==========================================================*/

void usage() {

cerr <<  "\n  Usage: AnalyseADST <ADSTfileName1> <ADSTfileName2> <Outputfilename>  \n"

<<  "         where <ADSTfileName> is the name of an ADST file, e.g. HybridRec.root \n"

<<  "         or a list of ADST files, e.g. \"HybridRec*.root\". \n "<<endl;

}



//##############################################################

int main(int argc, char **argv) 
{

/*if ( argc != 3 ) {
    usage();
    exit(1);
  }*/
 std::fstream outFile;
 outFile.open("output.dat", outFile.out);

 // At least 2 arguments are experected: the program and a list of ADST files
 if (argc < 2) {  
    std::cerr << "Usage: " << argv[0] << " [adst files]" << endl;
    return -1;
    }

 // output file names
 std::string inputFilename(argv[0]);
 std::string outputFilenameBase = inputFilename.substr(0, inputFilename.rfind("."));
  
 // output-definitions
// TFile* histoFile= new TFile(argv[2],"RECREATE");


 
// load adst file names    
 std::vector<std::string> dataFileNames;
 for (int i=1; i<argc; ++i) {
     dataFileNames.push_back(argv[i]);
     }
  
 /* Event */
 RecEvent* theRecEvent = new RecEvent();
 /* Detector */
 theGeometry = new DetectorGeometry();
 /* Info */
 FileInfo fileInfo;
 RecEventFile dataFile(dataFileNames);
 /* Get Event */
 dataFile.SetBuffers(&theRecEvent);
 /* Get Detector */
 dataFile.ReadDetectorGeometry(*theGeometry);
 /* Get Info */
 dataFile.ReadFileInfo(fileInfo); 

 // constructor to write events 
 RecEventFile writeFile_saturBetaFree("~/work/LDF/LDF_gaia/ADST_SaturatedFreebeta.root", RecEventFile::eWrite);
 RecEventFile writeFile_NosaturBetaFree("~/work/LDF/LDF_gaia/ADST_NoSaturatedFreebeta.root", RecEventFile::eWrite);
 // Get Event 
 writeFile_saturBetaFree.SetBuffers(&theRecEvent);
 writeFile_saturBetaFree.WriteDetectorGeometry(*theGeometry);
 writeFile_saturBetaFree.WriteFileInfo(fileInfo); 
 // Get Event 
 writeFile_NosaturBetaFree.SetBuffers(&theRecEvent);
 writeFile_NosaturBetaFree.WriteDetectorGeometry(*theGeometry);
 writeFile_NosaturBetaFree.WriteFileInfo(fileInfo); 




 /* References */ 
 SDEvent& theSdEvent = theRecEvent->GetSDEvent();
 
 book_histograms();
 myVectorGraphAsymmTrigSt           = createGraphAsymm(4);
 myVectorProfilesTrigSt             = createProfiles(4,16,18,20,"ProfTrig");
 myVectorGraphAsymm                  = createGraphAsymm(6);
 myVectorProfiles                    = createProfiles(10,0,800,15,"ProfS");
 myVectorProfBeta_EnBin              = createProfiles(4,1,1.4,10,"profBeta");
 myVectorProfBeta_EnBinNozenithCut   = createProfiles(4,1,1.4,10,"profBetaNoZenith");
 myVectorBetaS                       = createProfiles(10,1,1.4,10,"BetaS");
 myVectorHistosMultZenith            = createTH1F(8, 1, 1.4, 10, "hMultZenith");
 myVectorHistosMultEn                = createTH1F(8, 16, 18, 20, "hMultEn");
 myVectorHistos                      = createTH1F(4, 0.5, 19.5, 19, "hMult");
 myVectorRoptLogSTheta               = createTH1F(4, 200, 600, 20, "hRoptlogSTheta");//lo uso per ropt in bin di logS e theta
 myVectorRopt                        = createTH1F(3, 200, 600, 80, "hRopt");//lo uso per lo studio dell ropt con multiplocitá
 myVectorRopt_saturation             = createTH1F(3, 200, 600, 80, "hRoptSatur");//lo uso per studio della ropt e con saturati




 double thetaMaxDeg     = 0;
 double thetaMaxRad     = 0;
 int nHdEvents_750FD    = 0; 
 int nHdEvents_433FD    = 0; 
 int n250               = 0;
 int n450 	        = 0; 
 int n450_6T5 	        = 0; 
 int n433_750	        = 0; 
 int n433_750_bin0      = 0;
 int n433_750_bin1      = 0;
 int n433_750_bin2      = 0;
 int n433_750_bin3      = 0;
 int n433_750_bin4      = 0;
 int Ev6T5_433	        = 0;
 int nFdEvents          = 0;
 double minEn 	        = 18;
 double minEnCut        = 17;
 int nSelectedEvents 	= 0;
 int n433RecEvents	= 0;
 int N6T5Events		= 0;
 int nSaturatedEvents	= 0;
 int nSaturatedEvents6T5= 0;


TGraph *gEn_Err = new TGraph();
TGraph *gEn6T5_Err = new TGraph();

  double Northing = 6114576.09;
  double Easting  = 450992.16;
  double KT_Northing = 0;
  double KT_Easting  = 0;
  double Lety_Northing = 6115009.85 - Northing;
  double Lety_Easting  = 451002.08 - Easting;
  double Guili_Northing = 6114787.76 - Northing;  
  double Guili_Easting  = 451371.31 - Easting;
  double Pipi_Northing = 6114355.47 - Northing; 
  double Pipi_Easting  = 451376.50 - Easting;
  double Rosella_Northing = 6114128.70 - Northing;
  double Rosella_Easting  = 450997.78 - Easting;
  double Chichino_Northing = 6114356.58 - Northing;
  double Chichino_Easting  = 450622.16 - Easting;
  double Catherina_Northing = 6114794.00 - Northing;
  double Catherina_Easting  = 450620.00 - Easting;

std::vector<limits> vBounds(AngularBinCalculator());
  double Deg_433 = TMath::RadToDeg()*(theSdEvent.GetSdRecShower().GetZenith());


 // Read events - Loop over events
  while (dataFile.ReadNextEvent() == RecEventFile::eSuccess) 
        {    
        int eventId2 = theSdEvent.GetEventId();
//      std::string eventIdAuger = theRecEvent->GetAugerId();
        int eventIdAuger = theRecEvent->GetAugerIdNumber();
        std::string eventId = theRecEvent->GetEventId();
        std::string sizeLabel = theSdEvent.GetSdRecShower().GetShowerSizeLabel();
        const auto& sRecShower = theSdEvent.GetSdRecShower();
        const auto& stationVector = theSdEvent.GetStationVector();
        const auto& badStationVector = theSdEvent.GetBadStationVector();
 
       //std::string eventId(eventId.substr(6,12));
       //cout << "Processing event Auger=" << eventIdAuger <<" Sd" << eventId2<< " boh=" << eventId << endl;
        ++nEvents;
        

  bool trigger6T5 = theSdEvent.Is6T5();

     

//Infill 750 analysis -----------------
if(sizeLabel == "S450"){
  //Skip  events
  double thetaRad_750 = theSdEvent.GetSdRecShower().GetZenith();
  double thetaDeg_750 = TMath::RadToDeg()*(theSdEvent.GetSdRecShower().GetZenith());
  energyRec_750 = theSdEvent.GetSdRecShower().GetEnergy();  
  //----------- CUTs -----------
  if(thetaDeg_750>55) {
    continue;
    }
  if(!trigger6T5) continue;	  
  //if(log10(energyRec_750)<17.0)
  //  continue; 

double energy_err_750 = theSdEvent.GetSdRecShower().GetEnergyError();
const TVector3 &theShowerCoreRec_750 = theSdEvent.GetSdRecShower().GetCoreUTMCS();
double xCoreRec_750 = theShowerCoreRec_750.X(); //Easting
double yCoreRec_750 = theShowerCoreRec_750.Y(); //Northing
double rRec_750     = std::sqrt(xCoreRec_750*xCoreRec_750+yCoreRec_750*yCoreRec_750);
double dist_x_rel_750 = TMath::Abs(450992.16 - xCoreRec_750); //relativo a KT
double dist_y_rel_750 = TMath::Abs(6114576.09 - yCoreRec_750);
double Total_distRel_750 = std::sqrt(dist_x_rel_750*dist_x_rel_750+dist_y_rel_750*dist_y_rel_750); // in metri
    
/*if(!(pointInTriangle(KT_Northing, KT_Easting, Lety_Northing, Lety_Easting, Guili_Northing, Guili_Easting, dist_y_rel_750, dist_x_rel_750) || pointInTriangle(KT_Northing, KT_Easting, Guili_Northing, Guili_Easting, Pipi_Northing, Pipi_Easting, dist_y_rel_750, dist_x_rel_750) || pointInTriangle(KT_Northing, KT_Easting, Pipi_Northing, Pipi_Easting, Rosella_Northing, Rosella_Easting, dist_y_rel_750, dist_x_rel_750) || pointInTriangle(KT_Northing, KT_Easting, Rosella_Northing, Rosella_Easting, Chichino_Northing, Chichino_Easting, dist_y_rel_750, dist_x_rel_750) || pointInTriangle(KT_Northing, KT_Easting, Chichino_Northing, Chichino_Easting, Catherina_Northing, Catherina_Easting, dist_y_rel_750, dist_x_rel_750) || pointInTriangle(KT_Northing, KT_Easting, Catherina_Northing, Catherina_Easting, Lety_Northing, Lety_Easting, dist_y_rel_750, dist_x_rel_750)))  continue;
*/

// END CUTs


    if(trigger6T5 ==1) n450_6T5++;
    

    EnRec_750.push_back((energyRec_750 > 0)? log10(energyRec_750) : 0);
    hEnRec_750->Fill(log10(energyRec_750));


    const ESdRecLevel reclevel750 = theSdEvent.GetRecLevel();
    if(reclevel750>3){

      if(!mymap.count(eventId2)){  //Veo si existe la estructura para el evento, si existe la uso, sino la creo
        struct var_SD SD_info = { -1, -1, -1, -1, -1, -1, -1, -1};
        mymap.insert(std::make_pair(eventId2,SD_info) );
        }

    mymap[eventId2].S_opt750  = theSdEvent.GetSdRecShower().GetShowerSize(); 
    mymap[eventId2].Energy750 = theSdEvent.GetSdRecShower().GetEnergy();
    mymap[eventId2].Theta750  = TMath::RadToDeg()*(theSdEvent.GetSdRecShower().GetZenith()); 
    mymap[eventId2].S_refAngle750 = energyCalib(sizeLabel, theSdEvent.GetSdRecShower().GetShowerSize(), thetaRad_750);
   // mymap[eventId2].beta750   = ShapeModel(cos(thetaRad_750), theSdEvent.GetSdRecShower().GetShowerSize());
    }



  //  }// fine 6T5 
    	
       if(log10(energyRec_750)<minEn) minEn=log10(energyRec_750);
	




    if(energyRec_750>0){

          if(trigger6T5 ==1){
            En_Err_y_6T5_750.push_back(energy_err_750/energyRec_750);
            Total_distRel_x_6T5_750.push_back(Total_distRel_750);
            }
          En_Err_y_750.push_back(energy_err_750/energyRec_750);
          Total_distRel_x_750.push_back(Total_distRel_750);
          }


}// fine 750 analysis








//Infill 433 analysis -----------------
//if(sizeLabel == "S250"){
  //Skip  events
  double energyRec_433 = theSdEvent.GetSdRecShower().GetEnergy();  
  double thetaRad_433 = theSdEvent.GetSdRecShower().GetZenith();
  double thetaDeg_433 = TMath::RadToDeg()*(theSdEvent.GetSdRecShower().GetZenith());
  double beta = sRecShower.GetBeta(), betaError = sRecShower.GetBetaError();
  int multiplicity = theSdEvent.GetNumberOfCandidates();
  bool bin1 = false;
  bool bin2 = false;
  bool bin3 = false;
  double logEn = log10(energyRec_433); 
  //cout << "log en " << logEn << endl; 
  if (theSdEvent.IsSaturated())  nSaturatedEvents++;


        //Canditate stations
        int Nstations = 0;
        const std::vector<SdRecStation>& fStation = theSdEvent.GetStationVector();
        vector <double> Sb_vector_tot, Sb_vector_mu, Sb_vector_el, Sb_vector_S250;
        for(unsigned int iS=0;iS<fStation.size();iS++) {
           const SdRecStation& recStation= fStation[iS];
           int stationId = recStation.GetId();
           if(recStation.GetRejectionStatus() & eOffGrid)
            continue;
           if(recStation.GetRejectionStatus() & eNoTrigger)
            continue;
                Nstations++;

		}

          if(thetaRad_433>vBounds[0].min &&  thetaRad_433<vBounds[0].max) myVectorProfilesTrigSt[0]->Fill(log10(energyRec_433),Nstations);
          if(thetaRad_433>vBounds[1].min &&  thetaRad_433<vBounds[1].max) myVectorProfilesTrigSt[1]->Fill(log10(energyRec_433),Nstations);
          if(thetaRad_433>vBounds[2].min &&  thetaRad_433<vBounds[2].max) myVectorProfilesTrigSt[2]->Fill(log10(energyRec_433),Nstations);
          if(thetaRad_433>vBounds[3].min &&  thetaRad_433<vBounds[3].max) myVectorProfilesTrigSt[3]->Fill(log10(energyRec_433),Nstations);
	   hTheta->Fill(thetaDeg_433);
  	   double secTheta = 1/(cos(thetaRad_433));
  	   evtDistrib->Fill(log10(energyRec_433),secTheta);
  if(log10(energyRec_433)<=16.5) myVectorProfBeta_EnBinNozenithCut[0]->Fill(secTheta,beta);       
  if (log10(energyRec_433)>16.5 && log10(energyRec_433)<=17) myVectorProfBeta_EnBinNozenithCut[1]->Fill(secTheta,beta);
  if (log10(energyRec_433)>17 && log10(energyRec_433)<=17.5) myVectorProfBeta_EnBinNozenithCut[2]->Fill(secTheta,beta);
  if (log10(energyRec_433)>17.5 && log10(energyRec_433)<=18) myVectorProfBeta_EnBinNozenithCut[3]->Fill(secTheta,beta);

  // CUTS-------------------------------------------------------------------------------------------------------
  if (thetaDeg_433>45) continue;






/*  if (log10(energyRec_433)<16.6)// con questo e l`angolo < 45 ho full eff
    continue; 
*/
  if (trigger6T5!=1)
     continue;
	   hEnRec_433->Fill(log10(energyRec_433));
if(log10(energyRec_433)<=16.7) continue;    
N6T5Events++; 

double energy_err_433 = theSdEvent.GetSdRecShower().GetEnergyError();    
EnRec_433.push_back((energyRec_433 > 0)? log10(energyRec_433) : 0);
//hEnRec_433->Fill(log10(energyRec_433));


//free + fixed beta
  if (theSdEvent.IsSaturated()) {
     myVectorHistos[3]->Fill(multiplicity);
     } else {	
     myVectorHistos[2]->Fill(multiplicity);
     }
  double Ropt = sRecShower.GetROpt();
  if (multiplicity < 5 ) myVectorRopt[0]->Fill(Ropt);
  else if (multiplicity > 5 && multiplicity < 10 ) myVectorRopt[1]->Fill(Ropt);
  else if (multiplicity > 10 ) myVectorRopt[2]->Fill(Ropt);
cout << "qui 5  " << endl;
        
//if (multiplicity <=4) cout << beta << " " << ShapeModel(cos(thetaRad_433), theSdEvent.GetSdRecShower().GetShowerSize()) << endl;
 // if (betaError<=0.) //---------------------------------------------beta libero
  //   continue;
     myVectorRopt_saturation[2]->Fill(Ropt);     


/*
//Full efficiency zenith dependent cut
if ((thetaDeg_433>0 &&  thetaDeg_433<=24) && logEn >= 16.4) {
   //cout << " in bin 1 " << thetaDeg_433 << " " << logEn << endl;
   en_bin1.push_back(logEn);
   }
else if ((thetaDeg_433>24 && thetaDeg_433<=35) && logEn >= 16.5) {
   //cout << " in bin 2 " << thetaDeg_433 << " " << logEn << endl;
   en_bin2.push_back(logEn);
   }
else if ((thetaDeg_433>35 && thetaDeg_433<=45) && logEn >= 16.7) {
   //cout << " in bin 3 " << thetaDeg_433 << " " << logEn << endl;
   en_bin3.push_back(logEn);
   } 
else { continue; }
*/


//just free beta 
  if (theSdEvent.IsSaturated()) {
     myVectorRopt_saturation[0]->Fill(Ropt);     
     writeFile_saturBetaFree.WriteEvent();
     myVectorHistos[1]->Fill(multiplicity);
    if (multiplicity >=12 && multiplicity <=15 ) {
/*       if (!mapMult.count(multiplicity)){  //Veo si existe la estructura para el evento, si existe la uso, sino la creo
          struct multVar structMult >};
          mapMult.insert(std::make_pair(multiplicity,structMult));
          }
  */   mapMult["12_15"].energyVector.push_back(log10(energyRec_433));
       mapMult["12_15"].zenithVector.push_back(secTheta);
//     myVectorHistosMultEn[0]->Fill(log10(energyRec_433));
//     myVectorHistosMultZenith[0]->Fill(thetaDeg_433);
       } else if (multiplicity >=16 && multiplicity <=19) {
       mapMult["16_19"].energyVector.push_back(log10(energyRec_433));
       mapMult["16_19"].zenithVector.push_back(secTheta);
//       myVectorHistosMultEn[1]->Fill(log10(energyRec_433));
//       myVectorHistosMultZenith[1]->Fill(thetaDeg_433);
       }
    for (int i=12; i<=19; i++) {
        if (multiplicity == i) myVectorHistosMultEn[i-12]->Fill(log10(energyRec_433));
        if (multiplicity == i) myVectorHistosMultZenith[i-12]->Fill(secTheta);
        }
     nSaturatedEvents6T5++;
     }//-------------------------------------------------------------------- fine saturated

     writeFile_NosaturBetaFree.WriteEvent();
     if(theSdEvent.IsSaturated())
	continue;
	myVectorRopt_saturation[1]->Fill(Ropt);     


    // Energy (Rec)
//double energy_err_433 = theSdEvent.GetSdRecShower().GetEnergyError();    
const TVector3 &theShowerCoreRec_433 = theSdEvent.GetSdRecShower().GetCoreUTMCS();
double xCoreRec_433 = theShowerCoreRec_433.X();
double yCoreRec_433 = theShowerCoreRec_433.Y();
double rRec_433     = std::sqrt(xCoreRec_433*xCoreRec_433+yCoreRec_433*yCoreRec_433);
double dist_x_rel_433 = TMath::Abs(450992.16 - xCoreRec_433);
double dist_y_rel_433 = TMath::Abs(6114576.09 - yCoreRec_433);
double Total_distRel_433 = std::sqrt(dist_x_rel_433*dist_x_rel_433+dist_y_rel_433*dist_y_rel_433);
sec_vector.push_back(secTheta);
dist_vec.push_back(Total_distRel_433);
if (Total_distRel_433>800) {
   cout << "rejecting event too far from KT " << Total_distRel_433 << endl;
   continue;
   }





/*
if(!(pointInTriangle(KT_Northing, KT_Easting, Lety_Northing, Lety_Easting, Guili_Northing, Guili_Easting, dist_y_rel_433, dist_x_rel_433) || pointInTriangle(KT_Northing, KT_Easting, Guili_Northing, Guili_Easting, Pipi_Northing, Pipi_Easting, dist_y_rel_433, dist_x_rel_433) || pointInTriangle(KT_Northing, KT_Easting, Pipi_Northing, Pipi_Easting, Rosella_Northing, Rosella_Easting, dist_y_rel_433, dist_x_rel_433) || pointInTriangle(KT_Northing, KT_Easting, Rosella_Northing, Rosella_Easting, Chichino_Northing, Chichino_Easting, dist_y_rel_433, dist_x_rel_433) || pointInTriangle(KT_Northing, KT_Easting, Chichino_Northing, Chichino_Easting, Catherina_Northing, Catherina_Easting, dist_y_rel_433, dist_x_rel_433) || pointInTriangle(KT_Northing, KT_Easting, Catherina_Northing, Catherina_Easting, Lety_Northing, Lety_Easting, dist_y_rel_433, dist_x_rel_433))) {continue; }

*/
// END CUTs
gEn_Err->SetPoint(n250,Total_distRel_433, energy_err_433/energyRec_433);
if(trigger6T5 ==1) gEn6T5_Err->SetPoint(n250,Total_distRel_433, energy_err_433/energyRec_433);



    //EnRec_433.push_back(energyRec_433/1e18);
/*EnRec_433.push_back((energyRec_433 > 0)? log10(energyRec_433) : 0);

    hEnRec_433->Fill(log10(energyRec_433));
  */  string timeEvt = timeStamp(theSdEvent);    
   

     const ESdRecLevel reclevel433 = theSdEvent.GetRecLevel();
     if(reclevel433>3){

	n433RecEvents++;

      if(!mymap.count(eventId2)){  //Veo si existe la estructura para el evento, si existe la uso, sino la creo
        struct var_SD SD_info = { -1, -1, -1, -1, -1, -1, -1, -1, -1};
        mymap.insert(std::make_pair(eventId2,SD_info) );
        }
      mymap[eventId2].S_opt433  = theSdEvent.GetSdRecShower().GetShowerSize(); 
      mymap[eventId2].Energy433 = theSdEvent.GetSdRecShower().GetEnergy();
      mymap[eventId2].Theta433  = theSdEvent.GetSdRecShower().GetZenith(); 
     // mymap[eventId2].Theta433  = TMath::RadToDeg()*(theSdEvent.GetSdRecShower().GetZenith()); 
      mymap[eventId2].S_refAngle433  = energyCalib(sizeLabel, theSdEvent.GetSdRecShower().GetShowerSize(), thetaRad_433); 
      mymap[eventId2].beta433   = ShapeModel(cos(thetaRad_433), theSdEvent.GetSdRecShower().GetShowerSize());

      }

      if(!mymapCher.count(eventIdAuger)){  //Veo si existe la estructura para el evento, si existe la uso, sino la creo
        struct varCher Eninfo = { -1, -1};
        mymapCher.insert(std::make_pair(eventIdAuger,Eninfo) );
        }
      mymapCher[eventIdAuger].Energy433 = theSdEvent.GetSdRecShower().GetEnergy();
      

     if(mymap[eventId2].S_opt433 > 0 && mymap[eventId2].S_opt750 > 0) {
       nSelectedEvents++;
//     writeFile_SelectedEvents.WriteEvent();
       }
/*
    if(energyRec_433>0){

        if(trigger6T5 ==1){
        En_Err_y_6T5_433.push_back(energy_err_433/energyRec_433);
        Total_distRel_x_6T5_433.push_back(Total_distRel_433);
        Ev6T5_433++;        }
        En_Err_y_433.push_back(energy_err_433/energyRec_433);
        Total_distRel_x_433.push_back(Total_distRel_433);
}*/


   double S250         = theSdEvent.GetSdRecShower().GetShowerSize();
   double S250_err      = theSdEvent.GetSdRecShower().GetShowerSizeError();
   //const TVector3 &theShowerCoreRec = theSdEventSi.GetSdRecShower().GetCoreSiteCS();
   double xCoreRec_err_433  = theSdEvent.GetSdRecShower().GetCoreNorthingError();
   double yCoreRec_err_433  = theSdEvent.GetSdRecShower().GetCoreEastingError();
   double energyRec_err_433 = theSdEvent.GetSdRecShower().GetEnergyError();

if((thetaDeg_433>=29 && thetaDeg_433<=31) && (S250>=180 && S250<=220)) cout << "per Nico energia " << S250 << " " << thetaDeg_433 << " " << energyRec_433 << endl;

   //Propagazione degli errori
   double rRec_err_433 = std::sqrt(((xCoreRec_433*xCoreRec_433)/(rRec_433))*(xCoreRec_err_433*xCoreRec_err_433) + ((yCoreRec_433*yCoreRec_433)/(rRec_433))*(yCoreRec_err_433*yCoreRec_err_433));



    double deltaR_433 = rRec_433/rRec_err_433;
    double deltaS250 = S250/S250_err;
    double deltaEn_433 = energyRec_433/energyRec_err_433;

    deltaEn_v_433.push_back(deltaEn_433);
    deltaS250_v.push_back(deltaS250);
    deltaR_v_433.push_back(deltaR_433);
    enRec_v_433.push_back((energyRec_433 > 0)? log10(energyRec_433) : 0);
//    const std::vector<SdRecStation>& fStation = theSdEvent.GetStationVector();
    int nHighSaturatedSt = 0, nLowSaturatedSt = 0;
    const  SdRecStation* fHottestStation = &fStation[0];
    double maxSignal = 0;
    hRopt->Fill(Ropt);

   /* for (int i=3; i<=10; i++) {
        if (multiplicity == i) myVectorRopt[i-3]->Fill(Ropt);
        }
*/
//Ropt per differenti Shower size e bin di theta
if ((thetaDeg_433>0 &&  thetaDeg_433<=30) && (log10(S250) >= 1 && log10(S250)<=1.5)) myVectorRoptLogSTheta[0]->Fill(Ropt);
if ((thetaDeg_433>30 &&  thetaDeg_433<=45) && (log10(S250) >= 1 && log10(S250)<=1.5)) myVectorRoptLogSTheta[1]->Fill(Ropt);
if ((thetaDeg_433>0 &&  thetaDeg_433<=30) && (log10(S250) > 1.5 && log10(S250)<=2)) myVectorRoptLogSTheta[2]->Fill(Ropt);
if ((thetaDeg_433>30 &&  thetaDeg_433<=45) && (log10(S250) > 1.5 && log10(S250)<=2)) myVectorRoptLogSTheta[3]->Fill(Ropt);
           int counter=0;
		 for (std::vector<SdRecStation>::const_iterator sIt=fStation.begin(); sIt!=fStation.end(); ++sIt) {

                   double totalSignal = sIt->GetTotalSignal();
                   double totalSignalErr = sIt->GetTotalSignalError();
                   double rSP = sIt->GetSPDistance();
		   int StationId = sIt->GetId();
		   ldfPura->SetPoint(counter,rSP,totalSignal);
		   ldfPura->SetPointError(counter,0, totalSignalErr);
	           counter++;
cout << " totalSignal" <<totalSignal<< " rSP" <<rSP <<  " beta" << beta << " S250" <<S250 <<endl;
cout << " computeNKG" <<computeNKG(rSP,S250, beta) << endl;

		   ldfModel->SetPoint(counter,rSP,computeNKG(rSP,S250, beta));
 		   

			if(maxSignal < totalSignal) {
			  maxSignal = totalSignal;
			    fHottestStation = &(*sIt);
		  }


double sPDistance = fHottestStation->GetSPDistance();
//cout << "sPDistance" << sPDistance << " " << fHottestStation->GetId() << endl; 
hprofRopt.Fill(sPDistance, Ropt);
ldf2DRopt->Fill(sPDistance, Ropt);
	}

/*
           //commentato        if(log(energyRec)<=16.5) {cout <<"into if"<< endl; hprofEn0.Fill(rSP,ratio);}
                   if(log(energyRec_433)<=16.5) {cout <<"into if"<< endl; hprofEn0.Fill(rSP,ratio);}
                   if(log(energyRec_433)>16.5 && log(energyRec_433)<=17) {hprofEn1.Fill(rSP,ratio*1.5);}
                   if(log(energyRec_433)>17 && log(energyRec_433)<=17.5) {hprofEn2.Fill(rSP,ratio*2.5);}
                   if(log(energyRec_433)>17.5 && log(energyRec_433)<=18) {hprofEn3.Fill(rSP,ratio*3.5);}
//fino a qui
                   if(energyRec_433<=pow(10,16.5)) { hprofEn0.Fill(rSP,ratio);}
                   if(energyRec_433>pow(10,16.5) && energyRec_433<=pow(10,17)) {hprofEn1.Fill(rSP,ratio);}
                   if(energyRec_433>pow(10,17) && energyRec_433<=pow(10,17.5)) {hprofEn2.Fill(rSP,ratio);}
                   if(energyRec_433>pow(10,17.5) && energyRec_433<=pow(10,18)) {hprofEn3.Fill(rSP,ratio);}


                   if (thetaRad_433>vBounds[0].min &&  thetaRad_433<vBounds[0].max){
                      double signal0 = sIt->GetTotalSignal();
                      signal0_v.push_back(signal0/S250);
                      double rSP0 = sIt->GetSPDistance();
                      hprof0.Fill(rSP,ratio);
                      rSP0_v.push_back(rSP0);
                      }
                   if (thetaRad_433>vBounds[1].min &&  thetaRad_433<vBounds[1].max){
                      double signal1 = sIt->GetTotalSignal();
                      signal1_v.push_back(signal1/S250);
                      double rSP1 = sIt->GetSPDistance();
                      hprof1.Fill(rSP,ratio*1.5);
                      rSP1_v.push_back(rSP1);
                      }
                   if (thetaRad_433>vBounds[2].min &&  thetaRad_433<vBounds[2].max){
                      double signal2 = sIt->GetTotalSignal();
                      signal2_v.push_back(signal2/S250);
                      double rSP2 = sIt->GetSPDistance();
                      hprof2.Fill(rSP,ratio*2.5);
                      rSP2_v.push_back(rSP2);
                      }
                   if (thetaRad_433>vBounds[3].min &&  thetaRad_433<vBounds[3].max){
                      double signal3 = sIt->GetTotalSignal();
                      signal3_v.push_back(signal3/S250);
                      double rSP3 = sIt->GetSPDistance();
                      hprof3.Fill(rSP,ratio);
                      rSP3_v.push_back(rSP3);
                      }
*/
// beta study
cout << "qui 2  " << endl;

     TBits zero = theRecEvent->GetDetector().GetActiveStations();
     beta_vector.push_back(beta);
     ldf3D->Fill(secTheta,log10(S250), beta);
     ldf2D->Fill(secTheta,log10(S250));
     ldf2D_S250->Fill(log10(S250), beta);
hprof2d.Fill(secTheta,log10(S250),beta);
profBetaError.Fill(betaError, log10(S250)); 
//cout << beta << " " << betaError << " "<< log10(S250) << " " << log10(S250_err)<< " " << secTheta << endl;
vector_beta.push_back(beta);
vectorErr_beta.push_back(betaError);
vector_logS250.push_back(log10(S250));
vectorErr_logS250.push_back(log10(S250_err));
vector_secTheta.push_back(secTheta);

	   if(log10(energyRec_433)<=16.5) myVectorProfBeta_EnBin[0]->Fill(secTheta,beta);       
	   if (log10(energyRec_433)>16.5 && log10(energyRec_433)<=17) myVectorProfBeta_EnBin[1]->Fill(secTheta,beta);
	   if (log10(energyRec_433)>17 && log10(energyRec_433)<=17.5) myVectorProfBeta_EnBin[2]->Fill(secTheta,beta);
	   if (log10(energyRec_433)>17.5 && log10(energyRec_433)<=18) myVectorProfBeta_EnBin[3]->Fill(secTheta,beta);
//cout << " beta " << beta << " beta433 " << mymap[eventId2].beta433 << endl;       
if (log10(S250)<1.0)  {

    stationFunc(fStation, myVectorProfiles[0], S250, mymap[eventId2].beta433);
    myVectorBetaS[0]->Fill(secTheta,beta);
    }

if (log10(S250)>=1.0 && log10(S250)<1.1)   {
    stationFunc(fStation, myVectorProfiles[1], S250, mymap[eventId2].beta433);
    myVectorBetaS[1]->Fill(secTheta,beta);
    }

if (log10(S250)>=1.1 && log10(S250)<1.2)  {
    stationFunc(fStation, myVectorProfiles[2], S250, mymap[eventId2].beta433);
    myVectorBetaS[2]->Fill(secTheta,beta);
    }	

if (log10(S250)>=1.2 && log10(S250)<1.3)  {
   stationFunc(fStation, myVectorProfiles[3], S250, mymap[eventId2].beta433);
   myVectorBetaS[3]->Fill(secTheta,beta);
    }

if (log10(S250)>=1.3 && log10(S250)<1.4)  {
   stationFunc(fStation, myVectorProfiles[4], S250, mymap[eventId2].beta433);
   myVectorBetaS[4]->Fill(secTheta,beta);
   }
if (log10(S250)>=1.4 && log10(S250)<1.5)  {
   stationFunc(fStation, myVectorProfiles[5], S250, mymap[eventId2].beta433);
   myVectorBetaS[5]->Fill(secTheta,beta);
   }
if (log10(S250)>=1.5 && log10(S250)<1.6)  {
   stationFunc(fStation, myVectorProfiles[6], S250, mymap[eventId2].beta433);
   myVectorBetaS[6]->Fill(secTheta,beta);
   }
if (log10(S250)>=1.6 && log10(S250)<1.7)  {
   stationFunc(fStation, myVectorProfiles[7], S250, mymap[eventId2].beta433);
   myVectorBetaS[7]->Fill(secTheta,beta);
   }
if (log10(S250)>=1.7 && log10(S250)<1.8)  {
   stationFunc(fStation, myVectorProfiles[8], S250, mymap[eventId2].beta433);
   myVectorBetaS[8]->Fill(secTheta,beta);
   }

if (log10(S250)>=1.8)  {
   stationFunc(fStation, myVectorProfiles[9], S250, mymap[eventId2].beta433);
   myVectorBetaS[9]->Fill(secTheta,beta);
   }
/* 
if (log10(S250)>=1. && log10(S250)<1.2)  cout << " interval 1  " << log10(energyRec_433) << endl;
if (log10(S250)>=1.2 && log10(S250)<1.4)  cout << " interval 2  " << log10(energyRec_433) << endl;
if (log10(S250)>=1.4 && log10(S250)<1.6)  cout << " interval 3  " << log10(energyRec_433) << endl;
if (log10(S250)>=1.6 && log10(S250)<1.8)  cout << " interval 4  " << log10(energyRec_433) << endl;
if (log10(S250)>=1.8 && log10(S250)<2.0)  cout << " interval 5  " << log10(energyRec_433) << endl;
if (log10(S250)>=2.0 && log10(S250)<2.2)  cout << " interval 6  " << log10(energyRec_433) << endl;
if (log10(S250)>=2.2 && log10(S250)<2.4)  cout << " interval 7  " << log10(energyRec_433) << endl;
if (log10(S250)>=2.4 && log10(S250)<2.6)  cout << " interval 8  " << log10(energyRec_433) << endl;
*/




// prepare silent extraction
for (const auto& station : stationVector) {
zero[station.GetId()] = false;
}   
for (const auto& badStation : badStationVector){
zero[badStation.GetId()] = false;
}
int nSilent = 0;
for (unsigned int id = 0, n = zero.GetNbits(); id < n; ++id) {
if(zero[id]){
  nSilent++;
}
}
//int multiplicity = theSdEvent.GetNumberOfCandidates() + nSilent;
//hMult->Fill(multiplicity);	
myVectorHistos[0]->Fill(multiplicity);
cout << "qui 3  " << endl;





//} // fine size label 250 fine 433 analysis


/*
//------------- Ojos -------------------
/// Creo que jEye = 6 es HeCO (hay que revisarlo)
for(unsigned int jEye = 1; jEye <= 6; ++jEye) {
//     if(jEye != 6)
//         continue;
if(!FdCuts(theRecEvent, jEye))
continue;
//dalla funzione FdCuts l´evento esce gia con un reclevel10


FdInfoEye(theRecEvent, jEye, outFile, &Eye_info);

//!!!!!!!!!!!!!!-------  FD Part -----------!!!!!!!!!!!!!!
//-------- XMax --------
//commentato      double XmaxMC = theGenShower.GetXmax();
double DeltaXmax = Eye_info.FdInfo_xMax - XmaxMC;
//cout <<"DeltaXMax= " << DeltaXmax << endl;
h8->Fill(DeltaXmax);
// fine commento

//-------- Energia --------
//double DeltaE_Eye = 100*(Eye_info.FdInfo_eTot - energyMC)/energyMC;
double En_Fd = Eye_info.FdInfo_eTot;
hEn_Fd->Fill(log10(En_Fd));
EnRec_FD.push_back((En_Fd > 0)? log10(En_Fd) : 0);
//EnRec_FD.push_back(En_Fd/1e18);
nFdEvents++;
//cout << "eventIdAuger number " << eventIdAuger << endl;

//commentato  if(!mymapCher.count(eventId)){  //Veo si existe la estructura para el evento, si existe la uso, sino la creo
   struct varCher Eninfo = { -1, -1};
   mymap.insert(std::make_pair(eventId,Eninfo) );
   }
mymapCher[eventId2].EnergyFD = En_Fd;
}

//-------- Theta --------
double DeltaTheta_Eye = Eye_info.FdInfo_zenith - zenithSim;
h9->Fill((DeltaTheta_Eye*180)/TMath::Pi());
double DeltaTheta_Eye = Eye_info.FdInfo_zenith - zenithSim;
h9->Fill((DeltaTheta_Eye*180)/TMath::Pi());
double DeltaTheta_Eye2 = (DeltaTheta_Eye*180)/TMath::Pi();

//-------- Core position --------
double d_xEye = Eye_info.FdInfo_xCore - xCoreSim;
double d_yEye = Eye_info.FdInfo_yCore - yCoreSim;
double Delta_r = std::sqrt(d_xEye*d_xEye + d_yEye*d_yEye);
h10->Fill(Delta_r);

//-------- Alpha --------
double CosDeltaAlpha_Eye = (sin(Eye_info.FdInfo_zenith)*cos(Eye_info.FdInfo_azimuth))*(sin(zenithSim)*cos(azimutSim))+(sin(Eye_info.FdInfo_zenith)*sin(Eye_info.FdInfo_azimuth))*(sin(zenithSim)*sin(azimutSim))+(cos(Eye_info.FdInfo_zenith)*cos(zenithSim));
double DeltaAlpha_Eye = acos(CosDeltaAlpha_Eye);
double DeltaAlpha_Eye2 = (DeltaAlpha_Eye*180)/TMath::Pi();
h11->Fill((DeltaAlpha_Eye*180)/TMath::Pi());
//fine commento

} /// End for(jEye)
*/





}//end of while
/*for (int i=0; i<4; i++){
        binomialErrors(*(myVectorProfilesTrigSt[i]),myVectorGraphAsymmTrigSt[i]);
}

graphFunc(sec_vector, beta_vector, "beta");
graphErrorFunc(vector_secTheta,vector_logS250, vector_beta, vectorErr_logS250, vectorErr_beta, "graphErrLdf2D");
for (std::map<char *,multVar>::iterator it=mapMult.begin(); it!=mapMult.end(); ++it) {
//graphFunc(it->second.energyVector, it->second.zenithVector, strcat("graphMult_",to_string(it->first).c_str()) );     
//char *name = const_cast<char*>(("graphMult_" + to_string(it->first)).c_str());
graphFunc(it->second.energyVector, it->second.zenithVector, it->first);
}
*/
cout << "qui 4  " << endl;


cout << " Max ref distance from KT " << calculateMax(dist_vec) << " ";
cout << " min ref distance from KT " << calculateMin(dist_vec) << " ";
cout << " Max beta " << calculateMax(vector_beta) << " ";
cout << " min beta " << calculateMin(vector_beta) << " ";
cout << " Max logS250 " << calculateMax(vector_logS250) << " ";
cout << " min logS250 " << calculateMin(vector_logS250) << " ";
cout << " min bin 1 0-24 grad " << calculateMin(en_bin1) << " ";
cout << " min bin 2 24-35 grad " << calculateMin(en_bin2)<< " ";
cout << " min bin 3 35-45 grad " << calculateMin(en_bin3) << endl;

/*
graphFunc(rSP_v, signal_v, "LDF");
graphFunc(rSP0_v, signal0_v, "LDF0");
graphFunc(rSP1_v, signal1_v, "LDF1");
graphFunc(rSP2_v, signal2_v, "LDF2");
graphFunc(rSP3_v, signal3_v, "LDF3");


AngularBinCalculator();
cout << "vBounds 0 " << "min= "<< vBounds[0].min* 180/TMath::Pi()<< "----- max  " << vBounds[0].max* 180/TMath::Pi()<< endl;
cout << "vBounds 1 " << "min= "<< vBounds[1].min* 180/TMath::Pi()<< "----- max  " << vBounds[1].max* 180/TMath::Pi()<< endl;
cout << "vBounds 2 " << "min= "<< vBounds[2].min* 180/TMath::Pi()<< "----- max  " << vBounds[2].max* 180/TMath::Pi()<< endl;
cout << "vBounds 3 " << "min= "<< vBounds[3].min* 180/TMath::Pi()<< "----- max  " << vBounds[3].max* 180/TMath::Pi()<< endl;

for(std::map<int,var_SD>::iterator it=mymap.begin(); it!=mymap.end(); ++it) {
// it->first tiene el key del map
// it->second tiene el valor del mapi, en mi caso la struct
// Hago algo si tengo los dos valores
if(it->second.S_opt433 > 0 && it->second.S_opt750 > 0) {
n433_750++;
s250_y.push_back(it->second.S_opt433);
s450_y.push_back(it->second.S_opt750);

s250_tot.push_back(it->second.S_opt433);
s450_tot.push_back(it->second.S_opt750);
sRef35_tot.push_back(it->second.S_refAngle750);
sRef33_tot.push_back(it->second.S_refAngle433);

if(it->second.Theta433>vBounds[0].min &&  it->second.Theta433<vBounds[0].max){       
s250_bin0.push_back(it->second.S_opt433);
s450_bin0.push_back(it->second.S_opt750);
sRef33_bin0.push_back(it->second.S_refAngle433);
sRef35_bin0.push_back(it->second.S_refAngle750);
n433_750_bin0++;
zenithVector0.push_back(TMath::RadToDeg()*(it->second.Theta433));
}
if(it->second.Theta433>vBounds[1].min &&  it->second.Theta433<vBounds[1].max){       
s250_bin1.push_back(it->second.S_opt433);
s450_bin1.push_back(it->second.S_opt750);
sRef33_bin1.push_back(it->second.S_refAngle433);
sRef35_bin1.push_back(it->second.S_refAngle750);
n433_750_bin1++;
zenithVector1.push_back(TMath::RadToDeg()*(it->second.Theta433));
}
if(it->second.Theta433>vBounds[2].min &&  it->second.Theta433<vBounds[2].max){       
s250_bin2.push_back(it->second.S_opt433);
s450_bin2.push_back(it->second.S_opt750);
sRef33_bin2.push_back(it->second.S_refAngle433);
sRef35_bin2.push_back(it->second.S_refAngle750);
n433_750_bin2++;
zenithVector2.push_back(TMath::RadToDeg()*(it->second.Theta433));
}
if(it->second.Theta433>vBounds[3].min &&  it->second.Theta433<vBounds[3].max){       
s250_bin3.push_back(it->second.S_opt433);
s450_bin3.push_back(it->second.S_opt750);
sRef33_bin3.push_back(it->second.S_refAngle433);
sRef35_bin3.push_back(it->second.S_refAngle750);
n433_750_bin3++;
zenithVector3.push_back(TMath::RadToDeg()*(it->second.Theta433));
}

}
} 
*/




//Prova per vedere se il vettore viene riempito correttamente
/*int N = s450_y.size();
double x[N];

for (int i=0 ; i<N ; i++)
    {
     x[i]=s450_y[i];
     cout << "Sopt750" << ": x[" << x[i] << "]" << endl;
    }
*/
//intervalli di sin2theta
/*double maxTheta0 = zenithBin(zenithVector0);
double maxTheta1 = zenithBin(zenithVector1);
double maxTheta2 = zenithBin(zenithVector2);
double maxTheta3 = zenithBin(zenithVector3);
//double maxTheta4 = zenithBin(zenithVector4);
cout << "bin0 maxTheta é "<< maxTheta0 << endl;
cout << "bin1 maxTheta é "<< maxTheta1 << endl;
cout << "bin2 maxTheta é "<< maxTheta2 << endl;
cout << "bin3 maxTheta é "<< maxTheta3 << endl;
//cout << "bin4 maxTheta é "<< maxTheta4 << endl;
cout << "fd events= " << nFdEvents << endl;
cout << nSelectedEvents << " nSelectedEventsHd433-750" << endl;
cout << n433_750 << " Hd 433 && 750 events" << endl;
cout << n433_750_bin0 << " Hd 433 && 750 events bin0 theta" << endl;
cout << n433_750_bin1 << " Hd 433 && 750 events bin1 theta" << endl;
cout << n433_750_bin2 << " Hd 433 && 750 events bin2 theta" << endl;
cout << n433_750_bin3 << " Hd 433 && 750 events bin3 theta" << endl;
//cout << n433_750_bin4 << " Hd 433 && 750 events bin4 theta" << endl;
*/
cout << nEvents << " total events" << nSaturatedEvents << " Saturated Events " << endl;
cout << N6T5Events << "selected events over " << n433RecEvents << " rec 433 events" << endl;

//---------------Istogrammi e grafici-------------------------------------------
// graphTheta(EnRec_750, s450_y,  "Extimator450");
// graphFunc(EnRec_433, EnRec_750, "Energies433_750"); 
// graphTheta(EnRec_433, s250_y,  "Extimator250");
/*graphFunc(EnRec_433, EnRec_FD,  "EnergySD433_FD");
graphFunc(s250_y, s450_y,       "Extimators");
graphFunc(s250_tot, s450_tot,       "Extimators_tot");
graphFunc(s250_bin0, s450_bin0,     "Extimators_bin0"); 
graphFunc(s250_bin1, s450_bin1,     "Extimators_bin1"); 
graphFunc(s250_bin2, s450_bin2,     "Extimators_bin2"); 
graphFunc(s250_bin3, s450_bin3,     "Extimators_bin3"); 
graphFunc(sRef33_tot, sRef35_tot,   "ExtimatorsRef_tot");
graphFunc(sRef33_bin0, sRef35_bin0,   "ExtimatorsRef_bin0");
graphFunc(sRef33_bin1, sRef35_bin1,   "ExtimatorsRef_bin1");
graphFunc(sRef33_bin2, sRef35_bin2,   "ExtimatorsRef_bin2");
graphFunc(sRef33_bin3, sRef35_bin3,   "ExtimatorsRef_bin3");
//plotResidual(sRef33_tot, sRef35_tot, "ResidualS_Ref");
plotResidual(mymap, "ResidualS_Ref");

gEn_Err->Write("EnergyErr_433");
gEn6T5_Err->Write("EnergyErr6T5_433");
*/

write_histos(myVectorProfiles, myVectorBetaS, myVectorHistos, myVectorRopt, myVectorRopt_saturation, myVectorProfBeta_EnBin, myVectorProfBeta_EnBinNozenithCut, myVectorHistosMultEn, myVectorHistosMultZenith, myVectorRoptLogSTheta, myVectorGraphAsymmTrigSt);
outFile.close();
file->Close();

return 0;



} //fine main









//---------------Funzioni------------------------------------------- 
void graphErrorFunc(vector <double> var_x, vector <double> var_y, vector <double> var_z, vector <double> err_y, vector <double> err_z, char *name){


int N = 0;
if(var_x.size() < var_y.size())
{N = var_x.size();
}else{
N = var_y.size();}
if(var_z.size() < N) N = var_z.size();
double x[N];
double y[N];
double z[N];
double ex[N];
double ey[N];
double ez[N];
if(N<=0) return;
for (int i=0 ; i<N ; i++)
    {	
     x[i]=var_x[i];
     y[i]=var_y[i];
     z[i]=var_z[i];
     ex[i]=0;
     ey[i]=err_y[i];
     ez[i]=err_z[i];
//             cout << "graphErrorFunc  N:" << i << ": x[" << x[i] << "], y[" << y[i] << "]" << endl;
//cout << var_x.size() << " = size s250 vector" <<" ----" << var_y.size() << " = size s450 vector" << " " << name << endl; 
    }
TGraph2DErrors *f = new TGraph2DErrors(N, x, y, z, ex, ey, ez);
//f->GetXaxis()->SetRangeUser(16.0,18.5);
f->Write(name);
}

void graphFunc(vector <double> var_x, vector <double> var_y, char *name){



int N = 0;
if(var_x.size() > var_y.size())
{N = var_y.size();
}else{
N = var_x.size();}


double x[N];
double y[N];

for (int i=0 ; i<N ; i++)
    {	
     x[i]=var_x[i];
     y[i]=var_y[i];
     //cout << "graphFunc N:" << i << ": x[" << x[i] << "], y[" << y[i] << "]" << endl;
    }
//cout << var_x.size() << " = size s250 vector" <<" ----" << var_y.size() << " = size s450 vector" << " " << name << endl; 
TGraph *f = new TGraph(N,x,y);
//f->GetXaxis()->SetRangeUser(16.0,18.5);
f->Write(name);
}


void write_histos(std::vector<TProfile*> prof_vector, std::vector<TProfile*> prof_vector2, std::vector<TH1F*> h_vector, std::vector<TH1F*> vec_Ropt, std::vector<TH1F*> vec_RoptSatur, std::vector<TProfile*> prof_vector3, std::vector<TProfile*> prof_vector4, std::vector<TH1F*> vectorMultEn, std::vector<TH1F*> vectorMultZenith, std::vector<TH1F*> vec_RoptLogSTheta, std::vector<TGraphAsymmErrors*> myVectorGraphAsymmTrigSt) {

for (int i=0; i<4; i++){
        std::string nameTrigSt = "TrigSt" + to_string(i);
        myVectorGraphAsymmTrigSt[i]->Write(nameTrigSt.c_str());
}
hEnRec_433->Write();
hTheta->Write();
hRopt->Write();
hEnRec_750->Write();
hEn_Fd->Write();
pTriggerET4.Write();
pTriggerETot.Write();
hprofRopt.Write();
ldfPura->Write("ldfPura");
ldfModel->Write("ldfModel");
hprof0.Write();
hprof1.Write();
hprof2.Write();
hprof3.Write();
hprof2d.Write();
//profBetaError.Write();
//	evtDistrib->Write();
	ldf3D->Write();
	ldf2D->Write();
	ldf2D_S250->Write();
	ldf2DRopt->Write();
	for (int i=0; i<h_vector.size(); i++) {
	      h_vector[i]->Write();
	      }
	for (int i=0; i<prof_vector.size(); i++) {
	      prof_vector[i]->Write();
	      }
	for (int i=0; i<prof_vector2.size(); i++) {
	     prof_vector2[i]->Write();
	     }
	for (int i=0; i<vec_Ropt.size(); i++) {
	     vec_Ropt[i]->Write();
	     }
	for (int i=0; i<vec_RoptLogSTheta.size(); i++) {
	     vec_RoptLogSTheta[i]->Write();
	     }
	for (int i=0; i<vec_RoptSatur.size(); i++) {
	     vec_RoptSatur[i]->Write();
	     }
	for (int i=0; i<prof_vector3.size(); i++) {
	      prof_vector3[i]->Write();
	      }
	for (int i=0; i<prof_vector4.size(); i++) {
	      prof_vector4[i]->Write();
	      }
	for (int i=0; i<vectorMultEn.size(); i++) {
	      vectorMultEn[i]->Write();
	      }
	for (int i=0; i<vectorMultZenith.size(); i++) {
	      vectorMultZenith[i]->Write();
	      }

	}

	std::vector<TGraphAsymmErrors*> createGraphAsymm (int n) {
	  std::vector<TGraphAsymmErrors*> graph_vector;
	  for (int i=0; i<n; i++) {
	      TGraphAsymmErrors *gAsym = new TGraphAsymmErrors();
	      graph_vector.push_back(gAsym);
	      }
	return graph_vector;

	}

	std::vector<TProfile*> createProfiles (int n, double inf, double sup, int bin, std::string  name) {
	  
	  std::vector<TProfile*> prof_vector;
	  for (int i=0; i<n; i++) {
	      std::string nameNew = name + to_string(i);
	      TProfile *profS = new TProfile(nameNew.c_str(),"Profile",bin,inf,sup);
	      prof_vector.push_back(profS);
	      }
	return prof_vector;

	}


	std::vector<TH1F*> createTH1F (int n, double inf, double sup, int bin, std::string  name) {
	  std::vector<TH1F*> th1_vector;
	  for (int i=0; i<n; i++) {
	      std::string nameNew = name + to_string(i);
	      TH1F *h = new TH1F(nameNew.c_str(), nameNew.c_str(), bin, inf, sup);
	      th1_vector.push_back(h);
	      }
	return th1_vector;

	}





	void book_histograms() {

	  file          = new TFile("histogramsLDFRopt250FixedBeta2013_201820Bin.root","RECREATE");
	  hRopt         = new TH1F("Histo_Ropt", "", 50, 100, 600);
	  hEn_Fd        = new TH1F("Histo_EnergyFd", "",     30, 16.0, 18.0);
	  hEnRec_433    = new TH1F("Histo_EnergyRec433", "", 20, 16.0, 18.0);
	  hTheta        = new TH1F("Histo_Theta", "", 40, 0, 80);
	  hEnRec_750    = new TH1F("Histo_EnergyRec750", "", 30, 16.0, 18.0);
	  pTriggerET4   = TProfile("pTET4" ,"",15,16,18.5);
	  pTriggerETot  = TProfile("pTETot" ,"",15,16,18.5);
	  hprofEn0      = TProfile("hProfEn0","",50,0,1000);
	  hprofEn1      = TProfile("hProfEn1","",50,0,1000);
	  hprofEn2      = TProfile("hProfEn2","",50,0,1000);
	  hprofEn3      = TProfile("hProfEn3","",50,0,1000);
	  hprof         = TProfile("hprof" ,"",50,0,1000);
	  hprof0        = TProfile("hprof0","",50,0,1000);
	  hprof1        = TProfile("hprof1","",50,0,1000);
	  hprof2d       = TProfile2D("hprof2d","Profile S250vssecTheta",10,1,1.4,8,0.6,2.2,-3.8,-1.2);//uso per minuit e Root per parametri
  	  hprof2        = TProfile("hprof2","",50,0,1000);
	  hprof3        = TProfile("hprof3","",50,0,1000);
	  ldf3D	        = new TH3D("ldfpar3D", "",10, 1, 1.4, 8, 1, 2.6, 10, -4, -1.5);
	  ldf2D	        = new TH2D("ldfpar2D", "",10, 1, 1.4, 8, 1, 2.6);
//	  ldf2D_S250    = new TH2D("ldfpar2D250", "",8, 1, 2.6, 10, -3, -1.5);
	  ldf2D_S250    = new TH2D("ldfpar2D250", "",8, 0.6, 2.2, 50 , -2.7, -1.9);
          profBetaError = TProfile("betaErr","",12,0.9,2.1);
	  ldf2DRopt	= new TH2D("ThRopt2D", "",150, 0, 300,50 , 0, 500);
	  evtDistrib	= new TH2D("evtDistrib", "",300, 15,18,100,1,2);
	  ldfGraph2D    = new TGraph2DErrors(); 
	  ldfPura       = new TGraphErrors(); 
	  ldfModel      = new TGraphErrors(); 
	  hprofRopt     = TProfile("hprofRopt" ,"",60,0,300, "s");
}

double zenithBin(vector <double> zenithVector){
int N = zenithVector.size();
double maxTheta = 0;
for(int i=0; i<N; i++){
   if(zenithVector[i]>maxTheta) maxTheta=zenithVector[i];
   }
return maxTheta;
}


const double kPi = 3.1415926535897932384626;
bool FdCuts(const RecEvent *recEvent, const unsigned int eyeId) {

  /// Has eye ???
  if(!(*recEvent).HasEye(eyeId))
   return false;

  /// Minimum level of reconstruction.
  if((*recEvent).GetEye(eyeId).GetRecLevel() < eHasAxis)
   return false;
  /// Check for hybrid reconstruction.
  if(!(*recEvent).GetEye(eyeId).IsHybridEvent())
   return false;

  /// Geometry sanity cuts.
  const FdRecGeometry &recGeometry =
   (*recEvent).GetEye(eyeId).GetFdRecGeometry();

  /// No horizontal events.
  const double chi0 = recGeometry.GetChi0();
  if((chi0 <= 0.) || (chi0 >= kPi))
   return false;

  /// No laser or back-eye events.
  const double Rp = recGeometry.GetRp();
  if(Rp <= 0.)
   return false;

  /// Olny vertical events.
  const FdRecShower &recShower = (*recEvent).GetEye(eyeId).GetFdRecShower();

  const double zenith = recShower.GetZenith();
  if(zenith*TMath::RadToDeg() > 70.)
   return false;

  /// Minimum level of reconstruction.
  if((*recEvent).GetEye(eyeId).GetRecLevel() < eHasEnergy)
   return false;

const GenShower &genShower = (*recEvent).GetGenShower();


const std::vector<double> &xDepth = recShower.GetDepth();
double xDepth_Max = 0;
double xDepth_min = 10000;

           for(int j=0; j<xDepth.size(); j++) {
             if(xDepth[j] > xDepth_Max) {
               xDepth_Max = xDepth[j];
              }
             if(xDepth[j] < xDepth_min) {
               xDepth_min = xDepth[j];
              }

        }

        //cout << "Max= " << xDepth_Max << " " << "min= " << xDepth_min << endl;



  /// Energy treshold.
  /// const double eTotal = recShower.GetEnergy();
  /// if(eTotal < kEnergyMin)
  ///  return false;
  /// Xmax observed.
  double XmaxMC          = genShower.GetXmax();
  //double XmaxMCInterpolated = genShower.GetXmaxInterpolated();
 // double XmaxMC2 = genShower.GetXmaxGaisserHillas();
  const double xMax      = recShower.GetXmax();
  const double xTrackMin = recShower.GetXTrackMin();
  const double xTrackMax = recShower.GetXTrackMax();
 //cout << "xmaxMC2= " << XmaxMC2 << " xMaxMC=" << XmaxMC <<  endl;

  if((xMax < xDepth_min) || (xMax > xDepth_Max))
   return false;

  return true;

} /// End FdCuts()

//----------

void FdInfoEye(const RecEvent *recEvent, const unsigned int eyeId, std::fstream &outFile, struct var *Eye_var) {

  const FdRecShower  &recShower = (*recEvent).GetEye(eyeId).GetFdRecShower();

  unsigned int evt   = (*recEvent).GetSDEvent().GetEventId();
  unsigned int gps   = (*recEvent).GetSDEvent().GetGPSSecond();



  double zenith      = recShower.GetZenith();
  double zenith_err  = recShower.GetZenithError();

  double azimuth     = recShower.GetAzimuth();
  double azimuth_err = recShower.GetAzimuthError();

  double xCore       = recShower.GetCoreSiteCS().X();
  // double xCore       = recShower.GetCoreUTMCS().X();
  double xCore_err   = recShower.GetCoreEastingError();

  double yCore       = recShower.GetCoreSiteCS().Y();
  // double yCore       = recShower.GetCoreUTMCS().Y();
  double yCore_err   = recShower.GetCoreNorthingError();

  double eEM         = recShower.GetEcal();
  double eTot        = recShower.GetEnergy();
  double eTotErr     = recShower.GetEnergyError();

  double xMax        = recShower.GetXmax();
  double xMaxErr     = recShower.GetXmaxError();


  outFile << "  " << evt
          << "  " << gps
          << "  " << TMath::RadToDeg() * zenith
          << "  " << TMath::RadToDeg() * zenith_err
          << "  " << TMath::RadToDeg() * azimuth
          << "  " << TMath::RadToDeg() * azimuth_err
          << "  " << (unsigned int) xCore
          << "  " << xCore_err
          << "  " << (unsigned int) yCore
          << "  " << yCore_err
          << "  " << xMax
          << "  " << xMaxErr
          << "  " << eEM       / 1.e+18
          << "  " << eTot      / 1.e+18
          << "  " << eTotErr   / 1.e+18
          << endl;


  Eye_var->FdInfo_xMax = xMax;
  Eye_var->FdInfo_eTot = eTot;
  Eye_var->FdInfo_zenith = zenith;
  Eye_var->FdInfo_azimuth = azimuth;
  Eye_var->FdInfo_xCore = xCore;
  Eye_var->FdInfo_yCore = yCore;



} /// End FdInfoEye()







double energyCalib(std::string showerSizeLabel, double showerSize, double theta) { 

/*  <!-- Energy calibration:
     attenuation:
       x(theta) = cos^2(theta) - cos^2(38*degree)
       CIC(x) = 1 + attenuationPar1 * x + attenuationPar2 * x^2 +
                    attenuationPar3 * x^3
       S38 = S(rOpt) / CIC(x(theta))

     energy = energyS38Const * pow(S38, energyS38Slope)

  -->
*/
std::vector<double> cicParameters(3);
double cosRef = 0;
double refAngle_433, refAngle_750;


if(showerSizeLabel == "S250"){
/*
refAngle_433 = 33;     
cicParameters[0] = 2.244;
cicParameters[1] = -0.608;
cicParameters[2] = -1.921;
cosRef = cos(TMath::DegToRad()*(refAngle_433));
*/
//cic ropt 250 e old LDF
refAngle_433 = 30;     
cicParameters[0] =  1.92174;//cambia segno
cicParameters[1] =  -1.03300 ;
cicParameters[2] =  -3.25485;//cambia segno
cosRef = cos(TMath::DegToRad()*(refAngle_433));
//cic ropt 300 e new LDF
/*refAngle_433 = 30;     
cicParameters[0] =  1.76083;
cicParameters[1] =  -1.14986 ;
cicParameters[2] =  - 1.55693;
cosRef = cos(TMath::DegToRad()*(refAngle_433));
*/


}


if(showerSizeLabel == "S450"){

refAngle_750 = 35;
cicParameters[0] = 1.639737;
cicParameters[1] = -1.400067;
cicParameters[2] = -2.190876;
cosRef = cos(TMath::DegToRad()*(refAngle_750));
}


const double x = std::sqrt(cos(theta)) - std::sqrt(cosRef);
double cic = 0;
double s1000 = showerSize;
double  sRefAngle = 0;

 for (vector<double>::const_reverse_iterator pIt = cicParameters.rbegin(), end = cicParameters.rend(); pIt != end; ++pIt)
     cic = *pIt + x * cic;
     cic = 1 + x * cic;

    if(cic > 0) {
      sRefAngle = s1000 / cic;
      }

return sRefAngle;
 }


// funzioni per core in hexagon aeralet
  bool pointInTriangle(double x1, double y1, double x2, double y2, double x3, double y3, double px, double py)
  {
  double denominator = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
  double alpha = ((y2 - y3)*(px - x3) + (x3 - x2)*(py - y3)) / denominator;
  double beta  = ((y3 - y1)*(px - x3) + (x1 - x3)*(py - y3)) / denominator;
  double gamma = 1 - alpha - beta;

  if(0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1)
    return true;

  return false;
  }






void plotResidual(std::map<int,var_SD> mymap, char *name){

int n433_750New  = -1; 
TGraph *f    = new TGraph();
std::fstream file2;
file2.open("indiciEventiBad.txt", file2.out); //e´ per leggere e scrivere
// o usare std::ofstream file2; // ofstream e´ solo pre scrivere
//file2.open("indiciEventiBad.txt");

for(std::map<int,var_SD>::iterator it=mymap.begin(); it!=mymap.end(); ++it) {
 
// it->first tiene el key del map
 // it->second tiene el valor del mapi, en mi caso la struct
 // Hago algo si tengo los dos valores
//cout << "idEventoMap "<< it->first<< endl;

if(it->second.S_opt433 > 0 && it->second.S_opt750 > 0) {
     n433_750New++;
//cout << n433_750New << " " << it->second.S_refAngle433 <<" " <<  it->second.S_refAngle750 << endl; 

f->SetPoint(n433_750New, it->second.S_refAngle433, it->second.S_refAngle750);
}
}

//f->Write("fGraphProva");


//Fit
double var_yFit;
double res;
int N = 0;
TGraph *gRes = new TGraph();


 //TF1 *f2 = new TF1("f2",  "[0]+[1]*x", 0, 200);
   TF1 *f2 = new TF1("f2",  "pol1", 0, 200);
   f->Fit(f2);
   double p0 = f2->GetParameter(0);
   double p1 = f2->GetParameter(1);


for(std::map<int,var_SD>::iterator it=mymap.begin(); it!=mymap.end(); ++it) {
   if(it->second.S_opt433 > 0 && it->second.S_opt750 > 0) {
     N++;     

//var_yFit e´ l´S35 fit
var_yFit = p0 + p1*(it->second.S_refAngle433);
res = ((it->second.S_refAngle750) - var_yFit)/var_yFit;

if(res>1) {
cout << "res value "  << res << "Eventid=" << it->first << endl;
file2 << it->first << endl; 
}
gRes->SetPoint(N, it->second.S_refAngle433, res);    

     }
   }
//file2.close();

gRes->Write(name);



}








using namespace std;

std::vector<limits> AngularBinCalculator()
{
    /**********/
    /* Inputs */
    /******* **/

    //Angular limits in deg
    double theta_1 = 0;
    double theta_2 = 45;
  
    std::vector<limits> pippo;
    //Number of bins
    double nbins = 4;

    vector<double> vLowBounds(nbins,0);
    vector<double> vUpBounds(nbins,0);

    theta_1 *= TMath::Pi()/180;//per passare a rad
    theta_2 *= TMath::Pi()/180;

    double totalIntensity = pow(cos(theta_1),2) - pow(cos(theta_2),2);
    double binIntensity = totalIntensity/nbins;

    cout << "totalIntensity: " << totalIntensity << endl;
    cout << "binIntensity: " << binIntensity << endl;

    int i=0;
    limits intervals;

    while (i<nbins) {

        if (i==0)
            vLowBounds.at(0) = theta_1;
        else
            vLowBounds.at(i) = vUpBounds.at(i-1);

        if (i == nbins-1)
            vUpBounds.at(nbins-1) = theta_2;
        else
            vUpBounds.at(i) = TMath::ACos(sqrt(pow(cos(vLowBounds.at(i)),2) - binIntensity));

        //Printing out in my code format ..
        cout << "Bin " << i << ": (" << vLowBounds.at(i) * 180/TMath::Pi() << " ; " << vUpBounds.at(i) * 180/TMath::Pi() << ") deg" << endl;
        intervals.max =vUpBounds.at(i);
        intervals.min =vLowBounds.at(i);
        pippo.push_back(intervals);
        ++i;
    }

    for (int j=0; j<nbins; ++j)
        cout << "vThetaMin.push_back(" << vLowBounds.at(j) * 180/TMath::Pi() << ");" << endl;

    cout << endl;

    for (int j=0; j<nbins; ++j)
        cout << "vThetaMax.push_back(" << vUpBounds.at(j) * 180/TMath::Pi() << ");" << endl;

return pippo;

}


string timeStamp(SDEvent& theSdEvent)
{
        const unsigned int yymmdd = theSdEvent.GetYYMMDD();
        const unsigned int hhmmss = theSdEvent.GetHHMMSS();
        const unsigned int year   = yymmdd/10000+2000;
        const unsigned int month  = (yymmdd%10000)/100;
        const unsigned int day    =(yymmdd%100);
        const unsigned int hour   = (int)(hhmmss/10000.);
        const unsigned int minute = (hhmmss%10000)/100;
        const unsigned int second = (hhmmss%100);
        unsigned int       gps    = theSdEvent.GetGPSSecond();

        stringstream ss;
        ss << year << '/'<< month <<'/'<<  day << ' '<< hour << ':'<< minute << ':'<< second << ' '<< "gps=" <<gps;

        string x = ss.str();

        return x;
}


double computeNKG(double dist, double Sref, double beta){
//double computeNKG(double dist, double Sref, double beta, double gamma){

//  double sLDF = Sref*pow(dist/1000., beta)*pow((dist+700.)/1700., gamma);
  double sLDF = Sref*pow(dist/250., beta)*pow((dist+700.)/950., beta);
  return sLDF;
}

//the LDFFinder_infill uses NKGAS
//  <!-- type of LDF model -->
//  <ldfType> NKGAS </ldfType>
//shapeMOdel function from NKGASLDF.h  
double ShapeModel(const double cosTheta, const double showerSize)
    {
      const double lgSRef = std::log10(showerSize);
      const double secTheta = 1/cosTheta;
  //<!-- updated values from Multi event fit (dipl. thesis Alex S.) -->    
     /* const double a0 = -2.378716;//pars from ldf par for the infill (Alex Shultz.)
      const double a1 = -0.190785;
      const double b0 = 0.693927;
      const double b1 = 0.096167;
      const double c0 = 0.101887;
      const double c1 = -0.182174;
*/
      const double a0 = -4.550;//pars from 433 parametrization S.Messina from GAP-2014-094
      const double a1 = 1.160;
      const double b0 = 3.414;
      const double b1 = -1.946;
      const double c0 = -0.994;
      const double c1 = 0.736;
      double beta;

      beta = a0 + a1*lgSRef + secTheta*(b0 + b1*lgSRef + secTheta*(c0 + c1*lgSRef));

      return beta;
    }




void stationFunc(const std::vector<SdRecStation>& fStation, TProfile *hprof, double S250, double beta433) { 
	    for (std::vector<SdRecStation>::const_iterator sIt=fStation.begin(); sIt!=fStation.end(); ++sIt) {
                if (sIt->IsCandidate()) {

                   double signal = sIt->GetTotalSignal();
//		   if (signal<5) continue;
//	if (signal >200) continue;
                   double rSP = sIt->GetSPDistance();
		   double expectedSignal = computeNKG(rSP, S250, beta433);
	           if (expectedSignal<5) continue;
	 	   if (expectedSignal >200) continue;
//cout << " expectedsignal after cut "<< expectedSignal << "--- signal " << signal <<  endl;
                   double ratio = signal/S250;
                   signal_v.push_back(signal/S250);
                   rSP_v.push_back(rSP);
                   hprof->Fill(rSP,ratio);           
                   }
                 }

        }


void binomialErrors(TProfile prof,TGraphAsymmErrors*& g)
{
     for (int i=1; i<=prof.GetNbinsX(); ++i) {
         double eff = prof.GetBinCenter(i);
         double p = prof.GetBinContent(i);
         double n = prof.GetBinEntries(i);
         double x = prof.GetBinContent(i) * prof.GetBinEntries(i);
         double CL = 0.90;
         double alpha = 1-CL;

         double y1 = 0;
         double y2 = 0;

/*         if (p>0 && p<1) {
             y1 = ROOT::Math::beta_quantile(alpha/2,x,n-x+1);
             y2 = ROOT::Math::beta_quantile(1-alpha/2,x+1,n-x);
         }

  */       g->SetPoint(i-1,eff,p);

         if (y1 != 0 && y2 !=0)
             g->SetPointError(i-1,0,0,p-y1,y2-p);
         else
             g->SetPointError(i-1,0,0,0,0);

     }

}

double calculateMax(vector<double> en) {
double en_max = -4;
 for (int i=0; i<en.size(); i++) {
    if (en[i]>en_max) en_max = en[i];
    }
return en_max;
}

double calculateMin(vector<double> en) {
double en_min = 1000;
 for (int i=0; i<en.size(); i++) {
    if (en[i]<en_min) en_min = en[i];
    }
return en_min;
}


