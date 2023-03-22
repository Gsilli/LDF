#include <TTimeStamp.h>
#include <iterator>
#include "TProfile.h"
#include <TProfile2D.h>
#include "TFile.h"
#include <utl/ErrorLogger.h>
#include <utl/AugerException.h>
#include <utl/Point.h>
#include <utl/Vector.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AugerUnits.h>
#include <utl/Branch.h>
#include <fwk/CentralConfig.h>
#include "TGraph.h"
#include "TEfficiency.h"
#include "Rtypes.h"
#include <ctime>
#include <evt/Event.h>
#include <evt/ShowerRecData.h>
#include <evt/ShowerSRecData.h>
#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <det/Detector.h>
#include <sdet/SDetector.h>
#include <sdet/Station.h>
#include <map>


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
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <RecEventFile.h>
#include <RecEvent.h>
#include <FDEvent.h>
#include <SDEvent.h>
#include <DetectorGeometry.h>
#include <io/RootFile.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <TSystem.h>

using namespace std;
namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

double RPerp(const TVector3& axis, const TVector3& station){

  const double scal = axis*station;
  return sqrt(station*station - scal*scal);
}


DetectorGeometry* theGeometry;
using namespace sevt;
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

int main(int argc, char **argv)
{

std::string path(argv[1]);
std::size_t pos = path.find("_adst.root");
std::string resultado3 = path.substr(pos-16, 4);
std::fstream outFile;
outFile.open(("output" + resultado3 +".txt").c_str(), outFile.out);

std::vector<std::string> dataFileNames;
for (int i=1; i<argc; ++i) {

     std::string nomeADST(argv[i]);
     dataFileNames.push_back(argv[i]);
}

 RecEvent* theRecEvent = new RecEvent();
 theGeometry = new DetectorGeometry();
 FileInfo fileInfoSi;
 RecEventFile dataFile(dataFileNames);
 dataFile.SetBuffers(&theRecEvent);
 dataFile.ReadDetectorGeometry(*theGeometry);
 dataFile.ReadFileInfo(fileInfoSi);

 SDEvent& theSdEvent = theRecEvent->GetSDEvent();

  while (dataFile.ReadNextEvent() == RecEventFile::eSuccess)
        {
        int eventId2 = theSdEvent.GetEventId();
        std::string eventId = theRecEvent->GetEventId();
        std::string sizeLabel = theSdEvent.GetSdRecShower().GetShowerSizeLabel();
        cout << "Processing event Si" << eventId <<" " << eventId2<< endl;
        ++nEvents;

      const auto& sRecShower = theSdEvent.GetSdRecShower();
      const auto& stationVector = theSdEvent.GetStationVector();
      const auto& badStationVector = theSdEvent.GetBadStationVector();
//    const auto& detector = dataFile.Get<DetectorGeometry>("detectorGeometry");
                                                        
      TBits zero = theRecEvent->GetDetector().GetActiveStations();
      if (!theSdEvent.HasLDF()){
        continue;
      }
      const TVector3& showerCore = sRecShower.GetCoreSiteCS();
      const TVector3& showerAxis = sRecShower.GetAxisSiteCS();

      nEvents++;
      string augerId = theRecEvent->GetAugerId();
      if(augerId[0] != '1'){
        augerId.erase(augerId.size()-1);
      }

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

      unsigned int gps = theSdEvent.GetGPSSecond();
      int sdId = theSdEvent.GetEventId();
      double theta = sRecShower.GetZenith(), thetaError = sRecShower.GetZenithError();
      double phi = sRecShower.GetAzimuth(), phiError = sRecShower.GetAzimuthError();
      double xCore = sRecShower.GetCoreSiteCS()[0], yCore = sRecShower.GetCoreSiteCS()[1];
      double dxCore = sRecShower.GetCoreEastingError(), dyCore = sRecShower.GetCoreNorthingError();
      //int neighbors = theSdEvent.GetT5PostActiveNeighbors().size();
      //int position = theSdEvent.GetT5PostCoreTriangle();
      double energy = sRecShower.GetEnergy();

      double energyError = sRecShower.GetEnergyError();
      double showerSize = sRecShower.GetShowerSize(); //GetLDF().GetS1000()
      double showerSizeErrorSyst = sRecShower.GetShowerSizeSys();
      double showerSizeErrorStat = sRecShower.GetShowerSizeError();
      double beta = sRecShower.GetBeta(), betaError = sRecShower.GetBetaError();
      double gamma = sRecShower.GetGamma(), gammaError = sRecShower.GetGammaError();
      double ropt = sRecShower.GetROpt();
      int fitBeta = 0, fitGamma = 0;
      if(betaError){
        fitBeta = 1;
      } 
      if(gammaError){
        fitGamma = 1;
      }
      double angularReducedChi2 = sRecShower.GetAngleChi2();
      if(sRecShower.GetAngleNDoF()){
        angularReducedChi2 = sRecShower.GetAngleChi2() / sRecShower.GetAngleNDoF();
      }
      double ldfReducedChi2 = sRecShower.GetLDFChi2();
      if(sRecShower.GetLDFNdof()){
        ldfReducedChi2 = sRecShower.GetLDFChi2() / sRecShower.GetLDFNdof();
      }
      double globalReducedChi2 = 0.;
      int multiplicity = theSdEvent.GetNumberOfCandidates() + nSilent;
/*      if(showerSize < 100){
        continue;
      }*/


     outFile
      <<theSdEvent.Is6T5()<<" "<<theSdEvent.Is5T5()<<" "
      <<augerId<<" "<<sdId<<" "<<gps<<" "
      <<energy<<" "<<energyError<<" "
      <<showerSize<<" "<<showerSizeErrorStat<<" "<<showerSizeErrorSyst<<" "<<ropt<<" "
      //<<neighbors<<" "<<position<<" "
      <<fitBeta<<" "<<beta<<" "<<betaError<<" "<<fitGamma<<" "<<gamma<<" "<<gammaError<<" "
      <<theta<<" "<<thetaError<<" "<<phi<<" "<<phiError<<" "
      <<xCore<<" "<<dxCore<<" "<<yCore<<" "<<dyCore<<" "
      <<angularReducedChi2<<" "<<ldfReducedChi2<<" "<<globalReducedChi2<<" "
      <<multiplicity<<" ";

      double maxDistance = 0.;
      for (const auto& station : stationVector) {
        int stId = station.GetId();
        zero[stId] = false;
        if (!station.IsCandidate())
          continue;
        double stSignal = station.GetTotalSignal(), stSignalError = station.GetTotalSignalError();
        double stSignalRec = 0.;
        double stDist = station.GetSPDistance();
        if(stDist > maxDistance){
          maxDistance = stDist;
        }
        double stDistError = station.GetSPDistanceError();
        if (station.IsLowGainSaturated() && station.GetRecoveredSignal() > 0) {
          stSignalRec = -station.GetRecoveredSignal(); // minus for saturated
          stSignalError = pow(station.GetTotalSignalError(),2.) + pow(station.GetRecoveredSignalError(),2.);
        }
        outFile<<stId<<" "<<station.IsLowGainSaturated()<<" "<<stDist<<" "<<stDistError<<" "<<stSignal<<" "<<stSignalRec<<" "<<stSignalError<<" ";
        outFile<<station.IsT1Threshold()<<" "<<station.IsT2Threshold()<<" "<<station.IsToT()<<" "<<station.IsToTd()<<" "<<station.IsMoPS()<<" ";
      }

      // extraction of silent (zero-signal) stations information
      const double nCrowns = 2;
      const double gridSpacing = 1500;
      const double tolerance = 150;
      const double distanceCut = maxDistance + nCrowns*gridSpacing;
      for (unsigned int id = 0, n = zero.GetNbits(); id < n; ++id) {
        if(zero[id]) {
          double distance = RPerp(showerAxis, theGeometry->GetStationPosition(id) - showerCore);
          if((distance - distanceCut) < tolerance) {
            outFile<<id<<" "<<0<<" "<<distance<<" "<<100.<<" "<<0<<" "<<0<<" "<<0<<" ";
            outFile<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" ";
          }
        }
      }
      //
      outFile<<std::endl;

  cout << "nEvents = " << nEvents << endl;

  //return EXIT_SUCCESS;


















        }//fine while


outFile.close();

return 0;



} //fine main



