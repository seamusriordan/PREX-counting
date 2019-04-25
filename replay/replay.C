
#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#pragma cling load("../TreeSearch/libTreeSearch-GEM.so");
#pragma cling load("../PREX-counting/libprexCounting.so");
//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"

void replay(){
    PREXStand *prex= new PREXStand("prex", "Generic apparatus");
    MPDGEMTracker *gems = new MPDGEMTracker("gems", "Collection of GEMs in stand");

    prex->AddDetector(gems);

  //
  //  Steering script for Hall A analyzer demo
  //
  
  // Set up the equipment to be analyzed.
  
  // add the two spectrometers with the "standard" configuration
  // (VDC planes, S1, and S2)
  // Collect information about a easily modified random set of channels
  // (see DB_DIR/db_D.dat)
  /*
  THaApparatus* DECDAT = new THaDecData("D","Misc. Decoder Data");
  gHaApps->Add( DECDAT );
  */
  

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  gHaApps->Add(prex);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;
  
  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
//  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  //THaRun* run = new THaRun( "/home/chad/data/test_589.dat" );
  int argc=0;
  //std::cout<<"enter file name: "<<std::endl;
  //cin>>argc;
  argc = 409;
  THaRun* run = new THaRun(Form("test_%d.dat",argc));
  run->SetLastEvent(-1);

  run->SetDataRequired(0);
  run->SetDate(TDatime());

  analyzer->SetVerbosity(10);
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile(Form("test_%d.root",argc));
  // File to record cuts accounting information
  analyzer->SetSummaryFile("summary_example.log"); // optional
  
  //analyzer->SetCompressionLevel(0); // turn off compression
  analyzer->Process(run);     // start the actual analysis
}
