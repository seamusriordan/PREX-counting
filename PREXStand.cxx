//////////////////////////////////////////////////////////////////////////
//
// Bare-bones SBB BigBite spectrometer class 
// Used for testing.
//
//////////////////////////////////////////////////////////////////////////

#include "PREXStand.h"

using namespace std;

ClassImp(PREXStand)


//_____________________________________________________________________________
PREXStand::PREXStand( const char* name, const char* description ) :
  THaSpectrometer( name, description )
{
  // Constructor. Defines standard detectors

}

//_____________________________________________________________________________
PREXStand::~PREXStand()
{
  // Destructor
}

//_____________________________________________________________________________
  Int_t PREXStand::FindVertices( TClonesArray& /* tracks */ )
{
  // Reconstruct target coordinates for all tracks found.

  // TODO

  return 0;
}

//_____________________________________________________________________________
Int_t PREXStand::TrackCalc()
{
  // Additioal track calculations

  // TODO

  return 0;
}

