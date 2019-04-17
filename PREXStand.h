#ifndef ROOT_TreeSearch_PREXStand
#define ROOT_TreeSearch_PREXStand

#include "THaSpectrometer.h"

class PREXStand : public THaSpectrometer {

    public:
    PREXStand( const char *name, const char *description );
    virtual ~PREXStand();

    virtual Int_t FindVertices( TClonesArray& tracks );
    virtual Int_t TrackCalc();

    protected:
    ClassDef(PREXStand,0) // BigBite spectrometer
};

#endif

