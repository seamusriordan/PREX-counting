#include <iostream>
#include <string>

#include "MPDGEMTracker.h"
#include "MPDGEMPlane.h"

#include "Plane.h"

MPDGEMTracker::MPDGEMTracker( const char* name, const char* desc, THaApparatus* app ):
    TreeSearch::GEMTracker(name,desc,app) {
}

MPDGEMTracker::~MPDGEMTracker(){
    return;
}


TreeSearch::Plane* MPDGEMTracker::MakePlane( const char* name, const char* description,
        THaDetectorBase* parent  ) const {
    return new MPDGEMPlane( name, description, parent );
}

