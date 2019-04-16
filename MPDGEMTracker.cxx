#include <iostream>
#include <string>

#include "GEMTracker.h"
#include "GEMPlane.h"

GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app ):
    TreeSearch::GEMTracker(name,desc,app) {
}

GEMTracker::~GEMTracker(){
    return;
}


virtual TreeSearch::Plane* MakePlane( const char* name, const char* description = "",
        THaDetectorBase* parent = 0 ) const {
    return new GEMPlane( name, description, parent );
}

