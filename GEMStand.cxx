#include <iostream>
#include <string>

#include "THaTrackingDetector.h"
#include "THaRunBase.h"
#include "THaCrateMap.h"
#include "THaAnalysisObject.h"

#include "GEMStand.h"
#include "GEMPlane.h"

GEMStand::GEMStand( const char* name, const char* desc, THaApparatus* app ):
    THaTrackingDetector(name,desc,app) {

        fPlanes.clear();

        fCrateMap = 0;
}

GEMStand::~GEMStand(){
    return;
}


THaAnalysisObject::EStatus GEMStand::Init( const TDatime& date ){
    assert( fCrateMap == 0 );

    THaAnalysisObject::EStatus status = THaTrackingDetector::Init(date);

    if( status == kOK ){
        for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
            status = (*it)->Init(date);
            if( status != kOK ){
                return status;
            }
        }
    } else {
        return kInitError;
    }

    return kOK;
}


Int_t GEMStand::ReadDatabase( const TDatime& date ){
    std::cout << "[Reading GEMStand database]" << std::endl;

    fIsInit = kFALSE;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    Int_t err = ReadGeometry( file, date );
    if( err ) {
        fclose(file);
        return err;
    }

    std::string planeconfig;
    std::vector<Int_t> *cmap = new std::vector<Int_t>;

    DBRequest request[] = {
        { "planeconfig",       &planeconfig,       kString   },
        { "cratemap",          cmap,               kIntV     },
        {0}
    };

    Int_t status = kInitError;
    err = LoadDB( file, date, request, fPrefix );
    fclose(file);

    std::vector<std::string> planes = vsplit(planeconfig);
    if( planes.empty()) {
            Error("", "[GEMStand::ReadDatabase] No planes defined");
    }

    for (std::vector<std::string>::iterator it = planes.begin() ; it != planes.end(); ++it){
        fPlanes.push_back(new GEMPlane( (*it).c_str(), (*it).c_str(), this));
    }

    status = kOK;

    if( status != kOK )
        return status;

    fIsInit = kTRUE;
    
    return kOK;
}


Int_t GEMStand::Begin( THaRunBase* run ){
    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->Begin(run);
    }

    return 0;
}

void GEMStand::Clear( Option_t *opt ){
    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->Clear(opt);
    }

    return;
}

Int_t GEMStand::Decode(const THaEvData& evdata ){
//    std::cout << "[GEMStand::Decode]" << std::endl;

    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->Decode(evdata);
    }

    return 0;
}


Int_t GEMStand::End( THaRunBase* run ){
    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->End(run);
    }


    return 0;
}

void GEMStand::Print(const Option_t* opt) const {
    std::cout << "GEM Stand " << fName << " with " << fPlanes.size() << " planes defined:" << std::endl;
    /*
    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        std::cout << "\t"
        (*it)->Print(opt);
    }
    */
    for( unsigned int i = 0; i < fPlanes.size(); i++ ){
        fPlanes[i]->Print(opt);
    }

    return;
 }


void GEMStand::SetDebug( Int_t level ){
      THaTrackingDetector::SetDebug( level );
    for (std::vector<GEMPlane *>::iterator it = fPlanes.begin() ; it != fPlanes.end(); ++it){
        (*it)->SetDebug(level);
    }

    return;
}

Int_t GEMStand::DefineVariables( EMode mode ){
    if( mode == kDefine and fIsSetup ) return kOK;
    fIsSetup = ( mode == kDefine );
    RVarDef vars[] = {
//        { "trkstat", "Track reconstruction status",  "fTrkStat" },
        { 0 },
    };
    DefineVarsFromList( vars, mode );

    return 0;
}


Int_t GEMStand::CoarseTrack( TClonesArray& ){
    return 0;
}
Int_t GEMStand::FineTrack( TClonesArray& ){
    return 0;
}

