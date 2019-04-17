#ifndef MPDGEMTRACKER_H
#define MPDGEMTRACKER_H
#include <GEMTracker.h>

namespace TreeSearch {
    class Plane;
};

class MPDGEMTracker : public TreeSearch::GEMTracker {
    public:
        MPDGEMTracker( const char *name, const char *description = "",
                THaApparatus *app = 0 );

        virtual ~MPDGEMTracker();

    protected:
        virtual TreeSearch::Plane* MakePlane( const char* name, const char* description = "",
                THaDetectorBase* parent = 0 ) const;

    public:
        ClassDef(MPDGEMTracker,0)
};

#endif
