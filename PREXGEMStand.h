#ifndef PREXGEMSTAND_H
#define PREXGEMSTAND_H
#include <vector>
#include <THaTrackingDetector.h>

class THaRunBase;
class THaApparatus;
class THaEvData;
class PREXGEMPlane;
class THaCrateMap;

class PREXGEMStand : public THaTrackingDetector {
    public:
        PREXGEMStand( const char *name, const char *description = "",
                THaApparatus *app = 0 );

        virtual ~PREXGEMStand();

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual EStatus Init( const TDatime& date );

        virtual Int_t   ReadDatabase( const TDatime& date );

        virtual Int_t   CoarseTrack( TClonesArray& tracks );
        virtual Int_t   FineTrack( TClonesArray& tracks );
        virtual Int_t   DefineVariables( EMode mode = kDefine );
        virtual void    Print(const Option_t* opt) const;
        virtual void    SetDebug( Int_t level );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:
        std::vector <PREXGEMPlane *> fPlanes;

        THaCrateMap *fCrateMap;
        ClassDef(PREXGEMStand ,0)
};

#endif
