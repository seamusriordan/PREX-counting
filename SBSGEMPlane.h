#include "THaSubDetector.h"

class THaDetectorBase;
class THaEvData;
class THaRunBase;

struct mpdmap_t {
    UInt_t crate;
    UInt_t slot;
    UInt_t mpd_id;
    UInt_t gem_id;
    UInt_t adc_id;
    UInt_t i2c;
    UInt_t pos;
    UInt_t invert;
};

#define N_APV25_CHAN    128
#define N_MPD_TIME_SAMP 6
#define MPDMAP_ROW_SIZE 8
#define N_CMN_CHAN 100

const int APVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113, 25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59, 91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93, 125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127, 0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34, 66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100, 12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46, 78, 110, 22, 54, 86, 118, 30, 62, 94, 126};

class SBSGEMPlane : public THaSubDetector {
    public:

        SBSGEMPlane( const char *name, const char *description = "",
                THaDetectorBase* parent = 0 );

        virtual ~SBSGEMPlane();

        virtual void    Clear( Option_t* opt="" );
        virtual Int_t   Decode( const THaEvData& );
        virtual void    Print( Option_t* opt="" ) const;

        virtual Int_t   ReadDatabase(const TDatime& );
        virtual Int_t   DefineVariables( EMode mode );

        virtual Int_t   Begin( THaRunBase* r=0 );
        virtual Int_t   End( THaRunBase* r=0 );

    private:
        std::vector<mpdmap_t>    fMPDmap;
        std::vector<Int_t>       fChanMapData;

        Double_t fZeroSuppressRMS;
        Bool_t fZeroSuppress;

        // Output variables
        Int_t  fNch;
        Int_t *fStrip; // [fNch]
        Int_t *fadc[N_MPD_TIME_SAMP]; 
        // Being obnoxious so we match the stand alone more closely
        Int_t *fadc0; // [fNch]
        Int_t *fadc1; // [fNch]
        Int_t *fadc2; // [fNch]
        Int_t *fadc3; // [fNch]
        Int_t *fadc4; // [fNch]
        Int_t *fadc5; // [fNch]
	Int_t *fadc_sum; //[fNch]
        Double_t *fPedestal;
  	Int_t *fcommon_mode; //[fNch]
        Double_t *fRMS;
        Double_t trigger_time;
        Int_t ev_num;
  
        ClassDef(SBSGEMPlane,0)

};




