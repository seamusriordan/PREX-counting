#ifndef MPDModule_
#define MPDModule_

/////////////////////////////////////////////////////////////////////
//
//   MPDModule
//   This is the MPD INFN (GEM readout for SBS) module decoder
//   Note, it must inherit from VmeModule.
//
//   other steps
//   1. Register (add "DoRegister" call and include the header)
//   note: the number (4444) that is registered must appear in db_cratemap.dat
//   2. Add to namespace Decoder.h
//   3. Add to Makefile
//   4. Add to haDecode_LinkDef.h
//   5. Add line(s) for [crate,slot] in db_cratemap.dat
//
//   if the pre-compiler flag "LIKEV792" is defined, the decoding is
//   sort of like a V792 ... as an example.
//
/////////////////////////////////////////////////////////////////////


#define MAXHIT    2048

#include "VmeModule.h"

namespace Decoder {

  class MPDModule : public VmeModule {

  public:

    MPDModule() {};
    MPDModule(Int_t crate, Int_t slot);
    virtual ~MPDModule();


    virtual Int_t GetData(Int_t adc, Int_t sample, Int_t chan) const;
    // Stubs for virtual overloading
    virtual UInt_t GetData(Int_t)        const { return 0;}
    virtual Int_t  GetData(Int_t, Int_t) const { return 0;}
    virtual Int_t  GetData(Decoder::EModuleType, Int_t, Int_t)        const { return 0;}
    virtual Int_t  GetData(Decoder::EModuleType, Int_t, Int_t, Int_t) const { return 0;}

    virtual void Init();
    virtual void Clear(const Option_t *opt);
    virtual Int_t Decode(const UInt_t *p); // { return 0; };
    
    virtual Int_t LoadSlot(THaSlotData *sldat,  const UInt_t *evbuffer, Int_t pos, Int_t len);
    virtual Int_t LoadSlot(THaSlotData *sldat,  const UInt_t *evbuffer, const UInt_t *pstop);
//#endif

  private:

    // configuration parameters
    Int_t fAcqMode; // normal, zero suppression, histogram, synch ...
    Int_t fSamplePeriod; // 25 ns, 75 ns ...
    Int_t fNumSample; // number of sample / event
    
    Int_t fNumADC; // number of ADC fifos (number of front end cards served by the MPD)
    
    // current indices
    Int_t fIdxA; // ADC
    Int_t fIdxS; // Sample
    Int_t fIdxC; // Channel

    Int_t fIdxMPD; // MPD ID
  
    Int_t fCountS; // Sample Counter from electronics
    Int_t fCountW; // Word 

    Int_t fNumHits;

    std::vector<Int_t> fFrameHeader;  // Frame Header
    std::vector<Int_t> fFrameTrailer;  // Frame Trailer

    static TypeIter_t fgThisType;

    // linearization of the indeces 
    inline Int_t as2i(Int_t adc, Int_t sample) const {
      return adc*fNumSample + sample;
    };

    inline Int_t asc2i(Int_t adc, Int_t sample, Int_t chan) const {
      return adc*fNumSample*fNumChan + sample*fNumChan + chan;
    };
    
    ClassDef(MPDModule,0)  //  INFN MPD Module 

  };

}

#endif
