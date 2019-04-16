/////////////////////////////////////////////////////////////////////
//
//   MPDModule
//   This is the MPD module decoder; based on SkeletonModule
//   (https://github.com/JeffersonLab/analyzer)
//
/////////////////////////////////////////////////////////////////////

#include "MPDModule.h"
#include "THaSlotData.h"

using namespace std;

namespace Decoder {

  Module::TypeIter_t MPDModule::fgThisType =
    DoRegister( ModuleType( "Decoder::MPDModule" , 3561 ));

  MPDModule::MPDModule(Int_t crate, Int_t slot) : VmeModule(crate, slot) {
    fDebugFile=0;
    Init();
  }
  
  MPDModule::~MPDModule() {
    
  }
  
  void MPDModule::Init() { 
    Module::Init();
    fDebugFile=0;
    Clear("");
    fName = "MPD Module";
  }

  Int_t MPDModule::LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, const UInt_t *pstop) {
      return LoadSlot( sldat, evbuffer, 0, pstop-evbuffer );
  }
  
  Int_t MPDModule::LoadSlot(THaSlotData *sldat, const UInt_t* evbuffer, Int_t pos, Int_t len) {
      const UInt_t *p = &evbuffer[pos];
      UInt_t data;
      fWordsSeen = 0;

      Int_t status;

      // From stand alone decoder
      // We declare an effective channel from the MPD ID 
      // and ADC channel
      Int_t mpdID = -1;
      Int_t adcCh = -1;
      Int_t effCh = 0;

      UInt_t data_count = 0;


      // Variables for trigger times
      Int_t trigger_time;

      //  v5 decoder
      int ii,jj,kk,ll;
      int thesewords;

      jj =  0;


      while( jj < len ){
          thesewords = p[jj++] & 0xFFFFFF;
          // printf("===============================================================================\n");
          // printf("=    CRATE   %d   ======    SLOT   %d   =======================================\n", fCrate, fSlot);
          // printf("BLOCK HEADER  %06x\n", thesewords);
          // printf("Good? (0)       %x\n", (thesewords & 0xe00000) >> 21);
          // printf("Module ID       %d\n", (thesewords & 0x1F0000) >> 16 );
          // printf("EVENT_PER_BLOCK %d\n", (thesewords & 0x00FF00) >> 8 );
          // printf("BLOCK COUNT     %d\n", (thesewords & 0x0000FF) >> 0);

          if( (thesewords & 0xe00000) >> 21 != 0 ){
              fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK HEADER NOT FOUND\n", __LINE__);
              return -1;
          }

          mpdID = (thesewords & 0x1F0000) >> 16;

          int nevent = (thesewords & 0x00FF00) >> 8;

          for( ii = 0; ii < nevent; ii++ ){
              thesewords = p[jj++] & 0xFFFFFF;

              // printf("EVENT HEADER  %06x\n", thesewords);
              // printf("Good? (4)       %x\n", (thesewords & 0xF00000) >> 20);
              // printf("EVENT COUNT     %d\n", (thesewords & 0x0FFFFF) >> 0);
              if( (thesewords & 0xF00000) >> 20 != 0x4 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT HEADER NOT FOUND\n", __LINE__);
                  return -1;
              }
	      
	      Int_t event_num = (thesewords & 0x0FFFFF);
	      Int_t effCh_evnum = (mpdID)<<4;
	      status = sldat->loadData("adc",effCh_evnum,event_num,event_num);
	      if( status != SD_OK ) return -1;

              thesewords = p[jj++] & 0xFFFFFF;
              // printf("TRIGGER TIME 1%06x\n", thesewords);
              // printf("Good? (6)       %x\n", (thesewords & 0xF00000) >> 20);
              // printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
	      
	      trigger_time = (thesewords & 0x0FFFFF) >> 0;
	      Int_t effCh_trigger  = (mpdID)<< 5;
	      status = sldat->loadData("adc",effCh_trigger,trigger_time,trigger_time);
	      if( status != SD_OK ) return -1;
	      
              if( (thesewords & 0xF00000) >> 20 != 0x6 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 1 WORD NOT FOUND\n", __LINE__);
                  return -1;
              }

              thesewords = p[jj++] & 0xFFFFFF;
              // printf("TRIGGER TIME 2%06x\n", thesewords);
              // printf("Good? (7)       %x\n", (thesewords & 0xF00000) >> 20);
              // printf("COURSE TIME     %d\n", (thesewords & 0x0FFFFF) >> 0);
	      
	      trigger_time = (thesewords & 0x0FFFFF) >> 0;
	      status = sldat->loadData("adc",effCh_trigger,trigger_time,trigger_time);
	      if( status != SD_OK ) return -1;
	      
              if( (thesewords & 0xF00000) >> 20 != 0x7 ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] TRIGGER TIME 2 WORD NOT FOUND\n", __LINE__);
                  return -1;
              }


              kk = 0;
              while( ((p[jj] & 0xE00000) >> 21 ) == 0x4  ){
                  kk++;
                  //printf("\n[Starting sample %d]\n", kk);

                  thesewords = p[jj++] & 0x1FFFFF;

                  adcCh = thesewords & 0x00000F;

                  // printf("APV HEADER        %06x\n", thesewords);
                  // printf("Headergood? (0) %x\n", (thesewords & 0x1C0000) >> 18);
                  // printf("Baselineval     %x\n", (thesewords & 0x020000) >> 17);
                  // printf("APV HEADER      %x\n", (thesewords & 0x01FFF0) >> 4);
                  // printf("APV ID          %x\n", (thesewords & 0x00000F) >> 0);
                  if( (thesewords & 0x1C0000) >> 18 != 0x0 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA HEADER NOT FOUND\n", __LINE__);
                      return -1;
                  }

                  // Loop while still seeing reduced data
                  while( ((p[jj] & 0x180000) >> 19) == 0x1 ){
                      for( ll = 0; ll < 8; ll++ ){
                          if( ((p[jj] & 0x180000) >> 19) != 0x1 ){
                              break;
                          }
                          //                      printf("%08x  ", p[jj++]);
                          int x_data = p[jj++];

                          data =  x_data& 0x00FFF;
                          //ch   = (x_data& 0x7F000)>>12;
                          //printf("%3d %03x  ", ch, data);

                          // Otherwise we have data
                          effCh = (mpdID) << 8 | adcCh;

                          status = sldat->loadData("adc",effCh, data, data);
                          if( status != SD_OK ) return -1;

                          fWordsSeen++;
                          data_count++;
                      }
                      //printf("\n");
                  }

                  thesewords = p[jj++] & 0x1FFFFF;
                  //printf("APV TRAILER   %06x\n", thesewords);
                  //printf("Good? (8)       %x\n", (thesewords & 0x1E0000) >> 17);
                  //printf("Module ID       %x\n", (thesewords & 0x01F000) >> 12);
                  //printf("Sample Count    %x\n", (thesewords & 0x000F00) >>  8);
                  //printf("Frame Counter   %x\n", (thesewords & 0x0000FF) >>  0);
                  if( (thesewords & 0x1E0000) >> 17 != 0x8 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] APV TRAILER NOT FOUND\n", __LINE__);
                      return -1;
                  }

                  thesewords = p[jj++] & 0x1FFFFF;
                  //printf("TRAILER       %06x\n", thesewords);
                  //printf("Good? (3)       %x\n", (thesewords & 0x180000) >> 19);
                  //printf("Baseline val    %x\n", (thesewords & 0x07FF00) >>  8);
                  //printf("Word count      %x\n", (thesewords & 0x0000FF) >>  0);
                  if( (thesewords & 0x180000) >> 19 != 0x3 ){
                      fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] DATA TRAILER NOT FOUND\n", __LINE__);
                      return -1;
                  }

              }
              //printf("[ %d SAMPLES ]\n", kk);

              thesewords = p[jj++] & 0xFFFFFF;
              // printf("EVENT TRAILER %06x\n", thesewords);
              // printf("Good? (a)       %x\n", (thesewords & 0xF00000) >> 20);
              // printf("N WORDS IN EVT  %d\n", (thesewords & 0x0FFF00) >> 8);
              // printf("FINE TRIGGER T  %d\n", (thesewords & 0x0000FF) >> 0);
	      
	      trigger_time = (thesewords & 0x0000FF) >> 0;
	      status = sldat->loadData("adc",effCh_trigger,trigger_time,trigger_time);
	      if( status != SD_OK ) return -1;
		
              if( (thesewords & 0xF00000) >> 20 != 0xa ){
                  fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] EVENT TRAILER NOT FOUND\n", __LINE__);
                  return -1;
              }

          }

          // Filler words to end
          thesewords = p[jj++] & 0xFFFFFF;
          while( thesewords == 0xe00000 ){
              thesewords = p[jj++] & 0xFFFFFF;
          }

          //printf("BLOCK TRAILER %06x\n", thesewords);
          //printf("Good? (2)       %x\n", (thesewords & 0xF00000) >> 20);
          //printf("NWORDS IN BLOCK %d\n", (thesewords & 0x0FFFFF) >> 0);
          if( (thesewords & 0xF00000) >> 20 != 0x2 ){
              fprintf(stderr, "[ERROR  MPDModule::LoadSlot, line %d] BLOCK TRAILER NOT FOUND\n", __LINE__);
              return -1;
          }
      }

      return fWordsSeen;
  }
  
  Int_t MPDModule::GetData(Int_t adc, Int_t sample, Int_t chan) const {
    printf("MPD GET DATA\n");
    Int_t idx = asc2i(adc, sample, chan);
    if ((idx < 0 ) || (idx >= fNumChan*fNumSample*fNumADC)) { return 0; }
    return fData[idx];
  }
  
  void MPDModule::Clear(const Option_t *) {
    fNumHits = 0;
    for (Int_t i=0; i<fNumChan*fNumSample*fNumADC; i++) fData[i]=0;
    for (Int_t i=0; i<fNumADC*fNumSample; i++) { 
      fFrameHeader[i]=0;
      fFrameTrailer[i]=0;
    }
    
  }
  
  Int_t MPDModule::Decode(const UInt_t *) {
    
    
    return 0;
  }


}

ClassImp(Decoder::MPDModule)
