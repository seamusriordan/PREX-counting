#include <iostream>
#include "SBSGEMPlane.h"
#include "TDatime.h"
#include "THaEvData.h"

SBSGEMPlane::SBSGEMPlane( const char *name, const char *description,
    THaDetectorBase* parent ):
    THaSubDetector(name,description,parent),
    fNch(0),fStrip(NULL),fPedestal(NULL),fcommon_mode(NULL),
    trigger_time(-1),ev_num(-1)
{
    // FIXME:  To database
    fZeroSuppress    = kFALSE;
    fZeroSuppressRMS = 5.0;

    for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
      fadc[i] = NULL;
    }
    fadc_sum = NULL;
    return;
}

SBSGEMPlane::~SBSGEMPlane() {
    if( fStrip ){
        fadc0 = NULL;
        fadc1 = NULL;
        fadc2 = NULL;
        fadc3 = NULL;
        fadc4 = NULL;
        fadc5 = NULL;
	fadc_sum = NULL;
        for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
            delete fadc[i];
            fadc[i] = NULL;
        }
	delete fcommon_mode;
	fcommon_mode = NULL;
        delete fPedestal;
        fPedestal = NULL;
        delete fStrip;
        fStrip = NULL;
	
    }

    return;
}

Int_t SBSGEMPlane::ReadDatabase( const TDatime& date ){
    std::cout << "[SBSGEMPlane::ReadDatabase]" << std::endl;

    Int_t status;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    std::vector<Double_t> rawped;
    std::vector<Double_t> rawrms;

    const DBRequest request[] = {
        { "chanmap",        &fChanMapData,        kIntV, 0, 0},
        { "ped",            &rawped,        kDoubleV, 0, 1},
        { "rms",            &rawrms,        kDoubleV, 0, 1},
        {}
    };
    status = LoadDB( file, date, request, fPrefix );
    fclose(file);

    Int_t nentry = fChanMapData.size()/MPDMAP_ROW_SIZE;
    for( Int_t mapline = 0; mapline < nentry; mapline++ ){
        mpdmap_t thisdata;
        thisdata.crate  = fChanMapData[0+mapline*MPDMAP_ROW_SIZE];
        thisdata.slot   = fChanMapData[1+mapline*MPDMAP_ROW_SIZE];
        thisdata.mpd_id = fChanMapData[2+mapline*MPDMAP_ROW_SIZE];
        thisdata.gem_id = fChanMapData[3+mapline*MPDMAP_ROW_SIZE];
        thisdata.adc_id = fChanMapData[4+mapline*MPDMAP_ROW_SIZE];
        thisdata.i2c    = fChanMapData[5+mapline*MPDMAP_ROW_SIZE];
        thisdata.pos    = fChanMapData[6+mapline*MPDMAP_ROW_SIZE];
        thisdata.invert = fChanMapData[7+mapline*MPDMAP_ROW_SIZE];
        fMPDmap.push_back(thisdata);
    }

    std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

    // FIXME:  make sure to delete if already initialized
    fStrip    = new Int_t [N_APV25_CHAN*nentry];
    
    for( Int_t i = 0; i < N_MPD_TIME_SAMP; i++ ){
        fadc[i] = new Int_t [N_APV25_CHAN*nentry];
        for( Int_t j = 0; j < N_MPD_TIME_SAMP; j++ ){
            fadc[i][j] = 0.0;
        }
    }
    fadc0 = fadc[0];
    fadc1 = fadc[1];
    fadc2 = fadc[2];
    fadc3 = fadc[3];
    fadc4 = fadc[4];
    fadc5 = fadc[5];

    fadc_sum = new Int_t[N_APV25_CHAN*nentry];
    fcommon_mode = new Int_t[N_APV25_CHAN*nentry];

    fPedestal = new Double_t [N_APV25_CHAN*nentry];
    fRMS      = new Double_t [N_APV25_CHAN*nentry];

    for( Int_t i = 0; i < N_APV25_CHAN*nentry; i++ ){
      fPedestal[i] = 0.0;
      fRMS[i] = 0.0;
    }


    for( UInt_t i = 0; i < rawped.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        int idx = (int) rawped[i];
	
    	if( idx < N_APV25_CHAN*nentry ){
    		fPedestal[idx] = rawped[i+1];
    	} else {
		
    	    std::cout << "[SBSGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
    	}
    }

    for( UInt_t i = 0; i < rawrms.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        int idx = (int) rawrms[i];
	if( idx < N_APV25_CHAN*nentry ){
		fRMS[idx] = rawrms[i+1];
	} else {
	    std::cout << "[SBSGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
	}
    }


    return 0;
}

Int_t SBSGEMPlane::DefineVariables( EMode mode ) {
    if( mode == kDefine and fIsSetup ) return kOK;
      fIsSetup = ( mode == kDefine );

      RVarDef vars[] = {
          { "nch",   "Number of channels",   "fNch" },
          { "strip", "Strip number mapping", "fStrip" },
          { "adc0", "ADC sample", "fadc0" },
          { "adc1", "ADC sample", "fadc1" },
          { "adc2", "ADC sample", "fadc2" },
          { "adc3", "ADC sample", "fadc3" },
          { "adc4", "ADC sample", "fadc4" },
          { "adc5", "ADC sample", "fadc5" },
          { "adc_sum", "ADC samples sum", "fadc_sum" },
          { "common_mode", "Common Mode", "fcommon_mode" },
	  { "trigger_time", "Trigger Time", "trigger_time" },
	  { "ev_num","event counter","ev_num"},
          { 0 },
      };


      Int_t ret = DefineVarsFromList( vars, mode );

      if( ret != kOK )
          return ret;

      return kOK;

}

void    SBSGEMPlane::Clear( Option_t* opt){
    fNch = 0;
    return;
}

Int_t   SBSGEMPlane::Decode( const THaEvData& evdata ){
//    std::cout << "[SBSGEMPlane::Decode " << fName << "]" << std::endl;

    int i;

    fNch = 0;
    for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin() ; it != fMPDmap.end(); ++it){

        // Find channel for trigger time first
      Int_t effChan = it->mpd_id << 5 ;  // Channel reserved for trigger time
	ULong_t coarse_time1 = evdata.GetData(it->crate,it->slot,effChan,0);
	UInt_t coarse_time2 = evdata.GetData(it->crate,it->slot,effChan,1);
	UInt_t fine_time = evdata.GetData(it->crate,it->slot,effChan,2);
	trigger_time = ((coarse_time1<<20)|coarse_time2)+fine_time/6.0;
	
	effChan = it->mpd_id<<4;
	ev_num = evdata.GetData(it->crate,it->slot,effChan,0);
	
	// Start reading data sample
        effChan = it->mpd_id << 8 | it->adc_id;
        // Find channel for this crate/slot

        Int_t nchan = evdata.GetNumChan( it->crate, it->slot );

//        printf("nchan = %d\n", nchan );

        for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
            Int_t chan = evdata.GetNextChan( it->crate, it->slot, ichan );
            if( chan != effChan ) continue; // not part of this detector


            Int_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );

            //std::cout << fName << " MPD " << it->mpd_id << " ADC " << it->adc_id << " found " << nsamp << std::endl;
            //std::cout << nsamp << " samples detected (" << nsamp/N_APV25_CHAN <<  ")" << std::endl;

            assert( nsamp == N_APV25_CHAN*N_MPD_TIME_SAMP );
	    
	    Double_t arrADCSum[128]; // Copy of ADC sum for CMN Calculation
	    Int_t arrfNch[128]; // Copy of fNch for CMN Calculation
            for( Int_t strip = 0; strip < N_APV25_CHAN; ++strip ) {
                // data is packed like this
                // [ts1 of 128 chan] [ts2 of 128chan] ... [ts6 of 128chan]
                
                // Madness....   copy pasted from stand alone decoder
                // I bet there's a more elegant way to express this

                // Int_t RstripNb= 32*(strip%4)+ 8*(int)(strip/4)- 31*(int)(strip/16);
                // RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);
	      
	        // New: Use a pre-computed array from Danning to skip the above two steps.
		Int_t RstripNb = APVMAP[strip];
                RstripNb=RstripNb+(127-2*RstripNb)*it->invert;
                Int_t RstripPos = RstripNb + 128*it->pos;

/*
                if( it->adc_id == 10 ){
                std::cout << "ADC " << it->adc_id << " final strip pos: " << RstripPos << std::endl;
                }
*/

                fStrip[fNch] = RstripPos;
		
		fadc_sum[fNch] = 0;
		
                for( Int_t adc_samp = 0; adc_samp < N_MPD_TIME_SAMP; adc_samp++ ){
                    int isamp = adc_samp*N_APV25_CHAN + strip;

                    assert(isamp < nsamp);

                    fadc[adc_samp][fNch] =  evdata.GetData(it->crate, it->slot, chan, isamp) -
                                            fPedestal[RstripPos];
		    fadc_sum[fNch] += fadc[adc_samp][fNch];
                    assert( ((UInt_t) fNch) < fMPDmap.size()*N_APV25_CHAN ); 
		    // Note fMPDmap.size() equals to Number of APV Cards
                }
		// copy adc sum and its fNCH
		arrADCSum[strip] = fadc_sum[fNch];
		arrfNch[strip] = fNch;
                // Zero suppression
                if( !fZeroSuppress ||  
                      ( fRMS[RstripPos] > 0.0 && fabs(fadc[2][fNch])/fRMS[RstripPos] > fZeroSuppressRMS ) ){
                    fNch++;
                }
            }// End Strip Loop

	    // Common Mode Calculation starts after strip loop -- TY
	    // I calculate common mode noise within each APV, that is 128 channels
	    // Insertion sort arrADCSum
	    Double_t swap_buff;
	    for(Int_t i =1; i<N_APV25_CHAN;++i){
	      swap_buff = arrADCSum[i];
	      Int_t j = i-1;
	      while( j>=0 && arrADCSum[j]>swap_buff ){
		arrADCSum[j+1] = arrADCSum[j];
		j = j-1;
	      }
	      arrADCSum[j+1] = swap_buff;
	    }
	    // Average the channels with middle 1/3 signals
	    Double_t cm_noise = 0;
	    Int_t n_cmn = N_APV25_CHAN/3;
	    for(Int_t strip =n_cmn; strip<2*n_cmn;++strip){
	      if(arrADCSum[strip+1]<arrADCSum[strip]) // Unless N_CMN_CHAN !=128
		std::cout << "Sorting went Wrong ! " << std::endl;
	      cm_noise += arrADCSum[strip];
	    }
	    cm_noise = cm_noise/ n_cmn / N_MPD_TIME_SAMP; // averaged to each sample

	    // Write to fcommon_mode[fNch] array
	    for(Int_t strip=0; strip<N_APV25_CHAN;++strip){
	      fcommon_mode[ arrfNch[strip] ] = cm_noise;  
	      // Defined as the same for all channels in one APV
	    }

        }// End ichan loop: nchan = total APVs 

    }

//    std::cout << fName << " channels found  " << fNch << std::endl;

    return 0;
}

void    SBSGEMPlane::Print( Option_t* opt) const{
    return;
}

Int_t   SBSGEMPlane::Begin( THaRunBase* r){
    return 0;
}

Int_t   SBSGEMPlane::End( THaRunBase* r){
    return 0;
}

