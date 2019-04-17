#include <iostream>
#include "MPDGEMPlane.h"
#include "TDatime.h"
#include "THaEvData.h"
#include <vector>



MPDGEMPlane::MPDGEMPlane( const char *name, const char *description,
    THaDetectorBase* parent ):
    GEMPlane(name,description,parent),
//    nch(0),fStrip(NULL),fPedestal(NULL),fcommon_mode(NULL),
    trigger_time(-1),ev_num(-1)
{
    // FIXME:  To database
    fZeroSuppress    = kFALSE;
    fZeroSuppressRMS = 5.0;

    fMaxSamp = N_MPD_TIME_SAMP;

    for( UInt_t i = 0; i < fMaxSamp; i++ ){
      fADCForm[i] = NULL;
    }
    fADCSum = NULL;


    return;
}

MPDGEMPlane::~MPDGEMPlane() {
    if( fSigStrips.size() > 0 ){
        fADC0 = NULL;
        fADC1 = NULL;
        fADC2 = NULL;
        fADC3 = NULL;
        fADC4 = NULL;
        fADC5 = NULL;
	fADCSum = NULL;
        for( UInt_t i = 0; i < fMaxSamp; i++ ){
            delete fADCForm[i];
            fADCForm[i] = NULL;
        }
        /*
	delete fcommon_mode;
	fcommon_mode = NULL;
        delete fPedestal;
        fPedestal = NULL;
        delete fStrip;
        fStrip = NULL;
        */
	
    }

    return;
}

Int_t MPDGEMPlane::ReadDatabase( const TDatime& date ){
    std::cout << "[MPDGEMPlane::ReadDatabase]" << std::endl;

    // Read the database for the base class, quit if error
    Int_t status = ReadDatabaseCommon(date);
    if( status != kOK )
        return status;

    Int_t err;

    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    TString mapping;
    fMaxClusterSize = kMaxUInt;
    fMinAmpl   = 0.0;
    fSplitFrac = 0.0;
    fMapType   = kOneToOne;
    fMaxSamp   = 1;
    fChanMap.clear();
    fPed.clear();
    fAmplSigma = 0.36; // default, an educated guess


    std::vector<Double_t> rawped;
    std::vector<Double_t> rawrms;

    Int_t gbl = GetDBSearchLevel(fPrefix);

    const DBRequest request[] = {
        { "chanmap",        &fChanMapData,  kIntV,    0, 0},
        { "ped",            &rawped,        kDoubleV, 0, 1},
        { "rms",            &rawrms,        kDoubleV, 0, 1},

        { "strip.pos",      &fStart },
        { "strip.pitch",    &fPitch,          kDouble,  0, 0, gbl },
        { "maxclustsiz",    &fMaxClusterSize, kUInt,    0, 1, gbl },
        { "maxsamp",        &fMaxSamp,        kUInt,    0, 1, gbl },
        { "adc.min",        &fMinAmpl,        kDouble,  0, 1, gbl },
        { "split.frac",     &fSplitFrac,      kDouble,  0, 1, gbl },
        { "mapping",        &mapping,         kTString, 0, 1, gbl },

        {}
    };
    err = LoadDB( file, date, request, fPrefix );

    if( err ) return err;

    fclose(file);

    UInt_t nentry = fChanMapData.size()/MPDMAP_ROW_SIZE;
    for( UInt_t mapline = 0; mapline < nentry; mapline++ ){
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

    fNelem = N_APV25_CHAN*fMPDmap.size();

    SafeDelete(fADCraw);
    SafeDelete(fADC);
    SafeDelete(fHitTime);
    SafeDelete(fADCcor);
    SafeDelete(fGoodHit);

    std::cout << fName << " mapped to " << nentry << " APV25 chips" << std::endl;

    for( UInt_t i = 0; i < fMaxSamp; i++ ){
        fADCForm[i] = new Int_t [fNelem];
        for( UInt_t j = 0; j < fMaxSamp; j++ ){
            fADCForm[i][j] = 0.0;
        }
    }
    fADC0 = fADCForm[0];
    fADC1 = fADCForm[1];
    fADC2 = fADCForm[2];
    fADC3 = fADCForm[3];
    fADC4 = fADCForm[4];
    fADC5 = fADCForm[5];

    fADCSum = new Int_t[fNelem];
//    fcommon_mode = new Int_t[N_APV25_CHAN*nentry];

    fPed.clear();
    fRMS.clear();

    for( UInt_t i = 0; i < rawped.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        UInt_t idx = (UInt_t) rawped[i];
	
    	if( idx < (UInt_t) fNelem && fPed.size() == idx ){
    		fPed.push_back(rawped[i+1]);
    	} else {
		
    	    std::cout << "[MPDGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
    	}
    }

    for( UInt_t i = 0; i < rawrms.size(); i++ ){
        if( (i % 2) == 1 ) continue;
        UInt_t idx = (UInt_t) rawrms[i];
	if( idx < ((UInt_t) fNelem) && fRMS.size() == idx ){
		fRMS.push_back(rawrms[i+1]);
	} else {
	    std::cout << "[MPDGEMPlane::ReadDatabase]  WARNING: " << " strip " << idx  << " listed but not enough strips in cratemap" << std::endl;
	}
    }

    fADCraw = new Float_t[fNelem];
    fADC = new Float_t[fNelem];
    fHitTime = new Float_t[fNelem];
    fADCcor = new Float_t[fNelem];
    fGoodHit = new Byte_t[fNelem];
    fSigStrips.reserve(fNelem);
    fStripsSeen.resize(fNelem);




    return 0;
}

Int_t MPDGEMPlane::DefineVariables( EMode mode ) {
    if( mode == kDefine and fIsSetup ) return kOK;
      fIsSetup = ( mode == kDefine );

      RVarDef vars[] = {
          { "nrawstrips",     "nstrips with decoder data",        "fNrawStrips" },
          { "nhitstrips",     "nstrips > 0",                      "fNhitStrips" },
          { "nstrips",        "Num strips with hits > adc.min",   "GetNsigStrips()" },
          { "hitocc",         "strips > 0 / n_all_strips",        "fHitOcc" },
          { "occupancy",      "nstrips / n_all_strips",           "fOccupancy" },
          { "strip.adcraw",   "Raw strip ADC sum",                "fADCraw" },
          { "strip.adc",      "Deconvoluted strip ADC sum",       "fADC" },
          { "strip.adc_c",    "Pedestal-sub strip ADC sum",       "fADCcor" },
          { "strip.time",     "Leading time of strip signal (ns)","fHitTime" },
          { "strip.good",     "Good pulse shape on strip",        "fGoodHit" },
          { "nhits",          "Num hits (clusters of strips)",    "GetNhits()" },
          { "noise",          "Noise level (avg below adc.min)",  "fDnoise" },
          { "ncoords",        "Num fit coords",                   "GetNcoords()" },
          { "coord.pos",      "Position used in fit (m)",         "fFitCoords.TreeSearch::FitCoord.fPos" },
          { "coord.trkpos",   "Track pos from projection fit (m)","fFitCoords.TreeSearch::FitCoord.fTrackPos" },
          { "coord.trkslope", "Track slope from projection fit",  "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
          { "coord.resid",    "Residual of trkpos (m)",           "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
          { "coord.3Dpos",    "Crossing position of fitted 3D track (m)", "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
          { "coord.3Dresid",  "Residual of 3D trkpos (m)",        "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
          { "coord.3Dslope",  "Slope of fitted 3D track wrt projection",  "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },


          //{ "nch",   "Number of channels",   "nch" },
          { "strip.number", "Strip number mapping", "fSigStrips" },
          { "adc0", "ADC sample", "fADC0" },
          { "adc1", "ADC sample", "fADC1" },
          { "adc2", "ADC sample", "fADC2" },
          { "adc3", "ADC sample", "fADC3" },
          { "adc4", "ADC sample", "fADC4" },
          { "adc5", "ADC sample", "fADC5" },
          { "adc_sum", "ADC samples sum", "fADCSum" },
          //{ "common_mode", "Common Mode", "fcommon_mode" },
          { "trigger_time", "Trigger Time", "trigger_time" },
          { "ev_num","event counter","ev_num"},
          { 0 },
      };


      Int_t ret = DefineVarsFromList( vars, mode );

      if( ret != kOK )
          return ret;


      RVarDef nonmcvars[] = {
          { "hit.pos",  "Hit centroid (m)",      "fHits.TreeSearch::GEMHit.fPos" },
          { "hit.adc",  "Hit ADC sum",           "fHits.TreeSearch::GEMHit.fADCsum" },
          { "hit.size", "Num strips ",           "fHits.TreeSearch::GEMHit.fSize" },
          { "hit.type", "Hit analysis result",   "fHits.TreeSearch::GEMHit.fType" },
          { 0 }
      };
      ret = DefineVarsFromList( nonmcvars, mode );


      return kOK;

}

void MPDGEMPlane::Clear( Option_t* opt ){

    TreeSearch::GEMPlane::Clear(opt);
    return;
}

Int_t MPDGEMPlane::Decode( const THaEvData& evdata ){
//    std::cout << "[MPDGEMPlane::Decode " << fName << "]" << std::endl;

    Int_t nch = 0;
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


            UInt_t nsamp = evdata.GetNumHits( it->crate, it->slot, chan );

            //std::cout << fName << " MPD " << it->mpd_id << " ADC " << it->adc_id << " found " << nsamp << std::endl;
            //std::cout << nsamp << " samples detected (" << nsamp/N_APV25_CHAN <<  ")" << std::endl;

            assert( nsamp == N_APV25_CHAN*fMaxSamp );
	    
	    Double_t arrADCSum[128]; // Copy of ADC sum for CMN Calculation

            Int_t nchStartOfAPV = nch;
            for( Int_t strip = 0; strip < N_APV25_CHAN; ++strip ) {
                // data is packed like this
                // [ts1 of 128 chan] [ts2 of 128chan] ... [ts6 of 128chan]
	      
                Int_t RstripPos = GetRStripNumber( strip, it->pos, it->invert );

//                fStrip[nch] = RstripPos;
		
		fADCSum[nch] = 0;

                Vflt_t samples;
                samples.clear();
		
                for( UInt_t adc_samp = 0; adc_samp < fMaxSamp; adc_samp++ ){
                    UInt_t isamp = adc_samp*N_APV25_CHAN + strip;

                    assert(isamp < nsamp);

                    Int_t rawadc =  evdata.GetData(it->crate, it->slot, chan, isamp);
                    fADCForm[adc_samp][nch] = rawadc - fPed[RstripPos];
		    fADCSum[nch] += fADCForm[adc_samp][nch];
                    assert( nch < fNelem ); 

                    samples.push_back((Float_t) rawadc);
		    // Note fMPDmap.size() equals to Number of APV Cards
                }

                MPDStripData_t stripdata = ChargeDep(samples);

		// copy adc sum and its fNCH
		arrADCSum[strip] = fADCSum[nch];

                ++fNrawStrips;
                ++fNhitStrips;

                fADCraw[RstripPos]  = stripdata.adcraw;
                fADC[RstripPos]     = stripdata.adc;
                fHitTime[RstripPos] = stripdata.time;
                fGoodHit[RstripPos] = stripdata.pass;

                fADCcor[RstripPos]  = stripdata.adc;


                // Zero suppression
                if( !fZeroSuppress ||  
                       ( fRMS[RstripPos] > 0.0 && fabs(fADCForm[2][nch])/fRMS[RstripPos] > fZeroSuppressRMS) ){
                    fSigStrips.push_back(RstripPos);
                    nch++;
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
	    Int_t n_cutoff = 20 ;
	    Int_t n_cmn = N_APV25_CHAN - 2*n_cutoff;
	    for(Int_t strip =n_cutoff; strip< N_APV25_CHAN -n_cutoff ;++strip){
	      if(arrADCSum[strip+1]<arrADCSum[strip]) // Unless N_CMN_CHAN !=128
		std::cout << "Sorting went Wrong ! " << std::endl;
	      cm_noise += arrADCSum[strip];
	    }
	    cm_noise = cm_noise/ n_cmn / fMaxSamp; // averaged to each sample
            fDnoise = cm_noise;

            // Correct superclass' variables for common mode noise
	    for(Int_t strip=0; strip<N_APV25_CHAN;++strip){
                UInt_t RstripPos = GetRStripNumber( strip, it->pos, it->invert );
                fADCcor[RstripPos] -= fPed[RstripPos] + cm_noise;
	    }

            // Correct this class' variables for common mode noise
	    for(Int_t ch= nchStartOfAPV; ch < nch; ch++){
                for( UInt_t adc_samp = 0; adc_samp < fMaxSamp; adc_samp++ ){
                    fADCForm[adc_samp][ch] -= cm_noise;
                }
                fADCSum[nch] -= cm_noise*fMaxSamp;
            }

        }// End ichan loop: nchan = total APVs 

    }

    fHitOcc    = static_cast<Double_t>(fNhitStrips) / fNelem;
    fOccupancy = static_cast<Double_t>(GetNsigStrips()) / fNelem;

//    std::cout << fName << " channels found  " << nch << std::endl;

    return 0;
}

Int_t MPDGEMPlane::GetRStripNumber( UInt_t strip, UInt_t pos, UInt_t invert ){
    Int_t RstripNb = APVMAP[strip];
    RstripNb=RstripNb+(127-2*RstripNb)*invert;
    Int_t RstripPos = RstripNb + 128*pos;

    return RstripPos;
}

Int_t   MPDGEMPlane::Begin( THaRunBase* run ){
    TreeSearch::GEMPlane::Begin(run);
    return 0;
}

Int_t   MPDGEMPlane::End( THaRunBase* run ){
    TreeSearch::GEMPlane::End(run);
    return 0;
}

MPDStripData_t MPDGEMPlane::ChargeDep( const std::vector<Float_t>& amp ) {
  // Deconvolute signal given by samples in 'amp', return approximate integral.
  // Currently analyzes exactly 3 samples.
  // From Kalyan Allada
  // NIM A326, 112 (1993)

  //FIXME: from database, proper value for Tp
  const Float_t delta_t = 25.0; // time interval between samples (ns)
  const Float_t Tp      = 50.0; // RC filter time constant (ns)

  assert( amp.size() >= 3 );

  Float_t adcraw = delta_t*(amp[0]+amp[1]+amp[2]);

  // Weight factors calculated based on the response of the silicon microstrip
  // detector:
  // v(t) = (delta_t/Tp)*exp(-delta_t/Tp)
  // Need to update this for GEM detector response(?):
  // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
  // where A is the amplitude, t0 the begin of the rise, tau1 the time
  // parameter for the rising edge and tau2 the for the falling edge.

  Float_t x = delta_t/Tp;

  Float_t w1 = TMath::Exp(x-1)/x;
  Float_t w2 = -2*TMath::Exp(-1)/x;
  Float_t w3 = TMath::Exp(-x-1)/x;

  // Deconvoluted signal samples, assuming measurements of zero before the
  // leading edge
  Float_t sig[3] = { amp[0]*w1,
		     amp[1]*w1+amp[0]*w2,
		     amp[2]*w1+amp[1]*w2+amp[0]*w3 };

  Float_t adc    = delta_t*(sig[0]+sig[1]+sig[2]);
  Float_t time   = 0;     // TODO

  Bool_t pass;
  // Calculate ratios for 3 samples and check for bad signals
  if( amp[2] > 0 ) {
    Float_t r1 = amp[0]/amp[2];
    Float_t r2 = amp[1]/amp[2];
    pass = (r1 < 1.0 and r2 < 1.0 and r1 < r2);
  } else
    pass = false;

  return MPDStripData_t(adcraw,adc,time,pass);
}


