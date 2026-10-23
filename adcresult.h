//#include "TRoot.h"
//#include <string.h>


constexpr int MEAN = 0;
constexpr int VARI = 1;
constexpr int RESO = 2;
constexpr int STDDEV = 3;
constexpr int MPV = 4;
constexpr int FWHM = 5;
constexpr int n_chara = 6;


constexpr int VAL = 0;
constexpr int ERR = 1;

constexpr int gHV = 0;
constexpr int gBX = 1;
constexpr int gBY = 2;
constexpr int gBZ = 3;
constexpr int gLD = 4;
constexpr int gTHETA = 5;
constexpr int gPHI = 6;
constexpr int gALPHA = 7;
constexpr int gN_COND = 8;

class FitResultSPE
{
   public:
  //Setting
  Int_t    run;
  Int_t    pmtnum;
  Double_t theta;
  Double_t phi;
  Double_t distance;
  Float_t  nicf_x;
  Float_t  nicf_y;
  Float_t  nicf_z;
  Int_t    pmt_x;
  Int_t    pmt_y;
  Int_t    pmt_z;
  //std::string   serial;
  Int_t    ch;
  char serial[16];
  Float_t  HVS;
  Double_t unixtime;
  Int_t    LDA;
  Int_t    LDW;

  Double_t cond[gN_COND][2];
  Double_t qdc[n_chara][2];
  Double_t tdc[n_chara][2];

  //Result
  Double_t h_mean[2]      ;
  Double_t f_mean[2]      ;
  Double_t f2_mean[2]      ;
  Double_t h_variance[2]  ;
  Double_t f_variance[2]  ;
  Double_t f2_variance[2]  ;
  Double_t h_stddev[2]    ;
  Double_t f_stddev[2]    ;
  Double_t f2_stddev[2]    ;
  Double_t h_resolution[2]    ;
  Double_t f_resolution[2]    ;
  Double_t f2_resolution[2]    ;
  Double_t meanpredict      ;
  Double_t variancepredict  ;
  Double_t peakx ;
  Double_t peaky ;
  Double_t FWHMlow  ;
  Double_t FWHMhigh ;
  Double_t FWHM ;
  Double_t sigmalow  ;
  Double_t sigmahigh ;
  Double_t sigma           ;
  Double_t peak1pex        ;
  Double_t peak1pey        ;
  Double_t peak1pesigma    ;
  Double_t bssigma         ;
  Double_t peak1pexerr     ;
  Double_t peak1peyerr     ;
  Double_t peak1pesigmaerr ;
  Double_t valley          ;
  Double_t valleyx         ;
  Double_t peak            ;
  Double_t pv              ;
  Double_t area1pe         ;
  Double_t area1peall      ;
  Double_t areabs          ;
  Double_t areapredict     ;
  Double_t areasigmalow    ;
  Double_t areasigmahigh   ;
  Double_t areasigma       ;
  Double_t gainpeakx       ;
  Double_t gainpeak1pex    ;
  Double_t gainpeak1pexerr ;

  Double_t chi2;
  Double_t ndf;

  Double_t tmean      ;
  Double_t tvariance  ;
  Double_t tpeakx ;
  Double_t tpeaky ;
  Double_t tFWHMlow  ;
  Double_t tFWHMhigh ;
  Double_t tEdgelow  ;
  Double_t tEdgehigh ;
  Double_t tFWHM ;
  Double_t tgamma ;
  Double_t tmu ;
  Double_t tsigma ;
  Double_t tlambda ;
  Double_t ttau ;
  Double_t tgammaerr ;
  Double_t tmuerr ;
  Double_t tsigmaerr ;
  Double_t tlambdaerr ;
  Double_t ttauerr ;
  Double_t tchi2;
  Double_t tndf;

  Long64_t mppchn;
  Double_t pmthn;
  Double_t pmtftn;
  Double_t pmtfqn;
  Double_t mppcn[2];  // mppc's number of single and double photoelectron 
  Double_t mppcnpe;   //mppc's number of photoelectron
  Double_t pmtnpegain; //pmt's number of photoelectron \times gain
  Double_t pmtspegain; //pmt's number of photoelectron \times gain

  Double_t noiserate[2];
  
};

int result_fill(TString out_f_name, FitResultSPE res, const std::vector<TCanvas*> can = {0}){
    TFile* out_f = TFile::Open(out_f_name,"update");
    if (!out_f){
		std::cout<< "Error opening file" << std::endl;
		exit(-1);
    }
    if(can[0]){
	    for ( int i = 0; i<can.size(); i++ ){
		    can[i]->Write();
		    std::cout << "canvas[" << i << "] saved." << std::endl;
	    }
    }

    const char *adctrname = "adcr";
    FitResultSPE *pres = &res;
    TTree *adctr = 0;
    TBranch * b_res;
    TObject *tree_obj = out_f->Get(adctrname);
    if (tree_obj) {
       adctr = (TTree*)tree_obj;
       adctr->SetBranchAddress("res", &pres, &b_res );
    } else {
       adctr = new TTree(adctrname,"ADC SPE Result");
       adctr->Branch("res","FitResultSPE", &res);
    }
    adctr->Fill();
    adctr->Write(adctrname,TObject::kOverwrite);
    if (!tree_obj) { delete adctr; }
    out_f->Close();
    return 0;
}










// Below is trash
// Below is trash
// Below is trash
// Below is trash
// Below is trash
// Below is trash


class ValueWithError
{
public:
	Double_t val;
	Double_t err;
};
class DistChara
{
public:
  ValueWithError mean;
  ValueWithError variance;
  ValueWithError stddev;
  ValueWithError resolution;
  ValueWithError mpv;
  ValueWithError fwhm;
};


int adcresult_fill(TFile *rootf, FitResultSPE res)
{
    const char *adctrname = "adcr";
    //std::vector<std::string> serialstrings;
    //    serialstrings.push_back(str);

    TTree *adctr = 0;
    TBranch        *b_run            ;
    TBranch        *b_pmtnum         ;
    TBranch        *b_theta          ;
    TBranch        *b_phi            ;
    TBranch        *b_distance       ;
    TBranch        *b_nicf_x         ;
    TBranch        *b_nicf_y         ;
    TBranch        *b_nicf_z         ;
    TBranch        *b_pmt_x          ;
    TBranch        *b_pmt_y          ;
    TBranch        *b_pmt_z          ;
    TBranch        *b_ch             ;
    //TBranch        *b_nserial        ;
    TBranch        *b_serial         ;
    TBranch        *b_HVS            ;
    TBranch        *b_unixtime       ;
    TBranch        *b_LDA            ;
    TBranch        *b_LDW            ;
    TBranch        *b_h_mean         ;
    TBranch        *b_f_mean         ;
    TBranch        *b_f2_mean        ;
    TBranch        *b_h_variance     ;
    TBranch        *b_f_variance     ;
    TBranch        *b_f2_variance    ;
    TBranch        *b_h_stddev       ;
    TBranch        *b_f_stddev       ;
    TBranch        *b_f2_stddev      ;
    TBranch        *b_h_resolution   ;
    TBranch        *b_f_resolution   ;
    TBranch        *b_f2_resolution  ;
    TBranch        *b_meanpredict    ;
    TBranch        *b_variancepredict;
    TBranch        *b_peakx          ;
    TBranch        *b_peaky          ;
    TBranch        *b_FWHMlow        ;
    TBranch        *b_FWHMhigh       ;
    TBranch        *b_FWHM           ;
    TBranch        *b_sigmalow       ;
    TBranch        *b_sigmahigh      ;
    TBranch        *b_sigma          ;
    TBranch        *b_peak1pex       ;
    TBranch        *b_peak1pey       ;
    TBranch        *b_peak1pesigma   ;
    TBranch        *b_bssigma        ;
    TBranch        *b_peak1pexerr    ;
    TBranch        *b_peak1peyerr    ;
    TBranch        *b_peak1pesigmaerr;
    TBranch        *b_valley         ;
    TBranch        *b_valleyx        ;
    TBranch        *b_peak           ;
    TBranch        *b_pv             ;
    TBranch        *b_area1pe        ;
    TBranch        *b_area1peall     ;
    TBranch        *b_areabs         ;
    TBranch        *b_areapredict    ;
    TBranch        *b_areasigmalow   ;
    TBranch        *b_areasigmahigh  ;
    TBranch        *b_areasigma      ;
    TBranch        *b_gainpeakx      ;
    TBranch        *b_gainpeak1pex   ;
    TBranch        *b_gainpeak1pexerr;
    TBranch        *b_chi2;
    TBranch        *b_ndf;
    //TBranch        *b_length;
    //TBranch        *b_description;
    //TBranch        *b_start;
    //TBranch        *b_serial;
    //TBranch        *b_nhv;
    //TBranch        *b_hv;
    //TBranch        *b_hvi;
    //TBranch        *b_hvs;
    //TBranch        *b_nthr;
    //TBranch        *b_thr;
    //TBranch        *b_lsfactor;
    //TBranch        *b_lsamp;
    //TBranch        *b_lswidth;
    //TBranch        *b_nadc;
    //TBranch        *b_gaincorr;
    //TBranch        *b_pedestal;
    //TBranch        *b_ntdc;
    //TBranch        *b_delay;
    TBranch        *b_tmean           ;
    TBranch        *b_tvariance       ;
    TBranch        *b_tpeakx          ;
    TBranch        *b_tpeaky          ;
    TBranch        *b_tFWHMlow        ;
    TBranch        *b_tFWHMhigh       ;
    TBranch        *b_tEdgelow        ;
    TBranch        *b_tEdgehigh       ;
    TBranch        *b_tFWHM           ;
    TBranch        *b_tgamma          ;
    TBranch        *b_tmu             ;
    TBranch        *b_tsigma          ;
    TBranch        *b_tlambda         ;
    TBranch        *b_tgammaerr       ;
    TBranch        *b_tmuerr          ;
    TBranch        *b_tsigmaerr       ;
    TBranch        *b_tlambdaerr      ;
    TBranch        *b_tchi2           ;
    TBranch        *b_tndf            ;
    TBranch        *b_mppchn          ;
    TBranch        *b_pmthn           ;
    TBranch        *b_mppcnpe         ;
    TBranch        *b_pmtnpegain      ;
    TBranch        *b_pmtspegain      ;
    TBranch        *b_noiserate       ;

    if (rootf->Get(adctrname)) {
       adctr = (TTree*)(rootf->Get(adctrname));
       adctr->SetBranchAddress("run"            , &(res.run            ),   &b_run            );
       adctr->SetBranchAddress("pmtnum"         , &(res.pmtnum         ),   &b_pmtnum         );
       adctr->SetBranchAddress("theta"          , &(res.theta          ),   &b_theta          );
       adctr->SetBranchAddress("phi"            , &(res.phi            ),   &b_phi            );
       adctr->SetBranchAddress("distance"       , &(res.distance       ),   &b_distance       );
       adctr->SetBranchAddress("nicf_x"         , &(res.nicf_x         ),   &b_nicf_x         );
       adctr->SetBranchAddress("nicf_y"       	, &(res.nicf_y         ),   &b_nicf_y         );
       adctr->SetBranchAddress("nicf_z"       	, &(res.nicf_z         ),   &b_nicf_z         );
       adctr->SetBranchAddress("pmt_x"       	, &(res.pmt_x          ),   &b_pmt_x          );
       adctr->SetBranchAddress("pmt_y"       	, &(res.pmt_y          ),   &b_pmt_y          );
       adctr->SetBranchAddress("pmt_z"       	, &(res.pmt_z          ),   &b_pmt_z          );
       adctr->SetBranchAddress("ch"             , &(res.ch             ),   &b_ch             );
       adctr->SetBranchAddress("serial"         , &(res.serial         ),   &b_serial         );
       adctr->SetBranchAddress("HVS"            , &(res.HVS            ),   &b_HVS            );
       adctr->SetBranchAddress("unixtime"       , &(res.unixtime       ),   &b_unixtime       );
       adctr->SetBranchAddress("LDA"            , &(res.LDA            ),   &b_LDA            );
       adctr->SetBranchAddress("LDW"            , &(res.LDW            ),   &b_LDW            );
       adctr->SetBranchAddress("h_mean"         , &(res.h_mean         ),   &b_h_mean         );
       adctr->SetBranchAddress("f_mean"         , &(res.f_mean         ),   &b_f_mean         );
       adctr->SetBranchAddress("f2_mean"        , &(res.f2_mean        ),   &b_f2_mean        );
       adctr->SetBranchAddress("h_variance"     , &(res.h_variance     ),   &b_h_variance     );
       adctr->SetBranchAddress("f_variance"     , &(res.f_variance     ),   &b_f_variance     );
       adctr->SetBranchAddress("f2_variance"    , &(res.f2_variance    ),   &b_f2_variance    );
       adctr->SetBranchAddress("h_stddev"       , &(res.h_stddev       ),   &b_h_stddev       );
       adctr->SetBranchAddress("f_stddev"       , &(res.f_stddev       ),   &b_f_stddev       );
       adctr->SetBranchAddress("f2_stddev"      , &(res.f2_stddev      ),   &b_f2_stddev      );
       adctr->SetBranchAddress("h_resolution"   , &(res.h_resolution   ),   &b_h_resolution   );
       adctr->SetBranchAddress("f_resolution"   , &(res.f_resolution   ),   &b_f_resolution   );
       adctr->SetBranchAddress("f2_resolution"  , &(res.f2_resolution  ),   &b_f2_resolution  );
       adctr->SetBranchAddress("meanpredict"    , &(res.meanpredict    ),   &b_meanpredict    );
       adctr->SetBranchAddress("variancepredict", &(res.variancepredict),   &b_variancepredict);
       adctr->SetBranchAddress("peakx"          , &(res.peakx          ),   &b_peakx          );
       adctr->SetBranchAddress("peaky"          , &(res.peaky          ),   &b_peaky          );
       adctr->SetBranchAddress("FWHMlow"        , &(res.FWHMlow        ),   &b_FWHMlow        );
       adctr->SetBranchAddress("FWHMhigh"       , &(res.FWHMhigh       ),   &b_FWHMhigh       );
       adctr->SetBranchAddress("FWHM"           , &(res.FWHM           ),   &b_FWHM           );
       adctr->SetBranchAddress("sigmalow"       , &(res.sigmalow       ),   &b_sigmalow       );
       adctr->SetBranchAddress("sigmahigh"      , &(res.sigmahigh      ),   &b_sigmahigh      );
       adctr->SetBranchAddress("sigma"          , &(res.sigma          ),   &b_sigma          );
       adctr->SetBranchAddress("peak1pex"       , &(res.peak1pex       ),   &b_peak1pex       );
       adctr->SetBranchAddress("peak1pey"       , &(res.peak1pey       ),   &b_peak1pey       );
       adctr->SetBranchAddress("peak1pesigma"   , &(res.peak1pesigma   ),   &b_peak1pesigma   );
       adctr->SetBranchAddress("bssigma"        , &(res.bssigma        ),   &b_bssigma        );
       adctr->SetBranchAddress("peak1pexerr"    , &(res.peak1pexerr    ),   &b_peak1pexerr    );
       adctr->SetBranchAddress("peak1peyerr"    , &(res.peak1peyerr    ),   &b_peak1peyerr    );
       adctr->SetBranchAddress("peak1pesigmaerr", &(res.peak1pesigmaerr),   &b_peak1pesigmaerr);
       adctr->SetBranchAddress("valley"         , &(res.valley         ),   &b_valley         );
       adctr->SetBranchAddress("valleyx"        , &(res.valleyx        ),   &b_valleyx        );
       adctr->SetBranchAddress("peak"           , &(res.peak           ),   &b_peak           );
       adctr->SetBranchAddress("pv"             , &(res.pv             ),   &b_pv             );
       adctr->SetBranchAddress("area1pe"        , &(res.area1pe        ),   &b_area1pe        );
       adctr->SetBranchAddress("area1peall"     , &(res.area1peall     ),   &b_area1peall     );
       adctr->SetBranchAddress("areabs"         , &(res.areabs         ),   &b_areabs         );
       adctr->SetBranchAddress("areapredict"    , &(res.areapredict    ),   &b_areapredict    );
       adctr->SetBranchAddress("areasigmalow"   , &(res.areasigmalow   ),   &b_areasigmalow   );
       adctr->SetBranchAddress("areasigmahigh"  , &(res.areasigmahigh  ),   &b_areasigmahigh  );
       adctr->SetBranchAddress("areasigma"      , &(res.areasigma      ),   &b_areasigma      );
       adctr->SetBranchAddress("gainpeakx"      , &(res.gainpeakx      ),   &b_gainpeakx      );
       adctr->SetBranchAddress("gainpeak1pex"   , &(res.gainpeak1pex   ),   &b_gainpeak1pex   );
       adctr->SetBranchAddress("gainpeak1pexerr", &(res.gainpeak1pexerr),   &b_gainpeak1pexerr);
       adctr->SetBranchAddress("chi2"           , &(res.chi2           ),   &b_chi2           );
       adctr->SetBranchAddress("ndf"            , &(res.ndf            ),   &b_ndf            );
       adctr->SetBranchAddress("tmean"           , &(res.tmean          ),  &b_tmean           );
       adctr->SetBranchAddress("tvariance"       , &(res.tvariance      ),  &b_tvariance       );
       adctr->SetBranchAddress("tpeakx"          , &(res.tpeakx         ),  &b_tpeakx          );
       adctr->SetBranchAddress("tpeaky"          , &(res.tpeaky         ),  &b_tpeaky          );
       adctr->SetBranchAddress("tFWHMlow"        , &(res.tFWHMlow       ),  &b_tFWHMlow        );
       adctr->SetBranchAddress("tFWHMhigh"       , &(res.tFWHMhigh      ),  &b_tFWHMhigh       );
       adctr->SetBranchAddress("tEdgelow"        , &(res.tEdgelow       ),  &b_tEdgelow        );
       adctr->SetBranchAddress("tEdgehigh"       , &(res.tEdgehigh      ),  &b_tEdgehigh       );
       adctr->SetBranchAddress("tFWHM"           , &(res.tFWHM          ),  &b_tFWHM           );
       adctr->SetBranchAddress("tgamma"          , &(res.tgamma         ),  &b_tgamma          );
       adctr->SetBranchAddress("tmu"             , &(res.tmu            ),  &b_tmu             );
       adctr->SetBranchAddress("tsigma"          , &(res.tsigma         ),  &b_tsigma          );
       adctr->SetBranchAddress("tlambda"         , &(res.tlambda        ),  &b_tlambda         );
       adctr->SetBranchAddress("tgammaerr"       , &(res.tgammaerr      ),  &b_tgammaerr       );
       adctr->SetBranchAddress("tmuerr"          , &(res.tmuerr         ),  &b_tmuerr          );
       adctr->SetBranchAddress("tsigmaerr"       , &(res.tsigmaerr      ),  &b_tsigmaerr       );
       adctr->SetBranchAddress("tlambdaerr"      , &(res.tlambdaerr     ),  &b_tlambdaerr      );
       adctr->SetBranchAddress("tchi2"           , &(res.tchi2          ),  &b_tchi2           );
       adctr->SetBranchAddress("tndf"            , &(res.tndf           ),  &b_tndf            );
       adctr->SetBranchAddress("mppchn"          , &(res.mppchn         ),  &b_mppchn          );
       adctr->SetBranchAddress("pmthn"           , &(res.pmthn          ),  &b_pmthn           );
       adctr->SetBranchAddress("mppcnpe"         , &(res.mppcnpe        ),  &b_mppcnpe         );
       adctr->SetBranchAddress("pmtnpegain"      , &(res.pmtnpegain     ),  &b_pmtnpegain      );
       adctr->SetBranchAddress("pmtspegain"      , &(res.pmtspegain     ),  &b_pmtspegain      );
       adctr->SetBranchAddress("noiserate"       , &(res.noiserate      ),  &b_noiserate       );
    } else {
       adctr = new TTree(adctrname,"ADC SPE Result");
       adctr->Branch("run"            , &(res.run            ),"run/I"            );
       adctr->Branch("pmtnum"         , &(res.pmtnum         ),"pmtnum/I"         );
       adctr->Branch("theta"          , &(res.theta          ),"theta/D"          );
       adctr->Branch("phi"            , &(res.phi            ),"phi/D"            );
       adctr->Branch("distance"       , &(res.distance       ),"distance/D"       );
       adctr->Branch("nicf_x"         , &(res.nicf_x         ),"nicf_x/F"         );
       adctr->Branch("nicf_y"         , &(res.nicf_y         ),"nicf_y/F"         );
       adctr->Branch("nicf_z"         , &(res.nicf_z         ),"nicf_z/F"         );
       adctr->Branch("pmt_x"          , &(res.pmt_x          ),"pmt_x/I"          );
       adctr->Branch("pmt_y"          , &(res.pmt_y          ),"pmt_y/I"          );
       adctr->Branch("pmt_z"          , &(res.pmt_z          ),"pmt_z/I"          );
       adctr->Branch("ch"             , &(res.ch             ),"ch/I"             );
       adctr->Branch("serial"         , &(res.serial         ),"serial[16]/C"     );
       adctr->Branch("HVS"            , &(res.HVS            ),"HVS/F"            );
       adctr->Branch("unixtime"       , &(res.unixtime       ),"unixtime/D"       );
       adctr->Branch("LDA"            , &(res.LDA            ),"LDA/I"            );
       adctr->Branch("LDW"            , &(res.LDW            ),"LDW/I"            );
       adctr->Branch("h_mean"         , &(res.h_mean         ),"h_mean[2]/D"      );
       adctr->Branch("f_mean"         , &(res.f_mean         ),"f_mean[2]/D"      );
       adctr->Branch("f2_mean"        , &(res.f2_mean        ),"f2_mean[2]/D"     );
       adctr->Branch("h_variance"     , &(res.h_variance     ),"h_variance[2]/D"  );
       adctr->Branch("f_variance"     , &(res.f_variance     ),"f_variance[2]/D"  );
       adctr->Branch("f2_variance"    , &(res.f2_variance    ),"f2_variance[2]/D" );
       adctr->Branch("h_stddev"       , &(res.h_stddev       ),"h_stddev[2]/D"    );
       adctr->Branch("f_stddev"       , &(res.f_stddev       ),"f_stddev[2]/D"    );
       adctr->Branch("f2_stddev"      , &(res.f2_stddev      ),"f2_stddev[2]/D"   );
       adctr->Branch("h_resolution"   , &(res.h_resolution   ),"h_resolution[2]/D");
       adctr->Branch("f_resolution"   , &(res.f_resolution   ),"f_resolution[2]/D");
       adctr->Branch("f2_resolution"  , &(res.f2_resolution  ),"f2_resolution[2]/D");
       adctr->Branch("meanpredict"    , &(res.meanpredict    ),"meanpredict/D"    );
       adctr->Branch("variancepredict", &(res.variancepredict),"variancepredict/D");
       adctr->Branch("peakx"          , &(res.peakx          ),"peakx/D"          );
       adctr->Branch("peaky"          , &(res.peaky          ),"peaky/D"          );
       adctr->Branch("FWHMlow"        , &(res.FWHMlow        ),"FWHMlow/D"        );
       adctr->Branch("FWHMhigh"       , &(res.FWHMhigh       ),"FWHMhigh/D"       );
       adctr->Branch("FWHM"           , &(res.FWHM           ),"FWHM/D"           );
       adctr->Branch("sigmalow"       , &(res.sigmalow       ),"sigmalow/D"       );
       adctr->Branch("sigmahigh"      , &(res.sigmahigh      ),"sigmahigh/D"      );
       adctr->Branch("sigma"          , &(res.sigma          ),"sigma/D"          );
       adctr->Branch("peak1pex"       , &(res.peak1pex       ),"peak1pex/D"       );
       adctr->Branch("peak1pey"       , &(res.peak1pey       ),"peak1pey/D"       );
       adctr->Branch("peak1pesigma"   , &(res.peak1pesigma   ),"peak1pesigma/D"   );
       adctr->Branch("bssigma"        , &(res.bssigma        ),"bssigma/D"        );
       adctr->Branch("peak1pexerr"    , &(res.peak1pexerr    ),"peak1pexerr/D"    );
       adctr->Branch("peak1peyerr"    , &(res.peak1peyerr    ),"peak1peyerr/D"    );
       adctr->Branch("peak1pesigmaerr", &(res.peak1pesigmaerr),"peak1pesigmaerr/D");
       adctr->Branch("valley"         , &(res.valley         ),"valley/D"         );
       adctr->Branch("valleyx"        , &(res.valleyx        ),"valleyx/D"        );
       adctr->Branch("peak"           , &(res.peak           ),"peak/D"           );
       adctr->Branch("pv"             , &(res.pv             ),"pv/D"             );
       adctr->Branch("area1pe"        , &(res.area1pe        ),"area1pe/D"        );
       adctr->Branch("area1peall"     , &(res.area1peall     ),"area1peall/D"     );
       adctr->Branch("areabs"         , &(res.areabs         ),"areabs/D"         );
       adctr->Branch("areapredict"    , &(res.areapredict    ),"areapredict/D"    );
       adctr->Branch("areasigmalow"   , &(res.areasigmalow   ),"areasigmalow/D"   );
       adctr->Branch("areasigmahigh"  , &(res.areasigmahigh  ),"areasigmahigh/D"  );
       adctr->Branch("areasigma"      , &(res.areasigma      ),"areasigma/D"      );
       adctr->Branch("gainpeakx"      , &(res.gainpeakx      ),"gainpeakx/D"      );
       adctr->Branch("gainpeak1pex"   , &(res.gainpeak1pex   ),"gainpeak1pex/D"   );
       adctr->Branch("gainpeak1pexerr", &(res.gainpeak1pexerr),"gainpeak1pexerr/D");
       adctr->Branch("chi2"           , &(res.chi2           ),"chi2/D"           );
       adctr->Branch("ndf"            , &(res.ndf            ),"ndf/D"            );
       adctr->Branch("tmean"          , &(res.tmean          ),"tmean/D"           );
       adctr->Branch("tvariance"      , &(res.tvariance      ),"tvariance/D"       );
       adctr->Branch("tpeakx"         , &(res.tpeakx         ),"tpeakx/D"          );
       adctr->Branch("tpeaky"         , &(res.tpeaky         ),"tpeaky/D"          );
       adctr->Branch("tFWHMlow"       , &(res.tFWHMlow       ),"tFWHMlow/D"        );
       adctr->Branch("tFWHMhigh"      , &(res.tFWHMhigh      ),"tFWHMhigh/D"       );
       adctr->Branch("tEdgelow"       , &(res.tEdgelow       ),"tEdgelow/D"        );
       adctr->Branch("tEdgehigh"      , &(res.tEdgehigh      ),"tEdgehigh/D"       );
       adctr->Branch("tFWHM"          , &(res.tFWHM          ),"tFWHM/D"           );
       adctr->Branch("tgamma"         , &(res.tgamma         ),"tgamma/D"          );
       adctr->Branch("tmu"            , &(res.tmu            ),"tmu/D"             );
       adctr->Branch("tsigma"         , &(res.tsigma         ),"tsigma/D"          );
       adctr->Branch("tlambda"        , &(res.tlambda        ),"tlambda/D"         );
       adctr->Branch("tgammaerr"      , &(res.tgammaerr      ),"tgammaerr/D"       );
       adctr->Branch("tmuerr"         , &(res.tmuerr         ),"tmuerr/D"          );
       adctr->Branch("tsigmaerr"      , &(res.tsigmaerr      ),"tsigmaerr/D"       );
       adctr->Branch("tlambdaerr"     , &(res.tlambdaerr     ),"tlambdaerr/D"      );
       adctr->Branch("tchi2"          , &(res.tchi2          ),"tchi2/D"           );
       adctr->Branch("tndf"           , &(res.tndf           ),"tndf/D"            );
       adctr->Branch("mppchn"         , &(res.mppchn         ),"mppchn/L"          );
       adctr->Branch("pmthn"          , &(res.pmthn          ),"pmthn/L"           );
       adctr->Branch("mppcnpe"        , &(res.mppcnpe        ),"mppcnpe/D"         );
       adctr->Branch("pmtnpegain"     , &(res.pmtnpegain     ),"pmtnpegain/D"      );
       adctr->Branch("pmtspegain"     , &(res.pmtspegain     ),"pmtspegain/D"      );
       adctr->Branch("noiserate"      , &(res.noiserate      ),"noiserate[2]/D"       );
    }
    adctr->Fill();
    adctr->Write(adctrname,TObject::kOverwrite);
    rootf->Close();

    return 0;
}
