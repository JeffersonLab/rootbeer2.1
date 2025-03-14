//
//                  *** TSCtimeCorr.h ***
//--Description
//  Class for correcting SC time information.
//
//  SC counter# is addressed as scid=sector*100+paddle 
//  
// NOTE: if the code does not compile, check and correct the member names of the 
//       structures: gpid_t, part_t, evnt_t, tagr_t, tgpb_t ...
//         (depending on your local clasbanks.ddl and bankheader.h)
//
// Initialization:
// (a)  constructor: TSCtimeCorr(runno,SCIDcutChoice)  (see options below this description)
//           (if SCIDcutChoice not specified, defaults are used for runno from g9a and g9b;
//            to load your own SCIDcut arrays, use runno=0)
//     e.g. TSCtimeCorr *mySCtimeCorr = new TSCtimeCorr(runno,0,"G9b_Franz");
//
// (b)  initialize arrays for time correction (or counter turn off) and correction functions:
//     e.g. mySCtimeCorr->SCIDcut_init(runno)
//
// (c)  print information: which counters are turned off or time corrected
//     e.g. mySCtimeCorr->printInfo();
//
// (d)  load your own SCIDcuts (arrays with scid's) and (optional) your TimeOffsets
//     e.g. int mylist[]={ 101, 201, 301, 401, 501, 601, -1};  // (last entry must be -1 !!)
//          mySCtimeCorr->SCIDcut_arrayInit(mylist);
//
// (e)  correction functions (depending on energy deposit in SCpaddle) can be loaded
//
//
// Methods:
// (a)  check whether a SCpaddle is turned off: 
//     e.g. if( mySCtimeCorr->SCpaddleOff(scid) ) continue;
//
// (b)  get corrected SCtime for a track (row in GPID, PART, EVNT, or provide info)
//        if runno is specified, corrections will be initialized whenever runno changes
//     e.g. float timeSC = mySCtimeCorr->TimeCorr(GPID[row], runno); 
//
// (c)  get corrected Beta for a track (row in GPID, PART, EVNT, or provide Time info):
//        (TAGtime (Tphoton) and vertex pos. (Zvertex) can be specified, e.g. different tagged
//         photon than used in particle bank or MVRT[0].z instead of Zvertex from particle bank) 
//     e.g. float betaSC = mySCtimeCorr->BetaCorr(GPID[row]); 
//
// (d)  get time offset based on correction function (depends on dedx=energy deposit in paddle)
//     e.g. float deltaT = mySCtimeCorr->dedt_Proton(scid, dedx);
//
//
//  Note: this class has to be integrated into ROOT class container list when using under CINT:
//        (i)  generate dictionary SCtimeCorrDict.o 
//              $ROOTSYS/bin/rootcint -f SCtimeCorrDict.cxx -c TSCtimeCorr.h
//        (ii) compile the dictionary file
//              g++ $FLAGS -o SCtimeCorrDict.o -I. -I$ROOTSYS/include -c SCtimeCorrDict.cxx
//              (where $FLAGS= -DR__THREAD -fno-exceptions -fno-strict-aliasing -fPIC  or similar) 
// 
//--Author      F.Klein     Apr 2014
//--Update
//
//---------------------------------------------------------------------------
#ifndef __TSCtimeCorr__
#define __TSCtimeCorr__

#include <TObject.h>
#include "bankvars.h"

const int SC_MAX_COUNTERS = 342;       //all paddles:  id=(sector#-1)*57+paddle#
const int SC_COUNTER_OFF = -499;       //switch paddle off (in SCid_Tcorr[id])

// last entry in list must be SCcut_Nall
enum SCIDcut_List { SCIDcut_none, SCIDcutG9a_Steffen, SCIDcutG9a_Liam, 
		    SCIDcutG9b_Natalie, SCIDcutG9b_Aneta, SCIDcutG9b_Franz, SCIDcutG14b_Haiyun, SCIDcut_Nall};
const char *SCIDcut_Name[SCIDcut_Nall]={"none","G9a_Steffen","G9a_Liam",
					"G9b_Natalie","G9b_Aneta","G9b_Franz","G14b_Haiyun"};

//----------- class definition ----------------------------------------------------

class TSCtimeCorr {

 private:
 protected:
  static const float Clight = 29.97925;         // cm/ns
  int SCIDcutChoice;
  int last_runno;
  int SC_dEdTcorr[SC_MAX_COUNTERS];  // dE dependent time correction for SCpaddles (flag=1: proton, =2: pion; =4: use STdedx)
  float SCid_Tcorr[SC_MAX_COUNTERS];  // time correction for SCpaddles (switch off: SCid_Tcorr[isc] < -100)

  float dedtProt_yoff[SC_MAX_COUNTERS];
  float dedtProt_scal[SC_MAX_COUNTERS];
  float dedtProt_xoff[SC_MAX_COUNTERS];
  float dedtProt_xpow[SC_MAX_COUNTERS];
  float dedtProt_xmin[SC_MAX_COUNTERS];
  float dedtProt_xmax[SC_MAX_COUNTERS];

  float dedtPion_yoff[SC_MAX_COUNTERS];
  float dedtPion_scal[SC_MAX_COUNTERS];
  float dedtPion_xoff[SC_MAX_COUNTERS];
  float dedtPion_xpow[SC_MAX_COUNTERS];
  float dedtPion_xmin[SC_MAX_COUNTERS];
  float dedtPion_xmax[SC_MAX_COUNTERS];

 public:
  TSCtimeCorr(int runno=-1, int mySCcutChoice=-1, char const *NameSCcutChoice=NULL);
  virtual ~TSCtimeCorr();

  virtual void SCIDcut_init(int runno=0);
  virtual void SCIDcut_arrayInit(int *SCIDarray, float *TimeOffsets=NULL);
  virtual void printInfo(int runno=0);
  virtual int SCpaddleId(int scid);
  virtual int SCpaddleId(int sector, int paddle);

  virtual bool SCpaddleOff(int scid);
  virtual bool SCpaddleOff(int sector, int paddle);
  virtual bool SCpaddleOff(GPID_t gpid);
  virtual bool SCpaddleOff(PART_t part);
  virtual bool SCpaddleOff(EVNT_t evnt);

  virtual float TimeCorr(GPID_t gpid, int runno=0);
  virtual float TimeCorr(PART_t part, int runno=0);
  virtual float TimeCorr(EVNT_t evnt, int runno=0);
  virtual float TimeCorr(float SCtime, int scid=0, int runno=0, float mom=0.0, 
			 float dedx=0.0, float STdedx=0.0);

  virtual float BetaCorr(GPID_t gpid, int runno=0, float Tphoton=-500.0, float Zvertex=-500.0);
  virtual float BetaCorr(PART_t part, int runno=0, float Tphoton=-500.0, float Zvertex=-500.0);
  virtual float BetaCorr(EVNT_t evnt, int runno=0, float Tphoton=-500.0, float Zvertex=-500.0);
  virtual float BetaCorr(float SCtime, float SCpathlen, float TAGtime=0.0, float Zvertex=0.0, 
			 int scid=0, int runno=0, float mom=0.0, float dedx=0.0, float STdedx=0.0);

  virtual float dedt_Proton(int scid, float dedx, float Toff=0.0, float yscal=0.0, float xoff=0.0, 
			    float xpow=0.0, float xmin=2.0, float xmax=70.);
  virtual float dedt_Pion(int scid, float dedx, float Toff=0.0, float yscal=0.0, float xoff=0.0, 
			  float xpow=0.0, float xmin=2.0, float xmax=70.);

#ifdef __CINT__
  ClassDef(TSCtimeCorr,1);  //for integration into ROOT: has to be processed using rootcint
#endif
};

//------------------------------------------------------------------------
//---------- constructor and methods -------------------------------------

#ifdef __CINT__
ClassImp(TSCtimeCorr)  //for integration into ROOT: has to be processed using rootcint
#endif

TSCtimeCorr::TSCtimeCorr(int runno, int mySCcutChoice, char const * NameSCcutChoice) {

  last_runno = -2;
  if(mySCcutChoice >= SCIDcut_Nall) 
    SCIDcutChoice = -1; 
  else
    SCIDcutChoice = mySCcutChoice; 
  if(NameSCcutChoice) {
    if(SCIDcutChoice<=0) {
      if( !(strcmp(NameSCcutChoice,"g9a")) || !(strcmp(NameSCcutChoice,"G9a")) )
	SCIDcutChoice = SCIDcutG9a_Steffen;
      else if( !(strcmp(NameSCcutChoice,"g9b")) || !(strcmp(NameSCcutChoice,"G9b")) ) 
	SCIDcutChoice = SCIDcutG9b_Franz;
      else if( !(strcmp(NameSCcutChoice,"g14b")) || !(strcmp(NameSCcutChoice,"G14b")) ) 
	SCIDcutChoice = SCIDcutG14b_Haiyun;
      else {
	for(int i=0; i<SCIDcut_Nall; i++) {
	  if( !(strcmp(NameSCcutChoice,SCIDcut_Name[i])) ) {
	    SCIDcutChoice = i;
	    break;
	  }
	}
      }
    }
  }
  if(SCIDcutChoice < 0) {
    if((runno > 55479) && (runno < 56235))
      SCIDcutChoice = SCIDcutG9a_Steffen;
    else if((runno > 62200) && (runno < 63600))
      SCIDcutChoice = SCIDcutG9b_Franz;
    else if((runno > 67800) && (runno < 69650))
      SCIDcutChoice = SCIDcutG14b_Haiyun;
    else
      SCIDcutChoice = 0;
  }
  if(SCIDcutChoice) fprintf(stdout,"Using list %d (%s) for removing SCpaddles from analysis\n",SCIDcutChoice,SCIDcut_Name[SCIDcutChoice]);

  SCIDcut_init(runno);
}

//---------------------------------------------------------------------------

TSCtimeCorr::~TSCtimeCorr() {
}

//------------------ SCpaddleId --------------------------------------------

int TSCtimeCorr::SCpaddleId(int scid) {
  int isec = int(scid/100) -1;
  if((isec<0) || (isec>5)) return -1;
  int ipd  = (scid%100) -1;
  if((ipd<0) || (ipd>56)) return -1;
  return isec*57+ipd;
}

int TSCtimeCorr::SCpaddleId(int sector, int paddle) {
  if((sector<1) || (sector>6))  return -1;
  if((paddle<1) || (paddle>57)) return -1;
  return (sector-1)*57+paddle-1;
}

//------------------ SCpaddleOff --------------------------------------------

bool TSCtimeCorr::SCpaddleOff(int scid) {

  int isc = SCpaddleId(scid);
  if(isc<0) return true;
  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return true;
  return false;
}

bool TSCtimeCorr::SCpaddleOff(int sector, int paddle) {

  int isc = SCpaddleId(sector,paddle);
  if(isc<0) return true;
  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return true;
  return false;
}

//------------------ SCpaddleOff wrappers ---------------------------------

bool TSCtimeCorr::SCpaddleOff(GPID_t gpid) {

  int isc = SCpaddleId(gpid.sec*100+gpid.paddle);
  if(isc<0) return true;
  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return true;
  return false;
}

bool TSCtimeCorr::SCpaddleOff(PART_t part) {

  int row = part.trkid -1;
  if((row < 0) || (row >= TBID_NH)) return true; 
  int SCpaddle = -1;
  int sc_row = TBID[row].sc_id -1;
  if((sc_row>=0) && (SCRC_NS>0)) {
    for(int irec=0; irec<SCRC_NS; irec++) {
      if((SCRC_S[irec]==TBID[row].sec) && (sc_row<SCRC_NH[irec])) {
	SCpaddle = SCRC[irec][sc_row].id;
	break;
      }
    }
  }
  int isc = SCpaddleId(TBID[row].sec*100+SCpaddle);
  if(isc<0) return true;
  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return true;
  return false;
}

bool TSCtimeCorr::SCpaddleOff(EVNT_t evnt) {

  int sc_row = evnt.SCstat -1;
  if((sc_row < 0) || (sc_row >= SCPB_NH)) return true;
  int scid = int(SCPB[sc_row].ScPdHt/100);
  int isc = SCpaddleId(scid);
  if(isc<0) return true;
  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return true;
  return false;
}

//---------------------------------------------------------------------------
//----------- BetaCorr and TimeCorr with particle bank(row) as input -------------
// time and beta correction with particle bank as input
// Input parameters: 
//        particle bank row (structure): GPID[row], PART[row], EVNT[row]
//        runno (optional => if runno changed, correction functions are initialized)
//        Tphoton (optional - RF corrected TAGtime from PID scheme is used if not specified)
//        Zvertex (optional - Z position from particle bank is used if not specified)
//
// Output: corrected beta (as measured from SC pathlength, SC time, TAG time, particle Z-vertex)
// 
// e.g. float timeSC=TSCtimeCorr::TimeCorr(GPID[row]); 
// e.g. float betaSC=TSCtimeCorr::BetaCorr(GPID[row]); 
//
// NOTE: if it does not compile please check and correct the member names of these structures 
//         (they depend on your local clasbanks.ddl and bankheader.h)

//-------------- BetaCorr(GPID[row]) --------------------------------------------------

float TSCtimeCorr::BetaCorr(GPID_t gpid, int runno, float Tphoton, float Zvertex) {

  if(gpid.q == 0.0) return gpid.betam;  //?? or gpid.beta;

  if(Tphoton < -499) Tphoton = gpid.tpho;
  if(Zvertex < -499) Zvertex = gpid.z;
  return (gpid.sc_len/(TimeCorr(gpid,runno) - Tphoton - Zvertex/Clight)/Clight); 
}

//-------------- BetaCorr(PART[row]) ---------------------------------------------------

float TSCtimeCorr::BetaCorr(PART_t part, int runno, float Tphoton, float Zvertex) {

  int row = part.trkid -1;
  if((row < 0) || (row>=TBID_NH)) return 0.0; 
  if(part.q == 0.0) return TBID[row].beta;

  if(Tphoton < -499) {   //check for smallest time difference
    float diff_time = 5.;              //max. time difference considered
    for( int i=0; i<TAGR_NH; i++ ){
      if( ((TAGR[i].STAT&7)==7) && (fabs(TAGR[i].TPHO-TBID[row].vtime) < diff_time) ){
	Tphoton = TAGR[i].TPHO;
	diff_time = fabs(Tphoton - TBID[row].vtime);
      }
    }
  }

  //get path length from GPID or TDPL
  if(Zvertex < -499) Zvertex = part.z;
  if(GPID_NH) {  
    return (GPID[row].sc_len/(TimeCorr(part,runno) - Tphoton - Zvertex/Clight)/Clight); 
  }
  int tb_row = TBID[row].track -1;
  if( (tb_row>=0) && (TDPL_NS>0) && ((TBER_NH>0) || (TBTR_NH>0)) ) {
    int itrsec = (TBER_NH>0)? (TBER[tb_row].layinfo2>>8)&0xff : (TBTR[tb_row].itr_sec%100);
    int isec   = (TBER_NH>0)? (TBER[tb_row].layinfo2>>24)     : int(TBTR[tb_row].itr_sec/100);
    if( (itrsec>0) && (isec==TBID[row].sec) ) {
      for(int irec=0; irec<TDPL_NS; irec++) {
	if( (TDPL_S[irec]==isec) && (itrsec<=TDPL_NH[irec]) ) {
	  for(int i=0; i<TDPL_NH[irec]; i++) {
	    if( (TDPL[irec][i].trk_pln>tb_row*100+140) && (TDPL[irec][i].trk_pln<tb_row*100+145) && (TDPL[irec][i].x>-900) )
	      return (TDPL[irec][i].tlen/(TimeCorr(part,runno) - Tphoton - Zvertex/Clight)/Clight); 
	  }
	}
      }
    }
  }
  return TBID[row].beta;
}

//-------------- BetaCorr(EVNT[row])  ---------------------------------------------------

float TSCtimeCorr::BetaCorr(EVNT_t evnt, int runno, float Tphoton, float Zvertex) {

  int sc_row = evnt.SCstat -1;
  if((evnt.Charge==0) || (sc_row<0) || (sc_row>=SCPB_NH)) return evnt.Beta;

  if(Zvertex < -499) Zvertex = evnt.Z;
  if(Tphoton < -499) {
    for( int i=0; i<TGPB_NH; i++ ){
      if( TGPB[i].pointer<0 ){
	Tphoton = TGPB[i].Time;
	break;
      }
    }
  }
  return (SCPB[sc_row].Path/(TimeCorr(evnt,runno) - Tphoton - Zvertex/Clight)/Clight); 
}

//-------------- BetaCorr(float,float,float,float,int,int,float,float,float) -------------------

float TSCtimeCorr::BetaCorr(float SCtime, float SCpathlen, float TAGtime, float Zvertex, 
			    int scid, int runno, float mom, float dedx, float STdedx) {

  if(SCtime <= 0.01 || SCpathlen <= 0.01) return 0.0;
  return (SCpathlen/(TimeCorr(SCtime,scid,runno,mom,dedx,STdedx) - TAGtime - Zvertex/Clight)/Clight); 
}

//-------------- TimeCorr(GPID[row]) -------------------------------------------------------

float TSCtimeCorr::TimeCorr(GPID_t gpid, int runno) {

  if(gpid.q == 0.0) return gpid.sc_time;
  float STedep = 0.0;
  float mom = sqrt(pow(gpid.px,2)+pow(gpid.py,2)+pow(gpid.pz,2)); 
  int row = gpid.trkid -1;
  if((row >= 0) && (row < TBID_NH)) {
    int st_row = TBID[row].st_id -1;
    if( STRE_NS>0 ) {
      if( st_row>=0 && st_row < STRE_NH[0] ) 
	STedep = STRE[0][st_row].st_edep;
    }
  }
  int scid = gpid.sec*100 + gpid.paddle;
  return TimeCorr(gpid.sc_time, scid, runno, mom, gpid.dedx, STedep); 
}  

//-------------- TimeCorr(PART[row]) ---------------------------------------------------------

float TSCtimeCorr::TimeCorr(PART_t part, int runno) {

  int row = part.trkid -1;
  if((row < 0) || (row >= TBID_NH)) return -1.0;
  if(part.q == 0.0) return TBID[row].sc_time;
  int scid = 0;
  float STedep = 0.0;
  float SCedep = 0.0;
  float mom = sqrt(pow(part.px,2)+pow(part.py,2)+pow(part.pz,2)); 
  int st_row = TBID[row].st_id -1;
  if( STRE_NS>0 ) {
    if((st_row >= 0) && (st_row < STRE_NH[0])) 
      STedep = STRE[0][st_row].st_edep;
  }
  int sc_row = TBID[row].sc_id -1;
  if((sc_row >= 0) && (SCRC_NS > 0)) {
    for(int irec=0; irec<SCRC_NS; irec++) {
      if( (SCRC_S[irec]==TBID[row].sec) && (sc_row<SCRC_NH[irec]) ) {
	SCedep = SCRC[irec][sc_row].energy;
	scid = TBID[row].sec*100 + SCRC[irec][sc_row].id;
	break;
      }
    }
  }
  return TimeCorr(TBID[row].sc_time, scid, runno, mom, SCedep, STedep); 
}

//-------------- TimeCorr(EVNT[row]) -------------------------------------------------------

float TSCtimeCorr::TimeCorr(EVNT_t evnt, int runno) {

  if(evnt.Charge == 0.0) return 0.0;
  int scid = 0;
  float STedep = 0.0;
  float SCedep = 0.0;
  float SCtime =-1.0;
  int st_row = evnt.STstat -1;
  int sc_row = evnt.SCstat -1;
  if( (st_row>=0) && (st_row<STPB_NH) && (STRE_NS>0)) {
    int iptr=int(STPB[st_row].SThid/10000) -1;
    if((iptr>=0) && (iptr<STRE_NH[0])) 
      STedep = STRE[0][iptr].st_edep;
  }
  if((sc_row >= 0) && (sc_row < SCPB_NH)) {
    scid = int(SCPB[sc_row].ScPdHt/100);
    SCedep = SCPB[sc_row].Edep;
    SCtime = SCPB[sc_row].Time;
  }
  return TimeCorr(SCtime, scid, runno, evnt.Pmom, SCedep, STedep); 
}

//-------------- TimeCorr method -----------------------------------------------------
// Input:  SCtime: measured time from SC 
//         scid:   SC counter# (sector*100+paddle)
//         runno   (opt. => if runno changed, correction functions are initialized)
//         mom,SCdedx,STdedx (opt.: are used for 'timewalk' correction function)
//                           (output is SCtime+TimeOffset (TimeOffset from SCid_Tcorr list if not specified)
// Output: corrected SCtime
//
float TSCtimeCorr::TimeCorr(float SCtime, int scid, int runno, float mom, float SCdedx, float STdedx) {

  int isc = SCpaddleId(scid);
  if((isc<0) || (SCtime<=0.0)) return SCtime;

  if(runno > last_runno) SCIDcut_init(runno);

  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return -1.0;
  float Toffset = SCid_Tcorr[isc];
  if((mom==0.0) || (SCdedx==0.0)) return (SCtime + Toffset);
  if(!SC_dEdTcorr[isc]) return (SCtime + Toffset);

  if( (SC_dEdTcorr[isc]&4) ) {
    float STmxPiEloss=50.*pow((mom+0.2),2.2)+75.;
    if(STdedx > STmxPiEloss) 
      return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
    else
      return (SCtime + dedt_Pion(scid,SCdedx) + Toffset);
  }
  if(SCdedx < 4.0) return (SCtime + Toffset);
  float mxPiEloss=2.2/pow((mom+0.03),1.5) + 10.*pow((mom+1.),0.25);
  if(SCdedx < mxPiEloss) {
    return (SCtime + dedt_Pion(scid,SCdedx) + Toffset);
  }
  else {
    if(mom < 0.18) return (SCtime + Toffset);
    //skip differentiation of entries along dedx line for protons, deuterons 
    // if(mom < 0.41) { 
    // 	float maxEloss=380.*(mom-0.16);
    // 	float minEloss=250.*(mom-0.32);
    // 	if(SCdedx > minEloss && SCdedx < maxEloss)  
    // 	  return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
    // }
    // else {
    // 	float maxEloss=10./pow((mom-0.1),2) + 9.*pow((mom+1.),0.25);
    // 	float minEloss=2./pow((mom-0.06),2.3) + 6.*pow((mom+1.),0.25);
    // 	if(dedx > minEloss && dedx < maxEloss) 
    // 	  return (dedt_Proton(scid,dedx) + Toffset);
    // }
    return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
  }
}

//-------------- dedt_Proton (time correction depending on energy deposit) ----------------------------
//
// dE dependent correction:   dedt_Proton=-yoff-yscal*pow(xde-xoff),xpow)
// sc_id=100*sector+paddle;  xde=energy deposit in SC paddle
// lookup arrays filled in case that parameters (yoff,yscal,xoff,xpow,xmin,xmax) are specified
//
float TSCtimeCorr::dedt_Proton(int scid, float xde, float yoff, float yscal, float xoff, float xpow, float xmin, float xmax) {

  int isc=SCpaddleId(scid);
  if(isc<0) return 0.0;
  if(yscal!=0.0) {
    dedtProt_yoff[isc] = -yoff;
    dedtProt_scal[isc] = -yscal;
    dedtProt_xoff[isc] = xoff;
    dedtProt_xpow[isc] = xpow;
    dedtProt_xmin[isc] = xmin;
    dedtProt_xmax[isc] = xmax;
    SC_dEdTcorr[isc] |=1;
  }

  if(!(SC_dEdTcorr[isc]&1) || (xde <= 0.0)) return 0.0;
  if(dedtProt_xpow[isc] == 0.0) {
    if(xde < dedtProt_xmin[isc]) {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * (dedtProt_xmin[isc] - dedtProt_xoff[isc]);
    }
    else if(xde > dedtProt_xmax[isc]) {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * (dedtProt_xmax[isc] - dedtProt_xoff[isc]);
    }
    else {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * (xde - dedtProt_xoff[isc]);
    }
  }
  else {
    if(xde < dedtProt_xmin[isc]) {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * pow((dedtProt_xmin[isc] - dedtProt_xoff[isc]),dedtProt_xpow[isc]);
    }
    else if(xde > dedtProt_xmax[isc]) {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * pow((dedtProt_xmax[isc] - dedtProt_xoff[isc]),dedtProt_xpow[isc]);
    }
    else {
      return dedtProt_yoff[isc] + dedtProt_scal[isc] * pow((xde - dedtProt_xoff[isc]),dedtProt_xpow[isc]);
    }
  }
}

//-------------- dedt_Pion (time correction depending on energy deposit) ----------------------------
//
// dE dependent correction:   dedt_Pion=-yoff-yscal*pow(xde-xoff),xpow)
// sc_id=100*sector+paddle;  xde=energy deposit in SC paddle
// lookup arrays filled in case that parameters (yoff,yscal,xoff,xpow,xmin,xmax) are specified
//
float TSCtimeCorr::dedt_Pion(int scid, float xde, float yoff, float yscal, float xoff, float xpow, float xmin, float xmax) {

  int isc=SCpaddleId(scid);
  if(isc<0) return 0.0;
  if(yscal!=0.0) {
    dedtPion_yoff[isc] = -yoff;
    dedtPion_scal[isc] = -yscal;
    dedtPion_xoff[isc] = xoff;
    dedtPion_xpow[isc] = xpow;
    dedtPion_xmin[isc] = xmin;
    dedtPion_xmax[isc] = xmax;
    SC_dEdTcorr[isc] |=2;
  }

  if(!(SC_dEdTcorr[isc]&2) || (xde <= 0.0)) return 0.0;
  if(dedtPion_xpow[isc] == 0.0) {
    if(xde > dedtPion_xmax[isc]) {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * (dedtPion_xmin[isc] - dedtPion_xoff[isc]);
    }
    else if(xde > dedtPion_xmax[isc]) {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * (dedtPion_xmax[isc] - dedtPion_xoff[isc]);
    }
    else {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * (xde - dedtPion_xoff[isc]);
    }
  }
  else {
    if(xde < dedtPion_xmin[isc]) {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * pow((dedtPion_xmin[isc] - dedtPion_xoff[isc]),dedtPion_xpow[isc]);
    }
    else if(xde > dedtPion_xmax[isc]) {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * pow((dedtPion_xmax[isc] - dedtPion_xoff[isc]),dedtPion_xpow[isc]);
    }
    else {
      return dedtPion_yoff[isc] + dedtPion_scal[isc] * pow((xde - dedtPion_xoff[isc]),dedtPion_xpow[isc]);
    }
  }
}

//--------------- SCIDcut_init  ---------------------------------------------------------
//
// initialization of corrections for given run number (depending on SCIDcutChoice)
//
void TSCtimeCorr::SCIDcut_init(int runno) {
  
  if(runno <= last_runno) return;
  if((SCIDcutChoice > 0) || (last_runno < 0)) {
    memset(SCid_Tcorr, 0, sizeof(SCid_Tcorr));
    memset(SC_dEdTcorr,  0, sizeof(SC_dEdTcorr));
  }
  last_runno = runno;
  switch (SCIDcutChoice) {
  case SCIDcutG14b_Haiyun:
    {
      int badscid[]={116, 118, 123, 125, 127, 222, 230, 231, 305, 323, 402, 405, 428, 502, 505, 517, 
		     601, 607, 613, 614, -1};
      int i=-1;
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
      int scidno[]={124,  134,  144,  146,  147,  148,  149,  150,  151,  152, 
		    226,  229,  233,  236,  239,  244,  245,  246,  247,  248,  249,  250,  251,  252, 
		    331,  332,  334,  335,  344,  424,  429,  446,  447,  448,  449,  450,  451,
		    520,  525,  530,  535,  546,  547,  548,  549,  550,  551,  552,  553, 
		    605,  629,  630,  633,  642,  647,  649,  650,  651,  652, -1};
      float off[]={1.82, 3.56,-1.34,-1.65,-1.23,-2.11,-1.69,-1.86, -2.2,  3.4, 
		   -12.32,3.16,2.14,-6.68, 3.57,-14.23,-19.48,-2.0,-1.85,-1.7,-15.7, -1.4,-11.9, 1.1,
		   -10.1,1.89,-11.93,-12.0,1.27,-20.14,-20.28,-1.6,-1.7,-25.8,-30.77,-2.0, -2.3,
		   2.06,-19.72,-19.58,2.46,-1.16,-1.3,-1.39, 1.59, -9.5,-12.9, 19.0,  14.0,
		   -1.54,-20.38,-9.94,2.32, 2.15,-1.35,-7.9, -1.8, -2.2, -6.0};
      i=-1; 
      while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = off[i];
    }
    break;

  case SCIDcutG9a_Steffen:
    {
      int badscid[]={245, 423, 523, 552, 553, 644, 646, -1};
      int i=-1;
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
      if(runno >= 56164 && runno <= 56233) SCid_Tcorr[SCpaddleId(136)] = -2.0;
      if(runno >= 55521 && runno <= 55676) SCid_Tcorr[SCpaddleId(344)] = SC_COUNTER_OFF;
      if(runno >= 55669 && runno <= 55676) {
	SCid_Tcorr[SCpaddleId(443)] = -4.0;
	SCid_Tcorr[SCpaddleId(449)] = SCid_Tcorr[SCpaddleId(453)] = 2.0;
      }
      if(runno >= 55664 && runno <= 55668) SCid_Tcorr[SCpaddleId(649)] = 2.0;
      if(runno >= 55521 && runno <= 55676) SCid_Tcorr[SCpaddleId(655)] = -2.0;
    }
    break;

  case SCIDcutG9a_Liam:
    {
      int badscid[]={110, 111, 112, 113, 114, 115, 116, 156, 212, 244, 245, 330, 344, 
		     423, 433, 436, 522, 523, 529, 552, 553, 555, 601, 626, 630, 644, 646, 654, 656, -1};
      int i=-1;
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
    }
    break;
  case SCIDcutG9b_Natalie:
    {
      int badscid[]={224, 244, 245, 249, 251, 322, 338, 415, 448, 449,453, 522,550, 551,552, 553, 554, 555, 
		   612, 613, 649, 653, 656, -1};
      int i=-1;
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
    }
    break;
  case SCIDcutG9b_Aneta: 
  case SCIDcutG9b_Franz: 
    if(runno > 62717) { // && runno < 62987?) {
      int scidno[]= { 502, 405, 416, 323, 226, 633, 134, 239,
		      540, 244, 444, 245, 347, 448, 249, 449, 549, 649, 550, 251, 551, 152, 252, 
		      552, 652, 353, 453, 553, 653, 254, 354, 454, 554, 555, 355, 256, 356, 656, 457, 157,-1};
      float scoff[]={-1.1, 0.0, 2.3, 0.7,-12., 0.7, 0.0,  2.,
		     -0.5,-14., 0.0,-20.,-0.2,-24.,-15.,-29.5,-1.2, -5.,-10.,-11.,-11., 4., 2.,
		     19.2,-4.5,-14.5,-10., 14.,-14.,-0.5,-0.5,-0.5,-16.3,-26.5,-11.,10.,-27., 7.5,-1.2,-0.7};
      //554: 628:-16 653?b:0.45
      int i=-1;
      while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = scoff[i];
      if(runno>62827) {
	SCid_Tcorr[SCpaddleId(448)] -=0.5;
	SCid_Tcorr[SCpaddleId(549)] +=1.1;
	if(runno<62987) {
	  SCid_Tcorr[SCpaddleId(416)] -=2.0;
	  if(runno>62880)
	    SCid_Tcorr[SCpaddleId(633)] = -2.0; 
	}
      }
      if(runno>62880) {
	SCid_Tcorr[SCpaddleId(251)] +=0.8; 
	SCid_Tcorr[SCpaddleId(551)] -=0.8; 
      }
      if(runno>62986) {
	SCid_Tcorr[SCpaddleId(502)] +=1.5; //??+2.0
	SCid_Tcorr[SCpaddleId(124)] +=2.0;
	SCid_Tcorr[SCpaddleId(527)] -=2.0;
	SCid_Tcorr[SCpaddleId(633)] = 0.1;
	SCid_Tcorr[SCpaddleId(239)] +=2.0;
	SCid_Tcorr[SCpaddleId(444)] +=2.0;
	SCid_Tcorr[SCpaddleId(448)] -=0.3;
	SCid_Tcorr[SCpaddleId(249)] +=0.6;
	SCid_Tcorr[SCpaddleId(549)] +=1.1;
	SCid_Tcorr[SCpaddleId(256)] -=1.0;
	if(runno<63330) { //??
	  SCid_Tcorr[SCpaddleId(134)] +=1.0;
	  SCid_Tcorr[SCpaddleId(405)] -=0.5;
	}
	if(runno<63444) //??
	  SCid_Tcorr[SCpaddleId(236)] +=3.0;
	else
	  SCid_Tcorr[SCpaddleId(236)] +=2.5;
      }
      if(runno>63500) { //??  id=633??
	SCid_Tcorr[SCpaddleId(527)] = 0.0;
	SCid_Tcorr[SCpaddleId(134)] +=1.0;
	SCid_Tcorr[SCpaddleId(236)] = 0.0;
	SCid_Tcorr[SCpaddleId(239)] -=0.5;
	SCid_Tcorr[SCpaddleId(444)] +=0.5;
      }
      if(runno>63563) {
	SCid_Tcorr[SCpaddleId(617)] +=0.7;
	SCid_Tcorr[SCpaddleId(618)] +=1.0;
      }
    }
    else {  //runno<62718
      int scidno[]={ 301, 601,  502, 505, 416, 323, 523, 226, 326, 229, 429, 530, 236, 233, 237, 
		     633, 134,  540, 144, 244, 444, 245, 347, 448, 449, 249, 549, 649, 550, 551, 
		     552, 553,  554, 555, 251, 152, 652, 453, 353, 653, 254, 354, 155, 455,
		     355, 256, 356,  656, 157, 257, 457, 657, -1};
      float scoff[]={-0.1,-0.3,-1.0,-1.0, 2.5, 0.5,-2.0,-14., 0.2, 0.4,-0.3,-1.0,-0.5,-0.7,-0.2,
		     -0.4, 2.5,-0.8,-0.4,-14.,-0.2,-20.,-0.5,-24.5,-30.,-15.5,0.6,-5.0,-11.5,-12.2,
		     19.5, 14.,-16.,-26.5,-10.2,4.0,-5.,-9.0,-14.5,-14.,-0.5,-0.5,-0.4,-0.2,
		     -11.,  8.,-27.5,7.5,-0.8,-0.5,-0.5, -0.6}; 
      int i=-1; 
      while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = scoff[i];
      if(runno<62289) {
	SCid_Tcorr[SCpaddleId(344)]= 0.5;
	//	SCid_Tcorr[SCpaddleId(554)] =-14.0;  //??
      }
      if(runno>62288) {
	//	SCid_Tcorr[SCpaddleId(554)] =-14.0;  //??
	SCid_Tcorr[SCpaddleId(155)] -= 0.4;
      }
      if(runno>62604) {
	SCid_Tcorr[SCpaddleId(236)]= 1.0;
	SCid_Tcorr[SCpaddleId(140)]= -2.0;
	SCid_Tcorr[SCpaddleId(143)]= -2.0;
	SCid_Tcorr[SCpaddleId(147)] = -0.5;
      }
      if(runno>62400 && runno<62600) {
	if(runno<62500) { //??
	  SCid_Tcorr[SCpaddleId(237)] -= 0.5;
	}
	SCid_Tcorr[SCpaddleId(249)] += 0.7;
	SCid_Tcorr[SCpaddleId(251)] += 2.7;
	SCid_Tcorr[SCpaddleId(353)] += 1.5;
      }
      if(runno>62464) {
	SCid_Tcorr[SCpaddleId(237)]= -0.6;
      }
      if(runno>62564) {
	SCid_Tcorr[SCpaddleId(617)]= 0.7;
	SCid_Tcorr[SCpaddleId(618)]= 1.0;
	SCid_Tcorr[SCpaddleId(621)]= -1.3;
      }
    }
    break;
  }

  if(SCIDcutChoice != SCIDcutG9b_Franz) return;

  //set flag to use STdedx for time correction functions
  int scid_STedep[]={ 416, 123, 323, 326, 429, 629, 134, 234, 236, -1};
  int j=-1;
  while(scid_STedep[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep[j])]=4;
  if(runno<62604)
    SC_dEdTcorr[SCpaddleId(229)]=4;
  if(runno<62718) {
    int scid_STedep1[]={ 530, 331, 233, 444, 347, -1};
    j=-1;
    while(scid_STedep1[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep1[j])]=4;
    if(runno>62603) {
      SC_dEdTcorr[SCpaddleId(142)]=4;
      SC_dEdTcorr[SCpaddleId(144)]=4;
    }
  }
  else { // run>62717
    int scid_STedep2[]={ 334, 335, -1};
    j=-1;
    while(scid_STedep2[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep2[j])]=4;
  }
  //dE dependent corrections
  dedt_Proton(613,0.0, -0.221807,0.0046941,1.51377e-07,1.81129,2,70);
  dedt_Proton(614,0.0, -0.275289,0.0190288,2.4,1.5271,2.5,70); 
  if(runno>62288 && runno<62827) //?? 
    dedt_Proton(416,0.0, -0.822875,-0.0062907,1.6,1.8062,2,70.);
  dedt_Proton(116,0.0, -0.0746125,0.000757609,-5.6,2.09368,2,70);
  dedt_Proton(320,0.0, -0.162416,0.0249439,4,1.07305,4.1,70);
  dedt_Proton(321,0.0, 0.0138948,0.00277524,3.2,1.45354,3.3,70);
  dedt_Proton(222,0.0, -0.0201432,0.00626166,3.2,1.6035,3.3,70);
  dedt_Proton(522,0.0, 0.0100314,0.00951013,3.2,1.53726,3.3,70);
  dedt_Proton(622,0.0, -0.00240647,0.023751,4.8,1.31659,4.9,70);
  dedt_Proton(123,0.0, -0.117054,0.00850427,1.6,1.9268,2,70);
  dedt_Proton(223,0.0, 0.0979165,0.000503644,3.2,1.95479,3.3,70.);
  dedt_Proton(323,0.0, -0.981879,0.0139318,0.8,2.77031,2,70);
  dedt_Proton(423,0.0, -0.254501,0.0172978,4,1.1982,4.1,70.);
  dedt_Proton(623,0.0, 0.0332715,0.0187258,1.46817,1,2,70);
  if(runno<62828)
    dedt_Proton(124,0.0, -0.0298794,0.0141087,6.4,1.31286,6.5,70.);
  dedt_Proton(224,0.0, 0.050504,0.000523787,5.6,2.05784,5.7,70.);
  dedt_Proton(324,0.0, 0.023401,0.0024826,1.6,1.55345,2,70.);
  dedt_Proton(524,0.0, 0.070833,0.00282938,4,1.57203,4.1,70.);
  dedt_Proton(624,0.0, 0.029032,0.00166275,0.8,1.76111,2,70.);
  dedt_Proton(125,0.0, -0.161318,0.0165393,6.4,1.31379,6.5,70.);
  dedt_Proton(425,0.0, -0.376865,1.50562e-06,6.4,3.32393,6.5,70.);
  dedt_Proton(525,0.0, 0.139603,0.000841331,2.4,1.83562,2.5,70.);
  dedt_Proton(226,0.0, -0.0907017,0.00192274,4,1.51669,4.1,70.);
  dedt_Proton(127,0.0, -0.0581778,0.00416721,3.2,1.77332,3.3,70.);
  dedt_Proton(427,0.0, -0.216233,0.000160107,8.8,1.9781,8.9,70.);
  dedt_Proton(228,0.0, -0.351308,0.0352189,-1.21466,1,2,70.);
  dedt_Proton(528,0.0, 0.14166,0.000282926,-8.8,2.00977,2,70.);
  if(runno<62605) 
    dedt_Proton(229,0.0, -0.00994333,0.0201451,-2.27297,1,2,70.);
  else
    dedt_Proton(229,0.0, -0.816903,0.0570301,2.67814,1,2.77814,70.);
  dedt_Proton(529,0.0, -0.00371472,0.00138978,-6.21869,1.68372,2,70.);
  if(runno<62987)
    dedt_Proton(629,0.0, -1.13568,0.209839,-0.472198,1,2,70.);
  dedt_Proton(230,0.0, -0.335781,0.000142584,-14.6622,2.34831,2,70.);
  dedt_Proton(430,0.0, 0.254787,3.69113e-09,13.6,4.77487,13.7,70.);
  if(runno<62718)
    dedt_Proton(530,0.0, -0.168344,0.0432862,-1.21441,1,2,70.);
  dedt_Proton(630,0.0, -0.468713,0.00940293,-17.5995,1.29508,2,70.);
  dedt_Proton(331,0.0, 0.0302998,0.00628346,-4.78109,1.75404,2,70.);
  if(runno<62718)
    dedt_Proton(133,0.0, 0.338925,9.58889e-10,2.4,4.64712,2.5,70.);
  dedt_Proton(433,0.0, -3.0,1.6113,2.4,0.277499,2.5,70.);
  dedt_Proton(533,0.0, -0.0536817,0.000965426,0.8,1.68532,2,70.);
  dedt_Proton(234,0.0, 3.0,-2.473,0.8,0.1,2,70.);
  dedt_Proton(334,0.0, 2.74176,-1.99182,4.3469,0.116272,4.4469,70.);
  dedt_Proton(434,0.0, -1.75476,1.34208,-0.792445,0.101983,2,70.);
  dedt_Proton(335,0.0, -0.593243,0.225001,-2.4,0.391838,2,70.);
  dedt_Proton(337,0.0, -3,1.95507,4,0.204048,4.1,70.);
  dedt_Proton(538,0.0, -0.442077,0.0256601,0.8,1.17523,2,70.);
  dedt_Proton(239,0.0, 3.,-2.17317,6.36194,0.254568,6.46194,70.); //+pi
  dedt_Proton(439,0.0, 2.26499,-1.9095,1.6,0.1,2,70.);
  dedt_Proton(539,0.0, -0.199483,0.0144536,1.6,1.07901,2,70.);
  if(runno<62604 || runno>62717)
    dedt_Proton(142,0.0, -0.165495,0.0145086,-1.64358,0,2,70.);
  dedt_Proton(442,0.0, -0.2203,0.011732,1.59851,0,2,70.);
  dedt_Proton(542,0.0, -0.167378,0.0207299,6.87864,0,2,70.);
  dedt_Proton(144,0.0, 3.0,-1.40134,0.767807,0.264974,2,70.);
  if(runno<62718)
    dedt_Proton(445,0.0, 1.9195,-1.69332,4,0.1,4.1,70.);
  dedt_Proton(545,0.0, 1.4652,-1.21742,1.6,0.1,2,70.);
  dedt_Proton(246,0.0, -0.748499,0.0701757,2.4,1.0863,2.5,70.);
  dedt_Proton(446,0.0, -0.327215,0.00459207,1.6,1.65999,2,70.);
  dedt_Proton(546,0.0, -0.22831,0.00788455,4,1.48365,4.1,70.);
  dedt_Proton(247,0.0, -0.289359,0.0198042,3.2,1.21912,3.3,70.);
  dedt_Proton(547,0.0, -0.227409,0.0042275,2.4,1.62923,2.5,70.);
  dedt_Proton(148,0.0, -1.57372,0.238965,-1.6,0.750701,2,70.);
  dedt_Proton(248,0.0, -0.180264,0.00390619,1.6,1.64098,2,70.);
  dedt_Proton(348,0.0, -0.39176,0.0238184,1.6,1.21057,2,70.);
  dedt_Proton(149,0.0, -0.203037,0.0123528,3.2,1.31378,3.3,70.);
  dedt_Proton(649,0.0, 0.278938,0.0217861,0.8,1.26697,2,70.);
  dedt_Proton(250,0.0, -0.128981,0.00898725,5.6,1.50674,5.7,70.);
  dedt_Proton(350,0.0, -0.844932,0.0900319,0.842105,0,2,70.);
  dedt_Proton(650,0.0, -0.697032,0.0467138,-1.6905,0,2,70.);
  dedt_Proton(151,0.0, -0.296683,0.0178688,2.4,1.33207,2.5,70.);
  dedt_Proton(351,0.0, -0.313371,0.0227898,2.11182,0,2,70.);
  dedt_Proton(351,0.0, -0.149567,0.00331017,3.2,1.479,3.3,70.);
  dedt_Proton(451,0.0, -0.224098,0.00399342,1.6,1.6371,2,70.);
  dedt_Proton(651,0.0, -0.265437,0.0146552,3.2,1.30948,3.3,70.);
  dedt_Proton(153,0.0, -0.0329951,0.00968248,2.4,1.44213,2.5,70.);
  dedt_Proton(253,0.0, -0.276209,0.0176413,2.4,1.45461,2.5,70.);
  
  dedt_Pion(301,0.0, -0.0734455,0.0200517,9.98922,0,5.6,70);
  dedt_Pion(601,0.0, -0.208378,0.00170688,1.43809e-08,1.91284,10.4,70);
  dedt_Pion(502,0.0,  0.0195738,0.178761,6.72543,1.2528,6.8,70);
  dedt_Pion(305,0.0, -0.544275,0.0136702,0.8,1.63424,2,70.);
  dedt_Pion(505,0.0, -0.752772,0.221787,0.00338407,1,2.,70);
  dedt_Pion(604,0.0, -0.6591,0.0369152,4.8,1.62802,4.9,70.);
  if(runno > 62717)
    dedt_Pion(306,0.0, -0.34188,0.0462657,6.4,1.59913,6.5,70.);
  dedt_Pion(607,0.0, -0.50582,0.0025939,-3.99694,1.95829,2,70.);
  dedt_Pion(613,0.0, -0.593048,0.0773409,1.82273,0,4,70);
  dedt_Pion(614,0.0, -0.153899,0.100957,6.70603,0,4.8,70); 
  dedt_Pion(118,0.0, -0.0838246,0.0022169,1.6,1.7118,2,70.); //??
  dedt_Pion(122,0.0, -0.0175184,0.00049934,-7.99972,1.83505,2,70.);
  dedt_Pion(525,0.0, -0.779324,0.5,3.2,0.238471,3.3,70.);
  dedt_Pion(127,0.0, -0.277002,0.0187147,0.8,1.24111,2,70.);
  dedt_Pion(228,0.0, -0.102935,0.0157658,0.79733,0,2,70.);
  dedt_Pion(528,0.0, -0.830921,0.475982,-0.601109,0.232823,2,70.);
  if(runno>62604) 
    dedt_Pion(229,0.0, -0.85554,0.5,3.2,0.277431,3.3,70.);
  dedt_Pion(529,0.0, -0.023221,0.00184413,3.2,1.38569,3.3,70.);
  dedt_Pion(629,0.0, -0.376183,0.0186985,0.784218,0,2,70.);
  dedt_Pion(230,0.0, -0.10738,0.0241076,2.38848,0,2,70.);
  if(runno<62718) 
    dedt_Pion(530,0.0, 3.0,-2.24322,-0.16768,0.157661,2,70.);
  else 
    dedt_Pion(530,0.0, 3.0,-2.13406,0.8,0.133982,2,70.);
  dedt_Pion(630,0.0, -2.77573,2.34239,6.44143,0.114437,6.54143,70.);
  dedt_Pion(331,0.0, 0.0302998,0.00628346,-4.78109,1.75404,2,70.);
  if(runno<62718) 
    dedt_Pion(133,0.0, 0.338925,9.58889e-10,2.4,4.64712,2.5,70.);
  dedt_Pion(433,0.0, -3,1.6113,2.4,0.277499,2.5,70.);
  dedt_Pion(533,0.0, -0.0536817,0.000965426,0.8,1.68532,2,70.);
  dedt_Pion(234,0.0, 3.0,-2.473,0.8,0.1,2,70.);
  dedt_Pion(334,0.0, 2.74176,-1.99182,4.3469,0.116272,4.4469,70.);
  dedt_Pion(434,0.0, -1.75476,1.34208,-0.792445,0.101983,2,70.);
  dedt_Pion(335,0.0, -0.593243,0.225001,-2.4,0.391838,2,70.);
  dedt_Pion(337,0.0, -3,1.95507,4,0.204048,4.1,70.);
  dedt_Pion(538,0.0, -0.442077,0.0256601,0.8,1.17523,2,70.);
  dedt_Pion(439,0.0, 2.26499,-1.9095,1.6,0.1,2,70.);
  dedt_Pion(539,0.0, -0.199483,0.0144536,1.6,1.07901,2,70.);
  if(runno<62604 || runno>62717) 
    dedt_Pion(142,0.0, -0.165495,0.0145086,-1.64358,0,2,70.);
  dedt_Pion(442,0.0, -0.2203,0.011732,1.59851,0,2,70.);
  dedt_Pion(542,0.0, -0.167378,0.0207299,6.87864,0,2,70.);
  dedt_Pion(144,0.0, 3.0,-1.40134,0.767807,0.264974,2,70.);
  if(runno<62718) 
    dedt_Pion(445,0.0, 1.9195,-1.69332,4,0.1,4.1,70.);
  dedt_Pion(545,0.0, 1.4652,-1.21742,1.6,0.1,2,70.);
  dedt_Pion(246,0.0, -0.748499,0.0701757,2.4,1.0863,2.5,70.);
  dedt_Pion(446,0.0, -0.327215,0.00459207,1.6,1.65999,2,70.);
  dedt_Pion(546,0.0, -0.22831,0.00788455,4,1.48365,4.1,70.);
  dedt_Pion(247,0.0, -0.289359,0.0198042,3.2,1.21912,3.3,70.);
  dedt_Pion(547,0.0, -0.227409,0.0042275,2.4,1.62923,2.5,70.);
  dedt_Pion(148,0.0, -1.57372,0.238965,-1.6,0.750701,2,70.);
  dedt_Pion(248,0.0, -0.180264,0.00390619,1.6,1.64098,2,70.);
  dedt_Pion(348,0.0, -0.39176,0.0238184,1.6,1.21057,2,70.);
  dedt_Pion(149,0.0, -0.203037,0.0123528,3.2,1.31378,3.3,70.);
  dedt_Pion(649,0.0, 0.278938,0.0217861,0.8,1.26697,2,70.);
  dedt_Pion(250,0.0, -0.128981,0.00898725,5.6,1.50674,5.7,70.);
  dedt_Pion(350,0.0, -0.844932,0.0900319,0.842105,0,2,70.);
  dedt_Pion(450,0.0, -1.57119,0.0548546,-10.3995,1.04209,2,70.);
  dedt_Pion(650,0.0, -0.697032,0.0467138,-1.6905,0,2,70.);
  dedt_Pion(151,0.0, -0.296683,0.0178688,2.4,1.33207,2.5,70.);
  dedt_Pion(351,0.0, -0.149567,0.00331017,3.2,1.479,3.3,70.);
  dedt_Pion(451,0.0, -0.224098,0.00399342,1.6,1.6371,2,70.);
  if(runno>62717) {
    dedt_Pion(551,0.0, -0.647282,0.0368326,-7.99932,1.17183,2,70.);
    dedt_Pion(152,0.0, -0.481817,0.00244113,4,1.88254,4.1,70.);
    dedt_Pion(555,0.0, -1.72141,0.000759948,-7.19954,2.07456,2,70.);
  }
  dedt_Pion(651,0.0, -0.265437,0.0146552,3.2,1.30948,3.3,70.);
  dedt_Pion(153,0.0, -0.0329951,0.00968248,2.4,1.44213,2.5,70.);
  dedt_Pion(253,0.0, -0.276209,0.0176413,2.4,1.45461,2.5,70.);
}

//--------------- SCIDcut_arrayInit  -------------------------------------------------
//
// initialize SCID_Tcorr[] array:
// Input: array with SCpaddles (notation: sector*100+paddle) - LAST ENTRY MUST BE -1 !!
//        optional: array of same length with time offsets to be loaded.
//          if TimeOffsets[] is not specified, all counters in SCIDarray[] are turned off
//
void TSCtimeCorr::SCIDcut_arrayInit(int *SCIDarray, float *TimeOffsets) {

  int i=-1;
  if(TimeOffsets) {
    while(SCIDarray[++i] > 0) SCid_Tcorr[SCpaddleId(SCIDarray[i])] = SC_COUNTER_OFF;
  }
  else {
    while(SCIDarray[++i] > 0) SCid_Tcorr[SCpaddleId(SCIDarray[i])] = TimeOffsets[i];
  }
}

//-------------- printInfo(runno) -------------------------------------------------------

void TSCtimeCorr::printInfo(int runno) {
  
  if(runno>0)
    std::cout<<"SCtimeCorr:  SC paddles turned off for run "<<runno<<": ";
  else
    std::cout<<"SCtimeCorr:  SC paddles turned off: ";
  int j=0, k=0, l=0;
  for(int i=0; i<SC_MAX_COUNTERS; i++) {
    if(SCid_Tcorr[i] <= SC_COUNTER_OFF) {
      std::cout<<" "<<(int(i/57)+1)*100+(i%57)+1;
      j++;
    }
    else {
      if(SCid_Tcorr[i]!=0.0) k++;
      if((SC_dEdTcorr[i]&3)) l++;
    }
  }
  if(!j) 
    std::cout<<"  none"<<std::endl;
  else
    std::cout<<std::endl;
  if(k>0) {
    if(runno>0)
      std::cout<<"SCtimeCorr:  SC paddles time corrected for run "<<runno<<": ";
    else
      std::cout<<"SCtimeCorr:  SC paddles time corrected: ";
    for(int i=0; i<SC_MAX_COUNTERS; i++) {
      if((SCid_Tcorr[i] > SC_COUNTER_OFF) && (SCid_Tcorr[i] != 0.0)) 
	std::cout<<"  "<<(int(i/57)+1)*100+(i%57)+1<<" ("<<SCid_Tcorr[i]<<")";
    }
    std::cout<<std::endl;
  }
  if(l>0) {
    std::cout<<"  (SC timewalk correction function for proton:";
    for(int i=0; i<SC_MAX_COUNTERS; i++) {
      if( (SC_dEdTcorr[i]&1) ) { 
	std::cout<<"  "<<(int(i/57)+1)*100+(i%57)+1;
      }
    }
    std::cout<<")\n  (SC timewalk correction function for pion: ";
    for(int i=0; i<SC_MAX_COUNTERS; i++) {
      if( (SC_dEdTcorr[i]&2) ) { 
	std::cout<<"  "<<(int(i/57)+1)*100+(i%57)+1;
      }
    }
    std::cout<<")"<<std::endl;
  }
}

//---------------------------------------------------------------------------

#ifdef __CINT__ 
#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 
#pragma link C++ class TSCtimeCorr+; 
#endif 

#endif
