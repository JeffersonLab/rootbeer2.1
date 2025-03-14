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
//  for (a) and (b): target_center z-position can be added as optional last parameter
//                   (automatically set for default SCIDcutChoice)
//     e.g. TSCtimeCorr *mySCtimeCorr = new TSCtimeCorr(runno,0,"G9b_Franz",target_Zoffset);
//        or: mySCtimeCorr->SCIDcut_init(runno,target_Zoffset)  
//
// (c)  print information: which counters are turned off or time corrected
//     e.g. mySCtimeCorr->PrintInfo();  or optional: mySCtimeCorr->PrintInfo(runno,flag)
//       (with flag=1: print list of SC paddles that are turned off or time corrected;
//             flag=2: print SCPID beta cuts used for corrected particle id.
//             flag=3: print both)
//
// (d)  load your own SCIDcuts (arrays with scid's) and (optional) your TimeOffsets
//     e.g. int mylist[]={ 101, 201, 301, 401, 501, 601, -1};  // (last entry must be -1 !!)
//          mySCtimeCorr->SCIDcut_arrayInit(mylist);
//
// (e)  correction functions (depending on energy deposit in SCpaddle) can be loaded
//
// (f)  get particle-ID based on 2-dim cut in beta vs. momentum:
//    - define order of particle type checked in 2-dim cut
//       for e+,e-; pi+,pi-; K+,K-; proton,antiproton; deuteron; Helium3,triton
//     e.g. mySCtimeCorr->SCPID_SetSearchOrder(SCPID_Kaon,SCPID_Proton,SCPID_Pion,SCPID_Electron,SCPID_Deuteron,SCPID_Helium3);
//    - set parameters for the beta cut: 
//        beta_min < betaMeasured < beta_max  for p_min[itype] < momentum < p_max[itype]
//       with beta_min=momentum/sqrt(momentum^2+mass_min[itype]^2)+off_min[itype]  (dito for beta_max)
//     e.g. mySCtimeCorr->SCPID_CutInit(SCPID_Kaon,0,55,0.44,-0.04,0.04,0.25,1.7)
//        (note: add another variable to change the 'PDGmass' for beta_calc=mom/sqrt(mom^2+PDGmass^2)
//
// Methods:
// (a)  check whether a SCpaddle is turned off: 
//     e.g. if( mySCtimeCorr->SCpaddleOff(scid) ) continue;
//
// (b)  get corrected SCtime for a track (row in GPID, PART, EVNT, or provide info)
//        if runno is specified, corrections will be initialized whenever runno changes
//     e.g. double timeSC = mySCtimeCorr->TimeCorr(GPID[row], runno); 
//
// (c)  get corrected Beta for a track (row in GPID, PART, EVNT, or provide Time info):
//        (TAGtime (Tphoton) and vertex pos. (Zvertex) can be specified, e.g. different tagged
//         photon than used in particle bank or MVRT[0].z instead of Zvertex from particle bank) 
//     e.g. double betaSC = mySCtimeCorr->BetaCorr(GPID[row]); 
//
// (d)  get time offset based on correction function (depends on dedx=energy deposit in paddle)
//     e.g. double deltaT = mySCtimeCorr->dedt_Proton(scid, dedx);
//
// (e) get corrected particle id. from beta vs. momentum cuts  
//      optional input: betaMeasured; if betaMeasured=0: take betaMeasure=BetaCorr(...);
//      optional output: probabilities for all tested particle tpyes, based on |betaMeasured-betaCalc|
//     e.g. double myPIDprob[SCPID_ntypes];
//          int myPID = mySCtimeCorr->PIDcorr(GPID[row],myPIDprob) 
//            (myPID is Geant-ID for GPID,PART or (Lund,PDG) MCid for EVNT;
//             probabilities for other particle types in array myPIDprob[]  
//
//  (f) - PID conversion Geant-ID <-> MCid
//            e.g. int myGID = mySCtimeCorr->PIDconvert(myMCID,SCPIDin_MCID) <=from MCid to Geant-ID
//      - getPDGmass to calculate beta=mom/energy
//            e.g  double mass0 = mySCtimeCorr->getPDGmass(myGID,SCPIDin_GID)  <=input= Geant-ID
//
//  Note: this class has to be integrated into ROOT class container list when using under CINT:
//        (i)  generate dictionary SCtimeCorrDict.o 
//              $ROOTSYS/bin/rootcint -f SCtimeCorrDict.cxx -c TSCtimeCorr.h
//        (ii) compile the dictionary file
//              g++ $FLAGS -o SCtimeCorrDict.o -I. -I$ROOTSYS/include -c SCtimeCorrDict.cxx
//              (where $FLAGS= -DR__THREAD -fno-exceptions -fno-strict-aliasing -fPIC  or similar) 
// 
//--Author      F.Klein     Apr 2014
//--Update      F.Klein     Feb 2015
//
//---------------------------------------------------------------------------
#ifndef __TSCtimeCorr__
#define __TSCtimeCorr__

#include <TObject.h>
#include "bankvars.h"

const int SC_MAX_COUNTERS = 342;       //all paddles:  id=(sector#-1)*57+paddle#
const int SC_COUNTER_OFF = -500;       //switch paddle off (in SCid_Tcorr[id])
const int SC_DEFAULT_VALUE = -500;     //take times,positions from reconstructed banks
const int SC_CHECK_TPHO   = 500;       //check all entries in TAGR bank 

// last entry in list must be SCcut_Nall
enum SCIDcut_List { SCIDcut_none, SCIDcutG9a_Steffen, SCIDcutG9a_Liam, 
		    SCIDcutG9b_Natalie, SCIDcutG9b_Aneta, SCIDcutG9b_Franz, 
		    SCIDcutG14b_Haiyun, SCIDcutG14b_Franz, SCIDcut_Nall};
#ifndef __INCLUDE_ONLY__
const char *SCIDcut_Name[SCIDcut_Nall]={"none","G9a_Steffen","G9a_Liam",
					"G9b_Natalie","G9b_Aneta","G9b_Franz",
					"G14b_Haiyun","G14b_Franz"};
#endif

// SCPIDcuts for: [0]=e+-, [1]=pi+-, [2]=K+-, [3]=proton, [4]=deuteron, [5]=helium3 
enum SCPID_list {SCPID_electron, SCPID_pion, SCPID_kaon, SCPID_proton, SCPID_deuteron, SCPID_helium3, SCPID_ntypes };

//PID conversion: MCID (PDG,Lund), GID (geant3), myid (internal: SCPID_list), 
//                postype,negtype: from SCPID_list
enum SCPIDin_Flag {SCPIDin_myid=1, SCPIDin_GID, SCPIDin_MCID, SCPIDin_postype, SCPIDin_negtype};

enum MYPIDlist {mypid_Gamma ,mypid_Positron, mypid_Electron, mypid_MuonPlus, mypid_MuonMinus, mypid_Pion0, mypid_PionPlus, mypid_PionMinus, mypid_Kaon0long, mypid_KaonPlus, mypid_KaonMinus, mypid_Neutron, mypid_Proton, mypid_AntiProton, mypid_Deuteron, mypid_Triton, mypid_Helium3, mypid_Helium4, nMYPIDcodes};

#ifndef DEF_Clight
#define DEF_Clight
static const double Clight = 29.97925;         // (cm/ns)
#endif

//----------- class definition ----------------------------------------------------

class TSCtimeCorr {

 private:
  //arrays for particle Id. cuts for different particle types
  //[0]=e+-, [1]=pi+-, [2]=K+-, [3]=proton, [4]=deuteron, [5]=helium3            
  //default search order: SCPID_kaon,SCPID_proton,SCPID_pion,SCPID_electron,SCPID_deuteron,SCPID_helium3
  int PIDsearchorder[SCPID_ntypes];
  double PIDminProb;  //min.probability to assign particle id. in PIDcorr routine
  double PIDmass2[SCPID_ntypes];  //mass^2 (GeV^2)
  double beta_mass2_max[SCPID_ntypes]; //(GeV^2)
  double beta_mass2_min[SCPID_ntypes];
  double beta_off_max[SCPID_ntypes];
  double beta_off_min[SCPID_ntypes];
  double beta_p_max[SCPID_ntypes];
  double beta_p_min[SCPID_ntypes];
  double target_center;    //(cm) target offset (provided by user or default for run periods)

  int SCIDcutChoice;
  int last_runno;
  int SC_dEdTcorr[SC_MAX_COUNTERS];  //dE dependent time correction for SCpaddles 
  //(bit 0: for proton, 1: for pion; 2: use STdedx; 4: bad counter for PIDcorr; 
  //(bit 12-15: scale SCdedx by 2.8,2.4,1.5,1.3 for SC cut ranges in PIDcorr)
  //(bit 16-19: scale SCdedx by 0.9,0.8,0.7,0.6 for SC cut ranges in PIDcorr)
  //(bit 20-23: scale STdedx by 0.7,0.85,1.2,1.4 for SC cut ranges in PIDcorr)
  double SCid_Tcorr[SC_MAX_COUNTERS];  //(ns) time correction for SCpaddles (switch off: SCid_Tcorr[isc] < -100)

  double dedtProt_yoff[SC_MAX_COUNTERS];
  double dedtProt_scal[SC_MAX_COUNTERS];
  double dedtProt_xoff[SC_MAX_COUNTERS];
  double dedtProt_xpow[SC_MAX_COUNTERS];
  double dedtProt_xmin[SC_MAX_COUNTERS];
  double dedtProt_xmax[SC_MAX_COUNTERS];

  double dedtPion_yoff[SC_MAX_COUNTERS];
  double dedtPion_scal[SC_MAX_COUNTERS];
  double dedtPion_xoff[SC_MAX_COUNTERS];
  double dedtPion_xpow[SC_MAX_COUNTERS];
  double dedtPion_xmin[SC_MAX_COUNTERS];
  double dedtPion_xmax[SC_MAX_COUNTERS];
  double getSTedep(int st_row, int id_trksec=0);
  double scaleSCdedx(int scid);
  double scaleSTdedx(int scid);

 public:
  TSCtimeCorr(int runno=-1, int mySCcutChoice=-1, char *NameSCcutChoice=NULL, double tg_center=-500.0);
  virtual ~TSCtimeCorr();

  virtual void SCIDcut_init(int runno=0, double tg_center=-500.0);
  virtual void SCIDcut_arrayInit(int *SCIDarray, double *TimeOffsets=NULL);
  virtual void SCIDcut_arrayInit(int *SCIDarray, float *TimeOffsets);
  virtual void PrintInfo(int runno=0, int flag=1);
  virtual int SCpaddleId(int scid);
  virtual int SCpaddleId(int sector, int paddle);

  virtual bool SCpaddleOff(int scid);
  virtual bool SCpaddleOff(int sector, int paddle);
  virtual bool SCpaddleOff(GPID_t gpid);
  virtual bool SCpaddleOff(PART_t part);
  virtual bool SCpaddleOff(EVNT_t evnt);

  virtual double TimeCorr(GPID_t gpid, int runno=0);
  virtual double TimeCorr(PART_t part, int runno=0);
  virtual double TimeCorr(EVNT_t evnt, int runno=0);
  virtual double TimeCorr(double SCtime, int scid=0, int runno=0, double mom=0.0, 
			 double dedx=0.0, double STdedx=0.0);

  virtual double BetaCorr(GPID_t gpid, int runno=0, double Tphoton=-500.0, double Zvertex=-500.0);
  virtual double BetaCorr(GPID_t gpid, int *itagr, int runno=0, double Zvertex=-500.0);
  virtual double BetaCorr(PART_t part, int runno=0, double Tphoton=-500.0, double Zvertex=-500.0);
  virtual double BetaCorr(PART_t part, int *itagr, int runno=0, double Zvertex=-500.0);
  virtual double BetaCorr(EVNT_t evnt, int runno=0, double Tphoton=-500.0, double Zvertex=-500.0);
  virtual double BetaCorr(EVNT_t evnt, int *itagr, int runno=0, double Zvertex=-500.0);
  virtual double BetaCorr(double SCtime, double SCpathlen, double TAGtime=0.0, double Zvertex=0.0, 
			 int scid=0, int runno=0, double mom=0.0, double dedx=0.0, double STdedx=0.0);

  virtual int PIDcorr(GPID_t gpid, double *PIDprob=NULL, double mom_corr=0.0, double betaMeasured=0.0, double SCdedx=0.0, double STdedx=0.0);  //output: Gid   
  virtual int PIDcorr(PART_t part, double *PIDprob=NULL, double mom_corr=0.0, double betaMeasured=0.0, double SCdedx=0.0, double STdedx=0.0);  //output: Gid   
  virtual int PIDcorr(EVNT_t evnt, double *PIDprob=NULL, double mom_corr=0.0, double betaMeasured=0.0, double SCdedx=0.0, double STdedx=0.0);  //output: MCid 

  // Note: q_mom =charge*momentum (where charge=+1 or -1)
  virtual int PIDcorr(double q_mom, double betaMeasured, double *PIDprob=NULL, double SCdedx=0.0, double STdedx=0.0, int scid=0); //output: Gid 

  virtual void SCPID_SetSearchOrder(int j0=-1,int j1=-1,int j2=-1,int j3=-1,int j4=-1,int j5=-1);
  virtual void SCPID_CutInit(int itype=-1, double mass_min=-10., double mass_max=-10., double off_min=-10., double off_max=-10., double p_min=-10., double p_max=-10., double PIDmass=-10.);
  virtual void SCPID_SetMinProb(double minProb=0.85) { PIDminProb = minProb; };
  virtual double SCPID_GetMinProb(void) { return PIDminProb; };

  //PID conversion: Gid <-> MCid  (flag: see SCPID_convertFlag)
  virtual int PIDconvert(int id, int flag=SCPIDin_GID);
  virtual double getPDGmass(int id, int flag=SCPIDin_GID);
  virtual double getBetaCalc(double momentun, int id, int flag=SCPIDin_GID);

  virtual double dedt_Proton(int scid, double dedx, double Toff=0.0, double yscal=0.0, double xoff=0.0, 
			    double xpow=0.0, double xmin=2.0, double xmax=70.);
  virtual double dedt_Pion(int scid, double dedx, double Toff=0.0, double yscal=0.0, double xoff=0.0, 
			  double xpow=0.0, double xmin=2.0, double xmax=70.);

#ifdef __CINT__
  ClassDef(TSCtimeCorr,1);  //for integration into ROOT: has to be processed using rootcint
#endif
};

//------------------------------------------------------------------------
//---------- constructor and methods -------------------------------------
#ifndef __INCLUDE_ONLY__

#ifdef __CINT__
ClassImp(TSCtimeCorr)  //for integration into ROOT: has to be processed using rootcint
#endif

TSCtimeCorr::TSCtimeCorr(int runno, int mySCcutChoice, char* NameSCcutChoice, double tg_center) {

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
  if((SCIDcutChoice < 0) || (!SCIDcutChoice && !NameSCcutChoice)) {
    if((runno > 55479) && (runno < 56235))
      SCIDcutChoice = SCIDcutG9a_Steffen;
    else if((runno > 62200) && (runno < 63600))
      SCIDcutChoice = SCIDcutG9b_Franz;
    else if((runno > 67800) && (runno < 69650))
      SCIDcutChoice = SCIDcutG14b_Haiyun;
    else
      SCIDcutChoice = 0;
  }
  if(tg_center > -499.) {
    target_center=tg_center;
  }
  else {
    if((SCIDcutChoice==SCIDcutG14b_Haiyun)||(SCIDcutChoice==SCIDcutG14b_Franz)||((runno > 67800) && (runno < 69650)))
      target_center=-7.5;
  }
  if(SCIDcutChoice) fprintf(stdout,"Using list %d (%s) for removing SCpaddles from analysis\n",SCIDcutChoice,SCIDcut_Name[SCIDcutChoice]);

  SCIDcut_init(runno);

  //fill beta vs. momentum cut ranges with default values
  //[0]=e+-, [1]=pi+-, [2]=K+-, [3]=proton, [4]=deuteron, [5]=helium3            
  SCPID_SetMinProb();
  SCPID_SetSearchOrder( 2,3,1,0,4,5 );
  for(int i=0; i<SCPID_ntypes; i++) SCPID_CutInit(i);
}

//---------------------------------------------------------------------------

TSCtimeCorr::~TSCtimeCorr() {
}

//---------------------------------------------------------------------------
//-------------  private:  getSTedep(st_row, id_trksec)    ------------------
//  st_row, id_trksec starting from 1  (id_trksec = sector*100 + track# in this sector

double TSCtimeCorr::getSTedep(int st_row, int id_trksec) {
  if( !STRE_NS ) return 0.0;
  if( st_row>0 && st_row <= STRE_NH[STRE_NS-1] ) 
    return STRE[STRE_NS-1][st_row-1].st_edep;
  else if(id_trksec > 100) {
    int is100 = int(id_trksec/100)*100;
    for(int i=0; i<STRE_NH[STRE_NS-1]; i++) {
      int id = STRE[STRE_NS-1][i].ID;
      if( id > is100 && id < is100+99 ) {
	if(STRE[STRE_NS-1][i].Trk_no == (id_trksec%100)) 
	  return STRE[STRE_NS-1][i].st_edep;
      }
    }
  }
  return 0.0;
}

//-------------  private:  scaleSTdedx(scid)    ------------------
// scid<=0 means SCpaddleId already calculated: SCpaddleId=-scid 
double TSCtimeCorr::scaleSTdedx(int scid) {
  int isc = (scid>0)? SCpaddleId(scid) : -scid;
  if(isc < 0) return 1.0;
  double factor = 1.0;
  if(SC_dEdTcorr[isc]&0x0F00000) {
    if(SC_dEdTcorr[isc]&0x10000) factor *=0.7;
    if(SC_dEdTcorr[isc]&0x20000) factor *=0.85;
    if(SC_dEdTcorr[isc]&0x40000) factor *=1.2;
    if(SC_dEdTcorr[isc]&0x80000) factor *=1.4;
  }
  return factor;
}
//-------------  private:  scaleSCdedx(scid)    ------------------
// scid<=0 means SCpaddleId already calculated: SCpaddleId=-scid 
double TSCtimeCorr::scaleSCdedx(int scid) {
  int isc = (scid>0)? SCpaddleId(scid) : -scid;
  if(isc < 0) return 1.0;
  double factor = 1.0;
  if(SC_dEdTcorr[isc]&0x0FF000) {
    if(SC_dEdTcorr[isc]&0x01000) factor *=2.8;
    if(SC_dEdTcorr[isc]&0x02000) factor *=2.4;
    if(SC_dEdTcorr[isc]&0x04000) factor *=1.5;
    if(SC_dEdTcorr[isc]&0x08000) factor *=1.3;
    if(SC_dEdTcorr[isc]&0x10000) factor *=0.9;
    if(SC_dEdTcorr[isc]&0x20000) factor *=0.8;
    if(SC_dEdTcorr[isc]&0x40000) factor *=0.7;
    if(SC_dEdTcorr[isc]&0x80000) factor *=0.6;
  }
  return factor;
}

//---------------------------------------------------------------------------
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
//             (note: Zvertex is corrected for target_offset in the routine)
//
// Output: corrected beta (as measured from SC pathlength, SC time, TAG time, particle Z-vertex)
// 
// e.g. double timeSC=TSCtimeCorr::TimeCorr(GPID[row]); 
// e.g. double betaSC=TSCtimeCorr::BetaCorr(GPID[row]); 
//
// NOTE: if it does not compile please check and correct the member names of these structures 
//         (they depend on your local clasbanks.ddl and bankheader.h)

//-------------- BetaCorr(GPID[row]) --------------------------------------------------

double TSCtimeCorr::BetaCorr(GPID_t gpid, int runno, double Tphoton, double Zvertex) {

  if(gpid.q == 0.0) return gpid.betam; 
  if(Tphoton < -499 || Tphoton > 499) {
    int itagr = (Tphoton>499)? TAGR_NH : gpid.tagrid -1;
    return BetaCorr(gpid, &itagr, runno, Zvertex);
  }
  if(Zvertex < -499) Zvertex = gpid.z;
  return (gpid.sc_len/(TimeCorr(gpid,runno) - Tphoton - (Zvertex-target_center)/Clight)/Clight); 
}

//  itagr: input =0...TAGR_NH-1: take time for this TAGR entry;
//         input <0   take gpid value;    input >=TAGR_NH: find best TAGR entry
//         output: row number (0,..,TAGR_NH-1) of TAGR entry used to calculate Beta
double TSCtimeCorr::BetaCorr(GPID_t gpid, int *itagr, int runno, double Zvertex) {

  if(gpid.q == 0.0) {
    *itagr = gpid.tagrid -1;
    return gpid.betam; 
  }
  if(Zvertex < -499) Zvertex = gpid.z;
  double dtime = TimeCorr(gpid,runno) - (Zvertex-target_center)/Clight;
  double Tphoton = -500.;
  
  if(*itagr<0) {
    *itagr = gpid.tagrid -1;
    if(*itagr<0 && gpid.st_time < dtime && gpid.st_len > 1.) {
      return (gpid.sc_len - gpid.st_len)/(TimeCorr(gpid) - gpid.st_time)/Clight;
    }
  }
  if(*itagr >= 0 && *itagr < TAGR_NH) 
    Tphoton = TAGR[*itagr].TPHO; 
  else {
    *itagr = -1;
    double diff_time = 9.;              //max. time difference considered
    double mom2 = gpid.px*gpid.px+gpid.py*gpid.py+gpid.pz*gpid.pz;
    for( int j=0; j<SCPID_ntypes; j++) {                //check all particle types
      double betaX = sqrt(mom2)/sqrt(mom2+PIDmass2[j]);
      for( int i=0; i<TAGR_NH; i++ ){
	//	if( (TAGR[i].stat&7)==7 && fabs(dtime-gpid.sc_len/betaX/Clight-TAGR[i].TPHO) < diff_time){
	if(TAGR[i].STAT==15 && fabs(dtime-gpid.sc_len/betaX/Clight-TAGR[i].TPHO) < diff_time){
	  Tphoton = TAGR[i].TPHO;
	  *itagr = i;
	  diff_time = fabs(dtime - gpid.sc_len/betaX/Clight -Tphoton);
	}
      }
    }
  }
  double mybeta = gpid.sc_len/(dtime - Tphoton)/Clight; 
  if(mybeta < 0.05 || mybeta > 1.3) mybeta = -1.0;
  return mybeta;
}

//-------------- BetaCorr(PART[row]) ---------------------------------------------------

double TSCtimeCorr::BetaCorr(PART_t part, int runno, double Tphoton, double Zvertex) {

  int row = part.trkid -1;
  if((row < 0) || (row>=TBID_NH)) return -1.0; 
  if(part.q == 0.0) return TBID[row].beta;
  if(Tphoton < -499 || Tphoton > 499) {
    int itagr = TAGR_NH;
    return BetaCorr(part, &itagr, runno, Zvertex);
  }
  // fake entry in TAGR bank and correct afterwards (need path length to TOF!)
  int itagr = TAGR_NH-1;
  float lastTpho = TAGR[itagr].TPHO;
  TAGR[itagr].TPHO = Tphoton;
  double mybeta = BetaCorr(part, &itagr, runno, Zvertex);
  TAGR[itagr].TPHO = lastTpho;
  return mybeta;
}

double TSCtimeCorr::BetaCorr(PART_t part, int *itagr, int runno, double Zvertex) {

  int row = part.trkid -1;
  if((row < 0) || (row>=TBID_NH)) return -1.0;
  if(part.q == 0.0)  return TBID[row].beta;

  //get path length from GPID or TDPL
  if(Zvertex < -499) Zvertex = part.z;
  double dtime = TimeCorr(part,runno) - (Zvertex-target_center)/Clight;
  double pathlen = 0.0;
  if(GPID_NH) {  
    pathlen = GPID[row].sc_len;
  }
  else {
    int tb_row = TBID[row].track -1;
    if( (tb_row>=0) && (TDPL_NS>0) && ((TBER_NH>0) || (TBTR_NH>0)) ) {
      int itrsec = (TBER_NH>0)? (TBER[tb_row].layinfo2>>8)&0xff : (TBTR[tb_row].itr_sec%100);
      int isec   = (TBER_NH>0)? (TBER[tb_row].layinfo2>>24)     : int(TBTR[tb_row].itr_sec/100);
      if( (itrsec>0) && (isec==TBID[row].sec) ) {
	for(int irec=0; irec<TDPL_NS; irec++) {
	  if( (TDPL_S[irec]==isec) && (itrsec<=TDPL_NH[irec]/10) ) {
	    for(int i=0; i<TDPL_NH[irec]; i++) {
	      if( (TDPL[irec][i].trk_pln>tb_row*100+140) && (TDPL[irec][i].trk_pln<tb_row*100+145) && (TDPL[irec][i].x>-900) )
		pathlen = TDPL[irec][i].tlen;
	    }
	  }
	}
      }
    }
  }
  double Tphoton = -500.;
  if(*itagr>=0 && *itagr<TAGR_NH) {
    Tphoton = TAGR[*itagr].TPHO; 
  }
  else {              //check for smallest time difference
    *itagr = -1;
    double diff_time = 9.;              //max. time difference considered
    double mom2 = part.px*part.px+part.py*part.py+part.pz*part.pz;
    for( int j=0; j<SCPID_ntypes; j++) {
      double betaX = sqrt(mom2)/sqrt(mom2+PIDmass2[j]);
      for( int i=0; i<TAGR_NH; i++ ){
	if(TAGR[i].STAT==15 && fabs(dtime-pathlen/betaX/Clight-TAGR[i].TPHO) < diff_time){
	  Tphoton = TAGR[i].TPHO;
	  *itagr = i;
	  diff_time = fabs(dtime - pathlen/betaX/Clight - Tphoton);
	}
      }
    }
  }
  double mybeta = pathlen/(dtime - Tphoton)/Clight; 
  if(mybeta < 0.05 || mybeta > 1.3) mybeta = -1.0;
  return mybeta;
}

//-------------- BetaCorr(EVNT[row])  ---------------------------------------------------

double TSCtimeCorr::BetaCorr(EVNT_t evnt, int runno, double Tphoton, double Zvertex) {

  if((Tphoton < -499) || (Tphoton > 499)) {
    int itgpb = (Tphoton>499)? TGPB_NH : -1;
    return BetaCorr(evnt, &itgpb, runno, Zvertex);
  }
  int sc_row = evnt.SCstat -1;
  if((evnt.Charge==0) || (sc_row<0) || (sc_row>=SCPB_NH)) return evnt.Beta;

  if(Zvertex < -499) Zvertex = evnt.Z;
  return (SCPB[sc_row].Path/(TimeCorr(evnt,runno) - Tphoton - (Zvertex-target_center)/Clight)/Clight); 
}

double TSCtimeCorr::BetaCorr(EVNT_t evnt, int *itgpb, int runno, double Zvertex) {

  int sc_row = evnt.SCstat -1;
  if((evnt.Charge==0) || (sc_row<0) || (sc_row>=SCPB_NH)) return evnt.Beta;

  if(Zvertex < -499) Zvertex = evnt.Z;
  double dtime = TimeCorr(evnt,runno) - (Zvertex-target_center)/Clight;

  double Tphoton = -500.;
  if((*itgpb>=0) && (*itgpb<TGPB_NH)) {
    Tphoton = TGPB[*itgpb].Time;
  }
  else if(*itgpb < 0) {
    *itgpb = TGPB_NH;
    for( int i=0; i<TGPB_NH; i++ ){
      if( TGPB[i].pointer<0 ){
	Tphoton = TGPB[i].Time;
	*itgpb = i;
	break;
      }
    }
  }
  if(*itgpb >= TGPB_NH) {
    *itgpb = -1;
    double diff_time = 9.;              //max. time difference considered
    for( int j=0; j<SCPID_ntypes; j++) {
      double betaX = evnt.Pmom/sqrt(evnt.Pmom*evnt.Pmom+PIDmass2[j]);
      for( int i=0; i<TGPB_NH; i++ ){
	int tagstat = (TMath::Abs(TGPB[i].pointer))%100;
	if( (tagstat&7)==7 && fabs(dtime-SCPB[sc_row].Path/betaX/Clight-TGPB[i].Time) < diff_time){
	  Tphoton = TGPB[i].Time;
	  *itgpb = i;
	  diff_time = fabs(dtime - SCPB[sc_row].Path/betaX/Clight - Tphoton);
	}
      }
    }
  }
  double mybeta = SCPB[sc_row].Path/(dtime - Tphoton)/Clight; 
  if(mybeta < 0.05 || mybeta > 1.3) mybeta = -1.0;
  return mybeta;
}

//-------------- BetaCorr(double,double,double,double,int,int,double,double,double) -------------------

double TSCtimeCorr::BetaCorr(double SCtime, double SCpathlen, double TAGtime, double Zvertex, 
			    int scid, int runno, double mom, double dedx, double STdedx) {

  if(SCtime <= 0.01 || SCpathlen <= 0.01) return 0.0;
  return (SCpathlen/(TimeCorr(SCtime,scid,runno,mom,dedx,STdedx) - TAGtime - (Zvertex-target_center)/Clight)/Clight); 
}

//---------------------------------------------------------------------------
//-------------- TimeCorr(GPID[row]) -------------------------------------------------------

double TSCtimeCorr::TimeCorr(GPID_t gpid, int runno) {

  if(gpid.q == 0.0) return gpid.sc_time;
  double STedep = 0.0;
  double mom = sqrt(pow(gpid.px,2)+pow(gpid.py,2)+pow(gpid.pz,2)); 
  int row = gpid.trkid -1;
  if((row >= 0) && (row < TBID_NH)) {
    STedep = getSTedep(TBID[row].st_id);
  }
  int scid = gpid.sec*100 + gpid.paddle;
  return TimeCorr(gpid.sc_time, scid, runno, mom, gpid.dedx, STedep); 
}  

//-------------- TimeCorr(PART[row]) ---------------------------------------------------------

double TSCtimeCorr::TimeCorr(PART_t part, int runno) {

  int row = part.trkid -1;
  if((row < 0) || (row >= TBID_NH)) return -1.0;
  if(part.q == 0.0) return TBID[row].sc_time;
  int scid = 0;
  double STedep = getSTedep(TBID[row].st_id);
  double SCedep = 0.0;
  double mom = sqrt(pow(part.px,2)+pow(part.py,2)+pow(part.pz,2)); 
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

double TSCtimeCorr::TimeCorr(EVNT_t evnt, int runno) {

  if(evnt.Charge == 0.0) return 0.0;
  int scid = 0;
  double STedep = 0.0;
  double SCedep = 0.0;
  double SCtime =-1.0;
  int st_row = evnt.STstat -1;
  int sc_row = evnt.SCstat -1;
  if( (st_row>=0) && (st_row<STPB_NH)) {
    int is100 = (int(STPB[st_row].SThid/100)%100)*100;
    STedep = getSTedep( (STPB[st_row].SThid/10000), is100+STPB[st_row].charge);
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
double TSCtimeCorr::TimeCorr(double SCtime, int scid, int runno, double mom, double SCdedx, double STdedx) {

  int isc = SCpaddleId(scid);
  if((isc<0) || (SCtime<=0.0)) return SCtime;

  if(runno > last_runno) SCIDcut_init(runno);

  if(SCid_Tcorr[isc] <= SC_COUNTER_OFF) return -1.0;
  double Toffset = SCid_Tcorr[isc];
  if((mom==0.0) || (SCdedx==0.0)) return (SCtime + Toffset);
  if(!(SC_dEdTcorr[isc]&0x0FFF)) return (SCtime + Toffset);

  if( (SC_dEdTcorr[isc]&4) ) {
    STdedx *= scaleSTdedx(-isc);
    double STmxPiEloss=50./pow((mom+0.2),2.2)+75.;
    if(STdedx > STmxPiEloss) 
      return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
    else
      return (SCtime + dedt_Pion(scid,SCdedx) + Toffset);
  }
  double dEscale = scaleSCdedx(-isc);
  if(dEscale >= 1.0) {
    if(SCdedx < 4.0) return (SCtime + Toffset);
  }
  else {
    if(SCdedx*dEscale < 4.0) return (SCtime + Toffset);
  }
  double mxPiEloss=2.2/pow((mom+0.03),1.5) + 10.*pow((mom+1.),0.25);
  if(SCdedx*dEscale < mxPiEloss) {
    return (SCtime + dedt_Pion(scid,SCdedx) + Toffset);
  }
  else {
    if(mom < 0.18) return (SCtime + Toffset);
    //skip differentiation of entries along dedx line for protons, deuterons 
    // if(mom < 0.41) { 
    // 	double maxEloss=380.*(mom-0.16);
    // 	double minEloss=250.*(mom-0.32);
    // 	if(SCdedx > minEloss && SCdedx < maxEloss)  
    // 	  return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
    // }
    // else {
    // 	double maxEloss=10./pow((mom-0.1),2) + 9.*pow((mom+1.),0.25);
    // 	double minEloss=2./pow((mom-0.06),2.3) + 6.*pow((mom+1.),0.25);
    // 	if(dedx > minEloss && dedx < maxEloss) 
    // 	  return (dedt_Proton(scid,dedx) + Toffset);
    // }
    return (SCtime + dedt_Proton(scid,SCdedx) + Toffset);
  }
}

//-----------------------------------------------------------------------------------
//-------------- SCPID methods ------------------------------------------------------
//   calculate particle id. based on beta vs. momentum cuts:
//       beta_min < betaMeasured < beta_max  and  p_min[itype] < mom < p_max[itype]
//   with beta_min=mom/sqrt(mom^2+mass_min[itype]^2) + offset_min[itype]   (similar for beta_max);
//   particle types are checked in specified order (set via: SCPID_SetSearchOrder(....))
//
//------------  PIDcorr (output: GEANT id)   ----------------------------------------
//   q_mom = charge*momentum  (charge=+1 or -1)
//   optional output array: PIDprob[itype]: probability for this particle being of type 'itype'
//       (probability based on |betaMeas - betaCalc| and energy loss)

int TSCtimeCorr::PIDcorr(double q_mom, double betaMeas, double *PIDprob, double SCdedx, double STdedx, int scid) {

  bool calcProb = (PIDprob)? true : false;
  if(calcProb) for(int i=0; i<SCPID_ntypes; i++) PIDprob[i]=0.0;
  if(betaMeas<=0.0) return 0;
  double mom = (q_mom>0)? q_mom : -q_mom;
  double mom2 = mom*mom;

  //bad counter? reduce all probabilities by 0.2 ... or if SCpaddle off by 0.8
  double badctrOff = 0.0;
  double SCdecut = 5.0;
  int isc = SCpaddleId(scid);
  if(isc >=0) {
    if(SCid_Tcorr[isc]<=SC_COUNTER_OFF) badctrOff = 0.8;
    else if(SC_dEdTcorr[isc]&16)        badctrOff = 0.2;
    if(SC_dEdTcorr[isc]&4) SCdecut = 20.0;
    SCdedx *= scaleSCdedx(-isc);
    STdedx *= scaleSTdedx(-isc);
  }
  //guess using dedx for SC and ST
  int id_dedx = SCPID_ntypes;
  if(mom > 0.08 && SCdedx > 0.0) {
    double STmxpieloss = 50./pow(mom+0.2,2.2) +75.;
    id_dedx = -1;
    if(SCdedx < SCdecut) {
      if(STdedx > STmxpieloss-10.)
	id_dedx = SCPID_proton;
      else if(STdedx > 2.) 
	id_dedx = SCPID_pion;
    }
    else {
      double mxpieloss=2.2/pow((mom+0.03),1.5) + 10.*pow((mom+1.),0.25);
      if(SCdedx < mxpieloss && STdedx < STmxpieloss+10.) {
	if(STdedx==0.0 || STdedx > 8.) id_dedx = SCPID_pion;
      }
      else if(STdedx==0.0 || STdedx > STmxpieloss*0.8) {
	if(mom > 0.18 && mom < 0.41) {
	  double maxeloss=380.*(mom-0.16);
	  double mineloss=250.*(mom-0.32);
	  if(SCdedx > mineloss && SCdedx < maxeloss)
	    id_dedx = SCPID_proton;
	}
	else if(mom > 0.41) {
	  double maxeloss=10./pow((mom-0.1),2) + 9.*pow((mom+1.),0.25);
	  double mineloss=2./pow((mom-0.06),2.3) + 6.*pow((mom+1.),0.25);
	  if(SCdedx > maxeloss) 
	    id_dedx = SCPID_deuteron;
	  else if(SCdedx > mineloss) 
	    id_dedx = SCPID_proton;
	}
      }
    }
  }
  double dedxProb[SCPID_ntypes];
  memset(dedxProb , 0, sizeof(dedxProb));
  for(int i=0; i<SCPID_ntypes; i++) {
    int itype=PIDsearchorder[i];
    if(itype<0) continue;
    double betaCalc = mom/sqrt(mom2+PIDmass2[itype]);
    double betaMin = mom/sqrt(mom2+beta_mass2_min[itype]) + beta_off_min[itype];
    double betaMax = mom/sqrt(mom2+beta_mass2_max[itype]) + beta_off_max[itype];
    int fac_dedx = (id_dedx==itype)? 1.0 : 3.0; //factor for difference to calculated value;
    if(fac_dedx > 1.0) {
      if(itype==SCPID_kaon && mom>1.1 && (id_dedx==SCPID_pion || id_dedx==SCPID_proton)) fac_dedx=2.;
      if(itype==SCPID_helium3 && id_dedx==SCPID_deuteron) fac_dedx=2.;
      if(id_dedx==SCPID_ntypes) fac_dedx=1.5; //small reduction if SCdedx==0
    }
    if(betaMeas > betaMin && betaMeas < betaMax &&
       mom > beta_p_min[itype] && mom < beta_p_max[itype]) {
      dedxProb[itype] = 1.0 - pow(fac_dedx*(betaMeas - betaCalc),2) - badctrOff;
      if(calcProb && dedxProb[itype]>0.0) PIDprob[itype] = dedxProb[itype];
    }
    else if(calcProb) {
      double diffsum = (fac_dedx-1.0) * 0.05;
      if(betaMeas > betaMax) {
	diffsum += pow(5*(betaMax - betaCalc),2) + 3*(betaMeas - betaCalc);
      }
      else if(betaMeas < betaMin) {
	diffsum += pow(5*(betaCalc - betaMin),2) + 3*(betaMin - betaMeas);
      }
      if(mom > beta_p_max[itype]) {
	diffsum += (mom - beta_p_max[itype])/beta_p_max[itype];
      }
      else if(mom < beta_p_min[itype]) {
	diffsum += beta_p_min[itype] - mom;
      }
      PIDprob[itype] = 1.0 - diffsum - badctrOff;
      if(PIDprob[itype]<0) PIDprob[itype] = 0.0;
    }
  }
  int ifound = TMath::LocMax(SCPID_ntypes, dedxProb);
  if(dedxProb[ifound] < PIDminProb) return 0;
  //particle id
  if(q_mom > 0) 
    return PIDconvert(ifound, SCPIDin_postype+10*SCPIDin_GID);
  else
    return PIDconvert(ifound, SCPIDin_negtype+10*SCPIDin_GID);
}

//------------  PIDcorr wrapper for GPID (output: GEANT id)   ---------------------------------

int TSCtimeCorr::PIDcorr(GPID_t gpid, double *PIDprob, double mom_corr, double betaMeas, double SCdedx, double STdedx) {
  if(gpid.q == 0.0) {
    if(PIDprob) {
      for(int i=0; i<SCPID_ntypes; i++) PIDprob[i]=0.0;
    }
    return gpid.pid;
  }
  if(betaMeas<=0.0) betaMeas = BetaCorr(gpid);
  if(SCdedx==0.0)   SCdedx = gpid.dedx;
  if(STdedx==0.0) {
    int row = gpid.trkid -1;
    if((row >= 0) && (row < TBID_NH)) 
      STdedx = getSTedep(TBID[row].st_id);
  }
  double mom = (mom_corr>0.0)? mom_corr : sqrt(gpid.px*gpid.px+gpid.py*gpid.py+gpid.pz*gpid.pz);
  if(gpid.q<0) mom *= -1.;
  return PIDcorr(mom, betaMeas, PIDprob, SCdedx, STdedx, (gpid.sec*100+gpid.paddle));
}

//------------  PIDcorr wrapper for PART (output: GEANT id)   ---------------------------------

int TSCtimeCorr::PIDcorr(PART_t part, double *PIDprob, double mom_corr, double betaMeas, double SCdedx, double STdedx) {
  if(part.q == 0.0) {
    if(PIDprob) {
      for(int i=0; i<SCPID_ntypes; i++) PIDprob[i]=0.0;
    }
    return part.pid;
  }
  if(betaMeas<=0.0) betaMeas = BetaCorr(part);
  int scid = 0;
  int row = part.trkid -1;
  if((row >= 0) && (row < TBID_NH)) {
    if(STdedx==0.0)  STdedx = getSTedep(TBID[row].st_id);
    int sc_row = TBID[row].sc_id -1;
    if((sc_row>=0) && (SCRC_NS>0)) {
      for(int irec=0; irec<SCRC_NS; irec++) {
	if((SCRC_S[irec]==TBID[row].sec) && (sc_row<SCRC_NH[irec])) {
	  if(SCdedx==0.0)  SCdedx = SCRC[irec][sc_row].energy;
	  scid = TBID[row].sec*100 + SCRC[irec][sc_row].id;
	  break;
	}
      }
    }
  }
  double mom = (mom_corr>0.0)? mom_corr : sqrt(part.px*part.px+part.py*part.py+part.pz*part.pz);
  if(part.q<0) mom *= -1.;
  return PIDcorr(mom, betaMeas, PIDprob, SCdedx, STdedx, scid);
}

//------------  PIDcorr wrapper for EVNT (output: MCid (Lund,PDG))   ----------------------------

int TSCtimeCorr::PIDcorr(EVNT_t evnt, double *PIDprob, double mom_corr, double betaMeas, double SCdedx, double STdedx) {
  if(evnt.Charge == 0.0) {
    if(PIDprob) {
      for(int i=0; i<SCPID_ntypes; i++) PIDprob[i]=0.0;
    }
    return evnt.ID;
  }
  int scid = 0;
  if(betaMeas<=0.0) betaMeas = BetaCorr(evnt);
  if(evnt.SCstat>0 && evnt.SCstat<=SCPB_NH) { 
    if(SCdedx==0.0)  SCdedx = SCPB[evnt.SCstat-1].Edep;
    scid = int(SCPB[evnt.SCstat-1].ScPdHt/100);
  }
  if(STdedx==0.0 && evnt.STstat>0 && evnt.STstat<=STPB_NH) {
    int is100 = (int(STPB[evnt.STstat-1].SThid/100)%100)*100;
    STdedx = getSTedep( (STPB[evnt.STstat-1].SThid/10000), is100+STPB[evnt.STstat-1].charge);
  }
  int mycharge = (evnt.Charge<0)? -1 : 1;
  double mom = (mom_corr>0.0)? mycharge*mom_corr : mycharge*evnt.Pmom;
  int myid = PIDcorr(mom, betaMeas, PIDprob, SCdedx, STdedx, scid);
  return PIDconvert(myid);
}

//------------  SCPID_SetSearchOrder (sequence in which particle types are tested)  --------------
//       for SCPID_electron(e+,e-); SCPID_pion(pi+,pi-); SCPID_kaon (K+,K-); 
//           SCPID_proton (proton,antiproton); SCPID_deuteron; SCPID_helium3

void TSCtimeCorr::SCPID_SetSearchOrder(int j0, int j1, int j2, int j3, int j4, int j5) {

  for(int i=0; i<SCPID_ntypes; i++) PIDsearchorder[i]=-1;
  const char *clst[SCPID_ntypes]={"electron","kaon","pion","proton","deuteron","helium3"};
  int mylst[SCPID_ntypes];
  mylst[0]=j0; mylst[1]=j1; mylst[2]=j2; mylst[3]=j3; mylst[4]=j4; mylst[5]=j5;
  int k=0;
  for(int i=0; i<SCPID_ntypes; i++) {
    if(mylst[i]>=0 && mylst[i]<SCPID_ntypes) { 
      if(PIDsearchorder[mylst[i]]<0) 
	PIDsearchorder[mylst[i]]=i;
      else {
	std::cout<<"TSCtimeCorr::SCPID_SearchOrder: double entry for SCPID_type \""<<clst[mylst[i]]<<"\" discarded"<<std::endl;
	k++;
      }
    }
    else {
      mylst[i]=-1;
      k++;
    }
  }
  if(k>0) {
    std::cout<<"TSCtimeCorr::SCPID_SetSearchOrder: ";
    for(int i=0; i<SCPID_ntypes; i++) {
      if(mylst[i]>=0) std::cout<<clst[mylst[i]]<<", ";
    }
    std::cout<<std::endl;
  }
}

//------------  SCPID_CutInit for particle type 'itype'   ---------------------------------
//  set parameters for the beta vs momentum cut: 
//     beta_min < betaMeasured < beta_max  for p_min(itype) < mom < p_max(itype)
//     with beta_min=mom/sqrt(mom^2+mass_min(itype)^2)+off_min(itype)  (dito for beta_max)

void TSCtimeCorr::SCPID_CutInit(int itype, double mass_min, double mass_max, double off_min, double off_max, double p_min, double p_max, double mass) {
  
  if(itype<0 || itype>=SCPID_ntypes) {
    std::cout<<"TSCtimeCorr::SCPID_CutInit type out of range (0,"<<SCPID_ntypes<<"): "<<itype<<std::endl;
    return;
  }
  //[0]=e+-, [1]=pi+-, [2]=K+-, [3]=proton, [4]=deuteron, [5]=helium3            
  int myid[SCPID_ntypes]={mypid_Electron,mypid_PionPlus,mypid_KaonPlus,mypid_Proton,mypid_Deuteron,mypid_Helium3};
  double bmass_max[SCPID_ntypes]={0.0, 0.05,  0.44, 0.85, 1.7,  2.65}; 
  double bmass_min[SCPID_ntypes]={0.05, 0.20, 0.55, 1.05, 2.0,  3.1};  
  double boff_max[SCPID_ntypes]={0.05,  0.04, 0.04, 0.04, 0.03, 0.02}; 
  double boff_min[SCPID_ntypes]={0.02, -0.04,-0.04,-0.04,-0.03,-0.02}; 
  double bp_max[SCPID_ntypes]=  {0.5,  3.0,  1.7,  3.0,  2.5,  1.8};  
  double bp_min[SCPID_ntypes]=  {0.01, 0.07, 0.12, 0.25, 0.35, 0.45}; 

  if(mass<-9) mass = getPDGmass(myid[itype],SCPIDin_myid);
  PIDmass2[itype]      = fabs(mass)*mass;
  beta_mass2_max[itype] = (mass_max>-9)? fabs(mass_max)*mass_max : fabs(bmass_max[itype])*bmass_max[itype];
  beta_mass2_min[itype] = (mass_min>-9)? mass_min*mass_min : bmass_min[itype]*bmass_min[itype];
  beta_off_max[itype] = (off_max>-9)? off_max : boff_max[itype];
  beta_off_min[itype] = (off_min>-9)? off_min : boff_min[itype];
  beta_p_max[itype]   = (p_max>-9)? p_max : bp_max[itype];
  beta_p_min[itype]   = (p_min>-9)? p_min : bp_min[itype];
}


//------------  getPDGmass (flag=SCPIDin_myid,SCPIDin_MCID,SCPIDin_GID)  --------------------------

double TSCtimeCorr::getPDGmass(int id, int flag) {
  int myid=PIDconvert(id, flag+10*SCPIDin_myid); 
  if(myid<0 || myid>=nMYPIDcodes) return -1.0; 
  //                            [0]=g, [1]=e+, [2]=e-, [3]=mu+, [4]=mu-, [5]=pi0, [6]=pi+, [7]=pi-,
  const double Mass[nMYPIDcodes]={0.0, 0.000511, 0.000511, 0.10566, 0.10566, 0.13498, 0.13957, 0.13957,
				// [8]=K0l, [10]=K+, [11]=K-, [12]=neutr, [13]=prot, [14]=antiprot, 
				0.497672, 0.493677, 0.493677, 0.9395656, 0.938272, 0.938272,
				//[15]=deuteron, [16]=triton, [17]=He3, [18]=He4
				1.875613, 2.80925, 2.80923, 3.727417 }; 
  return Mass[myid];
}

//------------  getBetaCalc (flag=SCPIDin_myid,SCPIDin_MCID,SCPIDin_GID)  --------------------------

double TSCtimeCorr::getBetaCalc(double mom, int id, int flag) {
  double mass = getPDGmass(id, flag);
  if(mass < 0.0) return -1.0;
  return mom/sqrt(mom*mom + mass*mass);
}
  
//------------  PID conversion: GID <-> MCID   ----------------------------------------
// flag for input PID class: SCPIDin_MCID (Lund, PDG), SCPIDin_GID (geant3 Id)

int TSCtimeCorr::PIDconvert(int id, int flag) {

  // [0]=g, [1]=e+, [2]=e-, [3]=mu+, [4]=mu-, [5]=pi0, [6]=pi+, [7]=pi-, [8]=K0l, [9]=K+, [10]=K-,
  // [11]=neutron, [12]=proton, [13]=antiproton, [14]=deuteron, [15]=triton, [16]=He3, [17]=He4
  static const int GID[nMYPIDcodes]={  1,  3,  2,  6,   5,  7,   8,  9,  10,  11,  12,  13,   14,   15,  45, 46, 49, 47}; 
  static const int MCID[nMYPIDcodes]={22,-11, 11, -13, 13, 111, 211,-211, 130, 321,-321,2112,2212,-2212, 45, 46, 49, 47}; 

  int ival=0;
  if(flag < 0) return ival;
  int flg1 = flag/10;
  switch(flag%10) {
  case SCPIDin_GID:
    {
      int i=0;
      while(GID[i]!=id && i<nMYPIDcodes) i++;
      if(i<nMYPIDcodes) {
	if(flg1>=SCPIDin_postype)
	  ival = PIDconvert(i, SCPIDin_myid+10*flg1);
	else if(flg1==SCPIDin_myid) 
	  ival=i;
	else
	  ival=MCID[i];
      }
      else if(flg1>0 && flg1!=SCPIDin_MCID)
	ival=-1;
    }
    break;
  case SCPIDin_MCID:
    {
      int i=0;
      while(MCID[i]!=id && i<nMYPIDcodes) i++;
      if(i<nMYPIDcodes) {
	if(flg1>=SCPIDin_postype)
	  ival = PIDconvert(i, SCPIDin_myid+10*flg1);
	if(flg1==SCPIDin_myid) 
	  ival=i;
	else
	  ival=GID[i];
      }
      else if(flg1>0 && flg1!=SCPIDin_GID)
	ival=-1;
    }
    break;
  case SCPIDin_myid:
    {
      if(flg1>=SCPIDin_postype) {
	ival = -1;
	if(id==mypid_Positron || id==mypid_Electron)       ival=SCPID_electron;
	else if(id==mypid_PionPlus || id==mypid_PionMinus) ival=SCPID_pion;
	else if(id==mypid_KaonPlus || id==mypid_KaonMinus) ival=SCPID_kaon;
	else if(id==mypid_Proton || id==mypid_AntiProton)  ival=SCPID_proton;
	else if(id==mypid_Deuteron)                       ival=SCPID_deuteron;
	else if(id==mypid_Helium3 || id==mypid_Triton)     ival=SCPID_helium3;
      }      
      else if(id>=0 && id<nMYPIDcodes) {
	if(flg1==SCPIDin_MCID)
	  ival=MCID[id];
	else if(flg1==SCPIDin_myid)
	  ival=id;
	else
	  ival=GID[id];
      }
    }
    break;
  case SCPIDin_postype:
    {
      if(id==0)      ival = mypid_Positron;
      else if(id==1) ival = mypid_PionPlus;
      else if(id==2) ival = mypid_KaonPlus;
      else if(id==3) ival = mypid_Proton;
      else if(id==4) ival = mypid_Deuteron;
      else if(id==5) ival = mypid_Helium3;
    //    else if(id==6) ival = id_Helium4;
      else ival = -1;
      if(flg1>SCPIDin_myid) ival = PIDconvert(ival, SCPIDin_myid+10*flg1);
    }
    break;
  case SCPIDin_negtype:
    {
      if(id==0)      ival = mypid_Electron;
      else if(id==1) ival = mypid_PionMinus;
      else if(id==2) ival = mypid_KaonMinus;
      else if(id==3) ival = mypid_AntiProton;
      else ival = -1;
      if(flg1>SCPIDin_myid) ival = PIDconvert(ival, SCPIDin_myid+10*flg1);
    }
  }
  return ival;
}

//-------------- dedt_Proton (time correction depending on energy deposit) ----------------------------
//
// dE dependent correction:   dedt_Proton=-yoff-yscal*pow(xde-xoff),xpow)
// sc_id=100*sector+paddle;  xde=energy deposit in SC paddle
// lookup arrays filled in case that parameters (yoff,yscal,xoff,xpow,xmin,xmax) are specified
//
double TSCtimeCorr::dedt_Proton(int scid, double xde, double yoff, double yscal, double xoff, double xpow, double xmin, double xmax) {

  int isc=SCpaddleId(scid);
  if(isc<0) return 0.0;
  if(yoff!=0.0) {
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
double TSCtimeCorr::dedt_Pion(int scid, double xde, double yoff, double yscal, double xoff, double xpow, double xmin, double xmax) {

  int isc=SCpaddleId(scid);
  if(isc<0) return 0.0;
  if(yoff!=0.0) {
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
    if(xde < dedtPion_xmin[isc]) {
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
void TSCtimeCorr::SCIDcut_init(int runno, double tg_center) {
  
  if(runno <= last_runno) return;
  if((SCIDcutChoice > 0) || (last_runno < 0)) {
    memset(SCid_Tcorr, 0, sizeof(SCid_Tcorr));
    memset(SC_dEdTcorr,  0, sizeof(SC_dEdTcorr));
  }
  last_runno = runno;
  if(tg_center>-400) target_center=tg_center;
  int i=-1;

  switch (SCIDcutChoice) {
  case SCIDcutG14b_Franz:
    {
      int iscbad[]={123, 323, 405, 505, -1}; //not switched off but not good for PIDcorr (maybe 305? 601?)
      int j=-1;
      while(iscbad[++j] > 0) SC_dEdTcorr[SCpaddleId(iscbad[j])]=16; 
      int scid_STedep1[]={134, 228, 229, 233, 323, 332, 416, 439, 520, 523, 535, 549, 605, 639, 646, -1};
      //429,630 dE too high
      j=-1;
      while(scid_STedep1[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep1[j])]=4;
      int scidno[]={116,  118,  124,  130,  134,  135,  152, 
		    205,  226,  229,  232,  233,  236,  239,  244,  245,  249,  251,  252,  255,  256, 257,
		    307,  311,  316,  323,  331,  332,  334,  335,  344, 349, 351, 352, 353,  354,  355,  356, 
		    405,  416,  424,  429,  439,  440,  448,  449,  450,  453,  455,  456,
		    523,  525,  530,  535,  549,  550,  551,  552,  553,  554,  555,
		    605,  623,  629,  630,  633,  639,  642,  649,  652,  653,  656, 657, -1};
      double off[]={0.2,  0.0,  1.9, 0.15,  4.3,  0.8,  4.9,
		    0.2, -12.,  4.0, -0.1,  2.1, -5.9,  4.2,-14.2, -19.2, -14., -10., 2.4,  1.9,  9.3, -1.3,
		    -0.2,-0.4, -0.2,  2.3, -9.8,  2.0, -12., -11.3, 3.3, 2.0, 1.0, 1.0,-13.6,-0.4,-10.2,-25.0,
		    3.5, 0.25,-18.7,-20.6, 0.45, 1.9, -24.4,-29.1, -0.5, -8.1,  1.0, -0.56,
		    0., -19., -18.9, 1.8,  0.0, -9.0, -11.,  20.2, 15.0, -14.5, -26.,
		    -0.5,-0.1,-20.7, -10., 0.20,  1.0, 2.15, -5.6, -3.6, -8.2,  7.7, 1.4};
      i=-1; 
      while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = off[i];
      //in general: STdedx at ~60% of g9 values?
      for(int k=0; k<SC_MAX_COUNTERS; k++) 
	SC_dEdTcorr[k] |=0x0800000; //scale STdep x1.4

      SC_dEdTcorr[SCpaddleId(130)] =0x0c09004; //SCdep x3.12, STdep x1.68
      SC_dEdTcorr[SCpaddleId(135)] =0x0c0c000; //SCdep x2.0, STdep x1.68
      if(runno>68180) {
	SC_dEdTcorr[SCpaddleId(106)] =0x0c08000;
	if(runno<69050) SC_dEdTcorr[SCpaddleId(130)] =0x0c18004; //SCdep x1.12, STdep x1.68
	int idlist1[]={150, 229,  239, 256, 257, 307, 323,  331, 334, 335, 344, 349, 352, 354, 356,
		       405, 424, 429, 440, 444, 455, 456, 520, 523, 525, 530, 535, 549, 601, 629, 630, 633, 646, 649, 656, -1};
	double cor1[]={-0.2,-0.2,-0.25,-0.6, 0.5,0.22,-2.3,-0.23,-0.2,-0.5,-1.3,-2.2,-0.2, 0.5,-0.8,
		       -4.0,-1.3, 0.4,-1.7,0.35,-1.1,-0.8, 2.0,1.0,-0.8,-1.0, 0.5, 2.4,-3.7,0.45, 0.2, 2.3, 0.4,-0.4, 0.58};
	i=-1;
	while (idlist1[++i] > 0) SCid_Tcorr[SCpaddleId(idlist1[i])] += cor1[i];
	if(runno<68422) SCid_Tcorr[SCpaddleId(555)] +=0.5;
	if(runno>68421) SCid_Tcorr[SCpaddleId(605)] -=1.2;
	if(runno>68421) SCid_Tcorr[SCpaddleId(323)] -=10.;
      }
      if(runno > 69050) { 
	SC_dEdTcorr[SCpaddleId(135)] =0x0c09004; //SCdep x3.12, STdep x1.68
	SC_dEdTcorr[SCpaddleId(607)] =0x0820010; //0.8,1.4,noPID
	int iscbad2[]={116, 305, 306, 416, 639, -1};//not good for PIDcorr
	while(iscbad2[++j] > 0) SC_dEdTcorr[SCpaddleId(iscbad[j])]=16;
	SC_dEdTcorr[SCpaddleId(130)] |=4; //use STdedx
	SC_dEdTcorr[SCpaddleId(440)] |=4; //use STdedx
	SCid_Tcorr[SCpaddleId(605)] = SC_COUNTER_OFF;
	int idlist2[]={101, 102, 103,  104, 106, 107, 108, 109, 112, 113, 115, 116, 117, 118, 119, 121, 122, 125, 128,
		       130,  131, 132, 134, 135, 136, 137, 138,  139,  140,  142,  143,  144, 145, 147, 149,
		       201, 205,  210, 213, 214,  215,  218,  222, 224,  225, 228, 229, 
		       231,  233, 238,  242, 244, 245, 248, 249, 252, 253, 
		       301, 302,  305, 307,  310, 311, 313, 316, 317, 320, 321, 323, 326, 
		       331, 332,  334,  335, 336, 337,  338, 341, 342, 351, 352, 
		       402,  408, 411, 413, 415, 418, 424, 425, 427, 431, 434, 435,  436, 439,
		       440, 441,  443, 444,  446, 447,  448, 450, 453, 456, 
		       503,  508,  509, 517, 519, 520, 522, 523, 524, 525,  527,  530, 531,  532,  536,
		       540, 541, 542,  549,  550, 551,   552,  553,  556, 
		       601,  602,  603, 606, 607, 608, 609, 611, 613, 615,  620, 621,  622, 624, 626, 628,  
		       633, 636, 639, 640, 641, 642,  647,  648,   649, 651,  654, 656,  657, -1};
	double cor2[]={0.16,0.2,-0.15,0.25, 0.2,-0.2,0.15,0.15,0.24, 0.3, 0.3, 2.2,0.15,0.7, 0.2,-0.3, 0.3, 0.3,-0.25,
		       -0.25,0.25,-0.2,-0.5,-1.3,-0.5,-0.3,-0.25,-0.45,-0.57,-0.62,-0.35,-0.2,-0.45,-0.6,-0.1,
		       -0.2,-0.25,-0.2,-0.4,-0.25,-0.35,-0.15,0.2,-0.25, 0.2, 0.3,-4.7,
		       -0.25,-0.3,-0.35, 0.3, 0.4, 0.1, 0.2,-0.1,0.25,0.25,
		       -0.2,-0.35,0.5,-0.35,-0.2,0.55,-0.2,0.25,0.1,-0.3,-0.4,10.0,0.43, 
		       -0.7,-0.7,-0.45,0.23,-0.33,-0.2,-0.6,-0.35,-0.3,-0.8,-0.6,
		       -1.45, 0.7,0.15, 0.1,-0.2,-0.2,0.25,-0.3,0.55,-0.2,-0.3,-0.25,-0.2,-0.42,
		       -0.7,-0.47,-0.4,-0.27, 0.2,-0.35, 0.7, 0.4,-0.25, 0.5,
		       -0.47,0.17,-0.2,-2.2,-0.3, 0.0, 0.0,-0.5,-0.7,-0.44,-0.66, 2.7,-0.15,-0.24,-0.17,
		       0.23,0.35,-0.52,-0.6,-0.25,-0.45, -20.,-15., -0.2, 
		       3.7,-0.2,-0.25,-0.2,-0.4,-0.3,-0.2,0.76,0.46,-0.12,-0.3,-0.25,-0.2, 0.2,-0.4,-0.15,
		       -1.0,-0.4,-1.9,0.17,-0.2,-0.25,-0.35,-0.55,0.32,-0.36,0.2, 0.0,-1.8};
	//332v: -2.0?//601v: or-0.7//654v: 0.72//656v:or-8.1
	i=-1;
	while (idlist2[++i] > 0) SCid_Tcorr[SCpaddleId(idlist2[i])] += cor2[i];
      }
      if(runno<69050)
	dedt_Proton(106,0.0,0.339766,-0.687186,17.6,0.1,17.7,70.);
      if(runno>68180) {
	if(runno<69050) {
	  dedt_Proton(110,0.0,0.0248533,0.00912621,10.4,1.25583,10.5,70.);
	  dedt_Proton(116,0.0,0.5,0.0133108,3.2,1.73798,3.3,70.);
	  dedt_Proton(117,0.0,3.0,-1.89163,0.8,0.19818,2,70.);
	  dedt_Proton(118,0.0,-0.8,0.1455341,3.79579,1.1444,3.89579,70.);
	  dedt_Proton(124,0.0,0.2,0.0086,14.4,1.11711,14.5,70.);
	}
	else { //i.e. runno>69050
	  dedt_Proton(109,0.0,0.2,0.0);
	  dedt_Proton(110,0.0,-1.36921,0.294143,-6.45066,0.483435,2,70.);
	  dedt_Proton(116,0.0,-1.96,0.195,-8.76348,0.906899,2,70.);
	  dedt_Proton(117,0.0,-0.381815,0.185048,8.32603,0.404365,8.42603,70.);
	  dedt_Proton(118,0.0,0.76,0.0048,3.6,1.75,3.7,70.);
	  dedt_Proton(120,0.0,0.077067,4.96195e-06,-13.3709,2.74608,2,70.);
	  dedt_Proton(122,0.0,0.05,0.004,5,1.50913,5.1,70.);
	  dedt_Proton(124,0.0,1.0,0.0);
	  dedt_Proton(125,0.0,-0.1,0.0165825,9.6,1.21733,9.7,70.);
	  dedt_Proton(126,0.0,0.055176,5.42093e-06,-4.22615,2.8609,2,70.);
	  dedt_Proton(127,0.0,0.1,0.00189935,12.8,1.70156,12.9,70.);
	}
      }
      else { //i.e. runno<68181
	dedt_Proton(116,0.0,-0.5,0.0133108,3.2,1.73798,3.3,70.);
	dedt_Proton(118,0.0,-1.016988,0.0855341,3.79579,1.1444,3.89579,70.);
	dedt_Proton(124,0.0,0.5,0.0086,14.4,1.11711,14.5,70.);
      }
      if(runno<69050) {
	dedt_Proton(125,0.0,0.2,0.038516,15.2,1.21981,15.3,70.);
	dedt_Proton(127,0.0,-2.99279,1.80896,2.85506,0.354969,2.95506,70.); 
	dedt_Proton(128,0.0,3.,-2.25722,8.4,0.109756,10.5,70.);
      }
      dedt_Proton(129,0.0,0.299897,3.98324e-05,17.6,2.42225,17.7,70.);
      if(runno<68422) dedt_Proton(130,0.0,1.1,0.0);
      dedt_Proton(132,0.0,-0.065937,9.77536e-06,14.423,2.82131,14.523,70.);
      dedt_Proton(133,0.0,0.765581,0.000393981,12,2.11545,12.1,70.); 
      if(runno>68180) 
	dedt_Proton(134,0.0,3.0,-1.21415,6.4,0.233603,6.5,70.);
      dedt_Proton(135,0.0,1.0,6.9079e-07,0.4,4.62707,6.5,70.);
      dedt_Proton(204,0.0,0.164403,0.00752692,9.6,1.24344,9.7,70.);
      if(runno<69050) {
	dedt_Proton(205,0.0,3,-2.36302,1.6,0.103254,2,70.);
	dedt_Proton(206,0.0,-0.000924096,0.0120375,1.6,1.19284,2,70.);
	dedt_Proton(207,0.0,-0.128709,0.0335274,-3.2,0.814941,2,70.);
      }
      else {
	dedt_Proton(205,0.0,0.4,0.0);
	dedt_Proton(206,0.0,-0.355259,0.0101376,-11.199,1.12939,2,70.);
	dedt_Proton(207,0.0,-1.08443,0.25426,-10.8204,0.480466,2,70.);
	dedt_Proton(208,0.0,-0.401746,0.0183789,-9.6,1.10099,2,70.);
	dedt_Proton(209,0.0,-1.46175,1.3182,11.2,0.1,11.3,70.);
     }
      if(runno>68180) { 
	dedt_Proton(211,0.0,-0.0430456,0.00683376,-7.47659,1.14556,2,70.);
	if(runno<69050) {
	  dedt_Proton(209,0.0,-0.029333,0.149366,9.6,0.384276,9.7,70.);
	  dedt_Proton(213,0.0,0.152925,0.0111915,11.2836,1.07715,11.3836,70.);
	  dedt_Proton(214,0.0,0.294864,0.0017581,12.8,1.57622,12.9,70.);
	  dedt_Proton(215,0.0,-0.0143694,0.0127837,-11.9995,0.956487,2,70.);
	  dedt_Proton(216,0.0,0.342308,5.25159e-05,12,2.18391,12.1,70.);
	  dedt_Proton(217,0.0,-0.467679,0.361379,-0.658485,0.304473,2,70.);
	  dedt_Proton(218,0.0,0.425061,0.00575405,12,1.19093,12.1,70.);
	  dedt_Proton(219,0.0,0.390069,0.00305787,13.6,1.43754,13.7,70.);
	  dedt_Proton(220,0.0,-0.256008,0.0634041,-12.8,0.694375,2,70.);
	  dedt_Proton(221,0.0,0.430322,0.00114637,13.6,1.53328,13.7,70.);
	}
	else { // i.e. runno>69050
	  dedt_Proton(212,0.0,-0.341629,0.00322355,-13.4985,1.3837,2,70.);
	  dedt_Proton(214,0.0,-0.1423572,0.00825564,10.6392,1.17887,10.7392,70.);
	  dedt_Proton(216,0.0,-0.935355,0.644466,9.59757,0.247819,9.69757,70.);
	  dedt_Proton(217,0.0,-0.179364,0.000466564,-12.7962,1.8431,2,70.);
	  dedt_Proton(218,0.0,-0.405134,0.0146889,-12.7999,1.0613,2,70.);
	  dedt_Proton(219,0.0,0.187068,0.000883086,16.8,1.75139,16.9,70.);
	  dedt_Proton(220,0.0,0.141895,0.00199367,13.6,1.66842,13.7,70.);
	  dedt_Proton(221,0.0,0.180322,0.00114637,13.6,1.53328,13.7,70.);
	}
      }
      dedt_Proton(223,0.0,0.275846,7.36779e-05,0.8,2.20614,2,70.);
      if(runno<69050) { 
	dedt_Proton(222,0.0,-1.06897,0.109409,-1.59891,0.99796,2,70.);
	dedt_Proton(224,0.0,0.28,0.00192978,2.4,1.49997,2.5,70.);
      }
      else { 
      	dedt_Proton(222,0.0,-0.50293,0.14,-9.59827,0.893405,2,70.);
	dedt_Proton(224,0.0,-0.056118,0.01,10,0.0,10,70.);
      }
      if(runno>68180) { 
	dedt_Proton(225,0.0,0.412086,0.000426269,16,1.67714,16.1,70.);
	if(runno<69050) 
	  dedt_Proton(226,0.0,0.199293,0.0119975,-16.6788,0.749856,2,70.);
	else
	  dedt_Proton(226,0.0,12.0,0.0);
      }
      dedt_Proton(227,0.0,0.435195,1.51339e-07,3.16513,3.46945,3.26513,70.);
      dedt_Proton(228,0.0,3,-2.00159,2.4,0.161822,2.5,70.);
      if(runno<69050) {
	dedt_Proton(229,0.0,0.6,0.021764,-3.98111,1.24822,2,70.);
	dedt_Proton(230,0.0,0.03501,0.025188,12,1.35777,12.1,70.);
	dedt_Proton(231,0.0,0.22,0.00193528,11.2,1.92824,11.3,70.); 
      }
      else
	dedt_Proton(229,0.0,-2.2,0.0);
      dedt_Proton(232,0.0,0.179841,1.17926e-06,17.6,3.1882,17.7,70.);
      dedt_Proton(234,0.0,2.7,-2.473,0.8,0.1,2,70.);
      dedt_Proton(235,0.0,1.0,5.35193e-07,23.8657,3.57663,23.9657,70.);
      dedt_Proton(246,0.0,-0.289359,0.0198042,3.2,1.21912,3.3,70.);
      dedt_Proton(247,0.0,-0.289359,0.0198042,3.2,1.21912,3.3,70.);
      dedt_Proton(248,0.0,-0.180264,0.00390619,1.6,1.64098,2,70.);
      dedt_Proton(250,0.0,-0.128981,0.00898725,5.6,1.50674,5.7,70.);
      if(runno>69050) {
	dedt_Proton(301,0.0,-0.11562,0.0052724,3.2,1.62533,3.3,70.);
	dedt_Proton(304,0.0,-0.19932,0.0409465,5.43405,1.12163,5.53405,70.);
	dedt_Proton(305,0.0,-2.2,0.16016,-1.5998,0.985968,2,70.);
	dedt_Proton(306,0.0,2.0,0.0);
      }
      else {
	dedt_Proton(305,0.0,-1.21446,0.16016,-1.5998,0.985968,2,70.);
	dedt_Proton(306,0.0,-0.258798,0.00709846,-2.39999,1.60919,2,70.);
      }
      if(runno>68180) {
	dedt_Proton(308,0.0,0.235914,0.000321488,22.4,2.07878,22.5,70.);
	dedt_Proton(312,0.0,-0.0842454,0.0462802,-3.40217,0.661343,2,70.);
      }
      dedt_Proton(309,0.0,0.184351,0.0105604,19.2,1.25885,19.3,70.);
      if(runno<69050) 
	dedt_Proton(313,0.0,0.0210356,0.00559967,-6.427,1.21496,2,70.);
      else
	dedt_Proton(313,0.0,-0.536469,0.496381,12.8,0.100001,12.9,70.);
      if(runno>68180 && runno<69050) {
	dedt_Proton(314,0.0,0.184943,0.00227853,12.8,1.57594,12.9,70.);
	dedt_Proton(315,0.0,-0.184435,0.0335196,-14.3999,0.799217,2,70.);
	dedt_Proton(318,0.0,0.1,0.000822703,13.6,1.68086,13.7,70.);
	dedt_Proton(319,0.0,0.389325,0.000280599,12.8,1.91346,12.9,70.);
	dedt_Proton(320,0.0,0.374887,0.00169897,12.8,1.49971,12.9,70.);
	dedt_Proton(322,0.0,0.335077,0.0191201,14.3997,1.23334,14.4997,70.);
      }
      else if(runno>69050) {
	dedt_Proton(314,0.0,-0.0968975,9.08344e-05,2.05684,2.32854,2.15684,70.);
	dedt_Proton(315,0.0,-0.102969,0.119238,12.8704,0.430942,12.9704,70.);
	dedt_Proton(316,0.0,-0.27439,0.247299,2.2814,0.443767,2.3814,70.);
	dedt_Proton(317,0.0,-0.609198,0.454691,9.59999,0.320583,9.69999,70.);
	dedt_Proton(318,0.0,-0.1128,0.006,4.3488,1.26802,4.4488,70.);
	dedt_Proton(319,0.0, 0.10124,0.000193167,14.4,1.98199,14.5,70.);
	dedt_Proton(320,0.0,-0.823519,0.176302,1.49159,0.503761,2,70.);
	dedt_Proton(322,0.0,-0.62,0.0023,-16.3883,1.58865,2,70.);
      }
      if(runno<69050) 
	dedt_Proton(321,0.0,0.1,0.00277524,3.2,1.45354,3.3,70);
      if(runno<68181) 
	dedt_Proton(324,0.0,0.29,0.00203975,14.4,1.69673,14.5,70.);
      else
	dedt_Proton(324,0.0,0.22,0.0);
      dedt_Proton(325,0.0,0.186464,4.95448e-06,4,2.96137,4.1,70.);
      dedt_Proton(346,0.0,-0.390492,0.0133108,3.2,1.33798,3.3,70.);
      if(runno<69050) 
	dedt_Proton(402,0.0,-2.91515,1.62134,0.499403,0.28756,2,70.);
      else {
	dedt_Proton(402,0.0,-0.75,0.0);
	dedt_Proton(405,0.0,-5.3,3.51378e-05,-1.23233,4.28107,2,70.);
	dedt_Proton(408,0.0, 0.6,0.0);
	dedt_Proton(416,0.0,-1.87687,0.56184,6.23077,0.667429,6.33077,70.);
	dedt_Proton(417,0.0,-0.25724,0.00601316,-12.1981,1.16911,2,70.);
	dedt_Proton(418,0.0,-0.772378,0.106651,-12.7946,0.584051,2,70.);
      }
      if(runno>68180 && runno<69050) 
	dedt_Proton(422,0.0,0.04238,0.00598497,13.6,1.5067,13.7,70.);
      if(runno<69050)
	dedt_Proton(423,0.0,0.192459,0.00439352,6.4,1.68668,6.5,70.);
      else {
	dedt_Proton(422,0.0,-0.2449111,0.000254904,-15.9997,2.00524,2,70.);
	dedt_Proton(423,0.0,-0.215851,0.000704815,-4.65782,1.92677,2,70.);
      }
     if(runno<68180) 
	dedt_Proton(424,0.0,-0.279176,0.036537,-2.4,0.910531,2,70.);
      if(runno<69050) 
	dedt_Proton(428,0.0,0.810444,0.0596767,15.2,1.20456,15.3,70.);
      else 
	dedt_Proton(428,0.0,0.108183,0.00249084,17.6,1.61283,17.7,70.);
      dedt_Proton(451,0.0,-0.224098,0.00399342,1.6,1.6371,2,70.);
      if(runno<69050) 
	dedt_Proton(502,0.0,-0.322486,0.0567768,5.6,1.33922,5.7,70.);
      else {
	dedt_Proton(502,0.0,0.0823802,8.3227e-05,-9.6,2.45129,2,70.);
	dedt_Proton(503,0.0,-0.21,0.0);
      }
      dedt_Proton(505,0.0,-2.73895,0.0974445,8.38352,1.3496,8.48352,70.);
      //	if(runno<69050) dedt_Proton(506,0.0,-0.253216,0.0249845,0.8,1.45487,2,70.);
      if(runno<68181) {
	dedt_Proton(517,0.0,-0.95,0.00996271,1.6,1.52977,2,70.);
	dedt_Proton(524,0.0,-1.7,0.470896,-2.82986,0.509576,2,70.);
      }
      else { //i.e. runno>68180
	if(runno>69050) {
	  dedt_Proton(517,0.0,-1.2,0.0);
	  dedt_Proton(518,0.0,-0.126859,0.071055,13.4883,0.546451,13.5883,70.);
	  dedt_Proton(520,0.0, 0.75,0,0); 
	  dedt_Proton(522,0.0,-2.98647,0.274102,-11.9534,0.793694,2,70.);
	  dedt_Proton(524,0.0,-2.2,0.470896,-2.82986,0.509576,2,70.);
	}
	else {
	  dedt_Proton(517,0.0,-1.9,0.300084,4.85728,0.796582,4.95728,70.);
	  dedt_Proton(522,0.0,0.159836,0.00609667,14.4,1.59971,14.5,70.);
	  dedt_Proton(524,0.0,-1.2,0.470896,-2.82986,0.509576,2,70.);
	}
      }
      dedt_Proton(519,0.0,-0.0643748,0.000269477,-1.6,1.93514,2,70.);
      if(runno<69050) 
	dedt_Proton(520,0.0,2.7,-1.82306,1.6,0.1,2.,70.);
      if(runno>69050) 
	dedt_Proton(521,0.0,-0.189384,0.0192006,2.75081,1.09703,2.85081,70.);
      dedt_Proton(523,0.0,0.9,-0.001,-0.8,1.,2.,70.);
      if(runno<68188) 
	dedt_Proton(529,0.0,-0.35,0.00138978,-6.21869,1.68372,2,70.);
      dedt_Proton(531,0.0,0.38,9.80106e-08,12.8,3.99115,12.9,70.);
      dedt_Proton(533,0.0,0.05,0.000965426,0.8,1.68532,2,70.);
      dedt_Proton(535,0.0,-0.3,0.0);
      if(runno<69050) 
	dedt_Proton(601,0.0,-1.81064,0.0142246,-15.1993,1.37558,2,70.);
      else
	dedt_Proton(601,0.0,1.3,0.0);//601v: no corr
      dedt_Proton(602,0.0,0.02,0.000726521,9.6,1.84287,9.7,70.);
      dedt_Proton(604,0.0,0.1205,0.000357013,2.4,2.40479,2.5,70.);
      if(runno<69050) {
	dedt_Proton(607,0.0,-0.332136,0.0513749,1.6,1.23628,2,70.);
	dedt_Proton(613,0.0,-2.1,0.269007,3.07733,1.04175,3.17733,70.);
	dedt_Proton(614,0.0,-2.0,0.180144,2.17871,1.10602,2.27871,70.);
      }
      else { //i.e. runno>69050
	dedt_Proton(611,0.0,0.32,0.0);
	dedt_Proton(612,0.0,-0.742536,0.0253143,-11.9989,0.97483,2,70.);
	dedt_Proton(613,0.0,0.1,0.063185,1.92354,1.1064,2,70.);
	dedt_Proton(614,0.0,-0.35,0.198917,5.812,0.58774,5.912,70.);
	dedt_Proton(615,0.0,-0.93882,0.107693,-12.0107,0.607528,2,70.);
	dedt_Proton(616,0.0,-0.605878,0.0268689,-12.6021,0.880928,2,70.);
	dedt_Proton(617,0.0,-0.5844,0.134058,0.79737,0.554473,2.89737,70.);
	dedt_Proton(618,0.0,-0.2641416,0.000332653,-14.4,1.9085,2,70.);
	dedt_Proton(619,0.0,0.016689,0.000108355,2.10432,2.25211,2.20432,70.);
	dedt_Proton(620,0.0,-0.4,0.00322193,-14.3999,1.34424,2,70.);
	dedt_Proton(621,0.0,-0.067087,0.000479573,12.501,2.01732,12.601,70.);
      }
      if(runno>68180 && runno<69050) {
	dedt_Proton(612,0.0,0.0131479,0.00501798,-11.1954,1.10857,2,70.);
	dedt_Proton(615,0.0,-0.902466,0.152379,-12.2732,0.553558,2,70.);
	dedt_Proton(616,0.0,-0.554589,0.0945652,-12.7345,0.618711,2,70.);
	dedt_Proton(617,0.0,-0.41005,0.0275619,-10.9866,0.924217,2,70.);
	dedt_Proton(618,0.0,-0.173394,0.0017677,-14.4,1.53019,2,70.);
	dedt_Proton(619,0.0,0.232655,0.00248399,11.2,1.46466,11.3,70.);
	dedt_Proton(621,0.0,0.159603,0.00235898,15.2,1.57374,15.3,70.);
      }
      dedt_Proton(622,0.0,0.120788,0.000157589,4.8,2.11884,4.9,70.);
      if(runno<68181) {
	dedt_Proton(623,0.0,0.2,0.0187258,1.46817,1,2,70);
	dedt_Proton(624,0.0,0.3,0.00464152,0.8,1.37912,2,70.);
      }
      else { //i.e. runno>68180
	if(runno<68422)
	  dedt_Proton(623,0.0,0.8,0.0187258,1.46817,1,2,70);
	else
	  dedt_Proton(623,0.0,-0.3,0.0187258,1.46817,1,2,70);
	dedt_Proton(624,0.0,0.01,0.00464152,0.8,1.37912,2,70.);
      }
      dedt_Proton(625,0.0,0.16,0.000942636,-2.39995,1.56677,2,70.);
      if(runno<69050) { 
	dedt_Proton(626,0.0,0.01,0.0393641,-3.12252,0.769709,2,70.);
	dedt_Proton(631,0.0,0.4,0.00819722,14.4,1.40512,14.5,70.);
      }
      dedt_Proton(627,0.0,0.3,8.0783e-06,3.99809,2.65945,4.09809,70.);
      dedt_Proton(628,0.0,0.5,9.36724e-06,11.2,2.64661,11.3,70.);
      dedt_Proton(632,0.0,0.652094,0.000197481,17.6,2.1167,17.7,70.);
      dedt_Proton(633,0.0,1.12797,0.00178548,17.5946,1.76987,17.6946,70.);
      if(runno>69050)
	dedt_Pion(639,0.0,-0.22,0.0);
      if(runno>68180 && runno<69050)
	dedt_Pion(106,0.0,-2.22898,3.2,0.142409,3.3,70.);
      dedt_Pion(128,0.0,2.09389,-1.6775,2.4,0.1,2.5,70.);
      dedt_Pion(133,0.0,-0.2,0.032,3.99413,0,2,70.);
      dedt_Pion(134,0.0,0.488595,-0.109491,6.56829,0,2,70.);
      dedt_Pion(135,0.0,-2.3,1.99156,0.65841,0.180996,2,70.);
      dedt_Pion(144,0.0,-0.208654,0.0288841,6.96856,0,2,70.);
      dedt_Pion(146,0.0,-0.508029,0.0278119,3.09441,1.16445,3.19441,70.);
      dedt_Pion(149,0.0,-0.203037,0.0123528,3.2,1.31378,3.3,70.);
      dedt_Pion(150,0.0,-0.553017,0.0473225,8.69626,0,2,70.);
      dedt_Pion(151,0.0,-0.344003,0.0508589,8.65295,0,2,70.);
      dedt_Pion(157,0.0,-0.262501,0.0490657,4.73149,0,2,70.);
      dedt_Pion(222,0.0,-0.282021,0.00580945,0.8,1.84357,2,70.);
      if(runno<69050) {
	dedt_Pion(230,0.0,-0.533388,0.0810058,3.56382,0,2,70.);
	dedt_Pion(231,0.0,-0.453431,0.0436925,0.882027,0,2,70.);
      }
      dedt_Pion(235,0.0,-0.120471,0.0118643,1.79287,0,2,70.);
      dedt_Pion(236,0.0,0.239422,1.17916e-06,-6.20444,3.14362,2,70.);
      dedt_Pion(246,0.0,-0.370072,0.0223605,5.6,1.26078,5.7,70.);
      dedt_Pion(247,0.0,-0.289359,0.0198042,3.2,1.21912,3.3,70.);
      dedt_Pion(248,0.0,-0.180264,0.00390619,1.6,1.64098,2,70.);
      dedt_Pion(249,0.0,-0.592665,0.035359,4.00887,0,2,70.);
      dedt_Pion(250,0.0,-0.128981,0.00898725,5.6,1.50674,5.7,70.);
      dedt_Pion(251,0.0,-0.452382,0.00766532,2.4,1.56103,2.5,70.);
      dedt_Pion(305,0.0,-1.10374,0.0223994,-2.4,1.58119,2,70.);
      dedt_Pion(346,0.0,-0.390492,0.0133108,3.2,1.33798,3.3,70.);
      dedt_Pion(347,0.0,-0.48464,0.0263562,-3.12613,1.01845,2,70.);
      dedt_Pion(348,0.0,-0.39176,0.0238184,1.6,1.21057,2,70.);
      dedt_Pion(350,0.0,-0.359504,0.0115457,1.6,1.36009,2,70.);
      if(runno<69050) 
	dedt_Pion(402,0.0,-3.0,1.80603,4.21093,0.287011,4.31093,70.);
      dedt_Pion(423,0.0,-0.371619,0.0121308,1.6,1.44795,2,70.);
      if(runno<69050) 
	dedt_Pion(428,0.0,-0.74231,0.0530278,2.4,1.24907,2.5,70.);
      else 
	dedt_Pion(428,0.0,-0.266846,0.000676469,2.4,1.86849,2.5,70.);
      dedt_Pion(447,0.0,-1.07251,0.0426144,-8.79932,1.02244,2,70.);
      dedt_Pion(451,0.0,-0.1,0.00399342,1.6,1.6371,2,70.);
      dedt_Pion(452,0.0,-0.930448,0.0190526,-8.67281,1.2508,2,70.);
      dedt_Pion(454,0.0,-0.2,0.00832288,-3.33961,1.40149,2,70.);
      dedt_Pion(456,0.0,-0.206692,2.22695e-06,-7.2,3.08648,2,70.);
      if(runno<69050) 
	dedt_Pion(502,0.0,-0.729432,0.0362269,2.4,1.44565,2.5,70.);
      else
	dedt_Pion(502,0.0,-0.00501296,0.0180926,2.4,1.32062,2.5,70.);
      if(runno<68180) 
	dedt_Pion(517,0.0,-1.0,0.0781069,1.72179,0,2,70.);
      else
	dedt_Pion(517,0.0,-0.6,0.0781069,1.72179,0,2,70.);
      dedt_Pion(520,0.0,2.5,-1.85405,1.6,0.1,2,70.);
      dedt_Pion(524,0.0,-3.0,2.10873,1.43705,0.159743,2,70.);
      dedt_Pion(539,0.0,-0.782373,0.13418,-3.2,0.612977,2,70.);
      dedt_Pion(546,0.0,-0.318555,0.00786176,1.6,1.42305,2,70.);
      dedt_Pion(547,0.0,-0.776792,0.0182674,-6.39944,1.22746,2,70.);
      dedt_Pion(556,0.0,-0.558956,0.0608207,4.62534,0,2,70.);
      if(runno<69050) { 
	dedt_Pion(601,0.0,-1.81064,0.0142246,-15.1993,1.37558,2,70.);
	dedt_Pion(607,0.0,-0.50582,0.0025939,-3.99694,1.95829,2,70.);
	dedt_Pion(613,0.0,-2.65,0.0520837,-4.78794,1.43585,2,70.);
	dedt_Pion(614,0.0,-1.90149,0.0312718,-3.1916,1.56516,2,70.);
      }
      else {
	dedt_Pion(613,0.0,-2.97535,0.167841,-4.79996,1.0891,2,70.);
	dedt_Pion(614,0.0,-0.0276313,0.0434743,2.4,1.13338,2.5,70.);
      }
      dedt_Pion(628,0.0,-1.0769,0.77901,-0.801587,0.132035,2,70.);
      if(runno<69050) {
	dedt_Pion(631,0.0,-0.656054,0.0500413,-3.2,0.960601,2,70.);
	dedt_Pion(632,0.0,-0.185571,0.0174655,4.75556,0,2,70.);
      }
      else
	dedt_Pion(639,0.0,-1.0,0.0);
      dedt_Pion(650,0.0,-0.194173,0.00266067,4,1.75718,4.1,70.);
      dedt_Pion(653,0.0,-0.39933,0.127391,4.2312,0,2,30.);
    }
    break;

  case SCIDcutG14b_Haiyun:
    {
      int badscid[]={116, 118, 123, 125, 127, 222, 230, 231, 305, 323, 402, 405, 428, 502, 505, 517, 601, 607, 613, 614, -1};
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
      int scidno[]={124,  134,  144,  146,  147,  148,  149,  150,  151,  152, 
		    226,  229,  233,  236,  239,  244,  245,  246,  247,  248,  249,  250,  251,  252, 
		    331,  332,  334,  335,  344,  
		    424,  429,  446,  447,  448,  449,  450,  451,
		    520,  525,  530,  535,  546,  547,  548,  549,  550,  551,  552,  553, 
		    605,  629,  630,  633,  642,  647,  649,  650,  651,  652, -1};
      //for silver data
      double offs[]={1.82, 3.56,-1.34,-1.65,-1.23,-2.11,-1.69,-1.86, -2.2,  3.4, 
		   -12.32,3.16,2.14,-6.68, 3.57,-14.23,-19.48,-2.0,-1.85,-1.7,-15.7, -1.4,-11.9, 1.1,
		   -10.1,1.89,-11.93,-11.99,1.27,
		   -20.14,-20.28,-1.6,-1.7,-25.8,-30.77,-2.0, -2.3,
		   2.06,-19.72,-19.58,2.46,-1.16,-1.3,-1.39, 1.59, -9.5,-12.9, 19.0,  14.0,
		   -1.54,-20.38,-9.94,2.32, 2.15,-1.35,-7.9, -1.8, -2.2, -6.0};
      //for gold data
      double offg[]={1.62, 2.24,-1.75,-1.45,-1.21,-2.0, -1.72,-1.96,-1.58, 2.12, 
		   -12.33, 2.25, 1.26,-6.42,3.57,-13.94,-19.48,-1.88,-1.86,-1.9,-15.34, -1.5,-11.33, 0.78,
		   -10.96, 0.85,-12.35,-11.99,1.92,
		   -19.66,-20.15,-0.92,-1.54,-25.54,-30.83,-1.57, -2.5,
		    2.17,-20.22,-16.84,2.295, -1.02, -1.46,-1.00, 2.08, -9.68,-12.74, 17.73, 14.0,
		    0.0, -20.42, -9.92, 0.84, 1.78,-1.25,-7.7, -1.8, -1.51, -6.0};
      i=-1; 
      if(runno < 69100) {
	while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = offs[i];
      }
      else {
	while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = offg[i];
      }
    }
    break;

  case SCIDcutG9a_Steffen:
    {
      int badscid[]={245, 423, 523, 552, 553, 644, 646, -1};
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
      SCid_Tcorr[SCpaddleId(555)] = -22.0;
      if(runno >= 56164 && runno <= 56233) SCid_Tcorr[SCpaddleId(136)] = -2.0;
      if(runno >= 55521 && runno <= 55676) SCid_Tcorr[SCpaddleId(344)] = SC_COUNTER_OFF;
      if(runno >= 55669 && runno <= 55676) {
	SCid_Tcorr[SCpaddleId(443)] = -4.0;
	SCid_Tcorr[SCpaddleId(449)] = SCid_Tcorr[SCpaddleId(453)] = 2.0;
      }
      if(runno >= 55664 && runno <= 55668) SCid_Tcorr[SCpaddleId(649)] = 2.0;
      if(runno >= 55521 && runno <= 55676) SCid_Tcorr[SCpaddleId(655)] = -2.0;
      SCid_Tcorr[SCpaddleId(656)] = 1.0;
    }
    break;

  case SCIDcutG9a_Liam:
    {
      int badscid[]={110, 111, 112, 113, 114, 115, 116, 156, 212, 244, 245, 330, 344, 
		     423, 433, 436, 522, 523, 529, 552, 553, 555, 601, 626, 630, 644, 646, 654, 656, -1};
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
    }
    break;
  case SCIDcutG9b_Natalie:
    {
      int badscid[]={224, 244, 245, 249, 251, 322, 338, 415, 448, 449,453, 522,550, 551,552, 553, 554, 555, 
		   612, 613, 649, 653, 656, -1};
      while (badscid[++i] > 0) SCid_Tcorr[SCpaddleId(badscid[i])] = SC_COUNTER_OFF;
    }
    break;
  case SCIDcutG9b_Aneta: 
  case SCIDcutG9b_Franz: 
    if(runno > 62717) { // && runno < 62987?) {
      int scidno[]= { 502, 405, 416, 323, 226, 633, 134, 239,
		      540, 244, 444, 245, 347, 448, 249, 449, 549, 649, 550, 251, 551, 152, 252, 
		      552, 652, 353, 453, 553, 653, 254, 354, 454, 554, 555, 355, 256, 356, 656, 457, 157,-1};
      double scoff[]={-1.1, 0.0, 2.3, 0.7,-12., 0.7, 0.0,  2.,
		     -0.5,-14., 0.0,-20.,-0.2,-24.,-15.,-29.5,-1.2, -5.,-10.,-11.,-11., 4., 2.,
		     19.2,-4.5,-14.5,-10., 14.,-14.,-0.5,-0.5,-0.5,-16.3,-26.5,-11.,10.,-27., 7.5,-1.2,-0.7};
      //554: 628:-16 653?b:0.45
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
      double scoff[]={-0.1,-0.3,-1.0,-1.0, 2.5, 0.5,-2.0,-14., 0.2, 0.4,-0.4,-1.0,-0.5,-0.7,-0.2,
		     -0.4, 2.5,-0.8,-0.4,-14.,-0.2,-20.,-0.5,-24.5,-30.,-15.5,0.6,-5.0,-11.5,-12.2,
		     19.5, 14.,-16.,-26.5,-10.2,4.0,-5.,-9.0,-14.5,-14.,-0.5,-0.5,-0.4,-0.2,
		     -11.,  8.,-27.5,7.5,-0.8,-0.5,-0.5, -0.6}; 
      while (scidno[++i] > 0) SCid_Tcorr[SCpaddleId(scidno[i])] = scoff[i];
      if(runno>62604) {
	int scnodiff[]={ 123, 124, 133, 134, 140, 141, 142, 143, 144, 147, 152, 229, 233, 236, 237, 245,
			 323, 356, 416, 530, 542, 552, 652, 653, 657,-1};
	double scdiff[]={0.2, 1.3,-0.15,0.4,-1.7,-0.7,-1.7,-2.0,-0.3,-0.5, 0.2, 1.2, 0.2, 1.6, 0.8, 0.5,
			 0.2,-0.4,-0.8, 0.6,-0.2,-0.5, 0.4, 0.3,0.22};
	i=0;
	while (scnodiff[++i] > 0) SCid_Tcorr[SCpaddleId(scnodiff[i])] += scdiff[i];
      }
      else {
	if(runno>62500) {
	  SCid_Tcorr[SCpaddleId(144)] += 0.6;
	  SCid_Tcorr[SCpaddleId(236)] += 0.4;
	}
	if(runno<62289) {
	  SCid_Tcorr[SCpaddleId(344)]= 0.5;
	  //	SCid_Tcorr[SCpaddleId(554)] =-14.0;  //??
	}
	else {
	  //	SCid_Tcorr[SCpaddleId(554)] =-14.0;  //??
	  if(runno<62500) SCid_Tcorr[SCpaddleId(155)] -= 0.4;
	}
	if(runno>62400) {
	  if(runno<62500)  //??
	    SCid_Tcorr[SCpaddleId(237)] -= 0.5;
	  SCid_Tcorr[SCpaddleId(249)] += 0.7;
	  SCid_Tcorr[SCpaddleId(251)] += 2.7;
	  SCid_Tcorr[SCpaddleId(353)] += 1.5;
	}
	if(runno>62564) { //??
	  SCid_Tcorr[SCpaddleId(617)]= 0.7;
	  SCid_Tcorr[SCpaddleId(618)]= 1.0;
	  SCid_Tcorr[SCpaddleId(621)]= -1.3;
	}
      }
    }
    break;
  }

  if(SCIDcutChoice != SCIDcutG9b_Franz) return;

  //set flag to use STdedx for time correction functions
  int scid_STedep[]={ 416, 323, 326, 429, 629, 134, 234, 236, -1}; //123
  int j=-1;
  while(scid_STedep[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep[j])]=4;
  SC_dEdTcorr[SCpaddleId(123)]|=0x04000;
  SC_dEdTcorr[SCpaddleId(144)]|=0x0c000;
  if(runno<62604)
    SC_dEdTcorr[SCpaddleId(229)]=4;
  if(runno<62718) {
    int scid_STedep1[]={ 523, 530, 331, 233, 416, 444, 347, -1};
    j=-1;
    while(scid_STedep1[++j] > 0) SC_dEdTcorr[SCpaddleId(scid_STedep1[j])]=4;
    if(runno>62604) {
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
  if(runno<62605) dedt_Proton(123,0.0, -0.117054,0.00850427,1.6,1.9268,2,70);
  dedt_Proton(223,0.0, 0.0979165,0.000503644,3.2,1.95479,3.3,70.);
  dedt_Proton(323,0.0, -0.981879,0.0139318,0.8,2.77031,2,70);
  dedt_Proton(423,0.0, -0.254501,0.0172978,4,1.1982,4.1,70.);
  dedt_Proton(623,0.0, 0.0332715,0.0187258,1.46817,1,2,70);
  if(runno<62500) // || (runno>62604 && runno<62828)) //???
    dedt_Proton(124,0.0, -0.0298794,0.0141087,6.4,1.31286,6.5,70.);
  else if(runno>62604 && runno<62828) 
    dedt_Proton(124,0.0, 0.043791,0.00353826,18.4,1.48102,18.5,70.);
  if(runno<62500 || runno>62604)
    dedt_Proton(224,0.0, 0.050504,0.000523787,5.6,2.05784,5.7,70.);
  dedt_Proton(324,0.0, 0.023401,0.0024826,1.6,1.55345,2,70.);
  dedt_Proton(524,0.0, 0.070833,0.00282938,4,1.57203,4.1,70.);
  dedt_Proton(624,0.0, 0.029032,0.00166275,0.8,1.76111,2,70.);
  if(runno>62604 && runno<62718)
    dedt_Proton(125,0.0, -0.06,0.0085,6.4,1.31379,6.5,70.);
  else
    dedt_Proton(125,0.0, -0.161318,0.0165393,6.4,1.31379,6.5,70.);
  dedt_Proton(425,0.0, -0.376865,1.50562e-06,6.4,3.32393,6.5,70.);
  dedt_Proton(525,0.0, 0.139603,0.000841331,2.4,1.83562,2.5,70.);
  dedt_Proton(226,0.0, -0.0907017,0.00192274,4,1.51669,4.1,70.);
  if(runno<62500 || runno>62717) 
    dedt_Proton(127,0.0, -0.0581778,0.00416721,3.2,1.77332,3.3,70.);
  if(runno>62604 && runno<62718) {
    dedt_Proton(326,0.0, -1.2,0.0);
    dedt_Proton(426,0.0, -0.4,0.0);
    dedt_Proton(127,0.0, -0.0615916,0.00978456,5.6,1.44812,5.7,70.);
    dedt_Proton(427,0.0, -0.2,0.0);
    dedt_Proton(228,0.0,  0.05,0.0052189,-1.2,1.2,2,70.); //??
  }
  else { 
    dedt_Proton(427,0.0, -0.216233,0.000160107,8.8,1.9781,8.9,70.);
    dedt_Proton(228,0.0, -0.351308,0.0352189,-1.21466,1,2,70.);
  }
  dedt_Proton(528,0.0, 0.14166,0.000282926,-8.8,2.00977,2,70.);
  if(runno<62605) 
    dedt_Proton(229,0.0, -0.00994333,0.0201451,-2.27297,1,2,70.);
  else
    dedt_Proton(229,0.0, -0.816903,0.0570301,2.67814,1,2.77814,70.);
  dedt_Proton(529,0.0, -0.00371472,0.00138978,-6.21869,1.68372,2,70.);
  if(runno<62987)
    dedt_Proton(629,0.0, -1.13568,0.209839,-0.472198,1,2,70.);
  if(runno>62604 && runno<62718) {
    dedt_Proton(230,0.0, 0.43052,0.00117186,15.2,1.70164,15.3,70.);
    dedt_Proton(429,0.0, -0.8,0.0);
    dedt_Proton(430,0.0, -0.5,0.0);
    dedt_Proton(431,0.0, -0.4,0.0);
    dedt_Proton(432,0.0, -0.5,0.0);
    dedt_Proton(433,0.0, 0.5,-3.2,2.4,0.277499,2.5,70.); //x4: -0.5,0.0)
    dedt_Proton(434,0.0, -0.5,0.0);
    dedt_Proton(435,0.0, -0.6,0.0);
    dedt_Proton(436,0.0, -0.5,0.0);
    dedt_Proton(437,0.0, -0.5,0.0);
    dedt_Proton(438,0.0, -0.7,0.0);
  }
  else {
    dedt_Proton(230,0.0, -0.335781,0.000142584,-14.6622,2.34831,2,70.);
    dedt_Proton(430,0.0, 0.254787,3.69113e-09,13.6,4.77487,13.7,70.);
    dedt_Proton(433,0.0, -3.0,1.6113,2.4,0.277499,2.5,70.);
    dedt_Proton(434,0.0, -1.75476,1.34208,-0.792445,0.101983,2,70.);
  }
  if(runno<62605) //x2: was 62718)
    dedt_Proton(530,0.0, -0.168344,0.0432862,-1.21441,1,2,70.);
  dedt_Proton(331,0.0, -0.1302998,0.00628346,-4.78109,1.75404,2,70.); //x2:0.0,0.0302998
  if(runno<62605) //x2: was 62718
    dedt_Proton(133,0.0, 0.338925,9.58889e-10,2.4,4.64712,2.5,70.);
  else if(runno<62718)
    dedt_Proton(233,0.0, -0.4,0.0);
  dedt_Proton(533,0.0, -0.0536817,0.000965426,0.8,1.68532,2,70.);
  if(runno>62500 && runno<62718) //?? 62607)
    dedt_Proton(134,0.0, 2.3,0.0);
  dedt_Proton(334,0.0, 2.74176,-1.99182,4.3469,0.116272,4.4469,70.);
  dedt_Proton(335,0.0, -0.593243,0.225001,-2.4,0.391838,2,70.);
  if(runno>62604 && runno<62718)
    dedt_Proton(236,0.0,-2.89233,2.77,0.378457,0.200822,2,70.);
  else
    dedt_Proton(236,0.0,-2.89233,2.26755,0.378457,0.200822,2,70.);
  dedt_Proton(337,0.0, -3,1.95507,4,0.204048,4.1,70.);
  dedt_Proton(538,0.0, -0.442077,0.0256601,0.8,1.17523,2,70.);
  if(runno>62604 && runno<62718) {
    dedt_Proton(239,0.0, 3.,-2.17317,6.36194,0.254568,6.46194,70.); //+pi
    dedt_Proton(439,0.0, -0.7,0.0);
  }
  else {
    dedt_Proton(239,0.0, 3.,-2.17317,6.36194,0.254568,6.46194,70.); //+pi
    dedt_Proton(439,0.0, 2.26499,-1.9095,1.6,0.1,2,70.);
  }
  dedt_Proton(539,0.0, -0.199483,0.0144536,1.6,1.07901,2,70.);
  if(runno<62604 || runno>62717)
    dedt_Proton(142,0.0, -0.165495,0.0145086,-1.64358,0,2,70.);
  else
    dedt_Proton(440,0.0, -1.4,0.0);
  dedt_Proton(442,0.0, -0.2203,0.011732,1.59851,0,2,70.);
  dedt_Proton(542,0.0, -0.167378,0.0207299,6.87864,0,2,70.);
  if(runno<62500 || runno>62604)
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
  if(runno>62604 && runno<62718) {
    dedt_Pion(323,0.0, -0.6,0.);
    dedt_Pion(523,0.0, -1.7,0.0); 
  }
  dedt_Pion(525,0.0, -0.779324,0.5,3.2,0.238471,3.3,70.);
  dedt_Pion(127,0.0, -0.277002,0.0187147,0.8,1.24111,2,70.);
  if(runno<62605) dedt_Pion(228,0.0, -0.102935,0.0157658,0.79733,0,2,70.);
  dedt_Pion(528,0.0, -0.830921,0.475982,-0.601109,0.232823,2,70.);
  if(runno>62604) 
    dedt_Pion(229,0.0, -0.45554,0.5,3.2,0.277431,3.3,70.);//<x3:-0.85554
  dedt_Pion(529,0.0, -0.023221,0.00184413,3.2,1.38569,3.3,70.);
  dedt_Pion(629,0.0, -0.376183,0.0186985,0.784218,0,2,70.);
  if(runno>62604) dedt_Pion(230,0.0, -0.10738,0.0241076,2.38848,0,2,70.);
  if(runno<62718) 
    dedt_Pion(530,0.0, 3.0,-2.24322,-0.16768,0.157661,2,70.);
  else if( runno>62717)
    dedt_Pion(530,0.0, 3.0,-2.13406,0.8,0.133982,2,70.);
  dedt_Pion(630,0.0, -2.77573,2.34239,6.44143,0.114437,6.54143,70.);
  dedt_Pion(331,0.0, 0.0302998,0.00628346,-4.78109,1.75404,2,70.);
  if(runno<62605) //x2: was 62718
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
  else
    dedt_Pion(239,0.0, -2.2,0.0);
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
  int ndat=-1;
  while(SCIDarray[++ndat] > 0) ndat++;
  double myTimeOffsets[ndat];
  for(int i=0; i<ndat; i++)
    myTimeOffsets[i] = TimeOffsets[i];
  SCIDcut_arrayInit(SCIDarray, myTimeOffsets);
}

void TSCtimeCorr::SCIDcut_arrayInit(int *SCIDarray, double *TimeOffsets) {

  int i=-1;
  if(TimeOffsets) {
    while(SCIDarray[++i] > 0) SCid_Tcorr[SCpaddleId(SCIDarray[i])] = SC_COUNTER_OFF;
  }
  else {
    while(SCIDarray[++i] > 0) SCid_Tcorr[SCpaddleId(SCIDarray[i])] = TimeOffsets[i];
  }
}

//-------------- PrintInfo(runno,flag) ----------------------------------------------------
// flag=1:  print list of SC paddles that are turned off and that are time corrected,
// flag=2:  print SCPID beta cuts used for corrected particle id.
// flag=3:  print both

void TSCtimeCorr::PrintInfo(int runno, int flag) {
  
  if(flag != 2) {
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
	std::cout.precision(3);
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
    if(target_center!=0.0) {
      std::cout.precision(3);
      std::cout<<" (target_center at z="<<target_center<<" cm)"<<std::endl;
    }
  }
  //
  if(flag > 1) {
    const char *clst[SCPID_ntypes]={"Electron","Pion","Kaon","Proton","Deuteron","Helium3 "};
    
    std::cout<<"\nBeta cuts for corrected particle ID (tested in the following order):"<<std::endl;
    std::cout<<"(type     PDGmass     mass_min off_min  mass_max off_max   p_min  p_max)"<<std::endl;
    for(int i=0; i<SCPID_ntypes; i++) {
      int it = PIDsearchorder[i];
      std::cout.width(8); 
      std::cout<<std::left<<clst[it]<<" (";
      std::cout.precision(5); 
      std::cout<<sqrt(TMath::Abs(PIDmass2[it]))<<" GeV)\t";
      std::cout.precision(2); std::cout.width(6);
      std::cout<<sqrt(beta_mass2_min[it])<<"  "<<beta_off_min[it]<<"\t";
      std::cout.precision(2); std::cout.width(6);
      std::cout<<TMath::Sign(sqrt(TMath::Abs(beta_mass2_max[it])),beta_mass2_max[it])<<"  "<<beta_off_max[it]<<"\t   ";
      std::cout.precision(2); std::cout.width(6);
      std::cout<<beta_p_min[it]<<"   "<<beta_p_max[it]<<std::endl;
    }
  }
  std::cout<<"(min. probability to assign particle ID in PIDcorr routine: "<<PIDminProb<<")"<<std::endl;
  std::cout<<std::endl;
}

//---------------------------------------------------------------------------

#ifdef __CINT__ 
#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 
#pragma link C++ class TSCtimeCorr+; 
#endif 

#endif

//end: ifndef __INCLUDE_ONLY__
#endif
