//
//  detector and particle classes used in rb_gtest.C:
//  class TDetectorHit: some detector info for ST, SC, EC
//  clas  TMyPart:  tracking and timing info copied from banks
//

#ifndef __TMYPART__
#define __TMYPART__

#include <iostream>
#include <vector>
#include <string>
#include <TObject.h>
#include "bankvars.h"
#include "TSCtimeCorr.h"

enum detector_type {det_st, det_dc, det_cc, det_sc, det_ec};

#ifndef ClightD
#define ClightD
static const double Clight = 29.97925;         // (cm/ns)
#endif

//************************************************************************************
//********************          class  TDetectorHit             **********************
//************************************************************************************
class TDetectorHit {
 private:
  enum { bSTRE, bSTR, bCC1, bSCRC, bSCR, bECHB, bSTPB, bCCPB, bSCPB, bECPB, bGPID, bPART, bEVNT, bNbanks};
  int ibank;          //100*record# + bank
  int irow;           //row-1 (from 0...SCRC_NH-1 etc) 
  int idet;           //det_type + 10*bank
  int itrk;           //TBER,TBTR track number (1...TBTR_NH or TBER_NH)
 public:
  int    Index;          //st,sc: 100*sector+paddle_id; ec: 100*sector+paddle_u
  int    Status;
  double Trec;           //time from reconstruction
  double Time;           //corrected time
  double Path;           //pathlength from track 'vertex'
  double Edep;           //deposited energy
  double PosDet;          //position on paddle (detector coord)
  TVector3 Position;      // position from tracking (extrapolation) (in sector coord)

  TDetectorHit();
  ~TDetectorHit() {};
  int setDetectorHit(int det_type, GPID_t gpid, double Tcorr=0.0);
  int setDetectorHit(int det_type, PART_t part, double Tcorr=0.0);
  int setDetectorHit(int det_type, EVNT_t evnt, double Tcorr=0.0);
  int setDetectorHit(int det_type, int ind, double trec, double path, double edep=0.0);
  int getTrackPlane(int det_type, int trkno);     //output: record#*100+tdpl_row; trkno=1,...,ntrk
  double CalculateBeta( TDetectorHit* det );
  int getDetectorType(char *cbank=NULL);
  int getDetectorBank(char *cbank, int *irec);    //output: row, bank name and recordno
#ifdef __CINT__
  ClassDef(TDetectorHit,1);  //for integration into ROOT: has to be processed using rootcint
#endif
};
// (methods at the end of this file)

//************************************************************************************
//********************          class  TMyPart                ************************
//************************************************************************************

class TMyPart {
 private:
  int sortTimeInfo(TSCtimeCorr *TSCcorr, int itg1, int itg2, double myProb2[SCPID_ntypes], int flag=0);
  int irow;   //particle info from row+100*bank (bank=1:GPID, bank=2:PART, bank=3:EVNT, irow=0: no bank) 
 public:
  TLorentzVector p4rec;
  TLorentzVector p4cor;
  TVector3 vert;
  TVector3 dir;
  double mom;
  double phisec;
  double betam;
  double beta;
  double beta2;
  int charge;
  int sector;
  int gid;        //geant ID from particle bank
  int myid;
  int myid2;
  int itype;    //always >=0 after initialization
  int itype2;   //always >=0 after initialization
  double PIDprob[SCPID_ntypes+1];
  double Tpho[SCPID_ntypes+1];
  int itagr[SCPID_ntypes+1];
  double beta_mtrk;
  int itagr_mtrk;
  TDetectorHit sc;
  TDetectorHit ec;
  TDetectorHit st;
  TMyPart(TVector3 mom, TVector3 vertex, double beta, int q);
  TMyPart(GPID_t gpid);
  TMyPart(PART_t part);
  TMyPart(EVNT_t evnt);
 TMyPart(void) :
  irow(0), mom(0.0), phisec(0.0), betam(0.0), beta(0.0), beta2(0.0), charge(0), sector(0), gid(0), myid(-1), myid2(-1), itype(-1), itype2(-1), beta_mtrk(0.0), itagr_mtrk(-1) { };
  ~TMyPart() {};
  int setTimeInfo(GPID_t gpid, TSCtimeCorr *TSCcorr=NULL, int itag=-1);
  int setTimeInfo(PART_t part, TSCtimeCorr *TSCcorr=NULL, int itag=-1);
  int setTimeInfo(EVNT_t evnt, TSCtimeCorr *TSCcorr=NULL, int itgpb=-1);
  TVector3* Intersect(TMyPart part2, TVector3 *dirIntersect=NULL, TVector3 *p1=NULL, TVector3 *p2=NULL);
  TVector3* Intersect(TVector3 vert2, TVector3 dir2, TVector3 *dirInterset=NULL, TVector3 *p1=NULL, TVector3 *p2=NULL);
  double getDOCA(TMyPart mp2, TVector3 *p1=NULL, TVector3 *p2=NULL);
  double getDOCA(TVector3 vert2, TVector3 dir2, TVector3 *p1=NULL, TVector3 *p2=NULL);
  double CalculateT0 (TDetectorHit* detA=NULL, TDetectorHit *detB=NULL);
  double CalculateBeta(TMyPart *mp2, TDetectorHit* det=NULL);
  double CalculateBeta(double Ttag, double Zoffset=0.0, TSCtimeCorr *TSCcorr=NULL, TDetectorHit* det=NULL);
  double CalculateMassSqr(double mybeta=0.0);
  double CalculateMass(double mybeta=0.0);
  int getCommonITAGR(int mypartid, TMyPart mpart2, int mpart2id=nMYPIDcodes);
  int getCommonITAGR(int mypartid, TMyPart &mpart2out, std::vector<TMyPart> &mpart2list, int mpart2id=nMYPIDcodes);
  friend ostream& operator<< (ostream&, const TMyPart);
#ifdef __CINT__
  ClassDef(TMyPart,1);  //for integration into ROOT: has to be processed using rootcint
#endif
};

//***************** TMyPart  constructor  *****************

TMyPart::TMyPart(TVector3 momentum, TVector3 vertex, double beta_meas, int q) {
  irow=0;
  charge = q;
  beta = betam = beta_meas;
  vert = vertex;
  mom = momentum.Mag();
  dir = momentum.Unit();
  if(betam < 0.05 || betam == 1.0) {
    if(charge) p4rec.SetVectM(momentum,0.000511);
    else       p4rec.SetVectM(momentum,0.0);
  }
  else {
    double mass2 = mom*mom*(1/betam/betam -1.0);
    p4rec.SetVectM(momentum,sqrt(TMath::Abs(mass2)));
  }
  double phirad = momentum.Phi()*TMath::RadToDeg() + 30.;
  if(phirad<0) phirad+=360.;
  sector = (int)phirad/60. +1;
  phisec = phirad - sector*60 + 30.;

  myid=myid2=-1;
  itype=itype2=SCPID_ntypes;
  beta_mtrk = beta2 = 0.0;
  for(int i=0; i<SCPID_ntypes; i++) {
    PIDprob[i] = 0.0;
    itagr[i] = -1;
    Tpho[i]  = -500.;
  }
  itagr_mtrk = -1;
}

//-------- GPID

TMyPart::TMyPart(GPID_t gpid) {
  irow = 100+gpid.trkid;
  charge = gpid.q;
  vert.SetXYZ(gpid.x,gpid.y,gpid.z);
  mom = sqrt(gpid.px*gpid.px + gpid.py*gpid.py + gpid.pz*gpid.pz);
  if(mom < 1.E-6)
    dir.SetXYZ(gpid.px, gpid.py, gpid.pz);
  else
    dir.SetXYZ(gpid.px/mom, gpid.py/mom, gpid.pz/mom);
  beta = betam = gpid.betam;
  if(betam < 0.05 || betam == 1.0) {
    if(charge) p4rec.SetXYZM(gpid.px,gpid.py,gpid.pz,0.000511);
    else       p4rec.SetXYZM(gpid.px,gpid.py,gpid.pz,0.0);
  }
  else {
    double mass2 = mom*mom*(1/betam/betam -1.0);
    p4rec.SetVectM( (dir*mom), sqrt(TMath::Abs(mass2)));
  }
  double phirad = dir.Phi()*TMath::RadToDeg() + 30.;  
  if(phirad<0) phirad+=360.;
  sector = (int)phirad/60. +1;
  phisec = phirad - sector*60 + 30.;
  gid = (beta<0.05)? 0 : gpid.pid;
  myid=myid2=-1;
  itype=itype2=SCPID_ntypes;
  beta_mtrk = beta2 = 0.0;
  for(int i=0; i<=SCPID_ntypes; i++) {
    PIDprob[i] = 0.0;
    itagr[i] = -1;
    Tpho[i]  = -500.;
  }
  itagr_mtrk = -1;
  sc.setDetectorHit(det_sc, gpid);
  st.setDetectorHit(det_st, gpid);
}

//------- PART

TMyPart::TMyPart(PART_t part) {
  irow = 200+part.trkid;
  charge = part.q;
  vert.SetXYZ(part.x,part.y,part.z);
  mom = sqrt(part.px*part.px + part.py*part.py + part.pz*part.pz);
  if(mom < 1.E-6)
    dir.SetXYZ(part.px, part.py, part.pz);
  else
    dir.SetXYZ(part.px/mom, part.py/mom, part.pz/mom);
  betam = 0.0;
  if(part.trkid>0 && part.trkid<=TBID_NH) betam = TBID[part.trkid-1].beta;
  if(betam < 0.05 || betam == 1.0) {
    if(charge) p4rec.SetXYZM(part.px,part.py,part.pz,0.000511);
    else       p4rec.SetXYZM(part.px,part.py,part.pz,0.0);
  }
  else {
    double mass2 = mom*mom*(1/betam/betam -1.0);
    p4rec.SetXYZM(part.px,part.py,part.pz,sqrt(TMath::Abs(mass2)));
  }
  beta = betam;
  double phirad = dir.Phi()*TMath::RadToDeg() + 30.;
  if(phirad<0) phirad+=360.;
  sector = (int)phirad/60. +1;
  phisec = phirad - sector*60 + 30.;
  gid = (beta<0.05)? 0 : part.pid;
  myid=myid2=-1;
  itype=itype2=SCPID_ntypes;
  beta_mtrk = beta2 = 0.0;
  for(int i=0; i<=SCPID_ntypes; i++) {
    PIDprob[i] = 0.0;
    itagr[i] = -1;
    Tpho[i]  = -500.;
  }
  itagr_mtrk = -1;
  sc.setDetectorHit(det_sc, part);
  st.setDetectorHit(det_st, part);
}

//-------- EVNT

TMyPart::TMyPart(EVNT_t evnt) {
  irow=300;
  charge = evnt.Charge;
  vert.SetXYZ(evnt.X,evnt.Y,evnt.Z);
  mom = evnt.Pmom;
  dir.SetXYZ(evnt.cx, evnt.cy, evnt.cz);
  beta = betam = evnt.Beta;
  if(betam < 0.05 || betam == 1.0) {
    if(charge) p4rec.SetVectM( (dir*mom), 0.000511);
    else       p4rec.SetVectM( (dir*mom), 0.0);
  }
  else {
    p4rec.SetVectM( (dir*mom), evnt.Mass);
  }
  double phirad = dir.Phi()*TMath::RadToDeg() + 30.;
  if(phirad<0) phirad+=360.;
  sector = (int)phirad/60. +1;
  phisec = phirad - sector*60 + 30.;
  gid=0;
  myid=myid2=-1;
  itype=itype2=SCPID_ntypes;
  beta_mtrk = beta2 = 0.0;
  for(int i=0; i<=SCPID_ntypes; i++) {
    PIDprob[i] = 0.0;
    itagr[i] = -1;
    Tpho[i]  = -500.;
  }
  itagr_mtrk = -1;
  sc.setDetectorHit(det_sc, evnt);
  st.setDetectorHit(det_st, evnt);
}

//*****************  CalculateT0 (using 1 or 2 detector times)  ***********

double TMyPart::CalculateT0( TDetectorHit* detA, TDetectorHit* detB) {

  TDetectorHit *mydet=detA;
  if(!mydet) {
    if(sc.Index>0)      mydet=&sc;
    else if(ec.Index>0) mydet=&ec;
    else if(st.Index>0) mydet=&st;
    if(!mydet) return -500.;
  }
  double myBeta = (detB)? mydet->CalculateBeta(detB) : beta;
  if(myBeta > 0.0) 
    return mydet->Time - mydet->Path / myBeta / Clight; 
  else
    return -500.0;
}

//****************  CalculateBeta (using 2 tracks)  *************

double TMyPart::CalculateBeta( TMyPart* mp2, TDetectorHit* det ) {

  TDetectorHit *mydet=det;
  if(!mydet) {
    if(sc.Index>0)      mydet=&sc;
    else if(ec.Index>0) mydet=&ec;
    else if(st.Index>0) mydet=&st;
    if(!mydet) return -1.;
  }
  return  mydet->Path / (mydet->Time - mp2->CalculateT0()) / Clight; 
}

//****************  CalculateBeta (using Tpho and detector hit)  *************

double TMyPart::CalculateBeta( double Ttag, double Zoffset, TSCtimeCorr *TSCcorr, TDetectorHit* det ) {

  TDetectorHit *mydet=det;
  if(!mydet) {
    if(sc.Index>0)      mydet=&sc;
    else if(ec.Index>0) mydet=&ec;
    else if(st.Index>0) mydet=&st;
    if(!mydet) return beta;
  }
  if(TSCcorr && mydet->getDetectorType()==det_sc) {
    return TSCcorr->BetaCorr(sc.Time, sc.Path, Ttag, vert.Z(), sc.Index, 0, mom, sc.Edep, st.Edep);
  }
  else {
    return mydet->Path / (mydet->Time - Ttag - (vert.Z() - Zoffset)/Clight) / Clight; 
  }
}

//**************   CalculateMass  *****************

double TMyPart::CalculateMassSqr( double mybeta ) {
  if(mybeta < 0.05) mybeta = beta;
  if(fabs(mybeta-1.0) < 1.E-5) 
    return 2.6E-7;   //M_Electron^2
  else
    return (mom*mom*(1./mybeta/mybeta -1.));
}

double TMyPart::CalculateMass( double mybeta ) {
  double mymass2 = CalculateMassSqr(mybeta);
  return (mymass2<0)? -sqrt(-mymass2) : sqrt(mymass2);
}

//***********************  print  *****************

ostream& operator<< (ostream& os, const TMyPart mp) {
  os.width(5);
  os << mp.myid;
  os.width(3);
  os << mp.charge;
  os.flags(std::ios::fixed);
  os.precision(4);
  os.width(8);
  os << mp.p4rec.X();
  os.width(8);
  os << mp.p4rec.Y();
  os.width(8);
  os << mp.p4rec.Z();
  os.width(8);
  os << mp.p4rec.E();
  return os;
}

//*******************  setTimeInfo  ***************************
//----------- GPID

int TMyPart::setTimeInfo(GPID_t gpid, TSCtimeCorr *TSCcorr, int itag) {

  if(!TSCcorr) return -2;
  sc.Time = TSCcorr->TimeCorr(gpid);
  beta    = TSCcorr->BetaCorr(gpid, &itag);
  double betast = betam;
  if((beta<0.05 || beta>1.15) && sc.Index>100 && st.Index>100) {
    betast = (sc.Path - st.Path)/(sc.Time - st.Time)/Clight;
    if(betam>0.05 && betam<1.1 && fabs(betast-betam)<fabs(beta-betam)) {
      beta = betast;
      itag = gpid.tagrid -1;
    }
  }
  int mygid = TSCcorr->PIDcorr(gpid, PIDprob, beta, sc.Edep, st.Edep);
  if(!mygid) mygid = TSCcorr->PIDcorr(gpid, PIDprob, betast, sc.Edep, st.Edep);
  myid = TSCcorr->PIDconvert(mygid, SCPIDin_GID+10*SCPIDin_myid);
  double myPIDprob[SCPID_ntypes];
  int itagX = -1;
  if(TAGR_NH > 1) {
    itagX = TAGR_NH;
    if(itag>=0) TAGR[itag].TPHO -=500;
    double mybeta  = TSCcorr->BetaCorr(gpid, &itagX);
    mygid = TSCcorr->PIDcorr(gpid, myPIDprob, mybeta, sc.Edep, st.Edep);
    myid2 = TSCcorr->PIDconvert(mygid, SCPIDin_GID+10*SCPIDin_myid);
    if(itag>=0) TAGR[itag].TPHO +=500;
  }
  return sortTimeInfo(TSCcorr, itag, itagX, myPIDprob);
}

//------------ PART

int TMyPart::setTimeInfo(PART_t part, TSCtimeCorr *TSCcorr, int itag) {

  if(!TSCcorr) return -2;
  sc.Time = TSCcorr->TimeCorr(part);
  beta  = TSCcorr->BetaCorr(part, &itag);
  double betast=betam;
  if((beta<0.05 || beta>1.15) && sc.Index>100 && st.Index>100) {
    betast = (sc.Path - st.Path)/(sc.Time - st.Time)/Clight;
    if(betam>0.05 && betam<1.1 && fabs(betast-betam)<fabs(beta-betam)) {
      beta = betast;
    }
  }
  if(itag < 0 && beta != betast) return -1;

  int mygid = TSCcorr->PIDcorr(part, PIDprob, beta, sc.Edep, st.Edep);
  if(!mygid) mygid = TSCcorr->PIDcorr(part, PIDprob, betast, sc.Edep, st.Edep);
  myid = TSCcorr->PIDconvert(mygid, SCPIDin_GID+10*SCPIDin_myid);
  double myPIDprob[SCPID_ntypes];
  int itagX = -1;
  if(TAGR_NH > 1) {
    itagX = TAGR_NH;
    if(itag>=0) TAGR[itag].TPHO -=500;
    double mybeta  = TSCcorr->BetaCorr(part, &itagX);
    mygid = TSCcorr->PIDcorr(part, myPIDprob, mybeta, sc.Edep, st.Edep);
    myid2 = TSCcorr->PIDconvert(mygid, SCPIDin_GID+10*SCPIDin_myid);
    if(itag>=0) TAGR[itag].TPHO +=500;
  }
  return sortTimeInfo(TSCcorr, itag, itagX, myPIDprob);
}

//----------- EVNT

int TMyPart::setTimeInfo(EVNT_t evnt, TSCtimeCorr *TSCcorr, int itgpb) {

  if(!TSCcorr) return -2;
  sc.Time = TSCcorr->TimeCorr(evnt);
  beta  = TSCcorr->BetaCorr(evnt, &itgpb);
  double betast=betam;
  if((beta<0.05 || beta>1.15) && sc.Index>100 && st.Index>100) {
    betast = (sc.Path - st.Path)/(sc.Time - st.Time)/Clight;
    if(betam>0.1 && betam<1.1 && fabs(betast-betam)<fabs(beta-betam)) {
      beta = betast;
    }
  }
  if(itgpb < 0 && beta != betast) return -1;

  int mygid = TSCcorr->PIDcorr(evnt, PIDprob, beta, sc.Edep, st.Edep);
  if(!mygid) mygid = TSCcorr->PIDcorr(evnt, PIDprob, betast, sc.Edep, st.Edep);
  myid = TSCcorr->PIDconvert(mygid, SCPIDin_MCID+10*SCPIDin_myid);
  double myPIDprob[SCPID_ntypes];
  int itgpb2 = -1;
  if(TGPB_NH > 1) {
    itgpb2 = TGPB_NH;
    if(itgpb>=0) TGPB[itgpb].Time -=500;
    double mybeta  = TSCcorr->BetaCorr(evnt, &itgpb2);
    mygid = TSCcorr->PIDcorr(evnt, myPIDprob, mybeta, sc.Edep, st.Edep);
    myid2 = TSCcorr->PIDconvert(mygid, SCPIDin_MCID+10*SCPIDin_myid);
    if(itgpb>=0) TGPB[itgpb].Time +=500;
  }
  return sortTimeInfo(TSCcorr, itgpb, itgpb2, myPIDprob, 1);
}

//***********************  sortTimeInfo  *****************

// flag=0: get info from TAGR bank (default);  flag=1: get info from TGPB bank
int TMyPart::sortTimeInfo(TSCtimeCorr *TSCcorr, int itg1, int itg2, double myProb[SCPID_ntypes], int flag) {

  int ntgmax = (TAGR_NH)? TAGR_NH : TGPB_NH;
  int *ptr_tagr = new int[ntgmax];
  double *myTime = new double[ntgmax];
  if(!flag) {
    for(int i=0; i<TAGR_NH; i++) {
      ptr_tagr[i] = i;
      myTime[i] = TAGR[i].TPHO;
    }
  }
  else if(!TAGR_NH) {
    for(int i=0; i<TGPB_NH; i++) {
      ptr_tagr[i] = i+100;
      myTime[i] = TGPB[i].Time;
    }
  }
  else {
    for(int i=0; i<TGPB_NH; i++) {
      int it = TMath::Abs(TGPB[i].pointer)/1000;
      ptr_tagr[i] = it;
      myTime[i] = TAGR[it].TPHO;
    }
  }
  int convflag = -1;
  if(charge>0)
    convflag = SCPIDin_postype+10*SCPIDin_myid;
  else if(charge<0)
    convflag = SCPIDin_negtype+10*SCPIDin_myid;
  int j1 = TMath::LocMax(SCPID_ntypes, PIDprob);
  int j2 = -1;
  double maxProb = 0.;
  for(int i=0; i<SCPID_ntypes; i++) {
    if(itg1 >= 0) {
      itagr[i] = ptr_tagr[itg1];
      Tpho[i]  = myTime[itg1];
    }
    else 
      Tpho[i]  = CalculateT0();
    if((i!=j1) && (PIDprob[i] > maxProb)) {
      maxProb = PIDprob[i];
      j2 = i;
    }
  }
  beta2 = beta;
  int k1=-1, i1=-1;
  if(itg2 >= 0) {
    k1 = TMath::LocMax(SCPID_ntypes, myProb);
    if( (j2 < 0) || (myProb[k1] > maxProb) ) i1 = k1+100;
  }
  if(i1 < 0) {
    if(maxProb>0.85)
      myid2 = TSCcorr->PIDconvert(j2, convflag);
    itype  = j1;
    if(j2>=0) itype2 = j2;
    if(itype<0 || itype2<0) {
      std::cerr<<"TMyPart::sortTimeInfo(1): itype="<<itype<<" itype2="<<itype2<<"  j1="<<j1<<" j2="<<j2<<" k1="<<k1<<" i1="<<i1<<std::endl;
      if(itype<0) itype=SCPID_ntypes;
      if(itype2<0) itype2=SCPID_ntypes;
    }
    for(int i=0; itg2>=0 && i<SCPID_ntypes; i++) {
      if(i!=j1 && i!=j2) {
	if(myProb[i] > PIDprob[i]) {
	  PIDprob[i]= myProb[i];
	  itagr[i]  = ptr_tagr[itg2];
	  Tpho[i]   = myTime[itg2];
	}
      }
    }
    delete myTime;
    delete ptr_tagr;
    return j1;
  }

  double betak = TSCcorr->BetaCorr(sc.Time, sc.Path, myTime[itg2], vert.Z(), sc.Index, 0, mom, sc.Edep, st.Edep);
  int k2 = -1;
  maxProb = 0.;
  for(int i=0; i<SCPID_ntypes; i++) {
    if((i!=k1) && (myProb[i] > maxProb)) {
      maxProb = myProb[i];
      k2 = i;
    }
  }
  int i2 = (PIDprob[j1]>=maxProb || k2<0)? j1 : k2+100;
  if(PIDprob[j1] >= myProb[k1]) {
    i1 = j1; 
    i2 = k1+100;
  }
  if(j1 == k1) {
    itype = j1;
    if(i1 == j1) {
      PIDprob[SCPID_ntypes]= myProb[k1];
      itagr[SCPID_ntypes]  = ptr_tagr[itg2];
      Tpho[SCPID_ntypes]   = myTime[itg2];
      beta2 = betak;
    }
    else {
      PIDprob[SCPID_ntypes]= PIDprob[j1];
      itagr[SCPID_ntypes]  = itagr[j1];
      Tpho[SCPID_ntypes]   = Tpho[j1];
      beta = betak;
      if(i2==k2+100) {
	itype2 = k2;
	beta2 = betak;
      }
    }
  }
  else {
    itype  = (i1%100);
    if(i1>99) beta = betak;
    itype2 = (i2%100); 
    if(i2>99) beta2 = betak;
  }
  if(i2==j1) 
    myid2 = myid;
  else 
    myid2 = TSCcorr->PIDconvert((i2%100), convflag);
  if(i1>99 && convflag>0) 
    myid = TSCcorr->PIDconvert((i1%100), convflag);
    if(itype<0 || itype2<0) {
      std::cerr<<"TMyPart::sortTimeInfo(2): itype="<<itype<<" itype2="<<itype2<<"  j1="<<j1<<" j2="<<j2<<" k1="<<k1<<" k2="<<k2<<" i1="<<i1<<" i2="<<i2<<std::endl;
      if(itype<0) itype=SCPID_ntypes;
      if(itype2<0) itype2=SCPID_ntypes;
    }
  for(int i=0; itg2>=0 && i<SCPID_ntypes; i++) {
    if(myProb[i] > PIDprob[i]) {
      PIDprob[i]= myProb[i];
      itagr[i]  = ptr_tagr[itg2];
      Tpho[i]   = myTime[itg2];
    }
  }
  delete myTime;
  delete ptr_tagr;
  return (i1%100);
}

//***********************  Intersect  *****************

TVector3* TMyPart::Intersect(TMyPart mp2, TVector3 *dirIntersect, TVector3 *p1, TVector3 *p2) {
  return Intersect(mp2.vert, mp2.dir, dirIntersect, p1, p2);
}

TVector3* TMyPart::Intersect(TVector3 vert2, TVector3 dir2, TVector3 *dirIntersect, TVector3 *p1, TVector3 *p2) {
  TVector3 *myP1=NULL, *myP2=NULL;
  if(p1) myP1 = p1;
  if(p2) myP2 = p2;
  getDOCA( vert2, dir2, myP1, myP2 );
  if(dirIntersect) {
    TVector3 temp = (dir + dir2) * 0.5;
    dirIntersect->SetXYZ(temp.X(), temp.Y(), temp.Z());
  }
  static TVector3 mP12((myP1->X()+myP2->X())*0.5, (myP1->Y()+myP2->Y())*0.5, (myP1->Z()+myP2->Z())*0.5); 
  //  TVector3 *myP12=&mP12;
  //  myP12->SetXYZ((myP1->X()+myP2->X())*0.5, (myP1->Y()+myP2->Y())*0.5, (myP1->Z()+myP2->Z())*0.5); 
  return &mP12;
}

double TMyPart::getDOCA( TMyPart mp2, TVector3 *p1, TVector3 *p2 ) {
  return getDOCA(mp2.vert, mp2.dir, p1, p2);
}

double TMyPart::getDOCA( TVector3 vert2, TVector3 dir2, TVector3 *p1, TVector3 *p2 ) {
  double cosPhi = dir*dir2;
  double denom = 1.0 - cosPhi * cosPhi;
  double t1 = 0.0;
  double t2;
  if (fabs (denom) < 1.E-10) {  // parallel? then no correction to vertex
    t2 = (vert2 - vert) * dir / cosPhi;
  } 
  else {
    double A = (vert2 - vert) * dir;
    double B = (vert2 - vert) * dir2;
    t1 = (A - cosPhi * B) / denom;
    t2 = (cosPhi * A - B) / denom;
  }

  if(p1 && p2) {
    *p1 = vert + t1 * dir; 
    *p2 = vert2 + t2 * dir2;
  }
  TVector3 pdiff(vert - vert2 + t1 * dir - t2 * dir2);
  return pdiff.Mag();
}

//----------- getCommonITAGR  ---------------------------------------------------------------
// check for same TAGR entry (row = 0,...,NTAGR-1) for two particles 
//  (myid<0: any particle (also those without part.id); myid=nMYPIDcodes: any particle with part.id)

int TMyPart::getCommonITAGR(int mypartid, TMyPart mpart2, int mpart2id) {
  int itag=-1;
  int itg1=0, itg2=0;
  if(mypartid<0) {
    if(itype >= 0)  itg1 |= BIT(itagr[itype]);
    if(itype2 >= 0) itg1 |= BIT(itagr[itype2]);
  }
  else if(mypartid >= nMYPIDcodes) {
    if(myid >= 0)  itg1 |= BIT(itagr[itype]);
    if(myid2 >= 0) itg1 |= BIT(itagr[itype2]);
  }
  else {
    if(myid == mypartid)  itg1 |= BIT(itagr[itype]);
    if(myid2 == mypartid) itg1 |= BIT(itagr[itype2]);
  }
  if(!itg1) return -1;

  if(mpart2id<0) {
    if(mpart2.itagr[mpart2.itype] >= 0)  itg2 |= BIT(mpart2.itagr[mpart2.itype]);
    if(mpart2.itagr[mpart2.itype2] >= 0) itg2 |= BIT(mpart2.itagr[mpart2.itype2]);
  }
  else if(mpart2id >= nMYPIDcodes) {
    if(mpart2.myid >= 0)  itg2 |= BIT(mpart2.itagr[mpart2.itype]);
    if(mpart2.myid2 >= 0) itg2 |= BIT(mpart2.itagr[mpart2.itype2]);
  }
  else {
    if(mpart2.myid == mpart2id)  itg2 |= BIT(mpart2.itagr[mpart2.itype]);
    if(mpart2.myid2 == mpart2id) itg2 |= BIT(mpart2.itagr[mpart2.itype2]);
  }
  int itg12 = (itg1&itg2);
  if(itg12) {
    int i=0;
    while( !(itg12&(1<<i)) ) i++;
    int j=30;
    while( !(itg12&(1<<j)) ) j--;
    if(i==j) itag = i;
    else     itag = i + 1000*j;
  }
  return itag;
}

//---------------------------------------------------------------------------

int TMyPart::getCommonITAGR(int mypartid, TMyPart &mpart2out, std::vector<TMyPart> &mpart2list, int mpart2id) {

  if(mpart2list.size() < 1) return -1;
  int itag=-1;
  int itg1=0;
  if(mypartid<0) {
    if(itype >= 0)  itg1 |= BIT(itagr[itype]);
    if(itype2 >= 0) itg1 |= BIT(itagr[itype2]);
  }
  else if(mypartid >= nMYPIDcodes) {
    if(myid >= 0)  itg1 |= BIT(itagr[itype]);
    if(myid2 >= 0) itg1 |= BIT(itagr[itype2]);
  }
  else {
    if(myid == mypartid)  itg1 |= BIT(itagr[itype]);
    if(myid2 == mypartid) itg1 |= BIT(itagr[itype2]);
  }
  if(!itg1) return -1;

  for(std::vector<TMyPart>::iterator jj=mpart2list.begin(); jj!=mpart2list.end() && itag<0; jj++) {
    int itg2 = 0;
    if(mpart2id < 0) {
      if((*jj).itagr[(*jj).itype] >= 0)  itg2 |= BIT((*jj).itagr[(*jj).itype]);
      if((*jj).itagr[(*jj).itype2] >= 0) itg2 |= BIT((*jj).itagr[(*jj).itype2]);
    }
    else if(mpart2id == nMYPIDcodes) {
      if((*jj).myid >= 0)  itg2 |= BIT((*jj).itagr[(*jj).itype]);
      if((*jj).myid2 >= 0) itg2 |= BIT((*jj).itagr[(*jj).itype2]);
    }
    else {
      if((*jj).myid == mpart2id)  itg2 |= BIT((*jj).itagr[(*jj).itype]);
      if((*jj).myid2 == mpart2id) itg2 |= BIT((*jj).itagr[(*jj).itype2]);
    }
    int itg12 = (itg1&itg2);
    if(itg12) {
      int i=0;
      while( !(itg12&(1<<i)) ) i++;
      itag = i;
      mpart2out = (*jj);
    }
  }
  return itag;
}


//********************************************************************************
//********************************************************************************
//---------------   TDetectorHit  constructor    ---------------------------------

TDetectorHit::TDetectorHit() : ibank(-1), irow(-1), idet(-1), itrk(0), Index(-1), Trec(0), Time(0), Path(0), Edep(0), PosDet(-1000.) {}

//---------------------------------------------------------------------------
int TDetectorHit::getDetectorType(char *cbank) {
  static const char *cbanks[]={"STRE","STR ","CC1 ","SCRC","SCR ","ECHB","STPB","CCPB","SCPB","ECPB","GPID","PART","EVNT"};
  if(idet < 0) return -1;
  if(cbank) strncpy(cbank,cbanks[(idet/10)],4);
  return idet%10;
}

int TDetectorHit::getDetectorBank(char *cbank, int *irec) {
  if(irow < 0 || ibank < 0 || !cbank || !irec) return -1;
  static const char *cbanks[]={"STRE","STR ","CC1 ","SCRC","SCR ","ECHB","STPB","CCPB","SCPB","ECPB","GPID","PART","EVNT"};
  strncpy(cbank,cbanks[(ibank%100)],4);
  *irec = ibank/100;
  return irow;
}

//------------ setDetectorHit  ---------------------------------------------------------
// copy info from GPID, TBID banks
// return code: -1= no info found;  0= all info found;  1= Edep info missing; 2= pos.&dir. info missing
int TDetectorHit::setDetectorHit(int det_type, GPID_t gpid, double Tcorr) {
  if(gpid.sec<1 || gpid.sec>6) return -1;
  if(gpid.trkid<=0 || gpid.trkid>TBID_NH) return -1; 
  idet = det_type + 10*bGPID;
  itrk = TBID[gpid.trkid-1].track;
  
  switch(det_type) {
  case det_sc:
    {
      Index = gpid.sec*100 + gpid.paddle;
      Trec  = gpid.sc_time;
      Time  = (Tcorr==0.0)? gpid.sc_time : Tcorr; 
      Path  = gpid.sc_len;
      Edep  = gpid.dedx;
      Status= gpid.sc_stat;
      
      irow = TBID[gpid.trkid-1].sc_id -1;
      int is = 0;
      if(irow < 0) return 2;
      if(SCRC_NS) {
	while( (is < SCRC_NS) && (SCRC_S[is] != gpid.sec) ) is++;
	if(is >= SCRC_NS) return 2;
	if(irow >= SCRC_NH[is]) return 2;
	ibank = is*100 + bSCRC;
	Position.SetXYZ(SCRC[is][irow].x,SCRC[is][irow].y,SCRC[is][irow].z);
	PosDet = SCRC[is][irow].y;
      }
      else if(SCR_NS) {
	while( (is < SCR_NS) && (SCR_S[is] != gpid.sec) ) is++;
	if(is >= SCR_NS) return 2;
	if(irow >= SCR_NH[is]) return 2;
	ibank = is*100 + bSCR;
	Position.SetXYZ(SCR[is][irow].x,SCR[is][irow].y,SCR[is][irow].z);
	PosDet = SCR[is][irow].y;
      }
      else 
	return 2;
    }
    break;
  case det_st:
    {
      Index = gpid.sec*100;
      Trec  = gpid.st_time;
      Time  = (Tcorr==0.0)? gpid.st_time : Tcorr; 
      Path  = gpid.st_len;
      Status = gpid.st_stat;

      irow = TBID[gpid.trkid-1].st_id -1;
      if(irow < 0) return 2;
      int ipd = -1;
      if(STRE_NS) {
	int irec = STRE_NS-1;                     //take the last record (i.e. filled after TBT)
	if(irow >= STRE_NH[irec]) return 2;
	ibank = irec*100 + bSTRE;
	ipd = (STRE[irec][irow].ID/100)*8+(STRE[irec][irow].ID%100)/5+STRE[irec][irow].ID%10-10;
	//	ipd = (STRE[irec][irow].ID/100-1)*8+(STRE[irec][irow].ID%100)/10*2-2+STRE[irec][irow].ID%10;
	//	ipd = (STRE[irec][irow].ID/100-1)*4+(STRE[irec][irow].ID%100)/10;
	if( (ipd < 1) || (ipd > 48) ) {
	  std::cout<<" STRE["<<irec<<"] stpd out of range "<<ipd<<" nhits="<<STRE_NH[irec]<<" STRE[0]:"<<STRE[0][irow].ID<<std::endl;
	  ipd = 0;
	}
	else {
	  Trec  = STRE[irec][irow].ST_TIME;
	  Path  = STRE[irec][irow].ST_L;
	  Edep  = STRE[irec][irow].st_edep;
	  PosDet = STRE[irec][irow].st_pos;
	  //	  if(itrk!=STRE[irec][irow].Trk_no)   //Trk_no=TBTR? - looks like trk no in sector
	  //	    std::cout<<"itrk="<<itrk<<"  STRE["<<irec<<"].Trk_no="<<STRE[irec][irow].Trk_no<<" [0]="<<STRE[0][irow].Trk_no<<std::endl;

	}
      }
      else if( irow < STR_NH ){
	ibank = bSTR;
	ipd = STR[irow].ID;
	if(ipd<1 || (ipd>12 && ipd <100)) {
	  std::cout<<" STR stpd out of range "<<ipd<<" "<<STR_NH<<" "<<STR[irow].ID<<std::endl;
	  ipd = 0;
	}
	else if(ipd > 100) {
	  //	  ipd = (STR[irow].ID/100-1)*4+(STR[irow].ID%100)/10;
	  ipd = (STR[irow].ID/100)*8+(STR[irow].ID%100)/5+STR[irow].ID%10-10;
	  if(ipd<0 || ipd>48) {
	    std::cout<<" STR stpd out of range "<<ipd<<" "<<STR_NH<<" "<<STR[irow].ID<<std::endl;
	    ipd=0;
	  }
	}
	if(ipd > 0) {
	  Trec  = STR[irow].ST_TIME;
	  Path  = STR[irow].ST_L;
	  PosDet= STR[irow].st_pos;
	  // if(itrk==STR[0][irow].Trk_no)   //Trk_no=HBTR?
	}
      }
      Index = gpid.sec*100 + ipd;
      int itdpl = getTrackPlane(det_st,itrk);
      if(itdpl >= 0) {
	int ir = itdpl/100;
	itdpl = (itdpl%100); 
	Position.SetXYZ(TDPL[ir][itdpl].x, TDPL[ir][itdpl].y, TDPL[ir][itdpl].z);
      }
    }
    break;
  case det_ec:
    // not done yet
    idet = det_ec;
  }
  return Index;
}

int TDetectorHit::setDetectorHit(int det_type, PART_t part, double Tcorr) {
  int tbrow = part.trkid -1;
  if((tbrow < 0) || (tbrow >= TBID_NH)) return -1;
  itrk = TBID[tbrow].track;
  idet = det_type + 10*bPART;

  switch(det_type) {
  case det_sc:
    {
      irow = TBID[tbrow].sc_id -1;
      if(irow < 0) return -1;
      int is = 0;
      if(SCRC_NS) {
	while( (is < SCRC_NS) && (SCRC_S[is] != TBID[tbrow].sec) ) is++;
	if(is >= SCRC_NS) return -1;
	if(irow >= SCRC_NH[is]) return -1;
	ibank = is*100 + bSCRC;
	Index = TBID[tbrow].sec*100 + SCRC[is][irow].id;
	Trec  = TBID[tbrow].sc_time;
	Time  = (Tcorr==0.0)? Trec : Tcorr; 
	Edep  = SCRC[is][irow].energy;
	Status= TBID[tbrow].sc_stat;
	Position.SetXYZ(SCRC[is][irow].x,SCRC[is][irow].y,SCRC[is][irow].z);
	PosDet = SCRC[is][irow].y;
      }
      else if(SCR_NS) {
	while( (is < SCR_NS) && (SCR_S[is] != TBID[tbrow].sec) ) is++;
	if(is >= SCR_NS) return -1;
	if(irow >= SCRC_NH[is]) return -1;
	ibank = is*100 + bSCR;
	Index = TBID[tbrow].sec*100 + SCRC[is][irow].id;
	Trec  = TBID[tbrow].sc_time;
	Time  = (Tcorr==0.0)? Trec : Tcorr; 
	Edep  = SCR[is][irow].energy;
	Status= TBID[tbrow].sc_stat;
	Position.SetXYZ(SCR[is][irow].x,SCR[is][irow].y,SCR[is][irow].z);
	PosDet = SCR[is][irow].y;
      }
      else 
	return 2;
      //get path length from GPID or TDPL
      if(GPID_NH) {  
	Path = GPID[tbrow].sc_len;
      }
      else {
	int itdpl = getTrackPlane(det_sc,itrk);
	if(itdpl >= 0) {
	  int ir = itdpl/100;
	  itdpl = (itdpl%100); 
	  Path = TDPL[ir][itdpl].tlen;
	}
      }
    }
    break;
  case det_st:
    {
      irow = TBID[tbrow].st_id -1;
      if(irow < 0) return -1;
      int ipd = -1;
      if(STRE_NS) {
	int irec = STRE_NS-1;
	if(irow >= STRE_NH[irec]) return -1;
	ibank = irec*100 + bSTRE;
	ipd = (STRE[irec][irow].ID/100)*8+(STRE[irec][irow].ID%100)/5+STRE[irec][irow].ID%10-10;
	//	ipd = (STRE[irec][irow].ID/100-1)*4+(STRE[irec][irow].ID%100)/10;
	if( (ipd < 1) || (ipd > 48) ) 
	  ipd=0;
	else {
	  Trec  = STRE[irec][irow].ST_TIME;
	  Path  = STRE[irec][irow].ST_L;
	  Edep  = STRE[irec][irow].st_edep;
	  PosDet= STRE[irec][irow].st_pos;
	    // if(itrk==STRE[irec][irow].Trk_no)   //Trk_no=HBTR? or TBTR? or trk in sec.?
	}
      }
      else if( irow < STR_NH ){
	ibank = bSTR;
	ipd = STR[irow].ID;
	if(ipd<1 || (ipd>12 && ipd<100)) 
	  ipd = 0;
	else if(ipd > 100) {
	  ipd = (STR[irow].ID/100)*8+(STR[irow].ID%100)/5+STR[irow].ID%10-10;
	//	  ipd = (STR[irow].ID/100-1)*4+(STR[irow].ID%100)/10;
	  if(ipd<0 || ipd>48) ipd=0;
	}
	if(ipd > 0) {
	  Trec  = STR[irow].ST_TIME;
	  Path  = STR[irow].ST_L;
	  PosDet= STR[irow].st_pos;
	  // if(itrk==STR[irow].Trk_no)   //Trk_no=HBTR?
	}
      }
      Index = TBID[tbrow].sec*100 + ipd;
      int itdpl = getTrackPlane(det_st,itrk);
      if(itdpl >= 0) {
	int ir = itdpl/100;
	itdpl = (itdpl%100); 
	Position.SetXYZ(TDPL[ir][itdpl].x, TDPL[ir][itdpl].y, TDPL[ir][itdpl].z);
      }
    }     
  }
  return Index;
}

int TDetectorHit::setDetectorHit(int det_type, EVNT_t evnt, double Tcorr) {
  idet = det_type + 10*bEVNT;
  // itrk=?

  switch (det_type) {
  case det_sc:
    {
      irow = evnt.SCstat -1;
      if(irow >= 0 && irow < SCPB_NH) {
	ibank = bSCPB;
	Index = int(SCPB[irow].ScPdHt/100);
	Edep = SCPB[irow].Edep;
	Trec = SCPB[irow].Time;
	Time = (Tcorr < 0.0)? Trec : Tcorr;
	Path = SCPB[irow].Path;
	Status = SCPB[irow].Status;
	int id = (SCPB[irow].ScPdHt%100) -1;
	int is = 0;
	if(SCRC_NS) {
	  while( (is < SCRC_NS) && (SCRC_S[is] != Index/100) ) is++;
	  if(is >= SCRC_NS) return 2;
	  if(id < 0 || id >= SCRC_NH[is]) return 2;
	  Position.SetXYZ(SCRC[is][id].x,SCRC[is][id].y,SCRC[is][id].z);
	  PosDet = SCRC[is][id].y;
	}
	else if (SCR_NS) {
	  while( (is < SCR_NS) && (SCR_S[is] != Index/100) ) is++;
	  if(is >= SCR_NS) return 2;
	  if(id < 0 || id >= SCR_NH[is]) return 2;
	  Position.SetXYZ(SCR[is][id].x,SCR[is][id].y,SCR[is][id].z);
	  PosDet = SCR[is][id].y;
	}
	else
	  return 2;
      }
    }
    break;
  case det_st:
    {
      irow = evnt.STstat -1;
      if(irow >= 0 && irow < STPB_NH) {
	ibank = bSTPB;
	int stid = STPB[irow].SThid%10000;
	Index = int(stid/100)*100+(stid/100)*8+(stid%100)/5+(stid%10)-10;
	Trec = STPB[irow].Time;
	Time = (Tcorr < 0.0)? Trec : Tcorr;
	Path = STPB[irow].Path;
	Status = STPB[irow].Status;
	if(STRE_NS) {
	  for(int i=0; i<STRE_NH[STRE_NS-1]; i++) {
	    if(STRE[STRE_NS-1][i].ID==STPB[irow].SThid) {
	      Edep  = STRE[STRE_NS-1][i].st_edep;
	      PosDet= STRE[STRE_NS-1][i].st_pos;
	      break;
	    }
	  }
	}
	else if( STR_NH > 0) {
	  for(int i=0; i<STR_NH; i++) {
	    if(STR[i].ID==STPB[irow].SThid) {
	      PosDet= STR[irow].st_pos;
	      break;
	    }
	  }
	}
      }
    }
  }
  return Index;
}

int TDetectorHit::setDetectorHit(int det_type, int ind, double time, double path, double edep){
  idet = det_type;
  Index= ind;
  Trec = time;
  Time = time;
  Path = path; 
  Edep = edep;
  return Index;
 }

//---------------------------------------------------------------------------

int TDetectorHit::getTrackPlane(int det_type, int trkno) {      //trkno=1,...,ntrk
  int trkpln = -1;
  if( (trkno <= 0) || !TDPL_NS) return -1;
  if( (trkno <= TBER_NH) || (trkno <= TBTR_NH) ) {
    int itrsec = (TBER_NH>0)? (TBER[trkno-1].layinfo2>>8)&0xff : (TBTR[trkno-1].itr_sec%100);
    int isec   = (TBER_NH>0)? (TBER[trkno-1].layinfo2>>24)     : int(TBTR[trkno-1].itr_sec/100);
    if((itrsec > 0) && (isec > 0) && (isec < 7)) {
      int pln_min=0, pln_max=0;
      if(det_type == det_st) {
	pln_min=2; pln_max=3;
      }
      else if(det_type == det_dc) {
	pln_min = 4; pln_max = 15;  //only region 1
      }
      else if(det_type == det_cc) {
	pln_min = pln_max = 40;
      }
      else if(det_type == det_sc) {
	pln_min = 41; pln_max = 44;
      }
      else if(det_type == det_ec) {
	pln_min = 45; pln_max = 46;
      }
      for(int irec=0; irec<TDPL_NS; irec++) {
	if( (TDPL_S[irec] == isec) && (itrsec <= TDPL_NH[irec]/10) ) {
	  for(int i=0; i<TDPL_NH[irec]; i++) {
	    if( (TDPL[irec][i].trk_pln>=trkno*100+pln_min) && (TDPL[irec][i].trk_pln<=trkno*100+pln_max) && (TDPL[irec][i].x<900) ) {
	      trkpln = irec*100 + i;
	      break;
	    }
	  }
	}
      } 
    }
  }
  return trkpln;
}

//---------- CalculateBeta (e.g. SC-ST) ------------------------------------------------------

double TDetectorHit::CalculateBeta( TDetectorHit* det ) {
  if(Time > det->Time && Path > det->Path) {
    return (Path - det->Path)/(Time - det->Time)/Clight;
  }
  else if(Time < det->Time && Path < det->Path) {
    return (det->Path - Path)/(det->Time - Time)/Clight;
  }
  else {
    return -1;
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#ifdef __CINT__ 
#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 
#pragma link C++ class TDetectorHit+; 
#pragma link C++ class TMyPart+; 
#endif 

#endif
