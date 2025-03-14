//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DST                                                                  //
//                                                                      //
// Class to handle ROOTBEER DST files within root                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef ROOTBEER_TDST
#define ROOTBEER_TDST


#include "TRootBeer.h"

#ifndef __CINT__
#include <TTree.h>
#include <TFile.h>
#endif

class TDST : public TRootBeer {

private:
 public:
  class TBranch* fServedBranches[MAX_BANK_TYPES][6]; 	//branches for served banks
  Int_t fBankType[MAX_BANK_TYPES];			//types of server banks
  Int_t fBankIndex[MAX_BANK_TYPES];			//indices of served banks
  Int_t fRowLength[MAX_BANK_TYPES];			//lengths rows in served banks 
  char fBankList[MAX_BANK_TYPES][8];			//names of served banks
  char fBankBranch[MAX_BANK_TYPES][8];			//names of main served branches 
  Int_t fRowTotal[MAX_BANK_TYPES];			//no of rows in each bank for event	
  class TFile *fDSTfile;				//input file
  class TTree *fBankTree;				//tree with bank data
  char *fBigBuff;					// huge data buffer
  char *fDataStartAddr;					// pointer within buffer
  char *fDataAddr;					// ditto
  Int_t fServedBankNo;					// no of banks being handled
  Int_t fEntry;						// counter for events
  Int_t fMaxEntry;					// max no of events to be done

  TDST(const char*, void*, int);			//class cnstructor
  Int_t GetEvent();					//main function to be called for each event
  Int_t StartServer();					//some init after class construction
  virtual ~TDST();					//class destructor
  ClassDef(TDST,0)
};

#endif
