//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DSTWriter                                                            //
//                                                                      //
// Class to handle writing BANKs in tree in a root file                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef ROOTBEER_TDST_WRITER
#define ROOTBEER_TDST_WRITER

#include "bankvars.h"

#include <TTree.h>
#include <TFile.h>
class TDSTWriter : public TObject {

 private:

 protected:

class TRootBeer *fRootbeer;		// the main object

TTree *bankTree;			// Tree to hold banks 

TFile *fRootFile;      			// root file device
char fBrName[100];			// branch name
char fBrFormat[100];			// branch format string
//int r;	

char *fBigBuff; 			// allocate a 100k buffer for multi instance banks
char *fBuffPtr;				// running pointer within bigBuff
int fUsedBankTot;	       		// total no of banks accessed
int fUsedBankIndex[MAX_BANK_TYPES];  	// indexes in bankAddress[]			
int fUsedBankRows[MAX_BANK_TYPES];   	// total no of rows for mulit instance banks
int fUsedBankType[MAX_BANK_TYPES];	// 0 for singles 1 for multis
char fUsedBankName[MAX_BANK_TYPES][5];	// names of used banks

char fDropBankName[MAX_BANK_TYPES][5];	// names of banks to be dropped (ie not written to DST )
int  fDropBankTotal;			// counter for no od dropped banks
int fDropBankFlag;

int fDSTInitFlag ;			// Flag if dst initislised

 public:
  TDSTWriter(TFile *);
  int Init(TRootBeer *);
  int WriteDST();
  int DropBank(char *);
  TTree *GetTree() const {return bankTree;}
  virtual ~TDSTWriter();
  ClassDef(TDSTWriter,0)
};

#endif
