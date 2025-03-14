//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Bos                                                                  //
//                                                                      //
// Class to handle BOS files within root                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef ROOTBEER_TBos
#define ROOTBEER_TBos


#include "TRootBeer.h"
#include "BosTypes.h"
#include "TThread.h"

class TBos : public TRootBeer {

private:
  TThread    	*fDataThread; 		//thread pointer
  TThread    	*fScannerThread; 	//thread pointer

 public:
  struct eventInf_t 	*fERing;
  Char_t     	fBankList[MAX_BANK_TYPES][5];
  Int_t	      	fBankType[MAX_BANK_TYPES];
  Char_t      	fDataFile[80];		// hold the name of the current data file 
  Int_t		fNanalysed;		// count no of analysed events
  Int_t       	fSwapFlag;		// flag to indicate byte swapping
  Int_t	       	fRingEvent;		// current ring event index
  Int_t	       	fPrevRingEvent;		// previous ring event index
  Int_t	       	fServedBankNo;		// count the no of banks to be served
  Int_t	     	fBufferStatus[2];	// is the buffer READ, SCANNED, or ANALYSED 
  Int_t	     	fBufferNread[2];	// no of bytes read into buffer 
  char 		*fRBuff[2];	       	// Two buffers to take records read
  Int_t		fKillScanner;		// Flag to end the theads
  Int_t		fKillDataServer;       	// Flag to end the theads
  Int_t		fIsScanner;		// Flag to say if thread running
  Int_t		fIsDataServer;       	// ditto
  
  TBos(const char*, void*, Int_t);	//class constructor
  Int_t GetEvent();			//main function called for each event
  static void WordSwap(char *, int);	//if written on big(little)endian and read on little(big)...
  static void DataSwap(char *);		//...endian, some byte swapping needs to be done
  Int_t StartServer();			//initialise some things after constructor called
  virtual ~TBos();			//class destructor
  
  static  void *DataServer(void *arg);			//some debugging print functions for bos files
  static  void *BankScanner(void *arg);
  static  void printRHead(struct RecHead_t *s);
  static  void printRSHead(struct recSegHead_t *s);
  static  void printDataSHead(struct dataSegHead_t *s);

  ClassDef(TBos,0)  //Bos structure
};

#endif
