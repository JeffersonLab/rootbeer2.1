//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootBeer                                                             //
//                                                                      //
// Base Class to handle BOS/ROOT_DST files within root                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef TROOTBEER
#define TROOTBEER


#include "TObject.h"
#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>


#include "RootBeerTypes.h"         // structs with pointers to bank variables


class TRootBeer : public TObject {

 private:
  
 public:
  struct addressBanks_t *fServedBankTab[MAX_BANK_TYPES];  //array of structs of all banks to be served
  struct addressBanks_t* fBankAddr;		//struct names and addresses of all bank variables
  Int_t fBankStatus[MAX_BANK_TYPES];		//flags to say which banks are being decoded
  Int_t	fBankTotal;				// total no of defined banks known about
  Int_t fError;                                 // flag for success / failure in class constructor
  
  TRootBeer(const char*, void*, int);
  TRootBeer();
  virtual   Int_t GetEvent();				//called for every event
  virtual   Int_t SetBankStatus(const char*,int);	//add /remove bank name to banks served
  virtual   Int_t SetBankStatus(char**,int);		//add / remove list of banks to banks served
  virtual   Int_t GetBankStatus(const char*);		//say whether a bank is being served or not
  virtual   Int_t ListServedBanks();			//Print list of served banks
  virtual   Int_t StartServer();			//initislise and start serving banks
  virtual ~TRootBeer();					//class destructor
  virtual   Int_t GetError(){return fError;};
  ClassDef(TRootBeer,0)
};

#endif
