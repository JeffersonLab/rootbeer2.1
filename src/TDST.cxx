//--Author	Ken Lvingston	Mar 2002                             
//--Update	Ken Lvingston	Sep 2005  
//		Moved completely to root                           

#include "TDST.h"


ClassImp(TDST)

//////////////////////////////////////////////////////////////////////////
//     
// TDST   
// 
// Class to handle root BANK DST files within root
//
// See RootBeer main clas
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TDST::TDST(const char *bosFile, void* addresses, int total)
{
// Create one TDST object and open the 1st data file
  fError=0;

  fBankAddr = (struct addressBanks_t*)addresses;	//all info about bank variables
  fBankTotal = total;					//total no of defined banks
  fDSTfile = new TFile(bosFile);			//open DST input file
  if(fDSTfile->IsZombie()){                             //if failed return error
    fError=1;
    return;
  }
  fBankTree = (TTree*)gDirectory->Get("bankTree");	//get the tree
  if(fBankTree==NULL){
    fprintf(stderr,"failed to open banktree \n");
    return;
  }
  fBigBuff = new char[100000];		//create a 100k buffer to hold the event
  fDataStartAddr = fBigBuff;     	//init the data pointer
  fDataAddr = fBigBuff;			//ditto
  fServedBankNo=0;			//init the no of banks being handled
  fEntry=0;
  fMaxEntry=(int)fBankTree->GetEntries(); //get the no of events in the DST file
  for(int n=0;n<fBankTotal;n++){      	//init the status of all banks to zero
      fBankStatus[n]=0;
  }
}	

//______________________________________________________________________________
TDST::~TDST()
{
  delete fBigBuff;
  delete fDSTfile;
  Clear();
}

//______________________________________________________________________________
int TDST::GetEvent()
{
  // Read the next event from the tree and make all the global BANK variables
  // (who's addresses are in fServedBankTab) point to the data

  void **sect;		//pointer to sector addresses
  int *nRows,nBank;	

  char *secAddr; 	//
  fDataAddr=fDataStartAddr;		//point running ptr to start of buffer
  for(int n =0;n <fServedBankNo;n++){	// loop round all required banks
    // Make all variables to point to the data for the current event
    if(fBankType[n]==MULTIPLE_BANK){		// ie multi instance bank
      for(int b=1;b<5;b++){			// get all the static branches
	fServedBranches[n][b]->GetEntry(fEntry);
      }
      fServedBranches[n][0]->SetAddress(fDataAddr); //set address for current bank
      fServedBranches[n][0]->GetEntry(fEntry);

      secAddr=fDataAddr;			//set the 1st sector address
      nBank = *fServedBankTab[n]->nBank;      	//get no of bank / sectors
      nRows = *fServedBankTab[n]->nRow;		//get array of rows per sector
      sect = (void**)*fServedBankTab[n]->address;//get array of sector pointers	
      for(int s=0;s<nBank;s++){		//loop over all sector banks
	sect[s] = secAddr;		//set address to data for that sector
	secAddr+=(nRows[s]*fRowLength[n]); //move ptr to next sector
      }

      fDataAddr+=(fRowTotal[n]*fRowLength[n]); //set ptr for reading next bank data
    }
    else{					// ie single instance bank ...
      fServedBranches[n][4]->GetEntry(fEntry);	// ... as for MULTI but simpler
      if(fServedBranches[n][5]) fServedBranches[n][5]->GetEntry(fEntry);
      fServedBranches[n][0]->SetAddress(fDataAddr);
      fServedBranches[n][0]->GetEntry(fEntry);
      *fServedBankTab[n]->address=fDataAddr;
      fDataAddr+=(fRowLength[n]*(*(int *)fServedBankTab[n]->nRow));      
    }
  }
  fEntry++;				//increment no of entries
  if(fEntry>=fMaxEntry) fEntry*=-1;	//return a -ve number if finished
  return(fEntry);
}

int TDST::StartServer(){
  //make up the list of served banks and associated variables
  //and get all the relevant branches set up to the relevan addresses
  char name[8];		//general name string	
  char brName[20];	//branch name string

  SetBankStatus("HEAD",1);	//force the reading of the HEAD bank even if not set by user
  for(int n=0;n<fBankTotal;n++){
    if(fBankStatus[n]==1){
      sprintf(name,"%.4s",fBankAddr[n].name);
      sscanf(name,"%s",name);	//in case less that 4 chars get it terminated properly
      
      fRowLength[fServedBankNo]=fBankAddr[n].length;		//save the length for later use
      if(fBankTree->GetBranchStatus(name)){			//check that bank is in the tree
	strcpy(fBankList[fServedBankNo], fBankAddr[n].name);	//copy the name
	sprintf(fBankBranch[fServedBankNo],"%.4s",fBankAddr[n].name);
	sscanf(fBankBranch[fServedBankNo],"%s",fBankBranch[fServedBankNo]);  //in case less that 4 chars get it terminated
	fBankTree->SetBranchStatus(fBankBranch[fServedBankNo],1);   	//set the branch status to on
	fServedBranches[fServedBankNo][0] = fBankTree->GetBranch(fBankBranch[fServedBankNo]); //save the branch
	
	fServedBankTab[fServedBankNo]=&fBankAddr[n];

	if(fBankAddr[n].nBank==NULL){	//ie a single instance / sector bank
	  fBankType[fServedBankNo]=0;	//Meaning single
 
	  sprintf(brName,"%s_NH",fBankBranch[fServedBankNo]); //branch for no of rows 
	  fBankTree->SetBranchAddress(brName,(int*)fServedBankTab[fServedBankNo]->nRow);
	  fServedBranches[fServedBankNo][4] = fBankTree->GetBranch(brName);

	  sprintf(brName,"%s_S",fBankBranch[fServedBankNo]);  //branch for sector no
	  //for backward compat. see if BANK_S exists, if not make it 0
	  if(fBankTree->GetBranchStatus(brName)){
	    fBankTree->SetBranchAddress(brName,(int*)fServedBankTab[fServedBankNo]->bankSec);
	    fServedBranches[fServedBankNo][5] = fBankTree->GetBranch(brName);
	  }
	  else{
	    fServedBranches[fServedBankNo][5]=NULL;
	    *(int*)fServedBankTab[fServedBankNo]->bankSec=0;  
	  }
	}
	else{				//ie. a multi sector / instance bank

	  *fBankAddr[n].address=fDataStartAddr;
	  fDataStartAddr+=(MAX_BANK_NUMBER*sizeof(void *));

	  sprintf(brName,"%s_NH",fBankBranch[fServedBankNo]);
	  fBankTree->SetBranchAddress(brName,fDataStartAddr);
	  *fServedBankTab[fServedBankNo]->nRow=(int *)fDataStartAddr;
	  fDataStartAddr+=(MAX_BANK_NUMBER*sizeof(int *));
	  fServedBranches[fServedBankNo][4] = fBankTree->GetBranch(brName);
  
	  sprintf(brName,"%s_S",fBankBranch[fServedBankNo]);
	  fBankTree->SetBranchAddress(brName,fDataStartAddr);
	  *fServedBankTab[fServedBankNo]->bankSec=(int *)fDataStartAddr;
	  fDataStartAddr+=(MAX_BANK_NUMBER*sizeof(int *));
	  fServedBranches[fServedBankNo][2] = fBankTree->GetBranch(brName);

  
	  fBankType[fServedBankNo]=1;	//Meaning multiple
  
	  sprintf(brName,"%s_T",fBankBranch[fServedBankNo]);
	  fBankTree->SetBranchAddress(brName,&fRowTotal[fServedBankNo]);
	  fServedBranches[fServedBankNo][3] = fBankTree->GetBranch(brName);

	  sprintf(brName,"%s_NS",fBankBranch[fServedBankNo]);
	  fBankTree->SetBranchAddress(brName,fServedBankTab[fServedBankNo]->nBank);
	  fServedBranches[fServedBankNo][1] = fBankTree->GetBranch(brName);

	}

	fBankIndex[fServedBankNo]=n;	//save the index of the bank
	fServedBankNo++;		//increment the no of server banks
      }
    }
  }
  return(0);
}

