//--Author	Ken Lvingston	Mar 2002                             
//--Update	Ken Lvingston	Sep 2005  
//		Moved completely to root                           

#include "TRootBeer.h"

ClassImp(TRootBeer)

//////////////////////////////////////////////////////////////////////////
//
// TRootBeer - Base Class for TBos and TDST    
//                         
// Class to handle CLAS BANK data files within root. This clas should not 
// be instantiated, rather used a base clas.
//
// The input files can either be standard BOS files, or 
// root files where the BANKs have been written in Tree format
// using the associated DSTWriter class
//
// This package is independent of the CLAS software.
//
// Philosophy and use
// BOS BANK data maps naturally to a C struct.
// A nice efficient way to handle the data for a BANK is to make 
// a struct type corresponding to each bank and point the struct to the 
// data.
// 
// Making classes corresponding to each bank is neither neccessary, appropriate
// nor efficient. ( ...In my opinion)
//
// In this scheme each BANK has some global variables which allow access
// to the BANK data in the simplest possible way. The RootBeer::GetEvent()
// function fills variables by doing a lot of shuffling of void** type things.
// Example of the variable are:
//			 
// 1) Single instance bank (eg EVNT):			
//	struct EVNT_t *EVNT		//1d array of structures	
//	int            EVNT_NH; 	//no of hits in the bank
//	int            EVNT_S; 		//sector no of bank not ususlly relevant for singles
//
//
// 2) Multiple instance bank (eg DC0):				
//	struct SCR_t    **SCR;		//2d array of structs
//	int              *SCR_NH;	//array of hits for each sector 
//	int              *SCR_S;	//aray of sector numbers
//	int               SCR_NS;	//number of sectors 
//								
//
// Here's an example of the main event loop accessing and filling histograms for these banks
//
//    while(event>=0){                                            // do all events
//      event = rootbeer->GetEvent();                             // get next event from the server
//
//      //eg of handling of single sector bank TAGR
//      for(int row=0;row<TAGR_NH;row++){                         // loop over rows in TAGR bank
//        tagr_hist->Fill(TAGR[row].T_id,TAGR[row].TPHO);         // fill 2d hist of time vs id
//      }
//
//      //eg of handling of multiple sector bank SCR
//      for(int sector=0; sector<SCR_NS;sector++){                // for all sectors
//        sect_no=SCR_S[sector];                                  // get the number of the sector
//        for(int nhit=0;nhit<SCR_NH[sector];nhit++){             // for all hits in that sector
//          scr_hist[sect_no]->Fill(SCR[sector][nhit].id);        // fill  hist
//        }
//      }
//      eventTot++;
//      if(eventTot%1000 == 0) fprintf(stdout,"done %d\r",eventTot);
//      if((nEvents>0)&&(eventTot >=nEvents)) break;              //break if nEvents done
//    }
//
//
//		
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
  TRootBeer::TRootBeer(const char *File, void* addresses, int total){
  fBankAddr = (struct addressBanks_t*)addresses;
  fBankTotal = total;
}	

//______________________________________________________________________________
TRootBeer::TRootBeer(){
}	

//______________________________________________________________________________
TRootBeer::~TRootBeer(){
}

//______________________________________________________________________________
int TRootBeer::GetEvent(){
  return(0);
}
      
int TRootBeer::GetBankStatus(const char *name){
  //returns the status of a bank (ie whether it's to be accessed (1) or not (0)

  int n;
  
  for(n=0;n<fBankTotal;n++){
    if(strncmp(name,fBankAddr[n].name,4)==0){
      break;
    }
  }
  if(n==fBankTotal){
    fprintf(stderr,"Don't know about %.4s bank, can't GetBankStatus(%.4s)\n",name,name);
    exit(0);
  }
  return(fBankStatus[n]);

}

int TRootBeer::SetBankStatus(char **names,int status){
  //Sets the status of a list of banks (ie whether they're to be accessed (1) or not (0)
  //names can be like this: 
  //char **names ={"all"};
  //or like this
  //char **names ={"SC  ","TAGR",...,...,"null"}; (where last strung must be "null")
  //status should be 1 or 0 to switch on or off
 
  int n,b=0;
  //if "all" set all bank status
  if(strncmp(names[0],"all",3)==0){
    for(int m=0;m<fBankTotal;m++){
      fBankStatus[m]= status;
    }
    return(0);
  }
  //else go through list of names until "null"
  while(strncmp(names[b],"null",4)!=0){
    for(n=0;n<fBankTotal;n++){
      if(strncmp(names[b],fBankAddr[n].name,4)==0){
	break;
      }
    }
    if(n==fBankTotal){
      fprintf(stderr,"Don't know about %s bank, in rootbeer->SetBankStatus(....,""%4s"",...)\n",names[b],names[b]);
      exit(0);
    }	
    else{
      fBankStatus[n]=status;
      b++;
    }
  }
  return(0);
}    

int TRootBeer::SetBankStatus(const char *name, int status){
  int m,n;

  if(strncmp(name,"all",3)==0){
    for(m=0;m<fBankTotal;m++){
      fBankStatus[m]= status;
    }
    return(0);
  }
  
  for(n=0;n<fBankTotal;n++){
    //fprintf(stderr,"%.4s\n",fBankAddr[n].name);
    if(strncmp(name,fBankAddr[n].name,4)==0){
      break;
    }
  }
  if(n==fBankTotal){
    fprintf(stderr,"Don't know about bank, can't SetBankStatus(%.4s, %d)\n",name,status);
    exit(0);
  }
  fBankStatus[n]=status;
  return(0);
}
int TRootBeer::ListServedBanks(){
  //print to stderr a list of the banks which are being served
  fprintf(stderr,"The following BANKs will be handled (if they are in the data file):\n");
  for(int m=0;m<fBankTotal;m++){
    if(fBankStatus[m]==1) fprintf(stderr,"\"%.4s\"\n",fBankAddr[m].name);
  }
  return(0);
}

int TRootBeer::StartServer(){
  return(0);
}

