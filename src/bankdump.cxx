#include "RootBeerUtil.h"

// ******** Put any required root headers here **************************
//                                                                      *
#include "TFile.h"              // for the output root file
//                                                                      *
// **********************************************************************
#define NUMBER          1
#define GROUP           2

void printUsage();
int processArgs(int, char**);
void bankdump(int, char *);

//stuff related to the input args for the code
int opts=4;
char *inargs[4]={"-h","-N","-G","null"};

int groups=0;				//no of banks to be dumped
char *group[MAX_BANK_TYPES]={"all"};	//names of banks to be dumped - default to all
int npfunc=0;				//no of print functions to be called
void (*pfunc[MAX_BANK_TYPES])();	//array of pointers to print functiond

class TRootBeer *rootbeer;         	// The main class

int 	nEventsTot=0;			// ie do all events in all files
int     event=0;                        // Event counter for current file
int	eventTot=0;			// Total event count
char  	InFileName[128];	

// this part is if it's being made as an executable - needs a main() function
#ifdef ROOTEXE
int main(int argc,char **argv){
  if((processArgs(argc,argv))!=0) return(-1);
  bankdump(nEventsTot,InFileName);
  return(0);
}
#endif

int processArgs(int argc,char **argv){
  // Process args 
  int m,n;
  if(argc==1){
    printUsage();
    return(-1);
  }
  if(strcmp(argv[1],"-h")==0){ 		// print usage then exit 
    printUsage();
    return(-1);
  }
  for(n=1;n<argc;n++){			// run though arguments 
    for( m=0;m<opts;m++){		// compare each with list of accepted args 	
      if( strncmp(argv[n],inargs[m],2) == 0 ) break;	// found a match 
    }
    switch(m){
    case(NUMBER):
      if(strlen(argv[n])>=2)		// ie something after "-N" 
      sscanf(argv[n]+2,"%d",&nEventsTot);
      break;
    case(GROUP):
      if(strlen(argv[n])>=2){		// ie something after "-G" 
	group[groups]=new char[5];	// character to hold bank name
	sprintf(group[groups],"    ");	// padd to 4 blanks
	strncpy(group[groups],argv[n]+2,strlen(argv[n]+2)); //copy the bank name
	groups++;			//increment the no of groups
      }
      break;

    default:
      if( argv[n][0]=='-'){
	fprintf(stderr,"\nError: \"%s\" is not a recognised option.\n",argv[n]);
	printUsage();
	return(-1);
      }
      break;
    }
  }
  sprintf(InFileName,argv[argc -1]);	//copy the file name
  fprintf(stderr,"file = %s\n",InFileName);
return (0);
}

void bankdump(int nEvents, char *file){
  int b=0;
  fprintf(stderr,"Sorting file - %s\n",file);
  if((rootbeer=createBeerObject(file))==NULL) return; 	// create rootbeer object
  if(groups){						// if groups > 0
    group[groups]=new char[5];				// append null as the last name
    sprintf(group[groups],"null");
  }
  rootbeer->SetBankStatus(group,ON);			// switch on all the banks in the list

  //now fill a list of pointers to print functions of banks to be dumped
  for(int b=0;b<nBankTypes;b++){			// loop over all bank types
    if(rootbeer->GetBankStatus(bankAddress[b].name)){  	// if status is ON
      pfunc[npfunc++]=bankAddress[b].pfunc;		// store ptr to print function 
    }
  }
  rootbeer->StartServer();                  	    	// start the server running
  
  event=0;						//zero the event counter for the file

  // ********************** main sort loop ******************************************************
  //                                                                                            *
  
  while(event>=0){                            	       	// do all events
    event = rootbeer->GetEvent();                	// get next event from the server
    
    for(b=0;b<npfunc;b++){		// loop over all print functions for banks to be dumped
      pfunc[b]();			// call the print function
    }
    eventTot++;				//increment the total no of events
    if(eventTot%1000 == 0) fprintf(stderr,"done %d\r",eventTot);//print message every 1000 
    if((nEvents>0)&&(eventTot >=nEvents)) break;		//break if nEvents done
  }	
    //                                                                                            *
    // ********************************************************************************************
  fprintf(stderr,"\nSorted %d events from file: %s\n",abs(event),file);
  delete rootbeer;              	//Call the class destructor
  
  return;
  
}

void printUsage(){
  fprintf(stderr,"\nusage:\n\n");
  fprintf(stderr,"bankdump -h   Print this message\n\n"); 
  fprintf(stderr,"bankdump [-Nevents] [-Gbank1 -Gbank2 ...] <dstfile|bosfile> \n");
}
