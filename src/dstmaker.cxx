//
// -*-*-*-*-*--*-*-*-*-*-*-*-*-*-* dstmaker -*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*-//
//			A part of the rootbeer package			      //
//                                                                            //
//		    Ken Livingstn (k.livingston@physics,gla,ac,uk)            //
//									      //
//	Version		Date			Comment			      //
//	1	 	20th Nov 2002					      //
//									      //
// -*-*-*-*-*--*-*-*-*-*-*-*-*-*-* dstmaker -*-*-*-*-*--*-*-*-*-*-*-*-*-*-*-*-//

// Comments
// 
// It makes a root DST from one or more bos files
// Edit the sections between the ********************'s to customize
//
// See other documents in the rootbeer package

#include "RootBeerUtil.h"
#include "DSTWriter.h"

class RootBeer *rootbeer;	// the main BANK file handler object
class DSTWriter *dstwriter;	// the writer object
class TFile *outFile;           // the file to put the results

int   nEventsTot=0;			// ie do all events in all files
int   event=0;		    	// Event counter
int   eventTot=0;	        // Total event count
int   dstEventCount=0;     	// Events written to DST
char  InFileName[128];	
char  OutFileName[128];	

RootBeer *createBeerObject();	// create instance of rootbeer Class.
int   dropBanks();		// drop banks (ie not write to DST MUST CUSTOMIZE)
int   makeCuts();		// cuts - (MUST CUSTOMIZE)
int   selectBanks();		// Choose banks of interest (MUST CUSTOMIZE)
int   closeUp();	        // Close files etc
int   getNextFile(char *, char *);
void  printUsage();
int   parseArgs(int,char **);
void  dstmaker(int, char *, char *);

#ifdef ROOTEXE
int main(int argc,char **argv){
  if((parseArgs(argc,argv))!=0) return(-1);
  dstmaker(nEventsTot,InFileName,OutFileName);
  return(0);
}
#endif

void dstmaker(int nEvents, char *file, char *outFileName){   // The main function

  char inFile[128];
  int fileNo=0;

  outFile=new TFile(outFileName,"recreate");  			//Open the output file
  dstwriter = new DSTWriter(outFile);				//Create the dstwriter
  dropBanks();							//Drop and banks you don't want in DST (CUSTOMIZE BELOW)

  while((fileNo=getNextFile(inFile,file))!=-1){  		// loop while files are still avialable 
    fprintf(stderr,"Sorting file - %s\n",inFile);
    if((rootbeer=createBeerObject(inFile))==NULL) return; 	// create rootbeer object
    if(selectBanks()!=0) return;             			// select the required banks (CUSTOMIZE BELOW)
    rootbeer->StartServer();                  			// start the server running
    dstwriter->Init(rootbeer);					// Init the writer


    // the main sort loop
    while(event>=0){                                    // do all events
      event = rootbeer->GetEvent();			// get the next event from the server

      if(makeCuts()==0){	        		// Check passed the cuts
	dstwriter->Write();                     	// Write the event in the Tree
	dstEventCount++;		     		// increment the count	
      }
      eventTot++;
      if(eventTot%1000 == 0) fprintf(stdout,"Sorted\t%d\tWrote\t%d\n",eventTot,dstEventCount);
      if((nEvents>0)&&(eventTot >=nEvents)) break;		//break if nEvents done
    } 
    //End the main sort loop
    delete rootbeer;					// delet the rootbeer object   
    fprintf(stdout,"Read %d events from %s\n",abs(event),inFile); //print the stats 
    if((nEvents>0)&&(eventTot >=nEvents)) break;		//break if nEvents done
  }

  closeUp();						// close files delete things etc
}

// *************************************************************************** 
//                         CUSTOMIZE THIS FUNCTION                           *
//                                                                           *
//                              BANK SELECTION                               * 
//                                                                           *
// *************************************************************************** 
//                                                                           *
int selectBanks(){

// First the ones we want to have available from the bank file
// to use - eg make cuts etc.
//rootbeer->SetBankStatus("all",1);
rootbeer->SetBankStatus("CSQL",1);
rootbeer->SetBankStatus("DC0 ",1);
rootbeer->SetBankStatus("TGPB",1);
rootbeer->SetBankStatus("HBER",1);
rootbeer->SetBankStatus("DC1 ",1);
rootbeer->SetBankStatus("VERT",1);
rootbeer->SetBankStatus("ST1 ",1);
rootbeer->SetBankStatus("DCPB",1);
rootbeer->SetBankStatus("TRKS",1);
rootbeer->SetBankStatus("HDPL",1);
rootbeer->SetBankStatus("HBTR",1);
rootbeer->SetBankStatus("CC  ",1);
rootbeer->SetBankStatus("TRPB",1);
rootbeer->SetBankStatus("TAGR",1);
rootbeer->SetBankStatus("TAGE",1);
rootbeer->SetBankStatus("HBLA",1);
rootbeer->SetBankStatus("HBID",1);
rootbeer->SetBankStatus("ECPO",1);
rootbeer->SetBankStatus("ECPB",1);
rootbeer->SetBankStatus("SCR ",1);
rootbeer->SetBankStatus("SCPB",1);
rootbeer->SetBankStatus("PART",1);
rootbeer->SetBankStatus("ECPC",1);
rootbeer->SetBankStatus("EC1 ",1);
rootbeer->SetBankStatus("TAGT",1);
rootbeer->SetBankStatus("EC  ",1);
rootbeer->SetBankStatus("TGBI",1);
rootbeer->SetBankStatus("TBER",1);
rootbeer->SetBankStatus("SC  ",1);
rootbeer->SetBankStatus("MTRK",1);
rootbeer->SetBankStatus("TAGI",1);
rootbeer->SetBankStatus("SCRC",1);
rootbeer->SetBankStatus("TRL1",1);
rootbeer->SetBankStatus("TDPL",1);
rootbeer->SetBankStatus("TBTR",1);
rootbeer->SetBankStatus("STR ",1);
rootbeer->SetBankStatus("HEVT",1);
rootbeer->SetBankStatus("ECHB",1);
rootbeer->SetBankStatus("DSTC",1);
rootbeer->SetBankStatus("CL01",1);
rootbeer->SetBankStatus("CALL",1);
rootbeer->SetBankStatus("TBLA",1);
rootbeer->SetBankStatus("TBID",1);
rootbeer->SetBankStatus("EVNT",1);
rootbeer->SetBankStatus("ST  ",1);
rootbeer->SetBankStatus("SC1 ",1);
rootbeer->SetBankStatus("HEAD",1);
rootbeer->SetBankStatus("ECPI",1);
rootbeer->SetBankStatus("DHCL",1);
rootbeer->SetBankStatus("STPB",1);

 return(0);
}

int dropBanks(){
  // Now the ones from the above lot which we DON'T WANT to write out 
  // - for example
  //dstwriter->DropBank("MVRT");
  dstwriter->DropBank("PART");

  return 0;
}
//                                                                           *
// ***************************************************************************
//                                                                           *
//                           END OF BANK SELECTION                           * 
//                                                                           *
// ***************************************************************************





// *************************************************************************** 
//                         CUSTOMIZE THIS FUNCTION                           *
//                                                                           *
//                         CUTS AND CONDITION FILTER                         * 
//                                                                           *
//             This function has to return 0 to get the event recorded.      *
// *************************************************************************** 
//                                                                           *

int makeCuts(){
  return(0);
  if((TAGR_NH>0)||(EPIC_NS)){
    return 0;
  }
  else return -1;
}
//                                                                           *
// *************************************************************************** 
//                                                                           *
//                     END OF CUTS AND CONDITION FILTER                      *
//                                                                           *
// *************************************************************************** 



int parseArgs(int argc,char **argv){
  switch(argc){
  case(1):				// 1 or 2 args not enough
  case(2):
    printUsage();
    return(-1);
    break;
  case(3):				// 3 args rbtest <infile> <outfile>
    if((strstr(argv[1],"-N"))!=NULL){
      printUsage();
      return(-1);
    }
    else{
      strcpy(InFileName,argv[1]);
      strcpy(OutFileName,argv[2]);
    }
    break;
  case(4):				// 4 args must be "rbtest -N# <infile> <outfile>" 
    if(strncmp(argv[1],"-N",2)!=0){  	//print usage then exit */
      printUsage();
      return(-1);
    }
    else{
      sscanf(argv[1]+2,"%d",&nEventsTot);	//get the required event no
      strcpy(InFileName,argv[2]);
      strcpy(OutFileName,argv[3]);
    }      
    break;
  default:
    printUsage();
    return(-1);
    break;
  }
  return(0);
}

void printUsage(){
  fprintf(stderr,"\nusage:\n\n");
  fprintf(stderr,"dstmaker -h   Print this message\n\n");
  fprintf(stderr,"dstmaker [-Nevents] <infile> <dstfile> \n");
  fprintf(stderr,"dstmaker [-Nevents] <-Lfilelist> <dstfile> \n");
}


int closeUp(){
  fprintf(stdout,"Wrote %d events to %s\n",dstEventCount,OutFileName);		
  outFile->Write();			// close the files
  outFile->Close();
  return 0;
}
