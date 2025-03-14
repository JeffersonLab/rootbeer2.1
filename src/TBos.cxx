//--Author    	Ken Lvingston 	Mar 2002                            
//--Update    	Ken Lvingston  	Sep 2005
//	      	Moved completely to root
//--Update    	Ken Livingston	Jan 2006
//		Fixed bug in byte swapping for events split between segments
//		Fixed bug in ~TBos() related to killing threads
//              Fixed another bug in ~TBos() related to killing threads (Aug 2017)
//		Fixed bug in GetEvent() TThread::Unlock() never got called after last event.
//


#include "TBos.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

ClassImp(TBos)

//////////////////////////////////////////////////////////////////////////
//                                                                   
// TBos
//           
// Class to handle BOS files within root
//
//
// There are 2 threads running: DataServer() and BankScanner()
//
// DataServer() reads large chunks of data into a pair of swinging buffers
//
// BankScanner() scans the data in the free buffer to get all the information 
// about the BANKs for an event. It fills an event ring buffer (of structs)
// with the appropriate information.
//
//
// GetEvent() takes the next event from the ring buffer and makes the global BANK variables
// (who's addresses are in fServedBankTab) point to the data
//
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TBos::TBos(const char *bosFile, void* addresses, int total) :
          fDataThread(NULL), fScannerThread(NULL)
{
  fError=0;
  int   path=-1;                                   // path for open()
// Create one TBos object and open the 1st data file
  
  fBankAddr = (struct addressBanks_t*)addresses;
  fBankTotal = total;

  strcpy(fDataFile,bosFile); 	// make local copy of the input data file name

  struct stat fileStat;
  if(stat(fDataFile,&fileStat) <0 || (path = open(fDataFile,O_RDONLY)) < 0){  
    fprintf(stderr,"Error: can't open file - %s\n",fDataFile);       // exit if can't open 
    fError=1;
    return;
  }
  if(fileStat.st_size < 400) {
    fprintf(stderr,"Error: empty file - %s\n",fDataFile);       // exit
    fError=1;
    return;
  }
  if(path!=-1) close(path);

  fERing =new eventInf_t[MAX_STORED_EVENTS]; // Create the event ring buffer
  //create the swinging data buffers
  fRBuff[0]= new char[MAX_BYTES_PER_BUFF_READ+SPLIT_BANK_SPACE];
  fRBuff[1]= new char[MAX_BYTES_PER_BUFF_READ+SPLIT_BANK_SPACE];
  fServedBankNo=0;		// zero counter for the no of banks to be served
  fRingEvent=0;			//init the ring event no
  fPrevRingEvent=-1;		//init the previous ring event no
  fNanalysed=0;			//init the no of events analysed
  fBufferStatus[0]=BUFFER_ANALYSED;	 	// set buffers as analysed and ready for raw data
  fBufferStatus[1]=BUFFER_ANALYSED;	     	// set buffers as analysed and ready for raw data
  if(!fDataThread){
    fDataThread=new TThread((void(*) (void *))&DataServer,(void*)this);
   }
  if(!fScannerThread) fScannerThread=new TThread((void(*) (void *))&BankScanner,(void*)this);
}	

//______________________________________________________________________________

TBos::~TBos()
{

  //ugly fudge to kill the threads if they haven't finished.
  if(fIsScanner)  fKillScanner=1;
  while(fKillScanner){    //wait tile both are dead
    usleep(100);
  }
  if(fIsDataServer) fKillDataServer=1;
  
  while(fKillDataServer){    //wait tile both are dead
    usleep(100);
  }

  if(fDataThread){			//delete threads
    TThread::Delete(fDataThread);
    //delete fDataThread;
    fDataThread=0;
  }
  if(fScannerThread){
    TThread::Delete(fScannerThread);
    //delete fScannerThread;
    fScannerThread=0;
  }
  if(fERing) delete fERing;		//delet other things
  if(fRBuff[0]) delete fRBuff[0];
  if(fRBuff[1]) delete fRBuff[1];
  TThread::CleanUp();
  Clear();
}

//______________________________________________________________________________
int TBos::GetEvent()
{
// GetEvent() takes the next event from the ring buffer and makes the global BANK variables
// (who's addresses are in fServedBankTab) point to the data

  int timeout=0;

  while(fERing[fRingEvent].status<EVENT_SCANNED){	//wait til next event is scanned
    gSystem->Sleep(10);
    //usleep(100);
    timeout++;
    if(timeout>1000){
      fprintf(stderr,"WARNING: timed out waiting for events. Ended file. Assumed it ended in some ugly way!\n");
      return (-1*fNanalysed);
    }
  }
  
  if(fPrevRingEvent>-1){				//if past 1st event
    TThread::Lock();					//lock out other threads before writing to data members
    if(fERing[fPrevRingEvent].buffer>=EVENT_LAST_IN_BUFFER){			       	//if last in buffer
      fBufferStatus[fERing[fPrevRingEvent].buffer-EVENT_LAST_IN_BUFFER]=BUFFER_ANALYSED;//mark buffer as analysed
    }
    if(fERing[fPrevRingEvent].status==EVENT_LAST){	//if last event
      TThread::UnLock();					//lock out other threads before writing to data members
      return (-1*fNanalysed);				//return (the -ve indicates last event)
    }
    else{
      fERing[fPrevRingEvent].status=EVENT_ANALYSED;	// flag the prev event as analysed
    }
    TThread::UnLock();					// unlock for data writing from threads
  }
  for(int n =0;n <fServedBankNo;n++){			// loop round all required banks

    // Make all bank variables to point to the data for the current event
    if(fBankType[n]==MULTIPLE_BANK){			// ie multi instance bank
      *fServedBankTab[n]->nBank = fERing[fRingEvent].BTab[n].number[MAX_BANK_NUMBER];
      if(fERing[fRingEvent].BTab[n].number[MAX_BANK_NUMBER]>0){ 
	*fServedBankTab[n]->address = fERing[fRingEvent].BTab[n].address;
	*fServedBankTab[n]->bankSec = fERing[fRingEvent].BTab[n].number;
	*fServedBankTab[n]->nRow = fERing[fRingEvent].BTab[n].nrow;
      }
    }
    else{						// ie single instance bank
	*(int *)fServedBankTab[n]->nRow = 0;
       if(fERing[fRingEvent].BTab[n].nrow[0]>0){
	*fServedBankTab[n]->address = fERing[fRingEvent].BTab[n].address[0];
	*(int *)fServedBankTab[n]->nRow = fERing[fRingEvent].BTab[n].nrow[0];
	*(int *)fServedBankTab[n]->bankSec = fERing[fRingEvent].BTab[n].number[0];
      }
    }
  } 
  fPrevRingEvent=fRingEvent;			//save prev event index
  fRingEvent++;					//move round the event ring
  fRingEvent%=MAX_STORED_EVENTS;

  fNanalysed++;					// inc the count of analysed events
  return(fNanalysed);
}


int TBos::StartServer(){
  //makes a list of served banks and sets up a related tables of structs with all addresses
  //of the Global BANK variables
  SetBankStatus("HEAD",1);	//force the reading of the HEAD bank even if not set by user
  for(int n=0;n<fBankTotal;n++){
    if(fBankStatus[n]==1){
      if(fBankAddr[n].nBank==NULL) fBankType[fServedBankNo]=0; 	//meaning single
      else fBankType[fServedBankNo]=1; 				//meaning multiple
      strncpy(fBankList[fServedBankNo], fBankAddr[n].name,4);	//copy the name
      fServedBankTab[fServedBankNo++]=&fBankAddr[n];		//copy the address of the struct
    }
  }

  //  TThread::SetCancelAsynchronous();
  //  TThread::SetCancelOn();
  //now create threads for DataServer and BankScanner
  fKillScanner=0;
  fKillDataServer=0;
  fIsScanner=1;
  fIsDataServer=1;
  fDataThread->Run();				//start the data server thread
  fScannerThread->Run();			//start the bank scanner thread
  return(0);
}


void *TBos::DataServer(void *arg){
  // Run as a thread. Waits until the current data buffer is analysed then fills
  // it with raw data and flags the buffer as READ_RAW. Flips to the next buffer
  // and does the same .... etc between the 2 raw data buffers buffers 

  int 	path;					// path for open() 
  int 	currentBuffer=0;		       	// nubmer of current buffer for reading 
  int 	lastBuffRead;				// flag 
  long	ntot=0;
  unsigned int recordLen;			// length of record in bytes 
  int readBytes;				// bytes per fileread 
  struct  RecHead_t RH;	      			// for first fecord in file, to get length 
  struct  RecHead_t *recH;    

  TBos *b=(TBos*)arg;				// fudge to allow this thread access to data members  

  //we open and read the header to figure out if btye swapping is needed
  if((path = open(b->fDataFile,O_RDONLY)) < 0){	       		// open the data file 
    fprintf(stderr,"Error: can't open file - %s\n",b->fDataFile);	// exit if can't open 
    TThread::Lock();
    b->fIsDataServer=0;
    TThread::UnLock();
    return NULL;
  }
  read(path,b->fRBuff[currentBuffer],sizeof(RecHead_t));	//read the 1st header
  recH=(RecHead_t*)b->fRBuff[currentBuffer];		

  RH.nrecwd=recH->nrecwd;			//figure out if byte swapping needed
  RH.nusedw=recH->nusedw;			//by looking at nrecwd
  WordSwap((char *)&RH,2);			//swapping and trying again
      
  if(abs(RH.nrecwd) < recH->nrecwd){
    b->fSwapFlag=SWAP;				//no swapping needed
    recordLen=WORD*RH.nrecwd;			//set record length in bytes
  }
  else{
    b->fSwapFlag=NO_SWAP;		     	//swapping needed
    recordLen=WORD*recH->nrecwd;		//set record length in bytes
  }
  
  close(path);					//close the file

  if((path = open(b->fDataFile,O_RDONLY)) < 0){	       			// reopen the data file 
    fprintf(stderr,"Error: can't open file - %s\n",b->fDataFile);	// exit if can't open 
    TThread::Lock();
    b->fIsDataServer=0;
    TThread::UnLock();
    return NULL;
  }
  

  readBytes=((int)MAX_BYTES_PER_BUFF_READ/recordLen)*recordLen;	//set no of bytes per file read
  
  for(;;){					// loop going 'til file finished 
    
    while(b->fBufferStatus[currentBuffer]!=BUFFER_ANALYSED){    //wait 'til current buffer is analysed
      usleep(1000);
      if(b->fKillDataServer){
	break;
      }
    }

    if(b->fKillDataServer){
      break;
    }

    if((b->fBufferNread[currentBuffer] = 
	read(path, b->fRBuff[currentBuffer], readBytes)) > 0){ 	// read file 'til done 
      ntot+=b->fBufferNread[currentBuffer];
      if(b->fBufferNread[currentBuffer]%recordLen!=0){ 		// check that int no of records read 
	fprintf(stderr,"WARNING last record incomplete in file %s\n",b->fDataFile);
	TThread::Lock();
	b->fBufferNread[currentBuffer]-=recordLen;
	b->fBufferStatus[currentBuffer]=BUFFER_LAST;	 	//flag buffer as last 
	TThread::UnLock();
	//	fprintf(stderr,"closed file\n");
      }
      else if(b->fBufferNread[currentBuffer]!=readBytes){
	TThread::Lock();
	b->fBufferStatus[currentBuffer]=BUFFER_LAST;		//flag buffer as last and ready for scanning
      //fprintf(stderr,"READ LAST BUFFER in file %s\n",b->fDataFile);
	TThread::UnLock();
      }
      else{
	TThread::Lock();
	b->fBufferStatus[currentBuffer]=BUFFER_READ_RAW;		//flag buffer as ready for scanning
	TThread::UnLock();
      }
      
      lastBuffRead=currentBuffer;      		// save the last buffer read 
      
      //fprintf(stderr,"Read %d to buffer %d\n",currentBuffer,b->fBufferNread[currentBuffer]);
      currentBuffer=(lastBuffRead+1)%2;		// ie flips between buffers 0 and 1 
    }
    else{
      break;	
    }
  }

  // all files processed, print any log stuff here 
  //fprintf(stderr,"Ended Neatly after %ld bytes \n",ntot);

  close(path);
  TThread::Lock();
  b->fKillDataServer=0;
  b->fIsDataServer=0;
  TThread::UnLock();
  return NULL;
}


void *TBos::BankScanner(void *arg){
  // scans the next data buffer for information on the relevan banks for each event
  // and fills an event ring (of structs)

  struct RecHead_t *RHead=NULL;				// -> Record header 
  struct recSegHead_t *RSHead=NULL;    			// -> Record segment header 
  struct dataSegHead_t *DataSHeadFirst=NULL;	       	// -> first data segment header  
  struct dataSegHead_t *DataSHead=NULL;			// ->Data segment header 
  char *buffEnd=NULL; 					// -> end of raw data buffer 
  char *recHead=NULL, *recEnd=NULL, *recUsedEnd=NULL;	// -> addresses of records  
  char *recSegHead=NULL, *recSegEnd=NULL;		// -> addresses of record segments 
  char *dataSegHead=NULL, *dataSegEnd=NULL;		// -> addresses of data segments 
  char *dataSegHeadFirst=NULL;	       		// -> addresses of first data segments if split 
  //  char *data, *dataEnd;			// -> addresses of data 
  char *prevBEnd=NULL;			// -> end of previous raw data buffer 
  int prevBuf=0;			// number of previous raw data buffer (0 or 1) 
  int curBuf=0;				// number of current raw data buffer (0 or 1) 
  int prevEvent=-1;
  int event=-1;	
  int nextEvent=-1;			// event no within array of processed events 
  int bank=0;				// bank counter 
  int prevBank=0;			// index (in desired bank list) of previous bank processed 
  int numCount=0;			// number of occurances of same bank type withing event 
  int n=0;				// the ubiquitous n etc
  int newevent=NEW_EVENT;		// flag that segment is assoriated with a new event 
  int savecurrentbank=0;		// flag that the bank status is 1 for the current bank
  int currentSegmentCode=0;		// save the code of the current data segment
  
  TBos *b=(TBos*)arg;				// fudge to allow this thread access to data members 

  for(;;){
    while(b->fBufferStatus[curBuf]<BUFFER_READ_RAW){		// wait for buffer 
      usleep(10);
      if(b->fKillScanner){
	    TThread::Lock();
	    b->fKillScanner=0;
	    b->fIsScanner=0;
	    TThread::UnLock();
	return NULL;
      }
    }
    recHead = b->fRBuff[curBuf];		      	// point to 1st record header 
    buffEnd = recHead+b->fBufferNread[curBuf];        	// point to end of buffer 
    
    while(recHead < buffEnd){				// do all records in buffer 

      if(b->fSwapFlag==SWAP){                         	// swap all the bytes if needed 
        WordSwap(recHead,2);
      }

      RHead = (struct RecHead_t *)recHead;		// point struct at record header, whose
							//   byte swapping done by data server 
							//   if needed 
      //printRHead(RHead);


      recEnd=recHead+(WORD*RHead->nrecwd);		// save end of record address 
    
      recUsedEnd=recHead+(WORD*RHead->nusedw);		// find address of data end 
      recSegHead=recHead+sizeof(RecHead_t);		// find address of 1st rec segment header 

      while(recSegHead < recUsedEnd){			// do all segments in record 
	RSHead = (struct recSegHead_t *)recSegHead;	// point struct at rec seg header 

	if(b->fSwapFlag==SWAP){		       		//swap all the bytes if needed
	  WordSwap(recSegHead,3);
	  WordSwap(recSegHead+20,6);
	}
	//printRSHead(RSHead);
	if((RSHead->segment_code == SEG_COMPLETE)||(RSHead->segment_code==SEG_FIRST)) {	
	  if(prevEvent>-1){					// if not the first event  
	    TThread::Lock();
	    b->fERing[prevEvent].status=EVENT_SCANNED;		// flag prev event is scanned 
	    TThread::UnLock();
	  }
	  
	  prevBank=0;						// reset the prevBank counter 
	  bank=0;						// reset the bank counter 
	  if(event>-1)prevEvent=event;		       		// save this as previous event
	  event++;						// move round the event ring
	  event%=MAX_STORED_EVENTS;

	  while(b->fERing[event].status!=EVENT_ANALYSED){	//wait 'til next space in ring is free
	    usleep(10);
	    if(b->fKillScanner){
	      TThread::Lock();
	      b->fKillScanner=0;
	      b->fIsScanner=0;
	      TThread::UnLock();
	      return NULL;
	    }
	  }

	  for(n=0;n<b->fServedBankNo;n++){
	    b->fERing[event].BTab[n].number[MAX_BANK_NUMBER]=0; // zero the counter 
	    b->fERing[event].BTab[n].nrow[0]=0;			// zero the no of rows 
	  }
	  newevent = NEW_EVENT;			       		// flag that it's a new event 

	  b->fERing[event].buffer=curBuf;
	  if(curBuf!=prevBuf){		       	// If now into a new buffer, flag event as last in buff
	    //printRSHead(RSHead);
	    b->fERing[event].buffer=prevBuf+EVENT_LAST_IN_BUFFER;
	    b->fBufferStatus[prevBuf]=BUFFER_SCANNED;
	    //fprintf(stderr,"event %d, buff %d\n",event,fERing[event].buffer);
	    prevBuf=curBuf;
	  }	
	}
	
	recSegEnd=recSegHead+sizeof(recSegHead_t)+(WORD*RSHead->nsegwd);  // point to end of rec seg head 
	dataSegHead=recSegHead+sizeof(recSegHead_t);		// point to start of first data seg header 
      
	while(dataSegHead < recSegEnd){				// do data segments in current record segment 
	  DataSHead = (struct dataSegHead_t *)dataSegHead;	// point struct at date seg header ..
								// this has the info about the BANK 
          if(b->fSwapFlag==SWAP){                         	//swap all the bytes if needed
          WordSwap(dataSegHead,1);
          WordSwap(dataSegHead+12,6);
	  }
	  currentSegmentCode=DataSHead->segment_code;		//make copy of the segment code

	  dataSegEnd=dataSegHead+WORD*(DataSHead->nwords+DataSHead->nsegwd_this); // save end of segment 
	  if(dataSegEnd>buffEnd){
	    break;				// break here if the record was incomplete 
	  }
	  //printDataSHead(DataSHead);
	  if(currentSegmentCode < SEG_MIDDLE){			//first or only data segment of BANK 
	    dataSegHeadFirst=dataSegHead;				//save address of first seg 
	    DataSHeadFirst=(struct dataSegHead_t *)dataSegHeadFirst;   	//struct->  first seg 
	    // first save pointers to all the banks for the event in the order that they appear. There's
	    //   no real overhead in doing this, since they need to be found to go through the record 
	    b->fERing[event].BHead[bank]=DataSHead;	// save the pointer to data segment  header 
	    b->fERing[event].BData[bank++]=dataSegHead+(WORD*DataSHead->nwords); // -> data 
	    b->fERing[event].nBanks=bank;		// save the no of the bank within the event 

	    savecurrentbank=0;				// assume that we don't save this bank

	    // Now save pointers to the actual data segments of ONLY the banks of interest. This is to allow
	    // efficient, rapid access to these. The array char *fBankList[] should contain the names of
	    // the banks of interest in the order which they appear in the recored 
	    for(n=0;n<b->fServedBankNo;n++){
	      if(strncmp(DataSHead->name,b->fBankList[n],4) == 0){  	// compare name with fBankList[] 
		savecurrentbank=1;					//flag that we save the current bank
		if((n==prevBank)&&(newevent==OLD_EVENT)) numCount++;	// if same as prev BANK and not 1st 
									// BANK in an event, inc numCount 
		else {				// otherwise 
		  numCount=0;			// zero numCount; since  1st  of this BANK in the event 
		  newevent=OLD_EVENT;		// and flag that ist not the first bank in the  event 
		}
		//printDataSHead(DataSHead);
		if(numCount>=MAX_BANK_NUMBER){
		  fprintf(stderr,"Warning %dth instance of %s bank not handled; increase MAX_BANK_NUMBER and recompile bankserver\n",numCount,b->fBankList[n]);
		}
		else{
		  b->fERing[event].BTab[n].number[numCount]=DataSHead->number; // Store the BANK number 
		  b->fERing[event].BTab[n].nrow[numCount]=DataSHead->nrow;     // store the no of rows 
		  b->fERing[event].BTab[n].number[MAX_BANK_NUMBER]++;	// Keep count of how many 
		  b->fERing[event].BTab[n].address[numCount]=dataSegHead+(WORD*DataSHead->nwords); // -> data 
		}
		prevBank=n;	// start the next name search at the position ot this BANK in fBankList[] 
		break;		// break since name was found 
	      }
	    }
	    
	  }
	  else{	       	// BANK was split between record segments. This happens when the first 
		       	// record segment in a record is a middle, or last part of an event.
		       	// We don't want to save any pointers to this, all we need to do is 
		       	// shuffle the appropriate chunk of data back up the buffer a bit until
		       	// it's contiguous with the first chunk of data for that BANK. It will 
		       	// overwrite the record, record segment, and data segment headers, but
		       	// nobody cares at this stage!
		       	// Or, in the case of the first BANK in a new raw buffer, it will go to the 
		       	// end of the previous raw buffer (which should have enough space to hold it 

	    if(savecurrentbank){				//if we're saving the current bank
	      DataSHeadFirst->nsegwd_this+=DataSHead->nsegwd_this; //add the no of words to 1st headed
	    
	      if(WORD*(DataSHead->nsegwd_this)>SPLIT_BANK_SPACE){	//some huge no of rows  -too big for buffer 
		printDataSHead(DataSHead);
		fprintf(stderr,"Too many rows in %s, setting hits = 0\n",DataSHead->name);
		b->fERing[event].BTab[n].nrow[numCount]=0;	// set the no of rows = 0 
	      }
	      else{
		memmove(prevBEnd,dataSegHead+(WORD*DataSHead->nwords), // move the data to end of prev 
			WORD*(DataSHead->nsegwd_this));
	      }
	    }
	  }
	  if((currentSegmentCode == SEG_COMPLETE)||(currentSegmentCode == SEG_LAST)){
	    // if it's last or only data segment fo the BANK 
	    if(b->fSwapFlag==SWAP){
	      DataSwap(dataSegHeadFirst);
	    }
	  }
	  dataSegHead=dataSegEnd; 	// set address of next data segment header 
	  prevBEnd=dataSegEnd;		// save end of this one 
	}
	recSegHead=recSegEnd;		// set address of next record segment header 
      }

      recHead=recEnd;			// set address of next record header 
    }
    if(b->fBufferStatus[curBuf]==BUFFER_LAST){
      //fprintf(stderr,"Found last buffer event = %d\n",event);
      b->fERing[prevEvent].status=EVENT_SCANNED;	// flag prev event as scanned
      b->fERing[event].status=EVENT_LAST;		// flag the event as scanned, and last
      nextEvent=event+1;				// also flag next event as last
      nextEvent%=MAX_STORED_EVENTS;	    		// If at end of event buff, goto start 
      while(b->fERing[nextEvent].status!=EVENT_ANALYSED){	//wait 'til next space in ring is free
	usleep(10);
      }
      b->fERing[nextEvent].status=EVENT_LAST;		// flag the event as scanned, and last
      TThread::Lock();
      b->fIsScanner=0;
      TThread::UnLock();
      return NULL;
    }
    else{
      curBuf++;
      curBuf%=2;
    }
  }
}

void TBos::WordSwap(char *word, int nwords){
  //swap words
  char temp[4];	
  int n;

  swab(word,word,4*nwords);
  for(n=0;n<nwords;n++){
    memmove(temp,word+(4*n),2);
    memmove(word+(4*n),word+(4*n)+2,2);
    memmove(word+(4*n)+2,temp,2);
  }
  //  fprintf(stdout,"wordswap\n");
  return;
}
void TBos::DataSwap(char *dataSegHead){

  // try to do the byte swapping efficiently 
  // if no A or B variables just swap all the words
  // B's are treated specially 
  // A's - ie characters don't need swapped, but we need to swap the 
  // any other variables in the rows 
  // at the moment, just specials for GPAR and EPIC 
  // need more specials if other such banks get defined 
  int row,nrow;
  char *word;
  char format[100];
  struct   dataSegHead_t *DSH=(dataSegHead_t *)dataSegHead;

  strncpy(format,DSH->format,4*(DSH->nwords)-36);
  strcpy(format+4*(DSH->nwords)-36,"");
  if((strchr(format,'A')==NULL)&&(strchr(format,'B')==NULL)){
    WordSwap(dataSegHead+(WORD*DSH->nwords),DSH->nsegwd_this);
  }
  else{
    if(strncmp(format+1,"16",2)==0)
      swab(dataSegHead+(WORD*DSH->nwords),dataSegHead+(WORD*DSH->nwords),
	   4*DSH->nsegwd_this);
    if(strncmp(format+1,"32",2)==0)
      WordSwap(dataSegHead+(WORD*DSH->nwords),DSH->nsegwd_this);
    //couple of specials for banks with char variables
    //need to swap the F,I's but not the chars ... a pain
    if(strncmp(DSH->name,"EPIC",4)==0){
      word=dataSegHead+(WORD*DSH->nwords);
      nrow=DSH->nrow;
      for(row=0;row<nrow;row++){
	WordSwap(word,1);
	word+=36;
      }
    }
    if(strncmp(DSH->name,"GPAR",4)==0){
      word=dataSegHead+(WORD*DSH->nwords);
      nrow=DSH->nrow;
      for(row=0;row<nrow;row++){
	WordSwap(word,2);
	word+=28;
      }
    }
  }
  //  fprintf(stdout,"dataswap\n");
  return;
}

void TBos::printRHead(struct RecHead_t *s){
  fprintf(stderr,"RecordHeader:\n\tnrecwd: \t%d\n\tnusedwd:\t%d\n",
	  s->nrecwd,
	  s->nusedw);
}
void TBos::printRSHead(struct recSegHead_t *s){
  fprintf(stderr,"RecordSegmentHeader:\n\tswapping:\t%d\n\tformat:\t\t%d\n\trecord_no:\t%d\n\tname:\t\t%.8s\n\tnumber1:\t%d\n\tnumber2:\t%d\n\tphysical_no:\t%d\n\tdataword:\t%d\n\tsegment_code:\t%d\n\tnsegwd:\t\t%d\n\n", 	s->swapping,
	  s->format,
	  s->record_no,
	  s->name,
	  s->number1,
	  s->number2,
	  s->physical_no,
	  s->dataword,
	  s->segment_code,
	  s->nsegwd);
}

void TBos::printDataSHead(struct dataSegHead_t *s){
  int len;
  char format[20];

  len = s->nwords - 9;
  sprintf(format,"%%.%ds\n",4*len);

  fprintf(stderr,"Data Segment Header:\n\tnwords:\t\t%d\n\tname:\t\t%.4s\n\tnumber\t\t%d\n\tncol:\t\t%d\n\tnrow:\t\t%d\n\tsegment_code:\t%d\n\tnsegwd_prev:\t%d\n\tnsegwd_this:\t%d\n\tformat:\t\t",
	  s->nwords,
	  s->name,
	  s->number,
	  s->ncol,
	  s->nrow,
	  s->segment_code,
	  s->nsegwd_prev,
	  s->nsegwd_this);
  fprintf(stderr,format,s->format);
}
