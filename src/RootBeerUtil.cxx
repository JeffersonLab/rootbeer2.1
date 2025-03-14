//some utility functions which do miscellaneous things
//relating to rootbeer and will be called from a standard executable (see rbtest.C)

#include <sys/stat.h>
#include <sys/types.h>
#include "RootBeerUtil.h"
#define	DONT_KNOW_YET	       -1
#define	ON_COMMAND_LINE		0 // source of data files
#define FROM_FILELIST		1
#define FILES_AVAILABLE		0 // file availability flag
#define	FILES_ALL_DONE		1

static int mode=DONT_KNOW_YET;
static FILE *fp_filelist=NULL;
static int file_no=0;

class TRootBeer *rootbeer=NULL;	//the main instance of the class

//************************* createBeerObject() ********************************
TRootBeer *createBeerObject(char *bankFileName){
  Int_t error=0;
  //work out whether it's a bos file or a root file and create the
  //appropriate object
  if((strstr(bankFileName,".root"))!=NULL){
    rootbeer = new TDST(bankFileName,bankAddress,nBankTypes);
  }
  else{
    rootbeer = new TBos(bankFileName,bankAddress,nBankTypes); 
  }
  error=rootbeer->GetError();
  if(error){
    fprintf(stderr,"WARNING, ROOT thinks there was some error opening file %s\n", bankFileName);
    fprintf(stderr,"If there are no further errors, it's probably safe to ignore this warning\n");
    //    delete rootbeer;
    //return NULL;
  }
  return rootbeer;
}
//******************************************************************************

//************************* getFileMode() **************************************
int getFileMode(char *filearg){
  // set mode for single file or filelist
  //if filelist, open the list file

  char flist[255];

  if(strncmp(filearg,"-L",2)==0){		      	// ie -L<filelist>
    if((fp_filelist = fopen(filearg+2,"r"))!=NULL){ 	// open <filelist>
      strcpy(flist,filearg+2);				// copy <filelist> to flist
      mode=FROM_FILELIST;				// set mode
    }
    else{						//if can't open <filelist>
      fprintf(stderr,"Error: can't open file - %s\n",filearg+2);
      exit(-1);
    }
  }
  else{						// just a single input file argument
    strcpy(flist,filearg);			// copy the infile arg
    mode=ON_COMMAND_LINE;           		// set mode 
  }
  return 0;
  
}
//*********************************************************************************

//************************* getNextFile() *****************************************
int getNextFile(char *nextfile, char *filearg){
  // get the next file
  char line[255];

  if(mode==DONT_KNOW_YET) getFileMode(filearg);

  if(mode==FROM_FILELIST){
    for(;;){			    			// scan lines from filelist
      if((fgets(line,255,fp_filelist)) == NULL){
	fclose(fp_filelist);
	return -1;
      }	
      if((line[0] == '*')||(line[0] == '#')){	 	// skip comment lines 
	continue;
      }
      else{
	sscanf(line,"%s",nextfile);	   		// put the name in
	struct stat fileStat;
	if(stat(nextfile,&fileStat) < 0) {
	  fprintf(stderr,"*** File not found: %s\n",nextfile);
	  continue;
	}
	//to prevent empty ROOT file (only ROOT header in file):
	if(fileStat.st_size < 400) {
	  fprintf(stderr,"*** File too short: %s\n",nextfile);
	  continue;
	}
	return ++file_no;
	break;
      }
    }
  }
  else{					// on command line
    strcpy(nextfile,filearg);		// copy the infile arg
    if(file_no>0) 
      return -1;                // only one file
    else {
      struct stat fileStat;
      if(stat(nextfile,&fileStat) < 0) return -1;
      //to prevent empty ROOT file (only ROOT header in file):
      if(fileStat.st_size < 400) return -1;
      return ++file_no;
    }
  }
  return 0;
}
//*********************************************************************************

//************************* resetNextFile() **************************************
int resetNextFile(){
  mode=DONT_KNOW_YET;
  fp_filelist=NULL;
  file_no=0;
  return(0);
}
//*********************************************************************************

//************************* getEpics() *******************************************
//function to return the value of an EPICS channel
//returns -1 if it fails to find the channel
//otherwise, returns the sector no and stores the value
int getEpics(char *name, float *value){
  for(int sector=0;sector<EPIC_NS;sector++){	//for all sectors
    for(int row=0;row<EPIC_NH[sector];row++){	//for all "hits" is epics channels	
      if(strncmp(name,EPIC[sector][row].char1,strlen(name))==0){ //check to see if it's the one you want
	*value=EPIC[sector][row].value;		// store the value
	return EPIC_S[sector];			// return the sector no (0 or 1 I think)
      }
    }
  }
  return -1;		//return -1 if the channel wasn't found
}
//*********************************************************************************

//************************* getBit() *******************************************
//function to return status of a bit within an integer word. (Eg for TRIGBIT)
Int_t getBit(int word,int bit){
  return 1&(word>>bit);
}
//*********************************************************************************
