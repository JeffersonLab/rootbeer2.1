//--Author    Ken Lvingston  	Mar 2002                            
//--Update    Ken Lvingston	Sep 2005  
//	      Moved completely to root                           
//

#include "TDSTWriter.h"


ClassImp(TDSTWriter)

//////////////////////////////////////////////////////////////////////////
//
// TDSTWriter
//
// Class to handle writing BANKs into a root tree
// 
// For each BANK there are some branches:
//
// For single instance banks (eg. TAGR)
// TAGR_NH/I		- An int saying how many rows
// TAGR_S/I		- An int giving the sector number of the bank
// TAGR[TAGR_NH][12]/S	- An array of shorts: bank length (in shorts) x (no of rows)
//
//
// For multiple instance bank (eg. SCR)
// SCR_T/I 		- An int saying how many rows summed over all sectors / instances 
// SCR_NS/I		- An int for the no of sectors / instances of the bank in the event
// SCR_NH[SCR_NS]/I	- An int array with the numbers of rows in each sector 
// SCR_S[SCR_NS]/I	- An int array with the sector number for each sector
// SCR[SCR_T][24]/S	- An array of shorts:  bank length (in shorts) x (total no of rows)
//
//
// single sector banks are contiguous in memory from rootbeer 
// However, multi sector banks are not. Hence, the data from multi sector banks has to
// be copied into a big buffer to form a contiguous thing
// 
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TDSTWriter::TDSTWriter(TFile *file){
  fRootFile=file;
  fDSTInitFlag=0;
  fDropBankFlag=0;
  fDropBankTotal=0;
}

TDSTWriter::~TDSTWriter()
{
  Clear();
}

int TDSTWriter::Init(TRootBeer *rb){
  char name[100];		    	// general name thing
  int nwords;			     	// bank word counter
  int match;			    	// index of matching bank 

  if(fDSTInitFlag==1) return(0);    	//check if already initialised

  fRootbeer=rb;				//->rootbeer calling class
  fUsedBankTot=0;			//init the used bank counter
  fBigBuff = new char[100000]; 		//allocate a 100k buffer for multi sector banks

  fRootFile->cd();			//get into the file working directory
  bankTree = new TTree("bankTree","All the data for the event");	//create the bank tree

  
  for(int n=0;n< nBankTypes;n++){ 	//loop round all known banks
    if(fRootbeer->GetBankStatus(bankAddress[n].name)==1){	//if we using the bank from the input file
      match=0;
      for(int drop=0;drop<fDropBankTotal;drop++){  	//check to see whether to drop it in the output
	if(strncmp(bankAddress[n].name,fDropBankName[drop],4)==0){
	  match=1;
	  fprintf(stdout,"Dropping %s bank from output file\n",fDropBankName[drop]);
	  break;
	}	
      }	
      if(match==1) continue;		//ie dropping it
      sprintf(name,"%.4s",bankAddress[n].name);
      sscanf(name,"%s",name);		//in case less that 4 chars get it terminated properly
      sprintf(fUsedBankName[fUsedBankTot],name);
     
      nwords=bankAddress[n].length/2; 	//get the size of the row in shorts(2byte)

      if(bankAddress[n].nBank==NULL){   	//if single sector bank
	fUsedBankIndex[fUsedBankTot]=n;
	fUsedBankRows[fUsedBankTot]=0;
	fUsedBankType[fUsedBankTot]=0;
       
	//create the required branches
       
	sprintf(fBrName,"%s_NH",name); 		//for total no of rows
	sprintf(fBrFormat,"%s_NH/I",name);	
	bankTree->Branch(fBrName,bankAddress[n].nRow,fBrFormat);

	sprintf(fBrName,"%s_S",name); 		//for sector no (not usually relevant for single sect bank)
	sprintf(fBrFormat,"%s_S/I",name);	
	bankTree->Branch(fBrName,bankAddress[n].bankSec,fBrFormat);
       
	sprintf(fBrName,"%s",name);		//for all the data
	sprintf(fBrFormat,"%s[%s_NH][%d]/S",name,name,nwords);
	bankTree->Branch(fBrName,bankAddress[n].address,fBrFormat);	
      }

      else{	//multi sector banks - I think 18 sectors is the maximum possible (ECPC)?
	
	fUsedBankIndex[fUsedBankTot]=n;
	fUsedBankRows[fUsedBankTot]=0;
	fUsedBankType[fUsedBankTot]=1;
	
	sprintf(fBrName,"%s_T",name);	//For total no of rows in all sectors
	sprintf(fBrFormat,"%s_T/I",name);
	bankTree->Branch(fBrName,&fUsedBankRows[fUsedBankTot],fBrFormat);
	sprintf(fBrName,"%s_NS",name);	//For no of instances of the bank in event
	sprintf(fBrFormat,"%s_NS/I",name);
	bankTree->Branch(fBrName,bankAddress[n].nBank,fBrFormat);
	sprintf(fBrName,"%s_NH",name);	//Array of nos of sector for each instance
	sprintf(fBrFormat,"%s_NH[%s_NS]/I",name,name);
	bankTree->Branch(fBrName,bankAddress[n].nRow,fBrFormat);
	sprintf(fBrName,"%s_S",name);	//Array of nos of sector for each instance
	sprintf(fBrFormat,"%s_S[%s_NS]/I",name,name);
	bankTree->Branch(fBrName,bankAddress[n].bankSec,fBrFormat);
       
	sprintf(fBrName,"%s",name);	//for all the data
	sprintf(fBrFormat,"%s[%s_T][%d]/S",name,name,nwords);
	bankTree->Branch(fBrName,bankAddress[n].address,fBrFormat);	
      }	
      fUsedBankTot++;
    }
  }
  fDSTInitFlag=1;	// flag that the dst is initialised
  return (0);
}

int TDSTWriter::WriteDST(){
  int bank;		//index of the bank
  int *bp;		//->bank data for each sector / instance in multi banks
  int bytes;		//no of bytes to move to put all sector rows contiguous in memory
  void **secAddress;	//array of pointers to data for each bank sector (multi banks) 
  char br[8];		//branch name (surely an inefficient way of doing this ... fix later)

  fRootFile->cd();			//make sure in correct dir / file in case parent has changed
  fBuffPtr = fBigBuff;			//point running ptr to start of buffer
  for(int n=0;n<fUsedBankTot;n++){	//loop over used banks	
    bank=fUsedBankIndex[n];		//save index of bank
    
    if(fUsedBankType[n]==0){		//ie single instance bank
      bankTree->SetBranchAddress(fUsedBankName[n],*bankAddress[bank].address);
    }
    else{			//ie multi instance bank
      bankTree->SetBranchAddress(fUsedBankName[n],fBuffPtr);
      
      sprintf(br,"%s_NH",fUsedBankName[n]);	
      bankTree->SetBranchAddress(br,bankAddress[bank].nRow[0]);
      sprintf(br,"%s_S",fUsedBankName[n]);
      bankTree->SetBranchAddress(br,bankAddress[bank].bankSec[0]);

      fUsedBankRows[n]=0;	//copy all rows in multi sects to one contiguous lot
      secAddress = (void **)(*bankAddress[bank].address);
      for(int s=0;s<bankAddress[bank].nBank[0];s++){
	bp = bankAddress[bank].nRow[0];
	bytes=bankAddress[bank].length*bp[s];
	memcpy(fBuffPtr,secAddress[s],bytes);
	fBuffPtr+=bytes;
	fUsedBankRows[n]+=bp[s];
      }
    }    
  }	
  bankTree->Fill(); //fill the tree
  return 0;
}
int TDSTWriter::DropBank(char *bankName){
  strcpy(fDropBankName[fDropBankTotal++],bankName);
  return 0;
}

