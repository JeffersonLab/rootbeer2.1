//
//                  *** TExpTable.h ***

//--Description
//  Class for handling lines from a simple text table

//--Author      K Livingston      8th March  2007
//--Update      K Livingston     10th Sept   2009
//
//  See source code for documentation.
//
//---------------------------------------------------------------------------
#ifndef __TExpTable_h__
#define __TExpTable_h__

#include "TObject.h"

//enumerate the data types
enum {EHInt,EHFloat,EHString};

struct VarNames{
  char* VarName;
  void* VarAddress;
};

class TExpTable{
 private:
 protected:
  void 		**fVariableAddresses;	//Array of addresses of variables to be updated from run table info
  Char_t 	**fVariableNames;	//Array of names of variables
  Char_t 	**fVariableFormats;	//Array of format strings for scanning in each variable
  Int_t 	*fVariableTypes;	//Array of Types
  Int_t 	fNVariables;		//No of variable columns in the run table
  Int_t 	fNFormats;		//No of format columns in the run table
  Char_t 	fConfigFileName[120];	//Name of the Run Table
  Bool_t 	fHaveVariables;		//flag that a Variables: line has alread been processed
  Bool_t 	fHaveTypes;		//flag that a Types: line has alread been processed
  Char_t	fSetupFileName[200];	//name of run table file
  FILE          *fp;                    //pointer the file
  Char_t        fTableName[200];        //name of the run table

 public:
  TExpTable( const char* SetupFile, const char * TableName);  	        //Constructor
  virtual int GetAllForValue(int value, char *variable1, char *variable2=NULL, int range=0);
  virtual int GetAllForValue(float value, char *variable1, char *variable2=NULL, float range=-1.0);
  virtual int GetAllForValue(char *value, char *variable1, int perfectmatch=0);
  //  virtual int GetAllInRangeFast(int value, int *low, int* high);

  virtual void *GetVarAddress(char *variable);
  virtual void *SetVarAddress(char *variable, void *address );

  virtual ~TExpTable();	     		//Destructor
  ClassDef(TExpTable,1)
};
#endif

