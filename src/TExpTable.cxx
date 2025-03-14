//----------------------------------------------------------------------------
//                        ROOTPHYS
//
//                  *** TExpTable.cxx ***
//--Description
//
// Tool for handling lines from a simple text table, where we need to search the
// table for a Row where a specific criteria is met by one of the values. 
//
// What it's intended for:
// It was designed to be used in analysis of data from Jlab, where we might want to analyse up to several
// hundred runs, but each run has different conditions and parameters which can be kept in a text file
// and read on the basis of the current run number. The whole text table file is opened and scanned 
// every time it is accessed. This is fine if it is relatively infrequent - eg . once per run.
// However, it has quite a lot of generality and can be used to extract values for a table in several ways.
// The details are in the function descriptions below.
//
// What it's NOT intended for.
// It is not intended to be accessed every event. Although possible, this would be very slow and inefficient.
// Such a table needs to be resident in memorry and accessed in a more efficient way.
// The aim is to produce another clas TEventTable for this purpose.
//
// A  TExpTable looks like this
// 
// #sample table  - mytable.dat
//
// Variables:   RunMin  RunMax 	Energy 	        Data                    Comment
// Types:         %d      %d      %f            %s                      %d
// #
// Row:            1      50     1.5            data/thing.dat          "Help me I'm drowning"
// Row:           51      90     2.3            data/thing1.dat         -
// Row:           91     105     4.6            data/thing3.dat         -
// .
// .
//
// The table can have as many rows columns as required.
// The column names are specified in the "Variables:" line, and must have no spaces.
// The column types are specified in the "Types:" line, and must be types must be %d, %s, %f 
// For entries in the column which ihas a string type, there must be no spaces, However,
// the last column can be a Comment/Description with (or without) spaces in quotes.
// Every row must have a entry in every column. (No Comment is just a '-' , see above).
// Lines beginning with a # or * are ignored.
//
//
//
// A Table hander is opened by creating a TExpTable class, variables must be defined to correspond
// to the required table columns, and a call to GetAllForValue() updates the variables with
// the values from the relevant line in the table.
//
// For example, in an analysis we might use the above table to set the beam energy depending on 
// the number of run being analysed. Here's some skeleton code.
//
// #include "TExpTable.h"
// .
// .
// Int BeamEnergy=0;  
// TExpTable *myTable = new TExpTable("mytable.dat","mt");
// ...
// myTable->SetVarAddress("Energy",&BeamEnergy);
// .
// for(int run=0;run<100;run++){
//   if(myTable->GetAllForValue(run,"RunMin","RunMax")!=1) continue;
//   AnalyseRun(run,BeamEnergy);
// }
// .     
// . 
// 
//
// Aternatively, the table may contain lines for specific runs eg.
//
// Variables:   Run     Energy 	        Data                    Comment
// Types:       %d      %f              %s                      %d
// #
// Row:         1       1.5             data/thing.dat          "Help me I'm drowning"
// Row:         2       2.3             data/thing1.dat         -
// Row:         3       4.6             data/thing3.dat         -
// .
// .
// .
//
// and the search for the line for the run would look like this:
// if(myTable->GetAllForValue(run,"Run")!=1) continue;
//
//
// 
// table2header
// This is a utility which makes a header file from a table. It is the easiest way to 
// make a table (particularly with many columns) accessible in analysis code.
// 
// eg. >table2vars mytable.dat mt produces a header file mt.h:
// 
// //---------------------------------- mt.h --------------------------------------------//
// #ifndef __mt_h__
// #define __mt_h__
//
// #include "TExpTable.h"
//
// // This header was created from the table mytable.dat on Fri Sep 11 09:58:40 IST 2009 by running the command:
// // /home/kl/rootphys/scripts/table2header mytable.dat mt
//
// // To use it, simple include the header
// // A globally available ExpTable class - mtTable is created
// // together with global variables relating to each column (mt_var) - defined below
// //
// // Calling the member function GetAllForValue(int value, char *variable1, char *variable2=NULL)
// // searches for a row in the table and updates all the global variables to values in that row
// // Usually it would be a Run number or a run range
// // Eg. mtTable->GetAllForValue(53772,"Run"); would find the row where Run=53772
// // and update all the variables to the values from that row.
// // Eg. mtTable->GetAllForValue(53772,"RunMin","RunMax"); would find the row where RunMin <= 53772 <= RunMax
// // and update all the variables to the values from that row.
//
//
// // Create a TExpTable with this file as the actual data table
// TExpTable *mtTable = new TExpTable("mytable.dat","mt");
//
// // Define variables for all the columns in the run table and set their addresses in mtTable 
//
// // An obscure way to define and set global table and variables that works in both
// // interpreted and compiled mode - void *axx are dummies
//
//
// Int_t   mt_Run;                         void *a2  = mtTable->SetVarAddress("Run",&mt_Run);
// Float_t mt_Energy;                      void *a3  = mtTable->SetVarAddress("Energy",&mt_Energy);
// Char_t  mt_Data[200];                   void *a4  = mtTable->SetVarAddress("Data",mt_Data);
// Int_t   mt_Comment;                     void *a5  = mtTable->SetVarAddress("Comment",&mt_Comment);
//
// 
// #endif
//
// //-----------------------------------------------------------------------------------//
//// To make the table accessible in analysis code simply include the file "mt.h"
//  
// Skeleton code like this:
//
// #include "mt.h"
// .
// .
// for(int run=0;run<100;run++){
//   if(mt->GetAllForValue(run,"RunMin","RunMax")!=1) continue;
//   AnalyseRun(run,mt_Energy);
// }
// .     
// .
//
//
//	
//--Author      K Livingston      8th March  2007
//--Update      K Livingston     10th Sept   2009
//--Update      K Livingston     30th Oct    2009
//
//
//---------------------------------------------------------------------------

#include "TExpTable.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

ClassImp(TExpTable)
	
  enum {EExpTableVariables,EExpTableTypes,EExpTableRow};

TExpTable::TExpTable( const char* SetupFile, const char* TableName ){
  //
  // Clas constructor
  // Setupfile is the file containing the table.
  // TableName is the name chosen for the table - this is just used to print out the 
  // information every time a line is gathered from the table.
  // If the table is made into a header using the script ./scripts/table2vars this all gets 
  // take care of.
  // Eg. >./table2vars mytable.dat mt
  // Will make the header mt.h
  // The variables will be mt_RunMin, mt_RunMax
  // and the following line will make TExpTable *mt as a global variable.
  // TExpTable *myTable = new TExpTable("mytable.dat","mt");
  // So to have access to mt in the analysis, it's merely neccessary to #include "mt.h"
  //
  char line[256];		//for lines from the setup file
  Char_t TableType[12];		//The type name for the line of the input file
  char *ptr;  			//General ptr
  char TempString[60];		//General string for scanning variables
  char FormatBase[256]="%*s";   //Format string for scanning variables, skip 1st field for starters		
  int n;
  char fullVarName[200];        //name of the actual variable  
  char qVarName[200];           //name of the column header (in quotes)" 


  fp = NULL;			//init general file ptr
  fNVariables = 0;
  fHaveTypes = kFALSE;
  fNFormats = 0;
  fHaveVariables = kFALSE;
  
  strcpy(fSetupFileName,SetupFile);	       	//save setup file name
  strcpy(fTableName,TableName);                 //save table name
  if((fp = fopen(SetupFile,"r")) == NULL){ 	//open setup file
    fprintf(stderr,"Error, TExpTable: Couldn't open file: %s\n",SetupFile);
    return;
  }
  //look thout file for variable names  and types and set up relevant arrays;
  for(;;){
    if(fHaveTypes&&fHaveVariables) break;	//already got the variable setup lines
    if(fgets(line,256,fp) == NULL) break;  		// check got a line from file
    if((line[0] == '*')||(line[0] == '#')) continue;	// ignore comment lines beginning #
    if((line[3] == '*')||(line[3] == '#')) continue;	// ignore comment lines beginning //# 
    if(sscanf(line,"%s",TableType)!=1) continue;
    //there are 3 possible valuse of TableType "Variables:" "Types:" "Row:"
    if(strstr(TableType,"Variables:")){
      ptr=line+strlen(TableType);		       	//skip over "Variables:" to scan rest of line
      while(sscanf(ptr,"%s",TempString)==1){		//scan rest of line, 1 at a time to count variables
	fNVariables++;
	ptr=strstr(ptr,TempString)+strlen(TempString);
      }
      fVariableNames=new Char_t*[fNVariables];		//create array to hold varaible names 
      fVariableAddresses = new void*[fNVariables];	//create array to hold variable addresses
      fVariableFormats = new Char_t*[fNVariables];	//array to hold format strings for scanning 
      fVariableTypes = new Int_t[fNVariables];	        //array to hold dayta types
 
      for(n=0;n<fNVariables;n++){
	fVariableAddresses[n]=NULL;
	fVariableFormats[n]=NULL;
      }
      ptr=line+strlen(TableType);		       	//skip over "Variables:" to scan rest of line      
      for(n=0;n<fNVariables;n++){			//scan 1 variable name at a time
	if(sscanf(ptr,"%s",TempString)!=1){
	  fprintf(stderr,"Error, TExpTable: Screwed up");
	  fclose(fp);
	  return;
	}
	else{
	  fVariableNames[n]= new Char_t[strlen(TempString)+1]; //  create a string to hold the name
	  sprintf(fVariableNames[n],"%s",TempString);		//save name
	  ptr=strstr(ptr,TempString)+strlen(TempString);
	}
      }
      fHaveVariables=kTRUE;
    }
    else if(strstr(TableType,"Types:")){
      if(!fHaveVariables){
	fprintf(stderr,"Error, TExpTable: reading Exp Table %s ",fSetupFileName );
	fclose(fp);
	return;
      }
      ptr=line+strlen(TableType);		       	//skip over "Types:" to scan rest of line
      for(n=0;n<fNVariables;n++){			//scan 1 variable name at a time
	if(sscanf(ptr,"%s",TempString)!=1){
	  fprintf(stderr,"Error, TExpTable: Screwed up reading Types in %s\n",fSetupFileName );
	  fclose(fp);
	  return;
	}
	else{
	  //check that it's an acceptable format
	  if(strcmp(TempString,"%d")==0) fVariableTypes[n]=EHInt;
	  else if (strcmp(TempString,"%f")==0) fVariableTypes[n]=EHFloat;
	  else if (strcmp(TempString,"%g")==0) fVariableTypes[n]=EHFloat;
	  else if (strcmp(TempString,"%s")==0) fVariableTypes[n]=EHString;
	  else {
	    fprintf(stderr,"Error, TExpTable:  %s is not a valid type (must be %%s %%d %%f, file: %s\n",TempString,fSetupFileName);
	    fclose(fp);
	    return;
	  }
	  fVariableFormats[n]= new Char_t[strlen(FormatBase)+4]; //  create a string to hold the name
	  sprintf(fVariableFormats[n],"%s",FormatBase);		 // save name
	  strcat(fVariableFormats[n],TempString);		 // save name
	  strcat(FormatBase,"%*s");			       
	  ptr=strstr(ptr,TempString)+strlen(TempString);
	  fNFormats++;
	}
      }
      if(fNFormats!=fNVariables){
	fprintf(stderr,"Error, TExpTable: reading Exp Table %s ",fSetupFileName );
	fclose(fp);
	return;
      }
      fHaveTypes=kTRUE;
    }
  }
  fclose(fp);

  for( n=0;n<fNVariables;n++){
    switch(fVariableTypes[n]){
    case EHInt:
      fVariableAddresses[n]=new Int_t;
      break;
    case EHFloat:
      fVariableAddresses[n]=new Float_t;
      break;
    case EHString:
      fVariableAddresses[n]=new Char_t[200];
      break;
    default:
      fprintf(stderr,"Error, TExpTable: reading Exp Table %s ",fSetupFileName );
      fclose(fp);
      return;
      break;
    }
  }

  

  //go though all the variable names requested
  fprintf(stdout,"\nTExpTable: Using table \"%s\" with variables and types like this:\n\n",fSetupFileName);
  fprintf(stdout,"--------------------------------------------------------\n");
  fprintf(stdout,"Col No.\tCol Title\t    Format    Variable Name\n");
  fprintf(stdout,"--------------------------------------------------------\n");
  for( n=0;n<fNVariables;n++){
    sprintf(fullVarName,"%s_%s",fTableName,fVariableNames[n]);
    sprintf(qVarName,"\"%s\"",fVariableNames[n]);


    fprintf(stdout,"%-8d%-20s%-10s%s\n", n, qVarName, strrchr(fVariableFormats[n],'%'), fullVarName);
  }
  fprintf(stdout,"--------------------------------------------------------\n\n");
}



int TExpTable::GetAllForValue(int value, char *variable1,  char *variable2, int range){
  //
  //Look for a line in the table where the value equals the value of variable1
  //or is in the range between variable1 and variable2 (inclusive)
  //and set all variables to be their values from that line in the table 
  //return -1 if error
  //return 1 if only single line found (success)
  //return 0 and print warning if multiple lines found
  //If multiple lines found, use values from the first one found.
  //
  FILE *fp;			//general file ptr
  Char_t line[256];		//for lines from the setup file
  Char_t goodline[256];		//for condition matched line from the setup file
  char *ptr;  			//General ptr
  int match=0;		        //default to not finding the matching line
  int column1=0,column2=0;     	//column indices of the variable(s) of interest 
  int n;
  int colValue1=0,colValue2=0;  //for value scanned from column
  char *quote;                  //pointer for quote
  char temp[200];
  char fullVarName[200];        //name of the actual variable  
  char qVarName[200];           //name of the column header (in quotes)" 

  //figure out which column(s) in the table we're looking at
  for( column1=0;column1<fNVariables;column1++){
    if((strcmp(variable1,fVariableNames[column1]))==0){
      break;
    }	
  }
  if(column1>=fNVariables){
    fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable1,fSetupFileName );
    return -1;
  }
  if(fVariableTypes[column1]!=EHInt){
    fprintf(stderr,"Error, TExpTable: can't compare int (= %d) with non-int variable %s\n",value,variable1);
    return -1;    
  }
  //if variable2 non null (ie we're looking in range var1 - var2)
  if(variable2){
    for( column2=0;column2<fNVariables;column2++){
      if((strcmp(variable2,fVariableNames[column2]))==0){
	break;
      }	
    }
    if(column2>=fNVariables){
      fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable2,fSetupFileName );
      return -1;
    }
    if(fVariableTypes[column2]!=EHInt){
      fprintf(stderr,"Error, TExpTable: can't compare int (= %d) with non-int variable %s\n",value,variable2);
      return -1;    
    }
  }
  // open the file and loop over all the rows
  if((fp = fopen(fSetupFileName,"r")) == NULL){ 	//open setup file
    fprintf(stderr,"Error, TExpTable: couldn't open file: %s\n",fSetupFileName);
    return -1;
  }
  for(;;){
    if(fgets(line,256,fp) == NULL) break;     	        // check got a line from file
    if((line[0] == '*')||(line[0] == '#')) continue;	// ignore comment lines beginning #
    if((line[3] == '*')||(line[3] == '#')) continue;	// ignore comment lines beginning //# 
    if(sscanf(line,"%s",temp)!=1) continue;             //blank line or something
    if(!strstr(temp,"Row:")) continue;                  //not a Row:

    if((sscanf(line,fVariableFormats[column1],&colValue1))!=1){
      fprintf(stderr,"Error, TExpTable: couldn't scan some line in: %s\n",fSetupFileName);
      fclose(fp);
      return -1;
    }
    if(variable2){	//if var2 exists read its value too
      if((sscanf(line,fVariableFormats[column2],&colValue2))!=1){
	fprintf(stderr,"Error, TExpTable: couldn't scan some line in: %s\n",fSetupFileName);
	fclose(fp);
	return -1;
      }
    }
    //if var2 exists we're looking in the range between var1 and var2
    if(((!variable2)&&(value>=colValue1-range)&&(value<=colValue1+range))||((variable2)&&(value>=colValue1)&&(value<=colValue2))){
      match++;

      if(match>1){
	if(variable2)fprintf(stderr,"Error, TExpTable: found 2 or more lines in %s with %s <= %d <= %s:",
			     fSetupFileName, variable1,value,variable2);
	else fprintf(stderr,"Error, TExpTable: found 2 or more lines in %s with %s == %d:\n",
		     fSetupFileName, variable1,value);
	fprintf(stderr,"Here are the 1st 2 lines:\n%s%s\nUsing the 1st one.\n",goodline,line);
	return 0;
      }
      strcpy(goodline,line);
    }
  }
  fclose(fp);	
  
  if(!match){
    if(variable2)fprintf(stderr,"Error, TExpTable: couldn't find line in %s with %s <= %d <= %s\n",
			 fSetupFileName, variable1,value,variable2);
    else fprintf(stderr,"Error, TExpTable: couldn't find line in %s with %s == %d\n",
			 fSetupFileName, variable1,value);
    return -1;
  }
  fprintf(stdout,"\nTExpTable: Loaded variables from \"%s\":\nColumn\t          Variable\t\t    Value\n",fSetupFileName);
  for(n=0;n<fNVariables;n++){			//scan 1 variable name at a time
    if(!fVariableAddresses[n]) continue;	//scip if address is null

    sprintf(fullVarName,"%s_%s",fTableName,fVariableNames[n]);
    sprintf(qVarName,"\"%s\"",fVariableNames[n]);

    //switch according to the last letter in the format string
    ptr=strrchr(fVariableFormats[n],'%')+1;
    switch(ptr[0]){
    case('d'):	//integer
      sscanf(goodline,fVariableFormats[n],(int *)fVariableAddresses[n]);
      fprintf(stdout,"%-18s%-26s%d\n",qVarName,fullVarName,*(int *)fVariableAddresses[n]);
      break;
    case('f'):
       sscanf(goodline,fVariableFormats[n],(float *)fVariableAddresses[n]);
       fprintf(stdout,"%-18s%-26s%f\n",qVarName,fullVarName,*(float *)fVariableAddresses[n]);
      break;
    case('s'):
      sscanf(goodline,fVariableFormats[n],(char *)fVariableAddresses[n]);
      if(n==fNVariables-1){ //last column can be a comment string containing white space between quotes
	if((quote=strchr(goodline,'"'))!=NULL){ //hunt for "
	  if((ptr=strrchr(quote+1,'"'))!=NULL){ //hunt for 2nd quote
	    if((ptr-quote)<200){ //if short enough string
	      strncpy((char *)fVariableAddresses[n],quote+1,(ptr-quote)-1); //copy bit within quotes
	      strcpy(ptr,""); //and finish the string properly
	    }
	  }
	}
      }
      fprintf(stdout,"%-18s%-26s%s\n",qVarName,fullVarName,(char *)fVariableAddresses[n]);
      break;
    default:
      break;
    }
  }
  fprintf(stdout,"\n");	
  return 0;
}

int TExpTable::GetAllForValue(float value, char *variable1,  char *variable2, float range){
  //
  //Look for a line in the table where the float value equals the value of variable1
  //or is in the range between variable1 and variable2 (inclusive)
  //and set all variables to be their values from that line in the table 
  //return -1 if error
  //return 1 if only single line found (success)
  //return 0 and print warning if multiple lines found
  //If multiple lines found, use values from the first one found.

  FILE *fp;			//general file ptr
  Char_t line[256];		//for lines from the setup file
  Char_t goodline[256];		//for condition matched line from the setup file
  char *ptr;  			//General ptr
  int match=0;		        //default to not finding the matching line
  int column1=0,column2=0;     	//column indices of the variable(s) of interest 
  int n;
  float colValue1=0,colValue2=0;  //for value scanned from column
  char *quote;                  //pointer for quote
  char temp[200];
  float limit=range;            //set the limit of the comparison to range
  char fullVarName[200];        //name of the actual variable  
  char qVarName[200];           //name of the column header (in quotes)" 

  //figure out which column(s) in the table we're looking at
  for( column1=0;column1<fNVariables;column1++){
    if((strcmp(variable1,fVariableNames[column1]))==0){
      break;
    }	
  }
  if(column1>=fNVariables){
    fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable1,fSetupFileName );
    return -1;
  }
  if(fVariableTypes[column1]!=EHFloat){
    fprintf(stderr,"Error, TExpTable: can't compare float (= %f) with non-float variable %s\n",value,variable1);
    return -1;    
  }
  //if variable2 non null (ie we're looking in range var1 - var2)
  if(variable2){
    for( column2=0;column2<fNVariables;column2++){
      if((strcmp(variable2,fVariableNames[column2]))==0){
	break;
      }	
    }
    if(column2>=fNVariables){
      fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable2,fSetupFileName );
      return -1;
    }
    if(fVariableTypes[column2]!=EHFloat){
      fprintf(stderr,"Error, TExpTable: can't compare int (= %f) with non-float variable %s\n",value,variable2);
      return -1;    
    }
  }
  //if no variable 2 we're going to allow some tolerance for searching +- range
  //but if range = <0 (default) well allow +- 1/1000 of the value
  else{ 
    if(range<0.0) limit=value/1000.0;
  }

  // open the file and loop over all the rows
  if((fp = fopen(fSetupFileName,"r")) == NULL){ 	//open setup file
    fprintf(stderr,"Error, TExpTable: couldn't open file: %s\n",fSetupFileName);
    return -1;
  }
  for(;;){
    if(fgets(line,256,fp) == NULL) break;     	        // check got a line from file
    if((line[0] == '*')||(line[0] == '#')) continue;	// ignore comment lines beginning #
    if((line[3] == '*')||(line[3] == '#')) continue;	// ignore comment lines beginning //# 
    if(sscanf(line,"%s",temp)!=1) continue;             //blank line or something
    if(!strstr(temp,"Row:")) continue;                  //not a Row:

    if((sscanf(line,fVariableFormats[column1],&colValue1))!=1){
      fprintf(stderr,"Error, TExpTable: couldn't scan some line in: %s\n",fSetupFileName);
      fclose(fp);
      return -1;
    }
    if(variable2){	//if var2 exists read its value too
      if((sscanf(line,fVariableFormats[column2],&colValue2))!=1){
	fprintf(stderr,"Error, TExpTable: couldn't scan some line in: %s\n",fSetupFileName);
	fclose(fp);
	return -1;
      }
    }
    //if var2 exists we're looking in the range between var1 and var2
    //if not we look at var1+-limit
    if(((!variable2)&&(value>=colValue1-limit)&&(value<=colValue1+limit))||((variable2)&&(value>=colValue1)&&(value<=colValue2))){
      match++;

      if(match>1){
	if(variable2)fprintf(stderr,"Error, TExpTable: found 2 or more lines in %s with %s <= %f <= %s:",
			     fSetupFileName, variable1,value,variable2);
	else fprintf(stderr,"Error, TExpTable: found 2 or more lines in %s with %s == %f:\n",
		     fSetupFileName, variable1,value);
	fprintf(stderr,"Here are the 1st 2 lines:\n%s%s\nUsing the 1st one.\n",goodline,line);
	return 0;
      }
      strcpy(goodline,line);
    }
  }
  fclose(fp);	
  
  if(!match){
    if(variable2)fprintf(stderr,"Error, TExpTable: couldn't find line in %s with %s <= %f <= %s\n",
			 fSetupFileName, variable1,value,variable2);
    else fprintf(stderr,"Error, TExpTable: couldn't find line in %s with %s == %f\n",
			 fSetupFileName, variable1,value);
    return -1;
  }
  fprintf(stdout,"\nTExpTable: Loaded variables from \"%s\":\nColumn\t          Variable\t      Value\n",fSetupFileName);
  for(n=0;n<fNVariables;n++){                   //scan 1 variable name at a time
    if(!fVariableAddresses[n]) continue;        //scip if address is null

    sprintf(fullVarName,"%s_%s",fTableName,fVariableNames[n]);
    sprintf(qVarName,"\"%s\"",fVariableNames[n]);

    //switch according to the last letter in the format string
    ptr=strrchr(fVariableFormats[n],'%')+1;
    switch(ptr[0]){
    case('d'):  //integer
      sscanf(goodline,fVariableFormats[n],(int *)fVariableAddresses[n]);
      fprintf(stdout,"%-18s%-20s%d\n",qVarName,fullVarName,*(int *)fVariableAddresses[n]);
      break;
    case('f'):
      sscanf(goodline,fVariableFormats[n],(float *)fVariableAddresses[n]);
      fprintf(stdout,"%-18s%-20s%f\n",qVarName,fullVarName,*(float *)fVariableAddresses[n]);
      break;
    case('s'):
      sscanf(goodline,fVariableFormats[n],(char *)fVariableAddresses[n]);
      if(n==fNVariables-1){ //last column can be a comment string containing white space between quotes
        if((quote=strchr(goodline,'"'))!=NULL){ //hunt for "
          if((ptr=strrchr(quote+1,'"'))!=NULL){ //hunt for 2nd quote
            if((ptr-quote)<200){ //if short enough string
              strncpy((char *)fVariableAddresses[n],quote+1,(ptr-quote)-1); //copy bit within quotes
              strcpy(ptr,""); //and finish the string properly
            }
          }
        }
      }
      fprintf(stdout,"%-18s%-20s%s\n",qVarName,fullVarName,(char *)fVariableAddresses[n]);
      break;
    default:
      break;
    }
  }
  fprintf(stdout,"\n");	
  return 0;
}

int TExpTable::GetAllForValue(char *value, char *variable1, int perfectmatch){
  //
  //Look for a line in the table where string value equals the variable string
  //If perfectmatch=0 variable must contain value as a substring
  //If perfectmatch=1 variable must equal value
  //return -1 if error
  //return 1 if only single line found (success)
  //return 0 and print warning if no lines or multiple lines found
  //If multiple lines found, use values from the first one found.
  //

  FILE *fp;			//general file ptr
  Char_t line[256];		//for lines from the setup file
  Char_t goodline[256];		//for condition matched line from the setup file
  char *ptr;  			//General ptr
  int match=0;		        //default to not finding the matching line
  int column1=0;           	//column indices of the variable(s) of interest 
  int n;
  char colValue1[200];           //for value scanned from column
  char *quote;                  //pointer for quote
  char temp[200];
  char fullVarName[200];        //name of the actual variable  
  char qVarName[200];           //name of the column header (in quotes)" 


  //figure out which column(s) in the table we're looking at
  for( column1=0;column1<fNVariables;column1++){
    if((strcmp(variable1,fVariableNames[column1]))==0){
      break;
    }	
  }
  if(column1>=fNVariables){
    fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable1,fSetupFileName );
    return -1;
  }
  if(fVariableTypes[column1]!=EHString){
    fprintf(stderr,"Error, TExpTable: can't compare string (= %s) with non-string variable %s\n",value,variable1);
    return -1;    
  }

  // open the file and loop over all the rows
  if((fp = fopen(fSetupFileName,"r")) == NULL){ 	//open setup file
    fprintf(stderr,"Error, TExpTable: couldn't open file: %s\n",fSetupFileName);
    return -1;
  }
  for(;;){
    if(fgets(line,256,fp) == NULL) break;     	        // check got a line from file
    if((line[0] == '*')||(line[0] == '#')) continue;	// ignore comment lines beginning #
    if((line[3] == '*')||(line[3] == '#')) continue;	// ignore comment lines beginning //# 
    if(sscanf(line,"%s",temp)!=1) continue;             //blank line or something
    if(!strstr(temp,"Row:")) continue;                  //not a Row:

    if((sscanf(line,fVariableFormats[column1],&colValue1))!=1){
      fprintf(stderr,"Error, TExpTable: couldn't scan some line in: %s\n",fSetupFileName);
      fclose(fp);
      return -1;
    }
    if(((!perfectmatch)&&(strstr(colValue1,value)))||((perfectmatch)&&(!strcmp(colValue1,value)))){
      match++;

      if(match>1){
	fprintf(stderr,"Error, TExpTable: found 2 or more lines in %s with %s == %s:\n",
		fSetupFileName, variable1,value);
	fprintf(stderr,"Here are the 2 lines:\n%s%s\nUsing the 1st one.\n",goodline,line);
	return 0;
      }
      strcpy(goodline,line);
    }
  }
  fclose(fp);	
  
  if(!match){
    fprintf(stderr,"Error, TExpTable: couldn't find line in %s with %s == %s\n",
	    fSetupFileName, variable1,value);
    return 0;
  }

  fprintf(stdout,"\nTExpTable: Loaded variables from \"%s\":\nColumn\t          Variable\t      Value\n",fSetupFileName);
  for(n=0;n<fNVariables;n++){                   //scan 1 variable name at a time
    if(!fVariableAddresses[n]) continue;        //scip if address is null

    sprintf(fullVarName,"%s_%s",fTableName,fVariableNames[n]);
    sprintf(qVarName,"\"%s\"",fVariableNames[n]);

    //switch according to the last letter in the format string
    ptr=strrchr(fVariableFormats[n],'%')+1;
    switch(ptr[0]){
    case('d'):  //integer
      sscanf(goodline,fVariableFormats[n],(int *)fVariableAddresses[n]);
      fprintf(stdout,"%-18s%-20s%d\n",qVarName,fullVarName,*(int *)fVariableAddresses[n]);
      break;
    case('f'):
      sscanf(goodline,fVariableFormats[n],(float *)fVariableAddresses[n]);
      fprintf(stdout,"%-18s%-20s%f\n",qVarName,fullVarName,*(float *)fVariableAddresses[n]);
      break;
    case('s'):
      sscanf(goodline,fVariableFormats[n],(char *)fVariableAddresses[n]);
      if(n==fNVariables-1){ //last column can be a comment string containing white space between quotes
        if((quote=strchr(goodline,'"'))!=NULL){ //hunt for "
          if((ptr=strrchr(quote+1,'"'))!=NULL){ //hunt for 2nd quote
            if((ptr-quote)<200){ //if short enough string
              strncpy((char *)fVariableAddresses[n],quote+1,(ptr-quote)-1); //copy bit within quotes
              strcpy(ptr,""); //and finish the string properly
            }
          }
        }
      }
      fprintf(stdout,"%-18s%-20s%s\n",qVarName,fullVarName,(char *)fVariableAddresses[n]);
      break;
    default:
      break;
    }
  }
  fprintf(stdout,"\n");	
  return 0;
}

void *TExpTable::GetVarAddress(char *variable){
  int column;
  //figure out which column(s) in the table we're looking at
  fprintf(stdout, "Getting address \n");
  for( column=0;column<fNVariables;column++){
    if((strcmp(variable,fVariableNames[column]))==0){
      break;
    }	
  }
  if(column>=fNVariables){
    fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable,fSetupFileName );
    return NULL;
  }
  
  else return(fVariableAddresses[column]);
}
  
void *TExpTable::SetVarAddress(char *variable, void *address){
  int column;
  //figure out which column(s) in the table we're looking at
  for( column=0;column<fNVariables;column++){
    if((strcmp(variable,fVariableNames[column]))==0){
      break;
    }	
  }
  if(column>=fNVariables){
    fprintf(stderr,"Error, TExpTable: couldn't find column \"%s\" in %s\n ",variable,fSetupFileName );
    return NULL;
  }
  
  
  else{
    fVariableAddresses[column]=address;
    return NULL;
  }
}
  
TExpTable::~TExpTable(){
  for(int n=0;n<fNVariables;n++){
    delete fVariableNames[n];
    delete fVariableFormats[n];
  }
  delete fVariableNames;
  delete fVariableAddresses;
  delete fVariableFormats;
}
