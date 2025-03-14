#ifndef __g14b_run_h__
#define __g14b_run_h__

#include "TExpTable.h"

// This header was created from the table g14b_allruns.dat on Fri Jun  9 11:19:43 EDT 2017 by running the command:
// table2vars g14b_allruns.dat g14b_run

// To use it, simple include the header
// A globally available ExpTable class - g14b_runTable is created
// together with global variables relating to each column (g14b_run_var) - defined below
//
// Calling the member function GetAllForValue(int value, char *variable1, char *variable2=NULL)
// searches for a row in the table and updates all the global variables to values in that row
// Usually it would be a Run number or a run range
// Eg. g14b_runTable->GetAllForValue(53772,"Run"); would find the row where Run=53772
// and update all the variables to the values from that row.
// Eg. g14b_runTable->GetAllForValue(53772,"RunMin","RunMax"); would find the row where RunMin <= 53772 <= RunMax
// and update all the variables to the values from that row.


// Create a TExpTable with this file as the actual data table
TExpTable *g14b_runTable = new TExpTable("/u/home/clasg14/rootbeer2.1/dat/g14b_allruns.dat","g14b_run");

// Define variables for all the columns in the run table and set their addresses in g14b_runTable 

// An obscure way to define and set global table and variables that works in both
// interpreted and compiled mode - void *axx are dummies


Int_t   g14b_run_Run;                   void *a2  = g14b_runTable->SetVarAddress("Run",&g14b_run_Run);
Int_t   g14b_run_Nfiles;                void *a3  = g14b_runTable->SetVarAddress("Nfiles",&g14b_run_Nfiles);
Int_t   g14b_run_Nevents;               void *a4  = g14b_runTable->SetVarAddress("Nevents",&g14b_run_Nevents);
Char_t  g14b_run_Period[200];           void *a5  = g14b_runTable->SetVarAddress("Period",g14b_run_Period);
Char_t  g14b_run_Tstart[200];           void *a6  = g14b_runTable->SetVarAddress("Tstart",g14b_run_Tstart);
Float_t g14b_run_Ebeam;                 void *a7  = g14b_runTable->SetVarAddress("Ebeam",&g14b_run_Ebeam);
Float_t g14b_run_Itorus;                void *a8  = g14b_runTable->SetVarAddress("Itorus",&g14b_run_Itorus);
Float_t g14b_run_Tlive;                 void *a9  = g14b_runTable->SetVarAddress("Tlive",&g14b_run_Tlive);
Float_t g14b_run_Pbeam;                 void *a10 = g14b_runTable->SetVarAddress("Pbeam",&g14b_run_Pbeam);
Float_t g14b_run_WienPhi;               void *a11 = g14b_runTable->SetVarAddress("WienPhi",&g14b_run_WienPhi);
Float_t g14b_run_WienH;                 void *a12 = g14b_runTable->SetVarAddress("WienH",&g14b_run_WienH);
Int_t   g14b_run_HWP;                   void *a13 = g14b_runTable->SetVarAddress("HWP",&g14b_run_HWP);
Float_t g14b_run_CohEdge;               void *a14 = g14b_runTable->SetVarAddress("CohEdge",&g14b_run_CohEdge);
Int_t   g14b_run_CohPlane;              void *a15 = g14b_runTable->SetVarAddress("CohPlane",&g14b_run_CohPlane);
Int_t   g14b_run_Radiator;              void *a16 = g14b_runTable->SetVarAddress("Radiator",&g14b_run_Radiator);
Float_t g14b_run_Hpol;                  void *a17 = g14b_runTable->SetVarAddress("Hpol",&g14b_run_Hpol);
Float_t g14b_run_Dpol;                  void *a18 = g14b_runTable->SetVarAddress("Dpol",&g14b_run_Dpol);
Char_t  g14b_run_Comment[200];          void *a19 = g14b_runTable->SetVarAddress("Comment",g14b_run_Comment);

#endif
