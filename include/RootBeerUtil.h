#ifndef ROOTBEERUTIL
#define ROOTBEERUTIL
#include "bankvars.h" // all the BANK_t structs and related global variables

#define OFF	0			//defs to set the status of banks
#define ON	1

TRootBeer *createBeerObject(char *);
int getFileMode(char *);
int getNextFile(char *, char *);
int resetNextFile();
int getEpics(char *, float *);
Int_t getBit(int word,int bit);
#endif
