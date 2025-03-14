#ifndef ROOTBEER_BOS_TYPES
#define ROOTBEER_BOS_TYPES

#define WORD			4      	// a word is 4 bytes on 32bit systems 
#define MAX_BYTES_PER_BUFF_READ 2500000	// size of raw data buffers 
#define MAX_STORED_EVENTS	1000 	// max no events per buff - size of ring buffer 

#define SPLIT_BANK_SPACE	20000	// space at end of data buff for event split between buffs 

#define SEG_COMPLETE     	0	// Segment contains full event 
#define SEG_FIRST		1	// Segment contains first part of event 
#define SEG_MIDDLE		2	// Segment contains a middle chunk of the event 
#define SEG_LAST		3	// Segment contains last part of the event 

#define NEW_EVENT		0	
#define OLD_EVENT		1
	      			
#define BUFFER_SCANNED		0
#define BUFFER_ANALYSED		1
#define BUFFER_READ_RAW		2
#define BUFFER_LAST		3

#define EVENT_ANALYSED		0
#define EVENT_SCANNED		1
#define EVENT_LAST		2

#define EVENT_LAST_IN_BUFFER	100


#define NO_SWAP			0
#define SWAP			1


// Following types are structs to map the various headers in an FPACK record
//   See the FPACK manual - Appendix A. For description of FPACK records format.
//   Only the fields used in the code are commented. 
//
//   Briefly, 
//   A Record contains a Header followed by 1 or more Record Segments.
//   ..A Record Segment contains a Record Segment Header followed by 1 or more Data Segments
//   ....A Data Segment contains a Data Segment Header followed by Data.
//   ....The Data Segment contains stuff relating to one BANK:
//   ....The Data Segment Header describes the BANK, and the data format
//   ....The data follows.
//   ....To make things difficult, a BANK can be split between records! 

typedef struct RecHead_t{	// Record Header - every record has this 

  int	nrecwd;			// No of words in the record (same for all records) 
  int	nusedw;			// No of words used for data 
} RecHead_t;

typedef struct recSegHead_t{	// The first Record Segment Header comes immediately after the 
				//   Record Header, Others may or may not follow 
  int swapping;			
  int format;			
  int record_no;
  char name[2*WORD];		
  int number1;
  int number2;
  int physical_no;
  int dataword;
  int segment_code;		// complete event,or first, middle or last rec segment 
  int nsegwd;			// No of words in the rec segment 
} recSegHead_t;

typedef struct dataSegHead_t{	// One of these always follows a Rec Seg Header,
				//   and others follow 
  int nwords;			// No of words in this data seg Header  
  char name[8];			// Name of BANK 
  int number;			// Number of BANK 
  int ncol;			// Number of columns in bank 
  int nrow;			// Number of rowsin BANK 
  int segment_code;		// complete BANK,or first, middle or last data segment 
  int nsegwd_prev;		// > 0 if second (or more) part of a split BANK  
  int nsegwd_this;		// no of datawords 
  char format[20];		// Describes the FPACK format of the data, immediately following 
} dataSegHead_t;


// Now some structs to hold bank data 


typedef struct bankTab_t {     		// For each instance of one BANKtype in an event, this holds 
					//   the BANK number and the address of the data segment header 
  int number[MAX_BANK_NUMBER+2];	// Number of the BANK 
  int nrow[MAX_BANK_NUMBER];		// number of rows in the bank 
  void *address[MAX_BANK_NUMBER];	// Address of the data 
} bankTab_t;



typedef struct eventInf_t {		// The information relating to a single event 

  int buffer;				// which buffer holds the data (0 or 1 ) 
  int status;      			// status not available, available, processed 
  int nBanks;				// number of banks in event 
  int dummy;				// for 8 byte alignment 
  struct bankTab_t 	BTab[MAX_BANK_TYPES]; 	// One table for each bank type 
  struct dataSegHead_t *BHead[MAX_BANK_TYPES];	// Addresses of all bank data seg headers 
  void  	       *BData[MAX_BANK_TYPES];	// Addresses of all bank data 
} eventInf_t;

#endif
