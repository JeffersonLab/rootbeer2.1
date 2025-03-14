#ifndef __G14POL_H__
#define __G14POL_H__
//some sample functions to load polarization tables and get the polarization from them

//enumerator for the column names in the pol tables
enum {
  E_ID    = 0,
  ENERGY  = 1,
  ENH     = 2,
  ENHERR  = 3,
  ENHFIT  = 4,
  PFIT    = 5,
  PCOR    = 6,
  PCORERR = 7,
  PSMOOTH = 8,
  EDGE    = 9
};

enum {
  PARA=0,
  PERP=1
};

Double_t polTable[2][1000][385][10];  //where its [plane][edge][E_id][field]
Int_t polTableN[2]={0,0};            //No of entries for para and perp
Char_t polFirstLines[2][1000][250];   //to keep the 1st lines if needed (info from files)

Int_t edgeEventLow[10000];            //hold the current table of edge positions for event ranges
Int_t edgeEventHigh[10000];
Double_t edgeEventEdge[10000];
Int_t edgeEventPlane[10000];
Int_t edgeEventN;
Int_t edgeIndex=0;
Int_t lastEdgeEvent=0;
Double_t lastCohEdge=0.0;
Int_t lastCohPlane=-1;



int LoadPolTable(int plane, Char_t *PolTableList){
  polTableN[plane]=0;//added by H.Lu to enable reload table list;
  FILE *fplist,*fpfile;              //file pointers for filelist and each file therein
  Char_t lline[250];                 //for reading in lines from filelist
  Char_t fline[250];                 //for reading in lines from file
  Char_t filename[250];              //file
  Int_t  fcount=0;                   //counter for no of files read in
  Int_t  chancount=0;                //counter for no of channels read in
  Int_t  eid=0;
  Double_t edge=0.0;                 
  
  if((fplist=fopen(PolTableList,"r"))==NULL){ //open filelist
    cerr << "Error Couldn't open file: " << PolTableList << endl;
    return -1;
  }

  fcount=0; 
  //for each file in the list
  while(fgets(lline,240,fplist) != NULL){
    if((lline[0] == '*')||(lline[0] == '#')) continue; //skip comments
    sscanf(lline,"%s",filename);                       //read in filename
     
    cout << "opening " << filename << "   " << endl;
    if((fpfile=fopen(filename,"r"))==NULL){              //open file
      cerr << "Error Couldn't open file: " << filename << endl;
      return -1;
    }
        
    fgets(polFirstLines[plane][polTableN[plane]],240,fpfile ); //save the 1st line

    //scan the bit after the last "_" in the filename to get the edge energy
    sscanf(strrchr(filename,'_')+1,"%lg",&polTable[plane][fcount][0][EDGE]);

    chancount=0;                                             //starting array at 1 for consistency with E_ID
    while((fgets(fline,240,fpfile)) != NULL){
      if((fline[0] == '*')||(fline[0] == '#')) continue;     //skip comments    
      sscanf(fline,"%d",&eid);                               //first get the E_ID
      sscanf(fline,"%*d%lg%lg%lg%lg%lg%lg%lg%lg",
	     &polTable[plane][fcount][eid][ENERGY],
	     &polTable[plane][fcount][eid][ENH],
	     &polTable[plane][fcount][eid][ENHERR],
	     &polTable[plane][fcount][eid][ENHFIT],
	     &polTable[plane][fcount][eid][PFIT],
	     &polTable[plane][fcount][eid][PCOR],
	     &polTable[plane][fcount][eid][PCORERR],
	     &polTable[plane][fcount][eid][PSMOOTH]);
      chancount++; 
    }
    fclose(fpfile); //close the file
    if(chancount!=384){
      cerr << "Should be 384 lines in " << filename << " - only got " << chancount << endl;
      return -1;
    }
    polTableN[plane]++;
    
    fcount++;
  }

  fclose(fplist);

  return(0);
}


Double_t GetPol(Int_t plane, Double_t edge, Int_t eid, Int_t poltype = PSMOOTH, Double_t lowThresh=0.2, Double_t highThresh=0.3){
  //get polarization based on eid and edge position

  Int_t eIndex=0;
  Double_t pol=-1.0;
  
  //Check edge in in range of tables
  if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;
  
  //find index
  for(eIndex=0;eIndex<polTableN[plane];eIndex++){
    if(polTable[plane][eIndex][0][EDGE]>=edge) break;
  }
  pol=polTable[plane][eIndex][eid][poltype];
  /*
  if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
  if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;
  */
  return pol;
}
Double_t GetPol(Int_t plane, Double_t edge, Double_t energy, Int_t poltype = PSMOOTH, Double_t lowThresh=0.2, Double_t highThresh=0.3){
  //get polarization based on eid and edge position

  Int_t eIndex=0;
  Double_t pol=-1.0;
  
  //Check edge in in range of tables
  cout<<"checking edge in range:"<<edge<<"?"<<polTable[plane][1][0][EDGE]<<","<<polTable[plane][polTableN[plane]-1][0][EDGE]<<endl;
  if((edge<polTable[plane][1][0][EDGE])||(edge>polTable[plane][polTableN[plane]-1][0][EDGE])) return -1.0;
  
  cout << "In range" << endl;

  //find index
  for(eIndex=0;eIndex<polTableN[plane];eIndex++){
    if(polTable[plane][eIndex][0][EDGE]>=edge) break;
  }
  cout << "Index = " << eIndex << endl;
  int eid;
  for(eid=1;eid<=384;eid++){
    cout<<polTable[plane][eIndex][eid][ENERGY]<<endl;
  }
  cout<<"eid="<<eid<<endl;
  pol=polTable[plane][eIndex][eid][poltype];
  /*
  if((polTable[plane][eIndex][0][ENERGY]<edge)&&(pol<lowThresh)) pol = -1.0;
  if((polTable[plane][eIndex][0][ENERGY]>edge)&&(pol<highThresh)) pol = -1.0;
  */
  return pol;
}

int LoadEdgeTable(Char_t *EdgeTable){
  FILE *fpfile;                      //file pointers for filelist and each file therein
  Char_t fline[250];                 //for reading in lines from file

  cout << "opening " << EdgeTable << "   " << endl;
  if((fpfile=fopen(EdgeTable,"r"))==NULL){              //open file
    cerr << "Error Couldn't open file: " << EdgeTable << endl;
    return -1;
  }


  edgeEventN=0;                      //initialize the counter
  edgeIndex=0;                       //and index for current table
  lastEdgeEvent=0;                   //etc
  lastCohEdge=0.0;
  lastCohPlane=-1;
  
  while((fgets(fline,240,fpfile)) != NULL){
    if((fline[0] == '*')||(fline[0] == '#')) continue;     //skip comments
    sscanf(fline,"%d%d%lg%d",        //scan in the tables
	   &edgeEventLow[edgeEventN],
	   &edgeEventHigh[edgeEventN],
	   &edgeEventEdge[edgeEventN],
	   &edgeEventPlane[edgeEventN]);
    edgeEventN++;                    //inc the counter fo the no of lines in the table
  }
  fclose(fpfile);
  return 0;
}

Int_t GetEdge(Int_t event, Double_t *edge, Int_t *plane){

  //check event no >= previous event
  
  if(event<lastEdgeEvent){
    //    cerr << "Warning: event (= " << event << ") is earlier than previous event (= " <<  lastEdgeEvent << ") - ignoring" << endl;
    return -1;
  }
  lastEdgeEvent=event;

  //keep going until we're in range or don't find a range 
  while((!((event >= edgeEventLow[edgeIndex])&&(event <= edgeEventHigh[edgeIndex])))&&(edgeIndex<edgeEventN)) edgeIndex++;
  if(edgeIndex>=edgeEventN){
    //    cerr << "Error: event (= " << event << ") is not in the range of this table" << endl;    
    edgeIndex=0; //reset to zero
    return -1;
  }
  *edge = edgeEventEdge[edgeIndex];
  *plane = edgeEventPlane[edgeIndex];
  return 0;
}
#endif
