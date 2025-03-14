#include <fstream>
#include <cmath>
#include "g14pcor.h"
#include <stdlib.h>
#include <algorithm>

char files_location[]="/home/clasg14/rootbeer2.1/g14pcor";


extern "C" { // 'dumb it down' for c and fortran
//_____________________________________________________________________________
  // Tagger Correction:
  /* C/C++ Function:
   * Parameters: eg         Measured photon energy
   *             e_beam     Beam energy
   *             eid        E-paddle id in the tagger
   *             runNumber  CLAS run number
   *
   * Returns: Corrected photon energy
   */
  double g14TaggerCor(double __eg,double __e_beam,int __eid);

  /* FORTRAN Subroutine:
   * Parameters: eg         Measured photon energy
   *             e_beam     Beam energy
   *             eid        E-paddle id in the tagger
   *             runNumber  CLAS run number
   *             egCor      Corrected photon energy
   */
  void g14TaggerCor_(double *__eg,double *__e_beam,int *__eid,
		     double *__egCor);
//_____________________________________________________________________________
  // Momentum Corrections:
  /* C/C++ Function:
   * Parameters: p3[3]    Measured 3-momentum [x,y,z]
   *             q        Charge (ex. -1 for pi-)
   *             p3cor[3] Corrected 3-momentum [x,y,z]
   */
  void g14MomCor(double __p3[3],int __q, int __runNo,double __p3cor[3]);

  /* FORTRAN Subroutine:
   * Parameters: p3[3]    Measured 3-momentum [x,y,z]
   *             q        Charge (ex. -1 for pi-)
   *             p3cor[3] Corrected 3-momentum [x,y,z]
   */
  void g14MomCor_(double __p3[3],int *__q,int *__runNo, double __p3cor[3]);

}

//_____________________________________________________________________________

// Function prototypes:
int g14GetSector(const double __p3[3]);

int g14GetBin(double __theta,double __phi,const double __phiBins_min[12],
	   const double __phiBins_max[12],const double __thetaBins_min[15],
	      const double __thetaBins_max[15]);			

void g14ReadBinsFile(double __phiBins_min[12],double __phiBins_max[12],
		  double __thetaBins_min[15],double __thetaBins_max[15]);

double g14GetLambdaTrack(const double __p3[3]);

double g14GetPhiTrack(const double __p3[3]);

void g14SetP3FromTrackingParameters(double __p,double __lambda,double __phi,
				 int __sector,double __p3[3]);

double g14GetSectorPhi(double __phi);

void Init_g14pcor(double __p_vals[2][6][180][4],int __num_pars_p[2][6][180],
		  double __phi_vals[2][6][180][4],
		  int __num_pars_phi[2][6][180],
		  double __lambda_vals[2][6][180][4],
		  int __num_pars_lambda[2][6][180] );

void Init_linear_pcor(double __p_vals[2][6][180][4],int __num_pars_p[2][6][180],
		  double __phi_vals[2][6][180][4],
		  int __num_pars_phi[2][6][180],
		  double __lambda_vals[2][6][180][4],
		  int __num_pars_lambda[2][6][180] );

//_____________________________________________________________________________

// FORTRAN wrappers:

void g14TaggerCor_(double *__eg,double *__e_beam,int *__eid,
		   double *__egCor){
  *__egCor = g14TaggerCor(*__eg,*__e_beam,*__eid);
}

void g14MomCor_(double __p3[3],int *__q, int *__runNo, double __p3cor[3]){
  g14MomCor(__p3,*__q,*__runNo,__p3cor);
}


//_____________________________________________________________________________

/// Tagger correction for the g14 run period
/** @param eg         Measured photon energy (GeV)
 *  @param e_beam     Beam energy (GeV)
 *  @param eid        E-paddle id in the tagger
 *  @param runNumber  CLAS run number
 *
 *  <b>Returns</b>    Corrected photon energy
 *
 *  Corrects the tagger measurement using parameters from tagger_cor.dat (to 
 *  find this file, the user's CLAS_PACK enviornment variable must be set).
 *  Next a run dependent beam energy correction is applied. The corrected
 *  energy is returned.
 *
 *  Note: If the correction parameter file can not be opened (most likely due
 *        to inproper enviornment set up), then -666 is returned.
 *
 */
double g14TaggerCor(double __eg,double __e_beam,int __eid){

  static bool first_call = true;
  static double eid_cor[767][4];
  //  static double offset[700];

  if(first_call){ // only do this on the 1st call to the function
    // read in the correction parameters
    char file[200];
    sprintf(file,"%s/tagger_cor_2281.dat",files_location);
    std::ifstream inFile_egcor(file);
    // check to make sure the file opened properly
    if(!inFile_egcor.is_open()){
      printf("Error!!!! <g14TaggerCor> Can NOT open tagger_cor_2281.dat. This ");
      printf("file should be located in %s. Check ",files_location);
      printf("to make sure your CLAS_PACK enviornmen variable is properly ");
      printf("set. \n");
      return -666.; // return nonsense
    }
      
    int eid,count = 0;
    double cor;
    while(inFile_egcor >> eid){
      inFile_egcor >> cor;
      eid_cor[eid-1][0] = cor;
      count++;
    }    
    inFile_egcor.close(); // close the input file 
    if(count != 767){
      printf("Error!!!! <g14TaggerCor> Read incorrect number of E-paddles ");
      printf("[%d instead of 767]. \n",count);
      return -666.;
    }
    
    
     sprintf(file,"%s/tagger_cor_2258.dat",files_location);
    inFile_egcor.open(file);
    // check to make sure the file opened properly
    if(!inFile_egcor.is_open()){
      printf("Error!!!! <g14TaggerCor> Can NOT open tagger_cor_2258.dat. This ");
      printf("file should be located in %s. Check ",files_location);
      printf("to make sure your CLAS_PACK enviornmen variable is properly ");
      printf("set. \n");
      return -666.; // return nonsense
    }
      
    eid=0;
    count = 0;
    
    while(inFile_egcor >> eid){
      inFile_egcor >> cor;
      eid_cor[eid-1][1] = cor;
      count++;
    }    
    inFile_egcor.close(); // close the input file 
    if(count != 767){
      printf("Error!!!! <g14TaggerCor> Read incorrect number of E-paddles ");
      printf("[%d instead of 767]. \n",count);
      return -666.;
    }
    
    
     sprintf(file,"%s/tagger_cor_2541.dat",files_location);
    inFile_egcor.open(file);
    // check to make sure the file opened properly
    if(!inFile_egcor.is_open()){
      printf("Error!!!! <g14TaggerCor> Can NOT open tagger_cor_2541.dat. This ");
      printf("file should be located in %s. Check ",files_location);
      printf("to make sure your CLAS_PACK enviornmen variable is properly ");
      printf("set. \n");
      return -666.; // return nonsense
    }
   
    eid=0;
    count = 0;
    
    while(inFile_egcor >> eid){
      inFile_egcor >> cor;
      eid_cor[eid-1][2] = cor;
      count++;
    }    
    inFile_egcor.close(); // close the input file 
    if(count != 767){
      printf("Error!!!! <g14TaggerCor> Read incorrect number of E-paddles ");
      printf("[%d instead of 767]. \n",count);
      return -666.;
    }
    
     sprintf(file,"%s/linear_correction/tagger_cor_5552.dat",files_location);
    inFile_egcor.open(file);
    // check to make sure the file opened properly
    if(!inFile_egcor.is_open()){
      printf("Error!!!! <g14TaggerCor> Can NOT open tagger_cor_5552.dat. This ");
      printf("file should be located in %s/linear_correction. Check ",files_location);
      printf("to make sure your CLAS_PACK enviornmen variable is properly ");
      printf("set. \n");
      return -666.; // return nonsense
    }
      
     eid=0;
    count = 0;
    
    while(inFile_egcor >> eid){
      inFile_egcor >> cor;
      eid_cor[eid-1][3] = cor;
      count++;
    }    
    inFile_egcor.close(); // close the input file 
    if(count != 767){
      printf("Error!!!! <g14TaggerCor> Read incorrect number of E-paddles ");
      printf("[%d instead of 767]. \n",count);
      return -666.;
    }

     first_call = false;
  }

  // check that eid is valid
  if(__eid < 1 || __eid > 767){
    printf("Warning <g14TaggerCor> E-paddle out of range [id = %d].",__eid);
    printf(" No correction will be applied. \n");
    return __eg;
  }

  double dEg=0;
  if(3.0<__e_beam && __e_beam<3.4) dEg= __e_beam*eid_cor[__eid-1][2]; //use gold2 correction for empty target
  if(2.260<__e_beam && __e_beam<2.290) dEg= __e_beam*eid_cor[__eid-1][0]; // tagger correction E_beam =2281
  if(__e_beam<2.260) dEg= __e_beam*eid_cor[__eid-1][1]; // tagger correction E_beam =2258 
  if(__e_beam>2.535 && __e_beam<3.000) dEg= __e_beam*eid_cor[__eid-1][2]; // tagger correction E_beam =2541
  if(__e_beam>5.550) dEg= __e_beam*eid_cor[__eid-1][3]; // tagger correction E_beam =5552
  
  
  return __eg + dEg; // return the corrected energy
  
  }
//_____________________________________________________________________________

/// Momentum corrections for g14.
/** @param p3  Measured 3-momentum (x,y,z)
 *  @param q   Charge in units of @a e (ex. -1 for pi-)
 *  @param p3cor Corrected 3-momentum
 *
 *  Sets p3cor to be the corrected momentum.
 */
void g14MomCor(double __p3[3],int __q,int __runNo, double __p3cor[3]){

  static bool firstCall = true;
  static double phiBins_min[12],phiBins_max[12];
  static double thetaBins_min[15],thetaBins_max[15];

 
 static double p_vals[2][6][180][4];
  static double phi_vals[2][6][180][4];
  static double lambda_vals[2][6][180][4];
  static int num_pars_p[2][6][180];
  static int num_pars_lambda[2][6][180];
  static int num_pars_phi[2][6][180];
 
  static double p_vals_lin[2][6][180][4];
  static double phi_vals_lin[2][6][180][4];
  static double lambda_vals_lin[2][6][180][4];
  static int num_pars_p_lin[2][6][180];
  static int num_pars_lambda_lin[2][6][180];
  static int num_pars_phi_lin[2][6][180];

  int bin;

  if(firstCall){ // on first call, do the following...
    // read in binning info
    g14ReadBinsFile(phiBins_min,phiBins_max,thetaBins_min,thetaBins_max);
    
    // read in momcor parameters
    Init_g14pcor(p_vals,num_pars_p,phi_vals,num_pars_phi,lambda_vals,num_pars_lambda);
    Init_linear_pcor(p_vals_lin,num_pars_p_lin,phi_vals_lin,num_pars_phi_lin,lambda_vals_lin,num_pars_lambda_lin);
        
    firstCall = false;
  }
  
  // Get the kinematic quantities we need
  
  double px = __p3[0],py = __p3[1],pz = __p3[2];
  double p = sqrt(px*px + py*py + pz*pz),p_cor;
  double theta_lab = atan2(sqrt(px*px + py*py),pz);
  double phi_lab = atan2(py,px);
  double lam = g14GetLambdaTrack(__p3),lam_cor;
  double phi = g14GetPhiTrack(__p3),phi_cor;
  double par[4];
  int sector,q_index;
  sector = g14GetSector(__p3);
 
  //by default this is momCor for negative torrus field. So the negative charge correction was obtained from pim(negative field) 
  //and proton with positive(field), 
  //and the positive charge correction was from the proton (negative field) and pim(positive field).
  //HOWEVER, silver 1 and 2 has opposite torrus field, therefore the correction for positive charge was from pim and negative 
  //charge correction was from proton. 
  if(68020<__runNo && __runNo<68180) //silver period 1 and 2 with positive torrus field
    {
      if(__q < 0) {q_index = 0;}
    else {q_index = 1;}
    }

  else
    {
      if(__q < 0) {q_index = 1;}//neg
      else {q_index = 0;} //pos

    }
  //------------------------------------------------------------//


  bin =g14GetBin(theta_lab,phi_lab,phiBins_min,phiBins_max,thetaBins_min,
		 thetaBins_max);

  if(bin < 0) {// no correction for this (theta,phi)
    for(int i = 0; i < 3; i++) __p3cor[i] = __p3[i];
    return;
  }
 
  //------------------------------------------------------------------------//
  if(69374<__runNo && __runNo<69637)
    { //linear run
       // correct |p|
  double p_scale = 1.0;
  
  for(int i = 0; i < 4; i++) par[i] = p_vals_lin[q_index][sector-1][bin][i];
  
  if(num_pars_p_lin[q_index][sector-1][bin] == 4)
    p_cor = p + p_scale*(par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_p_lin[q_index][sector-1][bin] == 2)
    p_cor = p + p_scale*(par[1]*p + par[0]);
  
  else if(num_pars_p_lin[q_index][sector-1][bin] == 3)
    p_cor = p + p_scale*(par[2]*p*p + par[1]*p + par[0]);
  
  else p_cor = p;
  
  
  // correct tracking angle lambda
  for(int i = 0; i < 4; i++) par[i] = lambda_vals_lin[q_index][sector-1][bin][i];
  
  if(num_pars_lambda_lin[q_index][sector-1][bin] == 4)
    lam_cor = lam + (par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_lambda_lin[q_index][sector-1][bin] == 3)
    lam_cor = lam + (par[2]*(p*p) + par[1]*p +par[0]);
  
  else if(num_pars_lambda_lin[q_index][sector-1][bin] == 2)
    lam_cor = lam + (par[1]*p + par[0]);
  
  else lam_cor = lam;
 
  

  // correct tracking angle phi
  for(int i = 0; i < 4; i++) par[i] = phi_vals_lin[q_index][sector-1][bin][i];
  
  if(num_pars_phi_lin[q_index][sector-1][bin] == 4)
    phi_cor = phi + (par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_phi_lin[q_index][sector-1][bin] == 2)
    phi_cor = phi + (par[1]*p + par[2]);
  
  else if(num_pars_phi_lin[q_index][sector-1][bin] == 3)
    phi_cor = phi + (par[2]*p*p + par[1]*p + par[0]);

  else phi_cor = phi;
    } //end linear run correction

 //------------------------------------------------------------------------//

  else{ //all the rest 
  // correct |p|
  double p_scale = 1.0;
  
  for(int i = 0; i < 4; i++) par[i] = p_vals[q_index][sector-1][bin][i];
  
  if(num_pars_p[q_index][sector-1][bin] == 4)
    p_cor = p + p_scale*(par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_p[q_index][sector-1][bin] == 2)
    p_cor = p + p_scale*(par[1]*p + par[0]);
  
  else if(num_pars_p[q_index][sector-1][bin] == 3)
    p_cor = p + p_scale*(par[2]*p*p + par[1]*p + par[0]);
  
  else p_cor = p;
  
  
  // correct tracking angle lambda
  for(int i = 0; i < 4; i++) par[i] = lambda_vals[q_index][sector-1][bin][i];
  
  if(num_pars_lambda[q_index][sector-1][bin] == 4)
    lam_cor = lam + (par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_lambda[q_index][sector-1][bin] == 3)
    lam_cor = lam + (par[2]*(p*p) + par[1]*p +par[0]);
  
  else if(num_pars_lambda[q_index][sector-1][bin] == 2)
    lam_cor = lam + (par[1]*p + par[0]);
  
  else lam_cor = lam;
 
  

  // correct tracking angle phi
  for(int i = 0; i < 4; i++) par[i] = phi_vals[q_index][sector-1][bin][i];
  
  if(num_pars_phi[q_index][sector-1][bin] == 4)
    phi_cor = phi + (par[3]*(p*p*p) + par[2]*(p*p) + par[1]*p + par[0]);
  
  else if(num_pars_phi[q_index][sector-1][bin] == 2)
    phi_cor = phi + (par[1]*p + par[2]);
  
  else if(num_pars_phi[q_index][sector-1][bin] == 3)
    phi_cor = phi + (par[2]*p*p + par[1]*p + par[0]);

  else phi_cor = phi;
  } //all the rest

 //------------------------------------------------------------------------//

  

  if(!(-0.05<(phi_cor-phi) && (phi_cor -phi)<0.05 )) {phi_cor=phi;}
  
  if(!(-0.05<(lam_cor-lam) && (lam_cor -lam)<0.05 )) {lam_cor=lam;}
  
  if(!(-0.03<(p_cor-p) && (p_cor -p)<0.03 )) {p_cor=p;}

  //-----------------------------------------------------------------------------//


  
    //addtional correction for silver 1 runs:
  if(68020<__runNo && __runNo<68093) 
    {
      if((q_index==0 &&  p<0.46) || (q_index==1 && p>1.0))
	{p_cor += 0.15*(p_cor-p);}
      
      if(q_index==0)//negative  
	{lam_cor += 0.015; phi_cor+= 0.0007; p_cor += -0.001;}
      else {lam_cor += -0.01; phi_cor+= -0.0007; p_cor += -0.002;}
    }

   //addtional correction for silver 2 runs:
  if(68093<__runNo && __runNo<68180) 
    {
      if((q_index==0 &&  p<0.46) || (q_index==1 && p>1.0))
	{p_cor += 0.15*(p_cor-p);}
      
      if(q_index==0)//negative  
	{lam_cor += 0.015; phi_cor+= 0.0007; p_cor += -0.0005;}
      else {lam_cor += -0.005; phi_cor+= -0.0005; p_cor += -0.001;}   
    }

 //-----------------------------------------------------------------------------//
  

  //addtional correction for silver 3 runs: 
  if(68187<__runNo && __runNo<68231)
    {
      if((q_index==1 &&  p<0.46) || (q_index==0 && p>1.0))
	{p_cor += 0.15*(p_cor-p);}
      
      if(q_index==1)//negative  
	{lam_cor+= 0.028;p_cor += -0.001;}
      else {lam_cor += -0.005;}
    }

 //addtional correction for silver 4 runs:
  if(68231<__runNo && __runNo<68310) 
    {
      if((q_index==1 &&  p<0.46) || (q_index==0 && p>1.0))
	{p_cor += 0.15*(p_cor-p);}
      
      if(q_index==1)//negative  
	{lam_cor+= -0.009;p_cor += -0.001;}
      else {lam_cor += +0.005;}
    }
  

//-----------------------------------------------------------------------------//

   //addtional correction for gold 2 runs:
   if(69193<__runNo && __runNo<69374)  
    {
     
      if(q_index==1)//negative  
	{lam_cor+= 0.025; phi_cor += -0.0005;}
      else {lam_cor += -0.007; phi_cor += -0.0005;}
    }


//-----------------------------------------------------------------------------//  

  //addtional correction for linear runs:
   if(69374<__runNo && __runNo<69637) 
    {
      p_cor += -0.003;
      if(q_index==1)//negative  
	{lam_cor+= 0.015; phi_cor +=0.0005;}
      else {lam_cor += -0.01; phi_cor +=0.0015;}
     
    }
      
//-----------------------------------------------------------------------------//
 
//addtional correction for silver 5 and empty runs:
   if((68340<__runNo && __runNo<68771) ||(68990<__runNo && __runNo<69038))  
    {
      if(q_index==1)//negative  
	lam_cor+= -0.007;
      else {lam_cor += 0.007;}
    }
  



  //----------------------RETURN VALUES------------------------------------------//

  // set corrected 3-momenta
  g14SetP3FromTrackingParameters(p_cor,lam_cor,phi_cor,sector,__p3cor);  
  
  
  }

//_____________________________________________________________________________

/// Reads the binning.dat file to get (theta,phi) bins
void g14ReadBinsFile(double __phiBins_min[12],double __phiBins_max[12],
		  double __thetaBins_min[15],double __thetaBins_max[15]){

  char file[200];
  sprintf(file,"%s/binning.dat",files_location);
  std::ifstream binsFile(file);
  int num_bins;

  if(!binsFile.is_open()){ // couldn't open the file
    printf("Error!!! <ReadBinsFile> Could NOT open binning.dat. ");
    printf("This function is called by g14MomCor and should be located in ");
    printf("%s, make sure your CLAS_PACK ",files_location);
    printf("enviornment variable is set properly. \n");
    abort();
  }

  std::string str;
  // skip over the commented lines (lines that start with #)
  while(binsFile >> str){
    if(str[0] != '#') break;
    else binsFile.ignore(1000,'\n');
  }

  // get phi bins
  num_bins = atoi(str.c_str());
  if(num_bins != 12){ 
    printf("Error!!! <ReadBinsFile> Read incorrect number of phi bins ");
    printf("[%d instead of 12].",num_bins);
    printf("This function is called by g14MomCor and should be located in ");
    printf("%s, make sure your CLAS_PACK ",files_location);
    printf("enviornment variable is set properly. \n");
    abort();
  }
  for(int i = 0; i < num_bins; i++) 
    binsFile >> __phiBins_min[i] >> __phiBins_max[i];
  
  // skip next block of commented lines
  while(binsFile >> str){
    if(str[0] != '#') break;
    else binsFile.ignore(1000,'\n');
  }

  // get theta bins
  num_bins = atoi(str.c_str());
  if(num_bins != 15){ 
    printf("Error!!! <ReadBinsFile> Read incorrect number of theta bins ");
    printf("[%d instead of 15].",num_bins);
    printf("This function is called by g14MomCor and should be located in ");
    printf("%s, make sure your CLAS_PACK ",files_location);
    printf("enviornment variable is set properly. \n");
    abort();
  }

  for(int i = 0; i < num_bins; i++) 
    binsFile >> __thetaBins_min[i] >> __thetaBins_max[i];

  binsFile.close(); // close the file
}
//_____________________________________________________________________________

/// Gets CLAS sector number from 3-momentum
int g14GetSector(const double __p3[3]){

  int sector = 0;
  double px = __p3[0],py = __p3[1];//,pz = __p3[2];
  double pi = 3.14159;
  double phi_lab = atan2(py,px);
  double phi = (180./pi)*phi_lab; // radians --> degrees

  if(std::abs(phi) <= 30.) sector = 1;
  else if(phi > 0.){
    if(phi <= 90.) sector = 2;
    else if(phi <= 150) sector = 3;
    else sector = 4;
  }
  else {
    // phi < 0
    if(std::abs(phi) <= 90.) sector = 6;
    else if(std::abs(phi) <= 150.) sector = 5;
    else sector = 4;
  }
  return sector;
}
//_____________________________________________________________________________

/// Get (theta,phi) bin
int g14GetBin(double __theta,double __phi,const double __phiBins_min[12],
	   const double __phiBins_max[12],const double __thetaBins_min[15],
	      const double __thetaBins_max[15]){

  double sec_phi = g14GetSectorPhi(__phi);
  sec_phi *= 180./3.14159; // convert to degrees
  double theta = __theta * 180./3.14159; // convert to degrees

  int phi_bin = -1;
  for(int i = 0; i < 12; i++){
    if(sec_phi >= __phiBins_min[i] && sec_phi < __phiBins_max[i]){
      phi_bin = i;
      break;
    }
  }
  if(phi_bin == -1) return -1;

  int theta_bin = -1;
  for(int i = 0; i < 15; i++){
    if(theta >= __thetaBins_min[i] && theta < __thetaBins_max[i]){
      theta_bin = i;
      break;
    }
  }
  if(theta_bin == -1) return -1;
 
  return 12*theta_bin + phi_bin;  
}
//_____________________________________________________________________________

/// Calculates the tracking angle \f$ \lambda \f$.
double g14GetLambdaTrack(const double __p3[3]){

  double lambda;
  double px = __p3[0],py = __p3[1],pz = __p3[2];
  double p_mag = sqrt(px*px + py*py + pz*pz);
  double x = px/p_mag,y = py/p_mag;
  
  double alpha = (3.14159/3.)*(g14GetSector(__p3)-1);
  lambda = asin(cos(alpha)*y - sin(alpha)*x);
  return lambda;
}
//_____________________________________________________________________________

/// Calculates the tracking angle \f$ \phi \f$.
double g14GetPhiTrack(const double __p3[3]){

  double phi;
  double px = __p3[0],py = __p3[1],pz = __p3[2];
  double z = pz/sqrt(px*px + py*py + pz*pz); // normalized z_lab
  double lambda = g14GetLambdaTrack(__p3);

  phi = acos(z/cos(lambda));
  return phi;
}
//_____________________________________________________________________________

void g14SetP3FromTrackingParameters(double __p,double __lam,double __phi,
				 int __sector,double __p3[3]){

  double alpha = (3.14159/3.)*(__sector - 1);
  __p3[0] = __p*(cos(__lam)*sin(__phi)*cos(alpha) - sin(__lam)*sin(alpha));
  __p3[1] = __p*(cos(__lam)*sin(__phi)*sin(alpha) + sin(__lam)*cos(alpha));
  __p3[2] = __p*cos(__lam)*cos(__phi);    
}
//_____________________________________________________________________________

double g14GetSectorPhi(double __phi){

  double sec_phi = __phi;
  int sign = 1;
  if(__phi < 0) sign = -1;
  if(std::abs(sec_phi) < 3.14159/6.) return __phi;
  else{
    sec_phi -= sign*3.14159/3.;
    return g14GetSectorPhi(sec_phi);
  }
}
//_____________________________________________________________________________

/// Reads in the momentum correction parameters
void Init_g14pcor(double __p_vals[2][6][180][4],int __num_pars_p[2][6][180],
		  double __phi_vals[2][6][180][4],
		  int __num_pars_phi[2][6][180],
		  double __lambda_vals[2][6][180][4],
		  int __num_pars_lambda[2][6][180]){
		 
  int bin=0,num_pars=0;
  std::ifstream *inFile = 0;
  double x[10];
  std::string file_base = files_location;
  file_base += "/pars_";
  char fileName[200];
  std::string charge_ext[2];
  charge_ext[0] = "pos";
  charge_ext[1] = "neg";

  for(int s = 1; s <= 6; s++){ // loop over sectors
    for(int qi = 0; qi < 2; qi++){ // loop over charges

      // get |p| parameters    
      sprintf(fileName,"%sp_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){ // check that file was successfully opened
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");	
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	
	__num_pars_p[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];
	for(int i = 0; i < 4; i++){
	  __p_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __p_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file
	
      // get lambda parameters
      sprintf(fileName,"%slambda_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");	
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	
	__num_pars_lambda[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];
	for(int i = 0; i < 4; i++){
	  __lambda_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __lambda_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file

      // get phi parameters
      sprintf(fileName,"%sphi_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	//if(num_pars !=4) printf("charged bin phi and num_pars %d %d %f \n",qi, bin, num_pars);
	__num_pars_phi[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];	  
	for(int i = 0; i < 4; i++){
	  __phi_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __phi_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file
     
     

    }
  }

}


void Init_linear_pcor (double __p_vals[2][6][180][4],int __num_pars_p[2][6][180],
		  double __phi_vals[2][6][180][4],
		  int __num_pars_phi[2][6][180],
		  double __lambda_vals[2][6][180][4],
		  int __num_pars_lambda[2][6][180]){
		 
  int bin=0,num_pars=0;
  std::ifstream *inFile = 0;
  double x[10];
  std::string file_base = files_location;
  file_base += "/linear_correction/pars_";
  char fileName[200];
  std::string charge_ext[2];
  charge_ext[0] = "pos";
  charge_ext[1] = "neg";

  for(int s = 1; s <= 6; s++){ // loop over sectors
    for(int qi = 0; qi < 2; qi++){ // loop over charges

      // get |p| parameters    
      sprintf(fileName,"%sp_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){ // check that file was successfully opened
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");	
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	
	__num_pars_p[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];
	for(int i = 0; i < 4; i++){
	  __p_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __p_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file
	
      // get lambda parameters
      sprintf(fileName,"%slambda_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");	
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	
	__num_pars_lambda[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];
	for(int i = 0; i < 4; i++){
	  __lambda_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __lambda_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file

      // get phi parameters
      sprintf(fileName,"%sphi_sector%d.%s",file_base.c_str(),s,
	      charge_ext[qi].c_str());
      inFile = new std::ifstream(fileName);
      if(!(inFile->is_open())){
	printf("Error!!!! <Init_g14pcor> Could NOT open %s. ",fileName);
	printf("This function is called by g14MomCor to read in the ");
	printf("correction parameters. This function looks for files in ");
	printf("%s, make sure your CLAS_PACK ",files_location);
	printf("enviornment variable is set correctly. \n");
	abort();
      }

      while(*inFile >> bin){
	*inFile >> num_pars;
	//if(num_pars !=4) printf("charged bin phi and num_pars %d %d %f \n",qi, bin, num_pars);
	__num_pars_phi[qi][s-1][bin] = num_pars;
	for(int i = 0; i < num_pars; i++) *inFile >> x[i];	  
	for(int i = 0; i < 4; i++){
	  __phi_vals[qi][s-1][bin][i] = 0.;
	  if(i < num_pars) __phi_vals[qi][s-1][bin][i] = x[i];
	}
      }
      delete inFile; inFile = 0; // close the file
     
     

    }
  }

}
