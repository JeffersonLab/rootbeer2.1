// g14pcor package header file -*- C++ -*-
//#ifndef _g14pcor_H
//#define _g14pcor_H

// For more detailed descriptions of the functions, see g14pcor.cc

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
//_____________________________________________________________________________
  // Fiducial Cuts:
  /* C/C++ Function:
   * Parameters: p3[3]  Measured 3-momentum [x,y,z]
   *             q      Charge (ex 1 for proton)
   *
   * Returns: true = good fiducial region, false otherwise
   */
  bool g14FiducialCut(double __p3[3],int __q);

  /* FORTRAN Subroutine:
   * Parameters: p3[3]   Measured 3-momentum [x,y,z]
   *             q       Charge (ex 1 for pi+)
   *             cut     Set to 1 if good fiducial region, 0 otherwise
   */
  void g14FiducialCut_(double __p3[3],int *__q,int *__cut);
//_____________________________________________________________________________
}

//#endif /* _g14pcor_H */
