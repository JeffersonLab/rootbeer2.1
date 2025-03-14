#define ELECTRON_MASS     0.00051100
#define MUON_MASS         0.10566
#define PI_CHARGED_MASS   0.13957
#define PI_ZERO_MASS      0.13498
#define KAON_CHARGED_MASS 0.49368
#define KAON_ZERO_MASS    0.49767
#define ETA_MASS          0.54745
#define PROTON_MASS       0.93827
#define NEUTRON_MASS      0.93957 
#define OMEGA_MASS        0.78194
#define RHO_MASS          0.770
#define ETA_PRIME_MASS    0.95778
#define PHI_MASS          1.019413
#define DEUTERON_MASS     1.8756
#define TRITON_MASS       2.1573
#define LAMBDA_MASS       1.115683
#define SIGMA_ZERO_MASS   1.192642

typedef enum {

  /*
   * These constants are defined to be
   * same as GEANT. see http://wwwcn.cern.ch/asdoc/geant/H2GEANTCONS300.html
   * for more details.
  */

  Unknown        =  0,
  Gamma          =  1,
  Positron       =  2,
  Electron       =  3,
  Neutrino       =  4,
  MuonPlus       =  5,
  MuonMinus      =  6,
  Pi0            =  7,
  PiPlus         =  8,
  PiMinus        =  9,
  KLong          = 10,
  KPlus          = 11,
  KMinus         = 12,
  Neutron        = 13,
  Proton         = 14,
  AntiProton     = 15,
  KShort         = 16,
  Eta            = 17,
  Lambda         = 18,
  SigmaPlus      = 19,
  Sigma0         = 20,
  SigmaMinus     = 21,
  Xi0            = 22,
  XiMinus        = 23,
  OmegaMinus     = 24,
  AntiNeutron    = 25,
  AntiLambda     = 26,
  AntiSigmaMinus = 27,
  AntiSigma0     = 28,
  AntiSigmaPlus  = 29,
  AntiXi0        = 30,
  AntiXiPlus     = 31,
  AntiOmegaPlus  = 32,

  Deuteron       = 45,
  Triton         = 49,

  Rho0           = 57,
  RhoPlus        = 58,
  RhoMinus       = 59,
  omega          = 60,
  EtaPrime       = 61,
  phiMeson       = 62

} Particle_t;

/* end file */
