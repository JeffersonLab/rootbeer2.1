#ifndef RBANA_CONSTANTS
#define RBANA_CONSTANTS

#define RADDEG   57.2957795
#define CLIGHT   29.97925      // cm/ns
#define HBAR    6.58212e-16    // GeV*ns 
#ifndef PI
#define PI      3.14159265359
#endif
#define TWOPI   6.28318530718
#define PIBY2   1.57079632680

#define M_ELECTRON 0.000511
#define M_PROTON   0.938272
#define M_NEUTRON  0.939565
#define M_DELTA    1.231 
#define M_LAMBDA   1.115683
#define M_SIGMAPLUS 1.18937
#define M_SIGMA0   1.192642
#define M_PION     0.139700
#define M_PI0      0.134900
#define M_KAON     0.493677
#define M_KA0      0.497672
#define M_OMEGA    0.78194
#define M_RHO      0.770 

#define G_DELTA    0.120 
#define G_LAMBDA   2.5e-15
#define G_SIGMAPLUS 8.21e-15
#define G_SIGMA0   8.9e-6
#define G_PION     2.53e-17
#define G_PI0      7.8e-9
#define G_KAON     5.315e-17
#define G_KA0S     7.352e-15
#define G_KA0L     1.27e-17
#define G_OMEGA    0.00826
#define G_RHO      0.150 

// Monte Carlo particle numbering convention (PDG200, sect.30)
#define MCID_PROTON   2212
#define MCID_NEUTRON  2112
#define MCID_KAPLUS   321 
#define MCID_KAMINUS  -321
#define MCID_PIPLUS   211 
#define MCID_PIMINUS  -211
#define MCID_PI0      111 
#define MCID_PHOTON   22 
#define MCID_ELECTRON 11 

// GEANT3 particle id
#define GID_PROTON   14
#define GID_NEUTRON  13
#define GID_KAPLUS   11
#define GID_KAMINUS  12
#define GID_PIPLUS   8 
#define GID_PIMINUS  9
#define GID_PI0      7
#define GID_PHOTON   1
#define GID_ELECTRON 2

#ifndef ABS
#define ABS(x)   ((x) < 0 ? -(x) : (x))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef POW2
#define POW2(x)  ((x)*(x))
#endif

#endif
