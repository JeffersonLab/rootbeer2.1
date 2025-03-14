#ifndef ROOTBEER_BANKHEADER
#define ROOTBEER_BANKHEADER

// WARNING - This file is generated from clasbanks.ddl by the command
//           Using the script bank2struct.gk
//	      Edit at your peril!!

typedef struct BEAM_t {		//  BEAM bank for mySQL 
	float	ENERGY;         //  Electron beam energy in GeV
	float	ITORUS;         //  Torus magnet current in Amps
	float	IMINI;          //  Mini-torus magnet current in Amps
	float	ITAG;           //  Tagger magnet current in Amps
} BEAM_t;
typedef struct BMPR_t {		//  BEAM parameters bank 
	float	Q_TOT;          //  Total Charge (Coulomb)
	float	Q_TOT_LV;       //  Total charge in Live Time (Coulomb) 
	float	TL_TU;          //  ( Time Live ) / ( Time Ungated )
	float	CURRENT;        //  Current in Ampers
} BMPR_t;
typedef struct BREP_t {		//  Begin run epics values 
	float	MBSY2C_energy;  //  beam energy
	float	IGT0I00BIASET;  //  thermionic gun
	float	IGL1I00DAC2;    //  B polarized gun
	float	SMRPOSA;        //  A slit position
	float	SMRPOSB;        //  B slit position
	float	SMRPOSC;        //  C slit position
	float	Harp;           //  harp
	float	hallb_sf_xy560; //  torus current
	float	MTSETI;         //  Mini set current
	float	MTIRBCK;        //  Mini current readback
	float	MTVRBCK;        //  Mini voltage readback
	float	TMSETI;         //  Tagger set current
	float	TMIRBCK;        //  Tagger current readback
	float	TMVRBCK;        //  Tagger voltage readback
	float	Cryo_pressure;  //  cryotarget pressure
	float	Cryo_temperature; //  cryotarget temperature
	float	Cryo_status;    //  cryotarget status
	float	VCG2C24;        //  upstream beam vacuum
	float	VCG2H01;        //  target vacuum
	float	scalerS2o;      //  Halo UP upstream
	float	scalerS3o;      //  Halo DOWN upstream
	float	scalerS4o;      //  Halo LEFT upstream
	float	scalerS5o;      //  Halo RIGHT upstream
	float	scalerS6o;      //  Halo UP downstream
	float	scalerS7o;      //  Halo DOWN downstream
	float	scalerS8o;      //  Halo LEFT downstream
	float	scalerS9o;      //  Halo RIGHT downstream
	float	IPM2H01_XPOS;   //  bpm 1 x
	float	IPM2H01_YPOS;   //  bpm 1 y
	float	IPM2H01;        //  bpm 1 current
	float	IPM2C24A_XPOS;  //  bpm 2 x
	float	IPM2C24A_YPOS;  //  bpm 2 y
	float	IPM2C24A;       //  bpm 2 current
	float	IPM2C22A_XPOS;  //  bpm 3 x
	float	IPM2C22A_YPOS;  //  bpm 3 y
	float	IPM2C22A;       //  bpm 3 current
} BREP_t;
typedef struct CALB_t {		//  Monhist fit results for mySQL 
	float	meanRFe;        //  RF offset for electrons (all sectors)
	float	sigmaRFe;       //  Time resolution for electrons (RF)
	float	sigmaRFh;       //  Time resolution for pions
	float	sigmaECt;       //  Time resolution of EC, tEC(e)-tSC(e) 
	float	SFECe;          //  Sampling fraction E_EC(e)/p(e)
	float	sigmaSF;        //  width of the sampling fraction
	float	ResSL1;         //  DC residuals in R1 (all sectors)
	float	ResSL2;         //  DC residuals in R2 (all sectors)
	float	ResSL3;         //  DC residuals in R3 (all sectors)
	float	ResSL4;         //  DC residuals in R1 (all sectors)
	float	ResSL5;         //  DC residuals in R2 (all sectors)
	float	ResSL6;         //  DC residuals in R3 (all sectors)
} CALB_t;
typedef struct CALL_t {		//  Catch-ALL event bank 
	unsigned short	ID;             //  catch-all element(RF,laser diode,etc)
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} CALL_t;
typedef struct CC01_t {		//  Cerenkov Counter hits bank 
	int	id;             //  the address of the hit detector element
	float	time;           //  time(ns) 
	float	n_pe;           //  number of photoelectrons
} CC01_t;
typedef struct CC1_t {		//  Large angle Cerenkov counters event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} CC1_t;
typedef struct CC_t {		//  Forward Cerenkov counters event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} CC_t;
typedef struct CCDI_t {		//  CC channel discriminator thresholds 
	int	crate;          //  CAMAC crate number
	int	slot;           //  slot
	int	mask;           //  mask
	int	threshold;      //  actual threshold value (mV)
	int	width;          //  actual width value
} CCDI_t;
typedef struct CCH_t {		//  GSIM Cerenkov hit info 
	int	tknum;          //  track number, 1000*istak+itra
	int	id;             //  track PDG id
	int	nhits;          //  number of CC hits per track
	int	sector;         //  sector number of track
	int	segment;        //  cc segment number of track
	int	pmom;           //  track momentum 
	float	xin;            //  x pos track entry into cerenkov
	float	yin;            //  y pos
	float	zin;            //  z pos
	float	xout;           //  x pos track exit from cerenkov
	float	yout;           //  y pos
	float	zout;           //  z pos
} CCH_t;
typedef struct CCMT_t {		//  Mean CC pretrigger thresholds 
	int	mean_high;      //  mean high threshold (mV)
	int	mean_lo;        //  mean lo threshold (mV)
} CCMT_t;
typedef struct CCPB_t {		//  CC hits involved in the event 
	int	ScSgHt;         //  100*sector + Cluster # in CCRC
	float	Nphe;           //  Number of photo-electrons*10
	float	Time;           //  Flight time relative to the evnt start time
	float	Path;           //  Path lenght from target (from tracking)
	float	Chi2CC;         //  Geometrical matching: angle between CC hit and
	int	Status;         //  Status word - now 10*(CC segment number) 
} CCPB_t;
typedef struct CCPE_t {		//  Translated CC pedestal bank, measured during the pedestal run. 
	int	ID;             //  the address of the hit detector element
	int	mean;           //  adc pedestal mean value (channel)
	int	sigma;          //  sigma of the pedestal distribution (channel)
} CCPE_t;
typedef struct CCRC_t {		//  Cerenkov reconstruction bank 
	int	nrsect;         //  Sector #
	int	nrsegm;         //  10 * Mean segment #
	int	nrsegm_p;       //  10 * Max segment # in the cluster
	int	nrsegm_m;       //  10 * Min segment # in the cluster
	int	nrphe;          //  Number of photoelectrons obtained
	int	nrtime;         //  TOF in channels (50psec/channel)
	int	nrthet;         //  Estimated angle Theta
	int	nrdthet;        //  Estimated error of angle Theta
	int	nrphy;          //  Phy index 
	int	nriec;          //  Estimated i-coordinate in EC 
	int	nrdiec;         //  Estimated error of i-coord. in EC
	int	nrstat;         //  Status  word (yet unclear) 
} CCRC_t;
typedef struct CCS_t {		//  Cerenkov scaler bank 
	int	CCS1;           //  scaler sector 1
	int	CCS2;           //  scaler sector 1
	int	CCS3;           //  scaler sector 1
	int	CCS4;           //  scaler sector 1
	int	CCS5;           //  scaler sector 1
	int	CCS6;           //  scaler sector 1
	int	CCS7;           //  scaler sector 1
	int	CCS8;           //  scaler sector 1
	int	CCS9;           //  scaler sector 1
	int	CCS10;          //  scaler sector 1
	int	CCS11;          //  scaler sector 1
	int	CCS12;          //  scaler sector 1
	int	CCS13;          //  scaler sector 1
	int	CCS14;          //  scaler sector 1
	int	CCS15;          //  scaler sector 1
	int	CCS16;          //  scaler sector 1
	int	CCS17;          //  scaler sector 2
	int	CCS18;          //  scaler sector 2
	int	CCS19;          //  scaler sector 2
	int	CCS20;          //  scaler sector 2
	int	CCS21;          //  scaler sector 2
	int	CCS22;          //  scaler sector 2
	int	CCS23;          //  scaler sector 2
	int	CCS24;          //  scaler sector 2
	int	CCS25;          //  scaler sector 2
	int	CCS26;          //  scaler sector 2
	int	CCS27;          //  scaler sector 2
	int	CCS28;          //  scaler sector 2
	int	CCS29;          //  scaler sector 2
	int	CCS30;          //  scaler sector 2
	int	CCS31;          //  scaler sector 2
	int	CCS32;          //  scaler sector 2
	int	CCS33;          //  scaler sector 3
	int	CCS34;          //  scaler sector 3
	int	CCS35;          //  scaler sector 3
	int	CCS36;          //  scaler sector 3
	int	CCS37;          //  scaler sector 3
	int	CCS38;          //  scaler sector 3
	int	CCS39;          //  scaler sector 3
	int	CCS40;          //  scaler sector 3
	int	CCS41;          //  scaler sector 3
	int	CCS42;          //  scaler sector 3
	int	CCS43;          //  scaler sector 3
	int	CCS44;          //  scaler sector 3
	int	CCS45;          //  scaler sector 3
	int	CCS46;          //  scaler sector 3
	int	CCS47;          //  scaler sector 3
	int	CCS48;          //  scaler sector 3
	int	CCS49;          //  scaler sector 4
	int	CCS50;          //  scaler sector 4
	int	CCS51;          //  scaler sector 4
	int	CCS52;          //  scaler sector 4
	int	CCS53;          //  scaler sector 4
	int	CCS54;          //  scaler sector 4
	int	CCS55;          //  scaler sector 4
	int	CCS56;          //  scaler sector 4
	int	CCS57;          //  scaler sector 4
	int	CCS58;          //  scaler sector 4
	int	CCS59;          //  scaler sector 4
	int	CCS60;          //  scaler sector 4
	int	CCS61;          //  scaler sector 4
	int	CCS62;          //  scaler sector 4
	int	CCS63;          //  scaler sector 4
	int	CCS64;          //  scaler sector 4
	int	CCS65;          //  scaler sector 5
	int	CCS66;          //  scaler sector 5
	int	CCS67;          //  scaler sector 5
	int	CCS68;          //  scaler sector 5
	int	CCS69;          //  scaler sector 5
	int	CCS70;          //  scaler sector 5
	int	CCS71;          //  scaler sector 5
	int	CCS72;          //  scaler sector 5
	int	CCS73;          //  scaler sector 5
	int	CCS74;          //  scaler sector 5
	int	CCS75;          //  scaler sector 5
	int	CCS76;          //  scaler sector 5
	int	CCS77;          //  scaler sector 5
	int	CCS78;          //  scaler sector 5
	int	CCS79;          //  scaler sector 5
	int	CCS80;          //  scaler sector 5
	int	CCS81;          //  scaler sector 6
	int	CCS82;          //  scaler sector 6
	int	CCS83;          //  scaler sector 6
	int	CCS84;          //  scaler sector 6
	int	CCS85;          //  scaler sector 6
	int	CCS86;          //  scaler sector 6
	int	CCS87;          //  scaler sector 6
	int	CCS88;          //  scaler sector 6
	int	CCS89;          //  scaler sector 6
	int	CCS90;          //  scaler sector 6
	int	CCS91;          //  scaler sector 6
	int	CCS92;          //  scaler sector 6
	int	CCS93;          //  scaler sector 6
	int	CCS94;          //  scaler sector 6
	int	CCS95;          //  scaler sector 6
	int	CCS96;          //  scaler sector 6
} CCS_t;
typedef struct CCT_t {		//  Forward Cerenkov counters event bank (TDC) 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
} CCT_t;
typedef struct CHI2_t {		// Kinematic Fitter Info Bank 
	float	chi2;           //  overall chi2 for event  
	float	cl;             //  percentage (confidence level)
	int	ndf;            //  number of effective degrees of freedom
	int	iter;           
} CHI2_t;
typedef struct CL01_t {		//  Catch-ALL event bank 
	int	ac_amp;         //  60 Hz AC amplitude (pedestal subtracted)
	int	fc_diode_amp;   //  Forward Carriage diode amplitude (ped sub.)
	float	fc_diode_t;     //  Forward Carriage diode time
	int	nc_diode_amp;   //  North Clamshell diode amplitude (ped sub.)
	float	nc_diode_t;     //  Forward  diode time
	int	sc_diode_amp;   //  Forward Carriage diode amplitude (ped sub.)
	float	sc_diode_t;     //  Forward Carriage diode time(ns)
	int	sf_diode_amp;   //  Forward Carriage diode amplitude (ped sub.)
	float	sf_diode_t;     //  Forward Carriage diode time(ns)
	float	rf1;            //  RF time 1 (ns)
	float	rf2;            //  RF time 2 (ns)
	float	rf;             //  GOOD RF time (ns)
} CL01_t;
typedef struct CLST_t {		//  
	int	Clust;          //  bit packed,  see: include/bosddl.h, clasCLST_t
} CLST_t;
typedef struct CPED_t {		//  Untranslated pedestal bank, measured during the pedestal run. 
	int	slot;           //  ADC slot
	int	channel;        //  ADC channel
	int	mean;           //  adc pedestal mean value (channel)
	float	sigma;          //  sigma of the pedestal distribution (channel)
	int	offset;         //  offset for sparsification threshold calculation
} CPED_t;
typedef struct CSQL_t {		//  Data bank for mySQL 
	unsigned int	EVID;           //  Event ID (number of triggers)
	unsigned int	NPROC;          //  Number of processed triggers
	float	CPU;            //  CPU used (sec) 
	float	FC;             //  Faraday Cup (K)
	float	FCG;            //  Faraday Cup Gated (K)
	float	TG;             //  Clock Gated
	float	IBEAM;          //  Beam current 
	unsigned int	NeS1;           //   Number of electrons in sect 1
	unsigned int	NeS2;           //   Number of electrons in sect 2
	unsigned int	NeS3;           //   Number of electrons in sect 3 
	unsigned int	NeS4;           //   Number of electrons in sect 4 
	unsigned int	NeS5;           //   Number of electrons in sect 5  
	unsigned int	NeS6;           //   Number of electrons in sect 6
	unsigned int	Nhb;            //  Number of HB 
	unsigned int	Ntb;            //  Number of TB
	unsigned int	Nprot;          //  Number of protons
	unsigned int	Npip;           //  number of pip
	unsigned int	Ndeut;          //  number of deutrons
	unsigned int	Nphot;          //  number of photons	
	unsigned int	Nelhp;          //  Number of electrons at pos. Helic.
	unsigned int	Nelhn;          //  Number of electrons at neg. helic.
} CSQL_t;
typedef struct DC0_t {		//  Drift chamber event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
} DC0_t;
typedef struct DC1_t {		//  Drift chamber hits bank 
	int	ID;             //  the address of the hit detector element
	float	time;           //  time(ns) 
} DC1_t;
typedef struct DCDW_t {		//  Drift Chamber Dead Wire bank 
	int	ID;             //  ID_wire 
	float	TIDLY;          //  TIme DeLaY (ns) 		 
	int	STAT;           //  wire status word 
} DCDW_t;
typedef struct DCGM_t {		//  DC geometry info by layer - in sector coordinate system 
	float	x_curve;        //  x center of curvature (cm)
	float	y_curve;        //  y center of curvature (cm)
	float	z_curve;        //  z center of curvature (cm)
	float	r_curve;        //  radius of curvature (cm)
	float	theta_start;    //  angle of first logical wire WRT the center of curvature (rad)
	float	d_theta;        //  delta theta between wires WRT center of curvature
	float	x_nmid;         //  x normal vector to the midplane
	float	y_nmid;         //  y normal vector to the midplane
	float	z_nmid;         //  z normal vector to the midplane
	float	theta_min;      //  theta of first physical wire (rad)
	float	theta_max;      //  theta of last physical wire (rad)
	int	min_wire;       //  minimum physical wire number
	int	max_wire;       //  maximum physical wire number
	float	stereo;         //  approximate stereo angle
	float	cell_size;      //  cell size (cm)
	float	x_norm;         //  x normal vector to the plane(region 1)
	float	y_norm;         //  y normal vector to the plane(region 1)
	float	z_norm;         //  z normal vector to the plane(region 1)
	float	p_dist;         //  distance of plane to origin(cm) (region 1)
	float	p_sep;          //  planar separation(cm) (region 1)
	int	max_cylw;       //  maximum cylindrical wire
} DCGM_t;
typedef struct DCGW_t {		//  DC geometry info by wire - in sector coordinate system 
	float	x_mid;          //  x at midplane of wire(cm)
	float	y_mid;          //  y at midplane of wire(cm)
	float	z_mid;          //  z at midplane of wire(cm)
	float	x_dir;          //  x direction cosine along wire (cm)
	float	y_dir;          //  y direction cosine along wire (cm)
	float	z_dir;          //  z direction cosine along wire (cm)
	float	w_len;          //  wire length from midplane to amplifier (cm)
	float	w_len_hv;       //  wire length from midplane to HV (cm)
} DCGW_t;
typedef struct DCH_t {		//  GSIM Drift chamber hits 
	float	x;              //  x at layer center
	float	y;              //  y at layer center
	float	z;              //  z at layer center
	float	step;           //  step size through layer
	float	dedx;           //  energy deposit in layer
	float	pmom;           //  track momentum at layer center 
	float	time;           //  time of hit at layer center
	float	cx;             //  track x dir cosine at layer center
	float	cy;             //  track y dir cosine at layer center
	float	cz;             //  track z dir cosine at layer center
	int	track;          //  track number
	int	id;             //  track PDG id
	int	layer;          //  layer number
} DCH_t;
typedef struct DCMN_t {		//  
	float	HBT_evt_1;      //  Hit Based tracks per event for Sector 1 
	float	TBT_evt_1;      //  Time Based tracks per event for Sector 1 
	float	HBT_evt_2;      //  Hit Based tracks per event for Sector 2
	float	TBT_evt_2;      //  Time Based tracks per event for Sector 2 
	float	HBT_evt_3;      //  Hit Based tracks per event for Sector 3 
	float	TBT_evt_3;      //  Time Based tracks per event for Sector 3 
	float	HBT_evt_4;      //  Hit Based tracks per event for Sector 4 
	float	TBT_evt_4;      //  Time Based tracks per event for Sector 4 
	float	HBT_evt_5;      //  Hit Based tracks per event for Sector 5 
	float	TBT_evt_5;      //  Time Based tracks per event for Sector 5 
	float	HBT_evt_6;      //  Hit Based tracks per event for Sector 6 
	float	TBT_evt_6;      //  Time Based tracks per event for Sector 6 
	float	HBT_evt_all;    //  Hit Based tracks per event for all sec.  
	float	TBT_evt_all;    //  Time Based tracks per event for all sec. 
	float	Res_s1_sl1;     //  Residual rms for sec 1, superlayer 1  
	float	Res_s1_sl2;     //  Residual rms for sec 1, superlayer 2  
	float	Res_s1_sl3;     //  Residual rms for sec 1, superlayer 3  
	float	Res_s1_sl4;     //  Residual rms for sec 1, superlayer 4  
	float	Res_s1_sl5;     //  Residual rms for sec 1, superlayer 5  
	float	Res_s1_sl6;     //  Residual rms for sec 1, superlayer 6  
	float	Res_s2_sl1;     //  Residual rms for sec 2, superlayer 1  
	float	Res_s2_sl2;     //  Residual rms for sec 2, superlayer 2  
	float	Res_s2_sl3;     //  Residual rms for sec 2, superlayer 3  
	float	Res_s2_sl4;     //  Residual rms for sec 2, superlayer 4  
	float	Res_s2_sl5;     //  Residual rms for sec 2, superlayer 5  
	float	Res_s2_sl6;     //  Residual rms for sec 2, superlayer 6  
	float	Res_s3_sl1;     //  Residual rms for sec 3, superlayer 1  
	float	Res_s3_sl2;     //  Residual rms for sec 3, superlayer 2  
	float	Res_s3_sl3;     //  Residual rms for sec 3, superlayer 3  
	float	Res_s3_sl4;     //  Residual rms for sec 3, superlayer 4  
	float	Res_s3_sl5;     //  Residual rms for sec 3, superlayer 5  
	float	Res_s3_sl6;     //  Residual rms for sec 3, superlayer 6  
	float	Res_s4_sl1;     //  Residual rms for sec 4, superlayer 1  
	float	Res_s4_sl2;     //  Residual rms for sec 4, superlayer 2  
	float	Res_s4_sl3;     //  Residual rms for sec 4, superlayer 3  
	float	Res_s4_sl4;     //  Residual rms for sec 4, superlayer 4  
	float	Res_s4_sl5;     //  Residual rms for sec 4, superlayer 5  
	float	Res_s4_sl6;     //  Residual rms for sec 4, superlayer 6  
	float	Res_s5_sl1;     //  Residual rms for sec 5, superlayer 1  
	float	Res_s5_sl2;     //  Residual rms for sec 5, superlayer 2  
	float	Res_s5_sl3;     //  Residual rms for sec 5, superlayer 3  
	float	Res_s5_sl4;     //  Residual rms for sec 5, superlayer 4  
	float	Res_s5_sl5;     //  Residual rms for sec 5, superlayer 5  
	float	Res_s5_sl6;     //  Residual rms for sec 5, superlayer 6  
	float	Res_s6_sl1;     //  Residual rms for sec 6, superlayer 1  
	float	Res_s6_sl2;     //  Residual rms for sec 6, superlayer 2  
	float	Res_s6_sl3;     //  Residual rms for sec 6, superlayer 3  
	float	Res_s6_sl4;     //  Residual rms for sec 6, superlayer 4  
	float	Res_s6_sl5;     //  Residual rms for sec 6, superlayer 5  
	float	Res_s6_sl6;     //  Residual rms for sec 6, superlayer 6  
} DCMN_t;
typedef struct DCPB_t {		//  DC tracks involved in the event 
	int	ScTr;           //  100*sector+track_ID in *BTR  
	float	x_SC;           //  x coordinate of track intersection with SC plane 
	float	y_SC;           //  y coordinate of track intersection with SC plane
	float	z_SC;           //  z coordinate of track intersection with SC plane
	float	CX_SC;          //  X dir cosine at (x_SC,y_SC,z_SC)
	float	CY_SC;          //  y dir cosine at (x_SC,y_SC,z_SC)
	float	CZ_SC;          //  z dir cosine at (x_SC,y_SC,z_SC)
	float	X_v;            //  vertex X after fiting to the beam position
	float	Y_v;            //  vertex Y after fiting to the beam position
	float	Z_v;            //  vertex Z after fiting to the beam position
	float	R_v;            //  distance from production vertex to the bemam. 
	float	Chi2;           //  Chisquare of track fitting
	int	Status;         //  Status word
} DCPB_t;
typedef struct DCV1_t {		//  parameters for T->D 
	float	TsR1;           //  Tsmear 
	float	V0R1;           //  drift velocity (slope)
	float	TmR1;           //  Maximum drift time (Tmax)
	float	sp1R1;          //  spare  
	float	sp2R1;          //  spare	 
} DCV1_t;
typedef struct DCV2_t {		//  parameters for T->D 
	float	Ts1R2;          //  Tsmear 
	float	V01R2;          //  drift velocity (slope)
	float	VA1R2;          //  drift velocity function parameter
	float	VB1R2;          //  drift velocity function parameter  
	float	Tm1R2;          //  Maximum drift time (Tmax) 		 
	float	TA1R2;          //  Tmax function parameter
	float	TB1R2;          //  Tmax function parameter
	float	Ts2R2;          //  Tsmear 
	float	V02R2;          //  drift velocity (slope)
	float	VA2R2;          //  drift velocity function parameter
	float	VB2R2;          //  drift velocity function parameter  
	float	Tm2R2;          //  Maximum drift time (Tmax) 		 
	float	TA2R2;          //  Tmax function parameter
	float	TB2R2;          //  Tmax function parameter
} DCV2_t;
typedef struct DCV3_t {		//  parameters for T->D 
	float	TsR3;           //  Tsmear 
	float	V0R3;           //  drift velocity (slope)
	float	TmR3;           //  Maximum drift time (Tmax)
	float	sp1R3;          //  spare  
	float	sp2R3;          //  spare	 
} DCV3_t;
typedef struct DDLY_t {		//  Drift chamber time DeLaY bank 
	int	ID;             //  ID_wire 
	float	TIDLY;          //  TIme DeLaY (ns) 		 
	int	STAT;           //  wire STATus = 100*C + 10*A + B
} DDLY_t;
typedef struct DGEO_t {		//  Drift chamber GEOmetry bank 
	int	ID_sec;         //  ID_sector 
	int	ID_reg;         //  ID_region 
	float	xpos;           //  x misalignment 		 
	float	ypos;           //  y misalignment 		 
	float	zpos;           //  z misalignment 		 
	float	sxpos;          //  sx sine of little x angle 		 
	float	sypos;          //  sy sine of little y angle 
	float	szpos;          //  sz sine of little z angle 
} DGEO_t;
typedef struct DHCL_t {		//  drift chamber hit cluster 
	int	SLY;            //  superlayer
	int	BTRK;           //  track_in_sector# (bit set)
	int	TRKS1;          //  combinations to track segments
	int	TRKS2;          //    "             "       "
	int	WIRE1;          //  1.wire# in this cluster in 1.layer
	int	BWIR1;          //  hits in this layer (starting from WIRE1) (bit st)
	int	WIRE2;          //  1.wire# in this cluster in 1.layer
	int	BWIR2;          //  hits in this layer (starting from WIRE1) (bit st)
	int	WIRE3;          //  1.wire# in this cluster in 1.layer
	int	BWIR3;          //  hits in this layer (starting from WIRE1) (bit st)
	int	WIRE4;          //  1.wire# in this cluster in 1.layer
	int	BWIR4;          //  hits in this layer (starting from WIRE1) (bit st)
	int	WIRE5;          //  1.wire# in this cluster in 1.layer
	int	BWIR5;          //  hits in this layer (starting from WIRE1) (bit st)
	int	WIRE6;          //  1.wire# in this cluster in 1.layer
	int	BWIR6;          //  hits in this layer (starting from WIRE1) (bit st)
} DHCL_t;
typedef struct DITM_t {		//  Discriminator calibration time 
	int	time;           //  time of discriminator calibration
} DITM_t;
typedef struct DOCA_t {		//  Drift chamber event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	DOCA;           //  doca from GSIM (micron)
} DOCA_t;
typedef struct DPCP_t {		//  Translated pair counter  
	int	ID;             //  the address of the hit detector element
	int	MN_mean;        //  adc pedestal mean value (channel)
	float	MN_sigma;       //  sigma of the pedestal distribution (channel)
	int	LT_mean;        //  adc pedestal mean value (channel)
	float	LT_sigma;       //  sigma of the pedestal distribution (channel)
	int	RB_mean;        //  adc pedestal mean value (channel)
	float	RB_sigma;       //  sigma of the pedestal distribution (channel)
	int	LB_mean;        //  adc pedestal mean value (channel)
	float	LB_sigma;       //  sigma of the pedestal distribution (channel)
	int	RT_mean;        //  adc pedestal mean value (channel)
	float	RT_sigma;       //  sigma of the pedestal distribution (channel)
	int	VT_mean;        //  adc pedestal mean value (channel)
	float	VT_sigma;       //  sigma of the pedestal distribution (channel)
} DPCP_t;
typedef struct DPSP_t {		//  Translated pair specrometer  
	int	ID;             //  the address of the hit detector element
	int	mean;           //  adc pedestal mean value (channel)
	float	sigma;          //  sigma of the pedestal distribution (channel)
} DPSP_t;
typedef struct DSPC_t {		//  Down Stream Pair Counter event bank 
	unsigned short	PCID;           //  Id compelled by DAQ (always 1)
	unsigned short	TDCPC;          //  tdc information ( scintillator)
	unsigned short	ADCMN;          //  adc information (main)
	unsigned short	ADCLT;          //  adc information (left top)
	unsigned short	ADCRB;          //  adc information (right bottom)
	unsigned short	ADCLB;          //  adc information (left bottom)
	unsigned short	ADCRT;          //  adc information (right top)
	unsigned short	ADCVE;          //  adc information (veto)
} DSPC_t;
typedef struct DSPS_t {		//  Down Stream Pair Spectrometer event bank 
	unsigned short	ID;             //  paddle ID (left = 1 to 4)( right = 5 to 8)
	unsigned short	TDC;            //  tdc information
	unsigned short	ADC;            //  adc information 
} DSPS_t;
typedef struct DSTC_t {		//  Down Stream Total Absorption Shower Counter event bank 
	unsigned short	TACID;          //  ID compelled by DAQ (always 1)
	unsigned short	TDCTAC;         //  tdc information (on sum)
	unsigned short	ADCLT;          //  adc information (left top)
	unsigned short	ADCRT;          //  adc information (right top)
	unsigned short	ADCLB;          //  adc information (left bottom)
	unsigned short	ADCRB;          //  adc information (right bottom)
	unsigned short	ADCSUM1;        //  adc information (sum scale1)
	unsigned short	ADCSUM2;        //  adc information (sum scale2)
	unsigned short	ADCSUM3;        //  adc information (sum scale3)
} DSTC_t;
typedef struct DTCP_t {		//  Translated TAC   
	int	ID;             //  the address of the hit detector element
	int	LT_mean;        //  adc pedestal mean value (channel)
	float	LT_sigma;       //  sigma of the pedestal distribution (channel)
	int	RT_mean;        //  adc pedestal mean value (channel)
	float	RT_sigma;       //  sigma of the pedestal distribution (channel)
	int	LB_mean;        //  adc pedestal mean value (channel)
	float	LB_sigma;       //  sigma of the pedestal distribution (channel)
	int	RB_mean;        //  adc pedestal mean value (channel)
	float	RB_sigma;       //  sigma of the pedestal distribution (channel)
	int	sum1_mean;      //  adc pedestal mean value (channel)
	float	sum1_sigma;     //  sigma of the pedestal distribution (channel)
	int	sum2_mean;      //  adc pedestal mean value (channel)
	float	sum2_sigma;     //  sigma of the pedestal distribution (channel)
	int	sum3_mean;      //  adc pedestal mean value (channel)
	float	sum3_sigma;     //  sigma of the pedestal distribution (channel)
} DTCP_t;
typedef struct DTRK_t {		//  Drift chamber tracks 
	float	X;              
	float	Y;              
	float	Z;              
} DTRK_t;
typedef struct EC01_t {		//  Scintillation counter hits bank 
	int	ID;             //  the address of the hit detector element
	float	time;           //  time for left paddle(ns) 
	float	energy;         //  energy in left paddle(MeV) 
} EC01_t;
typedef struct EC1_t {		//  Large angle calorimeter event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDCL;           //  tdc information (channels)
	unsigned short	ADCL;           //  adc information (channels)
	unsigned short	TDCR;           //  tdc information (channels)
	unsigned short	ADCR;           //  adc information (channels)
} EC1_t;
typedef struct EC1P_t {		//  
	float	N1x;            //  x component of outward normal to plane
	float	N1y;            //  y component of outward normal to plane
	float	N1z;            //  z component of outward normal to plane
	float	R1n;            //  distanse in cm from origin to plane
} EC1P_t;
typedef struct EC1R_t {		//  LA calorimeter result bank 
	float	E_tot;          //  cluster energy (sum of inner+outer)
	float	dE_tot;         //   error on the cluster energy
	float	t_tot;          //  time
	float	dt_tot;         //  error on the time
	float	x_m;            //  x-in CLAS frame
	float	y_m;            //  y-in CLAS frame
	float	z_m;            //  z-in CLAS frame
	float	dx_m;           //  error on x
	float	dy_m;           //  error on y
	float	dz_m;           //  error on z
	float	E_in;           //   cluster energy in inner layer
	float	t_in;           //  time from inner layer
	float	x_in;           //  lab coordinate , inner layer
	float	y_in;           //  lab coordinate , inner layer
	float	x_out;          //  lab coordinate , outer layer
	float	y_out;          //  lab coordinate , outer layer
	float	x2_in_l;        //  second moment of x inner left
	float	x2_in_r;        //  second moment of x inner right 
	float	y2_in_l;        //  second moment of y inner left
	float	y2_in_r;        //  second moment of y inner right
	float	x2_out_l;       //  second moment of x outer hit left
	float	x2_out_r;       //  second moment of x outer hit right
	float	y2_out_l;       //  second moment of y outer hit left
	float	y2_out_r;       //  second moment of y outer hit right
	int	i_in;           //  cluster center in X inner short layer
	int	j_in;           //  cluster center in Y inner long layer
	int	i_out;          //  cluster center in X outer layers
	int	j_out;          //  cluster center in Y outer layer
	float	a_in_xl;        //  energy sum in inner short left pmts
	float	a_in_xr;        //  energy sum in inner short right pmts
	float	a_in_yl;        //  energy sum in inner long  left pmts
	float	a_in_yr;        //  energy sum in inner long  right pmts
	float	a_out_xl;       //  energy sum in outer short left pmts
	float	a_out_xr;       //  energy sum in outer short right pmts
	float	a_out_yl;       //  energy sum in outer long left pmts
	float	a_out_yr;       //  energy sum in outer long right pmts
	float	t_in_xs;        //  tdc sum in inner short righ+left  for cluster center
	float	t_in_xd;        //  tdc dif in inner short right-left 
	float	t_in_ys;        //  tdc sum in inner long  righ+left
	float	t_in_yd;        //  tdc dif in inner long  right-left
	float	t_out_xs;       //  tdc sum in outer short righ+left 
	float	t_out_xd;       //  tdc dif in outer short right-left 
	float	t_out_ys;       //  tdc sum in outer long righ+left
	float	t_out_yd;       //  tdc dif in outer long right-left 
	int	ibl;            //  LAC  block number
	int	ncluster;       //  1000xNclust4+100xNclust3+10xNclust2+Nclust1 
	int	pmtfired;       //  Number of fired pmt (more than threshold) 
	float	z_in;           //  Z in CLAS frame
	float	z_out;          //  Z out in CLAS frame
	int	istat;          //  status word
} EC1R_t;
typedef struct ECCA_t {		//  Forward EC Calibration Bank 
	int	id;             //  PMT ID
	float	aPED;           //  ADC pedestals (channels)
	float	aSIG;           //  ADC pedestal variance (channels)
	float	aMIP;           //  ADC calibration from MIP (N.I.M.P/ch.)
	float	aMIPu;          //  aMIP Error
	float	aSHR;           //  ADC calibration from showers (GeV/ch.)
	float	aSHRu;          //  aSHR Error
	int	stat;           //  4 byte status word    
} ECCA_t;
typedef struct ECCL_t {		//  Forward EC Calibration Bank 
	int	id;             //  PMT ID
	float	lDB;            //  Stack atten length (database)(cm)
	float	lDBu;           //  lDB Error (cm)
	float	lMIP;           //  Stack atten length (MinIonPart)(cm) 
	float	lMIPu;          //  lMIP Error (cm)
	float	lSHR;           //  Stack atten length (showers)(cm)  
	float	lSHRu;          //  lSHR Error (cm) 
	int	stat;           //  4 byte status word    
} ECCL_t;
typedef struct ECCT_t {		//  Forward EC Calibration Bank 
	int	id;             //  PMT ID
	float	tOFF;           //  TDC offset (channels)
	float	tOFFu;          //  tOFF Error
	float	tGAIN;          //  TDC conversion gain (nS/channel)
	float	tGAINu;         //  tGAIN Error
	float	tW0;            //  Time walk constant (channels)
	float	tW0u;           //  tWOu Error
	float	tW1;            //  Time walk correction parameter 
	float	tW1u;           //  tW1u Error
	float	tVEF;           //  Effective velocity of light (cm/ns)
	float	tVEFu;          //  tVEFu Error
	int	stat;           //  4 byte status word    
} ECCT_t;
typedef struct EC_t {		//  Forward calorimeter event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} EC_t;
typedef struct ECDI_t {		//  EC channel discriminator thresholds 
	int	crate;          //  CAMAC crate number
	int	slot;           //  slot
	int	mask;           //  mask
	int	threshold;      //  actual threshold value (mV)
	int	width;          //  actual width value
} ECDI_t;
typedef struct ECG_t {		//  Calorimeter geometry bank 
	float	L0;             //  distance from the target 
	float	THE0;           //  angle between the beam and perpendicular
	float	YLOW;           
	float	YHI;            
	float	DYLOW;          
	float	DYHI;           
	float	LRTH;           //  thickness of the single layer
	float	TANGR;          
	int	SECTOR;         
	float	PHISEC;         
	float	X0OFF;          
	float	Y0OFF;          
	float	Z0OFF;          
	float	ROTM11;         
	float	ROTM12;         
	float	ROTM13;         
	float	ROTM21;         
	float	ROTM22;         
	float	ROTM23;         
	float	ROTM31;         
	float	ROTM32;         
	float	ROTM33;         
} ECG_t;
typedef struct ECHB_t {		//  Forward calorimeter result bank 
	int	Sect;           //  Sector number & Layer number
	float	E__hit;         //  energy found 
	float	dE_hit;         //  error on the energy found 
	float	t_hit;          //  time found 
	float	dt_hi;          //  error time found 
	float	i_hit;          //  sector rectangular coordinate 
	float	j_hit;          //  sector rectangular coordinate 
	float	di_hit;         //  sector rectangular coordinate error,
	float	dj_hit;         //  sector rectangular coordinate error,
	float	x_hit;          //  lab coordinate,
	float	y_hit;          //  lab coordinate,
	float	z_hit;          //  lab coordinate,
	float	dx_hit;         //  lab coordinate error,
	float	dy_hit;         //  lab coordinate error, 
	float	dz_hit;         //  lab coordinate error,
	float	u2_hit;         //  second moment of u  _hit pattern
	float	v2_hit;         //  second moment of v  _hit pattern
	float	w2_hit;         //  second moment of w  _hit pattern
	float	u3_hit;         //  third moment of u  _hit pattern 
	float	v3_hit;         //  third moment of v  _hit pattern 
	float	w3_hit;         //  third moment of w  _hit pattern 
	float	u4_hit;         //  forth moment of u  _hit pattern 
	float	v4_hit;         //  forth moment of v  _hit pattern 
	float	w4_hit;         //  forth moment of w  _hit pattern 
	float	centr_U;        //  peak position on U axis 
	float	centr_V;        //  peak position on V axis 
	float	centr_W;        //  peak position on W axis 
	float	path_U;         //  path length from hit position to U axis 
	float	path_V;         //  path length from hit position to V axis 
	float	path_W;         //  path length from hit position to W axis 
	int	Nstrp_U;        //  Number of U strips in the hit
	int	Nstrp_V;        //  Number of V strips in the hit
	int	Nstrp_W;        //  Number of W strips in the hit
	int	MatchID1;       //  Id of matched hit in the layer1
	float	CH21;           //  Quality measure of matching with layer1
	int	MatchID2;       //  Id of matched hit in the layer2
	float	CH22;           //  Quality measure of matching with layer2
	int	istat;          //  Number of hits & hit ID
} ECHB_t;
typedef struct ECH_t {		//  GSIM Calorimeter track hit info 
	float	x;              //  x of hit
	float	y;              //  y of hit
	float	z;              //  z of hit
	float	cx;             //  track x dir cosine
	float	cy;             //  track y dir cosine
	float	cz;             //  track z dir cosine
	float	pmom;           //  track momentum
	int	id;             //  track PDG id
	float	tof;            //  time of flight
	float	edepin;         //  deposit energy (inner part)
	float	edepout;        //  deposit energy (outer part)
} ECH_t;
typedef struct ECMT_t {		//  Mean EC pretrigger thresholds 
	int	in_high;        //  inner high threshold (mV)
	int	out_high;       //  outer high threshold (mV)
	int	total_high;     //  total high threshold (mV)
	int	in_lo;          //  inner lo threshold (mV)
	int	out_lo;         //  outer lo threshold (mV)
	int	total_lo;       //  total lo threshold (mV)
} ECMT_t;
typedef struct ECP1_t {		//  Translated LAC pedestal bank, measured during the pedestal run.  
	int	ID;             //  the address of the hit detector element
	int	mean_left;      //  left adc pedestal mean value (channel)
	float	sigma_left;     //  sigma of the pedestal distribution for left adc (channel)
	int	mean_right;     //  right adc pedestal mean value (channel)
	float	sigma_right;    //  sigma of the pedestal distribution for right adc (channel
} ECP1_t;
typedef struct ECPB_t {		//  EC hits involved in the event 
	int	ScHt;           //  100*sector+Whole_Hit_ID in ECHB 
	float	Etot;           //  Reconstructed total energy
	float	Ein;            //  Inner energy
	float	Eout;           //  Outer energy 
	float	Time;           //  Flight time relative to the evnt start time
	float	Path;           //  Path lenght from target
	float	X;              //  x coordinate of hit 
	float	Y;              //  y coordinate of hit
	float	Z;              //  z coordinate of hit
	float	M2_hit;         //  second moment of _hit pattern
	float	M3_hit;         //  third moment of  _hit pattern 
	float	M4_hit;         //  forth moment of  _hit pattern 
	int	InnStr;         //  10000*UI+100*VI+WI 
	int	OutStr;         //  10000*UO+100*VO+WO 
	float	Chi2EC;         //  Quality measure of geometrical matching
	int	Status;         //  Status word (not implemented yet)
} ECPB_t;
typedef struct ECPC_t {		//  EC Particle Calibration bank 
	int	ID;             //  layer(1-3) * 100 + strip ID
	float	TDC;            //  tdc (channels)
	float	ADC;            //  adc - pedestal (channels)
} ECPC_t;
typedef struct ECP_t {		//  
	float	N1x;            //  x component of outward normal to plane
	float	N1y;            //  y component of outward normal to plane
	float	N1z;            //  z component of outward normal to plane
	float	R1n;            //  distanse in cm from origin to plane
	float	PlW;            //  whole plane depth
	float	PlI;            //  inner plane depth
	float	PlO;            //  outer plane depth
} ECP_t;
typedef struct ECPE_t {		//  Translated EC pedestal bank, measured during the pedestal run. 
	int	ID;             //  the address of the hit detector element
	int	mean;           //  adc pedestal mean value (channel)
	float	sigma;          //  sigma of the pedestal distribution (channel)
} ECPE_t;
typedef struct ECPI_t {		//  Calorimeter pixels for DISPLAY only 
	int	ID;             //  sector number
	int	Layer;          //  layer number, note 1-inner, 2-outer, 3-whole
	int	HITID;          //  number of hits (first 16 bits) and hit ID (last 16 bits)
	float	iloc;           //  position of the hit in the local coordinate system
	float	jloc;           //  position of the hit in the local coordinate system
	float	diloc;          //  i width of the hit
	float	djloc;          //  j width of the hit 
	float	R;              //  radius of the shower
	float	E;              //  energy
} ECPI_t;
typedef struct ECPO_t {		//  Forward EC hit pointerw 
	int	UID1;           //  Pointer to the U-strips 1-18 for hit
	int	UID2;           //  Pointer to the U-strips 19-3 for hit
	int	VID1;           //  Pointer to the V-strips 1-18 for hit
	int	VID2;           //  Pointer to the V-strips 19-3 for hit
	int	WID1;           //  Pointer to the W-strips 1-18 for hit 
	int	WID2;           //  Pointer to the W-strips 19-36 for hit 
	int	SLH;            //  Sector*1000+Layer*100+Hit
} ECPO_t;
typedef struct ECRB_t {		//  Forward calorimeter result bank 
	float	E_in;           //  energy found for the inner layer
	float	E_out;          //  energy found for the outer layer
	float	dE_in;          //  error on the energy found for the inner layer
	float	dE_out;         //  error on the energy found for the outer layer
	float	t_in;           //  time found for the inner layer
	float	t_out;          //  time found for the outer layer
	float	dt_in;          //  error on the time found for the inner layer
	float	dt_out;         //  error on the time found for the outer layer
	float	i_in;           //  sector rectangular coordinate for the inner layer
	float	j_in;           //  sector rectangular coordinate for the inner layer
	float	i_out;          //  sector rectangular coordinate for the outer layer
	float	j_out;          //  sector rectangular coordinate for the outer layer
	float	di_in;          //  sector rectangular coordinate error, inner layer
	float	dj_in;          //  sector rectangular coordinate error, inner layer
	float	di_out;         //  sector rectangular coordinate error, outer layer
	float	dj_out;         //  sector rectangular coordinate error, outer layer
	float	x_in;           //  lab coordinate, inner layer
	float	y_in;           //  lab coordinate, inner layer
	float	z_in;           //  lab coordinate, inner layer
	float	x_out;          //  lab coordinate, outer layer
	float	y_out;          //  lab coordinate, outer layer
	float	z_out;          //  lab coordinate, outer layer
	float	dx_in;          //  lab coordinate error, inner layer
	float	dy_in;          //  lab coordinate error, inner layer
	float	dz_in;          //  lab coordinate error, inner layer
	float	dx_out;         //  lab coordinate error, outer layer
	float	dy_out;         //  lab coordinate error, outer layer
	float	dz_out;         //  lab coordinate error, outer layer
	float	u2_in;          //  second moment of u inner hit pattern
	float	v2_in;          //  second moment of v inner hit pattern
	float	w2_in;          //  second moment of w inner hit pattern
	float	u2_out;         //  second moment of u outer hit pattern
	float	v2_out;         //  second moment of v outer hit pattern
	float	w2_out;         //  second moment of w outer hit pattern
	float	u3_in;          //  third moment of u inner hit pattern 
	float	v3_in;          //  third moment of v inner hit pattern 
	float	w3_in;          //  third moment of w inner hit pattern 
	float	u3_out;         //  third moment of u outer hit pattern 
	float	v3_out;         //  third moment of v outer hit pattern 
	float	w3_out;         //  third moment of w outer hit pattern 
	float	i2;             //  second moment of overall shower, sector coordinates
	float	j2;             //  second moment of overall shower, sector coordinates
	float	i3;             //  third moment of overall shower, sector coordinates
	float	j3;             //  third moment of overall shower, sector coordinates
	float	spare1;         //  spare
	float	spare2;         //  spare
	float	spare3;         //  spare
	float	spare4;         //  spare
	int	istat;          //  status word
} ECRB_t;
typedef struct ECS_t {		//  Calorimeter scaler bank 
	int	ECS1;           //  scaler sector 1
	int	ECS2;           //  scaler sector 1
	int	ECS3;           //  scaler sector 1
	int	ECS4;           //  scaler sector 1
	int	ECS5;           //  scaler sector 1
	int	ECS6;           //  scaler sector 1
	int	ECS7;           //  scaler sector 1
	int	ECS8;           //  scaler sector 1
	int	ECS9;           //  scaler sector 1
	int	ECS10;          //  scaler sector 1
	int	ECS11;          //  scaler sector 1
	int	ECS12;          //  scaler sector 1
	int	ECS13;          //  scaler sector 1
	int	ECS14;          //  scaler sector 1
	int	ECS15;          //  scaler sector 1
	int	ECS16;          //  scaler sector 1
	int	ECS17;          //  scaler sector 2
	int	ECS18;          //  scaler sector 2
	int	ECS19;          //  scaler sector 2
	int	ECS20;          //  scaler sector 2
	int	ECS21;          //  scaler sector 2
	int	ECS22;          //  scaler sector 2
	int	ECS23;          //  scaler sector 2
	int	ECS24;          //  scaler sector 2
	int	ECS25;          //  scaler sector 2
	int	ECS26;          //  scaler sector 2
	int	ECS27;          //  scaler sector 2
	int	ECS28;          //  scaler sector 2
	int	ECS29;          //  scaler sector 2
	int	ECS30;          //  scaler sector 2
	int	ECS31;          //  scaler sector 2
	int	ECS32;          //  scaler sector 2
	int	ECS33;          //  scaler sector 3
	int	ECS34;          //  scaler sector 3
	int	ECS35;          //  scaler sector 3
	int	ECS36;          //  scaler sector 3
	int	ECS37;          //  scaler sector 3
	int	ECS38;          //  scaler sector 3
	int	ECS39;          //  scaler sector 3
	int	ECS40;          //  scaler sector 3
	int	ECS41;          //  scaler sector 3
	int	ECS42;          //  scaler sector 3
	int	ECS43;          //  scaler sector 3
	int	ECS44;          //  scaler sector 3
	int	ECS45;          //  scaler sector 3
	int	ECS46;          //  scaler sector 3
	int	ECS47;          //  scaler sector 3
	int	ECS48;          //  scaler sector 3
	int	ECS49;          //  scaler sector 4
	int	ECS50;          //  scaler sector 4
	int	ECS51;          //  scaler sector 4
	int	ECS52;          //  scaler sector 4
	int	ECS53;          //  scaler sector 4
	int	ECS54;          //  scaler sector 4
	int	ECS55;          //  scaler sector 4
	int	ECS56;          //  scaler sector 4
	int	ECS57;          //  scaler sector 4
	int	ECS58;          //  scaler sector 4
	int	ECS59;          //  scaler sector 4
	int	ECS60;          //  scaler sector 4
	int	ECS61;          //  scaler sector 4
	int	ECS62;          //  scaler sector 4
	int	ECS63;          //  scaler sector 4
	int	ECS64;          //  scaler sector 4
	int	ECS65;          //  scaler sector 5
	int	ECS66;          //  scaler sector 5
	int	ECS67;          //  scaler sector 5
	int	ECS68;          //  scaler sector 5
	int	ECS69;          //  scaler sector 5
	int	ECS70;          //  scaler sector 5
	int	ECS71;          //  scaler sector 5
	int	ECS72;          //  scaler sector 5
	int	ECS73;          //  scaler sector 5
	int	ECS74;          //  scaler sector 5
	int	ECS75;          //  scaler sector 5
	int	ECS76;          //  scaler sector 5
	int	ECS77;          //  scaler sector 5
	int	ECS78;          //  scaler sector 5
	int	ECS79;          //  scaler sector 5
	int	ECS80;          //  scaler sector 5
	int	ECS81;          //  scaler sector 6
	int	ECS82;          //  scaler sector 6
	int	ECS83;          //  scaler sector 6
	int	ECS84;          //  scaler sector 6
	int	ECS85;          //  scaler sector 6
	int	ECS86;          //  scaler sector 6
	int	ECS87;          //  scaler sector 6
	int	ECS88;          //  scaler sector 6
	int	ECS89;          //  scaler sector 6
	int	ECS90;          //  scaler sector 6
	int	ECS91;          //  scaler sector 6
	int	ECS92;          //  scaler sector 6
	int	ECS93;          //  scaler sector 6
	int	ECS94;          //  scaler sector 6
	int	ECS95;          //  scaler sector 6
	int	ECS96;          //  scaler sector 6
} ECS_t;
typedef struct ECT_t {		//  Forward calorimeter event bank (TDCs only) 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
} ECT_t;
typedef struct EID0_t {		//  Fast electron identification 
	int	jeid0sec;       //   Sector # for electron candidate
	int	jeid0cc;        //   Pointer to the hit in CCRC bank
	int	jeid0ec;        //   Pointer to the hit in ECHB bank
	int	jeid0sc;        //   Pointer to the hit in SCRX bank
} EID0_t;
typedef struct EPIC_t {		//  
	float	value;          //  value of the epics channel
	char	char1[32];      
} EPIC_t;
typedef struct EVNT_t {		//  RECSIS reconstructed event bank 
	int	ID;             //  Particle Data Group ID (from SEB) 
	float	Pmom;           //  momentum (from tracking)
	float	Mass;           //  mass squared (from SEB =p**2(1.-beta**2)/beta**2) 
	int	Charge;         //  charge (from tracking)
	float	Beta;           //  Particle velocity in the units of c (=R_trk/TOF/c)
	float	cx;             //  x dir cosine at track origin 
	float	cy;             //  y dir cosine at track origin
	float	cz;             //  z dir cosine at track origin
	float	X;              //  X coordinate of vertex (cm)
	float	Y;              //  Y coordinate of vertex (cm)
	float	Z;              //  Z coordinate of vertex (cm)
	int	DCstat;         //  Pointer to DCPB bank (=0 if DC is not involved) 
	int	CCstat;         //  Pointer to CCPB bank (=0 if CC is not involved) 
	int	SCstat;         //  Pointer to SCPB bank (=0 if SC is not involved)  
	int	ECstat;         //  Pointer to ECPB bank (=0 if EC is not involved) 
	int	LCstat;         //  Pointer to LCPB bank (=0 if LAC is not involved) 
	int	STstat;         //  Pointer to STPB bank (=0 if ST is not involved) 
	int	Status;         //  Status word (=0 out of time particle)
} EVNT_t;
typedef struct FBPM_t {		//  Fast BPM bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} FBPM_t;
typedef struct G1SL_t {		//  Gated1 T Photon Flux Scaler bank  
	int	G1T1;           
	int	G1T2;           
	int	G1T3;           
	int	G1T4;           
	int	G1T5;           
	int	G1T6;           
	int	G1T7;           
	int	G1T8;           
	int	G1T9;           
	int	G1T10;          
	int	G1T11;          
	int	G1T12;          
	int	G1T13;          
	int	G1T14;          
	int	G1T15;          
	int	G1T16;          
	int	G1T17;          
	int	G1T18;          
	int	G1T19;          
	int	G1T20;          
	int	G1T21;          
	int	G1T22;          
	int	G1T23;          
	int	G1T24;          
	int	G1T25;          
	int	G1T26;          
	int	G1T27;          
	int	G1T28;          
	int	G1T29;          
	int	G1T30;          
	int	G1T31;          
	int	G1T32;          
	int	G1T33;          
	int	G1T34;          
	int	G1T35;          
	int	G1T36;          
	int	G1T37;          
	int	G1T38;          
	int	G1T39;          
	int	G1T40;          
	int	G1T41;          
	int	G1T42;          
	int	G1T43;          
	int	G1T44;          
	int	G1T45;          
	int	G1T46;          
	int	G1T47;          
	int	G1T48;          
	int	G1T49;          
	int	G1T50;          
	int	G1T51;          
	int	G1T52;          
	int	G1T53;          
	int	G1T54;          
	int	G1T55;          
	int	G1T56;          
	int	G1T57;          
	int	G1T58;          
	int	G1T59;          
	int	G1T60;          
	int	G1T61;          
	int	G1T62;          
	int	G1T63;          
	int	G1T64;          
} G1SL_t;
typedef struct G2SL_t {		//  Gated2 T Photon Flux Scaler bank  
	int	G2T1;           
	int	G2T2;           
	int	G2T3;           
	int	G2T4;           
	int	G2T5;           
	int	G2T6;           
	int	G2T7;           
	int	G2T8;           
	int	G2T9;           
	int	G2T10;          
	int	G2T11;          
	int	G2T12;          
	int	G2T13;          
	int	G2T14;          
	int	G2T15;          
	int	G2T16;          
	int	G2T17;          
	int	G2T18;          
	int	G2T19;          
	int	G2T20;          
	int	G2T21;          
	int	G2T22;          
	int	G2T23;          
	int	G2T24;          
	int	G2T25;          
	int	G2T26;          
	int	G2T27;          
	int	G2T28;          
	int	G2T29;          
	int	G2T30;          
	int	G2T31;          
	int	G2T32;          
	int	G2T33;          
	int	G2T34;          
	int	G2T35;          
	int	G2T36;          
	int	G2T37;          
	int	G2T38;          
	int	G2T39;          
	int	G2T40;          
	int	G2T41;          
	int	G2T42;          
	int	G2T43;          
	int	G2T44;          
	int	G2T45;          
	int	G2T46;          
	int	G2T47;          
	int	G2T48;          
	int	G2T49;          
	int	G2T50;          
	int	G2T51;          
	int	G2T52;          
	int	G2T53;          
	int	G2T54;          
	int	G2T55;          
	int	G2T56;          
	int	G2T57;          
	int	G2T58;          
	int	G2T59;          
	int	G2T60;          
	int	G2T61;          
	int	G2T62;          
	int	G2T63;          
	int	G2T64;          
} G2SL_t;
typedef struct G3SL_t {		//  Gated 3 T Photon Flux Scaler bank  
	int	G3T1;           
	int	G3T2;           
	int	G3T3;           
	int	G3T4;           
	int	G3T5;           
	int	G3T6;           
	int	G3T7;           
	int	G3T8;           
	int	G3T9;           
	int	G3T10;          
	int	G3T11;          
	int	G3T12;          
	int	G3T13;          
	int	G3T14;          
	int	G3T15;          
	int	G3T16;          
	int	G3T17;          
	int	G3T18;          
	int	G3T19;          
	int	G3T20;          
	int	G3T21;          
	int	G3T22;          
	int	G3T23;          
	int	G3T24;          
	int	G3T25;          
	int	G3T26;          
	int	G3T27;          
	int	G3T28;          
	int	G3T29;          
	int	G3T30;          
	int	G3T31;          
	int	G3T32;          
	int	G3T33;          
	int	G3T34;          
	int	G3T35;          
	int	G3T36;          
	int	G3T37;          
	int	G3T38;          
	int	G3T39;          
	int	G3T40;          
	int	G3T41;          
	int	G3T42;          
	int	G3T43;          
	int	G3T44;          
	int	G3T45;          
	int	G3T46;          
	int	G3T47;          
	int	G3T48;          
	int	G3T49;          
	int	G3T50;          
	int	G3T51;          
	int	G3T52;          
	int	G3T53;          
	int	G3T54;          
	int	G3T55;          
	int	G3T56;          
	int	G3T57;          
	int	G3T58;          
	int	G3T59;          
	int	G3T60;          
	int	G3T61;          
	int	G3T62;          
	int	G3T63;          
	int	G3T64;          
} G3SL_t;
typedef struct G4SL_t {		//  Gated 4 T Photon Flux Scaler bank  
	int	G4T1;           
	int	G4T2;           
	int	G4T3;           
	int	G4T4;           
	int	G4T5;           
	int	G4T6;           
	int	G4T7;           
	int	G4T8;           
	int	G4T9;           
	int	G4T10;          
	int	G4T11;          
	int	G4T12;          
	int	G4T13;          
	int	G4T14;          
	int	G4T15;          
	int	G4T16;          
	int	G4T17;          
	int	G4T18;          
	int	G4T19;          
	int	G4T20;          
	int	G4T21;          
	int	G4T22;          
	int	G4T23;          
	int	G4T24;          
	int	G4T25;          
	int	G4T26;          
	int	G4T27;          
	int	G4T28;          
	int	G4T29;          
	int	G4T30;          
	int	G4T31;          
	int	G4T32;          
	int	G4T33;          
	int	G4T34;          
	int	G4T35;          
	int	G4T36;          
	int	G4T37;          
	int	G4T38;          
	int	G4T39;          
	int	G4T40;          
	int	G4T41;          
	int	G4T42;          
	int	G4T43;          
	int	G4T44;          
	int	G4T45;          
	int	G4T46;          
	int	G4T47;          
	int	G4T48;          
	int	G4T49;          
	int	G4T50;          
	int	G4T51;          
	int	G4T52;          
	int	G4T53;          
	int	G4T54;          
	int	G4T55;          
	int	G4T56;          
	int	G4T57;          
	int	G4T58;          
	int	G4T59;          
	int	G4T60;          
	int	G4T61;          
	int	G4T62;          
	int	G4T63;          
	int	G4T64;          
} G4SL_t;
typedef struct GPAR_t {		//   input parameter bank  
	float	fval;           //  float value
	int	ival;           //  int value 
	char	char1[20];      
} GPAR_t;
typedef struct GPID_t {		//  ST-SC PID bank 
	int	pid;            //  particle id (GEANT)
	float	x;              //  vector3_t vert; Vertex position {x,y,z} 
	float	y;              //  y
	float	z;              //  z
	float	E;              //  vector4_t p; Energy
	float	px;             //  momentum {x,y,z}
	float	py;             //   py
	float	pz;             //  pz
	int	q;              //  charge
	int	trkid;          //  index to TBID bank, counting from 1
	int	sec;            //  Sector track is in
	int	paddle;         //  Paddle hit
	float	dedx;           //  Energy deposited in TOF
	float	beta;           //  beta pmag/E
	int	sc_stat;        //  status of hit matching to SC: see sc.h
	float	sc_time;        //  SC calibrated time for this track (ns)
	float	sc_len;         //  track length [cm] from origin to SC
	int	st_stat;        //  ST status
	float	st_time;        //  ST calibrated time for this track (ns) 
	float	st_len;         //  track length [cm] from origin to ST
	float	mass;           //  mass from time-of-flight
	int	mass_ref;       //  0: mass calc from SC & TAG, 1: SC & ST; -1 neutral or no SC;2:from PART
	float	betam;          //  beta from time_of-flight
	float	epho;           //  closest photon energy (GeV)
	float	tpho;           //  Time of the photon after RF correction
	int	tagrid;         //  index to TAGR bank, counting from 1
	int	ngrf;           //  number of photons in the same RF bucket
	int	ppid;           //  pid as seen in PART bank
} GPID_t;
typedef struct GP_X_t {		//  
	float	value;          //  count rate per channel
} GP_X_t;
typedef struct GP_Y_t {		//  
	float	value;          //  count rate per channel
} GP_Y_t;
typedef struct HBER_t {		//  Hit Based Tracking ERror bank 
	float	q_over_p;       //  q/p
	float	lambda;         //  dip angle (pi/2 - theta)
	float	phi;            //  phi
	float	d0;             //  min.distance from (x=0,y=0,z=?)   [cm]
	float	z0;             //  z position of starting point  [cm]
	float	c11;            //  element C{1,1}
	float	c12;            //  element C{1,2}
	float	c13;            //  element C{1,3}
	float	c14;            //  element C{1,4}
	float	c15;            //  element C{1,5}
	float	c22;            //  element C{2,2}
	float	c23;            //  element C{2,3}
	float	c24;            //  element C{2,4}
	float	c25;            //  element C{2,5}
	float	c33;            //  element C{3,3}
	float	c34;            //  element C{3,4}
	float	c35;            //  element C{3,5}
	float	c44;            //  element C{4,4}
	float	c45;            //  element C{4,5}
	float	c55;            //  element C{5,5}
	float	chi2;           //  Chisquare for this Track
	int	layinfo1;       //  layerhit info
	int	layinfo2;       //  layerhit info&sector&track#in sector
} HBER_t;
typedef struct HBID_t {		//  Hit Based particle ID 
	int	track;          //  Track Candidate (ptr to HBTR)
	int	sec;            //  Sector track is in
	float	beta;           //  Beta of the track in units of c
	float	vtime;          //  vertex time of track
	int	sc_stat;        //  status of hit matching to SC: see sc.h
	int	sc_id;          //  Pointer to SCRC bank
	float	sc_time;        //  SC calibrated time for this track (ns)
	float	sc_qual;        //  quality of match for SC
	float	sc_vtime;       //  time at vertex for SC(ns)
	float	sc_beta;        //  Beta calculated from TOF from SC
	int	cc_stat;        //  status of hit matching to CC: see sc.h
	int	cc_id;          //  pointer to CC01 bank
	float	cc_time;        //  CC calibrated time for this track (ns)
	float	cc_qual;        //  quality of match for CC
	float	cc_vtime;       //  time at vertex for CC(ns)
	float	cc_beta;        //  Beta as calculated by the EC
	int	ec_stat;        //  status of hit matching to ec: see sc.h
	int	ec_id;          //  Pointer to ECHB bank
	float	ec_time;        //  EC calibrated time for this track (ns)
	float	ec_qual;        //  EC quality factor
	float	ec_vtime;       //  time at vertex for EC(ns)
	float	ec_beta;        //  Beta as calculated by the EC
	int	st_stat;        //  status of hit matching to ST
	int	st_id;          //  Pointer to STR bank
	float	st_time;        //  ST calibrated time for this track (ns)
	float	st_qual;        //  ST quality factor
	float	st_vtime;       //  time at vertex for ST(ns)
	float	st_beta;        //  Beta as calculated by the ST
	int	lac_stat;       //  status of hit matching to LAC
	int	lac_id;         //  Pointer to EC1R bank
	float	lac_time;       //  LAC calibrated time for this track (ns)
	float	lac_qual;       //  LAC quality factor
	float	lac_vtime;      //  time at vertex for LAC(ns)
	float	lac_beta;       //  Beta as calculated by the LAC
} HBID_t;
typedef struct HBLA_t {		//  Hit Based tracking LAyer bank 
	int	trk_pln;        //  (track_number) *100 + Trk_plane_number
	float	x;              //  z coord [cm]  for track in this plane
	float	y;              //  y coord [cm]  for track in this plane
	float	z;              //  z coord [cm]  for track in this plane
	float	Bx;             //  B-field in x [kG] at coord.{x,y,z}
	float	By;             //  B-field in y [kG] at coord.{x,y,z}
	float	Bz;             //  B-field in z [kG] at coord.{x,y,z}
	float	tlen;           //  track length [cm] from origin to this plane
	int	dc1;            //  Pointer to DC1 bank
	int	stat;           //  Status of the hit
	int	wire;           //  Wire number 
	float	dtime;          //  drift time [ns]  (not corrected for start ime)
	float	alpha;          //  track angle (relative to R of SL) [deg]
	float	wlen;           //  Wire length (hit pos. to preamp)  [cm]
	float	sgdoca;         //  sigma DOCA
	float	fitdoca;        //  Fitted DOCA
} HBLA_t;
typedef struct HBTB_t {		//  Hit Based to Time Based tracking 
	int	trk;            //  Track Candidate 
	int	sly;            //  Superlayer
	int	icl;            //  Cluster# that matched to the track
	int	isg;            //  Segment# in this cluster
	int	id1;            //  Ptr. to DC1 entry  for 1.Layer in SL
	int	id2;            //  Ptr. to DC1 entry  for 2.Layer in SL
	int	id3;            //  Ptr. to DC1 entry  for 3.Layer in SL
	int	id4;            //  Ptr. to DC1 entry  for 4.Layer in SL
	int	id5;            //  Ptr. to DC1 entry  for 5.Layer in SL
	int	id6;            //  Ptr. to DC1 entry  for 6.Layer in SL
} HBTB_t;
typedef struct HBTR_t {		//  Hit Based Tracking Result bank 
	float	x;              //  x
	float	y;              //  y 'Vertex' position {x,y,z}
	float	z;              //  z
	float	px;             //  Px
	float	py;             //  Py  momentum at 'vertex' {Px,Py,Pz}
	float	pz;             //  Pz
	float	q;              //  charge   (straight tracks: set q=0.)
	float	chi2;           //  Chisquare for this Track
	int	itr_sec;        //  Trackno_in_Sector + Sector*100
} HBTR_t;
typedef struct HCAL_t {		//  Header CALibration bank 
	int	VERSION;        //  Version number
	int	RUN_LOW;        //  LOW RUN number		 
	int	RUN_HIGH;       //  HIGH RUN number
	int	CATIME;         //  CAlibration TIME (unix time = 0 1/1/1970)		
	int	ROCCA;          //  32 bit ReadOut Controller CAlibration status
} HCAL_t;
typedef struct HDPL_t {		//  Hbt Detector PLane bank 
	int	trk_pln;        //  (track_number) *100 + plane_number
	float	x;              //  vector3_t pos;
	float	y;              //   (x, y, z coord. for track in on the plane)
	float	z;              
	float	cx;             //  vector3_t dir;
	float	cy;             //   direction cosines (x,y,z) for track at coord.{x,y,z}
	float	cz;             
	float	tlen;           //  track length [cm] from origin to this plane
} HDPL_t;
typedef struct HEAD_t {		//  Bank HEAD 
	int	VERSION;        //  Version Number
	int	NRUN;           //  Run Number (monotonically increasing)
	int	NEVENT;         //  Event Number (starting with 1 at run begin,
	int	TIME;           //  Event Time (UNIX time = seconds as of January 1,1970)
	int	TYPE;           //  Event Type (Defined by on-line system or MC run:
	int	ROC;            //             = 0        - sync status is OK
	int	EVTCLASS;       //  Event Classification from DAQ:
	int	TRIGBITS;       //  Level 1 Trigger Latch Word (16 bits)
} HEAD_t;
typedef struct HEVT_t {		// Reconstructed event, Head bank 
	int	ESTATUS;        //  Event Statuse after reconstruction
	int	NRUN;           //  Run Number (monotonically increasing)
	int	NEVENT;         //  Event Number in the run NRUN 
	int	TYPE;           //  Event Type (Data or MC)
	int	NPGP;           //  Number of final reconstructed particles*100 +
	int	TRGPRS;         //  Trigger type*10000 + 
	float	FC;             //  Faraday Cup (K)
	float	FCG;            //  Faraday Cup Gated (K)
	float	TG;             //  Clock Gated
	float	STT;            //  Event Start Time 
	float	RF1;            //  RF Time 
	float	RF2;            //  RF Time
	float	CON1;           //  Control Rates 
	float	CON2;           //  Control Rates  
	float	CON3;           //  Control Rates 
	int	PTIME;          //  Event Processing Time (UNIX time = seconds)
} HEVT_t;
typedef struct HLS_t {		//  Helicity scaler bank 
	unsigned int	S1;             //  10MHz clock
	unsigned int	S2;             //  OTR-1
	unsigned int	S3;             //  OTR-2
	unsigned int	S4;             //  SLM
	unsigned int	S5;             //  lvl-1 trigger rate
	unsigned int	S6;             //  L.R Moller coincidences
	unsigned int	S7;             //  L.R Moller accidentals
	unsigned int	S8;             //  F-CUP
	unsigned int	S9;             //  pmt-1
	unsigned int	S10;            //  pmt-2
	unsigned int	S11;            //  pmt-3
	unsigned int	S12;            //  pmt-4
	unsigned int	S13;            //  reserve
	unsigned int	S14;            //  reserve
	unsigned int	S15;            //  Helicity states accumulating counter
	unsigned int	S16;            //  HLS banks accumulating counter
} HLS_t;
typedef struct IC_t {		//  Inner calorimeter event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} IC_t;
typedef struct ICHB_t {		//  Inner crystal calorimeter reconstruction bank 
	float	Eclust;         //  Cluster energy 
	float	Eclmax;         //  Max energy in a crystal of the cluster 
	float	Tclust;         //  Reconstructed cluster time
	float	Tclmax;         //  Time in a crystal with Max enrgy
	float	xclust;         //  lab coordinate X,
	float	yclust;         //  lab coordinate Y,
	float	zclust;         //  lab coordinate Z,
	float	xclmax;         //  lab coordinate X,
	float	yclmax;         //  lab coordinate Y,
	float	zclmax;         //  lab coordinate Z,
	float	res1;           //  lab coordinate error,
	float	res2;           //  lab coordinate error,
	float	res3;           //  lab coordinate error,
	int	ncryst;         //  Number of crystal in the cluster
} ICHB_t;
typedef struct KFIT_t {		// Kinematic Fitter Bank 
	float	momenta_f;      //  fitted momentum parameter (GeV)   
	float	lambda_f;       //  fitted dip angle (radians) 
	float	phi_f;          //  fitted phi angle (radians)
	float	d0_f;           //  fitted d0 (cm)
	float	z0_f;           //  fitted z0 (cm)
	float	sigma11;        //    error element of known
	float	sigma12;        //    error element of known
	float	sigma13;        //    error element of known
	float	sigma14;        //    error element of known
	float	sigma15;        //    error element of known
	float	sigma22;        //    error element of known   
	float	sigma23;        //    error element of known
	float	sigma24;        //    error element of known
	float	sigma25;        //    error element of known
	float	sigma33;        //    error element of known 
	float	sigma34;        //    error element of known 
	float	sigma35;        //    error element of known 
	float	sigma44;        //    error element of known 
	float	sigma45;        //    error element of known 
	float	sigma55;        //    error element of known 
	float	cov11;          //    error element bewteen known and unknown
	float	cov12;          //    error element bewteen known and unknown
	float	cov13;          //    error element bewteen known and unknown
	float	cov14;          //    error element bewteen known and unknown
	float	cov15;          //    error element bewteen known and unknown
	float	cov21;          //    error element bewteen known and unknown
	float	cov22;          //    error element bewteen known and unknown
	float	cov23;          //    error element bewteen known and unknown
	float	cov24;          //    error element bewteen known and unknown
	float	cov25;          //    error element bewteen known and unknown
	float	cov31;          //    error element bewteen known and unknown
	float	cov32;          //    error element bewteen known and unknown
	float	cov33;          //    error element bewteen known and unknown
	float	cov34;          //    error element bewteen known and unknown
	float	cov35;          //    error element bewteen known and unknown
	float	cov41;          //    error element bewteen known and unknown
	float	cov42;          //    error element bewteen known and unknown
	float	cov43;          //    error element bewteen known and unknown
	float	cov44;          //    error element bewteen known and unknown
	float	cov45;          //    error element bewteen known and unknown
	float	cov51;          //    error element bewteen known and unknown
	float	cov52;          //    error element bewteen known and unknown
	float	cov53;          //    error element bewteen known and unknown
	float	cov54;          //    error element bewteen known and unknown
	float	cov55;          //    error element bewteen known and unknown
	float	chi2piece;      //  track contribution to chi2
} KFIT_t;
typedef struct L1PG_t {		//  Level 1 program 
	char	l1prog[4];      
} L1PG_t;
typedef struct L2H_t {		//  Level 2 hit bank 
	unsigned short	ID;             //  Superlayer + 256*(hit# in TDC)
	unsigned short	TDC;            //  tdc information (channels)
} L2H_t;
typedef struct L2S_t {		//  Level 2 report bank 
	unsigned short	ID;             //  ID+256*hit#
	unsigned short	TDC;            //  TDC value
} L2S_t;
typedef struct LASR_t {		//  events from laser calibration of the EC 
	int	ID;             //  identifier  0-1099  TOF 1100-10000 EC laser  
	int	stat1;          //  status words different for each id
	int	stat2;          //   "
	int	stat3;          //   "
} LASR_t;
typedef struct LCDI_t {		//  LAC channel discriminator thresholds 
	int	crate;          //  CAMAC crate number
	int	slot;           //  slot
	int	mask;           //  mask
	int	threshold;      //  actual threshold value (mV)
} LCDI_t;
typedef struct LCPB_t {		//  LAC hits involved in the event 
	int	ScHt;           //  100*sector+Hit_ID in EC1R
	float	Etot;           //  Reconstructed total energy
	float	Time;           //  Flight time relative to the evnt start time
	float	Path;           //  Path lenght from target
	float	X;              //  x coordinate of the hit
	float	Y;              //  y coordinate of the hit
	float	Z;              //  z coordinate of the hit
	float	Chi2LC;         //  Quality measure of geometrical matching
	int	Status;         //  Status word
	float	Ein;            //  Reconstructed energy in the inner part
} LCPB_t;
typedef struct LOWQ_t {		//  FET reconstruction 
	int	pid;            //  particle id (GEANT)
	float	x;              //  vector3_t vert; Vertex position {x,y,z}
	float	y;              //  y
	float	z;              //  z
	float	E;              //  vector4_t p; Energy
	float	px;             //  momentum {x,y,z}
	float	py;             //   py
	float	pz;             //  pz
	float	q;              //  charge
	float	qpid;           //  quality factor for the pid
	float	qtrk;           //  quality factor for the trk
	int	flags;          //  set of flags defining track (ie, BEAM)
} LOWQ_t;
typedef struct MCEV_t {		//  GSIM Monte Carlo event parameters 
	int	I1;             //  first geant random number seed for event
	int	i2;             //  second seed
} MCEV_t;
typedef struct MCHD_t {		//  MC (generator) HeaD Bank 
	int	NRUN;           //  Run Number
	int	NEVENT;         //  Event Number
	int	TIME;           //  Event Time (UNIX time)
	int	TYPE;           //  Event Type (MUST be negative)
	int	REACTYPE;       //  optional: user defined reaction type
	float	WEIGHT;         //  event weight
	float	W;              //  center_of_mass energy (inv.mass of hadronic states)
	float	Q2;             //  Q2 (photon virtuallity)
	float	E_PHOT;         //  energy of (real/virtual) photon
	float	PX_PHOT;        //  p_x: momentum (in x) of (real/virtual) photon
	float	PY_PHOT;        //  p_y: momentum (in y) of (real/virtual) photon
	float	PZ_PHOT;        //  p_z: momentum (in z) of (real/virtual) photon
	float	E_TARG;         //  energy of target particle
	float	PX_TARG;        //  p_x: momentum (in x) of target part.
	float	PY_TARG;        //  p_y: momentum (in y) of target part.
	float	PZ_TARG;        //  p_z: momentum (in z) of target part.
} MCHD_t;
typedef struct MCTK_t {		//  GSIM Monte Carlo track bank 
	float	cx;             //  x dir cosine at track origin
	float	cy;             //  y dir cosine
	float	cz;             //  z dir cosine
	float	pmom;           //  momentum
	float	mass;           //  mass
	float	charge;         //  charge
	int	id;             //  track Particle Data Group id
	int	flag;           //  track flag
	int	beg_vtx;        //  beginning vertex number 
	int	end_vtx;        //  ending vertex number
	int	parent;         //  parent track
} MCTK_t;
typedef struct MCVX_t {		//  GSIM Monte Carlo vertex parameters 
	float	x;              //  x of vertex
	float	y;              //  y
	float	z;              //  z
	float	tof;            //  secs of flight
	int	flag;           //  vertex flag
} MCVX_t;
typedef struct MS1_t {		//  Polarimeter data 
	unsigned short	ID;             //  the address of the struck microstrip/panel (4 panels)
	unsigned short	ADC;            //  adc information (channels)
} MS1_t;
typedef struct MTRK_t {		//  vertex Result bank 
	int	sect;           //  Sector Number for track
	int	trkl1;          //  track segment and track cluster for Superlayer 1 (lnk_t)
	int	trkl2;          //  track segment and track cluster for Superlayer 2
	int	trkl3;          //  track segment and track cluster for Superlayer 3
	int	trkl4;          //  track segment and track cluster for Superlayer 4
	int	trkl5;          //  track segment and track cluster for Superlayer 5
	int	trkl6;          //  track segment and track cluster for Superlayer 6
} MTRK_t;
typedef struct MVRT_t {		//  vertex Result bank 
	int	v_id;           //  info about track ids
	int	ntrk;           //  number of tracks used to make vertex
	float	x;              //  x vector3_t vert{x,y,z}
	float	y;              //  y 
	float	z;              //  z
	float	chi2;           //  chi2
	float	cxx;            //  Covariance matrix array element
	float	cxy;            //  Covariance matrix array element
	float	cxz;            //  Covariance matrix array element
	float	cyy;            //  Covariance matrix array element
	float	cyz;            //  Covariance matrix array element
	float	czz;            //  Covariance matrix array element
	int	stat;           //  status integer, not used yet
} MVRT_t;
typedef struct PART_t {		//  Hit Based Tracking Result bank 
	int	pid;            //  particle id (GEANT)
	float	x;              //  vector3_t vert; Vertex position {x,y,z} 
	float	y;              //  y
	float	z;              //  z
	float	E;              //  vector4_t p; Energy
	float	px;             //  momentum {x,y,z}
	float	py;             //   py
	float	pz;             //  pz
	float	q;              //  charge
	int	trkid;          //  index to TBID bank, counting from 1
	float	qpid;           //  quality factor for the pid
	float	qtrk;           //  quality factor for the trk
	int	flags;          //  set of flags defining track (ie, BEAM)
} PART_t;
typedef struct PCO_t {		//  down stream output bank 
	float	TIME;           //   TDC time 
	float	ELT;            //  energy deposit PC = left top
	float	ERB;            //  energy deposit PC = right bottom
	float	ELB;            //  energy deposit PC = left bottom 
	float	ERT;            //  energy deposit PC = right top
	float	EMAIN;          //  energy deposit PC = MAIN
	float	EVETO;          //  energy deposit PC = veto
	int	TID;            //  T id of the corresponding T
} PCO_t;
typedef struct PHTM_t {		//  Photon_mon bank for mySQL 
	unsigned int	EVID;           //  Event ID (number of triggers)
	unsigned int	Nsprot;         //  Number of single proton events
	unsigned int	Npip;           //  Number of single pi+ events
	unsigned int	Npim;           //  Number of single pi- events
	unsigned int	Nppippim;       //  Number of single proton pi+ pi- events
	unsigned int	Nppip;          //  Number of single proton pi+ events
	unsigned int	Npippim;        //  Number of single pi+ pi- events
	unsigned int	Nppim;          //  Number of single proton pi- events
	unsigned int	Nkp;            //  Number of single K+ events
	unsigned int	Npkp;           //  Number of single proton K+ events
	float	tag_pi_v;       //  Mean of (tagger_time - pion_vertex_time)
} PHTM_t;
typedef struct PID1_t {		//  
	unsigned int	NHBpos_1;       //  Number of pos Hit-Based tracks for sec=1
	unsigned int	NHBneg_1;       //  Number of neg Hit-Based tracks for sec=1
	unsigned int	NTBpos_1;       //  Number of pos Time-Based tracks  for sec=1
	unsigned int	NTBneg_1;       //  Number of neg Time-Based tracks for sec=1
	float	chi2pos_1;      //  Chi^2 for positive tracks for sec=1
	float	chi2neg_1;      //  Chi^2 for negative tracks for sec=1
	unsigned int	Nunknown_1;     //  Number of unknown particles for sec=1
	unsigned int	Ngamma_1;       //  Number of photons for sec=1
	unsigned int	Nelec_1;        //  Number of electrons for sec=1
	unsigned int	Npiplus_1;      //  Number of pi+ for sec=1
	unsigned int	Npiminus_1;     //  Number of pi- for sec=1
	unsigned int	Nprot_1;        //  Number of proton for sec=1
	unsigned int	Nkplus_1;       //  Number of K+ for sec=1
	unsigned int	Nkminus_1;      //  Number of K- for sec=1
	unsigned int	Ndeuteron_1;    //  Number of deuterons for sec=1
	unsigned int	Nneutron_1;     //  Number of neutrons for sec=1 
} PID1_t;
typedef struct PID2_t {		//  
	unsigned int	NHBpos_2;       //  Number of pos Hit-Based tracks for sec=2
	unsigned int	NHBneg_2;       //  Number of neg Hit-Based tracks for sec=2
	unsigned int	NTBpos_2;       //  Number of pos Time-Based tracks  for sec=2
	unsigned int	NTBneg_2;       //  Number of neg Time-Based tracks for sec=2
	float	chi2pos_2;      //  Chi^2 for positive tracks for sec=2
	float	chi2neg_2;      //  Chi^2 for negative tracks for sec=2
	unsigned int	Nunknown_2;     //  Number of unknown particles for sec=2
	unsigned int	Ngamma_2;       //  Number of photons for sec=2
	unsigned int	Nelec_2;        //  Number of electrons for sec=2
	unsigned int	Npiplus_2;      //  Number of pi+ for sec=2
	unsigned int	Npiminus_2;     //  Number of pi- for sec=2
	unsigned int	Nprot_2;        //  Number of proton for sec=2
	unsigned int	Nkplus_2;       //  Number of K+ for sec=2
	unsigned int	Nkminus_2;      //  Number of K- for sec=2
	unsigned int	Ndeuteron_2;    //  Number of deuterons for sec=2
	unsigned int	Nneutron_2;     //  Number of neutrons for sec=2
} PID2_t;
typedef struct PID3_t {		//  
	unsigned int	NHBpos_3;       //  Number of pos Hit-Based tracks for sec=3
	unsigned int	NHBneg_3;       //  Number of neg Hit-Based tracks for sec=3
	unsigned int	NTBpos_3;       //  Number of pos Time-Based tracks  for sec=3
	unsigned int	NTBneg_3;       //  Number of neg Time-Based tracks for sec=3
	float	chi2pos_3;      //  Chi^2 for positive tracks for sec=3
	float	chi2neg_3;      //  Chi^2 for negative tracks for sec=3
	unsigned int	Nunknown_3;     //  Number of unknown particles for sec=3
	unsigned int	Ngamma_3;       //  Number of photons for sec=3
	unsigned int	Nelec_3;        //  Number of electrons for sec=3
	unsigned int	Npiplus_3;      //  Number of pi+ for sec=3
	unsigned int	Npiminus_3;     //  Number of pi- for sec=3
	unsigned int	Nprot_3;        //  Number of proton for sec=3
	unsigned int	Nkplus_3;       //  Number of K+ for sec=3
	unsigned int	Nkminus_3;      //  Number of K- for sec=3
	unsigned int	Ndeuteron_3;    //  Number of deuterons for sec=3
	unsigned int	Nneutron_3;     //  Number of neutrons for sec=3
} PID3_t;
typedef struct PID4_t {		//  
	unsigned int	NHBpos_4;       //  Number of pos Hit-Based tracks for sec=4
	unsigned int	NHBneg_4;       //  Number of neg Hit-Based tracks for sec=4
	unsigned int	NTBpos_4;       //  Number of pos Time-Based tracks  for sec=4
	unsigned int	NTBneg_4;       //  Number of neg Time-Based tracks for sec=4
	float	chi2pos_4;      //  Chi^2 for positive tracks for sec=4
	float	chi2neg_4;      //  Chi^2 for negative tracks for sec=4
	unsigned int	Nunknown_4;     //  Number of unknown particles for sec=4
	unsigned int	Ngamma_4;       //  Number of photons for sec=4
	unsigned int	Nelec_4;        //  Number of electrons for sec=4
	unsigned int	Npiplus_4;      //  Number of pi+ for sec=4
	unsigned int	Npiminus_4;     //  Number of pi- for sec=4
	unsigned int	Nprot_4;        //  Number of proton for sec=4
	unsigned int	Nkplus_4;       //  Number of K+ for sec=4
	unsigned int	Nkminus_4;      //  Number of K- for sec=4
	unsigned int	Ndeuteron_4;    //  Number of deuterons for sec=4
	unsigned int	Nneutron_4;     //  Number of neutrons for sec=4   
} PID4_t;
typedef struct PID5_t {		//  
	unsigned int	NHBpos_5;       //  Number of pos Hit-Based tracks for sec=5
	unsigned int	NHBneg_5;       //  Number of neg Hit-Based tracks for sec=5
	unsigned int	NTBpos_5;       //  Number of pos Time-Based tracks  for sec=5
	unsigned int	NTBneg_5;       //  Number of neg Time-Based tracks for sec=5
	float	chi2pos_5;      //  Chi^2 for positive tracks for sec=5
	float	chi2neg_5;      //  Chi^2 for negative tracks for sec=5
	unsigned int	Nunknown_5;     //  Number of unknown particles for sec=5
	unsigned int	Ngamma_5;       //  Number of photons for sec=5
	unsigned int	Nelec_5;        //  Number of electrons for sec=5
	unsigned int	Npiplus_5;      //  Number of pi+ for sec=5
	unsigned int	Npiminus_5;     //  Number of pi- for sec=5
	unsigned int	Nprot_5;        //  Number of proton for sec=5
	unsigned int	Nkplus_5;       //  Number of K+ for sec=5
	unsigned int	Nkminus_5;      //  Number of K- for sec=5
	unsigned int	Ndeuteron_5;    //  Number of deuterons for sec=5
	unsigned int	Nneutron_5;     //  Number of neutrons for sec=5 
} PID5_t;
typedef struct PID6_t {		//  
	unsigned int	NHBpos_6;       //  Number of pos Hit-Based tracks for sec=6
	unsigned int	NHBneg_6;       //  Number of neg Hit-Based tracks for sec=6
	unsigned int	NTBpos_6;       //  Number of pos Time-Based tracks  for sec=6
	unsigned int	NTBneg_6;       //  Number of neg Time-Based tracks for sec=6
	float	chi2pos_6;      //  Chi^2 for positive tracks for sec=6
	float	chi2neg_6;      //  Chi^2 for negative tracks for sec=6
	unsigned int	Nunknown_6;     //  Number of unknown particles for sec=6
	unsigned int	Ngamma_6;       //  Number of photons for sec=6
	unsigned int	Nelec_6;        //  Number of electrons for sec=6
	unsigned int	Npiplus_6;      //  Number of pi+ for sec=6
	unsigned int	Npiminus_6;     //  Number of pi- for sec=6
	unsigned int	Nprot_6;        //  Number of proton for sec=6
	unsigned int	Nkplus_6;       //  Number of K+ for sec=6
	unsigned int	Nkminus_6;      //  Number of K- for sec=6
	unsigned int	Ndeuteron_6;    //  Number of deuterons for sec=6
	unsigned int	Nneutron_6;     //  Number of neutrons for sec=6 
} PID6_t;
typedef struct PIDT_t {		//  
	unsigned int	Nunknown;       //  Number of particles labelled unknown
	unsigned int	Ngamma;         //  Number of photons
	unsigned int	Nelec;          //  Number of electrons
	unsigned int	Npiplus;        //  Number of pi+ 
	unsigned int	Npiminus;       //  Number of pi-
	unsigned int	Nprot;          //  Number of proton
	unsigned int	Nkplus;         //  Number of K+
	unsigned int	Nkminus;        //  Number of K-
	unsigned int	Ndeuteron;      //  Number of deuterons
	unsigned int	Nneutron;       //  Number of neutrons
} PIDT_t;
typedef struct PRTM_t {		//  Pretrigger calibration time 
	int	time;           //  time of pretrig calibration
} PRTM_t;
typedef struct PSO_t {		//  down stream output bank 
	int	ID;             //  Counter Id ()
	float	TIME;           //   TDC time 
	float	ENER;           //  energy deposit
	int	TID;            //  T id of the corresponding T
} PSO_t;
typedef struct PTDB_t {		// Polarized Target Data Bank 
	int	PBeam;          //  beam polarization (% x 100)
	int	PB_DATE;        //  date of measurement
	int	TTYPE;          //  Target Type 
	int	B_TARG;         //  Target holding field (Tesla x100)
	int	PTarg;          //  Target polarization
	int	PT_time;        //  Time measured (according to PC)
	int	HE_Level;       //  Liquid Helium Level(% x 100)
	int	EIOF;           //  EIO microwave tube frequency (MHz)
	int	TTempH;         //  Target temp measured by He Cell (K x100)
	int	TTempC;         //  Target temp measured by Cernox (K x100)
	int	AnealT;         //  Date of last target anealing
} PTDB_t;
typedef struct RCST_t {		//  Roc status bank 
	int	down_count;     //  download count
	int	prestart_count; //  prestart count
	int	go_count;       //  go count
	int	trig_count;     //  total trigger count
	int	event_count;    //  phys event count
	int	sync_count;     //  force-sync event count
	int	run_trig_count; //  trigger count current run
	int	run_event_count; //  phys event count current run
	int	run_sync_count; //  force-sync event count current run
	int	pause_count;    //  coda pause count
	int	end_count;      //  end count
	int	illegal_count;  //  illegal count
	int	run_illegal_count; //  illegal count current run
	int	phys_sync_count; //  physics-sync event count
	int	run_phys_sync_count; //  physics-sync event count current run
} RCST_t;
typedef struct REF_t {		//  Reference signals for pipeline TDCs 
	unsigned short	ID;             //  roc id (high byte) & slot number (low byte)
	unsigned short	TDC;            //  tdc information (channels)
} REF_t;
typedef struct RF_t {		//  RF result bank  
	float	RF;             //  Best RF value in ns
	float	RF1;            //  RF1 in ns
	float	RF2;            //  RF2 in ns
} RF_t;
typedef struct RFT_t {		//  RF signal TDC (pipeline TDCs) 
	unsigned short	ID;             //  id
	unsigned short	TDC;            //  tdc information (channels)
} RFT_t;
typedef struct RGLK_t {		//  single region hits bank 
	int	IREGION;        //  region
	float	X;              //  X pos. of hit in CLAS (cm)
	float	Y;              //  Y pos. of hit in CLAS (cm)
	float	Z;              //  Z pos. of hit in CLAS (cm)
	float	THETA0;         //  polar angle of the link position (deg)
	float	PHI0;           //  azim. angle of the link position (deg)
	float	RTHETA;         //  polar angle of the link direction (deg)
	float	RPHI;           //  azim. angle of the link direction (deg)
	float	CHI2;           //  fit chi2        
	int	STATUS;         //  MINUIT fit status (from 0=bad to 3=ok)
} RGLK_t;
typedef struct RNLG_t {		//  Database archive bank  
	char	sql[4];         
} RNLG_t;
typedef struct RNPE_t {		//  Calibration run index 
	int	run;            //  calibration run number
	char	bank[4];        
} RNPE_t;
typedef struct RTSL_t {		//  Photon Flux Scaler bank  
	int	RAWT1;          
	int	RAWT2;          
	int	RAWT3;          
	int	RAWT4;          
	int	RAWT5;          
	int	RAWT6;          
	int	RAWT7;          
	int	RAWT8;          
	int	RAWT9;          
	int	RAWT10;         
	int	RAWT11;         
	int	RAWT12;         
	int	RAWT13;         
	int	RAWT14;         
	int	RAWT15;         
	int	RAWT16;         
	int	RAWT17;         
	int	RAWT18;         
	int	RAWT19;         
	int	RAWT20;         
	int	RAWT21;         
	int	RAWT22;         
	int	RAWT23;         
	int	RAWT24;         
	int	RAWT25;         
	int	RAWT26;         
	int	RAWT27;         
	int	RAWT28;         
	int	RAWT29;         
	int	RAWT30;         
	int	RAWT31;         
	int	RAWT32;         
	int	RAWT33;         
	int	RAWT34;         
	int	RAWT35;         
	int	RAWT36;         
	int	RAWT37;         
	int	RAWT38;         
	int	RAWT39;         
	int	RAWT40;         
	int	RAWT41;         
	int	RAWT42;         
	int	RAWT43;         
	int	RAWT44;         
	int	RAWT45;         
	int	RAWT46;         
	int	RAWT47;         
	int	RAWT48;         
	int	RAWT49;         
	int	RAWT50;         
	int	RAWT51;         
	int	RAWT52;         
	int	RAWT53;         
	int	RAWT54;         
	int	RAWT55;         
	int	RAWT56;         
	int	RAWT57;         
	int	RAWT58;         
	int	RAWT59;         
	int	RAWT60;         
	int	RAWT61;         
	int	RAWT62;         
	int	RAWT63;         
	int	RAWT64;         
} RTSL_t;
typedef struct RUNC_t {		//  Run Control Map information 
	int	runno;          //  Run number extracted from map
	float	beam_e;         //  beam energy in MeV
	int	b_first;        //  firsttime
	float	q_live;         //  Fcup * live time (1*E-10C)
	int	ql_first;       //  firsttime	
	float	q_all;          //  Fcup (1*E-10C)
	int	qa_first;       //  firsttime	
	float	tor_curr;       //  Torus Current (A)
	int	t_first;        //  firsttime	
	float	mtor_curr;      //  Minitorus Current (A)
	int	mt_first;       //  firsttime	
	float	tag_curr;       //  Tagger Current (A)
	int	tag_first;      //  firsttime	
} RUNC_t;
typedef struct S1ST_t {		//  Scaler bank  
	int	latch1_bit1_count; //  Count trigger bit 1  latched events
	int	latch1_bit2_count; //  Count trigger bit 2  latched events
	int	latch1_bit3_count; //  Count trigger bit 3  latched events
	int	latch1_bit4_count; //  Count trigger bit 4  latched events
	int	latch1_bit5_count; //  Count trigger bit 5  latched events
	int	latch1_bit6_count; //  Count trigger bit 6  latched events
	int	latch1_bit7_count; //  Count trigger bit 7  latched events
	int	latch1_bit8_count; //  Count trigger bit 8  latched events
	int	latch1_bit9_count; //  Count trigger bit 9  latched events
	int	latch1_bit10_count; //  Count trigger bit 10 latched events
	int	latch1_bit11_count; //  Count trigger bit 11 latched events
	int	latch1_bit12_count; //  Count trigger bit 12 latched events
	int	event_count;    //  Latched event count this run
	int	latch1_zero_count; //  Latch1 zero count (illegal)
	int	latch1_empty_count; //  Latch1 empty count (illegal)
	int	latch1_not_empty_count; //  Latch1 not empty on sync event (illegal)
	int	latch1_ok_count; //  Latch1 ok
	int	level2_sector1; //  Level2 sector1 count
	int	level2_sector2; //  Level2 sector2 count
	int	level2_sector3; //  Level2 sector3 count
	int	level2_sector4; //  Level2 sector4 count
	int	level2_sector5; //  Level2 sector5 count
	int	level2_sector6; //  Level2 sector6 count
	int	level2_pass;    //  Level2 pass count
	int	level2_fail;    //  Level2 fail count
	int	latch2_zero_count; //  Latch2 zero count (illegal)
	int	latch2_empty_count; //  Latch2 empty count (illegal)
	int	latch2_not_empty_count; //  Latch2 not empty on sync event (illegal)
	int	latch2_ok_count; //  Latch2 ok
	int	roc_13_count;   //  Roc code 13 count (zero latch)
	int	roc_15_count;   //  Roc code 15 count (illegal)
	int	l1l2_zero_count; //  (latch1==0)&&(latch2==0)
	int	l1zero_13_count; //  (latch1==0)&&(roc_code==13)
	int	l2zero_13_count; //  (latch2==0)&&(roc_code==13)
	int	l1l2zero_13_count; //  (latch1==0)&&(latch2==0)&&(roc_code==13)
} S1ST_t;
typedef struct SC1_t {		//  Scintillation counter hits bank 
	int	ID;             //  the address of the hit detector element
	float	time_l;         //  time for left paddle(ns) 
	float	energy_l;       //  energy in left paddle(MeV) 
	float	time_r;         //  time for right paddle(ns) 
	float	energy_r;       //  energy in right paddle(MeV) 
	float	dtime_l;        //  uncertainty in time for left paddle(ns) 
	float	denergy_l;      //  uncertainty in energy in left paddle(MeV) 
	float	dtime_r;        //  uncertainty in time for right paddle(ns) 
	float	denergy_r;      //  uncertainty in energy in right paddle(MeV) 
} SC1_t;
typedef struct SCC_t {		//  GSIM & SC Scintillator calibration info 
	int	id;             //  paddle id# 
	int	date;           //  UNIX univiersal time of calibraton (32 bits) 
	int	version;        //  sequential version# of calibration 
	int	status;         //  4 byte status word (see details) 
	float	TD0L;           //  (Left) gives 0 time at center for tube (ns) 
	float	TD0Lu;          //  (Left) uncertainty (ns) 
	float	TD0R;           //  (Right) gives 0 time at center for tube (ns) 
	float	TD0Ru;          //  (Right) uncertainty (ns) 
	float	TD1L;           //  (Left) TDC**1 coefficient (ns/ch) 
	float	TD1Lu;          //  (Left) uncertainty (ns/ch) 
	float	TD1R;           //  (Right) TDC**1 coefficient (ns/ch) 
	float	TD1Ru;          //  (Right) uncertainty (ns/ch) 
	float	TD2L;           //  (Left) TDC**2 coefficient 
	float	TD2Lu;          //  (Left) uncertainty (ns/ch) 
	float	TD2R;           //  (Right)TDC**2 coefficient 
	float	TD2Ru;          //  (Right) uncertainty (ns/ch) 
	float	TW0L;           //  (Left) time walk constant (ns) 
	float	TW0Lu;          //  (Left) uncertainty (ns) 
	float	TW0R;           //  (Right) time walk constant (ns) 
	float	TW0Ru;          //  (Right) uncertainty (ns) 
	float	TW1L;           //  (Left) time walk 1st factor 
	float	TW1Lu;          //  (Left) uncertainty 
	float	TW1R;           //  (Right) time walk 1st factor 
	float	TW1Ru;          //  (Right) uncertainty 
	float	TW2L;           //  (Left) time walk 2nd factor 
	float	TW2Lu;          //  (Left) uncertainty 
	float	TW2R;           //  (Right) time walk 2nd factor 
	float	TW2Ru;          //  (Right) uncertainty 
	float	ADPL;           //  (Left) ADC pedestal (ch) 
	float	ADPLu;          //  (Left) uncertainty 
	float	ADPR;           //  (Right) ADC pedestal (ch) 
	float	ADPRu;          //  (Right) uncertainty 
	float	M0L;            //  (Left) nmip adc channel 
	float	M0Lu;           //  (Left) uncertainty (ch) 
	float	M0R;            //  (Right)nmip adc channel 
	float	M0Ru;           //  (Right) uncertainty(ch) 
	float	VEFL;           //  (Left) effective velocity of light (cm/ns) 
	float	VEFLu;          //  (Left) uncertainty (cm/ns) 
	float	VEFR;           //  (Right) effective velocity of light (cm/ns) 
	float	VEFRu;          //  (Right) uncertainty (cm/ns) 
	float	ATNL;           //  (Left) attenuation length (cm) 
	float	ATNLu;          //  (Left) uncertainty (cm) 
	float	ATNR;           //  (Right) attenuation length (cm) 
	float	ATNRu;          //  (Right) uncertainty (cm) 
	float	TSIG0L;         //  (Left) 1st parameter of measured resolution (ns) 
	float	TSIG0R;         //  (Right) 1st parameter of measured resolution (ns) 
	float	TSIG1L;         //  (Left) 2nd parameter of measured resolution 
	float	TSIG1R;         //  (Right) 2nd parameter of measured resolution 
} SCC_t;
typedef struct SC_t {		//  Scintillation counter event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDCL;           //  tdc information (channels)
	unsigned short	ADCL;           //  adc information (channels)
	unsigned short	TDCR;           //  tdc information (channels)
	unsigned short	ADCR;           //  adc information (channels)
} SC_t;
typedef struct SCDI_t {		//  SC channel discriminator thresholds 
	int	crate;          //  CAMAC crate number
	int	slot;           //  slot
	int	mask;           //  mask
	int	threshold;      //  actual threshold value (mV)
	int	width;          //  actual width value
} SCDI_t;
typedef struct SCGD_t {		// SC Geometry of paddles for DC (in SS) - M.K. 
	int	ID;             //  paddle # inside the plane (1-23(max))
	float	NORM_Z;         //  Z(X') coordinate of unit normal to plane
	float	NORM_X;         //  Z(X') coordinate of unit normal to plane
	float	NORM_D;         //  distance to plane along unit vector
	float	ALON_Z;         //  half width in X direction (along id's)
	float	ALON_X;         //  half length of SC in Y direction
	float	BEG_PD;         //  begin point of the paddle along id
	float	END_PD;         //  end point of the paddle along id
	float	PDL_SH;         //  shift of the paddle in respect to the Mid Plane
	float	R_BEAM;         //  distance from the beam to the counte center in M.P.
} SCGD_t;
typedef struct SCG_t {		//  GSIM & SC Scintillator geometry info 
	int	id;             //   paddle id#
	int	panel;          //  panel number - to which plane (panel) the paddle belongs
	float	Xccw;           //  X center of CCW end of paddle in CLAS system
	float	Yccw;           //  Y center of CCW end 
	float	Zccw;           //  - Z center of CCW end
	float	Xcw;            //  X center of CW end 
	float	Ycw;            //  Y center of CW end
	float	Zcw;            // - Z center of CW end 
	float	WIDTH;          //  width (cm) (~|| CLAS theta) 
	float	THICK;          //  thickness (cm) (~|| CLAS r) 
	float	DELTA;          //  the difference in length of the higher-theta slab minus the lower-theta slab divided by 4 (cm) 
} SCG_t;
typedef struct SCH_t {		//  GSIM Scintillator hit info 
	float	x;              //  x of hit
	float	y;              //  y of hit
	float	z;              //  z of hit
	float	cx;             //  track x dir cosine
	float	cy;             //  track y dir cosine
	float	cz;             //  track z dir cosine
	float	pmom;           //  track momentum
	int	track;          //  track number
	int	id;             //  track PDG id
} SCH_t;
typedef struct SCMD_t {		// Sc reconstruction bank for DC M.K. 
	int	ID;             //  paddle # 
	int	STATUS;         //  general status information 
	float	X_NORM;         //  X of unit vector normal to paddle
	float	Y_NORM;         //  Y of unit vector normal to paddle
	float	D_NORM;         //  normal distance to paddle plane
	float	X_PANEL;        //  X of unit vector along the panel
	float	Y_PANEL;        //  Y of unit vector along the panel
	float	B_PAD;          //  Begin of paddle along the panel
	float	E_PAD;          //  End of paddle along the panel
	float	TOF;            //  time of flight (nS)
	float	EDP;            //  energy deposited (GeV) 
	float	P_TIME;         //  position i.r.t. mid.plane from time
	float	P_ATLN;         //  position i.r.t. mid.plane from atten.
} SCMD_t;
typedef struct SCMT_t {		//  Mean SC pretrigger thresholds 
	int	mean_thr;       //  mean threshold (mV)
} SCMT_t;
typedef struct SCMW_t {		//  Mean SC pretrigger threshold width 
	int	mean_width;     //  mean threshold width
} SCMW_t;
typedef struct SCPB_t {		//  EC hits involved in the event 
	int	ScPdHt;         //  10000*sector+100*SC_PD_ID+Hit_ID in SCR 
	float	Edep;           //  Deposited energy (dE/dX)
	float	Time;           //  measured time 
	float	Path;           //  Path lenght from target
	float	Chi2SC;         //  Quality measure of geometrical matching
	int	Status;         //  Status word (not defined yet)
} SCPB_t;
typedef struct SCP_t {		// Sc Scintillator panel geometry bank  
	float	N1x;            //  x component of outward normal to panel 1
	float	N1y;            //  y component of outward normal to panel 1
	float	N1z;            //  z component of outward normal to panel 1
	float	R1n;            //  min.dist.inner(cm)from origin to panel 1
	float	N2x;            //  x component of outward normal to panel 2
	float	N2y;            //  y component of outward normal to panel 2
	float	N2z;            //  z component of outward normal to panel 2
	float	R2n;            //  min.dist.inner(cm)from origin to panel 2
	float	N3x;            //  x component of outward normal to panel 3
	float	N3y;            //  y component of outward normal to panel 3
	float	N3z;            //  z component of outward normal to panel 3
	float	R3n;            //  min.dist.inner(cm)from origin to panel 3
	float	N4x;            //  x component of outward normal to panel 4
	float	N4y;            //  y component of outward normal to panel 4
	float	N4z;            //  z component of outward normal to panel 4
	float	R4n;            //  min.dist.inner(cm)from origin to panel 4
} SCP_t;
typedef struct SCPE_t {		// Translated SC pedestal bank, measured during the pedestal run.  
	int	ID;             //  the address of the hit detector element
	int	mean_left;      //  left adc pedestal mean value (channel)
	float	sigma_left;     //  sigma of the pedestal distribution for left adc (channel)
	int	mean_right;     //  right adc pedestal mean value (channel)
	float	sigma_right;    //  sigma of the pedestal distribution for right adc (channel
} SCPE_t;
typedef struct SCPS_t {		//  Sc plane geometry in SS 
	int	ID;             //  Plane #(1-6):1-23,24-34,35-39,40-42,43-46,47-48
	int	FIRSTPNUM;      //  The first paddle number in the plane
	int	LASTPNUM;       //  The last paddle number in the plane
	float	DISTANCE;       //  distance to the center of the plane
	float	BEGALONSHFT;    //  X-shift of Begin Plane poin along id of SC
	float	ENDALONSHFT;    //  X-shift of End of Plane poin along id of SC
	float	LENGTHSHFT;     //  Y-shift of Plane SYS along the SC length
	float	ALONGIDX;       //  XP.x direction for width in Sector System
	float	ALONGIDY;       //  XP.y direction for width in Sector System
	float	ALONGIDZ;       //  XP.z direction for width in Sector System
	float	SLENGTHX;       //  YP.x direction for length in Sector System
	float	SLENGTHY;       //  YP.y direction for length in Sector System
	float	SLENGTHZ;       //  YP.z direction for length in Sector System
	float	UNORMALX;       //  ZP.x direction for thickness in Sector System
	float	UNORMALY;       //  ZP.y direction for thickness in Sector System
	float	UNORMALZ;       //  ZP.z direction for thickness in Sector System
} SCPS_t;
typedef struct SCRC_t {		// Sc Scintillator reconstruction hit bank 
	int	id;             //  cluster id   
	float	energy;         //  cluster Energy (MeV) 
	float	denergy;        //  error in cluster energy (ns)
	float	time;           //  cluster (energy-weighted) time(ns) 
	float	dtime;          //  error in cluster time (ns)
	float	x;              //  x position in sector coordinate system 
	float	y;              //  y position in sector coordinate system 
	float	z;              //  z position in sector coordinate system 
	float	dx;             //  x error in sector coordinate system 
	float	dy;             //  y error in sector coordinate system 
	float	dz;             //  z error in sector coordinate system 
	int	status;         //  status word defined in sc.h
} SCRC_t;
typedef struct SCR_t {		// Sc Scintillator reconstruction hit bank 
	int	id;             //  paddle id  
	float	energy;         //  Energy (MeV) 
	float	time;           //  time(ns)
	float	x;              //  x position in sector coodinate system 
	float	y;              //  y position in sector coodinate system 
	float	z;              //  z position in sector coodinate system 
	float	dx;             //  x error in sector coodinate system 
	float	dy;             //  y error in sector coodinate system 
	float	dz;             //  z error in sector coodinate system 
	int	status;         //  status word defined in sc.h
	float	denergy;        //  uncertainty in Energy (MeV) 
	float	dtime;          //  uncertainty in time (ns) 
} SCR_t;
typedef struct SCS_t {		//  Time of Flight scaler bank 
	int	SCS1;           //  scaler sector 1
	int	SCS2;           //  scaler sector 1
	int	SCS3;           //  scaler sector 1
	int	SCS4;           //  scaler sector 1
	int	SCS5;           //  scaler sector 1
	int	SCS6;           //  scaler sector 1
	int	SCS7;           //  scaler sector 1
	int	SCS8;           //  scaler sector 1
	int	SCS9;           //  scaler sector 1
	int	SCS10;          //  scaler sector 1
	int	SCS11;          //  scaler sector 1
	int	SCS12;          //  scaler sector 1
	int	SCS13;          //  scaler sector 1
	int	SCS14;          //  scaler sector 1
	int	SCS15;          //  scaler sector 1
	int	SCS16;          //  scaler sector 1
	int	SCS17;          //  scaler sector 1
	int	SCS18;          //  scaler sector 1
	int	SCS19;          //  scaler sector 1
	int	SCS20;          //  scaler sector 1
	int	SCS21;          //  scaler sector 1
	int	SCS22;          //  scaler sector 1
	int	SCS23;          //  scaler sector 1
	int	SCS24;          //  scaler sector 1
	int	SCS25;          //  scaler sector 1
	int	SCS26;          //  scaler sector 1
	int	SCS27;          //  scaler sector 1
	int	SCS28;          //  scaler sector 1
	int	SCS29;          //  scaler sector 1
	int	SCS30;          //  scaler sector 1
	int	SCS31;          //  scaler sector 1
	int	SCS32;          //  scaler sector 1
	int	SCS33;          //  scaler sector 2
	int	SCS34;          //  scaler sector 2
	int	SCS35;          //  scaler sector 2
	int	SCS36;          //  scaler sector 2
	int	SCS37;          //  scaler sector 2
	int	SCS38;          //  scaler sector 2
	int	SCS39;          //  scaler sector 2
	int	SCS40;          //  scaler sector 2
	int	SCS41;          //  scaler sector 2
	int	SCS42;          //  scaler sector 2
	int	SCS43;          //  scaler sector 2
	int	SCS44;          //  scaler sector 2
	int	SCS45;          //  scaler sector 2
	int	SCS46;          //  scaler sector 2
	int	SCS47;          //  scaler sector 2
	int	SCS48;          //  scaler sector 2
	int	SCS49;          //  scaler sector 2
	int	SCS50;          //  scaler sector 2
	int	SCS51;          //  scaler sector 2
	int	SCS52;          //  scaler sector 2
	int	SCS53;          //  scaler sector 2
	int	SCS54;          //  scaler sector 2
	int	SCS55;          //  scaler sector 2
	int	SCS56;          //  scaler sector 2
	int	SCS57;          //  scaler sector 2
	int	SCS58;          //  scaler sector 2
	int	SCS59;          //  scaler sector 2
	int	SCS60;          //  scaler sector 2
	int	SCS61;          //  scaler sector 2
	int	SCS62;          //  scaler sector 2
	int	SCS63;          //  scaler sector 2
	int	SCS64;          //  scaler sector 2
	int	SCS65;          //  scaler sector 3
	int	SCS66;          //  scaler sector 3
	int	SCS67;          //  scaler sector 3
	int	SCS68;          //  scaler sector 3
	int	SCS69;          //  scaler sector 3
	int	SCS70;          //  scaler sector 3
	int	SCS71;          //  scaler sector 3
	int	SCS72;          //  scaler sector 3
	int	SCS73;          //  scaler sector 3
	int	SCS74;          //  scaler sector 3
	int	SCS75;          //  scaler sector 3
	int	SCS76;          //  scaler sector 3
	int	SCS77;          //  scaler sector 3
	int	SCS78;          //  scaler sector 3
	int	SCS79;          //  scaler sector 3
	int	SCS80;          //  scaler sector 3
	int	SCS81;          //  scaler sector 3
	int	SCS82;          //  scaler sector 3
	int	SCS83;          //  scaler sector 3
	int	SCS84;          //  scaler sector 3
	int	SCS85;          //  scaler sector 3
	int	SCS86;          //  scaler sector 3
	int	SCS87;          //  scaler sector 3
	int	SCS88;          //  scaler sector 3
	int	SCS89;          //  scaler sector 3
	int	SCS90;          //  scaler sector 3
	int	SCS91;          //  scaler sector 3
	int	SCS92;          //  scaler sector 3
	int	SCS93;          //  scaler sector 3
	int	SCS94;          //  scaler sector 3
	int	SCS95;          //  scaler sector 3
	int	SCS96;          //  scaler sector 3
	int	SCS97;          //  scaler sector 4
	int	SCS98;          //  scaler sector 4
	int	SCS99;          //  scaler sector 4
	int	SCS100;         //  scaler sector 4
	int	SCS101;         //  scaler sector 4
	int	SCS102;         //  scaler sector 4
	int	SCS103;         //  scaler sector 4
	int	SCS104;         //  scaler sector 4
	int	SCS105;         //  scaler sector 4
	int	SCS106;         //  scaler sector 4
	int	SCS107;         //  scaler sector 4
	int	SCS108;         //  scaler sector 4
	int	SCS109;         //  scaler sector 4
	int	SCS110;         //  scaler sector 4
	int	SCS111;         //  scaler sector 4
	int	SCS112;         //  scaler sector 4
	int	SCS113;         //  scaler sector 4
	int	SCS114;         //  scaler sector 4
	int	SCS115;         //  scaler sector 4
	int	SCS116;         //  scaler sector 4
	int	SCS117;         //  scaler sector 4
	int	SCS118;         //  scaler sector 4
	int	SCS119;         //  scaler sector 4
	int	SCS120;         //  scaler sector 4
	int	SCS121;         //  scaler sector 4
	int	SCS122;         //  scaler sector 4
	int	SCS123;         //  scaler sector 4
	int	SCS124;         //  scaler sector 4
	int	SCS125;         //  scaler sector 4
	int	SCS126;         //  scaler sector 4
	int	SCS127;         //  scaler sector 4
	int	SCS128;         //  scaler sector 4
	int	SCS129;         //  scaler sector 5
	int	SCS130;         //  scaler sector 5
	int	SCS131;         //  scaler sector 5
	int	SCS132;         //  scaler sector 5
	int	SCS133;         //  scaler sector 5
	int	SCS134;         //  scaler sector 5
	int	SCS135;         //  scaler sector 5
	int	SCS136;         //  scaler sector 5
	int	SCS137;         //  scaler sector 5
	int	SCS138;         //  scaler sector 5
	int	SCS139;         //  scaler sector 5
	int	SCS140;         //  scaler sector 5
	int	SCS141;         //  scaler sector 5
	int	SCS142;         //  scaler sector 5
	int	SCS143;         //  scaler sector 5
	int	SCS144;         //  scaler sector 5
	int	SCS145;         //  scaler sector 5
	int	SCS146;         //  scaler sector 5
	int	SCS147;         //  scaler sector 5
	int	SCS148;         //  scaler sector 5
	int	SCS149;         //  scaler sector 5
	int	SCS150;         //  scaler sector 5
	int	SCS151;         //  scaler sector 5
	int	SCS152;         //  scaler sector 5
	int	SCS153;         //  scaler sector 5
	int	SCS154;         //  scaler sector 5
	int	SCS155;         //  scaler sector 5
	int	SCS156;         //  scaler sector 5
	int	SCS157;         //  scaler sector 5
	int	SCS158;         //  scaler sector 5
	int	SCS159;         //  scaler sector 5
	int	SCS160;         //  scaler sector 5
	int	SCS161;         //  scaler sector 6
	int	SCS162;         //  scaler sector 6
	int	SCS163;         //  scaler sector 6
	int	SCS164;         //  scaler sector 6
	int	SCS165;         //  scaler sector 6
	int	SCS166;         //  scaler sector 6
	int	SCS167;         //  scaler sector 6
	int	SCS168;         //  scaler sector 6
	int	SCS169;         //  scaler sector 6
	int	SCS170;         //  scaler sector 6
	int	SCS171;         //  scaler sector 6
	int	SCS172;         //  scaler sector 6
	int	SCS173;         //  scaler sector 6
	int	SCS174;         //  scaler sector 6
	int	SCS175;         //  scaler sector 6
	int	SCS176;         //  scaler sector 6
	int	SCS177;         //  scaler sector 6
	int	SCS178;         //  scaler sector 6
	int	SCS179;         //  scaler sector 6
	int	SCS180;         //  scaler sector 6
	int	SCS181;         //  scaler sector 6
	int	SCS182;         //  scaler sector 6
	int	SCS183;         //  scaler sector 6
	int	SCS184;         //  scaler sector 6
	int	SCS185;         //  scaler sector 6
	int	SCS186;         //  scaler sector 6
	int	SCS187;         //  scaler sector 6
	int	SCS188;         //  scaler sector 6
	int	SCS189;         //  scaler sector 6
	int	SCS190;         //  scaler sector 6
	int	SCS191;         //  scaler sector 6
	int	SCS192;         //  scaler sector 6
} SCS_t;
typedef struct SCT_t {		//  Scintillation counter event bank (TDC) 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
} SCT_t;
typedef struct SGMP_t {		//  
	int	is;             //  superlayer number of the segment
	int	sgm1;           //  data structure see: include/bosddl.h, wire_t
	int	sgm2;           //  wire_t
	int	sgm3;           //  wire_t
	int	sgm4;           //  wire_t
	int	sgm5;           //  wire_t
	int	sgm6;           //  wire_t
} SGMP_t;
typedef struct SPAR_t {		//  Untranslated pedestal bank, read out from the hardware. 
	int	slot;           //  ADC slot
	int	channel;        //  ADC channel
	int	spar;           //  sparsification threshold (channel)
	int	pedmean;        //  pedestal mean value measured, using internal 500ns gate (channel)
} SPAR_t;
typedef struct SPIN_t {		//  spin vector of the particle 
	int	xspin;          //  x component of spin
	int	yspin;          //  y component of spin
	int	zspin;          //  z component of spin
} SPIN_t;
typedef struct ST1_t {		//  Start counter intermediate results bank 
	int	ID;             //  Pair id
	int	status;         //  status word
	float	TIME_1;         //  Time (ns) for side 1
	float	ADC_1;          //  Pedestal-subtracted adc for side 1
	float	TIME_2;         //  Time (ns) for side 2
	float	ADC_2;          //  Pedestal-subtracted ADC  for side 2
} ST1_t;
typedef struct ST_t {		//  Start counter event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels)
	unsigned short	ADC;            //  adc information (channels)
} ST_t;
typedef struct STG_t {		//  Start Counter Geometry bank 
	int	id;             //  sector id
	float	leg_x_max;      //  distance from beam axis to inner face of leg 
	float	leg_x_min;      
	float	leg_y_max;      
	float	leg_y_min;      
	float	leg_z_max;      
	float	let_z_min;      
	float	nose_x_max;     //   
	float	nose_x_min;     
	float	nose_y_max;     
	float	nose_y_min;     
	float	nose_z_max;     
	float	nose_z_min;     
	float	noseAngle;      
	float	lleg;           
	float	lnose;          
} STG_t;
typedef struct STH_t {		//  GSIM ST hit info 
	float	x;              //  x of hit
	float	y;              //  y of hit
	float	z;              //  z of hit
	float	cx;             //  track x dir cosine
	float	cy;             //  track y dir cosine
	float	cz;             //  track z dir cosine
	float	pmom;           //  track momentum
	int	track;          //  track number
	int	id;             //  track PDG id
	float	tof;            //  flight time
} STH_t;
typedef struct STN0_t {		//  NEW Start counter TDC bank 
	unsigned short	ID;             //  detector element
	unsigned short	TDC;            //  tdc information (channels)
} STN0_t;
typedef struct STN1_t {		//  NEW Start counter ADC bank 
	unsigned short	ID;             //  detector element 
	unsigned short	ADC;            //  adc information (channels)
} STN1_t;
typedef struct STPB_t {		//  ST hits involved in the event 
	int	SThid;          //  100*sector+Hit_ID  
	float	Time;           //  Flight time relative to the evnt start time
	float	Path;           //  Path lenght from target
	int	charge;         //  charge (get from tracking)
	int	Status;         //  Status word (not defined yet)
} STPB_t;
typedef struct STPE_t {		//  Translated start counter pedestal bank, measured during the pedestal run.  
	int	ID;             //  the address of the hit detector element
	int	mean;           //  adc pedestal mean value (channel)
	float	sigma;          //  sigma of the pedestal distribution (channel)
} STPE_t;
typedef struct STR_t {		//  Start counter results bank 
	int	ID;             //  sector
	int	Trk_no;         //  Pointer to track in HBTR
	float	ST_TIME;        //  flight time from ST (ns)
	float	ST_L;           //  flight path from ST (cm)
	float	st_pos;         //  position within the start counter
	int	status;         //  Status word
} STR_t;
typedef struct STS_t {		//  Start Counter scaler bank 
	int	STS1;           //  trigger rate sector 1-2
	int	STS2;           //  trigger rate sector 3-4
	int	STS3;           //  trigger rate sector 5-6
	int	STS4;           
	int	STS5;           
	int	STS6;           
	int	STS7;           
	int	STS8;           
	int	STS9;           
	int	STS10;          
	int	STS11;          
	int	STS12;          
	int	STS13;          
	int	STS14;          
	int	STS15;          
	int	STS16;          
} STS_t;
typedef struct STSN_t {		//  Start Counter scaler bank 
	int	STSN1;          
	int	STSN2;          
	int	STSN3;          
	int	STSN4;          
	int	STSN5;          
	int	STSN6;          
	int	STSN7;          
	int	STSN8;          
	int	STSN9;          
	int	STSN10;         
	int	STSN11;         
	int	STSN12;         
	int	STSN13;         
	int	STSN14;         
	int	STSN15;         
	int	STSN16;         
	int	STSN17;         
	int	STSN18;         
	int	STSN19;         
	int	STSN20;         
	int	STSN21;         
	int	STSN22;         
	int	STSN23;         
	int	STSN24;         
	int	STOR;           
	int	STMULT;         
	int	STANDMOR;       
	int	MULTANDMOR;     
	int	RES1;           
	int	RES2;           
	int	RES3;           
	int	RES4;           
} STSN_t;
typedef struct SYNC_t {		//  Physics sync event bank 
	unsigned short	ID;             //  slot number
	unsigned short	TDCL;           //  count of missing responses
	unsigned short	ADCL;           //  count of extra buffers
} SYNC_t;
typedef struct TACO_t {		//  down stream output bank 
	int	ID;             //  ID  1= TAC, 2 = USLG
	float	TIME;           //   TDC time 
	float	ELT;            //  energy deposit TAC = left top
	float	ERT;            //  energy deposit TAC = right top
	float	ELB;            //  energy deposit TAC = left bottom
	float	ERB;            //  energy deposit TAC = right bottom
	float	ESUM;           //  energy deposit TAC = sum1
	float	ESUM2;          //  energy deposit TAC = sum2
	float	ESUM3;          //  energy deposit TAC = sum3
	int	TID;            //  T id of the corresponding T
} TACO_t;
typedef struct TAGE_t {		//  
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information (channels) [multihit tdc]
} TAGE_t;
typedef struct TAGI_t {		//  Tagger Intermediary event bank 
	int	IDT;            //  T id
	float	TIMEL;          //  time information (Left counters channels)
	float	TIMER;          //  time information (Right counters channels)
	int	IDE;            //  E id
	float	TIMEE;          //  time information (E counters)
	float	TIMEMEAN;       //  time information (left/right mean value | after alignement)
	float	TRF;            //  time information (mean val - RF | after alignement)
	float	NEXTTIME;       //  time difference with the next T when in coincidence | after alignment)
} TAGI_t;
typedef struct TAGM_t {		//  Tagger result bank based on T-counter Multi-hit TDCs 
	float	energy;         //  Energy of the photon in GeV
	float	t;              //  T-counter time (ns)
	float	e_t;            //  E-counter time (ns)
	int	status;         //  Status (not yet used) 
	int	tid;            //  T channel Id
	int	eid;            //  E channel Id
} TAGM_t;
typedef struct TAGR_t {		//  Tagger result bank (AL-LYM-FR 9/29/1997 --- FYDW) 
	float	ERG;            //  Energy of the photon in GeV
	float	TTAG;           //  Time of the photon has reconstructed in the Tagger
	float	TPHO;           //  Time of the photon after RF correction
	int	STAT;           //  Status ( 7 or 15 are Good) other values have problems (see tag_process_TAGR.F) 
	int	T_id;           //  T counter Id
	int	E_id;           //  E counter Id
} TAGR_t;
typedef struct TAGT_t {		//  Tagger T countert event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDCL;           //  tdc information (Left counters channels)
	unsigned short	TDCR;           //  tdc information (Right counters channels)
} TAGT_t;
typedef struct TATL_t {		//  Tagger left T counter event bank 
	int	ID;             //  the address of the hit detector element
	int	TDC;            //  tdc information
} TATL_t;
typedef struct TATR_t {		//  Tagger right T counter event bank 
	int	ID;             //  the address of the hit detector element
	int	TDC;            //  tdc information
} TATR_t;
typedef struct TBER_t {		//  Time Based Tracking ERror bank 
	float	q_over_p;       //  q/p
	float	lambda;         //  dip angle (pi/2 - theta)
	float	phi;            //  phi
	float	d0;             //  min.distance from (x=0,y=0,z=?)  [cm]
	float	z0;             //  z position of starting point     [cm] 
	float	c11;            //  element C{1,1}
	float	c12;            //  element C{1,2}
	float	c13;            //  element C{1,3}
	float	c14;            //  element C{1,4}
	float	c15;            //  element C{1,5}
	float	c22;            //  element C{2,2}
	float	c23;            //  element C{2,3}
	float	c24;            //  element C{2,4}
	float	c25;            //  element C{2,5}
	float	c33;            //  element C{3,3}
	float	c34;            //  element C{3,4}
	float	c35;            //  element C{3,5}
	float	c44;            //  element C{4,4}
	float	c45;            //  element C{4,5}
	float	c55;            //  element C{5,5}
	float	chi2;           //  Chisquare for this Track
	int	layinfo1;       //  layerhit info
	int	layinfo2;       //  layerhit info&sector&track#in hber
} TBER_t;
typedef struct TBID_t {		//  Time Based particle ID 
	int	track;          //  Track (index to TBTR) or zero if neutral
	int	sec;            //  Sector track is in
	float	beta;           //  Beta of the track in units of c
	float	vtime;          //  vertex time of track
	int	sc_stat;        //  status of hit matching to SC: see sc.h
	int	sc_id;          //  Pointer to SCRC bank
	float	sc_time;        //  SC calibrated time for this track (ns)
	float	sc_qual;        //  quality of match for SC
	float	sc_vtime;       //  time at vertex for SC(ns)
	float	sc_beta;        //  Beta calculated from TOF from SC
	int	cc_stat;        //  status of hit matching to CC: see sc.h
	int	cc_id;          //  pointer to CC01 bank
	float	cc_time;        //  CC calibrated time for this track (ns)
	float	cc_qual;        //  quality of match for CC
	float	cc_vtime;       //  time at vertex for CC(ns)
	float	cc_beta;        //  Beta as calculated by the EC
	int	ec_stat;        //  status of hit matching to ec: see sc.h
	int	ec_id;          //  Pointer to ECHB bank
	float	ec_time;        //  EC calibrated time for this track (ns)
	float	ec_qual;        //  EC quality factor
	float	ec_vtime;       //  time at vertex for EC(ns)
	float	ec_beta;        //  Beta as calculated by the EC
	int	st_stat;        //  status of hit matching to ST
	int	st_id;          //  Pointer to STR bank
	float	st_time;        //  ST calibrated time for this track (ns)
	float	st_qual;        //  ST quality factor
	float	st_vtime;       //  time at vertex for ST(ns)
	float	st_beta;        //  Beta as calculated by the ST
	int	lac_stat;       //  status of hit matching to LAC
	int	lac_id;         //  Pointer to EC1R bank
	float	lac_time;       //  LAC calibrated time for this track (ns)
	float	lac_qual;       //  LAC quality factor
	float	lac_vtime;      //  time at vertex for LAC(ns)
	float	lac_beta;       //  Beta as calculated by the LAC
} TBID_t;
typedef struct TBLA_t {		//  Time Based tracking LAyer bank 
	int	trk_pln;        //  (track_number) *100 + Trk_plane_number
	float	x;              //  z coord [cm]  for track in this plane
	float	y;              //  y coord [cm]  for track in this plane
	float	z;              //  z coord [cm]  for track in this plane
	float	Bx;             //  B-field in x [kG] at coord.{x,y,z}
	float	By;             //  B-field in y [kG] at coord.{x,y,z}
	float	Bz;             //  B-field in z [kG] at coord.{x,y,z}
	float	tlen;           //  track length [cm] from origin to this plane
	int	dc1;            //  Pointer to DC1 bank
	int	stat;           //  Status of the hit
	int	wire;           //  Wire number 
	float	dtime;          //  drift time  [ns]
	float	alpha;          //  track angle (relative to R of SL) [deg]
	float	wlen;           //  Wire length (hit pos. to preamp)  [cm]
	float	sgdoca;         //  sigma DOCA  [cm]
	float	fitdoca;        //  Fitted DOCA [cm]
	float	calcdoca;       //  calculated DOCA (via dtime)  [cm]
} TBLA_t;
typedef struct TBTR_t {		//  Time Based Tracking Result bank 
	float	x;              //  x
	float	y;              //  y 'Vertex' position {x,y,z}
	float	z;              //  z
	float	px;             //  Px
	float	py;             //  Py  momentum at 'vertex' {Px,Py,Pz}
	float	pz;             //  Pz
	float	q;              //  charge   (straight tracks: set q=0.) 
	float	chi2;           //  Chisquare for this Track
	int	itr_sec;        //  Trackno_in_Sector + Sector*100
	int	itr_hbt;        //  Trackno in HBTR for this track
} TBTR_t;
typedef struct TCSB_t {		//  Toroidal Coordinate System Bank 
	float	xpos;           //  x misalignment of TCS in HCS 		 
	float	ypos;           //  y misalignment of TCS in HCS	 
	float	zpos;           //  z misalignment of TCS in HCS	 
	float	sxpos;          //  sx sine of x-axis misorientation of TCS vs HCS		 
	float	sypos;          //  sy sine of y-axis misorientation of TCS vs HCS
	float	szpos;          //  sz sine of z-axis misorientation of TCS vs HCS
} TCSB_t;
typedef struct TCT_t {		//  total absorbtion counter TDC (pipeline TDCs) 
	unsigned short	ID;             //  id
	unsigned short	TDC;            //  tdc information (channels)
} TCT_t;
typedef struct TDPL_t {		//  Tbt Detector PLane bank 
	int	trk_pln;        //  (track_number) *100 + plane_number
	float	x;              //  vector3_t pos;
	float	y;              //   (x, y, z coord. for track in on the plane)
	float	z;              
	float	cx;             //  vector3_t dir;
	float	cy;             //   direction cosines (x,y,z) for track at coord.{x,y,z}
	float	cz;             
	float	tlen;           //  track length [cm] from origin to this plane
} TDPL_t;
typedef struct TESC_t {		//  
	float	value;          //  count rate per E-counter channel
} TESC_t;
typedef struct TGBI_t {		//  Trigger bank 
	unsigned int	latch1;         //  level1 trigger latch word (16 bits)
	unsigned int	helicity_scaler; //  helicity interval count
	unsigned int	interrupt_time; //  interrupt time from microsec clock
	unsigned int	latch2;         //  level2 trigger latch word (16 bits)
	unsigned int	level3;         //  level3 trigger word (32 bits)
} TGBI_t;
typedef struct TGEO_t {		//  Target geometry bank 
	float	x;              //  X position of target
	float	y;              //  Y position of target
	float	z;              //  Z position of target
	float	radius;         //  radius of target
	float	length;         //  lenght of target
	int	material;       //  material of target: 0=empty; 1=Hydrogen; 
} TGEO_t;
typedef struct TGPB_t {		//  Tagger hits matched start counter 
	int	pointer;        //  1000*i+tagstat where "i" is the pointer to the TAGR bank. Negative for trigger photon  
	float	Time;           //  starttime_TAG at interaction point(ns)
	float	Energy;         //  photon energy(GeV)
	float	dt;             //  starttime_ST - starttime_TAG (ns)
} TGPB_t;
typedef struct TGS_t {		//  Photon Flux Scaler bank  
	int	RAWT1;          
	int	RAWT2;          
	int	RAWT3;          
	int	RAWT4;          
	int	RAWT5;          
	int	RAWT6;          
	int	RAWT7;          
	int	RAWT8;          
	int	RAWT9;          
	int	RAWT10;         
	int	RAWT11;         
	int	RAWT12;         
	int	RAWT13;         
	int	RAWT14;         
	int	RAWT15;         
	int	RAWT16;         
	int	RAWT17;         
	int	RAWT18;         
	int	RAWT19;         
	int	RAWT20;         
	int	RAWT21;         
	int	RAWT22;         
	int	RAWT23;         
	int	RAWT24;         
	int	RAWT25;         
	int	RAWT26;         
	int	RAWT27;         
	int	RAWT28;         
	int	RAWT29;         
	int	RAWT30;         
	int	RAWT31;         
	int	RAWT32;         
	int	RAWT33;         
	int	RAWT34;         
	int	RAWT35;         
	int	RAWT36;         
	int	RAWT37;         
	int	RAWT38;         
	int	RAWT39;         
	int	RAWT40;         
	int	RAWT41;         
	int	RAWT42;         
	int	RAWT43;         
	int	RAWT44;         
	int	RAWT45;         
	int	RAWT46;         
	int	RAWT47;         
	int	RAWT48;         
	int	RAWT49;         
	int	RAWT50;         
	int	RAWT51;         
	int	RAWT52;         
	int	RAWT53;         
	int	RAWT54;         
	int	RAWT55;         
	int	RAWT56;         
	int	RAWT57;         
	int	RAWT58;         
	int	RAWT59;         
	int	RAWT60;         
	int	RAWT61;         
	int	notused62;      
	int	notused63;      
	int	notused64;      
	int	BNK1T1;         
	int	BNK1T2;         
	int	BNK1T3;         
	int	BNK1T4;         
	int	BNK1T5;         
	int	BNK1T6;         
	int	BNK1T7;         
	int	BNK1T8;         
	int	BNK1T9;         
	int	BNK1T10;        
	int	BNK1T11;        
	int	BNK1T12;        
	int	BNK1T13;        
	int	BNK1T14;        
	int	BNK1T15;        
	int	BNK1T16;        
	int	BNK1T17;        
	int	BNK1T18;        
	int	BNK1T19;        
	int	BNK1T20;        
	int	BNK1T21;        
	int	BNK1T22;        
	int	BNK1T23;        
	int	BNK1T24;        
	int	BNK1T25;        
	int	BNK1T26;        
	int	BNK1T27;        
	int	BNK1T28;        
	int	BNK1T29;        
	int	BNK1T30;        
	int	BNK1T31;        
	int	BNK1T32;        
	int	BNK1T33;        
	int	BNK1T34;        
	int	BNK1T35;        
	int	BNK1T36;        
	int	BNK1T37;        
	int	BNK1T38;        
	int	BNK1T39;        
	int	BNK1T40;        
	int	BNK1T41;        
	int	BNK1T42;        
	int	BNK1T43;        
	int	BNK1T44;        
	int	BNK1T45;        
	int	BNK1T46;        
	int	BNK1T47;        
	int	BNK1T48;        
	int	BNK1T49;        
	int	BNK1T50;        
	int	BNK1T51;        
	int	BNK1T52;        
	int	BNK1T53;        
	int	BNK1T54;        
	int	BNK1T55;        
	int	BNK1T56;        
	int	BNK1T57;        
	int	BNK1T58;        
	int	BNK1T59;        
	int	BNK1T60;        
	int	BNK1T61;        
	int	notused126;     
	int	notused127;     
	int	notused128;     
	int	BNK2T1;         
	int	BNK2T2;         
	int	BNK2T3;         
	int	BNK2T4;         
	int	BNK2T5;         
	int	BNK2T6;         
	int	BNK2T7;         
	int	BNK2T8;         
	int	BNK2T9;         
	int	BNK2T10;        
	int	BNK2T11;        
	int	BNK2T12;        
	int	BNK2T13;        
	int	BNK2T14;        
	int	BNK2T15;        
	int	BNK2T16;        
	int	BNK2T17;        
	int	BNK2T18;        
	int	BNK2T19;        
	int	BNK2T20;        
	int	BNK2T21;        
	int	BNK2T22;        
	int	BNK2T23;        
	int	BNK2T24;        
	int	BNK2T25;        
	int	BNK2T26;        
	int	BNK2T27;        
	int	BNK2T28;        
	int	BNK2T29;        
	int	BNK2T30;        
	int	BNK2T31;        
	int	BNK2T32;        
	int	BNK2T33;        
	int	BNK2T34;        
	int	BNK2T35;        
	int	BNK2T36;        
	int	BNK2T37;        
	int	BNK2T38;        
	int	BNK2T39;        
	int	BNK2T40;        
	int	BNK2T41;        
	int	BNK2T42;        
	int	BNK2T43;        
	int	BNK2T44;        
	int	BNK2T45;        
	int	BNK2T46;        
	int	BNK2T47;        
	int	BNK2T48;        
	int	BNK2T49;        
	int	BNK2T50;        
	int	BNK2T51;        
	int	BNK2T52;        
	int	BNK2T53;        
	int	BNK2T54;        
	int	BNK2T55;        
	int	BNK2T56;        
	int	BNK2T57;        
	int	BNK2T58;        
	int	BNK2T59;        
	int	BNK2T60;        
	int	BNK2T61;        
	int	notused190;     
	int	notused191;     
	int	notused192;     
} TGS_t;
typedef struct TGTL_t {		//  Tagger T countert event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information
} TGTL_t;
typedef struct TGTR_t {		//  Tagger T countert event bank 
	unsigned short	ID;             //  the address of the hit detector element
	unsigned short	TDC;            //  tdc information
} TGTR_t;
typedef struct TLV1_t {		//  Level1 trigger bank 
	int	ID;             //  id = Detector * 256 + channel
	float	time;           //  time(ns) 
} TLV1_t;
typedef struct TRCF_t {		//  Trigger configuration 
	char	trigcfg[4];     
} TRCF_t;
typedef struct TRGS_t {		//  Scaler bank  
	int	CLOCK_UG;       //  Clock ungated
	int	FCUP_UG;        //  FCUP ungated
	int	Microsec;       //  Microsecond clock (will overflow during normal run)
	int	random_ug;      //  Random pulser ungated
	int	MOR_ST;         //  MOR.ST
	int	MOR_PC;         //  MOR.PC
	int	MOR_PS;         //  MOR.PS
	int	MOR_TAC;        //  MOR.TAC
	int	MOR;            //  Master OR
	int	PC;             //  Pair Counter
	int	PS;             //  Pair Spectrometer
	int	TAC;            //  Total Absorption Counter
	int	ST;             //  ST
	int	hel_sync;       //  Helicity sync signal
	int	clock_ug_2;     //  Duplicate of channel 1
	int	fcup_ug_2;      //  Duplicate of channel 2
	int	CLOCK_G1;       //  Clock with run gatei
	int	FCUP_G1;        //  FCUP with Run Gate
	int	NotUsed19;      //  Currently not used
	int	random_g1;      //  Random pulser run gated
	int	MOR_ST_rg;      //  MOR.ST  with run gate
	int	MOR_PC_rg;      //  MOR.PC with run gate
	int	MOR_PS_rg;      //  MOR.PS with run gate
	int	MOR_TAC_rg;     //  MOR.TAC with run gate
	int	MOR_rg;         //  MASTER_OR with run gate
	int	PC_rg;          //  PC with run gate
	int	PS_rg;          //  PS with run gate
	int	TAC_rg;         //  TAC with run gate
	int	ST_rg;          //  ST  with run gate
	int	random_g1_2;    //  duplicate of channel 20
	int	clock_g1_2;     //  duplicate of channel 17
	int	fcup_g1_2;      //  duplicate of channel 18
	int	CLOCK_G2;       //  CLOCK live gated
	int	FCUP_G2;        //  FCUP live gated
	int	trig_or_g2;     //  Trigger OR of 1-12 live gated
	int	random_g2;      //  Random pulser live gated
	int	NotUsed37;      //  Currently not used
	int	NotUsed38;      //  Currently not used
	int	NotUsed39;      //  Currently not used
	int	NotUsed40;      //  Currently not used
	int	MOR_lg;         //  MASTER_OR live gated
	int	NotUsed42;      //  Currently not used
	int	NotUsed43;      //  Currently not used
	int	NotUsed44;      //  Currently not used
	int	NotUsed45;      //  Currently not used
	int	random_g2_2;    //  duplicate of channel 36
	int	clock_g2_2;     //  duplicate of channel 33
	int	fcup_g2_2;      //  duplicate of channel 34
	int	trig1_ug;       //  Trigger 1 ungated, prescaled
	int	trig2_ug;       //  Trigger 2 ungated, prescaled
	int	trig3_ug;       //  Trigger 3 ungated, prescaled
	int	trig4_ug;       //  Trigger 4 ungated, prescaled
	int	trig5_ug;       //  Trigger 5 ungated, prescaled
	int	trig6_ug;       //  Trigger 6 ungated, prescaled
	int	trig7_ug;       //  Trigger 7 ungated, prescaled
	int	trig8_ug;       //  Trigger 8 ungated, prescaled
	int	trig9_ug;       //  Trigger 9 ungated, prescaled
	int	trig10_ug;      //  Trigger 10 ungated, prescaled
	int	trig11_ug;      //  Trigger 11 ungated, prescaled
	int	trig12_ug;      //  Trigger 12 ungated, prescaled
	int	trig_or_ug;     //  Trigger OR of 1-12 ungated,
	int	l1accept;       //  Level 1 accept
	int	notused63;      //  Currently not used
	int	notused64;      //  Currently not used
	int	l2fail;         //  Level2 fail
	int	l2pass;         //  Level2 pass
	int	l2start;        //  Level2 start
	int	l2clear;        //  Level2 clear
	int	l2accept;       //  Level2 accept
	int	l3accept;       //  Level3 accept
	int	notused71;      //  Currently not used
	int	notused72;      //  Currently not used
	int	l2sec1_g;       //  Level2 sec1 gated
	int	l2sec2_g;       //  Level2 sec2 gated
	int	l2sec3_g;       //  Level2 sec3 gated
	int	l2sec4_g;       //  Level2 sec4 gated
	int	l2sec5_g;       //  Level2 sec5 gated
	int	l2sec6_g;       //  Level2 sec6 gated
	int	l2_or_g;        //  OR level2 gated
	int	l2_ok_g;        //  Level 2 OK gated
	int	trig1_lg;       //  Trigger 1 live gated
	int	trig2_lg;       //  Trigger 2 live gated
	int	trig3_lg;       //  Trigger 3 live gated
	int	trig4_lg;       //  Trigger 4 live gated
	int	trig5_lg;       //  Trigger 5 live gated
	int	trig6_lg;       //  Trigger 6 live gated
	int	trig7_lg;       //  Trigger 7 live gated
	int	trig8_lg;       //  Trigger 8 live gated
	int	trig9_lg;       //  Trigger 9 live gated
	int	trig10_lg;      //  Trigger 10 live gated
	int	trig11_lg;      //  Trigger 11 live gated
	int	trig12_lg;      //  Trigger 12 live gated
	int	notused93;      //  not used
	int	notused94;      //  not used
	int	ignore95;       //  ignore
	int	ignore96;       //  ignore
} TRGS_t;
typedef struct TRKS_t {		//  TRacK Status Bank  
	int	trk_lev;        //  Tracknr + analysis_level*100
	float	beta;           //  Beta of the track in unit of c
	float	st_time;        //  Start time of the event (track: photon trigg)
	float	cc_time;        //  CC calibrated time for this track (ns)
	float	sc_time;        //  SC calibrated time for this track (ns)
	float	ec_time;        //  EC calibrated time for this track (ns)
	int	st_id;          //  ST hit id (ptr. to STx: photon trigg)
	int	cc_id;          //  CC hit id (ptr. to CCRC)
	int	sc_id;          //  SC hit id (ptr. to SCRW)
	int	ec_id;          //  EC Cluster-id (ptr. to ECHB)
} TRKS_t;
typedef struct TRL1_t {		//  Tracking Results for Layer 1 
	int	fit_flags;      //  flags used in track fitting
	int	hits_hbt;       //  hits used in HBT fit
	int	hits_tbt;       //  hits used in TBT fit
	float	x;              //  x
	float	y;              //  y   position in first DC layer 
	float	z;              //  z
	float	cx;             //  cx
	float	cy;             //  cy  dir.cosine at first DC layer
	float	cz;             //  cz
	float	tlen;           //  track length to this layer (starting
} TRL1_t;
typedef struct TRPB_t {		//  Tracking Results for Layer 1 
	int	hits_tbt;       //  hits used in TBT fit
	float	x;              //  x
	float	y;              //  y   position in first DC layer 
	float	z;              //  z
	float	cx;             //  cx
	float	cy;             //  cy  dir.cosine at first DC layer
	float	cz;             //  cz
	float	tlen;           //  track length to this layer (starting
} TRPB_t;
typedef struct TSPR_t {		//  Trigger supervisor program memory 
	unsigned int	mem;            //  memory value
} TSPR_t;
typedef struct TSRG_t {		//  Trigger supervisor registers 
	unsigned int	reg;            //  register value
} TSRG_t;
typedef struct UNUS_t {		//  Unused hits/tracks bank 
	int	NDCUN;          //  Number of unused DC tracks
	int	IDCUN;          //  Trek ID  
	int	NSCUN;          //  Number of unused SC hits
	int	ISCUN;          //  SC hit ID
	int	NCCUN;          //  Number of unused CC hits
	int	ICCUN;          //  CC hit ID
	int	NECUN;          //  Number of unused EC hits
	int	IECUN;          //  EC hit ID
} UNUS_t;
typedef struct VERT_t {		//  vertex Result bank 
	int	vertex;         //  vertex id
	int	trk1;           //  track #1
	int	trk2;           //  track #2
	float	x;              //  x vector3_t vert{x,y,z}
	float	y;              //  y 
	float	z;              //  z
	float	sepd;           //  seperation distance
} VERT_t;
typedef struct RC26_t {		//  Trigger bank 
	unsigned int	trigger_pattern; //  level1 trigger pattern from trigger supervisor
	unsigned int	gated_scaler;   //  gated scaler
	unsigned int	ungated_scaler; //  ungated scaler
} RC26_t;
typedef struct STRE_t {		//  extended Start counter results bank 
	int	ID;             //  sector
	int	Trk_no;         //  Pointer to track in HBTR
	float	ST_TIME;        //  flight time from ST (ns)
	float	ST_L;           //  flight path from ST (cm)
	float	st_pos;         //  position within the start counter
	int	status;         //  Status word
	float	st_edep;        //  energy deposit in ST
	float	st_corr;        //  energy correction factor for angle
} STRE_t;
#define	BEAM_ID	  0
#define	BMPR_ID	  1
#define	BREP_ID	  2
#define	CALB_ID	  3
#define	CALL_ID	  4
#define	CC01_ID	  5
#define	CC1_ID	  6
#define	CC_ID	  7
#define	CCDI_ID	  8
#define	CCH_ID	  9
#define	CCMT_ID	 10
#define	CCPB_ID	 11
#define	CCPE_ID	 12
#define	CCRC_ID	 13
#define	CCS_ID	 14
#define	CCT_ID	 15
#define	CHI2_ID	 16
#define	CL01_ID	 17
#define	CLST_ID	 18
#define	CPED_ID	 19
#define	CSQL_ID	 20
#define	DC0_ID	 21
#define	DC1_ID	 22
#define	DCDW_ID	 23
#define	DCGM_ID	 24
#define	DCGW_ID	 25
#define	DCH_ID	 26
#define	DCMN_ID	 27
#define	DCPB_ID	 28
#define	DCV1_ID	 29
#define	DCV2_ID	 30
#define	DCV3_ID	 31
#define	DDLY_ID	 32
#define	DGEO_ID	 33
#define	DHCL_ID	 34
#define	DITM_ID	 35
#define	DOCA_ID	 36
#define	DPCP_ID	 37
#define	DPSP_ID	 38
#define	DSPC_ID	 39
#define	DSPS_ID	 40
#define	DSTC_ID	 41
#define	DTCP_ID	 42
#define	DTRK_ID	 43
#define	EC01_ID	 44
#define	EC1_ID	 45
#define	EC1P_ID	 46
#define	EC1R_ID	 47
#define	ECCA_ID	 48
#define	ECCL_ID	 49
#define	ECCT_ID	 50
#define	EC_ID	 51
#define	ECDI_ID	 52
#define	ECG_ID	 53
#define	ECHB_ID	 54
#define	ECH_ID	 55
#define	ECMT_ID	 56
#define	ECP1_ID	 57
#define	ECPB_ID	 58
#define	ECPC_ID	 59
#define	ECP_ID	 60
#define	ECPE_ID	 61
#define	ECPI_ID	 62
#define	ECPO_ID	 63
#define	ECRB_ID	 64
#define	ECS_ID	 65
#define	ECT_ID	 66
#define	EID0_ID	 67
#define	EPIC_ID	 68
#define	EVNT_ID	 69
#define	FBPM_ID	 70
#define	G1SL_ID	 71
#define	G2SL_ID	 72
#define	G3SL_ID	 73
#define	G4SL_ID	 74
#define	GPAR_ID	 75
#define	GPID_ID	 76
#define	GP_X_ID	 77
#define	GP_Y_ID	 78
#define	HBER_ID	 79
#define	HBID_ID	 80
#define	HBLA_ID	 81
#define	HBTB_ID	 82
#define	HBTR_ID	 83
#define	HCAL_ID	 84
#define	HDPL_ID	 85
#define	HEAD_ID	 86
#define	HEVT_ID	 87
#define	HLS_ID	 88
#define	IC_ID	 89
#define	ICHB_ID	 90
#define	KFIT_ID	 91
#define	L1PG_ID	 92
#define	L2H_ID	 93
#define	L2S_ID	 94
#define	LASR_ID	 95
#define	LCDI_ID	 96
#define	LCPB_ID	 97
#define	LOWQ_ID	 98
#define	MCEV_ID	 99
#define	MCHD_ID	100
#define	MCTK_ID	101
#define	MCVX_ID	102
#define	MS1_ID	103
#define	MTRK_ID	104
#define	MVRT_ID	105
#define	PART_ID	106
#define	PCO_ID	107
#define	PHTM_ID	108
#define	PID1_ID	109
#define	PID2_ID	110
#define	PID3_ID	111
#define	PID4_ID	112
#define	PID5_ID	113
#define	PID6_ID	114
#define	PIDT_ID	115
#define	PRTM_ID	116
#define	PSO_ID	117
#define	PTDB_ID	118
#define	RCST_ID	119
#define	REF_ID	120
#define	RF_ID	121
#define	RFT_ID	122
#define	RGLK_ID	123
#define	RNLG_ID	124
#define	RNPE_ID	125
#define	RTSL_ID	126
#define	RUNC_ID	127
#define	S1ST_ID	128
#define	SC1_ID	129
#define	SCC_ID	130
#define	SC_ID	131
#define	SCDI_ID	132
#define	SCGD_ID	133
#define	SCG_ID	134
#define	SCH_ID	135
#define	SCMD_ID	136
#define	SCMT_ID	137
#define	SCMW_ID	138
#define	SCPB_ID	139
#define	SCP_ID	140
#define	SCPE_ID	141
#define	SCPS_ID	142
#define	SCRC_ID	143
#define	SCR_ID	144
#define	SCS_ID	145
#define	SCT_ID	146
#define	SGMP_ID	147
#define	SPAR_ID	148
#define	SPIN_ID	149
#define	ST1_ID	150
#define	ST_ID	151
#define	STG_ID	152
#define	STH_ID	153
#define	STN0_ID	154
#define	STN1_ID	155
#define	STPB_ID	156
#define	STPE_ID	157
#define	STR_ID	158
#define	STS_ID	159
#define	STSN_ID	160
#define	SYNC_ID	161
#define	TACO_ID	162
#define	TAGE_ID	163
#define	TAGI_ID	164
#define	TAGM_ID	165
#define	TAGR_ID	166
#define	TAGT_ID	167
#define	TATL_ID	168
#define	TATR_ID	169
#define	TBER_ID	170
#define	TBID_ID	171
#define	TBLA_ID	172
#define	TBTR_ID	173
#define	TCSB_ID	174
#define	TCT_ID	175
#define	TDPL_ID	176
#define	TESC_ID	177
#define	TGBI_ID	178
#define	TGEO_ID	179
#define	TGPB_ID	180
#define	TGS_ID	181
#define	TGTL_ID	182
#define	TGTR_ID	183
#define	TLV1_ID	184
#define	TRCF_ID	185
#define	TRGS_ID	186
#define	TRKS_ID	187
#define	TRL1_ID	188
#define	TRPB_ID	189
#define	TSPR_ID	190
#define	TSRG_ID	191
#define	UNUS_ID	192
#define	VERT_ID	193
#define	RC26_ID	194
#define	STRE_ID	195

#define	NO_OF_BANKS	196

#endif
