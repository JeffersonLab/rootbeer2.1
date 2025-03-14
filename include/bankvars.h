#ifndef ROOTBEER_BANKVARS
#define ROOTBEER_BANKVARS

// WARNING - This file is generated from clasbanks.ddl by the command
//           Using the script bank2struct.gk
//	      Edit at your peril!!

#include "bankheader.h"
#ifndef READBOS
#include "TBos.h"
#include "TDST.h"
#endif
#include <time.h>

// bank BEAM - multi sector bank
extern struct BEAM_t	**BEAM;
extern int 		 *BEAM_NH;
extern int 		 *BEAM_S;
extern int 		  BEAM_NS;
extern void BEAM_Print();

// bank BMPR - single sector bank
extern struct BMPR_t	 *BMPR;
extern int 		  BMPR_NH;
extern int 		  BMPR_S;
extern void BMPR_Print();

// bank BREP - single sector bank
extern struct BREP_t	 *BREP;
extern int 		  BREP_NH;
extern int 		  BREP_S;
extern void BREP_Print();

// bank CALB - multi sector bank
extern struct CALB_t	**CALB;
extern int 		 *CALB_NH;
extern int 		 *CALB_S;
extern int 		  CALB_NS;
extern void CALB_Print();

// bank CALL - single sector bank
extern struct CALL_t	 *CALL;
extern int 		  CALL_NH;
extern int 		  CALL_S;
extern void CALL_Print();

// bank CC01 - multi sector bank
extern struct CC01_t	**CC01;
extern int 		 *CC01_NH;
extern int 		 *CC01_S;
extern int 		  CC01_NS;
extern void CC01_Print();

// bank CC1 - multi sector bank
extern struct CC1_t	**CC1;
extern int 		 *CC1_NH;
extern int 		 *CC1_S;
extern int 		  CC1_NS;
extern void CC1_Print();

// bank CC - multi sector bank
extern struct CC_t	**CC;
extern int 		 *CC_NH;
extern int 		 *CC_S;
extern int 		  CC_NS;
extern void CC_Print();

// bank CCDI - multi sector bank
extern struct CCDI_t	**CCDI;
extern int 		 *CCDI_NH;
extern int 		 *CCDI_S;
extern int 		  CCDI_NS;
extern void CCDI_Print();

// bank CCH - multi sector bank
extern struct CCH_t	**CCH;
extern int 		 *CCH_NH;
extern int 		 *CCH_S;
extern int 		  CCH_NS;
extern void CCH_Print();

// bank CCMT - multi sector bank
extern struct CCMT_t	**CCMT;
extern int 		 *CCMT_NH;
extern int 		 *CCMT_S;
extern int 		  CCMT_NS;
extern void CCMT_Print();

// bank CCPB - single sector bank
extern struct CCPB_t	 *CCPB;
extern int 		  CCPB_NH;
extern int 		  CCPB_S;
extern void CCPB_Print();

// bank CCPE - multi sector bank
extern struct CCPE_t	**CCPE;
extern int 		 *CCPE_NH;
extern int 		 *CCPE_S;
extern int 		  CCPE_NS;
extern void CCPE_Print();

// bank CCRC - single sector bank
extern struct CCRC_t	 *CCRC;
extern int 		  CCRC_NH;
extern int 		  CCRC_S;
extern void CCRC_Print();

// bank CCS - multi sector bank
extern struct CCS_t	**CCS;
extern int 		 *CCS_NH;
extern int 		 *CCS_S;
extern int 		  CCS_NS;
extern void CCS_Print();

// bank CCT - multi sector bank
extern struct CCT_t	**CCT;
extern int 		 *CCT_NH;
extern int 		 *CCT_S;
extern int 		  CCT_NS;
extern void CCT_Print();

// bank CHI2 - multi sector bank
extern struct CHI2_t	**CHI2;
extern int 		 *CHI2_NH;
extern int 		 *CHI2_S;
extern int 		  CHI2_NS;
extern void CHI2_Print();

// bank CL01 - single sector bank
extern struct CL01_t	 *CL01;
extern int 		  CL01_NH;
extern int 		  CL01_S;
extern void CL01_Print();

// bank CLST - multi sector bank
extern struct CLST_t	**CLST;
extern int 		 *CLST_NH;
extern int 		 *CLST_S;
extern int 		  CLST_NS;
extern void CLST_Print();

// bank CPED - multi sector bank
extern struct CPED_t	**CPED;
extern int 		 *CPED_NH;
extern int 		 *CPED_S;
extern int 		  CPED_NS;
extern void CPED_Print();

// bank CSQL - multi sector bank
extern struct CSQL_t	**CSQL;
extern int 		 *CSQL_NH;
extern int 		 *CSQL_S;
extern int 		  CSQL_NS;
extern void CSQL_Print();

// bank DC0 - multi sector bank
extern struct DC0_t	**DC0;
extern int 		 *DC0_NH;
extern int 		 *DC0_S;
extern int 		  DC0_NS;
extern void DC0_Print();

// bank DC1 - multi sector bank
extern struct DC1_t	**DC1;
extern int 		 *DC1_NH;
extern int 		 *DC1_S;
extern int 		  DC1_NS;
extern void DC1_Print();

// bank DCDW - multi sector bank
extern struct DCDW_t	**DCDW;
extern int 		 *DCDW_NH;
extern int 		 *DCDW_S;
extern int 		  DCDW_NS;
extern void DCDW_Print();

// bank DCGM - multi sector bank
extern struct DCGM_t	**DCGM;
extern int 		 *DCGM_NH;
extern int 		 *DCGM_S;
extern int 		  DCGM_NS;
extern void DCGM_Print();

// bank DCGW - multi sector bank
extern struct DCGW_t	**DCGW;
extern int 		 *DCGW_NH;
extern int 		 *DCGW_S;
extern int 		  DCGW_NS;
extern void DCGW_Print();

// bank DCH - multi sector bank
extern struct DCH_t	**DCH;
extern int 		 *DCH_NH;
extern int 		 *DCH_S;
extern int 		  DCH_NS;
extern void DCH_Print();

// bank DCMN - multi sector bank
extern struct DCMN_t	**DCMN;
extern int 		 *DCMN_NH;
extern int 		 *DCMN_S;
extern int 		  DCMN_NS;
extern void DCMN_Print();

// bank DCPB - single sector bank
extern struct DCPB_t	 *DCPB;
extern int 		  DCPB_NH;
extern int 		  DCPB_S;
extern void DCPB_Print();

// bank DCV1 - multi sector bank
extern struct DCV1_t	**DCV1;
extern int 		 *DCV1_NH;
extern int 		 *DCV1_S;
extern int 		  DCV1_NS;
extern void DCV1_Print();

// bank DCV2 - multi sector bank
extern struct DCV2_t	**DCV2;
extern int 		 *DCV2_NH;
extern int 		 *DCV2_S;
extern int 		  DCV2_NS;
extern void DCV2_Print();

// bank DCV3 - multi sector bank
extern struct DCV3_t	**DCV3;
extern int 		 *DCV3_NH;
extern int 		 *DCV3_S;
extern int 		  DCV3_NS;
extern void DCV3_Print();

// bank DDLY - multi sector bank
extern struct DDLY_t	**DDLY;
extern int 		 *DDLY_NH;
extern int 		 *DDLY_S;
extern int 		  DDLY_NS;
extern void DDLY_Print();

// bank DGEO - multi sector bank
extern struct DGEO_t	**DGEO;
extern int 		 *DGEO_NH;
extern int 		 *DGEO_S;
extern int 		  DGEO_NS;
extern void DGEO_Print();

// bank DHCL - multi sector bank
extern struct DHCL_t	**DHCL;
extern int 		 *DHCL_NH;
extern int 		 *DHCL_S;
extern int 		  DHCL_NS;
extern void DHCL_Print();

// bank DITM - multi sector bank
extern struct DITM_t	**DITM;
extern int 		 *DITM_NH;
extern int 		 *DITM_S;
extern int 		  DITM_NS;
extern void DITM_Print();

// bank DOCA - multi sector bank
extern struct DOCA_t	**DOCA;
extern int 		 *DOCA_NH;
extern int 		 *DOCA_S;
extern int 		  DOCA_NS;
extern void DOCA_Print();

// bank DPCP - multi sector bank
extern struct DPCP_t	**DPCP;
extern int 		 *DPCP_NH;
extern int 		 *DPCP_S;
extern int 		  DPCP_NS;
extern void DPCP_Print();

// bank DPSP - multi sector bank
extern struct DPSP_t	**DPSP;
extern int 		 *DPSP_NH;
extern int 		 *DPSP_S;
extern int 		  DPSP_NS;
extern void DPSP_Print();

// bank DSPC - multi sector bank
extern struct DSPC_t	**DSPC;
extern int 		 *DSPC_NH;
extern int 		 *DSPC_S;
extern int 		  DSPC_NS;
extern void DSPC_Print();

// bank DSPS - multi sector bank
extern struct DSPS_t	**DSPS;
extern int 		 *DSPS_NH;
extern int 		 *DSPS_S;
extern int 		  DSPS_NS;
extern void DSPS_Print();

// bank DSTC - single sector bank
extern struct DSTC_t	 *DSTC;
extern int 		  DSTC_NH;
extern int 		  DSTC_S;
extern void DSTC_Print();

// bank DTCP - multi sector bank
extern struct DTCP_t	**DTCP;
extern int 		 *DTCP_NH;
extern int 		 *DTCP_S;
extern int 		  DTCP_NS;
extern void DTCP_Print();

// bank DTRK - multi sector bank
extern struct DTRK_t	**DTRK;
extern int 		 *DTRK_NH;
extern int 		 *DTRK_S;
extern int 		  DTRK_NS;
extern void DTRK_Print();

// bank EC01 - multi sector bank
extern struct EC01_t	**EC01;
extern int 		 *EC01_NH;
extern int 		 *EC01_S;
extern int 		  EC01_NS;
extern void EC01_Print();

// bank EC1 - multi sector bank
extern struct EC1_t	**EC1;
extern int 		 *EC1_NH;
extern int 		 *EC1_S;
extern int 		  EC1_NS;
extern void EC1_Print();

// bank EC1P - multi sector bank
extern struct EC1P_t	**EC1P;
extern int 		 *EC1P_NH;
extern int 		 *EC1P_S;
extern int 		  EC1P_NS;
extern void EC1P_Print();

// bank EC1R - multi sector bank
extern struct EC1R_t	**EC1R;
extern int 		 *EC1R_NH;
extern int 		 *EC1R_S;
extern int 		  EC1R_NS;
extern void EC1R_Print();

// bank ECCA - multi sector bank
extern struct ECCA_t	**ECCA;
extern int 		 *ECCA_NH;
extern int 		 *ECCA_S;
extern int 		  ECCA_NS;
extern void ECCA_Print();

// bank ECCL - multi sector bank
extern struct ECCL_t	**ECCL;
extern int 		 *ECCL_NH;
extern int 		 *ECCL_S;
extern int 		  ECCL_NS;
extern void ECCL_Print();

// bank ECCT - multi sector bank
extern struct ECCT_t	**ECCT;
extern int 		 *ECCT_NH;
extern int 		 *ECCT_S;
extern int 		  ECCT_NS;
extern void ECCT_Print();

// bank EC - multi sector bank
extern struct EC_t	**EC;
extern int 		 *EC_NH;
extern int 		 *EC_S;
extern int 		  EC_NS;
extern void EC_Print();

// bank ECDI - multi sector bank
extern struct ECDI_t	**ECDI;
extern int 		 *ECDI_NH;
extern int 		 *ECDI_S;
extern int 		  ECDI_NS;
extern void ECDI_Print();

// bank ECG - multi sector bank
extern struct ECG_t	**ECG;
extern int 		 *ECG_NH;
extern int 		 *ECG_S;
extern int 		  ECG_NS;
extern void ECG_Print();

// bank ECHB - single sector bank
extern struct ECHB_t	 *ECHB;
extern int 		  ECHB_NH;
extern int 		  ECHB_S;
extern void ECHB_Print();

// bank ECH - multi sector bank
extern struct ECH_t	**ECH;
extern int 		 *ECH_NH;
extern int 		 *ECH_S;
extern int 		  ECH_NS;
extern void ECH_Print();

// bank ECMT - multi sector bank
extern struct ECMT_t	**ECMT;
extern int 		 *ECMT_NH;
extern int 		 *ECMT_S;
extern int 		  ECMT_NS;
extern void ECMT_Print();

// bank ECP1 - multi sector bank
extern struct ECP1_t	**ECP1;
extern int 		 *ECP1_NH;
extern int 		 *ECP1_S;
extern int 		  ECP1_NS;
extern void ECP1_Print();

// bank ECPB - single sector bank
extern struct ECPB_t	 *ECPB;
extern int 		  ECPB_NH;
extern int 		  ECPB_S;
extern void ECPB_Print();

// bank ECPC - multi sector bank
extern struct ECPC_t	**ECPC;
extern int 		 *ECPC_NH;
extern int 		 *ECPC_S;
extern int 		  ECPC_NS;
extern void ECPC_Print();

// bank ECP - multi sector bank
extern struct ECP_t	**ECP;
extern int 		 *ECP_NH;
extern int 		 *ECP_S;
extern int 		  ECP_NS;
extern void ECP_Print();

// bank ECPE - multi sector bank
extern struct ECPE_t	**ECPE;
extern int 		 *ECPE_NH;
extern int 		 *ECPE_S;
extern int 		  ECPE_NS;
extern void ECPE_Print();

// bank ECPI - single sector bank
extern struct ECPI_t	 *ECPI;
extern int 		  ECPI_NH;
extern int 		  ECPI_S;
extern void ECPI_Print();

// bank ECPO - multi sector bank
extern struct ECPO_t	**ECPO;
extern int 		 *ECPO_NH;
extern int 		 *ECPO_S;
extern int 		  ECPO_NS;
extern void ECPO_Print();

// bank ECRB - multi sector bank
extern struct ECRB_t	**ECRB;
extern int 		 *ECRB_NH;
extern int 		 *ECRB_S;
extern int 		  ECRB_NS;
extern void ECRB_Print();

// bank ECS - multi sector bank
extern struct ECS_t	**ECS;
extern int 		 *ECS_NH;
extern int 		 *ECS_S;
extern int 		  ECS_NS;
extern void ECS_Print();

// bank ECT - multi sector bank
extern struct ECT_t	**ECT;
extern int 		 *ECT_NH;
extern int 		 *ECT_S;
extern int 		  ECT_NS;
extern void ECT_Print();

// bank EID0 - multi sector bank
extern struct EID0_t	**EID0;
extern int 		 *EID0_NH;
extern int 		 *EID0_S;
extern int 		  EID0_NS;
extern void EID0_Print();

// bank EPIC - multi sector bank
extern struct EPIC_t	**EPIC;
extern int 		 *EPIC_NH;
extern int 		 *EPIC_S;
extern int 		  EPIC_NS;
extern void EPIC_Print();

// bank EVNT - single sector bank
extern struct EVNT_t	 *EVNT;
extern int 		  EVNT_NH;
extern int 		  EVNT_S;
extern void EVNT_Print();

// bank FBPM - multi sector bank
extern struct FBPM_t	**FBPM;
extern int 		 *FBPM_NH;
extern int 		 *FBPM_S;
extern int 		  FBPM_NS;
extern void FBPM_Print();

// bank G1SL - multi sector bank
extern struct G1SL_t	**G1SL;
extern int 		 *G1SL_NH;
extern int 		 *G1SL_S;
extern int 		  G1SL_NS;
extern void G1SL_Print();

// bank G2SL - multi sector bank
extern struct G2SL_t	**G2SL;
extern int 		 *G2SL_NH;
extern int 		 *G2SL_S;
extern int 		  G2SL_NS;
extern void G2SL_Print();

// bank G3SL - multi sector bank
extern struct G3SL_t	**G3SL;
extern int 		 *G3SL_NH;
extern int 		 *G3SL_S;
extern int 		  G3SL_NS;
extern void G3SL_Print();

// bank G4SL - multi sector bank
extern struct G4SL_t	**G4SL;
extern int 		 *G4SL_NH;
extern int 		 *G4SL_S;
extern int 		  G4SL_NS;
extern void G4SL_Print();

// bank GPAR - multi sector bank
extern struct GPAR_t	**GPAR;
extern int 		 *GPAR_NH;
extern int 		 *GPAR_S;
extern int 		  GPAR_NS;
extern void GPAR_Print();

// bank GPID - single sector bank
extern struct GPID_t	 *GPID;
extern int 		  GPID_NH;
extern int 		  GPID_S;
extern void GPID_Print();

// bank GP_X - single sector bank
extern struct GP_X_t	 *GP_X;
extern int 		  GP_X_NH;
extern int 		  GP_X_S;
extern void GP_X_Print();

// bank GP_Y - single sector bank
extern struct GP_Y_t	 *GP_Y;
extern int 		  GP_Y_NH;
extern int 		  GP_Y_S;
extern void GP_Y_Print();

// bank HBER - single sector bank
extern struct HBER_t	 *HBER;
extern int 		  HBER_NH;
extern int 		  HBER_S;
extern void HBER_Print();

// bank HBID - single sector bank
extern struct HBID_t	 *HBID;
extern int 		  HBID_NH;
extern int 		  HBID_S;
extern void HBID_Print();

// bank HBLA - multi sector bank
extern struct HBLA_t	**HBLA;
extern int 		 *HBLA_NH;
extern int 		 *HBLA_S;
extern int 		  HBLA_NS;
extern void HBLA_Print();

// bank HBTB - multi sector bank
extern struct HBTB_t	**HBTB;
extern int 		 *HBTB_NH;
extern int 		 *HBTB_S;
extern int 		  HBTB_NS;
extern void HBTB_Print();

// bank HBTR - single sector bank
extern struct HBTR_t	 *HBTR;
extern int 		  HBTR_NH;
extern int 		  HBTR_S;
extern void HBTR_Print();

// bank HCAL - multi sector bank
extern struct HCAL_t	**HCAL;
extern int 		 *HCAL_NH;
extern int 		 *HCAL_S;
extern int 		  HCAL_NS;
extern void HCAL_Print();

// bank HDPL - multi sector bank
extern struct HDPL_t	**HDPL;
extern int 		 *HDPL_NH;
extern int 		 *HDPL_S;
extern int 		  HDPL_NS;
extern void HDPL_Print();

// bank HEAD - single sector bank
extern struct HEAD_t	 *HEAD;
extern int 		  HEAD_NH;
extern int 		  HEAD_S;
extern void HEAD_Print();

// bank HEVT - single sector bank
extern struct HEVT_t	 *HEVT;
extern int 		  HEVT_NH;
extern int 		  HEVT_S;
extern void HEVT_Print();

// bank HLS - multi sector bank
extern struct HLS_t	**HLS;
extern int 		 *HLS_NH;
extern int 		 *HLS_S;
extern int 		  HLS_NS;
extern void HLS_Print();

// bank IC - multi sector bank
extern struct IC_t	**IC;
extern int 		 *IC_NH;
extern int 		 *IC_S;
extern int 		  IC_NS;
extern void IC_Print();

// bank ICHB - multi sector bank
extern struct ICHB_t	**ICHB;
extern int 		 *ICHB_NH;
extern int 		 *ICHB_S;
extern int 		  ICHB_NS;
extern void ICHB_Print();

// bank KFIT - multi sector bank
extern struct KFIT_t	**KFIT;
extern int 		 *KFIT_NH;
extern int 		 *KFIT_S;
extern int 		  KFIT_NS;
extern void KFIT_Print();

// bank L1PG - multi sector bank
extern struct L1PG_t	**L1PG;
extern int 		 *L1PG_NH;
extern int 		 *L1PG_S;
extern int 		  L1PG_NS;
extern void L1PG_Print();

// bank L2H - multi sector bank
extern struct L2H_t	**L2H;
extern int 		 *L2H_NH;
extern int 		 *L2H_S;
extern int 		  L2H_NS;
extern void L2H_Print();

// bank L2S - multi sector bank
extern struct L2S_t	**L2S;
extern int 		 *L2S_NH;
extern int 		 *L2S_S;
extern int 		  L2S_NS;
extern void L2S_Print();

// bank LASR - multi sector bank
extern struct LASR_t	**LASR;
extern int 		 *LASR_NH;
extern int 		 *LASR_S;
extern int 		  LASR_NS;
extern void LASR_Print();

// bank LCDI - multi sector bank
extern struct LCDI_t	**LCDI;
extern int 		 *LCDI_NH;
extern int 		 *LCDI_S;
extern int 		  LCDI_NS;
extern void LCDI_Print();

// bank LCPB - single sector bank
extern struct LCPB_t	 *LCPB;
extern int 		  LCPB_NH;
extern int 		  LCPB_S;
extern void LCPB_Print();

// bank LOWQ - multi sector bank
extern struct LOWQ_t	**LOWQ;
extern int 		 *LOWQ_NH;
extern int 		 *LOWQ_S;
extern int 		  LOWQ_NS;
extern void LOWQ_Print();

// bank MCEV - multi sector bank
extern struct MCEV_t	**MCEV;
extern int 		 *MCEV_NH;
extern int 		 *MCEV_S;
extern int 		  MCEV_NS;
extern void MCEV_Print();

// bank MCHD - multi sector bank
extern struct MCHD_t	**MCHD;
extern int 		 *MCHD_NH;
extern int 		 *MCHD_S;
extern int 		  MCHD_NS;
extern void MCHD_Print();

// bank MCTK - multi sector bank
extern struct MCTK_t	**MCTK;
extern int 		 *MCTK_NH;
extern int 		 *MCTK_S;
extern int 		  MCTK_NS;
extern void MCTK_Print();

// bank MCVX - multi sector bank
extern struct MCVX_t	**MCVX;
extern int 		 *MCVX_NH;
extern int 		 *MCVX_S;
extern int 		  MCVX_NS;
extern void MCVX_Print();

// bank MS1 - multi sector bank
extern struct MS1_t	**MS1;
extern int 		 *MS1_NH;
extern int 		 *MS1_S;
extern int 		  MS1_NS;
extern void MS1_Print();

// bank MTRK - multi sector bank
extern struct MTRK_t	**MTRK;
extern int 		 *MTRK_NH;
extern int 		 *MTRK_S;
extern int 		  MTRK_NS;
extern void MTRK_Print();

// bank MVRT - single sector bank
extern struct MVRT_t	 *MVRT;
extern int 		  MVRT_NH;
extern int 		  MVRT_S;
extern void MVRT_Print();

// bank PART - single sector bank
extern struct PART_t	 *PART;
extern int 		  PART_NH;
extern int 		  PART_S;
extern void PART_Print();

// bank PCO - multi sector bank
extern struct PCO_t	**PCO;
extern int 		 *PCO_NH;
extern int 		 *PCO_S;
extern int 		  PCO_NS;
extern void PCO_Print();

// bank PHTM - multi sector bank
extern struct PHTM_t	**PHTM;
extern int 		 *PHTM_NH;
extern int 		 *PHTM_S;
extern int 		  PHTM_NS;
extern void PHTM_Print();

// bank PID1 - multi sector bank
extern struct PID1_t	**PID1;
extern int 		 *PID1_NH;
extern int 		 *PID1_S;
extern int 		  PID1_NS;
extern void PID1_Print();

// bank PID2 - multi sector bank
extern struct PID2_t	**PID2;
extern int 		 *PID2_NH;
extern int 		 *PID2_S;
extern int 		  PID2_NS;
extern void PID2_Print();

// bank PID3 - multi sector bank
extern struct PID3_t	**PID3;
extern int 		 *PID3_NH;
extern int 		 *PID3_S;
extern int 		  PID3_NS;
extern void PID3_Print();

// bank PID4 - multi sector bank
extern struct PID4_t	**PID4;
extern int 		 *PID4_NH;
extern int 		 *PID4_S;
extern int 		  PID4_NS;
extern void PID4_Print();

// bank PID5 - multi sector bank
extern struct PID5_t	**PID5;
extern int 		 *PID5_NH;
extern int 		 *PID5_S;
extern int 		  PID5_NS;
extern void PID5_Print();

// bank PID6 - multi sector bank
extern struct PID6_t	**PID6;
extern int 		 *PID6_NH;
extern int 		 *PID6_S;
extern int 		  PID6_NS;
extern void PID6_Print();

// bank PIDT - multi sector bank
extern struct PIDT_t	**PIDT;
extern int 		 *PIDT_NH;
extern int 		 *PIDT_S;
extern int 		  PIDT_NS;
extern void PIDT_Print();

// bank PRTM - multi sector bank
extern struct PRTM_t	**PRTM;
extern int 		 *PRTM_NH;
extern int 		 *PRTM_S;
extern int 		  PRTM_NS;
extern void PRTM_Print();

// bank PSO - multi sector bank
extern struct PSO_t	**PSO;
extern int 		 *PSO_NH;
extern int 		 *PSO_S;
extern int 		  PSO_NS;
extern void PSO_Print();

// bank PTDB - multi sector bank
extern struct PTDB_t	**PTDB;
extern int 		 *PTDB_NH;
extern int 		 *PTDB_S;
extern int 		  PTDB_NS;
extern void PTDB_Print();

// bank RCST - multi sector bank
extern struct RCST_t	**RCST;
extern int 		 *RCST_NH;
extern int 		 *RCST_S;
extern int 		  RCST_NS;
extern void RCST_Print();

// bank REF - single sector bank
extern struct REF_t	 *REF;
extern int 		  REF_NH;
extern int 		  REF_S;
extern void REF_Print();

// bank RF - multi sector bank
extern struct RF_t	**RF;
extern int 		 *RF_NH;
extern int 		 *RF_S;
extern int 		  RF_NS;
extern void RF_Print();

// bank RFT - single sector bank
extern struct RFT_t	 *RFT;
extern int 		  RFT_NH;
extern int 		  RFT_S;
extern void RFT_Print();

// bank RGLK - multi sector bank
extern struct RGLK_t	**RGLK;
extern int 		 *RGLK_NH;
extern int 		 *RGLK_S;
extern int 		  RGLK_NS;
extern void RGLK_Print();

// bank RNLG - multi sector bank
extern struct RNLG_t	**RNLG;
extern int 		 *RNLG_NH;
extern int 		 *RNLG_S;
extern int 		  RNLG_NS;
extern void RNLG_Print();

// bank RNPE - multi sector bank
extern struct RNPE_t	**RNPE;
extern int 		 *RNPE_NH;
extern int 		 *RNPE_S;
extern int 		  RNPE_NS;
extern void RNPE_Print();

// bank RTSL - multi sector bank
extern struct RTSL_t	**RTSL;
extern int 		 *RTSL_NH;
extern int 		 *RTSL_S;
extern int 		  RTSL_NS;
extern void RTSL_Print();

// bank RUNC - multi sector bank
extern struct RUNC_t	**RUNC;
extern int 		 *RUNC_NH;
extern int 		 *RUNC_S;
extern int 		  RUNC_NS;
extern void RUNC_Print();

// bank S1ST - multi sector bank
extern struct S1ST_t	**S1ST;
extern int 		 *S1ST_NH;
extern int 		 *S1ST_S;
extern int 		  S1ST_NS;
extern void S1ST_Print();

// bank SC1 - multi sector bank
extern struct SC1_t	**SC1;
extern int 		 *SC1_NH;
extern int 		 *SC1_S;
extern int 		  SC1_NS;
extern void SC1_Print();

// bank SCC - multi sector bank
extern struct SCC_t	**SCC;
extern int 		 *SCC_NH;
extern int 		 *SCC_S;
extern int 		  SCC_NS;
extern void SCC_Print();

// bank SC - multi sector bank
extern struct SC_t	**SC;
extern int 		 *SC_NH;
extern int 		 *SC_S;
extern int 		  SC_NS;
extern void SC_Print();

// bank SCDI - multi sector bank
extern struct SCDI_t	**SCDI;
extern int 		 *SCDI_NH;
extern int 		 *SCDI_S;
extern int 		  SCDI_NS;
extern void SCDI_Print();

// bank SCGD - multi sector bank
extern struct SCGD_t	**SCGD;
extern int 		 *SCGD_NH;
extern int 		 *SCGD_S;
extern int 		  SCGD_NS;
extern void SCGD_Print();

// bank SCG - multi sector bank
extern struct SCG_t	**SCG;
extern int 		 *SCG_NH;
extern int 		 *SCG_S;
extern int 		  SCG_NS;
extern void SCG_Print();

// bank SCH - multi sector bank
extern struct SCH_t	**SCH;
extern int 		 *SCH_NH;
extern int 		 *SCH_S;
extern int 		  SCH_NS;
extern void SCH_Print();

// bank SCMD - multi sector bank
extern struct SCMD_t	**SCMD;
extern int 		 *SCMD_NH;
extern int 		 *SCMD_S;
extern int 		  SCMD_NS;
extern void SCMD_Print();

// bank SCMT - multi sector bank
extern struct SCMT_t	**SCMT;
extern int 		 *SCMT_NH;
extern int 		 *SCMT_S;
extern int 		  SCMT_NS;
extern void SCMT_Print();

// bank SCMW - multi sector bank
extern struct SCMW_t	**SCMW;
extern int 		 *SCMW_NH;
extern int 		 *SCMW_S;
extern int 		  SCMW_NS;
extern void SCMW_Print();

// bank SCPB - single sector bank
extern struct SCPB_t	 *SCPB;
extern int 		  SCPB_NH;
extern int 		  SCPB_S;
extern void SCPB_Print();

// bank SCP - multi sector bank
extern struct SCP_t	**SCP;
extern int 		 *SCP_NH;
extern int 		 *SCP_S;
extern int 		  SCP_NS;
extern void SCP_Print();

// bank SCPE - multi sector bank
extern struct SCPE_t	**SCPE;
extern int 		 *SCPE_NH;
extern int 		 *SCPE_S;
extern int 		  SCPE_NS;
extern void SCPE_Print();

// bank SCPS - multi sector bank
extern struct SCPS_t	**SCPS;
extern int 		 *SCPS_NH;
extern int 		 *SCPS_S;
extern int 		  SCPS_NS;
extern void SCPS_Print();

// bank SCRC - multi sector bank
extern struct SCRC_t	**SCRC;
extern int 		 *SCRC_NH;
extern int 		 *SCRC_S;
extern int 		  SCRC_NS;
extern void SCRC_Print();

// bank SCR - multi sector bank
extern struct SCR_t	**SCR;
extern int 		 *SCR_NH;
extern int 		 *SCR_S;
extern int 		  SCR_NS;
extern void SCR_Print();

// bank SCS - multi sector bank
extern struct SCS_t	**SCS;
extern int 		 *SCS_NH;
extern int 		 *SCS_S;
extern int 		  SCS_NS;
extern void SCS_Print();

// bank SCT - multi sector bank
extern struct SCT_t	**SCT;
extern int 		 *SCT_NH;
extern int 		 *SCT_S;
extern int 		  SCT_NS;
extern void SCT_Print();

// bank SGMP - multi sector bank
extern struct SGMP_t	**SGMP;
extern int 		 *SGMP_NH;
extern int 		 *SGMP_S;
extern int 		  SGMP_NS;
extern void SGMP_Print();

// bank SPAR - multi sector bank
extern struct SPAR_t	**SPAR;
extern int 		 *SPAR_NH;
extern int 		 *SPAR_S;
extern int 		  SPAR_NS;
extern void SPAR_Print();

// bank SPIN - multi sector bank
extern struct SPIN_t	**SPIN;
extern int 		 *SPIN_NH;
extern int 		 *SPIN_S;
extern int 		  SPIN_NS;
extern void SPIN_Print();

// bank ST1 - multi sector bank
extern struct ST1_t	**ST1;
extern int 		 *ST1_NH;
extern int 		 *ST1_S;
extern int 		  ST1_NS;
extern void ST1_Print();

// bank ST - single sector bank
extern struct ST_t	 *ST;
extern int 		  ST_NH;
extern int 		  ST_S;
extern void ST_Print();

// bank STG - multi sector bank
extern struct STG_t	**STG;
extern int 		 *STG_NH;
extern int 		 *STG_S;
extern int 		  STG_NS;
extern void STG_Print();

// bank STH - multi sector bank
extern struct STH_t	**STH;
extern int 		 *STH_NH;
extern int 		 *STH_S;
extern int 		  STH_NS;
extern void STH_Print();

// bank STN0 - single sector bank
extern struct STN0_t	 *STN0;
extern int 		  STN0_NH;
extern int 		  STN0_S;
extern void STN0_Print();

// bank STN1 - single sector bank
extern struct STN1_t	 *STN1;
extern int 		  STN1_NH;
extern int 		  STN1_S;
extern void STN1_Print();

// bank STPB - single sector bank
extern struct STPB_t	 *STPB;
extern int 		  STPB_NH;
extern int 		  STPB_S;
extern void STPB_Print();

// bank STPE - multi sector bank
extern struct STPE_t	**STPE;
extern int 		 *STPE_NH;
extern int 		 *STPE_S;
extern int 		  STPE_NS;
extern void STPE_Print();

// bank STR - single sector bank
extern struct STR_t	 *STR;
extern int 		  STR_NH;
extern int 		  STR_S;
extern void STR_Print();

// bank STS - multi sector bank
extern struct STS_t	**STS;
extern int 		 *STS_NH;
extern int 		 *STS_S;
extern int 		  STS_NS;
extern void STS_Print();

// bank STSN - multi sector bank
extern struct STSN_t	**STSN;
extern int 		 *STSN_NH;
extern int 		 *STSN_S;
extern int 		  STSN_NS;
extern void STSN_Print();

// bank SYNC - multi sector bank
extern struct SYNC_t	**SYNC;
extern int 		 *SYNC_NH;
extern int 		 *SYNC_S;
extern int 		  SYNC_NS;
extern void SYNC_Print();

// bank TACO - multi sector bank
extern struct TACO_t	**TACO;
extern int 		 *TACO_NH;
extern int 		 *TACO_S;
extern int 		  TACO_NS;
extern void TACO_Print();

// bank TAGE - single sector bank
extern struct TAGE_t	 *TAGE;
extern int 		  TAGE_NH;
extern int 		  TAGE_S;
extern void TAGE_Print();

// bank TAGI - multi sector bank
extern struct TAGI_t	**TAGI;
extern int 		 *TAGI_NH;
extern int 		 *TAGI_S;
extern int 		  TAGI_NS;
extern void TAGI_Print();

// bank TAGM - single sector bank
extern struct TAGM_t	 *TAGM;
extern int 		  TAGM_NH;
extern int 		  TAGM_S;
extern void TAGM_Print();

// bank TAGR - single sector bank
extern struct TAGR_t	 *TAGR;
extern int 		  TAGR_NH;
extern int 		  TAGR_S;
extern void TAGR_Print();

// bank TAGT - single sector bank
extern struct TAGT_t	 *TAGT;
extern int 		  TAGT_NH;
extern int 		  TAGT_S;
extern void TAGT_Print();

// bank TATL - multi sector bank
extern struct TATL_t	**TATL;
extern int 		 *TATL_NH;
extern int 		 *TATL_S;
extern int 		  TATL_NS;
extern void TATL_Print();

// bank TATR - multi sector bank
extern struct TATR_t	**TATR;
extern int 		 *TATR_NH;
extern int 		 *TATR_S;
extern int 		  TATR_NS;
extern void TATR_Print();

// bank TBER - single sector bank
extern struct TBER_t	 *TBER;
extern int 		  TBER_NH;
extern int 		  TBER_S;
extern void TBER_Print();

// bank TBID - single sector bank
extern struct TBID_t	 *TBID;
extern int 		  TBID_NH;
extern int 		  TBID_S;
extern void TBID_Print();

// bank TBLA - multi sector bank
extern struct TBLA_t	**TBLA;
extern int 		 *TBLA_NH;
extern int 		 *TBLA_S;
extern int 		  TBLA_NS;
extern void TBLA_Print();

// bank TBTR - single sector bank
extern struct TBTR_t	 *TBTR;
extern int 		  TBTR_NH;
extern int 		  TBTR_S;
extern void TBTR_Print();

// bank TCSB - multi sector bank
extern struct TCSB_t	**TCSB;
extern int 		 *TCSB_NH;
extern int 		 *TCSB_S;
extern int 		  TCSB_NS;
extern void TCSB_Print();

// bank TCT - multi sector bank
extern struct TCT_t	**TCT;
extern int 		 *TCT_NH;
extern int 		 *TCT_S;
extern int 		  TCT_NS;
extern void TCT_Print();

// bank TDPL - multi sector bank
extern struct TDPL_t	**TDPL;
extern int 		 *TDPL_NH;
extern int 		 *TDPL_S;
extern int 		  TDPL_NS;
extern void TDPL_Print();

// bank TESC - single sector bank
extern struct TESC_t	 *TESC;
extern int 		  TESC_NH;
extern int 		  TESC_S;
extern void TESC_Print();

// bank TGBI - single sector bank
extern struct TGBI_t	 *TGBI;
extern int 		  TGBI_NH;
extern int 		  TGBI_S;
extern void TGBI_Print();

// bank TGEO - multi sector bank
extern struct TGEO_t	**TGEO;
extern int 		 *TGEO_NH;
extern int 		 *TGEO_S;
extern int 		  TGEO_NS;
extern void TGEO_Print();

// bank TGPB - single sector bank
extern struct TGPB_t	 *TGPB;
extern int 		  TGPB_NH;
extern int 		  TGPB_S;
extern void TGPB_Print();

// bank TGS - multi sector bank
extern struct TGS_t	**TGS;
extern int 		 *TGS_NH;
extern int 		 *TGS_S;
extern int 		  TGS_NS;
extern void TGS_Print();

// bank TGTL - single sector bank
extern struct TGTL_t	 *TGTL;
extern int 		  TGTL_NH;
extern int 		  TGTL_S;
extern void TGTL_Print();

// bank TGTR - single sector bank
extern struct TGTR_t	 *TGTR;
extern int 		  TGTR_NH;
extern int 		  TGTR_S;
extern void TGTR_Print();

// bank TLV1 - multi sector bank
extern struct TLV1_t	**TLV1;
extern int 		 *TLV1_NH;
extern int 		 *TLV1_S;
extern int 		  TLV1_NS;
extern void TLV1_Print();

// bank TRCF - multi sector bank
extern struct TRCF_t	**TRCF;
extern int 		 *TRCF_NH;
extern int 		 *TRCF_S;
extern int 		  TRCF_NS;
extern void TRCF_Print();

// bank TRGS - multi sector bank
extern struct TRGS_t	**TRGS;
extern int 		 *TRGS_NH;
extern int 		 *TRGS_S;
extern int 		  TRGS_NS;
extern void TRGS_Print();

// bank TRKS - multi sector bank
extern struct TRKS_t	**TRKS;
extern int 		 *TRKS_NH;
extern int 		 *TRKS_S;
extern int 		  TRKS_NS;
extern void TRKS_Print();

// bank TRL1 - multi sector bank
extern struct TRL1_t	**TRL1;
extern int 		 *TRL1_NH;
extern int 		 *TRL1_S;
extern int 		  TRL1_NS;
extern void TRL1_Print();

// bank TRPB - multi sector bank
extern struct TRPB_t	**TRPB;
extern int 		 *TRPB_NH;
extern int 		 *TRPB_S;
extern int 		  TRPB_NS;
extern void TRPB_Print();

// bank TSPR - multi sector bank
extern struct TSPR_t	**TSPR;
extern int 		 *TSPR_NH;
extern int 		 *TSPR_S;
extern int 		  TSPR_NS;
extern void TSPR_Print();

// bank TSRG - multi sector bank
extern struct TSRG_t	**TSRG;
extern int 		 *TSRG_NH;
extern int 		 *TSRG_S;
extern int 		  TSRG_NS;
extern void TSRG_Print();

// bank UNUS - multi sector bank
extern struct UNUS_t	**UNUS;
extern int 		 *UNUS_NH;
extern int 		 *UNUS_S;
extern int 		  UNUS_NS;
extern void UNUS_Print();

// bank VERT - multi sector bank
extern struct VERT_t	**VERT;
extern int 		 *VERT_NH;
extern int 		 *VERT_S;
extern int 		  VERT_NS;
extern void VERT_Print();

// bank RC26 - single sector bank
extern struct RC26_t	 *RC26;
extern int 		  RC26_NH;
extern int 		  RC26_S;
extern void RC26_Print();

// bank STRE - multi sector bank
extern struct STRE_t	**STRE;
extern int 		 *STRE_NH;
extern int 		 *STRE_S;
extern int 		  STRE_NS;
extern void STRE_Print();

extern int 		  nBankTypes;	// Total of different banks 
extern struct addressBanks_t *bankAddress;
#endif
