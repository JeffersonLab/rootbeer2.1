// WARNING - This file is generated from clasbanks.ddl by the command
//           Using the script bank2struct.gk
//	      Edit at your peril!!

#include "bankvars.h"

// bank BEAM - multi sector bank
struct BEAM_t	**BEAM;
int 		 *BEAM_NH;
int 		 *BEAM_S;
int 		  BEAM_NS;
void BEAM_Print(){
  int s,n;
  if(BEAM_NS){
    for(s=0;s<BEAM_NS;s++){
      if(BEAM_NH[s]){
        fprintf(stdout,"Group:\tBEAM\tSector:\t%d \tNhits:\t%d\n",BEAM_S[s],BEAM_NH[s]);
        for(n=0;n<BEAM_NH[s];n++){
          fprintf(stdout,"ENERGY: %8.3f\tITORUS: %8.3f\tIMINI: %8.3f\tITAG: %8.3f\t",
		BEAM[s][n].ENERGY,
		BEAM[s][n].ITORUS,
		BEAM[s][n].IMINI,
		BEAM[s][n].ITAG);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank BMPR - single sector bank
struct BMPR_t	 *BMPR;
int 		  BMPR_NH;
int 		  BMPR_S;
void BMPR_Print(){
  int n;
  if(BMPR_NH){
    fprintf(stdout,"Group:\tBMPR\tSector:\t%d \tNhits:\t%d\n",BMPR_S,BMPR_NH);
    for(n=0;n<BMPR_NH;n++){
      fprintf(stdout,"Q_TOT: %8.3f\tQ_TOT_LV: %8.3f\tTL_TU: %8.3f\tCURRENT: %8.3f\t",
		BMPR[n].Q_TOT,
		BMPR[n].Q_TOT_LV,
		BMPR[n].TL_TU,
		BMPR[n].CURRENT);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank BREP - single sector bank
struct BREP_t	 *BREP;
int 		  BREP_NH;
int 		  BREP_S;
void BREP_Print(){
  int n;
  if(BREP_NH){
    fprintf(stdout,"Group:\tBREP\tSector:\t%d \tNhits:\t%d\n",BREP_S,BREP_NH);
    for(n=0;n<BREP_NH;n++){
      fprintf(stdout,"MBSY2C_energy: %8.3f\tIGT0I00BIASET: %8.3f\tIGL1I00DAC2: %8.3f\tSMRPOSA: %8.3f\tSMRPOSB: %8.3f\tSMRPOSC: %8.3f\tHarp: %8.3f\thallb_sf_xy560: %8.3f\tMTSETI: %8.3f\tMTIRBCK: %8.3f\tMTVRBCK: %8.3f\tTMSETI: %8.3f\tTMIRBCK: %8.3f\tTMVRBCK: %8.3f\tCryo_pressure: %8.3f\tCryo_temperature: %8.3f\tCryo_status: %8.3f\tVCG2C24: %8.3f\tVCG2H01: %8.3f\tscalerS2o: %8.3f\tscalerS3o: %8.3f\tscalerS4o: %8.3f\tscalerS5o: %8.3f\tscalerS6o: %8.3f\tscalerS7o: %8.3f\tscalerS8o: %8.3f\tscalerS9o: %8.3f\tIPM2H01_XPOS: %8.3f\tIPM2H01_YPOS: %8.3f\tIPM2H01: %8.3f\tIPM2C24A_XPOS: %8.3f\tIPM2C24A_YPOS: %8.3f\tIPM2C24A: %8.3f\tIPM2C22A_XPOS: %8.3f\tIPM2C22A_YPOS: %8.3f\tIPM2C22A: %8.3f\t",
		BREP[n].MBSY2C_energy,
		BREP[n].IGT0I00BIASET,
		BREP[n].IGL1I00DAC2,
		BREP[n].SMRPOSA,
		BREP[n].SMRPOSB,
		BREP[n].SMRPOSC,
		BREP[n].Harp,
		BREP[n].hallb_sf_xy560,
		BREP[n].MTSETI,
		BREP[n].MTIRBCK,
		BREP[n].MTVRBCK,
		BREP[n].TMSETI,
		BREP[n].TMIRBCK,
		BREP[n].TMVRBCK,
		BREP[n].Cryo_pressure,
		BREP[n].Cryo_temperature,
		BREP[n].Cryo_status,
		BREP[n].VCG2C24,
		BREP[n].VCG2H01,
		BREP[n].scalerS2o,
		BREP[n].scalerS3o,
		BREP[n].scalerS4o,
		BREP[n].scalerS5o,
		BREP[n].scalerS6o,
		BREP[n].scalerS7o,
		BREP[n].scalerS8o,
		BREP[n].scalerS9o,
		BREP[n].IPM2H01_XPOS,
		BREP[n].IPM2H01_YPOS,
		BREP[n].IPM2H01,
		BREP[n].IPM2C24A_XPOS,
		BREP[n].IPM2C24A_YPOS,
		BREP[n].IPM2C24A,
		BREP[n].IPM2C22A_XPOS,
		BREP[n].IPM2C22A_YPOS,
		BREP[n].IPM2C22A);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank CALB - multi sector bank
struct CALB_t	**CALB;
int 		 *CALB_NH;
int 		 *CALB_S;
int 		  CALB_NS;
void CALB_Print(){
  int s,n;
  if(CALB_NS){
    for(s=0;s<CALB_NS;s++){
      if(CALB_NH[s]){
        fprintf(stdout,"Group:\tCALB\tSector:\t%d \tNhits:\t%d\n",CALB_S[s],CALB_NH[s]);
        for(n=0;n<CALB_NH[s];n++){
          fprintf(stdout,"meanRFe: %8.3f\tsigmaRFe: %8.3f\tsigmaRFh: %8.3f\tsigmaECt: %8.3f\tSFECe: %8.3f\tsigmaSF: %8.3f\tResSL1: %8.3f\tResSL2: %8.3f\tResSL3: %8.3f\tResSL4: %8.3f\tResSL5: %8.3f\tResSL6: %8.3f\t",
		CALB[s][n].meanRFe,
		CALB[s][n].sigmaRFe,
		CALB[s][n].sigmaRFh,
		CALB[s][n].sigmaECt,
		CALB[s][n].SFECe,
		CALB[s][n].sigmaSF,
		CALB[s][n].ResSL1,
		CALB[s][n].ResSL2,
		CALB[s][n].ResSL3,
		CALB[s][n].ResSL4,
		CALB[s][n].ResSL5,
		CALB[s][n].ResSL6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CALL - single sector bank
struct CALL_t	 *CALL;
int 		  CALL_NH;
int 		  CALL_S;
void CALL_Print(){
  int n;
  if(CALL_NH){
    fprintf(stdout,"Group:\tCALL\tSector:\t%d \tNhits:\t%d\n",CALL_S,CALL_NH);
    for(n=0;n<CALL_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		CALL[n].ID,
		CALL[n].TDC,
		CALL[n].ADC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank CC01 - multi sector bank
struct CC01_t	**CC01;
int 		 *CC01_NH;
int 		 *CC01_S;
int 		  CC01_NS;
void CC01_Print(){
  int s,n;
  if(CC01_NS){
    for(s=0;s<CC01_NS;s++){
      if(CC01_NH[s]){
        fprintf(stdout,"Group:\tCC01\tSector:\t%d \tNhits:\t%d\n",CC01_S[s],CC01_NH[s]);
        for(n=0;n<CC01_NH[s];n++){
          fprintf(stdout,"id: %d\ttime: %8.3f\tn_pe: %8.3f\t",
		CC01[s][n].id,
		CC01[s][n].time,
		CC01[s][n].n_pe);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CC1 - multi sector bank
struct CC1_t	**CC1;
int 		 *CC1_NH;
int 		 *CC1_S;
int 		  CC1_NS;
void CC1_Print(){
  int s,n;
  if(CC1_NS){
    for(s=0;s<CC1_NS;s++){
      if(CC1_NH[s]){
        fprintf(stdout,"Group:\tCC1\tSector:\t%d \tNhits:\t%d\n",CC1_S[s],CC1_NH[s]);
        for(n=0;n<CC1_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		CC1[s][n].ID,
		CC1[s][n].TDC,
		CC1[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CC - multi sector bank
struct CC_t	**CC;
int 		 *CC_NH;
int 		 *CC_S;
int 		  CC_NS;
void CC_Print(){
  int s,n;
  if(CC_NS){
    for(s=0;s<CC_NS;s++){
      if(CC_NH[s]){
        fprintf(stdout,"Group:\tCC\tSector:\t%d \tNhits:\t%d\n",CC_S[s],CC_NH[s]);
        for(n=0;n<CC_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		CC[s][n].ID,
		CC[s][n].TDC,
		CC[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCDI - multi sector bank
struct CCDI_t	**CCDI;
int 		 *CCDI_NH;
int 		 *CCDI_S;
int 		  CCDI_NS;
void CCDI_Print(){
  int s,n;
  if(CCDI_NS){
    for(s=0;s<CCDI_NS;s++){
      if(CCDI_NH[s]){
        fprintf(stdout,"Group:\tCCDI\tSector:\t%d \tNhits:\t%d\n",CCDI_S[s],CCDI_NH[s]);
        for(n=0;n<CCDI_NH[s];n++){
          fprintf(stdout,"crate: %d\tslot: %d\tmask: %d\tthreshold: %d\twidth: %d\t",
		CCDI[s][n].crate,
		CCDI[s][n].slot,
		CCDI[s][n].mask,
		CCDI[s][n].threshold,
		CCDI[s][n].width);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCH - multi sector bank
struct CCH_t	**CCH;
int 		 *CCH_NH;
int 		 *CCH_S;
int 		  CCH_NS;
void CCH_Print(){
  int s,n;
  if(CCH_NS){
    for(s=0;s<CCH_NS;s++){
      if(CCH_NH[s]){
        fprintf(stdout,"Group:\tCCH\tSector:\t%d \tNhits:\t%d\n",CCH_S[s],CCH_NH[s]);
        for(n=0;n<CCH_NH[s];n++){
          fprintf(stdout,"tknum: %d\tid: %d\tnhits: %d\tsector: %d\tsegment: %d\tpmom: %d\txin: %8.3f\tyin: %8.3f\tzin: %8.3f\txout: %8.3f\tyout: %8.3f\tzout: %8.3f\t",
		CCH[s][n].tknum,
		CCH[s][n].id,
		CCH[s][n].nhits,
		CCH[s][n].sector,
		CCH[s][n].segment,
		CCH[s][n].pmom,
		CCH[s][n].xin,
		CCH[s][n].yin,
		CCH[s][n].zin,
		CCH[s][n].xout,
		CCH[s][n].yout,
		CCH[s][n].zout);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCMT - multi sector bank
struct CCMT_t	**CCMT;
int 		 *CCMT_NH;
int 		 *CCMT_S;
int 		  CCMT_NS;
void CCMT_Print(){
  int s,n;
  if(CCMT_NS){
    for(s=0;s<CCMT_NS;s++){
      if(CCMT_NH[s]){
        fprintf(stdout,"Group:\tCCMT\tSector:\t%d \tNhits:\t%d\n",CCMT_S[s],CCMT_NH[s]);
        for(n=0;n<CCMT_NH[s];n++){
          fprintf(stdout,"mean_high: %d\tmean_lo: %d\t",
		CCMT[s][n].mean_high,
		CCMT[s][n].mean_lo);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCPB - single sector bank
struct CCPB_t	 *CCPB;
int 		  CCPB_NH;
int 		  CCPB_S;
void CCPB_Print(){
  int n;
  if(CCPB_NH){
    fprintf(stdout,"Group:\tCCPB\tSector:\t%d \tNhits:\t%d\n",CCPB_S,CCPB_NH);
    for(n=0;n<CCPB_NH;n++){
      fprintf(stdout,"ScSgHt: %d\tNphe: %8.3f\tTime: %8.3f\tPath: %8.3f\tChi2CC: %8.3f\tStatus: %d\t",
		CCPB[n].ScSgHt,
		CCPB[n].Nphe,
		CCPB[n].Time,
		CCPB[n].Path,
		CCPB[n].Chi2CC,
		CCPB[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank CCPE - multi sector bank
struct CCPE_t	**CCPE;
int 		 *CCPE_NH;
int 		 *CCPE_S;
int 		  CCPE_NS;
void CCPE_Print(){
  int s,n;
  if(CCPE_NS){
    for(s=0;s<CCPE_NS;s++){
      if(CCPE_NH[s]){
        fprintf(stdout,"Group:\tCCPE\tSector:\t%d \tNhits:\t%d\n",CCPE_S[s],CCPE_NH[s]);
        for(n=0;n<CCPE_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean: %d\tsigma: %d\t",
		CCPE[s][n].ID,
		CCPE[s][n].mean,
		CCPE[s][n].sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCRC - single sector bank
struct CCRC_t	 *CCRC;
int 		  CCRC_NH;
int 		  CCRC_S;
void CCRC_Print(){
  int n;
  if(CCRC_NH){
    fprintf(stdout,"Group:\tCCRC\tSector:\t%d \tNhits:\t%d\n",CCRC_S,CCRC_NH);
    for(n=0;n<CCRC_NH;n++){
      fprintf(stdout,"nrsect: %d\tnrsegm: %d\tnrsegm_p: %d\tnrsegm_m: %d\tnrphe: %d\tnrtime: %d\tnrthet: %d\tnrdthet: %d\tnrphy: %d\tnriec: %d\tnrdiec: %d\tnrstat: %d\t",
		CCRC[n].nrsect,
		CCRC[n].nrsegm,
		CCRC[n].nrsegm_p,
		CCRC[n].nrsegm_m,
		CCRC[n].nrphe,
		CCRC[n].nrtime,
		CCRC[n].nrthet,
		CCRC[n].nrdthet,
		CCRC[n].nrphy,
		CCRC[n].nriec,
		CCRC[n].nrdiec,
		CCRC[n].nrstat);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank CCS - multi sector bank
struct CCS_t	**CCS;
int 		 *CCS_NH;
int 		 *CCS_S;
int 		  CCS_NS;
void CCS_Print(){
  int s,n;
  if(CCS_NS){
    for(s=0;s<CCS_NS;s++){
      if(CCS_NH[s]){
        fprintf(stdout,"Group:\tCCS\tSector:\t%d \tNhits:\t%d\n",CCS_S[s],CCS_NH[s]);
        for(n=0;n<CCS_NH[s];n++){
          fprintf(stdout,"CCS1: %d\tCCS2: %d\tCCS3: %d\tCCS4: %d\tCCS5: %d\tCCS6: %d\tCCS7: %d\tCCS8: %d\tCCS9: %d\tCCS10: %d\tCCS11: %d\tCCS12: %d\tCCS13: %d\tCCS14: %d\tCCS15: %d\tCCS16: %d\tCCS17: %d\tCCS18: %d\tCCS19: %d\tCCS20: %d\tCCS21: %d\tCCS22: %d\tCCS23: %d\tCCS24: %d\tCCS25: %d\tCCS26: %d\tCCS27: %d\tCCS28: %d\tCCS29: %d\tCCS30: %d\tCCS31: %d\tCCS32: %d\tCCS33: %d\tCCS34: %d\tCCS35: %d\tCCS36: %d\tCCS37: %d\tCCS38: %d\tCCS39: %d\tCCS40: %d\tCCS41: %d\tCCS42: %d\tCCS43: %d\tCCS44: %d\tCCS45: %d\tCCS46: %d\tCCS47: %d\tCCS48: %d\tCCS49: %d\tCCS50: %d\tCCS51: %d\tCCS52: %d\tCCS53: %d\tCCS54: %d\tCCS55: %d\tCCS56: %d\tCCS57: %d\tCCS58: %d\tCCS59: %d\tCCS60: %d\tCCS61: %d\tCCS62: %d\tCCS63: %d\tCCS64: %d\tCCS65: %d\tCCS66: %d\tCCS67: %d\tCCS68: %d\tCCS69: %d\tCCS70: %d\tCCS71: %d\tCCS72: %d\tCCS73: %d\tCCS74: %d\tCCS75: %d\tCCS76: %d\tCCS77: %d\tCCS78: %d\tCCS79: %d\tCCS80: %d\tCCS81: %d\tCCS82: %d\tCCS83: %d\tCCS84: %d\tCCS85: %d\tCCS86: %d\tCCS87: %d\tCCS88: %d\tCCS89: %d\tCCS90: %d\tCCS91: %d\tCCS92: %d\tCCS93: %d\tCCS94: %d\tCCS95: %d\tCCS96: %d\t",
		CCS[s][n].CCS1,
		CCS[s][n].CCS2,
		CCS[s][n].CCS3,
		CCS[s][n].CCS4,
		CCS[s][n].CCS5,
		CCS[s][n].CCS6,
		CCS[s][n].CCS7,
		CCS[s][n].CCS8,
		CCS[s][n].CCS9,
		CCS[s][n].CCS10,
		CCS[s][n].CCS11,
		CCS[s][n].CCS12,
		CCS[s][n].CCS13,
		CCS[s][n].CCS14,
		CCS[s][n].CCS15,
		CCS[s][n].CCS16,
		CCS[s][n].CCS17,
		CCS[s][n].CCS18,
		CCS[s][n].CCS19,
		CCS[s][n].CCS20,
		CCS[s][n].CCS21,
		CCS[s][n].CCS22,
		CCS[s][n].CCS23,
		CCS[s][n].CCS24,
		CCS[s][n].CCS25,
		CCS[s][n].CCS26,
		CCS[s][n].CCS27,
		CCS[s][n].CCS28,
		CCS[s][n].CCS29,
		CCS[s][n].CCS30,
		CCS[s][n].CCS31,
		CCS[s][n].CCS32,
		CCS[s][n].CCS33,
		CCS[s][n].CCS34,
		CCS[s][n].CCS35,
		CCS[s][n].CCS36,
		CCS[s][n].CCS37,
		CCS[s][n].CCS38,
		CCS[s][n].CCS39,
		CCS[s][n].CCS40,
		CCS[s][n].CCS41,
		CCS[s][n].CCS42,
		CCS[s][n].CCS43,
		CCS[s][n].CCS44,
		CCS[s][n].CCS45,
		CCS[s][n].CCS46,
		CCS[s][n].CCS47,
		CCS[s][n].CCS48,
		CCS[s][n].CCS49,
		CCS[s][n].CCS50,
		CCS[s][n].CCS51,
		CCS[s][n].CCS52,
		CCS[s][n].CCS53,
		CCS[s][n].CCS54,
		CCS[s][n].CCS55,
		CCS[s][n].CCS56,
		CCS[s][n].CCS57,
		CCS[s][n].CCS58,
		CCS[s][n].CCS59,
		CCS[s][n].CCS60,
		CCS[s][n].CCS61,
		CCS[s][n].CCS62,
		CCS[s][n].CCS63,
		CCS[s][n].CCS64,
		CCS[s][n].CCS65,
		CCS[s][n].CCS66,
		CCS[s][n].CCS67,
		CCS[s][n].CCS68,
		CCS[s][n].CCS69,
		CCS[s][n].CCS70,
		CCS[s][n].CCS71,
		CCS[s][n].CCS72,
		CCS[s][n].CCS73,
		CCS[s][n].CCS74,
		CCS[s][n].CCS75,
		CCS[s][n].CCS76,
		CCS[s][n].CCS77,
		CCS[s][n].CCS78,
		CCS[s][n].CCS79,
		CCS[s][n].CCS80,
		CCS[s][n].CCS81,
		CCS[s][n].CCS82,
		CCS[s][n].CCS83,
		CCS[s][n].CCS84,
		CCS[s][n].CCS85,
		CCS[s][n].CCS86,
		CCS[s][n].CCS87,
		CCS[s][n].CCS88,
		CCS[s][n].CCS89,
		CCS[s][n].CCS90,
		CCS[s][n].CCS91,
		CCS[s][n].CCS92,
		CCS[s][n].CCS93,
		CCS[s][n].CCS94,
		CCS[s][n].CCS95,
		CCS[s][n].CCS96);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CCT - multi sector bank
struct CCT_t	**CCT;
int 		 *CCT_NH;
int 		 *CCT_S;
int 		  CCT_NS;
void CCT_Print(){
  int s,n;
  if(CCT_NS){
    for(s=0;s<CCT_NS;s++){
      if(CCT_NH[s]){
        fprintf(stdout,"Group:\tCCT\tSector:\t%d \tNhits:\t%d\n",CCT_S[s],CCT_NH[s]);
        for(n=0;n<CCT_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		CCT[s][n].ID,
		CCT[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CHI2 - multi sector bank
struct CHI2_t	**CHI2;
int 		 *CHI2_NH;
int 		 *CHI2_S;
int 		  CHI2_NS;
void CHI2_Print(){
  int s,n;
  if(CHI2_NS){
    for(s=0;s<CHI2_NS;s++){
      if(CHI2_NH[s]){
        fprintf(stdout,"Group:\tCHI2\tSector:\t%d \tNhits:\t%d\n",CHI2_S[s],CHI2_NH[s]);
        for(n=0;n<CHI2_NH[s];n++){
          fprintf(stdout,"chi2: %8.3f\tcl: %8.3f\tndf: %d\titer: %d\t",
		CHI2[s][n].chi2,
		CHI2[s][n].cl,
		CHI2[s][n].ndf,
		CHI2[s][n].iter);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CL01 - single sector bank
struct CL01_t	 *CL01;
int 		  CL01_NH;
int 		  CL01_S;
void CL01_Print(){
  int n;
  if(CL01_NH){
    fprintf(stdout,"Group:\tCL01\tSector:\t%d \tNhits:\t%d\n",CL01_S,CL01_NH);
    for(n=0;n<CL01_NH;n++){
      fprintf(stdout,"ac_amp: %d\tfc_diode_amp: %d\tfc_diode_t: %8.3f\tnc_diode_amp: %d\tnc_diode_t: %8.3f\tsc_diode_amp: %d\tsc_diode_t: %8.3f\tsf_diode_amp: %d\tsf_diode_t: %8.3f\trf1: %8.3f\trf2: %8.3f\trf: %8.3f\t",
		CL01[n].ac_amp,
		CL01[n].fc_diode_amp,
		CL01[n].fc_diode_t,
		CL01[n].nc_diode_amp,
		CL01[n].nc_diode_t,
		CL01[n].sc_diode_amp,
		CL01[n].sc_diode_t,
		CL01[n].sf_diode_amp,
		CL01[n].sf_diode_t,
		CL01[n].rf1,
		CL01[n].rf2,
		CL01[n].rf);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank CLST - multi sector bank
struct CLST_t	**CLST;
int 		 *CLST_NH;
int 		 *CLST_S;
int 		  CLST_NS;
void CLST_Print(){
  int s,n;
  if(CLST_NS){
    for(s=0;s<CLST_NS;s++){
      if(CLST_NH[s]){
        fprintf(stdout,"Group:\tCLST\tSector:\t%d \tNhits:\t%d\n",CLST_S[s],CLST_NH[s]);
        for(n=0;n<CLST_NH[s];n++){
          fprintf(stdout,"Clust: %d\t",
		CLST[s][n].Clust);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CPED - multi sector bank
struct CPED_t	**CPED;
int 		 *CPED_NH;
int 		 *CPED_S;
int 		  CPED_NS;
void CPED_Print(){
  int s,n;
  if(CPED_NS){
    for(s=0;s<CPED_NS;s++){
      if(CPED_NH[s]){
        fprintf(stdout,"Group:\tCPED\tSector:\t%d \tNhits:\t%d\n",CPED_S[s],CPED_NH[s]);
        for(n=0;n<CPED_NH[s];n++){
          fprintf(stdout,"slot: %d\tchannel: %d\tmean: %d\tsigma: %8.3f\toffset: %d\t",
		CPED[s][n].slot,
		CPED[s][n].channel,
		CPED[s][n].mean,
		CPED[s][n].sigma,
		CPED[s][n].offset);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank CSQL - multi sector bank
struct CSQL_t	**CSQL;
int 		 *CSQL_NH;
int 		 *CSQL_S;
int 		  CSQL_NS;
void CSQL_Print(){
  int s,n;
  if(CSQL_NS){
    for(s=0;s<CSQL_NS;s++){
      if(CSQL_NH[s]){
        fprintf(stdout,"Group:\tCSQL\tSector:\t%d \tNhits:\t%d\n",CSQL_S[s],CSQL_NH[s]);
        for(n=0;n<CSQL_NH[s];n++){
          fprintf(stdout,"EVID: %d\tNPROC: %d\tCPU: %8.3f\tFC: %8.3f\tFCG: %8.3f\tTG: %8.3f\tIBEAM: %8.3f\tNeS1: %d\tNeS2: %d\tNeS3: %d\tNeS4: %d\tNeS5: %d\tNeS6: %d\tNhb: %d\tNtb: %d\tNprot: %d\tNpip: %d\tNdeut: %d\tNphot: %d\tNelhp: %d\tNelhn: %d\t",
		CSQL[s][n].EVID,
		CSQL[s][n].NPROC,
		CSQL[s][n].CPU,
		CSQL[s][n].FC,
		CSQL[s][n].FCG,
		CSQL[s][n].TG,
		CSQL[s][n].IBEAM,
		CSQL[s][n].NeS1,
		CSQL[s][n].NeS2,
		CSQL[s][n].NeS3,
		CSQL[s][n].NeS4,
		CSQL[s][n].NeS5,
		CSQL[s][n].NeS6,
		CSQL[s][n].Nhb,
		CSQL[s][n].Ntb,
		CSQL[s][n].Nprot,
		CSQL[s][n].Npip,
		CSQL[s][n].Ndeut,
		CSQL[s][n].Nphot,
		CSQL[s][n].Nelhp,
		CSQL[s][n].Nelhn);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DC0 - multi sector bank
struct DC0_t	**DC0;
int 		 *DC0_NH;
int 		 *DC0_S;
int 		  DC0_NS;
void DC0_Print(){
  int s,n;
  if(DC0_NS){
    for(s=0;s<DC0_NS;s++){
      if(DC0_NH[s]){
        fprintf(stdout,"Group:\tDC0\tSector:\t%d \tNhits:\t%d\n",DC0_S[s],DC0_NH[s]);
        for(n=0;n<DC0_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		DC0[s][n].ID,
		DC0[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DC1 - multi sector bank
struct DC1_t	**DC1;
int 		 *DC1_NH;
int 		 *DC1_S;
int 		  DC1_NS;
void DC1_Print(){
  int s,n;
  if(DC1_NS){
    for(s=0;s<DC1_NS;s++){
      if(DC1_NH[s]){
        fprintf(stdout,"Group:\tDC1\tSector:\t%d \tNhits:\t%d\n",DC1_S[s],DC1_NH[s]);
        for(n=0;n<DC1_NH[s];n++){
          fprintf(stdout,"ID: %d\ttime: %8.3f\t",
		DC1[s][n].ID,
		DC1[s][n].time);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCDW - multi sector bank
struct DCDW_t	**DCDW;
int 		 *DCDW_NH;
int 		 *DCDW_S;
int 		  DCDW_NS;
void DCDW_Print(){
  int s,n;
  if(DCDW_NS){
    for(s=0;s<DCDW_NS;s++){
      if(DCDW_NH[s]){
        fprintf(stdout,"Group:\tDCDW\tSector:\t%d \tNhits:\t%d\n",DCDW_S[s],DCDW_NH[s]);
        for(n=0;n<DCDW_NH[s];n++){
          fprintf(stdout,"ID: %d\tTIDLY: %8.3f\tSTAT: %d\t",
		DCDW[s][n].ID,
		DCDW[s][n].TIDLY,
		DCDW[s][n].STAT);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCGM - multi sector bank
struct DCGM_t	**DCGM;
int 		 *DCGM_NH;
int 		 *DCGM_S;
int 		  DCGM_NS;
void DCGM_Print(){
  int s,n;
  if(DCGM_NS){
    for(s=0;s<DCGM_NS;s++){
      if(DCGM_NH[s]){
        fprintf(stdout,"Group:\tDCGM\tSector:\t%d \tNhits:\t%d\n",DCGM_S[s],DCGM_NH[s]);
        for(n=0;n<DCGM_NH[s];n++){
          fprintf(stdout,"x_curve: %8.3f\ty_curve: %8.3f\tz_curve: %8.3f\tr_curve: %8.3f\ttheta_start: %8.3f\td_theta: %8.3f\tx_nmid: %8.3f\ty_nmid: %8.3f\tz_nmid: %8.3f\ttheta_min: %8.3f\ttheta_max: %8.3f\tmin_wire: %d\tmax_wire: %d\tstereo: %8.3f\tcell_size: %8.3f\tx_norm: %8.3f\ty_norm: %8.3f\tz_norm: %8.3f\tp_dist: %8.3f\tp_sep: %8.3f\tmax_cylw: %d\t",
		DCGM[s][n].x_curve,
		DCGM[s][n].y_curve,
		DCGM[s][n].z_curve,
		DCGM[s][n].r_curve,
		DCGM[s][n].theta_start,
		DCGM[s][n].d_theta,
		DCGM[s][n].x_nmid,
		DCGM[s][n].y_nmid,
		DCGM[s][n].z_nmid,
		DCGM[s][n].theta_min,
		DCGM[s][n].theta_max,
		DCGM[s][n].min_wire,
		DCGM[s][n].max_wire,
		DCGM[s][n].stereo,
		DCGM[s][n].cell_size,
		DCGM[s][n].x_norm,
		DCGM[s][n].y_norm,
		DCGM[s][n].z_norm,
		DCGM[s][n].p_dist,
		DCGM[s][n].p_sep,
		DCGM[s][n].max_cylw);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCGW - multi sector bank
struct DCGW_t	**DCGW;
int 		 *DCGW_NH;
int 		 *DCGW_S;
int 		  DCGW_NS;
void DCGW_Print(){
  int s,n;
  if(DCGW_NS){
    for(s=0;s<DCGW_NS;s++){
      if(DCGW_NH[s]){
        fprintf(stdout,"Group:\tDCGW\tSector:\t%d \tNhits:\t%d\n",DCGW_S[s],DCGW_NH[s]);
        for(n=0;n<DCGW_NH[s];n++){
          fprintf(stdout,"x_mid: %8.3f\ty_mid: %8.3f\tz_mid: %8.3f\tx_dir: %8.3f\ty_dir: %8.3f\tz_dir: %8.3f\tw_len: %8.3f\tw_len_hv: %8.3f\t",
		DCGW[s][n].x_mid,
		DCGW[s][n].y_mid,
		DCGW[s][n].z_mid,
		DCGW[s][n].x_dir,
		DCGW[s][n].y_dir,
		DCGW[s][n].z_dir,
		DCGW[s][n].w_len,
		DCGW[s][n].w_len_hv);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCH - multi sector bank
struct DCH_t	**DCH;
int 		 *DCH_NH;
int 		 *DCH_S;
int 		  DCH_NS;
void DCH_Print(){
  int s,n;
  if(DCH_NS){
    for(s=0;s<DCH_NS;s++){
      if(DCH_NH[s]){
        fprintf(stdout,"Group:\tDCH\tSector:\t%d \tNhits:\t%d\n",DCH_S[s],DCH_NH[s]);
        for(n=0;n<DCH_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tstep: %8.3f\tdedx: %8.3f\tpmom: %8.3f\ttime: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\ttrack: %d\tid: %d\tlayer: %d\t",
		DCH[s][n].x,
		DCH[s][n].y,
		DCH[s][n].z,
		DCH[s][n].step,
		DCH[s][n].dedx,
		DCH[s][n].pmom,
		DCH[s][n].time,
		DCH[s][n].cx,
		DCH[s][n].cy,
		DCH[s][n].cz,
		DCH[s][n].track,
		DCH[s][n].id,
		DCH[s][n].layer);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCMN - multi sector bank
struct DCMN_t	**DCMN;
int 		 *DCMN_NH;
int 		 *DCMN_S;
int 		  DCMN_NS;
void DCMN_Print(){
  int s,n;
  if(DCMN_NS){
    for(s=0;s<DCMN_NS;s++){
      if(DCMN_NH[s]){
        fprintf(stdout,"Group:\tDCMN\tSector:\t%d \tNhits:\t%d\n",DCMN_S[s],DCMN_NH[s]);
        for(n=0;n<DCMN_NH[s];n++){
          fprintf(stdout,"HBT_evt_1: %8.3f\tTBT_evt_1: %8.3f\tHBT_evt_2: %8.3f\tTBT_evt_2: %8.3f\tHBT_evt_3: %8.3f\tTBT_evt_3: %8.3f\tHBT_evt_4: %8.3f\tTBT_evt_4: %8.3f\tHBT_evt_5: %8.3f\tTBT_evt_5: %8.3f\tHBT_evt_6: %8.3f\tTBT_evt_6: %8.3f\tHBT_evt_all: %8.3f\tTBT_evt_all: %8.3f\tRes_s1_sl1: %8.3f\tRes_s1_sl2: %8.3f\tRes_s1_sl3: %8.3f\tRes_s1_sl4: %8.3f\tRes_s1_sl5: %8.3f\tRes_s1_sl6: %8.3f\tRes_s2_sl1: %8.3f\tRes_s2_sl2: %8.3f\tRes_s2_sl3: %8.3f\tRes_s2_sl4: %8.3f\tRes_s2_sl5: %8.3f\tRes_s2_sl6: %8.3f\tRes_s3_sl1: %8.3f\tRes_s3_sl2: %8.3f\tRes_s3_sl3: %8.3f\tRes_s3_sl4: %8.3f\tRes_s3_sl5: %8.3f\tRes_s3_sl6: %8.3f\tRes_s4_sl1: %8.3f\tRes_s4_sl2: %8.3f\tRes_s4_sl3: %8.3f\tRes_s4_sl4: %8.3f\tRes_s4_sl5: %8.3f\tRes_s4_sl6: %8.3f\tRes_s5_sl1: %8.3f\tRes_s5_sl2: %8.3f\tRes_s5_sl3: %8.3f\tRes_s5_sl4: %8.3f\tRes_s5_sl5: %8.3f\tRes_s5_sl6: %8.3f\tRes_s6_sl1: %8.3f\tRes_s6_sl2: %8.3f\tRes_s6_sl3: %8.3f\tRes_s6_sl4: %8.3f\tRes_s6_sl5: %8.3f\tRes_s6_sl6: %8.3f\t",
		DCMN[s][n].HBT_evt_1,
		DCMN[s][n].TBT_evt_1,
		DCMN[s][n].HBT_evt_2,
		DCMN[s][n].TBT_evt_2,
		DCMN[s][n].HBT_evt_3,
		DCMN[s][n].TBT_evt_3,
		DCMN[s][n].HBT_evt_4,
		DCMN[s][n].TBT_evt_4,
		DCMN[s][n].HBT_evt_5,
		DCMN[s][n].TBT_evt_5,
		DCMN[s][n].HBT_evt_6,
		DCMN[s][n].TBT_evt_6,
		DCMN[s][n].HBT_evt_all,
		DCMN[s][n].TBT_evt_all,
		DCMN[s][n].Res_s1_sl1,
		DCMN[s][n].Res_s1_sl2,
		DCMN[s][n].Res_s1_sl3,
		DCMN[s][n].Res_s1_sl4,
		DCMN[s][n].Res_s1_sl5,
		DCMN[s][n].Res_s1_sl6,
		DCMN[s][n].Res_s2_sl1,
		DCMN[s][n].Res_s2_sl2,
		DCMN[s][n].Res_s2_sl3,
		DCMN[s][n].Res_s2_sl4,
		DCMN[s][n].Res_s2_sl5,
		DCMN[s][n].Res_s2_sl6,
		DCMN[s][n].Res_s3_sl1,
		DCMN[s][n].Res_s3_sl2,
		DCMN[s][n].Res_s3_sl3,
		DCMN[s][n].Res_s3_sl4,
		DCMN[s][n].Res_s3_sl5,
		DCMN[s][n].Res_s3_sl6,
		DCMN[s][n].Res_s4_sl1,
		DCMN[s][n].Res_s4_sl2,
		DCMN[s][n].Res_s4_sl3,
		DCMN[s][n].Res_s4_sl4,
		DCMN[s][n].Res_s4_sl5,
		DCMN[s][n].Res_s4_sl6,
		DCMN[s][n].Res_s5_sl1,
		DCMN[s][n].Res_s5_sl2,
		DCMN[s][n].Res_s5_sl3,
		DCMN[s][n].Res_s5_sl4,
		DCMN[s][n].Res_s5_sl5,
		DCMN[s][n].Res_s5_sl6,
		DCMN[s][n].Res_s6_sl1,
		DCMN[s][n].Res_s6_sl2,
		DCMN[s][n].Res_s6_sl3,
		DCMN[s][n].Res_s6_sl4,
		DCMN[s][n].Res_s6_sl5,
		DCMN[s][n].Res_s6_sl6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCPB - single sector bank
struct DCPB_t	 *DCPB;
int 		  DCPB_NH;
int 		  DCPB_S;
void DCPB_Print(){
  int n;
  if(DCPB_NH){
    fprintf(stdout,"Group:\tDCPB\tSector:\t%d \tNhits:\t%d\n",DCPB_S,DCPB_NH);
    for(n=0;n<DCPB_NH;n++){
      fprintf(stdout,"ScTr: %d\tx_SC: %8.3f\ty_SC: %8.3f\tz_SC: %8.3f\tCX_SC: %8.3f\tCY_SC: %8.3f\tCZ_SC: %8.3f\tX_v: %8.3f\tY_v: %8.3f\tZ_v: %8.3f\tR_v: %8.3f\tChi2: %8.3f\tStatus: %d\t",
		DCPB[n].ScTr,
		DCPB[n].x_SC,
		DCPB[n].y_SC,
		DCPB[n].z_SC,
		DCPB[n].CX_SC,
		DCPB[n].CY_SC,
		DCPB[n].CZ_SC,
		DCPB[n].X_v,
		DCPB[n].Y_v,
		DCPB[n].Z_v,
		DCPB[n].R_v,
		DCPB[n].Chi2,
		DCPB[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank DCV1 - multi sector bank
struct DCV1_t	**DCV1;
int 		 *DCV1_NH;
int 		 *DCV1_S;
int 		  DCV1_NS;
void DCV1_Print(){
  int s,n;
  if(DCV1_NS){
    for(s=0;s<DCV1_NS;s++){
      if(DCV1_NH[s]){
        fprintf(stdout,"Group:\tDCV1\tSector:\t%d \tNhits:\t%d\n",DCV1_S[s],DCV1_NH[s]);
        for(n=0;n<DCV1_NH[s];n++){
          fprintf(stdout,"TsR1: %8.3f\tV0R1: %8.3f\tTmR1: %8.3f\tsp1R1: %8.3f\tsp2R1: %8.3f\t",
		DCV1[s][n].TsR1,
		DCV1[s][n].V0R1,
		DCV1[s][n].TmR1,
		DCV1[s][n].sp1R1,
		DCV1[s][n].sp2R1);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCV2 - multi sector bank
struct DCV2_t	**DCV2;
int 		 *DCV2_NH;
int 		 *DCV2_S;
int 		  DCV2_NS;
void DCV2_Print(){
  int s,n;
  if(DCV2_NS){
    for(s=0;s<DCV2_NS;s++){
      if(DCV2_NH[s]){
        fprintf(stdout,"Group:\tDCV2\tSector:\t%d \tNhits:\t%d\n",DCV2_S[s],DCV2_NH[s]);
        for(n=0;n<DCV2_NH[s];n++){
          fprintf(stdout,"Ts1R2: %8.3f\tV01R2: %8.3f\tVA1R2: %8.3f\tVB1R2: %8.3f\tTm1R2: %8.3f\tTA1R2: %8.3f\tTB1R2: %8.3f\tTs2R2: %8.3f\tV02R2: %8.3f\tVA2R2: %8.3f\tVB2R2: %8.3f\tTm2R2: %8.3f\tTA2R2: %8.3f\tTB2R2: %8.3f\t",
		DCV2[s][n].Ts1R2,
		DCV2[s][n].V01R2,
		DCV2[s][n].VA1R2,
		DCV2[s][n].VB1R2,
		DCV2[s][n].Tm1R2,
		DCV2[s][n].TA1R2,
		DCV2[s][n].TB1R2,
		DCV2[s][n].Ts2R2,
		DCV2[s][n].V02R2,
		DCV2[s][n].VA2R2,
		DCV2[s][n].VB2R2,
		DCV2[s][n].Tm2R2,
		DCV2[s][n].TA2R2,
		DCV2[s][n].TB2R2);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DCV3 - multi sector bank
struct DCV3_t	**DCV3;
int 		 *DCV3_NH;
int 		 *DCV3_S;
int 		  DCV3_NS;
void DCV3_Print(){
  int s,n;
  if(DCV3_NS){
    for(s=0;s<DCV3_NS;s++){
      if(DCV3_NH[s]){
        fprintf(stdout,"Group:\tDCV3\tSector:\t%d \tNhits:\t%d\n",DCV3_S[s],DCV3_NH[s]);
        for(n=0;n<DCV3_NH[s];n++){
          fprintf(stdout,"TsR3: %8.3f\tV0R3: %8.3f\tTmR3: %8.3f\tsp1R3: %8.3f\tsp2R3: %8.3f\t",
		DCV3[s][n].TsR3,
		DCV3[s][n].V0R3,
		DCV3[s][n].TmR3,
		DCV3[s][n].sp1R3,
		DCV3[s][n].sp2R3);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DDLY - multi sector bank
struct DDLY_t	**DDLY;
int 		 *DDLY_NH;
int 		 *DDLY_S;
int 		  DDLY_NS;
void DDLY_Print(){
  int s,n;
  if(DDLY_NS){
    for(s=0;s<DDLY_NS;s++){
      if(DDLY_NH[s]){
        fprintf(stdout,"Group:\tDDLY\tSector:\t%d \tNhits:\t%d\n",DDLY_S[s],DDLY_NH[s]);
        for(n=0;n<DDLY_NH[s];n++){
          fprintf(stdout,"ID: %d\tTIDLY: %8.3f\tSTAT: %d\t",
		DDLY[s][n].ID,
		DDLY[s][n].TIDLY,
		DDLY[s][n].STAT);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DGEO - multi sector bank
struct DGEO_t	**DGEO;
int 		 *DGEO_NH;
int 		 *DGEO_S;
int 		  DGEO_NS;
void DGEO_Print(){
  int s,n;
  if(DGEO_NS){
    for(s=0;s<DGEO_NS;s++){
      if(DGEO_NH[s]){
        fprintf(stdout,"Group:\tDGEO\tSector:\t%d \tNhits:\t%d\n",DGEO_S[s],DGEO_NH[s]);
        for(n=0;n<DGEO_NH[s];n++){
          fprintf(stdout,"ID_sec: %d\tID_reg: %d\txpos: %8.3f\typos: %8.3f\tzpos: %8.3f\tsxpos: %8.3f\tsypos: %8.3f\tszpos: %8.3f\t",
		DGEO[s][n].ID_sec,
		DGEO[s][n].ID_reg,
		DGEO[s][n].xpos,
		DGEO[s][n].ypos,
		DGEO[s][n].zpos,
		DGEO[s][n].sxpos,
		DGEO[s][n].sypos,
		DGEO[s][n].szpos);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DHCL - multi sector bank
struct DHCL_t	**DHCL;
int 		 *DHCL_NH;
int 		 *DHCL_S;
int 		  DHCL_NS;
void DHCL_Print(){
  int s,n;
  if(DHCL_NS){
    for(s=0;s<DHCL_NS;s++){
      if(DHCL_NH[s]){
        fprintf(stdout,"Group:\tDHCL\tSector:\t%d \tNhits:\t%d\n",DHCL_S[s],DHCL_NH[s]);
        for(n=0;n<DHCL_NH[s];n++){
          fprintf(stdout,"SLY: %d\tBTRK: %d\tTRKS1: %d\tTRKS2: %d\tWIRE1: %d\tBWIR1: %d\tWIRE2: %d\tBWIR2: %d\tWIRE3: %d\tBWIR3: %d\tWIRE4: %d\tBWIR4: %d\tWIRE5: %d\tBWIR5: %d\tWIRE6: %d\tBWIR6: %d\t",
		DHCL[s][n].SLY,
		DHCL[s][n].BTRK,
		DHCL[s][n].TRKS1,
		DHCL[s][n].TRKS2,
		DHCL[s][n].WIRE1,
		DHCL[s][n].BWIR1,
		DHCL[s][n].WIRE2,
		DHCL[s][n].BWIR2,
		DHCL[s][n].WIRE3,
		DHCL[s][n].BWIR3,
		DHCL[s][n].WIRE4,
		DHCL[s][n].BWIR4,
		DHCL[s][n].WIRE5,
		DHCL[s][n].BWIR5,
		DHCL[s][n].WIRE6,
		DHCL[s][n].BWIR6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DITM - multi sector bank
struct DITM_t	**DITM;
int 		 *DITM_NH;
int 		 *DITM_S;
int 		  DITM_NS;
void DITM_Print(){
  int s,n;
  if(DITM_NS){
    for(s=0;s<DITM_NS;s++){
      if(DITM_NH[s]){
        fprintf(stdout,"Group:\tDITM\tSector:\t%d \tNhits:\t%d\n",DITM_S[s],DITM_NH[s]);
        for(n=0;n<DITM_NH[s];n++){
          fprintf(stdout,"time: %d\t",
		DITM[s][n].time);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DOCA - multi sector bank
struct DOCA_t	**DOCA;
int 		 *DOCA_NH;
int 		 *DOCA_S;
int 		  DOCA_NS;
void DOCA_Print(){
  int s,n;
  if(DOCA_NS){
    for(s=0;s<DOCA_NS;s++){
      if(DOCA_NH[s]){
        fprintf(stdout,"Group:\tDOCA\tSector:\t%d \tNhits:\t%d\n",DOCA_S[s],DOCA_NH[s]);
        for(n=0;n<DOCA_NH[s];n++){
          fprintf(stdout,"ID: %hd\tDOCA: %hd\t",
		DOCA[s][n].ID,
		DOCA[s][n].DOCA);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DPCP - multi sector bank
struct DPCP_t	**DPCP;
int 		 *DPCP_NH;
int 		 *DPCP_S;
int 		  DPCP_NS;
void DPCP_Print(){
  int s,n;
  if(DPCP_NS){
    for(s=0;s<DPCP_NS;s++){
      if(DPCP_NH[s]){
        fprintf(stdout,"Group:\tDPCP\tSector:\t%d \tNhits:\t%d\n",DPCP_S[s],DPCP_NH[s]);
        for(n=0;n<DPCP_NH[s];n++){
          fprintf(stdout,"ID: %d\tMN_mean: %d\tMN_sigma: %8.3f\tLT_mean: %d\tLT_sigma: %8.3f\tRB_mean: %d\tRB_sigma: %8.3f\tLB_mean: %d\tLB_sigma: %8.3f\tRT_mean: %d\tRT_sigma: %8.3f\tVT_mean: %d\tVT_sigma: %8.3f\t",
		DPCP[s][n].ID,
		DPCP[s][n].MN_mean,
		DPCP[s][n].MN_sigma,
		DPCP[s][n].LT_mean,
		DPCP[s][n].LT_sigma,
		DPCP[s][n].RB_mean,
		DPCP[s][n].RB_sigma,
		DPCP[s][n].LB_mean,
		DPCP[s][n].LB_sigma,
		DPCP[s][n].RT_mean,
		DPCP[s][n].RT_sigma,
		DPCP[s][n].VT_mean,
		DPCP[s][n].VT_sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DPSP - multi sector bank
struct DPSP_t	**DPSP;
int 		 *DPSP_NH;
int 		 *DPSP_S;
int 		  DPSP_NS;
void DPSP_Print(){
  int s,n;
  if(DPSP_NS){
    for(s=0;s<DPSP_NS;s++){
      if(DPSP_NH[s]){
        fprintf(stdout,"Group:\tDPSP\tSector:\t%d \tNhits:\t%d\n",DPSP_S[s],DPSP_NH[s]);
        for(n=0;n<DPSP_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean: %d\tsigma: %8.3f\t",
		DPSP[s][n].ID,
		DPSP[s][n].mean,
		DPSP[s][n].sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DSPC - multi sector bank
struct DSPC_t	**DSPC;
int 		 *DSPC_NH;
int 		 *DSPC_S;
int 		  DSPC_NS;
void DSPC_Print(){
  int s,n;
  if(DSPC_NS){
    for(s=0;s<DSPC_NS;s++){
      if(DSPC_NH[s]){
        fprintf(stdout,"Group:\tDSPC\tSector:\t%d \tNhits:\t%d\n",DSPC_S[s],DSPC_NH[s]);
        for(n=0;n<DSPC_NH[s];n++){
          fprintf(stdout,"PCID: %hd\tTDCPC: %hd\tADCMN: %hd\tADCLT: %hd\tADCRB: %hd\tADCLB: %hd\tADCRT: %hd\tADCVE: %hd\t",
		DSPC[s][n].PCID,
		DSPC[s][n].TDCPC,
		DSPC[s][n].ADCMN,
		DSPC[s][n].ADCLT,
		DSPC[s][n].ADCRB,
		DSPC[s][n].ADCLB,
		DSPC[s][n].ADCRT,
		DSPC[s][n].ADCVE);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DSPS - multi sector bank
struct DSPS_t	**DSPS;
int 		 *DSPS_NH;
int 		 *DSPS_S;
int 		  DSPS_NS;
void DSPS_Print(){
  int s,n;
  if(DSPS_NS){
    for(s=0;s<DSPS_NS;s++){
      if(DSPS_NH[s]){
        fprintf(stdout,"Group:\tDSPS\tSector:\t%d \tNhits:\t%d\n",DSPS_S[s],DSPS_NH[s]);
        for(n=0;n<DSPS_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		DSPS[s][n].ID,
		DSPS[s][n].TDC,
		DSPS[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DSTC - single sector bank
struct DSTC_t	 *DSTC;
int 		  DSTC_NH;
int 		  DSTC_S;
void DSTC_Print(){
  int n;
  if(DSTC_NH){
    fprintf(stdout,"Group:\tDSTC\tSector:\t%d \tNhits:\t%d\n",DSTC_S,DSTC_NH);
    for(n=0;n<DSTC_NH;n++){
      fprintf(stdout,"TACID: %hd\tTDCTAC: %hd\tADCLT: %hd\tADCRT: %hd\tADCLB: %hd\tADCRB: %hd\tADCSUM1: %hd\tADCSUM2: %hd\tADCSUM3: %hd\t",
		DSTC[n].TACID,
		DSTC[n].TDCTAC,
		DSTC[n].ADCLT,
		DSTC[n].ADCRT,
		DSTC[n].ADCLB,
		DSTC[n].ADCRB,
		DSTC[n].ADCSUM1,
		DSTC[n].ADCSUM2,
		DSTC[n].ADCSUM3);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank DTCP - multi sector bank
struct DTCP_t	**DTCP;
int 		 *DTCP_NH;
int 		 *DTCP_S;
int 		  DTCP_NS;
void DTCP_Print(){
  int s,n;
  if(DTCP_NS){
    for(s=0;s<DTCP_NS;s++){
      if(DTCP_NH[s]){
        fprintf(stdout,"Group:\tDTCP\tSector:\t%d \tNhits:\t%d\n",DTCP_S[s],DTCP_NH[s]);
        for(n=0;n<DTCP_NH[s];n++){
          fprintf(stdout,"ID: %d\tLT_mean: %d\tLT_sigma: %8.3f\tRT_mean: %d\tRT_sigma: %8.3f\tLB_mean: %d\tLB_sigma: %8.3f\tRB_mean: %d\tRB_sigma: %8.3f\tsum1_mean: %d\tsum1_sigma: %8.3f\tsum2_mean: %d\tsum2_sigma: %8.3f\tsum3_mean: %d\tsum3_sigma: %8.3f\t",
		DTCP[s][n].ID,
		DTCP[s][n].LT_mean,
		DTCP[s][n].LT_sigma,
		DTCP[s][n].RT_mean,
		DTCP[s][n].RT_sigma,
		DTCP[s][n].LB_mean,
		DTCP[s][n].LB_sigma,
		DTCP[s][n].RB_mean,
		DTCP[s][n].RB_sigma,
		DTCP[s][n].sum1_mean,
		DTCP[s][n].sum1_sigma,
		DTCP[s][n].sum2_mean,
		DTCP[s][n].sum2_sigma,
		DTCP[s][n].sum3_mean,
		DTCP[s][n].sum3_sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank DTRK - multi sector bank
struct DTRK_t	**DTRK;
int 		 *DTRK_NH;
int 		 *DTRK_S;
int 		  DTRK_NS;
void DTRK_Print(){
  int s,n;
  if(DTRK_NS){
    for(s=0;s<DTRK_NS;s++){
      if(DTRK_NH[s]){
        fprintf(stdout,"Group:\tDTRK\tSector:\t%d \tNhits:\t%d\n",DTRK_S[s],DTRK_NH[s]);
        for(n=0;n<DTRK_NH[s];n++){
          fprintf(stdout,"X: %8.3f\tY: %8.3f\tZ: %8.3f\t",
		DTRK[s][n].X,
		DTRK[s][n].Y,
		DTRK[s][n].Z);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EC01 - multi sector bank
struct EC01_t	**EC01;
int 		 *EC01_NH;
int 		 *EC01_S;
int 		  EC01_NS;
void EC01_Print(){
  int s,n;
  if(EC01_NS){
    for(s=0;s<EC01_NS;s++){
      if(EC01_NH[s]){
        fprintf(stdout,"Group:\tEC01\tSector:\t%d \tNhits:\t%d\n",EC01_S[s],EC01_NH[s]);
        for(n=0;n<EC01_NH[s];n++){
          fprintf(stdout,"ID: %d\ttime: %8.3f\tenergy: %8.3f\t",
		EC01[s][n].ID,
		EC01[s][n].time,
		EC01[s][n].energy);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EC1 - multi sector bank
struct EC1_t	**EC1;
int 		 *EC1_NH;
int 		 *EC1_S;
int 		  EC1_NS;
void EC1_Print(){
  int s,n;
  if(EC1_NS){
    for(s=0;s<EC1_NS;s++){
      if(EC1_NH[s]){
        fprintf(stdout,"Group:\tEC1\tSector:\t%d \tNhits:\t%d\n",EC1_S[s],EC1_NH[s]);
        for(n=0;n<EC1_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDCL: %hd\tADCL: %hd\tTDCR: %hd\tADCR: %hd\t",
		EC1[s][n].ID,
		EC1[s][n].TDCL,
		EC1[s][n].ADCL,
		EC1[s][n].TDCR,
		EC1[s][n].ADCR);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EC1P - multi sector bank
struct EC1P_t	**EC1P;
int 		 *EC1P_NH;
int 		 *EC1P_S;
int 		  EC1P_NS;
void EC1P_Print(){
  int s,n;
  if(EC1P_NS){
    for(s=0;s<EC1P_NS;s++){
      if(EC1P_NH[s]){
        fprintf(stdout,"Group:\tEC1P\tSector:\t%d \tNhits:\t%d\n",EC1P_S[s],EC1P_NH[s]);
        for(n=0;n<EC1P_NH[s];n++){
          fprintf(stdout,"N1x: %8.3f\tN1y: %8.3f\tN1z: %8.3f\tR1n: %8.3f\t",
		EC1P[s][n].N1x,
		EC1P[s][n].N1y,
		EC1P[s][n].N1z,
		EC1P[s][n].R1n);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EC1R - multi sector bank
struct EC1R_t	**EC1R;
int 		 *EC1R_NH;
int 		 *EC1R_S;
int 		  EC1R_NS;
void EC1R_Print(){
  int s,n;
  if(EC1R_NS){
    for(s=0;s<EC1R_NS;s++){
      if(EC1R_NH[s]){
        fprintf(stdout,"Group:\tEC1R\tSector:\t%d \tNhits:\t%d\n",EC1R_S[s],EC1R_NH[s]);
        for(n=0;n<EC1R_NH[s];n++){
          fprintf(stdout,"E_tot: %8.3f\tdE_tot: %8.3f\tt_tot: %8.3f\tdt_tot: %8.3f\tx_m: %8.3f\ty_m: %8.3f\tz_m: %8.3f\tdx_m: %8.3f\tdy_m: %8.3f\tdz_m: %8.3f\tE_in: %8.3f\tt_in: %8.3f\tx_in: %8.3f\ty_in: %8.3f\tx_out: %8.3f\ty_out: %8.3f\tx2_in_l: %8.3f\tx2_in_r: %8.3f\ty2_in_l: %8.3f\ty2_in_r: %8.3f\tx2_out_l: %8.3f\tx2_out_r: %8.3f\ty2_out_l: %8.3f\ty2_out_r: %8.3f\ti_in: %d\tj_in: %d\ti_out: %d\tj_out: %d\ta_in_xl: %8.3f\ta_in_xr: %8.3f\ta_in_yl: %8.3f\ta_in_yr: %8.3f\ta_out_xl: %8.3f\ta_out_xr: %8.3f\ta_out_yl: %8.3f\ta_out_yr: %8.3f\tt_in_xs: %8.3f\tt_in_xd: %8.3f\tt_in_ys: %8.3f\tt_in_yd: %8.3f\tt_out_xs: %8.3f\tt_out_xd: %8.3f\tt_out_ys: %8.3f\tt_out_yd: %8.3f\tibl: %d\tncluster: %d\tpmtfired: %d\tz_in: %8.3f\tz_out: %8.3f\tistat: %d\t",
		EC1R[s][n].E_tot,
		EC1R[s][n].dE_tot,
		EC1R[s][n].t_tot,
		EC1R[s][n].dt_tot,
		EC1R[s][n].x_m,
		EC1R[s][n].y_m,
		EC1R[s][n].z_m,
		EC1R[s][n].dx_m,
		EC1R[s][n].dy_m,
		EC1R[s][n].dz_m,
		EC1R[s][n].E_in,
		EC1R[s][n].t_in,
		EC1R[s][n].x_in,
		EC1R[s][n].y_in,
		EC1R[s][n].x_out,
		EC1R[s][n].y_out,
		EC1R[s][n].x2_in_l,
		EC1R[s][n].x2_in_r,
		EC1R[s][n].y2_in_l,
		EC1R[s][n].y2_in_r,
		EC1R[s][n].x2_out_l,
		EC1R[s][n].x2_out_r,
		EC1R[s][n].y2_out_l,
		EC1R[s][n].y2_out_r,
		EC1R[s][n].i_in,
		EC1R[s][n].j_in,
		EC1R[s][n].i_out,
		EC1R[s][n].j_out,
		EC1R[s][n].a_in_xl,
		EC1R[s][n].a_in_xr,
		EC1R[s][n].a_in_yl,
		EC1R[s][n].a_in_yr,
		EC1R[s][n].a_out_xl,
		EC1R[s][n].a_out_xr,
		EC1R[s][n].a_out_yl,
		EC1R[s][n].a_out_yr,
		EC1R[s][n].t_in_xs,
		EC1R[s][n].t_in_xd,
		EC1R[s][n].t_in_ys,
		EC1R[s][n].t_in_yd,
		EC1R[s][n].t_out_xs,
		EC1R[s][n].t_out_xd,
		EC1R[s][n].t_out_ys,
		EC1R[s][n].t_out_yd,
		EC1R[s][n].ibl,
		EC1R[s][n].ncluster,
		EC1R[s][n].pmtfired,
		EC1R[s][n].z_in,
		EC1R[s][n].z_out,
		EC1R[s][n].istat);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECCA - multi sector bank
struct ECCA_t	**ECCA;
int 		 *ECCA_NH;
int 		 *ECCA_S;
int 		  ECCA_NS;
void ECCA_Print(){
  int s,n;
  if(ECCA_NS){
    for(s=0;s<ECCA_NS;s++){
      if(ECCA_NH[s]){
        fprintf(stdout,"Group:\tECCA\tSector:\t%d \tNhits:\t%d\n",ECCA_S[s],ECCA_NH[s]);
        for(n=0;n<ECCA_NH[s];n++){
          fprintf(stdout,"id: %d\taPED: %8.3f\taSIG: %8.3f\taMIP: %8.3f\taMIPu: %8.3f\taSHR: %8.3f\taSHRu: %8.3f\tstat: %d\t",
		ECCA[s][n].id,
		ECCA[s][n].aPED,
		ECCA[s][n].aSIG,
		ECCA[s][n].aMIP,
		ECCA[s][n].aMIPu,
		ECCA[s][n].aSHR,
		ECCA[s][n].aSHRu,
		ECCA[s][n].stat);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECCL - multi sector bank
struct ECCL_t	**ECCL;
int 		 *ECCL_NH;
int 		 *ECCL_S;
int 		  ECCL_NS;
void ECCL_Print(){
  int s,n;
  if(ECCL_NS){
    for(s=0;s<ECCL_NS;s++){
      if(ECCL_NH[s]){
        fprintf(stdout,"Group:\tECCL\tSector:\t%d \tNhits:\t%d\n",ECCL_S[s],ECCL_NH[s]);
        for(n=0;n<ECCL_NH[s];n++){
          fprintf(stdout,"id: %d\tlDB: %8.3f\tlDBu: %8.3f\tlMIP: %8.3f\tlMIPu: %8.3f\tlSHR: %8.3f\tlSHRu: %8.3f\tstat: %d\t",
		ECCL[s][n].id,
		ECCL[s][n].lDB,
		ECCL[s][n].lDBu,
		ECCL[s][n].lMIP,
		ECCL[s][n].lMIPu,
		ECCL[s][n].lSHR,
		ECCL[s][n].lSHRu,
		ECCL[s][n].stat);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECCT - multi sector bank
struct ECCT_t	**ECCT;
int 		 *ECCT_NH;
int 		 *ECCT_S;
int 		  ECCT_NS;
void ECCT_Print(){
  int s,n;
  if(ECCT_NS){
    for(s=0;s<ECCT_NS;s++){
      if(ECCT_NH[s]){
        fprintf(stdout,"Group:\tECCT\tSector:\t%d \tNhits:\t%d\n",ECCT_S[s],ECCT_NH[s]);
        for(n=0;n<ECCT_NH[s];n++){
          fprintf(stdout,"id: %d\ttOFF: %8.3f\ttOFFu: %8.3f\ttGAIN: %8.3f\ttGAINu: %8.3f\ttW0: %8.3f\ttW0u: %8.3f\ttW1: %8.3f\ttW1u: %8.3f\ttVEF: %8.3f\ttVEFu: %8.3f\tstat: %d\t",
		ECCT[s][n].id,
		ECCT[s][n].tOFF,
		ECCT[s][n].tOFFu,
		ECCT[s][n].tGAIN,
		ECCT[s][n].tGAINu,
		ECCT[s][n].tW0,
		ECCT[s][n].tW0u,
		ECCT[s][n].tW1,
		ECCT[s][n].tW1u,
		ECCT[s][n].tVEF,
		ECCT[s][n].tVEFu,
		ECCT[s][n].stat);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EC - multi sector bank
struct EC_t	**EC;
int 		 *EC_NH;
int 		 *EC_S;
int 		  EC_NS;
void EC_Print(){
  int s,n;
  if(EC_NS){
    for(s=0;s<EC_NS;s++){
      if(EC_NH[s]){
        fprintf(stdout,"Group:\tEC\tSector:\t%d \tNhits:\t%d\n",EC_S[s],EC_NH[s]);
        for(n=0;n<EC_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		EC[s][n].ID,
		EC[s][n].TDC,
		EC[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECDI - multi sector bank
struct ECDI_t	**ECDI;
int 		 *ECDI_NH;
int 		 *ECDI_S;
int 		  ECDI_NS;
void ECDI_Print(){
  int s,n;
  if(ECDI_NS){
    for(s=0;s<ECDI_NS;s++){
      if(ECDI_NH[s]){
        fprintf(stdout,"Group:\tECDI\tSector:\t%d \tNhits:\t%d\n",ECDI_S[s],ECDI_NH[s]);
        for(n=0;n<ECDI_NH[s];n++){
          fprintf(stdout,"crate: %d\tslot: %d\tmask: %d\tthreshold: %d\twidth: %d\t",
		ECDI[s][n].crate,
		ECDI[s][n].slot,
		ECDI[s][n].mask,
		ECDI[s][n].threshold,
		ECDI[s][n].width);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECG - multi sector bank
struct ECG_t	**ECG;
int 		 *ECG_NH;
int 		 *ECG_S;
int 		  ECG_NS;
void ECG_Print(){
  int s,n;
  if(ECG_NS){
    for(s=0;s<ECG_NS;s++){
      if(ECG_NH[s]){
        fprintf(stdout,"Group:\tECG\tSector:\t%d \tNhits:\t%d\n",ECG_S[s],ECG_NH[s]);
        for(n=0;n<ECG_NH[s];n++){
          fprintf(stdout,"L0: %8.3f\tTHE0: %8.3f\tYLOW: %8.3f\tYHI: %8.3f\tDYLOW: %8.3f\tDYHI: %8.3f\tLRTH: %8.3f\tTANGR: %8.3f\tSECTOR: %d\tPHISEC: %8.3f\tX0OFF: %8.3f\tY0OFF: %8.3f\tZ0OFF: %8.3f\tROTM11: %8.3f\tROTM12: %8.3f\tROTM13: %8.3f\tROTM21: %8.3f\tROTM22: %8.3f\tROTM23: %8.3f\tROTM31: %8.3f\tROTM32: %8.3f\tROTM33: %8.3f\t",
		ECG[s][n].L0,
		ECG[s][n].THE0,
		ECG[s][n].YLOW,
		ECG[s][n].YHI,
		ECG[s][n].DYLOW,
		ECG[s][n].DYHI,
		ECG[s][n].LRTH,
		ECG[s][n].TANGR,
		ECG[s][n].SECTOR,
		ECG[s][n].PHISEC,
		ECG[s][n].X0OFF,
		ECG[s][n].Y0OFF,
		ECG[s][n].Z0OFF,
		ECG[s][n].ROTM11,
		ECG[s][n].ROTM12,
		ECG[s][n].ROTM13,
		ECG[s][n].ROTM21,
		ECG[s][n].ROTM22,
		ECG[s][n].ROTM23,
		ECG[s][n].ROTM31,
		ECG[s][n].ROTM32,
		ECG[s][n].ROTM33);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECHB - single sector bank
struct ECHB_t	 *ECHB;
int 		  ECHB_NH;
int 		  ECHB_S;
void ECHB_Print(){
  int n;
  if(ECHB_NH){
    fprintf(stdout,"Group:\tECHB\tSector:\t%d \tNhits:\t%d\n",ECHB_S,ECHB_NH);
    for(n=0;n<ECHB_NH;n++){
      fprintf(stdout,"Sect: %d\tE__hit: %8.3f\tdE_hit: %8.3f\tt_hit: %8.3f\tdt_hi: %8.3f\ti_hit: %8.3f\tj_hit: %8.3f\tdi_hit: %8.3f\tdj_hit: %8.3f\tx_hit: %8.3f\ty_hit: %8.3f\tz_hit: %8.3f\tdx_hit: %8.3f\tdy_hit: %8.3f\tdz_hit: %8.3f\tu2_hit: %8.3f\tv2_hit: %8.3f\tw2_hit: %8.3f\tu3_hit: %8.3f\tv3_hit: %8.3f\tw3_hit: %8.3f\tu4_hit: %8.3f\tv4_hit: %8.3f\tw4_hit: %8.3f\tcentr_U: %8.3f\tcentr_V: %8.3f\tcentr_W: %8.3f\tpath_U: %8.3f\tpath_V: %8.3f\tpath_W: %8.3f\tNstrp_U: %d\tNstrp_V: %d\tNstrp_W: %d\tMatchID1: %d\tCH21: %8.3f\tMatchID2: %d\tCH22: %8.3f\tistat: %d\t",
		ECHB[n].Sect,
		ECHB[n].E__hit,
		ECHB[n].dE_hit,
		ECHB[n].t_hit,
		ECHB[n].dt_hi,
		ECHB[n].i_hit,
		ECHB[n].j_hit,
		ECHB[n].di_hit,
		ECHB[n].dj_hit,
		ECHB[n].x_hit,
		ECHB[n].y_hit,
		ECHB[n].z_hit,
		ECHB[n].dx_hit,
		ECHB[n].dy_hit,
		ECHB[n].dz_hit,
		ECHB[n].u2_hit,
		ECHB[n].v2_hit,
		ECHB[n].w2_hit,
		ECHB[n].u3_hit,
		ECHB[n].v3_hit,
		ECHB[n].w3_hit,
		ECHB[n].u4_hit,
		ECHB[n].v4_hit,
		ECHB[n].w4_hit,
		ECHB[n].centr_U,
		ECHB[n].centr_V,
		ECHB[n].centr_W,
		ECHB[n].path_U,
		ECHB[n].path_V,
		ECHB[n].path_W,
		ECHB[n].Nstrp_U,
		ECHB[n].Nstrp_V,
		ECHB[n].Nstrp_W,
		ECHB[n].MatchID1,
		ECHB[n].CH21,
		ECHB[n].MatchID2,
		ECHB[n].CH22,
		ECHB[n].istat);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank ECH - multi sector bank
struct ECH_t	**ECH;
int 		 *ECH_NH;
int 		 *ECH_S;
int 		  ECH_NS;
void ECH_Print(){
  int s,n;
  if(ECH_NS){
    for(s=0;s<ECH_NS;s++){
      if(ECH_NH[s]){
        fprintf(stdout,"Group:\tECH\tSector:\t%d \tNhits:\t%d\n",ECH_S[s],ECH_NH[s]);
        for(n=0;n<ECH_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\tpmom: %8.3f\tid: %d\ttof: %8.3f\tedepin: %8.3f\tedepout: %8.3f\t",
		ECH[s][n].x,
		ECH[s][n].y,
		ECH[s][n].z,
		ECH[s][n].cx,
		ECH[s][n].cy,
		ECH[s][n].cz,
		ECH[s][n].pmom,
		ECH[s][n].id,
		ECH[s][n].tof,
		ECH[s][n].edepin,
		ECH[s][n].edepout);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECMT - multi sector bank
struct ECMT_t	**ECMT;
int 		 *ECMT_NH;
int 		 *ECMT_S;
int 		  ECMT_NS;
void ECMT_Print(){
  int s,n;
  if(ECMT_NS){
    for(s=0;s<ECMT_NS;s++){
      if(ECMT_NH[s]){
        fprintf(stdout,"Group:\tECMT\tSector:\t%d \tNhits:\t%d\n",ECMT_S[s],ECMT_NH[s]);
        for(n=0;n<ECMT_NH[s];n++){
          fprintf(stdout,"in_high: %d\tout_high: %d\ttotal_high: %d\tin_lo: %d\tout_lo: %d\ttotal_lo: %d\t",
		ECMT[s][n].in_high,
		ECMT[s][n].out_high,
		ECMT[s][n].total_high,
		ECMT[s][n].in_lo,
		ECMT[s][n].out_lo,
		ECMT[s][n].total_lo);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECP1 - multi sector bank
struct ECP1_t	**ECP1;
int 		 *ECP1_NH;
int 		 *ECP1_S;
int 		  ECP1_NS;
void ECP1_Print(){
  int s,n;
  if(ECP1_NS){
    for(s=0;s<ECP1_NS;s++){
      if(ECP1_NH[s]){
        fprintf(stdout,"Group:\tECP1\tSector:\t%d \tNhits:\t%d\n",ECP1_S[s],ECP1_NH[s]);
        for(n=0;n<ECP1_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean_left: %d\tsigma_left: %8.3f\tmean_right: %d\tsigma_right: %8.3f\t",
		ECP1[s][n].ID,
		ECP1[s][n].mean_left,
		ECP1[s][n].sigma_left,
		ECP1[s][n].mean_right,
		ECP1[s][n].sigma_right);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECPB - single sector bank
struct ECPB_t	 *ECPB;
int 		  ECPB_NH;
int 		  ECPB_S;
void ECPB_Print(){
  int n;
  if(ECPB_NH){
    fprintf(stdout,"Group:\tECPB\tSector:\t%d \tNhits:\t%d\n",ECPB_S,ECPB_NH);
    for(n=0;n<ECPB_NH;n++){
      fprintf(stdout,"ScHt: %d\tEtot: %8.3f\tEin: %8.3f\tEout: %8.3f\tTime: %8.3f\tPath: %8.3f\tX: %8.3f\tY: %8.3f\tZ: %8.3f\tM2_hit: %8.3f\tM3_hit: %8.3f\tM4_hit: %8.3f\tInnStr: %d\tOutStr: %d\tChi2EC: %8.3f\tStatus: %d\t",
		ECPB[n].ScHt,
		ECPB[n].Etot,
		ECPB[n].Ein,
		ECPB[n].Eout,
		ECPB[n].Time,
		ECPB[n].Path,
		ECPB[n].X,
		ECPB[n].Y,
		ECPB[n].Z,
		ECPB[n].M2_hit,
		ECPB[n].M3_hit,
		ECPB[n].M4_hit,
		ECPB[n].InnStr,
		ECPB[n].OutStr,
		ECPB[n].Chi2EC,
		ECPB[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank ECPC - multi sector bank
struct ECPC_t	**ECPC;
int 		 *ECPC_NH;
int 		 *ECPC_S;
int 		  ECPC_NS;
void ECPC_Print(){
  int s,n;
  if(ECPC_NS){
    for(s=0;s<ECPC_NS;s++){
      if(ECPC_NH[s]){
        fprintf(stdout,"Group:\tECPC\tSector:\t%d \tNhits:\t%d\n",ECPC_S[s],ECPC_NH[s]);
        for(n=0;n<ECPC_NH[s];n++){
          fprintf(stdout,"ID: %d\tTDC: %8.3f\tADC: %8.3f\t",
		ECPC[s][n].ID,
		ECPC[s][n].TDC,
		ECPC[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECP - multi sector bank
struct ECP_t	**ECP;
int 		 *ECP_NH;
int 		 *ECP_S;
int 		  ECP_NS;
void ECP_Print(){
  int s,n;
  if(ECP_NS){
    for(s=0;s<ECP_NS;s++){
      if(ECP_NH[s]){
        fprintf(stdout,"Group:\tECP\tSector:\t%d \tNhits:\t%d\n",ECP_S[s],ECP_NH[s]);
        for(n=0;n<ECP_NH[s];n++){
          fprintf(stdout,"N1x: %8.3f\tN1y: %8.3f\tN1z: %8.3f\tR1n: %8.3f\tPlW: %8.3f\tPlI: %8.3f\tPlO: %8.3f\t",
		ECP[s][n].N1x,
		ECP[s][n].N1y,
		ECP[s][n].N1z,
		ECP[s][n].R1n,
		ECP[s][n].PlW,
		ECP[s][n].PlI,
		ECP[s][n].PlO);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECPE - multi sector bank
struct ECPE_t	**ECPE;
int 		 *ECPE_NH;
int 		 *ECPE_S;
int 		  ECPE_NS;
void ECPE_Print(){
  int s,n;
  if(ECPE_NS){
    for(s=0;s<ECPE_NS;s++){
      if(ECPE_NH[s]){
        fprintf(stdout,"Group:\tECPE\tSector:\t%d \tNhits:\t%d\n",ECPE_S[s],ECPE_NH[s]);
        for(n=0;n<ECPE_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean: %d\tsigma: %8.3f\t",
		ECPE[s][n].ID,
		ECPE[s][n].mean,
		ECPE[s][n].sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECPI - single sector bank
struct ECPI_t	 *ECPI;
int 		  ECPI_NH;
int 		  ECPI_S;
void ECPI_Print(){
  int n;
  if(ECPI_NH){
    fprintf(stdout,"Group:\tECPI\tSector:\t%d \tNhits:\t%d\n",ECPI_S,ECPI_NH);
    for(n=0;n<ECPI_NH;n++){
      fprintf(stdout,"ID: %d\tLayer: %d\tHITID: %d\tiloc: %8.3f\tjloc: %8.3f\tdiloc: %8.3f\tdjloc: %8.3f\tR: %8.3f\tE: %8.3f\t",
		ECPI[n].ID,
		ECPI[n].Layer,
		ECPI[n].HITID,
		ECPI[n].iloc,
		ECPI[n].jloc,
		ECPI[n].diloc,
		ECPI[n].djloc,
		ECPI[n].R,
		ECPI[n].E);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank ECPO - multi sector bank
struct ECPO_t	**ECPO;
int 		 *ECPO_NH;
int 		 *ECPO_S;
int 		  ECPO_NS;
void ECPO_Print(){
  int s,n;
  if(ECPO_NS){
    for(s=0;s<ECPO_NS;s++){
      if(ECPO_NH[s]){
        fprintf(stdout,"Group:\tECPO\tSector:\t%d \tNhits:\t%d\n",ECPO_S[s],ECPO_NH[s]);
        for(n=0;n<ECPO_NH[s];n++){
          fprintf(stdout,"UID1: %d\tUID2: %d\tVID1: %d\tVID2: %d\tWID1: %d\tWID2: %d\tSLH: %d\t",
		ECPO[s][n].UID1,
		ECPO[s][n].UID2,
		ECPO[s][n].VID1,
		ECPO[s][n].VID2,
		ECPO[s][n].WID1,
		ECPO[s][n].WID2,
		ECPO[s][n].SLH);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECRB - multi sector bank
struct ECRB_t	**ECRB;
int 		 *ECRB_NH;
int 		 *ECRB_S;
int 		  ECRB_NS;
void ECRB_Print(){
  int s,n;
  if(ECRB_NS){
    for(s=0;s<ECRB_NS;s++){
      if(ECRB_NH[s]){
        fprintf(stdout,"Group:\tECRB\tSector:\t%d \tNhits:\t%d\n",ECRB_S[s],ECRB_NH[s]);
        for(n=0;n<ECRB_NH[s];n++){
          fprintf(stdout,"E_in: %8.3f\tE_out: %8.3f\tdE_in: %8.3f\tdE_out: %8.3f\tt_in: %8.3f\tt_out: %8.3f\tdt_in: %8.3f\tdt_out: %8.3f\ti_in: %8.3f\tj_in: %8.3f\ti_out: %8.3f\tj_out: %8.3f\tdi_in: %8.3f\tdj_in: %8.3f\tdi_out: %8.3f\tdj_out: %8.3f\tx_in: %8.3f\ty_in: %8.3f\tz_in: %8.3f\tx_out: %8.3f\ty_out: %8.3f\tz_out: %8.3f\tdx_in: %8.3f\tdy_in: %8.3f\tdz_in: %8.3f\tdx_out: %8.3f\tdy_out: %8.3f\tdz_out: %8.3f\tu2_in: %8.3f\tv2_in: %8.3f\tw2_in: %8.3f\tu2_out: %8.3f\tv2_out: %8.3f\tw2_out: %8.3f\tu3_in: %8.3f\tv3_in: %8.3f\tw3_in: %8.3f\tu3_out: %8.3f\tv3_out: %8.3f\tw3_out: %8.3f\ti2: %8.3f\tj2: %8.3f\ti3: %8.3f\tj3: %8.3f\tspare1: %8.3f\tspare2: %8.3f\tspare3: %8.3f\tspare4: %8.3f\tistat: %d\t",
		ECRB[s][n].E_in,
		ECRB[s][n].E_out,
		ECRB[s][n].dE_in,
		ECRB[s][n].dE_out,
		ECRB[s][n].t_in,
		ECRB[s][n].t_out,
		ECRB[s][n].dt_in,
		ECRB[s][n].dt_out,
		ECRB[s][n].i_in,
		ECRB[s][n].j_in,
		ECRB[s][n].i_out,
		ECRB[s][n].j_out,
		ECRB[s][n].di_in,
		ECRB[s][n].dj_in,
		ECRB[s][n].di_out,
		ECRB[s][n].dj_out,
		ECRB[s][n].x_in,
		ECRB[s][n].y_in,
		ECRB[s][n].z_in,
		ECRB[s][n].x_out,
		ECRB[s][n].y_out,
		ECRB[s][n].z_out,
		ECRB[s][n].dx_in,
		ECRB[s][n].dy_in,
		ECRB[s][n].dz_in,
		ECRB[s][n].dx_out,
		ECRB[s][n].dy_out,
		ECRB[s][n].dz_out,
		ECRB[s][n].u2_in,
		ECRB[s][n].v2_in,
		ECRB[s][n].w2_in,
		ECRB[s][n].u2_out,
		ECRB[s][n].v2_out,
		ECRB[s][n].w2_out,
		ECRB[s][n].u3_in,
		ECRB[s][n].v3_in,
		ECRB[s][n].w3_in,
		ECRB[s][n].u3_out,
		ECRB[s][n].v3_out,
		ECRB[s][n].w3_out,
		ECRB[s][n].i2,
		ECRB[s][n].j2,
		ECRB[s][n].i3,
		ECRB[s][n].j3,
		ECRB[s][n].spare1,
		ECRB[s][n].spare2,
		ECRB[s][n].spare3,
		ECRB[s][n].spare4,
		ECRB[s][n].istat);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECS - multi sector bank
struct ECS_t	**ECS;
int 		 *ECS_NH;
int 		 *ECS_S;
int 		  ECS_NS;
void ECS_Print(){
  int s,n;
  if(ECS_NS){
    for(s=0;s<ECS_NS;s++){
      if(ECS_NH[s]){
        fprintf(stdout,"Group:\tECS\tSector:\t%d \tNhits:\t%d\n",ECS_S[s],ECS_NH[s]);
        for(n=0;n<ECS_NH[s];n++){
          fprintf(stdout,"ECS1: %d\tECS2: %d\tECS3: %d\tECS4: %d\tECS5: %d\tECS6: %d\tECS7: %d\tECS8: %d\tECS9: %d\tECS10: %d\tECS11: %d\tECS12: %d\tECS13: %d\tECS14: %d\tECS15: %d\tECS16: %d\tECS17: %d\tECS18: %d\tECS19: %d\tECS20: %d\tECS21: %d\tECS22: %d\tECS23: %d\tECS24: %d\tECS25: %d\tECS26: %d\tECS27: %d\tECS28: %d\tECS29: %d\tECS30: %d\tECS31: %d\tECS32: %d\tECS33: %d\tECS34: %d\tECS35: %d\tECS36: %d\tECS37: %d\tECS38: %d\tECS39: %d\tECS40: %d\tECS41: %d\tECS42: %d\tECS43: %d\tECS44: %d\tECS45: %d\tECS46: %d\tECS47: %d\tECS48: %d\tECS49: %d\tECS50: %d\tECS51: %d\tECS52: %d\tECS53: %d\tECS54: %d\tECS55: %d\tECS56: %d\tECS57: %d\tECS58: %d\tECS59: %d\tECS60: %d\tECS61: %d\tECS62: %d\tECS63: %d\tECS64: %d\tECS65: %d\tECS66: %d\tECS67: %d\tECS68: %d\tECS69: %d\tECS70: %d\tECS71: %d\tECS72: %d\tECS73: %d\tECS74: %d\tECS75: %d\tECS76: %d\tECS77: %d\tECS78: %d\tECS79: %d\tECS80: %d\tECS81: %d\tECS82: %d\tECS83: %d\tECS84: %d\tECS85: %d\tECS86: %d\tECS87: %d\tECS88: %d\tECS89: %d\tECS90: %d\tECS91: %d\tECS92: %d\tECS93: %d\tECS94: %d\tECS95: %d\tECS96: %d\t",
		ECS[s][n].ECS1,
		ECS[s][n].ECS2,
		ECS[s][n].ECS3,
		ECS[s][n].ECS4,
		ECS[s][n].ECS5,
		ECS[s][n].ECS6,
		ECS[s][n].ECS7,
		ECS[s][n].ECS8,
		ECS[s][n].ECS9,
		ECS[s][n].ECS10,
		ECS[s][n].ECS11,
		ECS[s][n].ECS12,
		ECS[s][n].ECS13,
		ECS[s][n].ECS14,
		ECS[s][n].ECS15,
		ECS[s][n].ECS16,
		ECS[s][n].ECS17,
		ECS[s][n].ECS18,
		ECS[s][n].ECS19,
		ECS[s][n].ECS20,
		ECS[s][n].ECS21,
		ECS[s][n].ECS22,
		ECS[s][n].ECS23,
		ECS[s][n].ECS24,
		ECS[s][n].ECS25,
		ECS[s][n].ECS26,
		ECS[s][n].ECS27,
		ECS[s][n].ECS28,
		ECS[s][n].ECS29,
		ECS[s][n].ECS30,
		ECS[s][n].ECS31,
		ECS[s][n].ECS32,
		ECS[s][n].ECS33,
		ECS[s][n].ECS34,
		ECS[s][n].ECS35,
		ECS[s][n].ECS36,
		ECS[s][n].ECS37,
		ECS[s][n].ECS38,
		ECS[s][n].ECS39,
		ECS[s][n].ECS40,
		ECS[s][n].ECS41,
		ECS[s][n].ECS42,
		ECS[s][n].ECS43,
		ECS[s][n].ECS44,
		ECS[s][n].ECS45,
		ECS[s][n].ECS46,
		ECS[s][n].ECS47,
		ECS[s][n].ECS48,
		ECS[s][n].ECS49,
		ECS[s][n].ECS50,
		ECS[s][n].ECS51,
		ECS[s][n].ECS52,
		ECS[s][n].ECS53,
		ECS[s][n].ECS54,
		ECS[s][n].ECS55,
		ECS[s][n].ECS56,
		ECS[s][n].ECS57,
		ECS[s][n].ECS58,
		ECS[s][n].ECS59,
		ECS[s][n].ECS60,
		ECS[s][n].ECS61,
		ECS[s][n].ECS62,
		ECS[s][n].ECS63,
		ECS[s][n].ECS64,
		ECS[s][n].ECS65,
		ECS[s][n].ECS66,
		ECS[s][n].ECS67,
		ECS[s][n].ECS68,
		ECS[s][n].ECS69,
		ECS[s][n].ECS70,
		ECS[s][n].ECS71,
		ECS[s][n].ECS72,
		ECS[s][n].ECS73,
		ECS[s][n].ECS74,
		ECS[s][n].ECS75,
		ECS[s][n].ECS76,
		ECS[s][n].ECS77,
		ECS[s][n].ECS78,
		ECS[s][n].ECS79,
		ECS[s][n].ECS80,
		ECS[s][n].ECS81,
		ECS[s][n].ECS82,
		ECS[s][n].ECS83,
		ECS[s][n].ECS84,
		ECS[s][n].ECS85,
		ECS[s][n].ECS86,
		ECS[s][n].ECS87,
		ECS[s][n].ECS88,
		ECS[s][n].ECS89,
		ECS[s][n].ECS90,
		ECS[s][n].ECS91,
		ECS[s][n].ECS92,
		ECS[s][n].ECS93,
		ECS[s][n].ECS94,
		ECS[s][n].ECS95,
		ECS[s][n].ECS96);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ECT - multi sector bank
struct ECT_t	**ECT;
int 		 *ECT_NH;
int 		 *ECT_S;
int 		  ECT_NS;
void ECT_Print(){
  int s,n;
  if(ECT_NS){
    for(s=0;s<ECT_NS;s++){
      if(ECT_NH[s]){
        fprintf(stdout,"Group:\tECT\tSector:\t%d \tNhits:\t%d\n",ECT_S[s],ECT_NH[s]);
        for(n=0;n<ECT_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		ECT[s][n].ID,
		ECT[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EID0 - multi sector bank
struct EID0_t	**EID0;
int 		 *EID0_NH;
int 		 *EID0_S;
int 		  EID0_NS;
void EID0_Print(){
  int s,n;
  if(EID0_NS){
    for(s=0;s<EID0_NS;s++){
      if(EID0_NH[s]){
        fprintf(stdout,"Group:\tEID0\tSector:\t%d \tNhits:\t%d\n",EID0_S[s],EID0_NH[s]);
        for(n=0;n<EID0_NH[s];n++){
          fprintf(stdout,"jeid0sec: %d\tjeid0cc: %d\tjeid0ec: %d\tjeid0sc: %d\t",
		EID0[s][n].jeid0sec,
		EID0[s][n].jeid0cc,
		EID0[s][n].jeid0ec,
		EID0[s][n].jeid0sc);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EPIC - multi sector bank
struct EPIC_t	**EPIC;
int 		 *EPIC_NH;
int 		 *EPIC_S;
int 		  EPIC_NS;
void EPIC_Print(){
  int s,n;
  if(EPIC_NS){
    for(s=0;s<EPIC_NS;s++){
      if(EPIC_NH[s]){
        fprintf(stdout,"Group:\tEPIC\tSector:\t%d \tNhits:\t%d\n",EPIC_S[s],EPIC_NH[s]);
        for(n=0;n<EPIC_NH[s];n++){
          fprintf(stdout,"value: %8.3f\tchar1: %.32s\t",
		EPIC[s][n].value,
		EPIC[s][n].char1);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank EVNT - single sector bank
struct EVNT_t	 *EVNT;
int 		  EVNT_NH;
int 		  EVNT_S;
void EVNT_Print(){
  int n;
  if(EVNT_NH){
    fprintf(stdout,"Group:\tEVNT\tSector:\t%d \tNhits:\t%d\n",EVNT_S,EVNT_NH);
    for(n=0;n<EVNT_NH;n++){
      fprintf(stdout,"ID: %d\tPmom: %8.3f\tMass: %8.3f\tCharge: %d\tBeta: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\tX: %8.3f\tY: %8.3f\tZ: %8.3f\tDCstat: %d\tCCstat: %d\tSCstat: %d\tECstat: %d\tLCstat: %d\tSTstat: %d\tStatus: %d\t",
		EVNT[n].ID,
		EVNT[n].Pmom,
		EVNT[n].Mass,
		EVNT[n].Charge,
		EVNT[n].Beta,
		EVNT[n].cx,
		EVNT[n].cy,
		EVNT[n].cz,
		EVNT[n].X,
		EVNT[n].Y,
		EVNT[n].Z,
		EVNT[n].DCstat,
		EVNT[n].CCstat,
		EVNT[n].SCstat,
		EVNT[n].ECstat,
		EVNT[n].LCstat,
		EVNT[n].STstat,
		EVNT[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank FBPM - multi sector bank
struct FBPM_t	**FBPM;
int 		 *FBPM_NH;
int 		 *FBPM_S;
int 		  FBPM_NS;
void FBPM_Print(){
  int s,n;
  if(FBPM_NS){
    for(s=0;s<FBPM_NS;s++){
      if(FBPM_NH[s]){
        fprintf(stdout,"Group:\tFBPM\tSector:\t%d \tNhits:\t%d\n",FBPM_S[s],FBPM_NH[s]);
        for(n=0;n<FBPM_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		FBPM[s][n].ID,
		FBPM[s][n].TDC,
		FBPM[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank G1SL - multi sector bank
struct G1SL_t	**G1SL;
int 		 *G1SL_NH;
int 		 *G1SL_S;
int 		  G1SL_NS;
void G1SL_Print(){
  int s,n;
  if(G1SL_NS){
    for(s=0;s<G1SL_NS;s++){
      if(G1SL_NH[s]){
        fprintf(stdout,"Group:\tG1SL\tSector:\t%d \tNhits:\t%d\n",G1SL_S[s],G1SL_NH[s]);
        for(n=0;n<G1SL_NH[s];n++){
          fprintf(stdout,"G1T1: %d\tG1T2: %d\tG1T3: %d\tG1T4: %d\tG1T5: %d\tG1T6: %d\tG1T7: %d\tG1T8: %d\tG1T9: %d\tG1T10: %d\tG1T11: %d\tG1T12: %d\tG1T13: %d\tG1T14: %d\tG1T15: %d\tG1T16: %d\tG1T17: %d\tG1T18: %d\tG1T19: %d\tG1T20: %d\tG1T21: %d\tG1T22: %d\tG1T23: %d\tG1T24: %d\tG1T25: %d\tG1T26: %d\tG1T27: %d\tG1T28: %d\tG1T29: %d\tG1T30: %d\tG1T31: %d\tG1T32: %d\tG1T33: %d\tG1T34: %d\tG1T35: %d\tG1T36: %d\tG1T37: %d\tG1T38: %d\tG1T39: %d\tG1T40: %d\tG1T41: %d\tG1T42: %d\tG1T43: %d\tG1T44: %d\tG1T45: %d\tG1T46: %d\tG1T47: %d\tG1T48: %d\tG1T49: %d\tG1T50: %d\tG1T51: %d\tG1T52: %d\tG1T53: %d\tG1T54: %d\tG1T55: %d\tG1T56: %d\tG1T57: %d\tG1T58: %d\tG1T59: %d\tG1T60: %d\tG1T61: %d\tG1T62: %d\tG1T63: %d\tG1T64: %d\t",
		G1SL[s][n].G1T1,
		G1SL[s][n].G1T2,
		G1SL[s][n].G1T3,
		G1SL[s][n].G1T4,
		G1SL[s][n].G1T5,
		G1SL[s][n].G1T6,
		G1SL[s][n].G1T7,
		G1SL[s][n].G1T8,
		G1SL[s][n].G1T9,
		G1SL[s][n].G1T10,
		G1SL[s][n].G1T11,
		G1SL[s][n].G1T12,
		G1SL[s][n].G1T13,
		G1SL[s][n].G1T14,
		G1SL[s][n].G1T15,
		G1SL[s][n].G1T16,
		G1SL[s][n].G1T17,
		G1SL[s][n].G1T18,
		G1SL[s][n].G1T19,
		G1SL[s][n].G1T20,
		G1SL[s][n].G1T21,
		G1SL[s][n].G1T22,
		G1SL[s][n].G1T23,
		G1SL[s][n].G1T24,
		G1SL[s][n].G1T25,
		G1SL[s][n].G1T26,
		G1SL[s][n].G1T27,
		G1SL[s][n].G1T28,
		G1SL[s][n].G1T29,
		G1SL[s][n].G1T30,
		G1SL[s][n].G1T31,
		G1SL[s][n].G1T32,
		G1SL[s][n].G1T33,
		G1SL[s][n].G1T34,
		G1SL[s][n].G1T35,
		G1SL[s][n].G1T36,
		G1SL[s][n].G1T37,
		G1SL[s][n].G1T38,
		G1SL[s][n].G1T39,
		G1SL[s][n].G1T40,
		G1SL[s][n].G1T41,
		G1SL[s][n].G1T42,
		G1SL[s][n].G1T43,
		G1SL[s][n].G1T44,
		G1SL[s][n].G1T45,
		G1SL[s][n].G1T46,
		G1SL[s][n].G1T47,
		G1SL[s][n].G1T48,
		G1SL[s][n].G1T49,
		G1SL[s][n].G1T50,
		G1SL[s][n].G1T51,
		G1SL[s][n].G1T52,
		G1SL[s][n].G1T53,
		G1SL[s][n].G1T54,
		G1SL[s][n].G1T55,
		G1SL[s][n].G1T56,
		G1SL[s][n].G1T57,
		G1SL[s][n].G1T58,
		G1SL[s][n].G1T59,
		G1SL[s][n].G1T60,
		G1SL[s][n].G1T61,
		G1SL[s][n].G1T62,
		G1SL[s][n].G1T63,
		G1SL[s][n].G1T64);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank G2SL - multi sector bank
struct G2SL_t	**G2SL;
int 		 *G2SL_NH;
int 		 *G2SL_S;
int 		  G2SL_NS;
void G2SL_Print(){
  int s,n;
  if(G2SL_NS){
    for(s=0;s<G2SL_NS;s++){
      if(G2SL_NH[s]){
        fprintf(stdout,"Group:\tG2SL\tSector:\t%d \tNhits:\t%d\n",G2SL_S[s],G2SL_NH[s]);
        for(n=0;n<G2SL_NH[s];n++){
          fprintf(stdout,"G2T1: %d\tG2T2: %d\tG2T3: %d\tG2T4: %d\tG2T5: %d\tG2T6: %d\tG2T7: %d\tG2T8: %d\tG2T9: %d\tG2T10: %d\tG2T11: %d\tG2T12: %d\tG2T13: %d\tG2T14: %d\tG2T15: %d\tG2T16: %d\tG2T17: %d\tG2T18: %d\tG2T19: %d\tG2T20: %d\tG2T21: %d\tG2T22: %d\tG2T23: %d\tG2T24: %d\tG2T25: %d\tG2T26: %d\tG2T27: %d\tG2T28: %d\tG2T29: %d\tG2T30: %d\tG2T31: %d\tG2T32: %d\tG2T33: %d\tG2T34: %d\tG2T35: %d\tG2T36: %d\tG2T37: %d\tG2T38: %d\tG2T39: %d\tG2T40: %d\tG2T41: %d\tG2T42: %d\tG2T43: %d\tG2T44: %d\tG2T45: %d\tG2T46: %d\tG2T47: %d\tG2T48: %d\tG2T49: %d\tG2T50: %d\tG2T51: %d\tG2T52: %d\tG2T53: %d\tG2T54: %d\tG2T55: %d\tG2T56: %d\tG2T57: %d\tG2T58: %d\tG2T59: %d\tG2T60: %d\tG2T61: %d\tG2T62: %d\tG2T63: %d\tG2T64: %d\t",
		G2SL[s][n].G2T1,
		G2SL[s][n].G2T2,
		G2SL[s][n].G2T3,
		G2SL[s][n].G2T4,
		G2SL[s][n].G2T5,
		G2SL[s][n].G2T6,
		G2SL[s][n].G2T7,
		G2SL[s][n].G2T8,
		G2SL[s][n].G2T9,
		G2SL[s][n].G2T10,
		G2SL[s][n].G2T11,
		G2SL[s][n].G2T12,
		G2SL[s][n].G2T13,
		G2SL[s][n].G2T14,
		G2SL[s][n].G2T15,
		G2SL[s][n].G2T16,
		G2SL[s][n].G2T17,
		G2SL[s][n].G2T18,
		G2SL[s][n].G2T19,
		G2SL[s][n].G2T20,
		G2SL[s][n].G2T21,
		G2SL[s][n].G2T22,
		G2SL[s][n].G2T23,
		G2SL[s][n].G2T24,
		G2SL[s][n].G2T25,
		G2SL[s][n].G2T26,
		G2SL[s][n].G2T27,
		G2SL[s][n].G2T28,
		G2SL[s][n].G2T29,
		G2SL[s][n].G2T30,
		G2SL[s][n].G2T31,
		G2SL[s][n].G2T32,
		G2SL[s][n].G2T33,
		G2SL[s][n].G2T34,
		G2SL[s][n].G2T35,
		G2SL[s][n].G2T36,
		G2SL[s][n].G2T37,
		G2SL[s][n].G2T38,
		G2SL[s][n].G2T39,
		G2SL[s][n].G2T40,
		G2SL[s][n].G2T41,
		G2SL[s][n].G2T42,
		G2SL[s][n].G2T43,
		G2SL[s][n].G2T44,
		G2SL[s][n].G2T45,
		G2SL[s][n].G2T46,
		G2SL[s][n].G2T47,
		G2SL[s][n].G2T48,
		G2SL[s][n].G2T49,
		G2SL[s][n].G2T50,
		G2SL[s][n].G2T51,
		G2SL[s][n].G2T52,
		G2SL[s][n].G2T53,
		G2SL[s][n].G2T54,
		G2SL[s][n].G2T55,
		G2SL[s][n].G2T56,
		G2SL[s][n].G2T57,
		G2SL[s][n].G2T58,
		G2SL[s][n].G2T59,
		G2SL[s][n].G2T60,
		G2SL[s][n].G2T61,
		G2SL[s][n].G2T62,
		G2SL[s][n].G2T63,
		G2SL[s][n].G2T64);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank G3SL - multi sector bank
struct G3SL_t	**G3SL;
int 		 *G3SL_NH;
int 		 *G3SL_S;
int 		  G3SL_NS;
void G3SL_Print(){
  int s,n;
  if(G3SL_NS){
    for(s=0;s<G3SL_NS;s++){
      if(G3SL_NH[s]){
        fprintf(stdout,"Group:\tG3SL\tSector:\t%d \tNhits:\t%d\n",G3SL_S[s],G3SL_NH[s]);
        for(n=0;n<G3SL_NH[s];n++){
          fprintf(stdout,"G3T1: %d\tG3T2: %d\tG3T3: %d\tG3T4: %d\tG3T5: %d\tG3T6: %d\tG3T7: %d\tG3T8: %d\tG3T9: %d\tG3T10: %d\tG3T11: %d\tG3T12: %d\tG3T13: %d\tG3T14: %d\tG3T15: %d\tG3T16: %d\tG3T17: %d\tG3T18: %d\tG3T19: %d\tG3T20: %d\tG3T21: %d\tG3T22: %d\tG3T23: %d\tG3T24: %d\tG3T25: %d\tG3T26: %d\tG3T27: %d\tG3T28: %d\tG3T29: %d\tG3T30: %d\tG3T31: %d\tG3T32: %d\tG3T33: %d\tG3T34: %d\tG3T35: %d\tG3T36: %d\tG3T37: %d\tG3T38: %d\tG3T39: %d\tG3T40: %d\tG3T41: %d\tG3T42: %d\tG3T43: %d\tG3T44: %d\tG3T45: %d\tG3T46: %d\tG3T47: %d\tG3T48: %d\tG3T49: %d\tG3T50: %d\tG3T51: %d\tG3T52: %d\tG3T53: %d\tG3T54: %d\tG3T55: %d\tG3T56: %d\tG3T57: %d\tG3T58: %d\tG3T59: %d\tG3T60: %d\tG3T61: %d\tG3T62: %d\tG3T63: %d\tG3T64: %d\t",
		G3SL[s][n].G3T1,
		G3SL[s][n].G3T2,
		G3SL[s][n].G3T3,
		G3SL[s][n].G3T4,
		G3SL[s][n].G3T5,
		G3SL[s][n].G3T6,
		G3SL[s][n].G3T7,
		G3SL[s][n].G3T8,
		G3SL[s][n].G3T9,
		G3SL[s][n].G3T10,
		G3SL[s][n].G3T11,
		G3SL[s][n].G3T12,
		G3SL[s][n].G3T13,
		G3SL[s][n].G3T14,
		G3SL[s][n].G3T15,
		G3SL[s][n].G3T16,
		G3SL[s][n].G3T17,
		G3SL[s][n].G3T18,
		G3SL[s][n].G3T19,
		G3SL[s][n].G3T20,
		G3SL[s][n].G3T21,
		G3SL[s][n].G3T22,
		G3SL[s][n].G3T23,
		G3SL[s][n].G3T24,
		G3SL[s][n].G3T25,
		G3SL[s][n].G3T26,
		G3SL[s][n].G3T27,
		G3SL[s][n].G3T28,
		G3SL[s][n].G3T29,
		G3SL[s][n].G3T30,
		G3SL[s][n].G3T31,
		G3SL[s][n].G3T32,
		G3SL[s][n].G3T33,
		G3SL[s][n].G3T34,
		G3SL[s][n].G3T35,
		G3SL[s][n].G3T36,
		G3SL[s][n].G3T37,
		G3SL[s][n].G3T38,
		G3SL[s][n].G3T39,
		G3SL[s][n].G3T40,
		G3SL[s][n].G3T41,
		G3SL[s][n].G3T42,
		G3SL[s][n].G3T43,
		G3SL[s][n].G3T44,
		G3SL[s][n].G3T45,
		G3SL[s][n].G3T46,
		G3SL[s][n].G3T47,
		G3SL[s][n].G3T48,
		G3SL[s][n].G3T49,
		G3SL[s][n].G3T50,
		G3SL[s][n].G3T51,
		G3SL[s][n].G3T52,
		G3SL[s][n].G3T53,
		G3SL[s][n].G3T54,
		G3SL[s][n].G3T55,
		G3SL[s][n].G3T56,
		G3SL[s][n].G3T57,
		G3SL[s][n].G3T58,
		G3SL[s][n].G3T59,
		G3SL[s][n].G3T60,
		G3SL[s][n].G3T61,
		G3SL[s][n].G3T62,
		G3SL[s][n].G3T63,
		G3SL[s][n].G3T64);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank G4SL - multi sector bank
struct G4SL_t	**G4SL;
int 		 *G4SL_NH;
int 		 *G4SL_S;
int 		  G4SL_NS;
void G4SL_Print(){
  int s,n;
  if(G4SL_NS){
    for(s=0;s<G4SL_NS;s++){
      if(G4SL_NH[s]){
        fprintf(stdout,"Group:\tG4SL\tSector:\t%d \tNhits:\t%d\n",G4SL_S[s],G4SL_NH[s]);
        for(n=0;n<G4SL_NH[s];n++){
          fprintf(stdout,"G4T1: %d\tG4T2: %d\tG4T3: %d\tG4T4: %d\tG4T5: %d\tG4T6: %d\tG4T7: %d\tG4T8: %d\tG4T9: %d\tG4T10: %d\tG4T11: %d\tG4T12: %d\tG4T13: %d\tG4T14: %d\tG4T15: %d\tG4T16: %d\tG4T17: %d\tG4T18: %d\tG4T19: %d\tG4T20: %d\tG4T21: %d\tG4T22: %d\tG4T23: %d\tG4T24: %d\tG4T25: %d\tG4T26: %d\tG4T27: %d\tG4T28: %d\tG4T29: %d\tG4T30: %d\tG4T31: %d\tG4T32: %d\tG4T33: %d\tG4T34: %d\tG4T35: %d\tG4T36: %d\tG4T37: %d\tG4T38: %d\tG4T39: %d\tG4T40: %d\tG4T41: %d\tG4T42: %d\tG4T43: %d\tG4T44: %d\tG4T45: %d\tG4T46: %d\tG4T47: %d\tG4T48: %d\tG4T49: %d\tG4T50: %d\tG4T51: %d\tG4T52: %d\tG4T53: %d\tG4T54: %d\tG4T55: %d\tG4T56: %d\tG4T57: %d\tG4T58: %d\tG4T59: %d\tG4T60: %d\tG4T61: %d\tG4T62: %d\tG4T63: %d\tG4T64: %d\t",
		G4SL[s][n].G4T1,
		G4SL[s][n].G4T2,
		G4SL[s][n].G4T3,
		G4SL[s][n].G4T4,
		G4SL[s][n].G4T5,
		G4SL[s][n].G4T6,
		G4SL[s][n].G4T7,
		G4SL[s][n].G4T8,
		G4SL[s][n].G4T9,
		G4SL[s][n].G4T10,
		G4SL[s][n].G4T11,
		G4SL[s][n].G4T12,
		G4SL[s][n].G4T13,
		G4SL[s][n].G4T14,
		G4SL[s][n].G4T15,
		G4SL[s][n].G4T16,
		G4SL[s][n].G4T17,
		G4SL[s][n].G4T18,
		G4SL[s][n].G4T19,
		G4SL[s][n].G4T20,
		G4SL[s][n].G4T21,
		G4SL[s][n].G4T22,
		G4SL[s][n].G4T23,
		G4SL[s][n].G4T24,
		G4SL[s][n].G4T25,
		G4SL[s][n].G4T26,
		G4SL[s][n].G4T27,
		G4SL[s][n].G4T28,
		G4SL[s][n].G4T29,
		G4SL[s][n].G4T30,
		G4SL[s][n].G4T31,
		G4SL[s][n].G4T32,
		G4SL[s][n].G4T33,
		G4SL[s][n].G4T34,
		G4SL[s][n].G4T35,
		G4SL[s][n].G4T36,
		G4SL[s][n].G4T37,
		G4SL[s][n].G4T38,
		G4SL[s][n].G4T39,
		G4SL[s][n].G4T40,
		G4SL[s][n].G4T41,
		G4SL[s][n].G4T42,
		G4SL[s][n].G4T43,
		G4SL[s][n].G4T44,
		G4SL[s][n].G4T45,
		G4SL[s][n].G4T46,
		G4SL[s][n].G4T47,
		G4SL[s][n].G4T48,
		G4SL[s][n].G4T49,
		G4SL[s][n].G4T50,
		G4SL[s][n].G4T51,
		G4SL[s][n].G4T52,
		G4SL[s][n].G4T53,
		G4SL[s][n].G4T54,
		G4SL[s][n].G4T55,
		G4SL[s][n].G4T56,
		G4SL[s][n].G4T57,
		G4SL[s][n].G4T58,
		G4SL[s][n].G4T59,
		G4SL[s][n].G4T60,
		G4SL[s][n].G4T61,
		G4SL[s][n].G4T62,
		G4SL[s][n].G4T63,
		G4SL[s][n].G4T64);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank GPAR - multi sector bank
struct GPAR_t	**GPAR;
int 		 *GPAR_NH;
int 		 *GPAR_S;
int 		  GPAR_NS;
void GPAR_Print(){
  int s,n;
  if(GPAR_NS){
    for(s=0;s<GPAR_NS;s++){
      if(GPAR_NH[s]){
        fprintf(stdout,"Group:\tGPAR\tSector:\t%d \tNhits:\t%d\n",GPAR_S[s],GPAR_NH[s]);
        for(n=0;n<GPAR_NH[s];n++){
          fprintf(stdout,"fval: %8.3f\tival: %d\tchar1: %.20s\t",
		GPAR[s][n].fval,
		GPAR[s][n].ival,
		GPAR[s][n].char1);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank GPID - single sector bank
struct GPID_t	 *GPID;
int 		  GPID_NH;
int 		  GPID_S;
void GPID_Print(){
  int n;
  if(GPID_NH){
    fprintf(stdout,"Group:\tGPID\tSector:\t%d \tNhits:\t%d\n",GPID_S,GPID_NH);
    for(n=0;n<GPID_NH;n++){
      fprintf(stdout,"pid: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tE: %8.3f\tpx: %8.3f\tpy: %8.3f\tpz: %8.3f\tq: %d\ttrkid: %d\tsec: %d\tpaddle: %d\tdedx: %8.3f\tbeta: %8.3f\tsc_stat: %d\tsc_time: %8.3f\tsc_len: %8.3f\tst_stat: %d\tst_time: %8.3f\tst_len: %8.3f\tmass: %8.3f\tmass_ref: %d\tbetam: %8.3f\tepho: %8.3f\ttpho: %8.3f\ttagrid: %d\tngrf: %d\tppid: %d\t",
		GPID[n].pid,
		GPID[n].x,
		GPID[n].y,
		GPID[n].z,
		GPID[n].E,
		GPID[n].px,
		GPID[n].py,
		GPID[n].pz,
		GPID[n].q,
		GPID[n].trkid,
		GPID[n].sec,
		GPID[n].paddle,
		GPID[n].dedx,
		GPID[n].beta,
		GPID[n].sc_stat,
		GPID[n].sc_time,
		GPID[n].sc_len,
		GPID[n].st_stat,
		GPID[n].st_time,
		GPID[n].st_len,
		GPID[n].mass,
		GPID[n].mass_ref,
		GPID[n].betam,
		GPID[n].epho,
		GPID[n].tpho,
		GPID[n].tagrid,
		GPID[n].ngrf,
		GPID[n].ppid);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank GP_X - single sector bank
struct GP_X_t	 *GP_X;
int 		  GP_X_NH;
int 		  GP_X_S;
void GP_X_Print(){
  int n;
  if(GP_X_NH){
    fprintf(stdout,"Group:\tGP_X\tSector:\t%d \tNhits:\t%d\n",GP_X_S,GP_X_NH);
    for(n=0;n<GP_X_NH;n++){
      fprintf(stdout,"value: %8.3f\t",
		GP_X[n].value);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank GP_Y - single sector bank
struct GP_Y_t	 *GP_Y;
int 		  GP_Y_NH;
int 		  GP_Y_S;
void GP_Y_Print(){
  int n;
  if(GP_Y_NH){
    fprintf(stdout,"Group:\tGP_Y\tSector:\t%d \tNhits:\t%d\n",GP_Y_S,GP_Y_NH);
    for(n=0;n<GP_Y_NH;n++){
      fprintf(stdout,"value: %8.3f\t",
		GP_Y[n].value);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HBER - single sector bank
struct HBER_t	 *HBER;
int 		  HBER_NH;
int 		  HBER_S;
void HBER_Print(){
  int n;
  if(HBER_NH){
    fprintf(stdout,"Group:\tHBER\tSector:\t%d \tNhits:\t%d\n",HBER_S,HBER_NH);
    for(n=0;n<HBER_NH;n++){
      fprintf(stdout,"q_over_p: %8.3f\tlambda: %8.3f\tphi: %8.3f\td0: %8.3f\tz0: %8.3f\tc11: %8.3f\tc12: %8.3f\tc13: %8.3f\tc14: %8.3f\tc15: %8.3f\tc22: %8.3f\tc23: %8.3f\tc24: %8.3f\tc25: %8.3f\tc33: %8.3f\tc34: %8.3f\tc35: %8.3f\tc44: %8.3f\tc45: %8.3f\tc55: %8.3f\tchi2: %8.3f\tlayinfo1: %d\tlayinfo2: %d\t",
		HBER[n].q_over_p,
		HBER[n].lambda,
		HBER[n].phi,
		HBER[n].d0,
		HBER[n].z0,
		HBER[n].c11,
		HBER[n].c12,
		HBER[n].c13,
		HBER[n].c14,
		HBER[n].c15,
		HBER[n].c22,
		HBER[n].c23,
		HBER[n].c24,
		HBER[n].c25,
		HBER[n].c33,
		HBER[n].c34,
		HBER[n].c35,
		HBER[n].c44,
		HBER[n].c45,
		HBER[n].c55,
		HBER[n].chi2,
		HBER[n].layinfo1,
		HBER[n].layinfo2);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HBID - single sector bank
struct HBID_t	 *HBID;
int 		  HBID_NH;
int 		  HBID_S;
void HBID_Print(){
  int n;
  if(HBID_NH){
    fprintf(stdout,"Group:\tHBID\tSector:\t%d \tNhits:\t%d\n",HBID_S,HBID_NH);
    for(n=0;n<HBID_NH;n++){
      fprintf(stdout,"track: %d\tsec: %d\tbeta: %8.3f\tvtime: %8.3f\tsc_stat: %d\tsc_id: %d\tsc_time: %8.3f\tsc_qual: %8.3f\tsc_vtime: %8.3f\tsc_beta: %8.3f\tcc_stat: %d\tcc_id: %d\tcc_time: %8.3f\tcc_qual: %8.3f\tcc_vtime: %8.3f\tcc_beta: %8.3f\tec_stat: %d\tec_id: %d\tec_time: %8.3f\tec_qual: %8.3f\tec_vtime: %8.3f\tec_beta: %8.3f\tst_stat: %d\tst_id: %d\tst_time: %8.3f\tst_qual: %8.3f\tst_vtime: %8.3f\tst_beta: %8.3f\tlac_stat: %d\tlac_id: %d\tlac_time: %8.3f\tlac_qual: %8.3f\tlac_vtime: %8.3f\tlac_beta: %8.3f\t",
		HBID[n].track,
		HBID[n].sec,
		HBID[n].beta,
		HBID[n].vtime,
		HBID[n].sc_stat,
		HBID[n].sc_id,
		HBID[n].sc_time,
		HBID[n].sc_qual,
		HBID[n].sc_vtime,
		HBID[n].sc_beta,
		HBID[n].cc_stat,
		HBID[n].cc_id,
		HBID[n].cc_time,
		HBID[n].cc_qual,
		HBID[n].cc_vtime,
		HBID[n].cc_beta,
		HBID[n].ec_stat,
		HBID[n].ec_id,
		HBID[n].ec_time,
		HBID[n].ec_qual,
		HBID[n].ec_vtime,
		HBID[n].ec_beta,
		HBID[n].st_stat,
		HBID[n].st_id,
		HBID[n].st_time,
		HBID[n].st_qual,
		HBID[n].st_vtime,
		HBID[n].st_beta,
		HBID[n].lac_stat,
		HBID[n].lac_id,
		HBID[n].lac_time,
		HBID[n].lac_qual,
		HBID[n].lac_vtime,
		HBID[n].lac_beta);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HBLA - multi sector bank
struct HBLA_t	**HBLA;
int 		 *HBLA_NH;
int 		 *HBLA_S;
int 		  HBLA_NS;
void HBLA_Print(){
  int s,n;
  if(HBLA_NS){
    for(s=0;s<HBLA_NS;s++){
      if(HBLA_NH[s]){
        fprintf(stdout,"Group:\tHBLA\tSector:\t%d \tNhits:\t%d\n",HBLA_S[s],HBLA_NH[s]);
        for(n=0;n<HBLA_NH[s];n++){
          fprintf(stdout,"trk_pln: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tBx: %8.3f\tBy: %8.3f\tBz: %8.3f\ttlen: %8.3f\tdc1: %d\tstat: %d\twire: %d\tdtime: %8.3f\talpha: %8.3f\twlen: %8.3f\tsgdoca: %8.3f\tfitdoca: %8.3f\t",
		HBLA[s][n].trk_pln,
		HBLA[s][n].x,
		HBLA[s][n].y,
		HBLA[s][n].z,
		HBLA[s][n].Bx,
		HBLA[s][n].By,
		HBLA[s][n].Bz,
		HBLA[s][n].tlen,
		HBLA[s][n].dc1,
		HBLA[s][n].stat,
		HBLA[s][n].wire,
		HBLA[s][n].dtime,
		HBLA[s][n].alpha,
		HBLA[s][n].wlen,
		HBLA[s][n].sgdoca,
		HBLA[s][n].fitdoca);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank HBTB - multi sector bank
struct HBTB_t	**HBTB;
int 		 *HBTB_NH;
int 		 *HBTB_S;
int 		  HBTB_NS;
void HBTB_Print(){
  int s,n;
  if(HBTB_NS){
    for(s=0;s<HBTB_NS;s++){
      if(HBTB_NH[s]){
        fprintf(stdout,"Group:\tHBTB\tSector:\t%d \tNhits:\t%d\n",HBTB_S[s],HBTB_NH[s]);
        for(n=0;n<HBTB_NH[s];n++){
          fprintf(stdout,"trk: %d\tsly: %d\ticl: %d\tisg: %d\tid1: %d\tid2: %d\tid3: %d\tid4: %d\tid5: %d\tid6: %d\t",
		HBTB[s][n].trk,
		HBTB[s][n].sly,
		HBTB[s][n].icl,
		HBTB[s][n].isg,
		HBTB[s][n].id1,
		HBTB[s][n].id2,
		HBTB[s][n].id3,
		HBTB[s][n].id4,
		HBTB[s][n].id5,
		HBTB[s][n].id6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank HBTR - single sector bank
struct HBTR_t	 *HBTR;
int 		  HBTR_NH;
int 		  HBTR_S;
void HBTR_Print(){
  int n;
  if(HBTR_NH){
    fprintf(stdout,"Group:\tHBTR\tSector:\t%d \tNhits:\t%d\n",HBTR_S,HBTR_NH);
    for(n=0;n<HBTR_NH;n++){
      fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tpx: %8.3f\tpy: %8.3f\tpz: %8.3f\tq: %8.3f\tchi2: %8.3f\titr_sec: %d\t",
		HBTR[n].x,
		HBTR[n].y,
		HBTR[n].z,
		HBTR[n].px,
		HBTR[n].py,
		HBTR[n].pz,
		HBTR[n].q,
		HBTR[n].chi2,
		HBTR[n].itr_sec);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HCAL - multi sector bank
struct HCAL_t	**HCAL;
int 		 *HCAL_NH;
int 		 *HCAL_S;
int 		  HCAL_NS;
void HCAL_Print(){
  int s,n;
  if(HCAL_NS){
    for(s=0;s<HCAL_NS;s++){
      if(HCAL_NH[s]){
        fprintf(stdout,"Group:\tHCAL\tSector:\t%d \tNhits:\t%d\n",HCAL_S[s],HCAL_NH[s]);
        for(n=0;n<HCAL_NH[s];n++){
          fprintf(stdout,"VERSION: %d\tRUN_LOW: %d\tRUN_HIGH: %d\tCATIME: %d\tROCCA: %d\t",
		HCAL[s][n].VERSION,
		HCAL[s][n].RUN_LOW,
		HCAL[s][n].RUN_HIGH,
		HCAL[s][n].CATIME,
		HCAL[s][n].ROCCA);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank HDPL - multi sector bank
struct HDPL_t	**HDPL;
int 		 *HDPL_NH;
int 		 *HDPL_S;
int 		  HDPL_NS;
void HDPL_Print(){
  int s,n;
  if(HDPL_NS){
    for(s=0;s<HDPL_NS;s++){
      if(HDPL_NH[s]){
        fprintf(stdout,"Group:\tHDPL\tSector:\t%d \tNhits:\t%d\n",HDPL_S[s],HDPL_NH[s]);
        for(n=0;n<HDPL_NH[s];n++){
          fprintf(stdout,"trk_pln: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\ttlen: %8.3f\t",
		HDPL[s][n].trk_pln,
		HDPL[s][n].x,
		HDPL[s][n].y,
		HDPL[s][n].z,
		HDPL[s][n].cx,
		HDPL[s][n].cy,
		HDPL[s][n].cz,
		HDPL[s][n].tlen);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank HEAD - single sector bank
struct HEAD_t	 *HEAD;
int 		  HEAD_NH;
int 		  HEAD_S;
void HEAD_Print(){
  int n;
  if(HEAD_NH){
    fprintf(stdout,"Group:\tHEAD\tSector:\t%d \tNhits:\t%d\n",HEAD_S,HEAD_NH);
    for(n=0;n<HEAD_NH;n++){
      fprintf(stdout,"VERSION: %d\tNRUN: %d\tNEVENT: %d\tTIME: %24.24s\tTYPE: %d\tROC: %d\tEVTCLASS: %d\tTRIGBITS: 0x%4.4X\t",
		HEAD[n].VERSION,
		HEAD[n].NRUN,
		HEAD[n].NEVENT,
		ctime((const long*)&(HEAD[n].TIME)),
		HEAD[n].TYPE,
		HEAD[n].ROC,
		HEAD[n].EVTCLASS,
		HEAD[n].TRIGBITS);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HEVT - single sector bank
struct HEVT_t	 *HEVT;
int 		  HEVT_NH;
int 		  HEVT_S;
void HEVT_Print(){
  int n;
  if(HEVT_NH){
    fprintf(stdout,"Group:\tHEVT\tSector:\t%d \tNhits:\t%d\n",HEVT_S,HEVT_NH);
    for(n=0;n<HEVT_NH;n++){
      fprintf(stdout,"ESTATUS: %d\tNRUN: %d\tNEVENT: %d\tTYPE: %d\tNPGP: %d\tTRGPRS: %d\tFC: %8.3f\tFCG: %8.3f\tTG: %8.3f\tSTT: %8.3f\tRF1: %8.3f\tRF2: %8.3f\tCON1: %8.3f\tCON2: %8.3f\tCON3: %8.3f\tPTIME: %d\t",
		HEVT[n].ESTATUS,
		HEVT[n].NRUN,
		HEVT[n].NEVENT,
		HEVT[n].TYPE,
		HEVT[n].NPGP,
		HEVT[n].TRGPRS,
		HEVT[n].FC,
		HEVT[n].FCG,
		HEVT[n].TG,
		HEVT[n].STT,
		HEVT[n].RF1,
		HEVT[n].RF2,
		HEVT[n].CON1,
		HEVT[n].CON2,
		HEVT[n].CON3,
		HEVT[n].PTIME);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank HLS - multi sector bank
struct HLS_t	**HLS;
int 		 *HLS_NH;
int 		 *HLS_S;
int 		  HLS_NS;
void HLS_Print(){
  int s,n;
  if(HLS_NS){
    for(s=0;s<HLS_NS;s++){
      if(HLS_NH[s]){
        fprintf(stdout,"Group:\tHLS\tSector:\t%d \tNhits:\t%d\n",HLS_S[s],HLS_NH[s]);
        for(n=0;n<HLS_NH[s];n++){
          fprintf(stdout,"S1: %d\tS2: %d\tS3: %d\tS4: %d\tS5: %d\tS6: %d\tS7: %d\tS8: %d\tS9: %d\tS10: %d\tS11: %d\tS12: %d\tS13: %d\tS14: %d\tS15: %d\tS16: %d\t",
		HLS[s][n].S1,
		HLS[s][n].S2,
		HLS[s][n].S3,
		HLS[s][n].S4,
		HLS[s][n].S5,
		HLS[s][n].S6,
		HLS[s][n].S7,
		HLS[s][n].S8,
		HLS[s][n].S9,
		HLS[s][n].S10,
		HLS[s][n].S11,
		HLS[s][n].S12,
		HLS[s][n].S13,
		HLS[s][n].S14,
		HLS[s][n].S15,
		HLS[s][n].S16);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank IC - multi sector bank
struct IC_t	**IC;
int 		 *IC_NH;
int 		 *IC_S;
int 		  IC_NS;
void IC_Print(){
  int s,n;
  if(IC_NS){
    for(s=0;s<IC_NS;s++){
      if(IC_NH[s]){
        fprintf(stdout,"Group:\tIC\tSector:\t%d \tNhits:\t%d\n",IC_S[s],IC_NH[s]);
        for(n=0;n<IC_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		IC[s][n].ID,
		IC[s][n].TDC,
		IC[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ICHB - multi sector bank
struct ICHB_t	**ICHB;
int 		 *ICHB_NH;
int 		 *ICHB_S;
int 		  ICHB_NS;
void ICHB_Print(){
  int s,n;
  if(ICHB_NS){
    for(s=0;s<ICHB_NS;s++){
      if(ICHB_NH[s]){
        fprintf(stdout,"Group:\tICHB\tSector:\t%d \tNhits:\t%d\n",ICHB_S[s],ICHB_NH[s]);
        for(n=0;n<ICHB_NH[s];n++){
          fprintf(stdout,"Eclust: %8.3f\tEclmax: %8.3f\tTclust: %8.3f\tTclmax: %8.3f\txclust: %8.3f\tyclust: %8.3f\tzclust: %8.3f\txclmax: %8.3f\tyclmax: %8.3f\tzclmax: %8.3f\tres1: %8.3f\tres2: %8.3f\tres3: %8.3f\tncryst: %d\t",
		ICHB[s][n].Eclust,
		ICHB[s][n].Eclmax,
		ICHB[s][n].Tclust,
		ICHB[s][n].Tclmax,
		ICHB[s][n].xclust,
		ICHB[s][n].yclust,
		ICHB[s][n].zclust,
		ICHB[s][n].xclmax,
		ICHB[s][n].yclmax,
		ICHB[s][n].zclmax,
		ICHB[s][n].res1,
		ICHB[s][n].res2,
		ICHB[s][n].res3,
		ICHB[s][n].ncryst);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank KFIT - multi sector bank
struct KFIT_t	**KFIT;
int 		 *KFIT_NH;
int 		 *KFIT_S;
int 		  KFIT_NS;
void KFIT_Print(){
  int s,n;
  if(KFIT_NS){
    for(s=0;s<KFIT_NS;s++){
      if(KFIT_NH[s]){
        fprintf(stdout,"Group:\tKFIT\tSector:\t%d \tNhits:\t%d\n",KFIT_S[s],KFIT_NH[s]);
        for(n=0;n<KFIT_NH[s];n++){
          fprintf(stdout,"momenta_f: %8.3f\tlambda_f: %8.3f\tphi_f: %8.3f\td0_f: %8.3f\tz0_f: %8.3f\tsigma11: %8.3f\tsigma12: %8.3f\tsigma13: %8.3f\tsigma14: %8.3f\tsigma15: %8.3f\tsigma22: %8.3f\tsigma23: %8.3f\tsigma24: %8.3f\tsigma25: %8.3f\tsigma33: %8.3f\tsigma34: %8.3f\tsigma35: %8.3f\tsigma44: %8.3f\tsigma45: %8.3f\tsigma55: %8.3f\tcov11: %8.3f\tcov12: %8.3f\tcov13: %8.3f\tcov14: %8.3f\tcov15: %8.3f\tcov21: %8.3f\tcov22: %8.3f\tcov23: %8.3f\tcov24: %8.3f\tcov25: %8.3f\tcov31: %8.3f\tcov32: %8.3f\tcov33: %8.3f\tcov34: %8.3f\tcov35: %8.3f\tcov41: %8.3f\tcov42: %8.3f\tcov43: %8.3f\tcov44: %8.3f\tcov45: %8.3f\tcov51: %8.3f\tcov52: %8.3f\tcov53: %8.3f\tcov54: %8.3f\tcov55: %8.3f\tchi2piece: %8.3f\t",
		KFIT[s][n].momenta_f,
		KFIT[s][n].lambda_f,
		KFIT[s][n].phi_f,
		KFIT[s][n].d0_f,
		KFIT[s][n].z0_f,
		KFIT[s][n].sigma11,
		KFIT[s][n].sigma12,
		KFIT[s][n].sigma13,
		KFIT[s][n].sigma14,
		KFIT[s][n].sigma15,
		KFIT[s][n].sigma22,
		KFIT[s][n].sigma23,
		KFIT[s][n].sigma24,
		KFIT[s][n].sigma25,
		KFIT[s][n].sigma33,
		KFIT[s][n].sigma34,
		KFIT[s][n].sigma35,
		KFIT[s][n].sigma44,
		KFIT[s][n].sigma45,
		KFIT[s][n].sigma55,
		KFIT[s][n].cov11,
		KFIT[s][n].cov12,
		KFIT[s][n].cov13,
		KFIT[s][n].cov14,
		KFIT[s][n].cov15,
		KFIT[s][n].cov21,
		KFIT[s][n].cov22,
		KFIT[s][n].cov23,
		KFIT[s][n].cov24,
		KFIT[s][n].cov25,
		KFIT[s][n].cov31,
		KFIT[s][n].cov32,
		KFIT[s][n].cov33,
		KFIT[s][n].cov34,
		KFIT[s][n].cov35,
		KFIT[s][n].cov41,
		KFIT[s][n].cov42,
		KFIT[s][n].cov43,
		KFIT[s][n].cov44,
		KFIT[s][n].cov45,
		KFIT[s][n].cov51,
		KFIT[s][n].cov52,
		KFIT[s][n].cov53,
		KFIT[s][n].cov54,
		KFIT[s][n].cov55,
		KFIT[s][n].chi2piece);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank L1PG - multi sector bank
struct L1PG_t	**L1PG;
int 		 *L1PG_NH;
int 		 *L1PG_S;
int 		  L1PG_NS;
void L1PG_Print(){
  int s,n;
  if(L1PG_NS){
    for(s=0;s<L1PG_NS;s++){
      if(L1PG_NH[s]){
        fprintf(stdout,"Group:\tL1PG\tSector:\t%d \tNhits:\t%d\n",L1PG_S[s],L1PG_NH[s]);
        for(n=0;n<L1PG_NH[s];n++){
          fprintf(stdout,"l1prog: %.4s\t",
		L1PG[s][n].l1prog);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank L2H - multi sector bank
struct L2H_t	**L2H;
int 		 *L2H_NH;
int 		 *L2H_S;
int 		  L2H_NS;
void L2H_Print(){
  int s,n;
  if(L2H_NS){
    for(s=0;s<L2H_NS;s++){
      if(L2H_NH[s]){
        fprintf(stdout,"Group:\tL2H\tSector:\t%d \tNhits:\t%d\n",L2H_S[s],L2H_NH[s]);
        for(n=0;n<L2H_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		L2H[s][n].ID,
		L2H[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank L2S - multi sector bank
struct L2S_t	**L2S;
int 		 *L2S_NH;
int 		 *L2S_S;
int 		  L2S_NS;
void L2S_Print(){
  int s,n;
  if(L2S_NS){
    for(s=0;s<L2S_NS;s++){
      if(L2S_NH[s]){
        fprintf(stdout,"Group:\tL2S\tSector:\t%d \tNhits:\t%d\n",L2S_S[s],L2S_NH[s]);
        for(n=0;n<L2S_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		L2S[s][n].ID,
		L2S[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank LASR - multi sector bank
struct LASR_t	**LASR;
int 		 *LASR_NH;
int 		 *LASR_S;
int 		  LASR_NS;
void LASR_Print(){
  int s,n;
  if(LASR_NS){
    for(s=0;s<LASR_NS;s++){
      if(LASR_NH[s]){
        fprintf(stdout,"Group:\tLASR\tSector:\t%d \tNhits:\t%d\n",LASR_S[s],LASR_NH[s]);
        for(n=0;n<LASR_NH[s];n++){
          fprintf(stdout,"ID: %d\tstat1: %d\tstat2: %d\tstat3: %d\t",
		LASR[s][n].ID,
		LASR[s][n].stat1,
		LASR[s][n].stat2,
		LASR[s][n].stat3);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank LCDI - multi sector bank
struct LCDI_t	**LCDI;
int 		 *LCDI_NH;
int 		 *LCDI_S;
int 		  LCDI_NS;
void LCDI_Print(){
  int s,n;
  if(LCDI_NS){
    for(s=0;s<LCDI_NS;s++){
      if(LCDI_NH[s]){
        fprintf(stdout,"Group:\tLCDI\tSector:\t%d \tNhits:\t%d\n",LCDI_S[s],LCDI_NH[s]);
        for(n=0;n<LCDI_NH[s];n++){
          fprintf(stdout,"crate: %d\tslot: %d\tmask: %d\tthreshold: %d\t",
		LCDI[s][n].crate,
		LCDI[s][n].slot,
		LCDI[s][n].mask,
		LCDI[s][n].threshold);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank LCPB - single sector bank
struct LCPB_t	 *LCPB;
int 		  LCPB_NH;
int 		  LCPB_S;
void LCPB_Print(){
  int n;
  if(LCPB_NH){
    fprintf(stdout,"Group:\tLCPB\tSector:\t%d \tNhits:\t%d\n",LCPB_S,LCPB_NH);
    for(n=0;n<LCPB_NH;n++){
      fprintf(stdout,"ScHt: %d\tEtot: %8.3f\tTime: %8.3f\tPath: %8.3f\tX: %8.3f\tY: %8.3f\tZ: %8.3f\tChi2LC: %8.3f\tStatus: %d\tEin: %8.3f\t",
		LCPB[n].ScHt,
		LCPB[n].Etot,
		LCPB[n].Time,
		LCPB[n].Path,
		LCPB[n].X,
		LCPB[n].Y,
		LCPB[n].Z,
		LCPB[n].Chi2LC,
		LCPB[n].Status,
		LCPB[n].Ein);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank LOWQ - multi sector bank
struct LOWQ_t	**LOWQ;
int 		 *LOWQ_NH;
int 		 *LOWQ_S;
int 		  LOWQ_NS;
void LOWQ_Print(){
  int s,n;
  if(LOWQ_NS){
    for(s=0;s<LOWQ_NS;s++){
      if(LOWQ_NH[s]){
        fprintf(stdout,"Group:\tLOWQ\tSector:\t%d \tNhits:\t%d\n",LOWQ_S[s],LOWQ_NH[s]);
        for(n=0;n<LOWQ_NH[s];n++){
          fprintf(stdout,"pid: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tE: %8.3f\tpx: %8.3f\tpy: %8.3f\tpz: %8.3f\tq: %8.3f\tqpid: %8.3f\tqtrk: %8.3f\tflags: %d\t",
		LOWQ[s][n].pid,
		LOWQ[s][n].x,
		LOWQ[s][n].y,
		LOWQ[s][n].z,
		LOWQ[s][n].E,
		LOWQ[s][n].px,
		LOWQ[s][n].py,
		LOWQ[s][n].pz,
		LOWQ[s][n].q,
		LOWQ[s][n].qpid,
		LOWQ[s][n].qtrk,
		LOWQ[s][n].flags);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MCEV - multi sector bank
struct MCEV_t	**MCEV;
int 		 *MCEV_NH;
int 		 *MCEV_S;
int 		  MCEV_NS;
void MCEV_Print(){
  int s,n;
  if(MCEV_NS){
    for(s=0;s<MCEV_NS;s++){
      if(MCEV_NH[s]){
        fprintf(stdout,"Group:\tMCEV\tSector:\t%d \tNhits:\t%d\n",MCEV_S[s],MCEV_NH[s]);
        for(n=0;n<MCEV_NH[s];n++){
          fprintf(stdout,"I1: %d\ti2: %d\t",
		MCEV[s][n].I1,
		MCEV[s][n].i2);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MCHD - multi sector bank
struct MCHD_t	**MCHD;
int 		 *MCHD_NH;
int 		 *MCHD_S;
int 		  MCHD_NS;
void MCHD_Print(){
  int s,n;
  if(MCHD_NS){
    for(s=0;s<MCHD_NS;s++){
      if(MCHD_NH[s]){
        fprintf(stdout,"Group:\tMCHD\tSector:\t%d \tNhits:\t%d\n",MCHD_S[s],MCHD_NH[s]);
        for(n=0;n<MCHD_NH[s];n++){
          fprintf(stdout,"NRUN: %d\tNEVENT: %d\tTIME: %d\tTYPE: %d\tREACTYPE: %d\tWEIGHT: %8.3f\tW: %8.3f\tQ2: %8.3f\tE_PHOT: %8.3f\tPX_PHOT: %8.3f\tPY_PHOT: %8.3f\tPZ_PHOT: %8.3f\tE_TARG: %8.3f\tPX_TARG: %8.3f\tPY_TARG: %8.3f\tPZ_TARG: %8.3f\t",
		MCHD[s][n].NRUN,
		MCHD[s][n].NEVENT,
		MCHD[s][n].TIME,
		MCHD[s][n].TYPE,
		MCHD[s][n].REACTYPE,
		MCHD[s][n].WEIGHT,
		MCHD[s][n].W,
		MCHD[s][n].Q2,
		MCHD[s][n].E_PHOT,
		MCHD[s][n].PX_PHOT,
		MCHD[s][n].PY_PHOT,
		MCHD[s][n].PZ_PHOT,
		MCHD[s][n].E_TARG,
		MCHD[s][n].PX_TARG,
		MCHD[s][n].PY_TARG,
		MCHD[s][n].PZ_TARG);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MCTK - multi sector bank
struct MCTK_t	**MCTK;
int 		 *MCTK_NH;
int 		 *MCTK_S;
int 		  MCTK_NS;
void MCTK_Print(){
  int s,n;
  if(MCTK_NS){
    for(s=0;s<MCTK_NS;s++){
      if(MCTK_NH[s]){
        fprintf(stdout,"Group:\tMCTK\tSector:\t%d \tNhits:\t%d\n",MCTK_S[s],MCTK_NH[s]);
        for(n=0;n<MCTK_NH[s];n++){
          fprintf(stdout,"cx: %8.3f\tcy: %8.3f\tcz: %8.3f\tpmom: %8.3f\tmass: %8.3f\tcharge: %8.3f\tid: %d\tflag: %d\tbeg_vtx: %d\tend_vtx: %d\tparent: %d\t",
		MCTK[s][n].cx,
		MCTK[s][n].cy,
		MCTK[s][n].cz,
		MCTK[s][n].pmom,
		MCTK[s][n].mass,
		MCTK[s][n].charge,
		MCTK[s][n].id,
		MCTK[s][n].flag,
		MCTK[s][n].beg_vtx,
		MCTK[s][n].end_vtx,
		MCTK[s][n].parent);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MCVX - multi sector bank
struct MCVX_t	**MCVX;
int 		 *MCVX_NH;
int 		 *MCVX_S;
int 		  MCVX_NS;
void MCVX_Print(){
  int s,n;
  if(MCVX_NS){
    for(s=0;s<MCVX_NS;s++){
      if(MCVX_NH[s]){
        fprintf(stdout,"Group:\tMCVX\tSector:\t%d \tNhits:\t%d\n",MCVX_S[s],MCVX_NH[s]);
        for(n=0;n<MCVX_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\ttof: %8.3f\tflag: %d\t",
		MCVX[s][n].x,
		MCVX[s][n].y,
		MCVX[s][n].z,
		MCVX[s][n].tof,
		MCVX[s][n].flag);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MS1 - multi sector bank
struct MS1_t	**MS1;
int 		 *MS1_NH;
int 		 *MS1_S;
int 		  MS1_NS;
void MS1_Print(){
  int s,n;
  if(MS1_NS){
    for(s=0;s<MS1_NS;s++){
      if(MS1_NH[s]){
        fprintf(stdout,"Group:\tMS1\tSector:\t%d \tNhits:\t%d\n",MS1_S[s],MS1_NH[s]);
        for(n=0;n<MS1_NH[s];n++){
          fprintf(stdout,"ID: %hd\tADC: %hd\t",
		MS1[s][n].ID,
		MS1[s][n].ADC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MTRK - multi sector bank
struct MTRK_t	**MTRK;
int 		 *MTRK_NH;
int 		 *MTRK_S;
int 		  MTRK_NS;
void MTRK_Print(){
  int s,n;
  if(MTRK_NS){
    for(s=0;s<MTRK_NS;s++){
      if(MTRK_NH[s]){
        fprintf(stdout,"Group:\tMTRK\tSector:\t%d \tNhits:\t%d\n",MTRK_S[s],MTRK_NH[s]);
        for(n=0;n<MTRK_NH[s];n++){
          fprintf(stdout,"sect: %d\ttrkl1: %d\ttrkl2: %d\ttrkl3: %d\ttrkl4: %d\ttrkl5: %d\ttrkl6: %d\t",
		MTRK[s][n].sect,
		MTRK[s][n].trkl1,
		MTRK[s][n].trkl2,
		MTRK[s][n].trkl3,
		MTRK[s][n].trkl4,
		MTRK[s][n].trkl5,
		MTRK[s][n].trkl6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank MVRT - single sector bank
struct MVRT_t	 *MVRT;
int 		  MVRT_NH;
int 		  MVRT_S;
void MVRT_Print(){
  int n;
  if(MVRT_NH){
    fprintf(stdout,"Group:\tMVRT\tSector:\t%d \tNhits:\t%d\n",MVRT_S,MVRT_NH);
    for(n=0;n<MVRT_NH;n++){
      fprintf(stdout,"v_id: %d\tntrk: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tchi2: %8.3f\tcxx: %8.3f\tcxy: %8.3f\tcxz: %8.3f\tcyy: %8.3f\tcyz: %8.3f\tczz: %8.3f\tstat: %d\t",
		MVRT[n].v_id,
		MVRT[n].ntrk,
		MVRT[n].x,
		MVRT[n].y,
		MVRT[n].z,
		MVRT[n].chi2,
		MVRT[n].cxx,
		MVRT[n].cxy,
		MVRT[n].cxz,
		MVRT[n].cyy,
		MVRT[n].cyz,
		MVRT[n].czz,
		MVRT[n].stat);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank PART - single sector bank
struct PART_t	 *PART;
int 		  PART_NH;
int 		  PART_S;
void PART_Print(){
  int n;
  if(PART_NH){
    fprintf(stdout,"Group:\tPART\tSector:\t%d \tNhits:\t%d\n",PART_S,PART_NH);
    for(n=0;n<PART_NH;n++){
      fprintf(stdout,"pid: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tE: %8.3f\tpx: %8.3f\tpy: %8.3f\tpz: %8.3f\tq: %8.3f\ttrkid: %d\tqpid: %8.3f\tqtrk: %8.3f\tflags: %d\t",
		PART[n].pid,
		PART[n].x,
		PART[n].y,
		PART[n].z,
		PART[n].E,
		PART[n].px,
		PART[n].py,
		PART[n].pz,
		PART[n].q,
		PART[n].trkid,
		PART[n].qpid,
		PART[n].qtrk,
		PART[n].flags);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank PCO - multi sector bank
struct PCO_t	**PCO;
int 		 *PCO_NH;
int 		 *PCO_S;
int 		  PCO_NS;
void PCO_Print(){
  int s,n;
  if(PCO_NS){
    for(s=0;s<PCO_NS;s++){
      if(PCO_NH[s]){
        fprintf(stdout,"Group:\tPCO\tSector:\t%d \tNhits:\t%d\n",PCO_S[s],PCO_NH[s]);
        for(n=0;n<PCO_NH[s];n++){
          fprintf(stdout,"TIME: %8.3f\tELT: %8.3f\tERB: %8.3f\tELB: %8.3f\tERT: %8.3f\tEMAIN: %8.3f\tEVETO: %8.3f\tTID: %d\t",
		PCO[s][n].TIME,
		PCO[s][n].ELT,
		PCO[s][n].ERB,
		PCO[s][n].ELB,
		PCO[s][n].ERT,
		PCO[s][n].EMAIN,
		PCO[s][n].EVETO,
		PCO[s][n].TID);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PHTM - multi sector bank
struct PHTM_t	**PHTM;
int 		 *PHTM_NH;
int 		 *PHTM_S;
int 		  PHTM_NS;
void PHTM_Print(){
  int s,n;
  if(PHTM_NS){
    for(s=0;s<PHTM_NS;s++){
      if(PHTM_NH[s]){
        fprintf(stdout,"Group:\tPHTM\tSector:\t%d \tNhits:\t%d\n",PHTM_S[s],PHTM_NH[s]);
        for(n=0;n<PHTM_NH[s];n++){
          fprintf(stdout,"EVID: %d\tNsprot: %d\tNpip: %d\tNpim: %d\tNppippim: %d\tNppip: %d\tNpippim: %d\tNppim: %d\tNkp: %d\tNpkp: %d\ttag_pi_v: %8.3f\t",
		PHTM[s][n].EVID,
		PHTM[s][n].Nsprot,
		PHTM[s][n].Npip,
		PHTM[s][n].Npim,
		PHTM[s][n].Nppippim,
		PHTM[s][n].Nppip,
		PHTM[s][n].Npippim,
		PHTM[s][n].Nppim,
		PHTM[s][n].Nkp,
		PHTM[s][n].Npkp,
		PHTM[s][n].tag_pi_v);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID1 - multi sector bank
struct PID1_t	**PID1;
int 		 *PID1_NH;
int 		 *PID1_S;
int 		  PID1_NS;
void PID1_Print(){
  int s,n;
  if(PID1_NS){
    for(s=0;s<PID1_NS;s++){
      if(PID1_NH[s]){
        fprintf(stdout,"Group:\tPID1\tSector:\t%d \tNhits:\t%d\n",PID1_S[s],PID1_NH[s]);
        for(n=0;n<PID1_NH[s];n++){
          fprintf(stdout,"NHBpos_1: %d\tNHBneg_1: %d\tNTBpos_1: %d\tNTBneg_1: %d\tchi2pos_1: %8.3f\tchi2neg_1: %8.3f\tNunknown_1: %d\tNgamma_1: %d\tNelec_1: %d\tNpiplus_1: %d\tNpiminus_1: %d\tNprot_1: %d\tNkplus_1: %d\tNkminus_1: %d\tNdeuteron_1: %d\tNneutron_1: %d\t",
		PID1[s][n].NHBpos_1,
		PID1[s][n].NHBneg_1,
		PID1[s][n].NTBpos_1,
		PID1[s][n].NTBneg_1,
		PID1[s][n].chi2pos_1,
		PID1[s][n].chi2neg_1,
		PID1[s][n].Nunknown_1,
		PID1[s][n].Ngamma_1,
		PID1[s][n].Nelec_1,
		PID1[s][n].Npiplus_1,
		PID1[s][n].Npiminus_1,
		PID1[s][n].Nprot_1,
		PID1[s][n].Nkplus_1,
		PID1[s][n].Nkminus_1,
		PID1[s][n].Ndeuteron_1,
		PID1[s][n].Nneutron_1);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID2 - multi sector bank
struct PID2_t	**PID2;
int 		 *PID2_NH;
int 		 *PID2_S;
int 		  PID2_NS;
void PID2_Print(){
  int s,n;
  if(PID2_NS){
    for(s=0;s<PID2_NS;s++){
      if(PID2_NH[s]){
        fprintf(stdout,"Group:\tPID2\tSector:\t%d \tNhits:\t%d\n",PID2_S[s],PID2_NH[s]);
        for(n=0;n<PID2_NH[s];n++){
          fprintf(stdout,"NHBpos_2: %d\tNHBneg_2: %d\tNTBpos_2: %d\tNTBneg_2: %d\tchi2pos_2: %8.3f\tchi2neg_2: %8.3f\tNunknown_2: %d\tNgamma_2: %d\tNelec_2: %d\tNpiplus_2: %d\tNpiminus_2: %d\tNprot_2: %d\tNkplus_2: %d\tNkminus_2: %d\tNdeuteron_2: %d\tNneutron_2: %d\t",
		PID2[s][n].NHBpos_2,
		PID2[s][n].NHBneg_2,
		PID2[s][n].NTBpos_2,
		PID2[s][n].NTBneg_2,
		PID2[s][n].chi2pos_2,
		PID2[s][n].chi2neg_2,
		PID2[s][n].Nunknown_2,
		PID2[s][n].Ngamma_2,
		PID2[s][n].Nelec_2,
		PID2[s][n].Npiplus_2,
		PID2[s][n].Npiminus_2,
		PID2[s][n].Nprot_2,
		PID2[s][n].Nkplus_2,
		PID2[s][n].Nkminus_2,
		PID2[s][n].Ndeuteron_2,
		PID2[s][n].Nneutron_2);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID3 - multi sector bank
struct PID3_t	**PID3;
int 		 *PID3_NH;
int 		 *PID3_S;
int 		  PID3_NS;
void PID3_Print(){
  int s,n;
  if(PID3_NS){
    for(s=0;s<PID3_NS;s++){
      if(PID3_NH[s]){
        fprintf(stdout,"Group:\tPID3\tSector:\t%d \tNhits:\t%d\n",PID3_S[s],PID3_NH[s]);
        for(n=0;n<PID3_NH[s];n++){
          fprintf(stdout,"NHBpos_3: %d\tNHBneg_3: %d\tNTBpos_3: %d\tNTBneg_3: %d\tchi2pos_3: %8.3f\tchi2neg_3: %8.3f\tNunknown_3: %d\tNgamma_3: %d\tNelec_3: %d\tNpiplus_3: %d\tNpiminus_3: %d\tNprot_3: %d\tNkplus_3: %d\tNkminus_3: %d\tNdeuteron_3: %d\tNneutron_3: %d\t",
		PID3[s][n].NHBpos_3,
		PID3[s][n].NHBneg_3,
		PID3[s][n].NTBpos_3,
		PID3[s][n].NTBneg_3,
		PID3[s][n].chi2pos_3,
		PID3[s][n].chi2neg_3,
		PID3[s][n].Nunknown_3,
		PID3[s][n].Ngamma_3,
		PID3[s][n].Nelec_3,
		PID3[s][n].Npiplus_3,
		PID3[s][n].Npiminus_3,
		PID3[s][n].Nprot_3,
		PID3[s][n].Nkplus_3,
		PID3[s][n].Nkminus_3,
		PID3[s][n].Ndeuteron_3,
		PID3[s][n].Nneutron_3);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID4 - multi sector bank
struct PID4_t	**PID4;
int 		 *PID4_NH;
int 		 *PID4_S;
int 		  PID4_NS;
void PID4_Print(){
  int s,n;
  if(PID4_NS){
    for(s=0;s<PID4_NS;s++){
      if(PID4_NH[s]){
        fprintf(stdout,"Group:\tPID4\tSector:\t%d \tNhits:\t%d\n",PID4_S[s],PID4_NH[s]);
        for(n=0;n<PID4_NH[s];n++){
          fprintf(stdout,"NHBpos_4: %d\tNHBneg_4: %d\tNTBpos_4: %d\tNTBneg_4: %d\tchi2pos_4: %8.3f\tchi2neg_4: %8.3f\tNunknown_4: %d\tNgamma_4: %d\tNelec_4: %d\tNpiplus_4: %d\tNpiminus_4: %d\tNprot_4: %d\tNkplus_4: %d\tNkminus_4: %d\tNdeuteron_4: %d\tNneutron_4: %d\t",
		PID4[s][n].NHBpos_4,
		PID4[s][n].NHBneg_4,
		PID4[s][n].NTBpos_4,
		PID4[s][n].NTBneg_4,
		PID4[s][n].chi2pos_4,
		PID4[s][n].chi2neg_4,
		PID4[s][n].Nunknown_4,
		PID4[s][n].Ngamma_4,
		PID4[s][n].Nelec_4,
		PID4[s][n].Npiplus_4,
		PID4[s][n].Npiminus_4,
		PID4[s][n].Nprot_4,
		PID4[s][n].Nkplus_4,
		PID4[s][n].Nkminus_4,
		PID4[s][n].Ndeuteron_4,
		PID4[s][n].Nneutron_4);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID5 - multi sector bank
struct PID5_t	**PID5;
int 		 *PID5_NH;
int 		 *PID5_S;
int 		  PID5_NS;
void PID5_Print(){
  int s,n;
  if(PID5_NS){
    for(s=0;s<PID5_NS;s++){
      if(PID5_NH[s]){
        fprintf(stdout,"Group:\tPID5\tSector:\t%d \tNhits:\t%d\n",PID5_S[s],PID5_NH[s]);
        for(n=0;n<PID5_NH[s];n++){
          fprintf(stdout,"NHBpos_5: %d\tNHBneg_5: %d\tNTBpos_5: %d\tNTBneg_5: %d\tchi2pos_5: %8.3f\tchi2neg_5: %8.3f\tNunknown_5: %d\tNgamma_5: %d\tNelec_5: %d\tNpiplus_5: %d\tNpiminus_5: %d\tNprot_5: %d\tNkplus_5: %d\tNkminus_5: %d\tNdeuteron_5: %d\tNneutron_5: %d\t",
		PID5[s][n].NHBpos_5,
		PID5[s][n].NHBneg_5,
		PID5[s][n].NTBpos_5,
		PID5[s][n].NTBneg_5,
		PID5[s][n].chi2pos_5,
		PID5[s][n].chi2neg_5,
		PID5[s][n].Nunknown_5,
		PID5[s][n].Ngamma_5,
		PID5[s][n].Nelec_5,
		PID5[s][n].Npiplus_5,
		PID5[s][n].Npiminus_5,
		PID5[s][n].Nprot_5,
		PID5[s][n].Nkplus_5,
		PID5[s][n].Nkminus_5,
		PID5[s][n].Ndeuteron_5,
		PID5[s][n].Nneutron_5);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PID6 - multi sector bank
struct PID6_t	**PID6;
int 		 *PID6_NH;
int 		 *PID6_S;
int 		  PID6_NS;
void PID6_Print(){
  int s,n;
  if(PID6_NS){
    for(s=0;s<PID6_NS;s++){
      if(PID6_NH[s]){
        fprintf(stdout,"Group:\tPID6\tSector:\t%d \tNhits:\t%d\n",PID6_S[s],PID6_NH[s]);
        for(n=0;n<PID6_NH[s];n++){
          fprintf(stdout,"NHBpos_6: %d\tNHBneg_6: %d\tNTBpos_6: %d\tNTBneg_6: %d\tchi2pos_6: %8.3f\tchi2neg_6: %8.3f\tNunknown_6: %d\tNgamma_6: %d\tNelec_6: %d\tNpiplus_6: %d\tNpiminus_6: %d\tNprot_6: %d\tNkplus_6: %d\tNkminus_6: %d\tNdeuteron_6: %d\tNneutron_6: %d\t",
		PID6[s][n].NHBpos_6,
		PID6[s][n].NHBneg_6,
		PID6[s][n].NTBpos_6,
		PID6[s][n].NTBneg_6,
		PID6[s][n].chi2pos_6,
		PID6[s][n].chi2neg_6,
		PID6[s][n].Nunknown_6,
		PID6[s][n].Ngamma_6,
		PID6[s][n].Nelec_6,
		PID6[s][n].Npiplus_6,
		PID6[s][n].Npiminus_6,
		PID6[s][n].Nprot_6,
		PID6[s][n].Nkplus_6,
		PID6[s][n].Nkminus_6,
		PID6[s][n].Ndeuteron_6,
		PID6[s][n].Nneutron_6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PIDT - multi sector bank
struct PIDT_t	**PIDT;
int 		 *PIDT_NH;
int 		 *PIDT_S;
int 		  PIDT_NS;
void PIDT_Print(){
  int s,n;
  if(PIDT_NS){
    for(s=0;s<PIDT_NS;s++){
      if(PIDT_NH[s]){
        fprintf(stdout,"Group:\tPIDT\tSector:\t%d \tNhits:\t%d\n",PIDT_S[s],PIDT_NH[s]);
        for(n=0;n<PIDT_NH[s];n++){
          fprintf(stdout,"Nunknown: %d\tNgamma: %d\tNelec: %d\tNpiplus: %d\tNpiminus: %d\tNprot: %d\tNkplus: %d\tNkminus: %d\tNdeuteron: %d\tNneutron: %d\t",
		PIDT[s][n].Nunknown,
		PIDT[s][n].Ngamma,
		PIDT[s][n].Nelec,
		PIDT[s][n].Npiplus,
		PIDT[s][n].Npiminus,
		PIDT[s][n].Nprot,
		PIDT[s][n].Nkplus,
		PIDT[s][n].Nkminus,
		PIDT[s][n].Ndeuteron,
		PIDT[s][n].Nneutron);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PRTM - multi sector bank
struct PRTM_t	**PRTM;
int 		 *PRTM_NH;
int 		 *PRTM_S;
int 		  PRTM_NS;
void PRTM_Print(){
  int s,n;
  if(PRTM_NS){
    for(s=0;s<PRTM_NS;s++){
      if(PRTM_NH[s]){
        fprintf(stdout,"Group:\tPRTM\tSector:\t%d \tNhits:\t%d\n",PRTM_S[s],PRTM_NH[s]);
        for(n=0;n<PRTM_NH[s];n++){
          fprintf(stdout,"time: %d\t",
		PRTM[s][n].time);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PSO - multi sector bank
struct PSO_t	**PSO;
int 		 *PSO_NH;
int 		 *PSO_S;
int 		  PSO_NS;
void PSO_Print(){
  int s,n;
  if(PSO_NS){
    for(s=0;s<PSO_NS;s++){
      if(PSO_NH[s]){
        fprintf(stdout,"Group:\tPSO\tSector:\t%d \tNhits:\t%d\n",PSO_S[s],PSO_NH[s]);
        for(n=0;n<PSO_NH[s];n++){
          fprintf(stdout,"ID: %d\tTIME: %8.3f\tENER: %8.3f\tTID: %d\t",
		PSO[s][n].ID,
		PSO[s][n].TIME,
		PSO[s][n].ENER,
		PSO[s][n].TID);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank PTDB - multi sector bank
struct PTDB_t	**PTDB;
int 		 *PTDB_NH;
int 		 *PTDB_S;
int 		  PTDB_NS;
void PTDB_Print(){
  int s,n;
  if(PTDB_NS){
    for(s=0;s<PTDB_NS;s++){
      if(PTDB_NH[s]){
        fprintf(stdout,"Group:\tPTDB\tSector:\t%d \tNhits:\t%d\n",PTDB_S[s],PTDB_NH[s]);
        for(n=0;n<PTDB_NH[s];n++){
          fprintf(stdout,"PBeam: %d\tPB_DATE: %d\tTTYPE: %d\tB_TARG: %d\tPTarg: %d\tPT_time: %d\tHE_Level: %d\tEIOF: %d\tTTempH: %d\tTTempC: %d\tAnealT: %d\t",
		PTDB[s][n].PBeam,
		PTDB[s][n].PB_DATE,
		PTDB[s][n].TTYPE,
		PTDB[s][n].B_TARG,
		PTDB[s][n].PTarg,
		PTDB[s][n].PT_time,
		PTDB[s][n].HE_Level,
		PTDB[s][n].EIOF,
		PTDB[s][n].TTempH,
		PTDB[s][n].TTempC,
		PTDB[s][n].AnealT);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RCST - multi sector bank
struct RCST_t	**RCST;
int 		 *RCST_NH;
int 		 *RCST_S;
int 		  RCST_NS;
void RCST_Print(){
  int s,n;
  if(RCST_NS){
    for(s=0;s<RCST_NS;s++){
      if(RCST_NH[s]){
        fprintf(stdout,"Group:\tRCST\tSector:\t%d \tNhits:\t%d\n",RCST_S[s],RCST_NH[s]);
        for(n=0;n<RCST_NH[s];n++){
          fprintf(stdout,"down_count: %d\tprestart_count: %d\tgo_count: %d\ttrig_count: %d\tevent_count: %d\tsync_count: %d\trun_trig_count: %d\trun_event_count: %d\trun_sync_count: %d\tpause_count: %d\tend_count: %d\tillegal_count: %d\trun_illegal_count: %d\tphys_sync_count: %d\trun_phys_sync_count: %d\t",
		RCST[s][n].down_count,
		RCST[s][n].prestart_count,
		RCST[s][n].go_count,
		RCST[s][n].trig_count,
		RCST[s][n].event_count,
		RCST[s][n].sync_count,
		RCST[s][n].run_trig_count,
		RCST[s][n].run_event_count,
		RCST[s][n].run_sync_count,
		RCST[s][n].pause_count,
		RCST[s][n].end_count,
		RCST[s][n].illegal_count,
		RCST[s][n].run_illegal_count,
		RCST[s][n].phys_sync_count,
		RCST[s][n].run_phys_sync_count);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank REF - single sector bank
struct REF_t	 *REF;
int 		  REF_NH;
int 		  REF_S;
void REF_Print(){
  int n;
  if(REF_NH){
    fprintf(stdout,"Group:\tREF\tSector:\t%d \tNhits:\t%d\n",REF_S,REF_NH);
    for(n=0;n<REF_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		REF[n].ID,
		REF[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank RF - multi sector bank
struct RF_t	**RF;
int 		 *RF_NH;
int 		 *RF_S;
int 		  RF_NS;
void RF_Print(){
  int s,n;
  if(RF_NS){
    for(s=0;s<RF_NS;s++){
      if(RF_NH[s]){
        fprintf(stdout,"Group:\tRF\tSector:\t%d \tNhits:\t%d\n",RF_S[s],RF_NH[s]);
        for(n=0;n<RF_NH[s];n++){
          fprintf(stdout,"RF: %8.3f\tRF1: %8.3f\tRF2: %8.3f\t",
		RF[s][n].RF,
		RF[s][n].RF1,
		RF[s][n].RF2);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RFT - single sector bank
struct RFT_t	 *RFT;
int 		  RFT_NH;
int 		  RFT_S;
void RFT_Print(){
  int n;
  if(RFT_NH){
    fprintf(stdout,"Group:\tRFT\tSector:\t%d \tNhits:\t%d\n",RFT_S,RFT_NH);
    for(n=0;n<RFT_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		RFT[n].ID,
		RFT[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank RGLK - multi sector bank
struct RGLK_t	**RGLK;
int 		 *RGLK_NH;
int 		 *RGLK_S;
int 		  RGLK_NS;
void RGLK_Print(){
  int s,n;
  if(RGLK_NS){
    for(s=0;s<RGLK_NS;s++){
      if(RGLK_NH[s]){
        fprintf(stdout,"Group:\tRGLK\tSector:\t%d \tNhits:\t%d\n",RGLK_S[s],RGLK_NH[s]);
        for(n=0;n<RGLK_NH[s];n++){
          fprintf(stdout,"IREGION: %d\tX: %8.3f\tY: %8.3f\tZ: %8.3f\tTHETA0: %8.3f\tPHI0: %8.3f\tRTHETA: %8.3f\tRPHI: %8.3f\tCHI2: %8.3f\tSTATUS: %d\t",
		RGLK[s][n].IREGION,
		RGLK[s][n].X,
		RGLK[s][n].Y,
		RGLK[s][n].Z,
		RGLK[s][n].THETA0,
		RGLK[s][n].PHI0,
		RGLK[s][n].RTHETA,
		RGLK[s][n].RPHI,
		RGLK[s][n].CHI2,
		RGLK[s][n].STATUS);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RNLG - multi sector bank
struct RNLG_t	**RNLG;
int 		 *RNLG_NH;
int 		 *RNLG_S;
int 		  RNLG_NS;
void RNLG_Print(){
  int s,n;
  if(RNLG_NS){
    for(s=0;s<RNLG_NS;s++){
      if(RNLG_NH[s]){
        fprintf(stdout,"Group:\tRNLG\tSector:\t%d \tNhits:\t%d\n",RNLG_S[s],RNLG_NH[s]);
        for(n=0;n<RNLG_NH[s];n++){
          fprintf(stdout,"sql: %.4s\t",
		RNLG[s][n].sql);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RNPE - multi sector bank
struct RNPE_t	**RNPE;
int 		 *RNPE_NH;
int 		 *RNPE_S;
int 		  RNPE_NS;
void RNPE_Print(){
  int s,n;
  if(RNPE_NS){
    for(s=0;s<RNPE_NS;s++){
      if(RNPE_NH[s]){
        fprintf(stdout,"Group:\tRNPE\tSector:\t%d \tNhits:\t%d\n",RNPE_S[s],RNPE_NH[s]);
        for(n=0;n<RNPE_NH[s];n++){
          fprintf(stdout,"run: %d\tbank: %.4s\t",
		RNPE[s][n].run,
		RNPE[s][n].bank);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RTSL - multi sector bank
struct RTSL_t	**RTSL;
int 		 *RTSL_NH;
int 		 *RTSL_S;
int 		  RTSL_NS;
void RTSL_Print(){
  int s,n;
  if(RTSL_NS){
    for(s=0;s<RTSL_NS;s++){
      if(RTSL_NH[s]){
        fprintf(stdout,"Group:\tRTSL\tSector:\t%d \tNhits:\t%d\n",RTSL_S[s],RTSL_NH[s]);
        for(n=0;n<RTSL_NH[s];n++){
          fprintf(stdout,"RAWT1: %d\tRAWT2: %d\tRAWT3: %d\tRAWT4: %d\tRAWT5: %d\tRAWT6: %d\tRAWT7: %d\tRAWT8: %d\tRAWT9: %d\tRAWT10: %d\tRAWT11: %d\tRAWT12: %d\tRAWT13: %d\tRAWT14: %d\tRAWT15: %d\tRAWT16: %d\tRAWT17: %d\tRAWT18: %d\tRAWT19: %d\tRAWT20: %d\tRAWT21: %d\tRAWT22: %d\tRAWT23: %d\tRAWT24: %d\tRAWT25: %d\tRAWT26: %d\tRAWT27: %d\tRAWT28: %d\tRAWT29: %d\tRAWT30: %d\tRAWT31: %d\tRAWT32: %d\tRAWT33: %d\tRAWT34: %d\tRAWT35: %d\tRAWT36: %d\tRAWT37: %d\tRAWT38: %d\tRAWT39: %d\tRAWT40: %d\tRAWT41: %d\tRAWT42: %d\tRAWT43: %d\tRAWT44: %d\tRAWT45: %d\tRAWT46: %d\tRAWT47: %d\tRAWT48: %d\tRAWT49: %d\tRAWT50: %d\tRAWT51: %d\tRAWT52: %d\tRAWT53: %d\tRAWT54: %d\tRAWT55: %d\tRAWT56: %d\tRAWT57: %d\tRAWT58: %d\tRAWT59: %d\tRAWT60: %d\tRAWT61: %d\tRAWT62: %d\tRAWT63: %d\tRAWT64: %d\t",
		RTSL[s][n].RAWT1,
		RTSL[s][n].RAWT2,
		RTSL[s][n].RAWT3,
		RTSL[s][n].RAWT4,
		RTSL[s][n].RAWT5,
		RTSL[s][n].RAWT6,
		RTSL[s][n].RAWT7,
		RTSL[s][n].RAWT8,
		RTSL[s][n].RAWT9,
		RTSL[s][n].RAWT10,
		RTSL[s][n].RAWT11,
		RTSL[s][n].RAWT12,
		RTSL[s][n].RAWT13,
		RTSL[s][n].RAWT14,
		RTSL[s][n].RAWT15,
		RTSL[s][n].RAWT16,
		RTSL[s][n].RAWT17,
		RTSL[s][n].RAWT18,
		RTSL[s][n].RAWT19,
		RTSL[s][n].RAWT20,
		RTSL[s][n].RAWT21,
		RTSL[s][n].RAWT22,
		RTSL[s][n].RAWT23,
		RTSL[s][n].RAWT24,
		RTSL[s][n].RAWT25,
		RTSL[s][n].RAWT26,
		RTSL[s][n].RAWT27,
		RTSL[s][n].RAWT28,
		RTSL[s][n].RAWT29,
		RTSL[s][n].RAWT30,
		RTSL[s][n].RAWT31,
		RTSL[s][n].RAWT32,
		RTSL[s][n].RAWT33,
		RTSL[s][n].RAWT34,
		RTSL[s][n].RAWT35,
		RTSL[s][n].RAWT36,
		RTSL[s][n].RAWT37,
		RTSL[s][n].RAWT38,
		RTSL[s][n].RAWT39,
		RTSL[s][n].RAWT40,
		RTSL[s][n].RAWT41,
		RTSL[s][n].RAWT42,
		RTSL[s][n].RAWT43,
		RTSL[s][n].RAWT44,
		RTSL[s][n].RAWT45,
		RTSL[s][n].RAWT46,
		RTSL[s][n].RAWT47,
		RTSL[s][n].RAWT48,
		RTSL[s][n].RAWT49,
		RTSL[s][n].RAWT50,
		RTSL[s][n].RAWT51,
		RTSL[s][n].RAWT52,
		RTSL[s][n].RAWT53,
		RTSL[s][n].RAWT54,
		RTSL[s][n].RAWT55,
		RTSL[s][n].RAWT56,
		RTSL[s][n].RAWT57,
		RTSL[s][n].RAWT58,
		RTSL[s][n].RAWT59,
		RTSL[s][n].RAWT60,
		RTSL[s][n].RAWT61,
		RTSL[s][n].RAWT62,
		RTSL[s][n].RAWT63,
		RTSL[s][n].RAWT64);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RUNC - multi sector bank
struct RUNC_t	**RUNC;
int 		 *RUNC_NH;
int 		 *RUNC_S;
int 		  RUNC_NS;
void RUNC_Print(){
  int s,n;
  if(RUNC_NS){
    for(s=0;s<RUNC_NS;s++){
      if(RUNC_NH[s]){
        fprintf(stdout,"Group:\tRUNC\tSector:\t%d \tNhits:\t%d\n",RUNC_S[s],RUNC_NH[s]);
        for(n=0;n<RUNC_NH[s];n++){
          fprintf(stdout,"runno: %d\tbeam_e: %8.3f\tb_first: %d\tq_live: %8.3f\tql_first: %d\tq_all: %8.3f\tqa_first: %d\ttor_curr: %8.3f\tt_first: %d\tmtor_curr: %8.3f\tmt_first: %d\ttag_curr: %8.3f\ttag_first: %d\t",
		RUNC[s][n].runno,
		RUNC[s][n].beam_e,
		RUNC[s][n].b_first,
		RUNC[s][n].q_live,
		RUNC[s][n].ql_first,
		RUNC[s][n].q_all,
		RUNC[s][n].qa_first,
		RUNC[s][n].tor_curr,
		RUNC[s][n].t_first,
		RUNC[s][n].mtor_curr,
		RUNC[s][n].mt_first,
		RUNC[s][n].tag_curr,
		RUNC[s][n].tag_first);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank S1ST - multi sector bank
struct S1ST_t	**S1ST;
int 		 *S1ST_NH;
int 		 *S1ST_S;
int 		  S1ST_NS;
void S1ST_Print(){
  int s,n;
  if(S1ST_NS){
    for(s=0;s<S1ST_NS;s++){
      if(S1ST_NH[s]){
        fprintf(stdout,"Group:\tS1ST\tSector:\t%d \tNhits:\t%d\n",S1ST_S[s],S1ST_NH[s]);
        for(n=0;n<S1ST_NH[s];n++){
          fprintf(stdout,"latch1_bit1_count: %d\tlatch1_bit2_count: %d\tlatch1_bit3_count: %d\tlatch1_bit4_count: %d\tlatch1_bit5_count: %d\tlatch1_bit6_count: %d\tlatch1_bit7_count: %d\tlatch1_bit8_count: %d\tlatch1_bit9_count: %d\tlatch1_bit10_count: %d\tlatch1_bit11_count: %d\tlatch1_bit12_count: %d\tevent_count: %d\tlatch1_zero_count: %d\tlatch1_empty_count: %d\tlatch1_not_empty_count: %d\tlatch1_ok_count: %d\tlevel2_sector1: %d\tlevel2_sector2: %d\tlevel2_sector3: %d\tlevel2_sector4: %d\tlevel2_sector5: %d\tlevel2_sector6: %d\tlevel2_pass: %d\tlevel2_fail: %d\tlatch2_zero_count: %d\tlatch2_empty_count: %d\tlatch2_not_empty_count: %d\tlatch2_ok_count: %d\troc_13_count: %d\troc_15_count: %d\tl1l2_zero_count: %d\tl1zero_13_count: %d\tl2zero_13_count: %d\tl1l2zero_13_count: %d\t",
		S1ST[s][n].latch1_bit1_count,
		S1ST[s][n].latch1_bit2_count,
		S1ST[s][n].latch1_bit3_count,
		S1ST[s][n].latch1_bit4_count,
		S1ST[s][n].latch1_bit5_count,
		S1ST[s][n].latch1_bit6_count,
		S1ST[s][n].latch1_bit7_count,
		S1ST[s][n].latch1_bit8_count,
		S1ST[s][n].latch1_bit9_count,
		S1ST[s][n].latch1_bit10_count,
		S1ST[s][n].latch1_bit11_count,
		S1ST[s][n].latch1_bit12_count,
		S1ST[s][n].event_count,
		S1ST[s][n].latch1_zero_count,
		S1ST[s][n].latch1_empty_count,
		S1ST[s][n].latch1_not_empty_count,
		S1ST[s][n].latch1_ok_count,
		S1ST[s][n].level2_sector1,
		S1ST[s][n].level2_sector2,
		S1ST[s][n].level2_sector3,
		S1ST[s][n].level2_sector4,
		S1ST[s][n].level2_sector5,
		S1ST[s][n].level2_sector6,
		S1ST[s][n].level2_pass,
		S1ST[s][n].level2_fail,
		S1ST[s][n].latch2_zero_count,
		S1ST[s][n].latch2_empty_count,
		S1ST[s][n].latch2_not_empty_count,
		S1ST[s][n].latch2_ok_count,
		S1ST[s][n].roc_13_count,
		S1ST[s][n].roc_15_count,
		S1ST[s][n].l1l2_zero_count,
		S1ST[s][n].l1zero_13_count,
		S1ST[s][n].l2zero_13_count,
		S1ST[s][n].l1l2zero_13_count);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SC1 - multi sector bank
struct SC1_t	**SC1;
int 		 *SC1_NH;
int 		 *SC1_S;
int 		  SC1_NS;
void SC1_Print(){
  int s,n;
  if(SC1_NS){
    for(s=0;s<SC1_NS;s++){
      if(SC1_NH[s]){
        fprintf(stdout,"Group:\tSC1\tSector:\t%d \tNhits:\t%d\n",SC1_S[s],SC1_NH[s]);
        for(n=0;n<SC1_NH[s];n++){
          fprintf(stdout,"ID: %d\ttime_l: %8.3f\tenergy_l: %8.3f\ttime_r: %8.3f\tenergy_r: %8.3f\tdtime_l: %8.3f\tdenergy_l: %8.3f\tdtime_r: %8.3f\tdenergy_r: %8.3f\t",
		SC1[s][n].ID,
		SC1[s][n].time_l,
		SC1[s][n].energy_l,
		SC1[s][n].time_r,
		SC1[s][n].energy_r,
		SC1[s][n].dtime_l,
		SC1[s][n].denergy_l,
		SC1[s][n].dtime_r,
		SC1[s][n].denergy_r);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCC - multi sector bank
struct SCC_t	**SCC;
int 		 *SCC_NH;
int 		 *SCC_S;
int 		  SCC_NS;
void SCC_Print(){
  int s,n;
  if(SCC_NS){
    for(s=0;s<SCC_NS;s++){
      if(SCC_NH[s]){
        fprintf(stdout,"Group:\tSCC\tSector:\t%d \tNhits:\t%d\n",SCC_S[s],SCC_NH[s]);
        for(n=0;n<SCC_NH[s];n++){
          fprintf(stdout,"id: %d\tdate: %d\tversion: %d\tstatus: %d\tTD0L: %8.3f\tTD0Lu: %8.3f\tTD0R: %8.3f\tTD0Ru: %8.3f\tTD1L: %8.3f\tTD1Lu: %8.3f\tTD1R: %8.3f\tTD1Ru: %8.3f\tTD2L: %8.3f\tTD2Lu: %8.3f\tTD2R: %8.3f\tTD2Ru: %8.3f\tTW0L: %8.3f\tTW0Lu: %8.3f\tTW0R: %8.3f\tTW0Ru: %8.3f\tTW1L: %8.3f\tTW1Lu: %8.3f\tTW1R: %8.3f\tTW1Ru: %8.3f\tTW2L: %8.3f\tTW2Lu: %8.3f\tTW2R: %8.3f\tTW2Ru: %8.3f\tADPL: %8.3f\tADPLu: %8.3f\tADPR: %8.3f\tADPRu: %8.3f\tM0L: %8.3f\tM0Lu: %8.3f\tM0R: %8.3f\tM0Ru: %8.3f\tVEFL: %8.3f\tVEFLu: %8.3f\tVEFR: %8.3f\tVEFRu: %8.3f\tATNL: %8.3f\tATNLu: %8.3f\tATNR: %8.3f\tATNRu: %8.3f\tTSIG0L: %8.3f\tTSIG0R: %8.3f\tTSIG1L: %8.3f\tTSIG1R: %8.3f\t",
		SCC[s][n].id,
		SCC[s][n].date,
		SCC[s][n].version,
		SCC[s][n].status,
		SCC[s][n].TD0L,
		SCC[s][n].TD0Lu,
		SCC[s][n].TD0R,
		SCC[s][n].TD0Ru,
		SCC[s][n].TD1L,
		SCC[s][n].TD1Lu,
		SCC[s][n].TD1R,
		SCC[s][n].TD1Ru,
		SCC[s][n].TD2L,
		SCC[s][n].TD2Lu,
		SCC[s][n].TD2R,
		SCC[s][n].TD2Ru,
		SCC[s][n].TW0L,
		SCC[s][n].TW0Lu,
		SCC[s][n].TW0R,
		SCC[s][n].TW0Ru,
		SCC[s][n].TW1L,
		SCC[s][n].TW1Lu,
		SCC[s][n].TW1R,
		SCC[s][n].TW1Ru,
		SCC[s][n].TW2L,
		SCC[s][n].TW2Lu,
		SCC[s][n].TW2R,
		SCC[s][n].TW2Ru,
		SCC[s][n].ADPL,
		SCC[s][n].ADPLu,
		SCC[s][n].ADPR,
		SCC[s][n].ADPRu,
		SCC[s][n].M0L,
		SCC[s][n].M0Lu,
		SCC[s][n].M0R,
		SCC[s][n].M0Ru,
		SCC[s][n].VEFL,
		SCC[s][n].VEFLu,
		SCC[s][n].VEFR,
		SCC[s][n].VEFRu,
		SCC[s][n].ATNL,
		SCC[s][n].ATNLu,
		SCC[s][n].ATNR,
		SCC[s][n].ATNRu,
		SCC[s][n].TSIG0L,
		SCC[s][n].TSIG0R,
		SCC[s][n].TSIG1L,
		SCC[s][n].TSIG1R);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SC - multi sector bank
struct SC_t	**SC;
int 		 *SC_NH;
int 		 *SC_S;
int 		  SC_NS;
void SC_Print(){
  int s,n;
  if(SC_NS){
    for(s=0;s<SC_NS;s++){
      if(SC_NH[s]){
        fprintf(stdout,"Group:\tSC\tSector:\t%d \tNhits:\t%d\n",SC_S[s],SC_NH[s]);
        for(n=0;n<SC_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDCL: %hd\tADCL: %hd\tTDCR: %hd\tADCR: %hd\t",
		SC[s][n].ID,
		SC[s][n].TDCL,
		SC[s][n].ADCL,
		SC[s][n].TDCR,
		SC[s][n].ADCR);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCDI - multi sector bank
struct SCDI_t	**SCDI;
int 		 *SCDI_NH;
int 		 *SCDI_S;
int 		  SCDI_NS;
void SCDI_Print(){
  int s,n;
  if(SCDI_NS){
    for(s=0;s<SCDI_NS;s++){
      if(SCDI_NH[s]){
        fprintf(stdout,"Group:\tSCDI\tSector:\t%d \tNhits:\t%d\n",SCDI_S[s],SCDI_NH[s]);
        for(n=0;n<SCDI_NH[s];n++){
          fprintf(stdout,"crate: %d\tslot: %d\tmask: %d\tthreshold: %d\twidth: %d\t",
		SCDI[s][n].crate,
		SCDI[s][n].slot,
		SCDI[s][n].mask,
		SCDI[s][n].threshold,
		SCDI[s][n].width);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCGD - multi sector bank
struct SCGD_t	**SCGD;
int 		 *SCGD_NH;
int 		 *SCGD_S;
int 		  SCGD_NS;
void SCGD_Print(){
  int s,n;
  if(SCGD_NS){
    for(s=0;s<SCGD_NS;s++){
      if(SCGD_NH[s]){
        fprintf(stdout,"Group:\tSCGD\tSector:\t%d \tNhits:\t%d\n",SCGD_S[s],SCGD_NH[s]);
        for(n=0;n<SCGD_NH[s];n++){
          fprintf(stdout,"ID: %d\tNORM_Z: %8.3f\tNORM_X: %8.3f\tNORM_D: %8.3f\tALON_Z: %8.3f\tALON_X: %8.3f\tBEG_PD: %8.3f\tEND_PD: %8.3f\tPDL_SH: %8.3f\tR_BEAM: %8.3f\t",
		SCGD[s][n].ID,
		SCGD[s][n].NORM_Z,
		SCGD[s][n].NORM_X,
		SCGD[s][n].NORM_D,
		SCGD[s][n].ALON_Z,
		SCGD[s][n].ALON_X,
		SCGD[s][n].BEG_PD,
		SCGD[s][n].END_PD,
		SCGD[s][n].PDL_SH,
		SCGD[s][n].R_BEAM);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCG - multi sector bank
struct SCG_t	**SCG;
int 		 *SCG_NH;
int 		 *SCG_S;
int 		  SCG_NS;
void SCG_Print(){
  int s,n;
  if(SCG_NS){
    for(s=0;s<SCG_NS;s++){
      if(SCG_NH[s]){
        fprintf(stdout,"Group:\tSCG\tSector:\t%d \tNhits:\t%d\n",SCG_S[s],SCG_NH[s]);
        for(n=0;n<SCG_NH[s];n++){
          fprintf(stdout,"id: %d\tpanel: %d\tXccw: %8.3f\tYccw: %8.3f\tZccw: %8.3f\tXcw: %8.3f\tYcw: %8.3f\tZcw: %8.3f\tWIDTH: %8.3f\tTHICK: %8.3f\tDELTA: %8.3f\t",
		SCG[s][n].id,
		SCG[s][n].panel,
		SCG[s][n].Xccw,
		SCG[s][n].Yccw,
		SCG[s][n].Zccw,
		SCG[s][n].Xcw,
		SCG[s][n].Ycw,
		SCG[s][n].Zcw,
		SCG[s][n].WIDTH,
		SCG[s][n].THICK,
		SCG[s][n].DELTA);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCH - multi sector bank
struct SCH_t	**SCH;
int 		 *SCH_NH;
int 		 *SCH_S;
int 		  SCH_NS;
void SCH_Print(){
  int s,n;
  if(SCH_NS){
    for(s=0;s<SCH_NS;s++){
      if(SCH_NH[s]){
        fprintf(stdout,"Group:\tSCH\tSector:\t%d \tNhits:\t%d\n",SCH_S[s],SCH_NH[s]);
        for(n=0;n<SCH_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\tpmom: %8.3f\ttrack: %d\tid: %d\t",
		SCH[s][n].x,
		SCH[s][n].y,
		SCH[s][n].z,
		SCH[s][n].cx,
		SCH[s][n].cy,
		SCH[s][n].cz,
		SCH[s][n].pmom,
		SCH[s][n].track,
		SCH[s][n].id);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCMD - multi sector bank
struct SCMD_t	**SCMD;
int 		 *SCMD_NH;
int 		 *SCMD_S;
int 		  SCMD_NS;
void SCMD_Print(){
  int s,n;
  if(SCMD_NS){
    for(s=0;s<SCMD_NS;s++){
      if(SCMD_NH[s]){
        fprintf(stdout,"Group:\tSCMD\tSector:\t%d \tNhits:\t%d\n",SCMD_S[s],SCMD_NH[s]);
        for(n=0;n<SCMD_NH[s];n++){
          fprintf(stdout,"ID: %d\tSTATUS: %d\tX_NORM: %8.3f\tY_NORM: %8.3f\tD_NORM: %8.3f\tX_PANEL: %8.3f\tY_PANEL: %8.3f\tB_PAD: %8.3f\tE_PAD: %8.3f\tTOF: %8.3f\tEDP: %8.3f\tP_TIME: %8.3f\tP_ATLN: %8.3f\t",
		SCMD[s][n].ID,
		SCMD[s][n].STATUS,
		SCMD[s][n].X_NORM,
		SCMD[s][n].Y_NORM,
		SCMD[s][n].D_NORM,
		SCMD[s][n].X_PANEL,
		SCMD[s][n].Y_PANEL,
		SCMD[s][n].B_PAD,
		SCMD[s][n].E_PAD,
		SCMD[s][n].TOF,
		SCMD[s][n].EDP,
		SCMD[s][n].P_TIME,
		SCMD[s][n].P_ATLN);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCMT - multi sector bank
struct SCMT_t	**SCMT;
int 		 *SCMT_NH;
int 		 *SCMT_S;
int 		  SCMT_NS;
void SCMT_Print(){
  int s,n;
  if(SCMT_NS){
    for(s=0;s<SCMT_NS;s++){
      if(SCMT_NH[s]){
        fprintf(stdout,"Group:\tSCMT\tSector:\t%d \tNhits:\t%d\n",SCMT_S[s],SCMT_NH[s]);
        for(n=0;n<SCMT_NH[s];n++){
          fprintf(stdout,"mean_thr: %d\t",
		SCMT[s][n].mean_thr);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCMW - multi sector bank
struct SCMW_t	**SCMW;
int 		 *SCMW_NH;
int 		 *SCMW_S;
int 		  SCMW_NS;
void SCMW_Print(){
  int s,n;
  if(SCMW_NS){
    for(s=0;s<SCMW_NS;s++){
      if(SCMW_NH[s]){
        fprintf(stdout,"Group:\tSCMW\tSector:\t%d \tNhits:\t%d\n",SCMW_S[s],SCMW_NH[s]);
        for(n=0;n<SCMW_NH[s];n++){
          fprintf(stdout,"mean_width: %d\t",
		SCMW[s][n].mean_width);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCPB - single sector bank
struct SCPB_t	 *SCPB;
int 		  SCPB_NH;
int 		  SCPB_S;
void SCPB_Print(){
  int n;
  if(SCPB_NH){
    fprintf(stdout,"Group:\tSCPB\tSector:\t%d \tNhits:\t%d\n",SCPB_S,SCPB_NH);
    for(n=0;n<SCPB_NH;n++){
      fprintf(stdout,"ScPdHt: %d\tEdep: %8.3f\tTime: %8.3f\tPath: %8.3f\tChi2SC: %8.3f\tStatus: %d\t",
		SCPB[n].ScPdHt,
		SCPB[n].Edep,
		SCPB[n].Time,
		SCPB[n].Path,
		SCPB[n].Chi2SC,
		SCPB[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank SCP - multi sector bank
struct SCP_t	**SCP;
int 		 *SCP_NH;
int 		 *SCP_S;
int 		  SCP_NS;
void SCP_Print(){
  int s,n;
  if(SCP_NS){
    for(s=0;s<SCP_NS;s++){
      if(SCP_NH[s]){
        fprintf(stdout,"Group:\tSCP\tSector:\t%d \tNhits:\t%d\n",SCP_S[s],SCP_NH[s]);
        for(n=0;n<SCP_NH[s];n++){
          fprintf(stdout,"N1x: %8.3f\tN1y: %8.3f\tN1z: %8.3f\tR1n: %8.3f\tN2x: %8.3f\tN2y: %8.3f\tN2z: %8.3f\tR2n: %8.3f\tN3x: %8.3f\tN3y: %8.3f\tN3z: %8.3f\tR3n: %8.3f\tN4x: %8.3f\tN4y: %8.3f\tN4z: %8.3f\tR4n: %8.3f\t",
		SCP[s][n].N1x,
		SCP[s][n].N1y,
		SCP[s][n].N1z,
		SCP[s][n].R1n,
		SCP[s][n].N2x,
		SCP[s][n].N2y,
		SCP[s][n].N2z,
		SCP[s][n].R2n,
		SCP[s][n].N3x,
		SCP[s][n].N3y,
		SCP[s][n].N3z,
		SCP[s][n].R3n,
		SCP[s][n].N4x,
		SCP[s][n].N4y,
		SCP[s][n].N4z,
		SCP[s][n].R4n);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCPE - multi sector bank
struct SCPE_t	**SCPE;
int 		 *SCPE_NH;
int 		 *SCPE_S;
int 		  SCPE_NS;
void SCPE_Print(){
  int s,n;
  if(SCPE_NS){
    for(s=0;s<SCPE_NS;s++){
      if(SCPE_NH[s]){
        fprintf(stdout,"Group:\tSCPE\tSector:\t%d \tNhits:\t%d\n",SCPE_S[s],SCPE_NH[s]);
        for(n=0;n<SCPE_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean_left: %d\tsigma_left: %8.3f\tmean_right: %d\tsigma_right: %8.3f\t",
		SCPE[s][n].ID,
		SCPE[s][n].mean_left,
		SCPE[s][n].sigma_left,
		SCPE[s][n].mean_right,
		SCPE[s][n].sigma_right);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCPS - multi sector bank
struct SCPS_t	**SCPS;
int 		 *SCPS_NH;
int 		 *SCPS_S;
int 		  SCPS_NS;
void SCPS_Print(){
  int s,n;
  if(SCPS_NS){
    for(s=0;s<SCPS_NS;s++){
      if(SCPS_NH[s]){
        fprintf(stdout,"Group:\tSCPS\tSector:\t%d \tNhits:\t%d\n",SCPS_S[s],SCPS_NH[s]);
        for(n=0;n<SCPS_NH[s];n++){
          fprintf(stdout,"ID: %d\tFIRSTPNUM: %d\tLASTPNUM: %d\tDISTANCE: %8.3f\tBEGALONSHFT: %8.3f\tENDALONSHFT: %8.3f\tLENGTHSHFT: %8.3f\tALONGIDX: %8.3f\tALONGIDY: %8.3f\tALONGIDZ: %8.3f\tSLENGTHX: %8.3f\tSLENGTHY: %8.3f\tSLENGTHZ: %8.3f\tUNORMALX: %8.3f\tUNORMALY: %8.3f\tUNORMALZ: %8.3f\t",
		SCPS[s][n].ID,
		SCPS[s][n].FIRSTPNUM,
		SCPS[s][n].LASTPNUM,
		SCPS[s][n].DISTANCE,
		SCPS[s][n].BEGALONSHFT,
		SCPS[s][n].ENDALONSHFT,
		SCPS[s][n].LENGTHSHFT,
		SCPS[s][n].ALONGIDX,
		SCPS[s][n].ALONGIDY,
		SCPS[s][n].ALONGIDZ,
		SCPS[s][n].SLENGTHX,
		SCPS[s][n].SLENGTHY,
		SCPS[s][n].SLENGTHZ,
		SCPS[s][n].UNORMALX,
		SCPS[s][n].UNORMALY,
		SCPS[s][n].UNORMALZ);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCRC - multi sector bank
struct SCRC_t	**SCRC;
int 		 *SCRC_NH;
int 		 *SCRC_S;
int 		  SCRC_NS;
void SCRC_Print(){
  int s,n;
  if(SCRC_NS){
    for(s=0;s<SCRC_NS;s++){
      if(SCRC_NH[s]){
        fprintf(stdout,"Group:\tSCRC\tSector:\t%d \tNhits:\t%d\n",SCRC_S[s],SCRC_NH[s]);
        for(n=0;n<SCRC_NH[s];n++){
          fprintf(stdout,"id: %d\tenergy: %8.3f\tdenergy: %8.3f\ttime: %8.3f\tdtime: %8.3f\tx: %8.3f\ty: %8.3f\tz: %8.3f\tdx: %8.3f\tdy: %8.3f\tdz: %8.3f\tstatus: %d\t",
		SCRC[s][n].id,
		SCRC[s][n].energy,
		SCRC[s][n].denergy,
		SCRC[s][n].time,
		SCRC[s][n].dtime,
		SCRC[s][n].x,
		SCRC[s][n].y,
		SCRC[s][n].z,
		SCRC[s][n].dx,
		SCRC[s][n].dy,
		SCRC[s][n].dz,
		SCRC[s][n].status);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCR - multi sector bank
struct SCR_t	**SCR;
int 		 *SCR_NH;
int 		 *SCR_S;
int 		  SCR_NS;
void SCR_Print(){
  int s,n;
  if(SCR_NS){
    for(s=0;s<SCR_NS;s++){
      if(SCR_NH[s]){
        fprintf(stdout,"Group:\tSCR\tSector:\t%d \tNhits:\t%d\n",SCR_S[s],SCR_NH[s]);
        for(n=0;n<SCR_NH[s];n++){
          fprintf(stdout,"id: %d\tenergy: %8.3f\ttime: %8.3f\tx: %8.3f\ty: %8.3f\tz: %8.3f\tdx: %8.3f\tdy: %8.3f\tdz: %8.3f\tstatus: %d\tdenergy: %8.3f\tdtime: %8.3f\t",
		SCR[s][n].id,
		SCR[s][n].energy,
		SCR[s][n].time,
		SCR[s][n].x,
		SCR[s][n].y,
		SCR[s][n].z,
		SCR[s][n].dx,
		SCR[s][n].dy,
		SCR[s][n].dz,
		SCR[s][n].status,
		SCR[s][n].denergy,
		SCR[s][n].dtime);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCS - multi sector bank
struct SCS_t	**SCS;
int 		 *SCS_NH;
int 		 *SCS_S;
int 		  SCS_NS;
void SCS_Print(){
  int s,n;
  if(SCS_NS){
    for(s=0;s<SCS_NS;s++){
      if(SCS_NH[s]){
        fprintf(stdout,"Group:\tSCS\tSector:\t%d \tNhits:\t%d\n",SCS_S[s],SCS_NH[s]);
        for(n=0;n<SCS_NH[s];n++){
          fprintf(stdout,"SCS1: %d\tSCS2: %d\tSCS3: %d\tSCS4: %d\tSCS5: %d\tSCS6: %d\tSCS7: %d\tSCS8: %d\tSCS9: %d\tSCS10: %d\tSCS11: %d\tSCS12: %d\tSCS13: %d\tSCS14: %d\tSCS15: %d\tSCS16: %d\tSCS17: %d\tSCS18: %d\tSCS19: %d\tSCS20: %d\tSCS21: %d\tSCS22: %d\tSCS23: %d\tSCS24: %d\tSCS25: %d\tSCS26: %d\tSCS27: %d\tSCS28: %d\tSCS29: %d\tSCS30: %d\tSCS31: %d\tSCS32: %d\tSCS33: %d\tSCS34: %d\tSCS35: %d\tSCS36: %d\tSCS37: %d\tSCS38: %d\tSCS39: %d\tSCS40: %d\tSCS41: %d\tSCS42: %d\tSCS43: %d\tSCS44: %d\tSCS45: %d\tSCS46: %d\tSCS47: %d\tSCS48: %d\tSCS49: %d\tSCS50: %d\tSCS51: %d\tSCS52: %d\tSCS53: %d\tSCS54: %d\tSCS55: %d\tSCS56: %d\tSCS57: %d\tSCS58: %d\tSCS59: %d\tSCS60: %d\tSCS61: %d\tSCS62: %d\tSCS63: %d\tSCS64: %d\tSCS65: %d\tSCS66: %d\tSCS67: %d\tSCS68: %d\tSCS69: %d\tSCS70: %d\tSCS71: %d\tSCS72: %d\tSCS73: %d\tSCS74: %d\tSCS75: %d\tSCS76: %d\tSCS77: %d\tSCS78: %d\tSCS79: %d\tSCS80: %d\tSCS81: %d\tSCS82: %d\tSCS83: %d\tSCS84: %d\tSCS85: %d\tSCS86: %d\tSCS87: %d\tSCS88: %d\tSCS89: %d\tSCS90: %d\tSCS91: %d\tSCS92: %d\tSCS93: %d\tSCS94: %d\tSCS95: %d\tSCS96: %d\tSCS97: %d\tSCS98: %d\tSCS99: %d\tSCS100: %d\tSCS101: %d\tSCS102: %d\tSCS103: %d\tSCS104: %d\tSCS105: %d\tSCS106: %d\tSCS107: %d\tSCS108: %d\tSCS109: %d\tSCS110: %d\tSCS111: %d\tSCS112: %d\tSCS113: %d\tSCS114: %d\tSCS115: %d\tSCS116: %d\tSCS117: %d\tSCS118: %d\tSCS119: %d\tSCS120: %d\tSCS121: %d\tSCS122: %d\tSCS123: %d\tSCS124: %d\tSCS125: %d\tSCS126: %d\tSCS127: %d\tSCS128: %d\tSCS129: %d\tSCS130: %d\tSCS131: %d\tSCS132: %d\tSCS133: %d\tSCS134: %d\tSCS135: %d\tSCS136: %d\tSCS137: %d\tSCS138: %d\tSCS139: %d\tSCS140: %d\tSCS141: %d\tSCS142: %d\tSCS143: %d\tSCS144: %d\tSCS145: %d\tSCS146: %d\tSCS147: %d\tSCS148: %d\tSCS149: %d\tSCS150: %d\tSCS151: %d\tSCS152: %d\tSCS153: %d\tSCS154: %d\tSCS155: %d\tSCS156: %d\tSCS157: %d\tSCS158: %d\tSCS159: %d\tSCS160: %d\tSCS161: %d\tSCS162: %d\tSCS163: %d\tSCS164: %d\tSCS165: %d\tSCS166: %d\tSCS167: %d\tSCS168: %d\tSCS169: %d\tSCS170: %d\tSCS171: %d\tSCS172: %d\tSCS173: %d\tSCS174: %d\tSCS175: %d\tSCS176: %d\tSCS177: %d\tSCS178: %d\tSCS179: %d\tSCS180: %d\tSCS181: %d\tSCS182: %d\tSCS183: %d\tSCS184: %d\tSCS185: %d\tSCS186: %d\tSCS187: %d\tSCS188: %d\tSCS189: %d\tSCS190: %d\tSCS191: %d\tSCS192: %d\t",
		SCS[s][n].SCS1,
		SCS[s][n].SCS2,
		SCS[s][n].SCS3,
		SCS[s][n].SCS4,
		SCS[s][n].SCS5,
		SCS[s][n].SCS6,
		SCS[s][n].SCS7,
		SCS[s][n].SCS8,
		SCS[s][n].SCS9,
		SCS[s][n].SCS10,
		SCS[s][n].SCS11,
		SCS[s][n].SCS12,
		SCS[s][n].SCS13,
		SCS[s][n].SCS14,
		SCS[s][n].SCS15,
		SCS[s][n].SCS16,
		SCS[s][n].SCS17,
		SCS[s][n].SCS18,
		SCS[s][n].SCS19,
		SCS[s][n].SCS20,
		SCS[s][n].SCS21,
		SCS[s][n].SCS22,
		SCS[s][n].SCS23,
		SCS[s][n].SCS24,
		SCS[s][n].SCS25,
		SCS[s][n].SCS26,
		SCS[s][n].SCS27,
		SCS[s][n].SCS28,
		SCS[s][n].SCS29,
		SCS[s][n].SCS30,
		SCS[s][n].SCS31,
		SCS[s][n].SCS32,
		SCS[s][n].SCS33,
		SCS[s][n].SCS34,
		SCS[s][n].SCS35,
		SCS[s][n].SCS36,
		SCS[s][n].SCS37,
		SCS[s][n].SCS38,
		SCS[s][n].SCS39,
		SCS[s][n].SCS40,
		SCS[s][n].SCS41,
		SCS[s][n].SCS42,
		SCS[s][n].SCS43,
		SCS[s][n].SCS44,
		SCS[s][n].SCS45,
		SCS[s][n].SCS46,
		SCS[s][n].SCS47,
		SCS[s][n].SCS48,
		SCS[s][n].SCS49,
		SCS[s][n].SCS50,
		SCS[s][n].SCS51,
		SCS[s][n].SCS52,
		SCS[s][n].SCS53,
		SCS[s][n].SCS54,
		SCS[s][n].SCS55,
		SCS[s][n].SCS56,
		SCS[s][n].SCS57,
		SCS[s][n].SCS58,
		SCS[s][n].SCS59,
		SCS[s][n].SCS60,
		SCS[s][n].SCS61,
		SCS[s][n].SCS62,
		SCS[s][n].SCS63,
		SCS[s][n].SCS64,
		SCS[s][n].SCS65,
		SCS[s][n].SCS66,
		SCS[s][n].SCS67,
		SCS[s][n].SCS68,
		SCS[s][n].SCS69,
		SCS[s][n].SCS70,
		SCS[s][n].SCS71,
		SCS[s][n].SCS72,
		SCS[s][n].SCS73,
		SCS[s][n].SCS74,
		SCS[s][n].SCS75,
		SCS[s][n].SCS76,
		SCS[s][n].SCS77,
		SCS[s][n].SCS78,
		SCS[s][n].SCS79,
		SCS[s][n].SCS80,
		SCS[s][n].SCS81,
		SCS[s][n].SCS82,
		SCS[s][n].SCS83,
		SCS[s][n].SCS84,
		SCS[s][n].SCS85,
		SCS[s][n].SCS86,
		SCS[s][n].SCS87,
		SCS[s][n].SCS88,
		SCS[s][n].SCS89,
		SCS[s][n].SCS90,
		SCS[s][n].SCS91,
		SCS[s][n].SCS92,
		SCS[s][n].SCS93,
		SCS[s][n].SCS94,
		SCS[s][n].SCS95,
		SCS[s][n].SCS96,
		SCS[s][n].SCS97,
		SCS[s][n].SCS98,
		SCS[s][n].SCS99,
		SCS[s][n].SCS100,
		SCS[s][n].SCS101,
		SCS[s][n].SCS102,
		SCS[s][n].SCS103,
		SCS[s][n].SCS104,
		SCS[s][n].SCS105,
		SCS[s][n].SCS106,
		SCS[s][n].SCS107,
		SCS[s][n].SCS108,
		SCS[s][n].SCS109,
		SCS[s][n].SCS110,
		SCS[s][n].SCS111,
		SCS[s][n].SCS112,
		SCS[s][n].SCS113,
		SCS[s][n].SCS114,
		SCS[s][n].SCS115,
		SCS[s][n].SCS116,
		SCS[s][n].SCS117,
		SCS[s][n].SCS118,
		SCS[s][n].SCS119,
		SCS[s][n].SCS120,
		SCS[s][n].SCS121,
		SCS[s][n].SCS122,
		SCS[s][n].SCS123,
		SCS[s][n].SCS124,
		SCS[s][n].SCS125,
		SCS[s][n].SCS126,
		SCS[s][n].SCS127,
		SCS[s][n].SCS128,
		SCS[s][n].SCS129,
		SCS[s][n].SCS130,
		SCS[s][n].SCS131,
		SCS[s][n].SCS132,
		SCS[s][n].SCS133,
		SCS[s][n].SCS134,
		SCS[s][n].SCS135,
		SCS[s][n].SCS136,
		SCS[s][n].SCS137,
		SCS[s][n].SCS138,
		SCS[s][n].SCS139,
		SCS[s][n].SCS140,
		SCS[s][n].SCS141,
		SCS[s][n].SCS142,
		SCS[s][n].SCS143,
		SCS[s][n].SCS144,
		SCS[s][n].SCS145,
		SCS[s][n].SCS146,
		SCS[s][n].SCS147,
		SCS[s][n].SCS148,
		SCS[s][n].SCS149,
		SCS[s][n].SCS150,
		SCS[s][n].SCS151,
		SCS[s][n].SCS152,
		SCS[s][n].SCS153,
		SCS[s][n].SCS154,
		SCS[s][n].SCS155,
		SCS[s][n].SCS156,
		SCS[s][n].SCS157,
		SCS[s][n].SCS158,
		SCS[s][n].SCS159,
		SCS[s][n].SCS160,
		SCS[s][n].SCS161,
		SCS[s][n].SCS162,
		SCS[s][n].SCS163,
		SCS[s][n].SCS164,
		SCS[s][n].SCS165,
		SCS[s][n].SCS166,
		SCS[s][n].SCS167,
		SCS[s][n].SCS168,
		SCS[s][n].SCS169,
		SCS[s][n].SCS170,
		SCS[s][n].SCS171,
		SCS[s][n].SCS172,
		SCS[s][n].SCS173,
		SCS[s][n].SCS174,
		SCS[s][n].SCS175,
		SCS[s][n].SCS176,
		SCS[s][n].SCS177,
		SCS[s][n].SCS178,
		SCS[s][n].SCS179,
		SCS[s][n].SCS180,
		SCS[s][n].SCS181,
		SCS[s][n].SCS182,
		SCS[s][n].SCS183,
		SCS[s][n].SCS184,
		SCS[s][n].SCS185,
		SCS[s][n].SCS186,
		SCS[s][n].SCS187,
		SCS[s][n].SCS188,
		SCS[s][n].SCS189,
		SCS[s][n].SCS190,
		SCS[s][n].SCS191,
		SCS[s][n].SCS192);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SCT - multi sector bank
struct SCT_t	**SCT;
int 		 *SCT_NH;
int 		 *SCT_S;
int 		  SCT_NS;
void SCT_Print(){
  int s,n;
  if(SCT_NS){
    for(s=0;s<SCT_NS;s++){
      if(SCT_NH[s]){
        fprintf(stdout,"Group:\tSCT\tSector:\t%d \tNhits:\t%d\n",SCT_S[s],SCT_NH[s]);
        for(n=0;n<SCT_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		SCT[s][n].ID,
		SCT[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SGMP - multi sector bank
struct SGMP_t	**SGMP;
int 		 *SGMP_NH;
int 		 *SGMP_S;
int 		  SGMP_NS;
void SGMP_Print(){
  int s,n;
  if(SGMP_NS){
    for(s=0;s<SGMP_NS;s++){
      if(SGMP_NH[s]){
        fprintf(stdout,"Group:\tSGMP\tSector:\t%d \tNhits:\t%d\n",SGMP_S[s],SGMP_NH[s]);
        for(n=0;n<SGMP_NH[s];n++){
          fprintf(stdout,"is: %d\tsgm1: %d\tsgm2: %d\tsgm3: %d\tsgm4: %d\tsgm5: %d\tsgm6: %d\t",
		SGMP[s][n].is,
		SGMP[s][n].sgm1,
		SGMP[s][n].sgm2,
		SGMP[s][n].sgm3,
		SGMP[s][n].sgm4,
		SGMP[s][n].sgm5,
		SGMP[s][n].sgm6);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SPAR - multi sector bank
struct SPAR_t	**SPAR;
int 		 *SPAR_NH;
int 		 *SPAR_S;
int 		  SPAR_NS;
void SPAR_Print(){
  int s,n;
  if(SPAR_NS){
    for(s=0;s<SPAR_NS;s++){
      if(SPAR_NH[s]){
        fprintf(stdout,"Group:\tSPAR\tSector:\t%d \tNhits:\t%d\n",SPAR_S[s],SPAR_NH[s]);
        for(n=0;n<SPAR_NH[s];n++){
          fprintf(stdout,"slot: %d\tchannel: %d\tspar: %d\tpedmean: %d\t",
		SPAR[s][n].slot,
		SPAR[s][n].channel,
		SPAR[s][n].spar,
		SPAR[s][n].pedmean);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SPIN - multi sector bank
struct SPIN_t	**SPIN;
int 		 *SPIN_NH;
int 		 *SPIN_S;
int 		  SPIN_NS;
void SPIN_Print(){
  int s,n;
  if(SPIN_NS){
    for(s=0;s<SPIN_NS;s++){
      if(SPIN_NH[s]){
        fprintf(stdout,"Group:\tSPIN\tSector:\t%d \tNhits:\t%d\n",SPIN_S[s],SPIN_NH[s]);
        for(n=0;n<SPIN_NH[s];n++){
          fprintf(stdout,"xspin: %d\tyspin: %d\tzspin: %d\t",
		SPIN[s][n].xspin,
		SPIN[s][n].yspin,
		SPIN[s][n].zspin);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ST1 - multi sector bank
struct ST1_t	**ST1;
int 		 *ST1_NH;
int 		 *ST1_S;
int 		  ST1_NS;
void ST1_Print(){
  int s,n;
  if(ST1_NS){
    for(s=0;s<ST1_NS;s++){
      if(ST1_NH[s]){
        fprintf(stdout,"Group:\tST1\tSector:\t%d \tNhits:\t%d\n",ST1_S[s],ST1_NH[s]);
        for(n=0;n<ST1_NH[s];n++){
          fprintf(stdout,"ID: %d\tstatus: %d\tTIME_1: %8.3f\tADC_1: %8.3f\tTIME_2: %8.3f\tADC_2: %8.3f\t",
		ST1[s][n].ID,
		ST1[s][n].status,
		ST1[s][n].TIME_1,
		ST1[s][n].ADC_1,
		ST1[s][n].TIME_2,
		ST1[s][n].ADC_2);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank ST - single sector bank
struct ST_t	 *ST;
int 		  ST_NH;
int 		  ST_S;
void ST_Print(){
  int n;
  if(ST_NH){
    fprintf(stdout,"Group:\tST\tSector:\t%d \tNhits:\t%d\n",ST_S,ST_NH);
    for(n=0;n<ST_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\tADC: %hd\t",
		ST[n].ID,
		ST[n].TDC,
		ST[n].ADC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STG - multi sector bank
struct STG_t	**STG;
int 		 *STG_NH;
int 		 *STG_S;
int 		  STG_NS;
void STG_Print(){
  int s,n;
  if(STG_NS){
    for(s=0;s<STG_NS;s++){
      if(STG_NH[s]){
        fprintf(stdout,"Group:\tSTG\tSector:\t%d \tNhits:\t%d\n",STG_S[s],STG_NH[s]);
        for(n=0;n<STG_NH[s];n++){
          fprintf(stdout,"id: %d\tleg_x_max: %8.3f\tleg_x_min: %8.3f\tleg_y_max: %8.3f\tleg_y_min: %8.3f\tleg_z_max: %8.3f\tlet_z_min: %8.3f\tnose_x_max: %8.3f\tnose_x_min: %8.3f\tnose_y_max: %8.3f\tnose_y_min: %8.3f\tnose_z_max: %8.3f\tnose_z_min: %8.3f\tnoseAngle: %8.3f\tlleg: %8.3f\tlnose: %8.3f\t",
		STG[s][n].id,
		STG[s][n].leg_x_max,
		STG[s][n].leg_x_min,
		STG[s][n].leg_y_max,
		STG[s][n].leg_y_min,
		STG[s][n].leg_z_max,
		STG[s][n].let_z_min,
		STG[s][n].nose_x_max,
		STG[s][n].nose_x_min,
		STG[s][n].nose_y_max,
		STG[s][n].nose_y_min,
		STG[s][n].nose_z_max,
		STG[s][n].nose_z_min,
		STG[s][n].noseAngle,
		STG[s][n].lleg,
		STG[s][n].lnose);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank STH - multi sector bank
struct STH_t	**STH;
int 		 *STH_NH;
int 		 *STH_S;
int 		  STH_NS;
void STH_Print(){
  int s,n;
  if(STH_NS){
    for(s=0;s<STH_NS;s++){
      if(STH_NH[s]){
        fprintf(stdout,"Group:\tSTH\tSector:\t%d \tNhits:\t%d\n",STH_S[s],STH_NH[s]);
        for(n=0;n<STH_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\tpmom: %8.3f\ttrack: %d\tid: %d\ttof: %8.3f\t",
		STH[s][n].x,
		STH[s][n].y,
		STH[s][n].z,
		STH[s][n].cx,
		STH[s][n].cy,
		STH[s][n].cz,
		STH[s][n].pmom,
		STH[s][n].track,
		STH[s][n].id,
		STH[s][n].tof);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank STN0 - single sector bank
struct STN0_t	 *STN0;
int 		  STN0_NH;
int 		  STN0_S;
void STN0_Print(){
  int n;
  if(STN0_NH){
    fprintf(stdout,"Group:\tSTN0\tSector:\t%d \tNhits:\t%d\n",STN0_S,STN0_NH);
    for(n=0;n<STN0_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		STN0[n].ID,
		STN0[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STN1 - single sector bank
struct STN1_t	 *STN1;
int 		  STN1_NH;
int 		  STN1_S;
void STN1_Print(){
  int n;
  if(STN1_NH){
    fprintf(stdout,"Group:\tSTN1\tSector:\t%d \tNhits:\t%d\n",STN1_S,STN1_NH);
    for(n=0;n<STN1_NH;n++){
      fprintf(stdout,"ID: %hd\tADC: %hd\t",
		STN1[n].ID,
		STN1[n].ADC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STPB - single sector bank
struct STPB_t	 *STPB;
int 		  STPB_NH;
int 		  STPB_S;
void STPB_Print(){
  int n;
  if(STPB_NH){
    fprintf(stdout,"Group:\tSTPB\tSector:\t%d \tNhits:\t%d\n",STPB_S,STPB_NH);
    for(n=0;n<STPB_NH;n++){
      fprintf(stdout,"SThid: %d\tTime: %8.3f\tPath: %8.3f\tcharge: %d\tStatus: %d\t",
		STPB[n].SThid,
		STPB[n].Time,
		STPB[n].Path,
		STPB[n].charge,
		STPB[n].Status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STPE - multi sector bank
struct STPE_t	**STPE;
int 		 *STPE_NH;
int 		 *STPE_S;
int 		  STPE_NS;
void STPE_Print(){
  int s,n;
  if(STPE_NS){
    for(s=0;s<STPE_NS;s++){
      if(STPE_NH[s]){
        fprintf(stdout,"Group:\tSTPE\tSector:\t%d \tNhits:\t%d\n",STPE_S[s],STPE_NH[s]);
        for(n=0;n<STPE_NH[s];n++){
          fprintf(stdout,"ID: %d\tmean: %d\tsigma: %8.3f\t",
		STPE[s][n].ID,
		STPE[s][n].mean,
		STPE[s][n].sigma);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank STR - single sector bank
struct STR_t	 *STR;
int 		  STR_NH;
int 		  STR_S;
void STR_Print(){
  int n;
  if(STR_NH){
    fprintf(stdout,"Group:\tSTR\tSector:\t%d \tNhits:\t%d\n",STR_S,STR_NH);
    for(n=0;n<STR_NH;n++){
      fprintf(stdout,"ID: %d\tTrk_no: %d\tST_TIME: %8.3f\tST_L: %8.3f\tst_pos: %8.3f\tstatus: %d\t",
		STR[n].ID,
		STR[n].Trk_no,
		STR[n].ST_TIME,
		STR[n].ST_L,
		STR[n].st_pos,
		STR[n].status);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STS - multi sector bank
struct STS_t	**STS;
int 		 *STS_NH;
int 		 *STS_S;
int 		  STS_NS;
void STS_Print(){
  int s,n;
  if(STS_NS){
    for(s=0;s<STS_NS;s++){
      if(STS_NH[s]){
        fprintf(stdout,"Group:\tSTS\tSector:\t%d \tNhits:\t%d\n",STS_S[s],STS_NH[s]);
        for(n=0;n<STS_NH[s];n++){
          fprintf(stdout,"STS1: %d\tSTS2: %d\tSTS3: %d\tSTS4: %d\tSTS5: %d\tSTS6: %d\tSTS7: %d\tSTS8: %d\tSTS9: %d\tSTS10: %d\tSTS11: %d\tSTS12: %d\tSTS13: %d\tSTS14: %d\tSTS15: %d\tSTS16: %d\t",
		STS[s][n].STS1,
		STS[s][n].STS2,
		STS[s][n].STS3,
		STS[s][n].STS4,
		STS[s][n].STS5,
		STS[s][n].STS6,
		STS[s][n].STS7,
		STS[s][n].STS8,
		STS[s][n].STS9,
		STS[s][n].STS10,
		STS[s][n].STS11,
		STS[s][n].STS12,
		STS[s][n].STS13,
		STS[s][n].STS14,
		STS[s][n].STS15,
		STS[s][n].STS16);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank STSN - multi sector bank
struct STSN_t	**STSN;
int 		 *STSN_NH;
int 		 *STSN_S;
int 		  STSN_NS;
void STSN_Print(){
  int s,n;
  if(STSN_NS){
    for(s=0;s<STSN_NS;s++){
      if(STSN_NH[s]){
        fprintf(stdout,"Group:\tSTSN\tSector:\t%d \tNhits:\t%d\n",STSN_S[s],STSN_NH[s]);
        for(n=0;n<STSN_NH[s];n++){
          fprintf(stdout,"STSN1: %d\tSTSN2: %d\tSTSN3: %d\tSTSN4: %d\tSTSN5: %d\tSTSN6: %d\tSTSN7: %d\tSTSN8: %d\tSTSN9: %d\tSTSN10: %d\tSTSN11: %d\tSTSN12: %d\tSTSN13: %d\tSTSN14: %d\tSTSN15: %d\tSTSN16: %d\tSTSN17: %d\tSTSN18: %d\tSTSN19: %d\tSTSN20: %d\tSTSN21: %d\tSTSN22: %d\tSTSN23: %d\tSTSN24: %d\tSTOR: %d\tSTMULT: %d\tSTANDMOR: %d\tMULTANDMOR: %d\tRES1: %d\tRES2: %d\tRES3: %d\tRES4: %d\t",
		STSN[s][n].STSN1,
		STSN[s][n].STSN2,
		STSN[s][n].STSN3,
		STSN[s][n].STSN4,
		STSN[s][n].STSN5,
		STSN[s][n].STSN6,
		STSN[s][n].STSN7,
		STSN[s][n].STSN8,
		STSN[s][n].STSN9,
		STSN[s][n].STSN10,
		STSN[s][n].STSN11,
		STSN[s][n].STSN12,
		STSN[s][n].STSN13,
		STSN[s][n].STSN14,
		STSN[s][n].STSN15,
		STSN[s][n].STSN16,
		STSN[s][n].STSN17,
		STSN[s][n].STSN18,
		STSN[s][n].STSN19,
		STSN[s][n].STSN20,
		STSN[s][n].STSN21,
		STSN[s][n].STSN22,
		STSN[s][n].STSN23,
		STSN[s][n].STSN24,
		STSN[s][n].STOR,
		STSN[s][n].STMULT,
		STSN[s][n].STANDMOR,
		STSN[s][n].MULTANDMOR,
		STSN[s][n].RES1,
		STSN[s][n].RES2,
		STSN[s][n].RES3,
		STSN[s][n].RES4);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank SYNC - multi sector bank
struct SYNC_t	**SYNC;
int 		 *SYNC_NH;
int 		 *SYNC_S;
int 		  SYNC_NS;
void SYNC_Print(){
  int s,n;
  if(SYNC_NS){
    for(s=0;s<SYNC_NS;s++){
      if(SYNC_NH[s]){
        fprintf(stdout,"Group:\tSYNC\tSector:\t%d \tNhits:\t%d\n",SYNC_S[s],SYNC_NH[s]);
        for(n=0;n<SYNC_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDCL: %hd\tADCL: %hd\t",
		SYNC[s][n].ID,
		SYNC[s][n].TDCL,
		SYNC[s][n].ADCL);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TACO - multi sector bank
struct TACO_t	**TACO;
int 		 *TACO_NH;
int 		 *TACO_S;
int 		  TACO_NS;
void TACO_Print(){
  int s,n;
  if(TACO_NS){
    for(s=0;s<TACO_NS;s++){
      if(TACO_NH[s]){
        fprintf(stdout,"Group:\tTACO\tSector:\t%d \tNhits:\t%d\n",TACO_S[s],TACO_NH[s]);
        for(n=0;n<TACO_NH[s];n++){
          fprintf(stdout,"ID: %d\tTIME: %8.3f\tELT: %8.3f\tERT: %8.3f\tELB: %8.3f\tERB: %8.3f\tESUM: %8.3f\tESUM2: %8.3f\tESUM3: %8.3f\tTID: %d\t",
		TACO[s][n].ID,
		TACO[s][n].TIME,
		TACO[s][n].ELT,
		TACO[s][n].ERT,
		TACO[s][n].ELB,
		TACO[s][n].ERB,
		TACO[s][n].ESUM,
		TACO[s][n].ESUM2,
		TACO[s][n].ESUM3,
		TACO[s][n].TID);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TAGE - single sector bank
struct TAGE_t	 *TAGE;
int 		  TAGE_NH;
int 		  TAGE_S;
void TAGE_Print(){
  int n;
  if(TAGE_NH){
    fprintf(stdout,"Group:\tTAGE\tSector:\t%d \tNhits:\t%d\n",TAGE_S,TAGE_NH);
    for(n=0;n<TAGE_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		TAGE[n].ID,
		TAGE[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TAGI - multi sector bank
struct TAGI_t	**TAGI;
int 		 *TAGI_NH;
int 		 *TAGI_S;
int 		  TAGI_NS;
void TAGI_Print(){
  int s,n;
  if(TAGI_NS){
    for(s=0;s<TAGI_NS;s++){
      if(TAGI_NH[s]){
        fprintf(stdout,"Group:\tTAGI\tSector:\t%d \tNhits:\t%d\n",TAGI_S[s],TAGI_NH[s]);
        for(n=0;n<TAGI_NH[s];n++){
          fprintf(stdout,"IDT: %d\tTIMEL: %8.3f\tTIMER: %8.3f\tIDE: %d\tTIMEE: %8.3f\tTIMEMEAN: %8.3f\tTRF: %8.3f\tNEXTTIME: %8.3f\t",
		TAGI[s][n].IDT,
		TAGI[s][n].TIMEL,
		TAGI[s][n].TIMER,
		TAGI[s][n].IDE,
		TAGI[s][n].TIMEE,
		TAGI[s][n].TIMEMEAN,
		TAGI[s][n].TRF,
		TAGI[s][n].NEXTTIME);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TAGM - single sector bank
struct TAGM_t	 *TAGM;
int 		  TAGM_NH;
int 		  TAGM_S;
void TAGM_Print(){
  int n;
  if(TAGM_NH){
    fprintf(stdout,"Group:\tTAGM\tSector:\t%d \tNhits:\t%d\n",TAGM_S,TAGM_NH);
    for(n=0;n<TAGM_NH;n++){
      fprintf(stdout,"energy: %8.3f\tt: %8.3f\te_t: %8.3f\tstatus: %d\ttid: %d\teid: %d\t",
		TAGM[n].energy,
		TAGM[n].t,
		TAGM[n].e_t,
		TAGM[n].status,
		TAGM[n].tid,
		TAGM[n].eid);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TAGR - single sector bank
struct TAGR_t	 *TAGR;
int 		  TAGR_NH;
int 		  TAGR_S;
void TAGR_Print(){
  int n;
  if(TAGR_NH){
    fprintf(stdout,"Group:\tTAGR\tSector:\t%d \tNhits:\t%d\n",TAGR_S,TAGR_NH);
    for(n=0;n<TAGR_NH;n++){
      fprintf(stdout,"ERG: %8.3f\tTTAG: %8.3f\tTPHO: %8.3f\tSTAT: %d\tT_id: %d\tE_id: %d\t",
		TAGR[n].ERG,
		TAGR[n].TTAG,
		TAGR[n].TPHO,
		TAGR[n].STAT,
		TAGR[n].T_id,
		TAGR[n].E_id);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TAGT - single sector bank
struct TAGT_t	 *TAGT;
int 		  TAGT_NH;
int 		  TAGT_S;
void TAGT_Print(){
  int n;
  if(TAGT_NH){
    fprintf(stdout,"Group:\tTAGT\tSector:\t%d \tNhits:\t%d\n",TAGT_S,TAGT_NH);
    for(n=0;n<TAGT_NH;n++){
      fprintf(stdout,"ID: %hd\tTDCL: %hd\tTDCR: %hd\t",
		TAGT[n].ID,
		TAGT[n].TDCL,
		TAGT[n].TDCR);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TATL - multi sector bank
struct TATL_t	**TATL;
int 		 *TATL_NH;
int 		 *TATL_S;
int 		  TATL_NS;
void TATL_Print(){
  int s,n;
  if(TATL_NS){
    for(s=0;s<TATL_NS;s++){
      if(TATL_NH[s]){
        fprintf(stdout,"Group:\tTATL\tSector:\t%d \tNhits:\t%d\n",TATL_S[s],TATL_NH[s]);
        for(n=0;n<TATL_NH[s];n++){
          fprintf(stdout,"ID: %d\tTDC: %d\t",
		TATL[s][n].ID,
		TATL[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TATR - multi sector bank
struct TATR_t	**TATR;
int 		 *TATR_NH;
int 		 *TATR_S;
int 		  TATR_NS;
void TATR_Print(){
  int s,n;
  if(TATR_NS){
    for(s=0;s<TATR_NS;s++){
      if(TATR_NH[s]){
        fprintf(stdout,"Group:\tTATR\tSector:\t%d \tNhits:\t%d\n",TATR_S[s],TATR_NH[s]);
        for(n=0;n<TATR_NH[s];n++){
          fprintf(stdout,"ID: %d\tTDC: %d\t",
		TATR[s][n].ID,
		TATR[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TBER - single sector bank
struct TBER_t	 *TBER;
int 		  TBER_NH;
int 		  TBER_S;
void TBER_Print(){
  int n;
  if(TBER_NH){
    fprintf(stdout,"Group:\tTBER\tSector:\t%d \tNhits:\t%d\n",TBER_S,TBER_NH);
    for(n=0;n<TBER_NH;n++){
      fprintf(stdout,"q_over_p: %8.3f\tlambda: %8.3f\tphi: %8.3f\td0: %8.3f\tz0: %8.3f\tc11: %8.3f\tc12: %8.3f\tc13: %8.3f\tc14: %8.3f\tc15: %8.3f\tc22: %8.3f\tc23: %8.3f\tc24: %8.3f\tc25: %8.3f\tc33: %8.3f\tc34: %8.3f\tc35: %8.3f\tc44: %8.3f\tc45: %8.3f\tc55: %8.3f\tchi2: %8.3f\tlayinfo1: %d\tlayinfo2: %d\t",
		TBER[n].q_over_p,
		TBER[n].lambda,
		TBER[n].phi,
		TBER[n].d0,
		TBER[n].z0,
		TBER[n].c11,
		TBER[n].c12,
		TBER[n].c13,
		TBER[n].c14,
		TBER[n].c15,
		TBER[n].c22,
		TBER[n].c23,
		TBER[n].c24,
		TBER[n].c25,
		TBER[n].c33,
		TBER[n].c34,
		TBER[n].c35,
		TBER[n].c44,
		TBER[n].c45,
		TBER[n].c55,
		TBER[n].chi2,
		TBER[n].layinfo1,
		TBER[n].layinfo2);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TBID - single sector bank
struct TBID_t	 *TBID;
int 		  TBID_NH;
int 		  TBID_S;
void TBID_Print(){
  int n;
  if(TBID_NH){
    fprintf(stdout,"Group:\tTBID\tSector:\t%d \tNhits:\t%d\n",TBID_S,TBID_NH);
    for(n=0;n<TBID_NH;n++){
      fprintf(stdout,"track: %d\tsec: %d\tbeta: %8.3f\tvtime: %8.3f\tsc_stat: %d\tsc_id: %d\tsc_time: %8.3f\tsc_qual: %8.3f\tsc_vtime: %8.3f\tsc_beta: %8.3f\tcc_stat: %d\tcc_id: %d\tcc_time: %8.3f\tcc_qual: %8.3f\tcc_vtime: %8.3f\tcc_beta: %8.3f\tec_stat: %d\tec_id: %d\tec_time: %8.3f\tec_qual: %8.3f\tec_vtime: %8.3f\tec_beta: %8.3f\tst_stat: %d\tst_id: %d\tst_time: %8.3f\tst_qual: %8.3f\tst_vtime: %8.3f\tst_beta: %8.3f\tlac_stat: %d\tlac_id: %d\tlac_time: %8.3f\tlac_qual: %8.3f\tlac_vtime: %8.3f\tlac_beta: %8.3f\t",
		TBID[n].track,
		TBID[n].sec,
		TBID[n].beta,
		TBID[n].vtime,
		TBID[n].sc_stat,
		TBID[n].sc_id,
		TBID[n].sc_time,
		TBID[n].sc_qual,
		TBID[n].sc_vtime,
		TBID[n].sc_beta,
		TBID[n].cc_stat,
		TBID[n].cc_id,
		TBID[n].cc_time,
		TBID[n].cc_qual,
		TBID[n].cc_vtime,
		TBID[n].cc_beta,
		TBID[n].ec_stat,
		TBID[n].ec_id,
		TBID[n].ec_time,
		TBID[n].ec_qual,
		TBID[n].ec_vtime,
		TBID[n].ec_beta,
		TBID[n].st_stat,
		TBID[n].st_id,
		TBID[n].st_time,
		TBID[n].st_qual,
		TBID[n].st_vtime,
		TBID[n].st_beta,
		TBID[n].lac_stat,
		TBID[n].lac_id,
		TBID[n].lac_time,
		TBID[n].lac_qual,
		TBID[n].lac_vtime,
		TBID[n].lac_beta);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TBLA - multi sector bank
struct TBLA_t	**TBLA;
int 		 *TBLA_NH;
int 		 *TBLA_S;
int 		  TBLA_NS;
void TBLA_Print(){
  int s,n;
  if(TBLA_NS){
    for(s=0;s<TBLA_NS;s++){
      if(TBLA_NH[s]){
        fprintf(stdout,"Group:\tTBLA\tSector:\t%d \tNhits:\t%d\n",TBLA_S[s],TBLA_NH[s]);
        for(n=0;n<TBLA_NH[s];n++){
          fprintf(stdout,"trk_pln: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tBx: %8.3f\tBy: %8.3f\tBz: %8.3f\ttlen: %8.3f\tdc1: %d\tstat: %d\twire: %d\tdtime: %8.3f\talpha: %8.3f\twlen: %8.3f\tsgdoca: %8.3f\tfitdoca: %8.3f\tcalcdoca: %8.3f\t",
		TBLA[s][n].trk_pln,
		TBLA[s][n].x,
		TBLA[s][n].y,
		TBLA[s][n].z,
		TBLA[s][n].Bx,
		TBLA[s][n].By,
		TBLA[s][n].Bz,
		TBLA[s][n].tlen,
		TBLA[s][n].dc1,
		TBLA[s][n].stat,
		TBLA[s][n].wire,
		TBLA[s][n].dtime,
		TBLA[s][n].alpha,
		TBLA[s][n].wlen,
		TBLA[s][n].sgdoca,
		TBLA[s][n].fitdoca,
		TBLA[s][n].calcdoca);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TBTR - single sector bank
struct TBTR_t	 *TBTR;
int 		  TBTR_NH;
int 		  TBTR_S;
void TBTR_Print(){
  int n;
  if(TBTR_NH){
    fprintf(stdout,"Group:\tTBTR\tSector:\t%d \tNhits:\t%d\n",TBTR_S,TBTR_NH);
    for(n=0;n<TBTR_NH;n++){
      fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tpx: %8.3f\tpy: %8.3f\tpz: %8.3f\tq: %8.3f\tchi2: %8.3f\titr_sec: %d\titr_hbt: %d\t",
		TBTR[n].x,
		TBTR[n].y,
		TBTR[n].z,
		TBTR[n].px,
		TBTR[n].py,
		TBTR[n].pz,
		TBTR[n].q,
		TBTR[n].chi2,
		TBTR[n].itr_sec,
		TBTR[n].itr_hbt);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TCSB - multi sector bank
struct TCSB_t	**TCSB;
int 		 *TCSB_NH;
int 		 *TCSB_S;
int 		  TCSB_NS;
void TCSB_Print(){
  int s,n;
  if(TCSB_NS){
    for(s=0;s<TCSB_NS;s++){
      if(TCSB_NH[s]){
        fprintf(stdout,"Group:\tTCSB\tSector:\t%d \tNhits:\t%d\n",TCSB_S[s],TCSB_NH[s]);
        for(n=0;n<TCSB_NH[s];n++){
          fprintf(stdout,"xpos: %8.3f\typos: %8.3f\tzpos: %8.3f\tsxpos: %8.3f\tsypos: %8.3f\tszpos: %8.3f\t",
		TCSB[s][n].xpos,
		TCSB[s][n].ypos,
		TCSB[s][n].zpos,
		TCSB[s][n].sxpos,
		TCSB[s][n].sypos,
		TCSB[s][n].szpos);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TCT - multi sector bank
struct TCT_t	**TCT;
int 		 *TCT_NH;
int 		 *TCT_S;
int 		  TCT_NS;
void TCT_Print(){
  int s,n;
  if(TCT_NS){
    for(s=0;s<TCT_NS;s++){
      if(TCT_NH[s]){
        fprintf(stdout,"Group:\tTCT\tSector:\t%d \tNhits:\t%d\n",TCT_S[s],TCT_NH[s]);
        for(n=0;n<TCT_NH[s];n++){
          fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		TCT[s][n].ID,
		TCT[s][n].TDC);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TDPL - multi sector bank
struct TDPL_t	**TDPL;
int 		 *TDPL_NH;
int 		 *TDPL_S;
int 		  TDPL_NS;
void TDPL_Print(){
  int s,n;
  if(TDPL_NS){
    for(s=0;s<TDPL_NS;s++){
      if(TDPL_NH[s]){
        fprintf(stdout,"Group:\tTDPL\tSector:\t%d \tNhits:\t%d\n",TDPL_S[s],TDPL_NH[s]);
        for(n=0;n<TDPL_NH[s];n++){
          fprintf(stdout,"trk_pln: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\ttlen: %8.3f\t",
		TDPL[s][n].trk_pln,
		TDPL[s][n].x,
		TDPL[s][n].y,
		TDPL[s][n].z,
		TDPL[s][n].cx,
		TDPL[s][n].cy,
		TDPL[s][n].cz,
		TDPL[s][n].tlen);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TESC - single sector bank
struct TESC_t	 *TESC;
int 		  TESC_NH;
int 		  TESC_S;
void TESC_Print(){
  int n;
  if(TESC_NH){
    fprintf(stdout,"Group:\tTESC\tSector:\t%d \tNhits:\t%d\n",TESC_S,TESC_NH);
    for(n=0;n<TESC_NH;n++){
      fprintf(stdout,"value: %8.3f\t",
		TESC[n].value);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TGBI - single sector bank
struct TGBI_t	 *TGBI;
int 		  TGBI_NH;
int 		  TGBI_S;
void TGBI_Print(){
  int n;
  if(TGBI_NH){
    fprintf(stdout,"Group:\tTGBI\tSector:\t%d \tNhits:\t%d\n",TGBI_S,TGBI_NH);
    for(n=0;n<TGBI_NH;n++){
      fprintf(stdout,"latch1: %d\thelicity_scaler: %d\tinterrupt_time: %d\tlatch2: %d\tlevel3: %d\t",
		TGBI[n].latch1,
		TGBI[n].helicity_scaler,
		TGBI[n].interrupt_time,
		TGBI[n].latch2,
		TGBI[n].level3);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TGEO - multi sector bank
struct TGEO_t	**TGEO;
int 		 *TGEO_NH;
int 		 *TGEO_S;
int 		  TGEO_NS;
void TGEO_Print(){
  int s,n;
  if(TGEO_NS){
    for(s=0;s<TGEO_NS;s++){
      if(TGEO_NH[s]){
        fprintf(stdout,"Group:\tTGEO\tSector:\t%d \tNhits:\t%d\n",TGEO_S[s],TGEO_NH[s]);
        for(n=0;n<TGEO_NH[s];n++){
          fprintf(stdout,"x: %8.3f\ty: %8.3f\tz: %8.3f\tradius: %8.3f\tlength: %8.3f\tmaterial: %d\t",
		TGEO[s][n].x,
		TGEO[s][n].y,
		TGEO[s][n].z,
		TGEO[s][n].radius,
		TGEO[s][n].length,
		TGEO[s][n].material);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TGPB - single sector bank
struct TGPB_t	 *TGPB;
int 		  TGPB_NH;
int 		  TGPB_S;
void TGPB_Print(){
  int n;
  if(TGPB_NH){
    fprintf(stdout,"Group:\tTGPB\tSector:\t%d \tNhits:\t%d\n",TGPB_S,TGPB_NH);
    for(n=0;n<TGPB_NH;n++){
      fprintf(stdout,"pointer: %d\tTime: %8.3f\tEnergy: %8.3f\tdt: %8.3f\t",
		TGPB[n].pointer,
		TGPB[n].Time,
		TGPB[n].Energy,
		TGPB[n].dt);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TGS - multi sector bank
struct TGS_t	**TGS;
int 		 *TGS_NH;
int 		 *TGS_S;
int 		  TGS_NS;
void TGS_Print(){
  int s,n;
  if(TGS_NS){
    for(s=0;s<TGS_NS;s++){
      if(TGS_NH[s]){
        fprintf(stdout,"Group:\tTGS\tSector:\t%d \tNhits:\t%d\n",TGS_S[s],TGS_NH[s]);
        for(n=0;n<TGS_NH[s];n++){
          fprintf(stdout,"RAWT1: %d\tRAWT2: %d\tRAWT3: %d\tRAWT4: %d\tRAWT5: %d\tRAWT6: %d\tRAWT7: %d\tRAWT8: %d\tRAWT9: %d\tRAWT10: %d\tRAWT11: %d\tRAWT12: %d\tRAWT13: %d\tRAWT14: %d\tRAWT15: %d\tRAWT16: %d\tRAWT17: %d\tRAWT18: %d\tRAWT19: %d\tRAWT20: %d\tRAWT21: %d\tRAWT22: %d\tRAWT23: %d\tRAWT24: %d\tRAWT25: %d\tRAWT26: %d\tRAWT27: %d\tRAWT28: %d\tRAWT29: %d\tRAWT30: %d\tRAWT31: %d\tRAWT32: %d\tRAWT33: %d\tRAWT34: %d\tRAWT35: %d\tRAWT36: %d\tRAWT37: %d\tRAWT38: %d\tRAWT39: %d\tRAWT40: %d\tRAWT41: %d\tRAWT42: %d\tRAWT43: %d\tRAWT44: %d\tRAWT45: %d\tRAWT46: %d\tRAWT47: %d\tRAWT48: %d\tRAWT49: %d\tRAWT50: %d\tRAWT51: %d\tRAWT52: %d\tRAWT53: %d\tRAWT54: %d\tRAWT55: %d\tRAWT56: %d\tRAWT57: %d\tRAWT58: %d\tRAWT59: %d\tRAWT60: %d\tRAWT61: %d\tnotused62: %d\tnotused63: %d\tnotused64: %d\tBNK1T1: %d\tBNK1T2: %d\tBNK1T3: %d\tBNK1T4: %d\tBNK1T5: %d\tBNK1T6: %d\tBNK1T7: %d\tBNK1T8: %d\tBNK1T9: %d\tBNK1T10: %d\tBNK1T11: %d\tBNK1T12: %d\tBNK1T13: %d\tBNK1T14: %d\tBNK1T15: %d\tBNK1T16: %d\tBNK1T17: %d\tBNK1T18: %d\tBNK1T19: %d\tBNK1T20: %d\tBNK1T21: %d\tBNK1T22: %d\tBNK1T23: %d\tBNK1T24: %d\tBNK1T25: %d\tBNK1T26: %d\tBNK1T27: %d\tBNK1T28: %d\tBNK1T29: %d\tBNK1T30: %d\tBNK1T31: %d\tBNK1T32: %d\tBNK1T33: %d\tBNK1T34: %d\tBNK1T35: %d\tBNK1T36: %d\tBNK1T37: %d\tBNK1T38: %d\tBNK1T39: %d\tBNK1T40: %d\tBNK1T41: %d\tBNK1T42: %d\tBNK1T43: %d\tBNK1T44: %d\tBNK1T45: %d\tBNK1T46: %d\tBNK1T47: %d\tBNK1T48: %d\tBNK1T49: %d\tBNK1T50: %d\tBNK1T51: %d\tBNK1T52: %d\tBNK1T53: %d\tBNK1T54: %d\tBNK1T55: %d\tBNK1T56: %d\tBNK1T57: %d\tBNK1T58: %d\tBNK1T59: %d\tBNK1T60: %d\tBNK1T61: %d\tnotused126: %d\tnotused127: %d\tnotused128: %d\tBNK2T1: %d\tBNK2T2: %d\tBNK2T3: %d\tBNK2T4: %d\tBNK2T5: %d\tBNK2T6: %d\tBNK2T7: %d\tBNK2T8: %d\tBNK2T9: %d\tBNK2T10: %d\tBNK2T11: %d\tBNK2T12: %d\tBNK2T13: %d\tBNK2T14: %d\tBNK2T15: %d\tBNK2T16: %d\tBNK2T17: %d\tBNK2T18: %d\tBNK2T19: %d\tBNK2T20: %d\tBNK2T21: %d\tBNK2T22: %d\tBNK2T23: %d\tBNK2T24: %d\tBNK2T25: %d\tBNK2T26: %d\tBNK2T27: %d\tBNK2T28: %d\tBNK2T29: %d\tBNK2T30: %d\tBNK2T31: %d\tBNK2T32: %d\tBNK2T33: %d\tBNK2T34: %d\tBNK2T35: %d\tBNK2T36: %d\tBNK2T37: %d\tBNK2T38: %d\tBNK2T39: %d\tBNK2T40: %d\tBNK2T41: %d\tBNK2T42: %d\tBNK2T43: %d\tBNK2T44: %d\tBNK2T45: %d\tBNK2T46: %d\tBNK2T47: %d\tBNK2T48: %d\tBNK2T49: %d\tBNK2T50: %d\tBNK2T51: %d\tBNK2T52: %d\tBNK2T53: %d\tBNK2T54: %d\tBNK2T55: %d\tBNK2T56: %d\tBNK2T57: %d\tBNK2T58: %d\tBNK2T59: %d\tBNK2T60: %d\tBNK2T61: %d\tnotused190: %d\tnotused191: %d\tnotused192: %d\t",
		TGS[s][n].RAWT1,
		TGS[s][n].RAWT2,
		TGS[s][n].RAWT3,
		TGS[s][n].RAWT4,
		TGS[s][n].RAWT5,
		TGS[s][n].RAWT6,
		TGS[s][n].RAWT7,
		TGS[s][n].RAWT8,
		TGS[s][n].RAWT9,
		TGS[s][n].RAWT10,
		TGS[s][n].RAWT11,
		TGS[s][n].RAWT12,
		TGS[s][n].RAWT13,
		TGS[s][n].RAWT14,
		TGS[s][n].RAWT15,
		TGS[s][n].RAWT16,
		TGS[s][n].RAWT17,
		TGS[s][n].RAWT18,
		TGS[s][n].RAWT19,
		TGS[s][n].RAWT20,
		TGS[s][n].RAWT21,
		TGS[s][n].RAWT22,
		TGS[s][n].RAWT23,
		TGS[s][n].RAWT24,
		TGS[s][n].RAWT25,
		TGS[s][n].RAWT26,
		TGS[s][n].RAWT27,
		TGS[s][n].RAWT28,
		TGS[s][n].RAWT29,
		TGS[s][n].RAWT30,
		TGS[s][n].RAWT31,
		TGS[s][n].RAWT32,
		TGS[s][n].RAWT33,
		TGS[s][n].RAWT34,
		TGS[s][n].RAWT35,
		TGS[s][n].RAWT36,
		TGS[s][n].RAWT37,
		TGS[s][n].RAWT38,
		TGS[s][n].RAWT39,
		TGS[s][n].RAWT40,
		TGS[s][n].RAWT41,
		TGS[s][n].RAWT42,
		TGS[s][n].RAWT43,
		TGS[s][n].RAWT44,
		TGS[s][n].RAWT45,
		TGS[s][n].RAWT46,
		TGS[s][n].RAWT47,
		TGS[s][n].RAWT48,
		TGS[s][n].RAWT49,
		TGS[s][n].RAWT50,
		TGS[s][n].RAWT51,
		TGS[s][n].RAWT52,
		TGS[s][n].RAWT53,
		TGS[s][n].RAWT54,
		TGS[s][n].RAWT55,
		TGS[s][n].RAWT56,
		TGS[s][n].RAWT57,
		TGS[s][n].RAWT58,
		TGS[s][n].RAWT59,
		TGS[s][n].RAWT60,
		TGS[s][n].RAWT61,
		TGS[s][n].notused62,
		TGS[s][n].notused63,
		TGS[s][n].notused64,
		TGS[s][n].BNK1T1,
		TGS[s][n].BNK1T2,
		TGS[s][n].BNK1T3,
		TGS[s][n].BNK1T4,
		TGS[s][n].BNK1T5,
		TGS[s][n].BNK1T6,
		TGS[s][n].BNK1T7,
		TGS[s][n].BNK1T8,
		TGS[s][n].BNK1T9,
		TGS[s][n].BNK1T10,
		TGS[s][n].BNK1T11,
		TGS[s][n].BNK1T12,
		TGS[s][n].BNK1T13,
		TGS[s][n].BNK1T14,
		TGS[s][n].BNK1T15,
		TGS[s][n].BNK1T16,
		TGS[s][n].BNK1T17,
		TGS[s][n].BNK1T18,
		TGS[s][n].BNK1T19,
		TGS[s][n].BNK1T20,
		TGS[s][n].BNK1T21,
		TGS[s][n].BNK1T22,
		TGS[s][n].BNK1T23,
		TGS[s][n].BNK1T24,
		TGS[s][n].BNK1T25,
		TGS[s][n].BNK1T26,
		TGS[s][n].BNK1T27,
		TGS[s][n].BNK1T28,
		TGS[s][n].BNK1T29,
		TGS[s][n].BNK1T30,
		TGS[s][n].BNK1T31,
		TGS[s][n].BNK1T32,
		TGS[s][n].BNK1T33,
		TGS[s][n].BNK1T34,
		TGS[s][n].BNK1T35,
		TGS[s][n].BNK1T36,
		TGS[s][n].BNK1T37,
		TGS[s][n].BNK1T38,
		TGS[s][n].BNK1T39,
		TGS[s][n].BNK1T40,
		TGS[s][n].BNK1T41,
		TGS[s][n].BNK1T42,
		TGS[s][n].BNK1T43,
		TGS[s][n].BNK1T44,
		TGS[s][n].BNK1T45,
		TGS[s][n].BNK1T46,
		TGS[s][n].BNK1T47,
		TGS[s][n].BNK1T48,
		TGS[s][n].BNK1T49,
		TGS[s][n].BNK1T50,
		TGS[s][n].BNK1T51,
		TGS[s][n].BNK1T52,
		TGS[s][n].BNK1T53,
		TGS[s][n].BNK1T54,
		TGS[s][n].BNK1T55,
		TGS[s][n].BNK1T56,
		TGS[s][n].BNK1T57,
		TGS[s][n].BNK1T58,
		TGS[s][n].BNK1T59,
		TGS[s][n].BNK1T60,
		TGS[s][n].BNK1T61,
		TGS[s][n].notused126,
		TGS[s][n].notused127,
		TGS[s][n].notused128,
		TGS[s][n].BNK2T1,
		TGS[s][n].BNK2T2,
		TGS[s][n].BNK2T3,
		TGS[s][n].BNK2T4,
		TGS[s][n].BNK2T5,
		TGS[s][n].BNK2T6,
		TGS[s][n].BNK2T7,
		TGS[s][n].BNK2T8,
		TGS[s][n].BNK2T9,
		TGS[s][n].BNK2T10,
		TGS[s][n].BNK2T11,
		TGS[s][n].BNK2T12,
		TGS[s][n].BNK2T13,
		TGS[s][n].BNK2T14,
		TGS[s][n].BNK2T15,
		TGS[s][n].BNK2T16,
		TGS[s][n].BNK2T17,
		TGS[s][n].BNK2T18,
		TGS[s][n].BNK2T19,
		TGS[s][n].BNK2T20,
		TGS[s][n].BNK2T21,
		TGS[s][n].BNK2T22,
		TGS[s][n].BNK2T23,
		TGS[s][n].BNK2T24,
		TGS[s][n].BNK2T25,
		TGS[s][n].BNK2T26,
		TGS[s][n].BNK2T27,
		TGS[s][n].BNK2T28,
		TGS[s][n].BNK2T29,
		TGS[s][n].BNK2T30,
		TGS[s][n].BNK2T31,
		TGS[s][n].BNK2T32,
		TGS[s][n].BNK2T33,
		TGS[s][n].BNK2T34,
		TGS[s][n].BNK2T35,
		TGS[s][n].BNK2T36,
		TGS[s][n].BNK2T37,
		TGS[s][n].BNK2T38,
		TGS[s][n].BNK2T39,
		TGS[s][n].BNK2T40,
		TGS[s][n].BNK2T41,
		TGS[s][n].BNK2T42,
		TGS[s][n].BNK2T43,
		TGS[s][n].BNK2T44,
		TGS[s][n].BNK2T45,
		TGS[s][n].BNK2T46,
		TGS[s][n].BNK2T47,
		TGS[s][n].BNK2T48,
		TGS[s][n].BNK2T49,
		TGS[s][n].BNK2T50,
		TGS[s][n].BNK2T51,
		TGS[s][n].BNK2T52,
		TGS[s][n].BNK2T53,
		TGS[s][n].BNK2T54,
		TGS[s][n].BNK2T55,
		TGS[s][n].BNK2T56,
		TGS[s][n].BNK2T57,
		TGS[s][n].BNK2T58,
		TGS[s][n].BNK2T59,
		TGS[s][n].BNK2T60,
		TGS[s][n].BNK2T61,
		TGS[s][n].notused190,
		TGS[s][n].notused191,
		TGS[s][n].notused192);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TGTL - single sector bank
struct TGTL_t	 *TGTL;
int 		  TGTL_NH;
int 		  TGTL_S;
void TGTL_Print(){
  int n;
  if(TGTL_NH){
    fprintf(stdout,"Group:\tTGTL\tSector:\t%d \tNhits:\t%d\n",TGTL_S,TGTL_NH);
    for(n=0;n<TGTL_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		TGTL[n].ID,
		TGTL[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TGTR - single sector bank
struct TGTR_t	 *TGTR;
int 		  TGTR_NH;
int 		  TGTR_S;
void TGTR_Print(){
  int n;
  if(TGTR_NH){
    fprintf(stdout,"Group:\tTGTR\tSector:\t%d \tNhits:\t%d\n",TGTR_S,TGTR_NH);
    for(n=0;n<TGTR_NH;n++){
      fprintf(stdout,"ID: %hd\tTDC: %hd\t",
		TGTR[n].ID,
		TGTR[n].TDC);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank TLV1 - multi sector bank
struct TLV1_t	**TLV1;
int 		 *TLV1_NH;
int 		 *TLV1_S;
int 		  TLV1_NS;
void TLV1_Print(){
  int s,n;
  if(TLV1_NS){
    for(s=0;s<TLV1_NS;s++){
      if(TLV1_NH[s]){
        fprintf(stdout,"Group:\tTLV1\tSector:\t%d \tNhits:\t%d\n",TLV1_S[s],TLV1_NH[s]);
        for(n=0;n<TLV1_NH[s];n++){
          fprintf(stdout,"ID: %d\ttime: %8.3f\t",
		TLV1[s][n].ID,
		TLV1[s][n].time);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TRCF - multi sector bank
struct TRCF_t	**TRCF;
int 		 *TRCF_NH;
int 		 *TRCF_S;
int 		  TRCF_NS;
void TRCF_Print(){
  int s,n;
  if(TRCF_NS){
    for(s=0;s<TRCF_NS;s++){
      if(TRCF_NH[s]){
        fprintf(stdout,"Group:\tTRCF\tSector:\t%d \tNhits:\t%d\n",TRCF_S[s],TRCF_NH[s]);
        for(n=0;n<TRCF_NH[s];n++){
          fprintf(stdout,"trigcfg: %.4s\t",
		TRCF[s][n].trigcfg);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TRGS - multi sector bank
struct TRGS_t	**TRGS;
int 		 *TRGS_NH;
int 		 *TRGS_S;
int 		  TRGS_NS;
void TRGS_Print(){
  int s,n;
  if(TRGS_NS){
    for(s=0;s<TRGS_NS;s++){
      if(TRGS_NH[s]){
        fprintf(stdout,"Group:\tTRGS\tSector:\t%d \tNhits:\t%d\n",TRGS_S[s],TRGS_NH[s]);
        for(n=0;n<TRGS_NH[s];n++){
          fprintf(stdout,"CLOCK_UG: %d\tFCUP_UG: %d\tMicrosec: %d\trandom_ug: %d\tMOR_ST: %d\tMOR_PC: %d\tMOR_PS: %d\tMOR_TAC: %d\tMOR: %d\tPC: %d\tPS: %d\tTAC: %d\tST: %d\thel_sync: %d\tclock_ug_2: %d\tfcup_ug_2: %d\tCLOCK_G1: %d\tFCUP_G1: %d\tNotUsed19: %d\trandom_g1: %d\tMOR_ST_rg: %d\tMOR_PC_rg: %d\tMOR_PS_rg: %d\tMOR_TAC_rg: %d\tMOR_rg: %d\tPC_rg: %d\tPS_rg: %d\tTAC_rg: %d\tST_rg: %d\trandom_g1_2: %d\tclock_g1_2: %d\tfcup_g1_2: %d\tCLOCK_G2: %d\tFCUP_G2: %d\ttrig_or_g2: %d\trandom_g2: %d\tNotUsed37: %d\tNotUsed38: %d\tNotUsed39: %d\tNotUsed40: %d\tMOR_lg: %d\tNotUsed42: %d\tNotUsed43: %d\tNotUsed44: %d\tNotUsed45: %d\trandom_g2_2: %d\tclock_g2_2: %d\tfcup_g2_2: %d\ttrig1_ug: %d\ttrig2_ug: %d\ttrig3_ug: %d\ttrig4_ug: %d\ttrig5_ug: %d\ttrig6_ug: %d\ttrig7_ug: %d\ttrig8_ug: %d\ttrig9_ug: %d\ttrig10_ug: %d\ttrig11_ug: %d\ttrig12_ug: %d\ttrig_or_ug: %d\tl1accept: %d\tnotused63: %d\tnotused64: %d\tl2fail: %d\tl2pass: %d\tl2start: %d\tl2clear: %d\tl2accept: %d\tl3accept: %d\tnotused71: %d\tnotused72: %d\tl2sec1_g: %d\tl2sec2_g: %d\tl2sec3_g: %d\tl2sec4_g: %d\tl2sec5_g: %d\tl2sec6_g: %d\tl2_or_g: %d\tl2_ok_g: %d\ttrig1_lg: %d\ttrig2_lg: %d\ttrig3_lg: %d\ttrig4_lg: %d\ttrig5_lg: %d\ttrig6_lg: %d\ttrig7_lg: %d\ttrig8_lg: %d\ttrig9_lg: %d\ttrig10_lg: %d\ttrig11_lg: %d\ttrig12_lg: %d\tnotused93: %d\tnotused94: %d\tignore95: %d\tignore96: %d\t",
		TRGS[s][n].CLOCK_UG,
		TRGS[s][n].FCUP_UG,
		TRGS[s][n].Microsec,
		TRGS[s][n].random_ug,
		TRGS[s][n].MOR_ST,
		TRGS[s][n].MOR_PC,
		TRGS[s][n].MOR_PS,
		TRGS[s][n].MOR_TAC,
		TRGS[s][n].MOR,
		TRGS[s][n].PC,
		TRGS[s][n].PS,
		TRGS[s][n].TAC,
		TRGS[s][n].ST,
		TRGS[s][n].hel_sync,
		TRGS[s][n].clock_ug_2,
		TRGS[s][n].fcup_ug_2,
		TRGS[s][n].CLOCK_G1,
		TRGS[s][n].FCUP_G1,
		TRGS[s][n].NotUsed19,
		TRGS[s][n].random_g1,
		TRGS[s][n].MOR_ST_rg,
		TRGS[s][n].MOR_PC_rg,
		TRGS[s][n].MOR_PS_rg,
		TRGS[s][n].MOR_TAC_rg,
		TRGS[s][n].MOR_rg,
		TRGS[s][n].PC_rg,
		TRGS[s][n].PS_rg,
		TRGS[s][n].TAC_rg,
		TRGS[s][n].ST_rg,
		TRGS[s][n].random_g1_2,
		TRGS[s][n].clock_g1_2,
		TRGS[s][n].fcup_g1_2,
		TRGS[s][n].CLOCK_G2,
		TRGS[s][n].FCUP_G2,
		TRGS[s][n].trig_or_g2,
		TRGS[s][n].random_g2,
		TRGS[s][n].NotUsed37,
		TRGS[s][n].NotUsed38,
		TRGS[s][n].NotUsed39,
		TRGS[s][n].NotUsed40,
		TRGS[s][n].MOR_lg,
		TRGS[s][n].NotUsed42,
		TRGS[s][n].NotUsed43,
		TRGS[s][n].NotUsed44,
		TRGS[s][n].NotUsed45,
		TRGS[s][n].random_g2_2,
		TRGS[s][n].clock_g2_2,
		TRGS[s][n].fcup_g2_2,
		TRGS[s][n].trig1_ug,
		TRGS[s][n].trig2_ug,
		TRGS[s][n].trig3_ug,
		TRGS[s][n].trig4_ug,
		TRGS[s][n].trig5_ug,
		TRGS[s][n].trig6_ug,
		TRGS[s][n].trig7_ug,
		TRGS[s][n].trig8_ug,
		TRGS[s][n].trig9_ug,
		TRGS[s][n].trig10_ug,
		TRGS[s][n].trig11_ug,
		TRGS[s][n].trig12_ug,
		TRGS[s][n].trig_or_ug,
		TRGS[s][n].l1accept,
		TRGS[s][n].notused63,
		TRGS[s][n].notused64,
		TRGS[s][n].l2fail,
		TRGS[s][n].l2pass,
		TRGS[s][n].l2start,
		TRGS[s][n].l2clear,
		TRGS[s][n].l2accept,
		TRGS[s][n].l3accept,
		TRGS[s][n].notused71,
		TRGS[s][n].notused72,
		TRGS[s][n].l2sec1_g,
		TRGS[s][n].l2sec2_g,
		TRGS[s][n].l2sec3_g,
		TRGS[s][n].l2sec4_g,
		TRGS[s][n].l2sec5_g,
		TRGS[s][n].l2sec6_g,
		TRGS[s][n].l2_or_g,
		TRGS[s][n].l2_ok_g,
		TRGS[s][n].trig1_lg,
		TRGS[s][n].trig2_lg,
		TRGS[s][n].trig3_lg,
		TRGS[s][n].trig4_lg,
		TRGS[s][n].trig5_lg,
		TRGS[s][n].trig6_lg,
		TRGS[s][n].trig7_lg,
		TRGS[s][n].trig8_lg,
		TRGS[s][n].trig9_lg,
		TRGS[s][n].trig10_lg,
		TRGS[s][n].trig11_lg,
		TRGS[s][n].trig12_lg,
		TRGS[s][n].notused93,
		TRGS[s][n].notused94,
		TRGS[s][n].ignore95,
		TRGS[s][n].ignore96);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TRKS - multi sector bank
struct TRKS_t	**TRKS;
int 		 *TRKS_NH;
int 		 *TRKS_S;
int 		  TRKS_NS;
void TRKS_Print(){
  int s,n;
  if(TRKS_NS){
    for(s=0;s<TRKS_NS;s++){
      if(TRKS_NH[s]){
        fprintf(stdout,"Group:\tTRKS\tSector:\t%d \tNhits:\t%d\n",TRKS_S[s],TRKS_NH[s]);
        for(n=0;n<TRKS_NH[s];n++){
          fprintf(stdout,"trk_lev: %d\tbeta: %8.3f\tst_time: %8.3f\tcc_time: %8.3f\tsc_time: %8.3f\tec_time: %8.3f\tst_id: %d\tcc_id: %d\tsc_id: %d\tec_id: %d\t",
		TRKS[s][n].trk_lev,
		TRKS[s][n].beta,
		TRKS[s][n].st_time,
		TRKS[s][n].cc_time,
		TRKS[s][n].sc_time,
		TRKS[s][n].ec_time,
		TRKS[s][n].st_id,
		TRKS[s][n].cc_id,
		TRKS[s][n].sc_id,
		TRKS[s][n].ec_id);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TRL1 - multi sector bank
struct TRL1_t	**TRL1;
int 		 *TRL1_NH;
int 		 *TRL1_S;
int 		  TRL1_NS;
void TRL1_Print(){
  int s,n;
  if(TRL1_NS){
    for(s=0;s<TRL1_NS;s++){
      if(TRL1_NH[s]){
        fprintf(stdout,"Group:\tTRL1\tSector:\t%d \tNhits:\t%d\n",TRL1_S[s],TRL1_NH[s]);
        for(n=0;n<TRL1_NH[s];n++){
          fprintf(stdout,"fit_flags: %d\thits_hbt: %d\thits_tbt: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\ttlen: %8.3f\t",
		TRL1[s][n].fit_flags,
		TRL1[s][n].hits_hbt,
		TRL1[s][n].hits_tbt,
		TRL1[s][n].x,
		TRL1[s][n].y,
		TRL1[s][n].z,
		TRL1[s][n].cx,
		TRL1[s][n].cy,
		TRL1[s][n].cz,
		TRL1[s][n].tlen);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TRPB - multi sector bank
struct TRPB_t	**TRPB;
int 		 *TRPB_NH;
int 		 *TRPB_S;
int 		  TRPB_NS;
void TRPB_Print(){
  int s,n;
  if(TRPB_NS){
    for(s=0;s<TRPB_NS;s++){
      if(TRPB_NH[s]){
        fprintf(stdout,"Group:\tTRPB\tSector:\t%d \tNhits:\t%d\n",TRPB_S[s],TRPB_NH[s]);
        for(n=0;n<TRPB_NH[s];n++){
          fprintf(stdout,"hits_tbt: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tcx: %8.3f\tcy: %8.3f\tcz: %8.3f\ttlen: %8.3f\t",
		TRPB[s][n].hits_tbt,
		TRPB[s][n].x,
		TRPB[s][n].y,
		TRPB[s][n].z,
		TRPB[s][n].cx,
		TRPB[s][n].cy,
		TRPB[s][n].cz,
		TRPB[s][n].tlen);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TSPR - multi sector bank
struct TSPR_t	**TSPR;
int 		 *TSPR_NH;
int 		 *TSPR_S;
int 		  TSPR_NS;
void TSPR_Print(){
  int s,n;
  if(TSPR_NS){
    for(s=0;s<TSPR_NS;s++){
      if(TSPR_NH[s]){
        fprintf(stdout,"Group:\tTSPR\tSector:\t%d \tNhits:\t%d\n",TSPR_S[s],TSPR_NH[s]);
        for(n=0;n<TSPR_NH[s];n++){
          fprintf(stdout,"mem: %d\t",
		TSPR[s][n].mem);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank TSRG - multi sector bank
struct TSRG_t	**TSRG;
int 		 *TSRG_NH;
int 		 *TSRG_S;
int 		  TSRG_NS;
void TSRG_Print(){
  int s,n;
  if(TSRG_NS){
    for(s=0;s<TSRG_NS;s++){
      if(TSRG_NH[s]){
        fprintf(stdout,"Group:\tTSRG\tSector:\t%d \tNhits:\t%d\n",TSRG_S[s],TSRG_NH[s]);
        for(n=0;n<TSRG_NH[s];n++){
          fprintf(stdout,"reg: %d\t",
		TSRG[s][n].reg);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank UNUS - multi sector bank
struct UNUS_t	**UNUS;
int 		 *UNUS_NH;
int 		 *UNUS_S;
int 		  UNUS_NS;
void UNUS_Print(){
  int s,n;
  if(UNUS_NS){
    for(s=0;s<UNUS_NS;s++){
      if(UNUS_NH[s]){
        fprintf(stdout,"Group:\tUNUS\tSector:\t%d \tNhits:\t%d\n",UNUS_S[s],UNUS_NH[s]);
        for(n=0;n<UNUS_NH[s];n++){
          fprintf(stdout,"NDCUN: %d\tIDCUN: %d\tNSCUN: %d\tISCUN: %d\tNCCUN: %d\tICCUN: %d\tNECUN: %d\tIECUN: %d\t",
		UNUS[s][n].NDCUN,
		UNUS[s][n].IDCUN,
		UNUS[s][n].NSCUN,
		UNUS[s][n].ISCUN,
		UNUS[s][n].NCCUN,
		UNUS[s][n].ICCUN,
		UNUS[s][n].NECUN,
		UNUS[s][n].IECUN);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank VERT - multi sector bank
struct VERT_t	**VERT;
int 		 *VERT_NH;
int 		 *VERT_S;
int 		  VERT_NS;
void VERT_Print(){
  int s,n;
  if(VERT_NS){
    for(s=0;s<VERT_NS;s++){
      if(VERT_NH[s]){
        fprintf(stdout,"Group:\tVERT\tSector:\t%d \tNhits:\t%d\n",VERT_S[s],VERT_NH[s]);
        for(n=0;n<VERT_NH[s];n++){
          fprintf(stdout,"vertex: %d\ttrk1: %d\ttrk2: %d\tx: %8.3f\ty: %8.3f\tz: %8.3f\tsepd: %8.3f\t",
		VERT[s][n].vertex,
		VERT[s][n].trk1,
		VERT[s][n].trk2,
		VERT[s][n].x,
		VERT[s][n].y,
		VERT[s][n].z,
		VERT[s][n].sepd);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

// bank RC26 - single sector bank
struct RC26_t	 *RC26;
int 		  RC26_NH;
int 		  RC26_S;
void RC26_Print(){
  int n;
  if(RC26_NH){
    fprintf(stdout,"Group:\tRC26\tSector:\t%d \tNhits:\t%d\n",RC26_S,RC26_NH);
    for(n=0;n<RC26_NH;n++){
      fprintf(stdout,"trigger_pattern: %d\tgated_scaler: %d\tungated_scaler: %d\t",
		RC26[n].trigger_pattern,
		RC26[n].gated_scaler,
		RC26[n].ungated_scaler);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }
}

// bank STRE - multi sector bank
struct STRE_t	**STRE;
int 		 *STRE_NH;
int 		 *STRE_S;
int 		  STRE_NS;
void STRE_Print(){
  int s,n;
  if(STRE_NS){
    for(s=0;s<STRE_NS;s++){
      if(STRE_NH[s]){
        fprintf(stdout,"Group:\tSTRE\tSector:\t%d \tNhits:\t%d\n",STRE_S[s],STRE_NH[s]);
        for(n=0;n<STRE_NH[s];n++){
          fprintf(stdout,"ID: %d\tTrk_no: %d\tST_TIME: %8.3f\tST_L: %8.3f\tst_pos: %8.3f\tstatus: %d\tst_edep: %8.3f\tst_corr: %8.3f\t",
		STRE[s][n].ID,
		STRE[s][n].Trk_no,
		STRE[s][n].ST_TIME,
		STRE[s][n].ST_L,
		STRE[s][n].st_pos,
		STRE[s][n].status,
		STRE[s][n].st_edep,
		STRE[s][n].st_corr);
          fprintf(stdout,"\n");
        }
        fprintf(stdout,"\n");
      }
    }
  }
}

int 		  nBankTypes = 196;	// Total of different banks 

// Make a table with the bank names, and addresses of the structs and variables
// associated with them
struct addressBanks_t bankAdd[196] = {
	{"BEAM", (void**)&BEAM, &BEAM_NS, &BEAM_NH, &BEAM_S, 16, BEAM_Print},
	{"BMPR", (void**)&BMPR, NULL, (int**)&BMPR_NH, (int**)&BMPR_S, 16, BMPR_Print},
	{"BREP", (void**)&BREP, NULL, (int**)&BREP_NH, (int**)&BREP_S, 144, BREP_Print},
	{"CALB", (void**)&CALB, &CALB_NS, &CALB_NH, &CALB_S, 48, CALB_Print},
	{"CALL", (void**)&CALL, NULL, (int**)&CALL_NH, (int**)&CALL_S, 6, CALL_Print},
	{"CC01", (void**)&CC01, &CC01_NS, &CC01_NH, &CC01_S, 12, CC01_Print},
	{"CC1 ", (void**)&CC1, &CC1_NS, &CC1_NH, &CC1_S, 6, CC1_Print},
	{"CC  ", (void**)&CC, &CC_NS, &CC_NH, &CC_S, 6, CC_Print},
	{"CCDI", (void**)&CCDI, &CCDI_NS, &CCDI_NH, &CCDI_S, 20, CCDI_Print},
	{"CCH ", (void**)&CCH, &CCH_NS, &CCH_NH, &CCH_S, 48, CCH_Print},
	{"CCMT", (void**)&CCMT, &CCMT_NS, &CCMT_NH, &CCMT_S, 8, CCMT_Print},
	{"CCPB", (void**)&CCPB, NULL, (int**)&CCPB_NH, (int**)&CCPB_S, 24, CCPB_Print},
	{"CCPE", (void**)&CCPE, &CCPE_NS, &CCPE_NH, &CCPE_S, 12, CCPE_Print},
	{"CCRC", (void**)&CCRC, NULL, (int**)&CCRC_NH, (int**)&CCRC_S, 48, CCRC_Print},
	{"CCS ", (void**)&CCS, &CCS_NS, &CCS_NH, &CCS_S, 384, CCS_Print},
	{"CCT ", (void**)&CCT, &CCT_NS, &CCT_NH, &CCT_S, 4, CCT_Print},
	{"CHI2", (void**)&CHI2, &CHI2_NS, &CHI2_NH, &CHI2_S, 16, CHI2_Print},
	{"CL01", (void**)&CL01, NULL, (int**)&CL01_NH, (int**)&CL01_S, 48, CL01_Print},
	{"CLST", (void**)&CLST, &CLST_NS, &CLST_NH, &CLST_S, 4, CLST_Print},
	{"CPED", (void**)&CPED, &CPED_NS, &CPED_NH, &CPED_S, 20, CPED_Print},
	{"CSQL", (void**)&CSQL, &CSQL_NS, &CSQL_NH, &CSQL_S, 84, CSQL_Print},
	{"DC0 ", (void**)&DC0, &DC0_NS, &DC0_NH, &DC0_S, 4, DC0_Print},
	{"DC1 ", (void**)&DC1, &DC1_NS, &DC1_NH, &DC1_S, 8, DC1_Print},
	{"DCDW", (void**)&DCDW, &DCDW_NS, &DCDW_NH, &DCDW_S, 12, DCDW_Print},
	{"DCGM", (void**)&DCGM, &DCGM_NS, &DCGM_NH, &DCGM_S, 84, DCGM_Print},
	{"DCGW", (void**)&DCGW, &DCGW_NS, &DCGW_NH, &DCGW_S, 32, DCGW_Print},
	{"DCH ", (void**)&DCH, &DCH_NS, &DCH_NH, &DCH_S, 52, DCH_Print},
	{"DCMN", (void**)&DCMN, &DCMN_NS, &DCMN_NH, &DCMN_S, 200, DCMN_Print},
	{"DCPB", (void**)&DCPB, NULL, (int**)&DCPB_NH, (int**)&DCPB_S, 52, DCPB_Print},
	{"DCV1", (void**)&DCV1, &DCV1_NS, &DCV1_NH, &DCV1_S, 20, DCV1_Print},
	{"DCV2", (void**)&DCV2, &DCV2_NS, &DCV2_NH, &DCV2_S, 56, DCV2_Print},
	{"DCV3", (void**)&DCV3, &DCV3_NS, &DCV3_NH, &DCV3_S, 20, DCV3_Print},
	{"DDLY", (void**)&DDLY, &DDLY_NS, &DDLY_NH, &DDLY_S, 12, DDLY_Print},
	{"DGEO", (void**)&DGEO, &DGEO_NS, &DGEO_NH, &DGEO_S, 32, DGEO_Print},
	{"DHCL", (void**)&DHCL, &DHCL_NS, &DHCL_NH, &DHCL_S, 64, DHCL_Print},
	{"DITM", (void**)&DITM, &DITM_NS, &DITM_NH, &DITM_S, 4, DITM_Print},
	{"DOCA", (void**)&DOCA, &DOCA_NS, &DOCA_NH, &DOCA_S, 4, DOCA_Print},
	{"DPCP", (void**)&DPCP, &DPCP_NS, &DPCP_NH, &DPCP_S, 52, DPCP_Print},
	{"DPSP", (void**)&DPSP, &DPSP_NS, &DPSP_NH, &DPSP_S, 12, DPSP_Print},
	{"DSPC", (void**)&DSPC, &DSPC_NS, &DSPC_NH, &DSPC_S, 16, DSPC_Print},
	{"DSPS", (void**)&DSPS, &DSPS_NS, &DSPS_NH, &DSPS_S, 6, DSPS_Print},
	{"DSTC", (void**)&DSTC, NULL, (int**)&DSTC_NH, (int**)&DSTC_S, 18, DSTC_Print},
	{"DTCP", (void**)&DTCP, &DTCP_NS, &DTCP_NH, &DTCP_S, 60, DTCP_Print},
	{"DTRK", (void**)&DTRK, &DTRK_NS, &DTRK_NH, &DTRK_S, 12, DTRK_Print},
	{"EC01", (void**)&EC01, &EC01_NS, &EC01_NH, &EC01_S, 12, EC01_Print},
	{"EC1 ", (void**)&EC1, &EC1_NS, &EC1_NH, &EC1_S, 10, EC1_Print},
	{"EC1P", (void**)&EC1P, &EC1P_NS, &EC1P_NH, &EC1P_S, 16, EC1P_Print},
	{"EC1R", (void**)&EC1R, &EC1R_NS, &EC1R_NH, &EC1R_S, 200, EC1R_Print},
	{"ECCA", (void**)&ECCA, &ECCA_NS, &ECCA_NH, &ECCA_S, 32, ECCA_Print},
	{"ECCL", (void**)&ECCL, &ECCL_NS, &ECCL_NH, &ECCL_S, 32, ECCL_Print},
	{"ECCT", (void**)&ECCT, &ECCT_NS, &ECCT_NH, &ECCT_S, 48, ECCT_Print},
	{"EC  ", (void**)&EC, &EC_NS, &EC_NH, &EC_S, 6, EC_Print},
	{"ECDI", (void**)&ECDI, &ECDI_NS, &ECDI_NH, &ECDI_S, 20, ECDI_Print},
	{"ECG ", (void**)&ECG, &ECG_NS, &ECG_NH, &ECG_S, 88, ECG_Print},
	{"ECHB", (void**)&ECHB, NULL, (int**)&ECHB_NH, (int**)&ECHB_S, 152, ECHB_Print},
	{"ECH ", (void**)&ECH, &ECH_NS, &ECH_NH, &ECH_S, 44, ECH_Print},
	{"ECMT", (void**)&ECMT, &ECMT_NS, &ECMT_NH, &ECMT_S, 24, ECMT_Print},
	{"ECP1", (void**)&ECP1, &ECP1_NS, &ECP1_NH, &ECP1_S, 20, ECP1_Print},
	{"ECPB", (void**)&ECPB, NULL, (int**)&ECPB_NH, (int**)&ECPB_S, 64, ECPB_Print},
	{"ECPC", (void**)&ECPC, &ECPC_NS, &ECPC_NH, &ECPC_S, 12, ECPC_Print},
	{"ECP ", (void**)&ECP, &ECP_NS, &ECP_NH, &ECP_S, 28, ECP_Print},
	{"ECPE", (void**)&ECPE, &ECPE_NS, &ECPE_NH, &ECPE_S, 12, ECPE_Print},
	{"ECPI", (void**)&ECPI, NULL, (int**)&ECPI_NH, (int**)&ECPI_S, 36, ECPI_Print},
	{"ECPO", (void**)&ECPO, &ECPO_NS, &ECPO_NH, &ECPO_S, 28, ECPO_Print},
	{"ECRB", (void**)&ECRB, &ECRB_NS, &ECRB_NH, &ECRB_S, 196, ECRB_Print},
	{"ECS ", (void**)&ECS, &ECS_NS, &ECS_NH, &ECS_S, 384, ECS_Print},
	{"ECT ", (void**)&ECT, &ECT_NS, &ECT_NH, &ECT_S, 4, ECT_Print},
	{"EID0", (void**)&EID0, &EID0_NS, &EID0_NH, &EID0_S, 16, EID0_Print},
	{"EPIC", (void**)&EPIC, &EPIC_NS, &EPIC_NH, &EPIC_S, 36, EPIC_Print},
	{"EVNT", (void**)&EVNT, NULL, (int**)&EVNT_NH, (int**)&EVNT_S, 72, EVNT_Print},
	{"FBPM", (void**)&FBPM, &FBPM_NS, &FBPM_NH, &FBPM_S, 6, FBPM_Print},
	{"G1SL", (void**)&G1SL, &G1SL_NS, &G1SL_NH, &G1SL_S, 256, G1SL_Print},
	{"G2SL", (void**)&G2SL, &G2SL_NS, &G2SL_NH, &G2SL_S, 256, G2SL_Print},
	{"G3SL", (void**)&G3SL, &G3SL_NS, &G3SL_NH, &G3SL_S, 256, G3SL_Print},
	{"G4SL", (void**)&G4SL, &G4SL_NS, &G4SL_NH, &G4SL_S, 256, G4SL_Print},
	{"GPAR", (void**)&GPAR, &GPAR_NS, &GPAR_NH, &GPAR_S, 28, GPAR_Print},
	{"GPID", (void**)&GPID, NULL, (int**)&GPID_NH, (int**)&GPID_S, 112, GPID_Print},
	{"GP_X", (void**)&GP_X, NULL, (int**)&GP_X_NH, (int**)&GP_X_S, 4, GP_X_Print},
	{"GP_Y", (void**)&GP_Y, NULL, (int**)&GP_Y_NH, (int**)&GP_Y_S, 4, GP_Y_Print},
	{"HBER", (void**)&HBER, NULL, (int**)&HBER_NH, (int**)&HBER_S, 92, HBER_Print},
	{"HBID", (void**)&HBID, NULL, (int**)&HBID_NH, (int**)&HBID_S, 136, HBID_Print},
	{"HBLA", (void**)&HBLA, &HBLA_NS, &HBLA_NH, &HBLA_S, 64, HBLA_Print},
	{"HBTB", (void**)&HBTB, &HBTB_NS, &HBTB_NH, &HBTB_S, 40, HBTB_Print},
	{"HBTR", (void**)&HBTR, NULL, (int**)&HBTR_NH, (int**)&HBTR_S, 36, HBTR_Print},
	{"HCAL", (void**)&HCAL, &HCAL_NS, &HCAL_NH, &HCAL_S, 20, HCAL_Print},
	{"HDPL", (void**)&HDPL, &HDPL_NS, &HDPL_NH, &HDPL_S, 32, HDPL_Print},
	{"HEAD", (void**)&HEAD, NULL, (int**)&HEAD_NH, (int**)&HEAD_S, 32, HEAD_Print},
	{"HEVT", (void**)&HEVT, NULL, (int**)&HEVT_NH, (int**)&HEVT_S, 64, HEVT_Print},
	{"HLS ", (void**)&HLS, &HLS_NS, &HLS_NH, &HLS_S, 64, HLS_Print},
	{"IC  ", (void**)&IC, &IC_NS, &IC_NH, &IC_S, 6, IC_Print},
	{"ICHB", (void**)&ICHB, &ICHB_NS, &ICHB_NH, &ICHB_S, 56, ICHB_Print},
	{"KFIT", (void**)&KFIT, &KFIT_NS, &KFIT_NH, &KFIT_S, 184, KFIT_Print},
	{"L1PG", (void**)&L1PG, &L1PG_NS, &L1PG_NH, &L1PG_S, 4, L1PG_Print},
	{"L2H ", (void**)&L2H, &L2H_NS, &L2H_NH, &L2H_S, 4, L2H_Print},
	{"L2S ", (void**)&L2S, &L2S_NS, &L2S_NH, &L2S_S, 4, L2S_Print},
	{"LASR", (void**)&LASR, &LASR_NS, &LASR_NH, &LASR_S, 16, LASR_Print},
	{"LCDI", (void**)&LCDI, &LCDI_NS, &LCDI_NH, &LCDI_S, 16, LCDI_Print},
	{"LCPB", (void**)&LCPB, NULL, (int**)&LCPB_NH, (int**)&LCPB_S, 40, LCPB_Print},
	{"LOWQ", (void**)&LOWQ, &LOWQ_NS, &LOWQ_NH, &LOWQ_S, 48, LOWQ_Print},
	{"MCEV", (void**)&MCEV, &MCEV_NS, &MCEV_NH, &MCEV_S, 8, MCEV_Print},
	{"MCHD", (void**)&MCHD, &MCHD_NS, &MCHD_NH, &MCHD_S, 64, MCHD_Print},
	{"MCTK", (void**)&MCTK, &MCTK_NS, &MCTK_NH, &MCTK_S, 44, MCTK_Print},
	{"MCVX", (void**)&MCVX, &MCVX_NS, &MCVX_NH, &MCVX_S, 20, MCVX_Print},
	{"MS1 ", (void**)&MS1, &MS1_NS, &MS1_NH, &MS1_S, 4, MS1_Print},
	{"MTRK", (void**)&MTRK, &MTRK_NS, &MTRK_NH, &MTRK_S, 28, MTRK_Print},
	{"MVRT", (void**)&MVRT, NULL, (int**)&MVRT_NH, (int**)&MVRT_S, 52, MVRT_Print},
	{"PART", (void**)&PART, NULL, (int**)&PART_NH, (int**)&PART_S, 52, PART_Print},
	{"PCO ", (void**)&PCO, &PCO_NS, &PCO_NH, &PCO_S, 32, PCO_Print},
	{"PHTM", (void**)&PHTM, &PHTM_NS, &PHTM_NH, &PHTM_S, 44, PHTM_Print},
	{"PID1", (void**)&PID1, &PID1_NS, &PID1_NH, &PID1_S, 64, PID1_Print},
	{"PID2", (void**)&PID2, &PID2_NS, &PID2_NH, &PID2_S, 64, PID2_Print},
	{"PID3", (void**)&PID3, &PID3_NS, &PID3_NH, &PID3_S, 64, PID3_Print},
	{"PID4", (void**)&PID4, &PID4_NS, &PID4_NH, &PID4_S, 64, PID4_Print},
	{"PID5", (void**)&PID5, &PID5_NS, &PID5_NH, &PID5_S, 64, PID5_Print},
	{"PID6", (void**)&PID6, &PID6_NS, &PID6_NH, &PID6_S, 64, PID6_Print},
	{"PIDT", (void**)&PIDT, &PIDT_NS, &PIDT_NH, &PIDT_S, 40, PIDT_Print},
	{"PRTM", (void**)&PRTM, &PRTM_NS, &PRTM_NH, &PRTM_S, 4, PRTM_Print},
	{"PSO ", (void**)&PSO, &PSO_NS, &PSO_NH, &PSO_S, 16, PSO_Print},
	{"PTDB", (void**)&PTDB, &PTDB_NS, &PTDB_NH, &PTDB_S, 44, PTDB_Print},
	{"RCST", (void**)&RCST, &RCST_NS, &RCST_NH, &RCST_S, 60, RCST_Print},
	{"REF ", (void**)&REF, NULL, (int**)&REF_NH, (int**)&REF_S, 4, REF_Print},
	{"RF  ", (void**)&RF, &RF_NS, &RF_NH, &RF_S, 12, RF_Print},
	{"RFT ", (void**)&RFT, NULL, (int**)&RFT_NH, (int**)&RFT_S, 4, RFT_Print},
	{"RGLK", (void**)&RGLK, &RGLK_NS, &RGLK_NH, &RGLK_S, 40, RGLK_Print},
	{"RNLG", (void**)&RNLG, &RNLG_NS, &RNLG_NH, &RNLG_S, 4, RNLG_Print},
	{"RNPE", (void**)&RNPE, &RNPE_NS, &RNPE_NH, &RNPE_S, 8, RNPE_Print},
	{"RTSL", (void**)&RTSL, &RTSL_NS, &RTSL_NH, &RTSL_S, 256, RTSL_Print},
	{"RUNC", (void**)&RUNC, &RUNC_NS, &RUNC_NH, &RUNC_S, 52, RUNC_Print},
	{"S1ST", (void**)&S1ST, &S1ST_NS, &S1ST_NH, &S1ST_S, 140, S1ST_Print},
	{"SC1 ", (void**)&SC1, &SC1_NS, &SC1_NH, &SC1_S, 36, SC1_Print},
	{"SCC ", (void**)&SCC, &SCC_NS, &SCC_NH, &SCC_S, 192, SCC_Print},
	{"SC  ", (void**)&SC, &SC_NS, &SC_NH, &SC_S, 10, SC_Print},
	{"SCDI", (void**)&SCDI, &SCDI_NS, &SCDI_NH, &SCDI_S, 20, SCDI_Print},
	{"SCGD", (void**)&SCGD, &SCGD_NS, &SCGD_NH, &SCGD_S, 40, SCGD_Print},
	{"SCG ", (void**)&SCG, &SCG_NS, &SCG_NH, &SCG_S, 44, SCG_Print},
	{"SCH ", (void**)&SCH, &SCH_NS, &SCH_NH, &SCH_S, 36, SCH_Print},
	{"SCMD", (void**)&SCMD, &SCMD_NS, &SCMD_NH, &SCMD_S, 52, SCMD_Print},
	{"SCMT", (void**)&SCMT, &SCMT_NS, &SCMT_NH, &SCMT_S, 4, SCMT_Print},
	{"SCMW", (void**)&SCMW, &SCMW_NS, &SCMW_NH, &SCMW_S, 4, SCMW_Print},
	{"SCPB", (void**)&SCPB, NULL, (int**)&SCPB_NH, (int**)&SCPB_S, 24, SCPB_Print},
	{"SCP ", (void**)&SCP, &SCP_NS, &SCP_NH, &SCP_S, 64, SCP_Print},
	{"SCPE", (void**)&SCPE, &SCPE_NS, &SCPE_NH, &SCPE_S, 20, SCPE_Print},
	{"SCPS", (void**)&SCPS, &SCPS_NS, &SCPS_NH, &SCPS_S, 64, SCPS_Print},
	{"SCRC", (void**)&SCRC, &SCRC_NS, &SCRC_NH, &SCRC_S, 48, SCRC_Print},
	{"SCR ", (void**)&SCR, &SCR_NS, &SCR_NH, &SCR_S, 48, SCR_Print},
	{"SCS ", (void**)&SCS, &SCS_NS, &SCS_NH, &SCS_S, 768, SCS_Print},
	{"SCT ", (void**)&SCT, &SCT_NS, &SCT_NH, &SCT_S, 4, SCT_Print},
	{"SGMP", (void**)&SGMP, &SGMP_NS, &SGMP_NH, &SGMP_S, 28, SGMP_Print},
	{"SPAR", (void**)&SPAR, &SPAR_NS, &SPAR_NH, &SPAR_S, 16, SPAR_Print},
	{"SPIN", (void**)&SPIN, &SPIN_NS, &SPIN_NH, &SPIN_S, 12, SPIN_Print},
	{"ST1 ", (void**)&ST1, &ST1_NS, &ST1_NH, &ST1_S, 24, ST1_Print},
	{"ST  ", (void**)&ST, NULL, (int**)&ST_NH, (int**)&ST_S, 6, ST_Print},
	{"STG ", (void**)&STG, &STG_NS, &STG_NH, &STG_S, 64, STG_Print},
	{"STH ", (void**)&STH, &STH_NS, &STH_NH, &STH_S, 40, STH_Print},
	{"STN0", (void**)&STN0, NULL, (int**)&STN0_NH, (int**)&STN0_S, 4, STN0_Print},
	{"STN1", (void**)&STN1, NULL, (int**)&STN1_NH, (int**)&STN1_S, 4, STN1_Print},
	{"STPB", (void**)&STPB, NULL, (int**)&STPB_NH, (int**)&STPB_S, 20, STPB_Print},
	{"STPE", (void**)&STPE, &STPE_NS, &STPE_NH, &STPE_S, 12, STPE_Print},
	{"STR ", (void**)&STR, NULL, (int**)&STR_NH, (int**)&STR_S, 24, STR_Print},
	{"STS ", (void**)&STS, &STS_NS, &STS_NH, &STS_S, 64, STS_Print},
	{"STSN", (void**)&STSN, &STSN_NS, &STSN_NH, &STSN_S, 128, STSN_Print},
	{"SYNC", (void**)&SYNC, &SYNC_NS, &SYNC_NH, &SYNC_S, 6, SYNC_Print},
	{"TACO", (void**)&TACO, &TACO_NS, &TACO_NH, &TACO_S, 40, TACO_Print},
	{"TAGE", (void**)&TAGE, NULL, (int**)&TAGE_NH, (int**)&TAGE_S, 4, TAGE_Print},
	{"TAGI", (void**)&TAGI, &TAGI_NS, &TAGI_NH, &TAGI_S, 32, TAGI_Print},
	{"TAGM", (void**)&TAGM, NULL, (int**)&TAGM_NH, (int**)&TAGM_S, 24, TAGM_Print},
	{"TAGR", (void**)&TAGR, NULL, (int**)&TAGR_NH, (int**)&TAGR_S, 24, TAGR_Print},
	{"TAGT", (void**)&TAGT, NULL, (int**)&TAGT_NH, (int**)&TAGT_S, 6, TAGT_Print},
	{"TATL", (void**)&TATL, &TATL_NS, &TATL_NH, &TATL_S, 8, TATL_Print},
	{"TATR", (void**)&TATR, &TATR_NS, &TATR_NH, &TATR_S, 8, TATR_Print},
	{"TBER", (void**)&TBER, NULL, (int**)&TBER_NH, (int**)&TBER_S, 92, TBER_Print},
	{"TBID", (void**)&TBID, NULL, (int**)&TBID_NH, (int**)&TBID_S, 136, TBID_Print},
	{"TBLA", (void**)&TBLA, &TBLA_NS, &TBLA_NH, &TBLA_S, 68, TBLA_Print},
	{"TBTR", (void**)&TBTR, NULL, (int**)&TBTR_NH, (int**)&TBTR_S, 40, TBTR_Print},
	{"TCSB", (void**)&TCSB, &TCSB_NS, &TCSB_NH, &TCSB_S, 24, TCSB_Print},
	{"TCT ", (void**)&TCT, &TCT_NS, &TCT_NH, &TCT_S, 4, TCT_Print},
	{"TDPL", (void**)&TDPL, &TDPL_NS, &TDPL_NH, &TDPL_S, 32, TDPL_Print},
	{"TESC", (void**)&TESC, NULL, (int**)&TESC_NH, (int**)&TESC_S, 4, TESC_Print},
	{"TGBI", (void**)&TGBI, NULL, (int**)&TGBI_NH, (int**)&TGBI_S, 20, TGBI_Print},
	{"TGEO", (void**)&TGEO, &TGEO_NS, &TGEO_NH, &TGEO_S, 24, TGEO_Print},
	{"TGPB", (void**)&TGPB, NULL, (int**)&TGPB_NH, (int**)&TGPB_S, 16, TGPB_Print},
	{"TGS ", (void**)&TGS, &TGS_NS, &TGS_NH, &TGS_S, 768, TGS_Print},
	{"TGTL", (void**)&TGTL, NULL, (int**)&TGTL_NH, (int**)&TGTL_S, 4, TGTL_Print},
	{"TGTR", (void**)&TGTR, NULL, (int**)&TGTR_NH, (int**)&TGTR_S, 4, TGTR_Print},
	{"TLV1", (void**)&TLV1, &TLV1_NS, &TLV1_NH, &TLV1_S, 8, TLV1_Print},
	{"TRCF", (void**)&TRCF, &TRCF_NS, &TRCF_NH, &TRCF_S, 4, TRCF_Print},
	{"TRGS", (void**)&TRGS, &TRGS_NS, &TRGS_NH, &TRGS_S, 384, TRGS_Print},
	{"TRKS", (void**)&TRKS, &TRKS_NS, &TRKS_NH, &TRKS_S, 40, TRKS_Print},
	{"TRL1", (void**)&TRL1, &TRL1_NS, &TRL1_NH, &TRL1_S, 40, TRL1_Print},
	{"TRPB", (void**)&TRPB, &TRPB_NS, &TRPB_NH, &TRPB_S, 32, TRPB_Print},
	{"TSPR", (void**)&TSPR, &TSPR_NS, &TSPR_NH, &TSPR_S, 4, TSPR_Print},
	{"TSRG", (void**)&TSRG, &TSRG_NS, &TSRG_NH, &TSRG_S, 4, TSRG_Print},
	{"UNUS", (void**)&UNUS, &UNUS_NS, &UNUS_NH, &UNUS_S, 32, UNUS_Print},
	{"VERT", (void**)&VERT, &VERT_NS, &VERT_NH, &VERT_S, 28, VERT_Print},
	{"RC26", (void**)&RC26, NULL, (int**)&RC26_NH, (int**)&RC26_S, 12, RC26_Print},
	{"STRE", (void**)&STRE, &STRE_NS, &STRE_NH, &STRE_S, 32, STRE_Print},
};
struct addressBanks_t *bankAddress=bankAdd;
