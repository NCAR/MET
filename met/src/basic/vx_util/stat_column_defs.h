// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

#ifndef  __STAT_COLUMN_DEFS_H__
#define  __STAT_COLUMN_DEFS_H__

////////////////////////////////////////////////////////////////////////

static const char * hdr_columns [] = {
   "VERSION",        "MODEL",
   "DESC",           "FCST_LEAD",
   "FCST_VALID_BEG", "FCST_VALID_END",
   "OBS_LEAD",
   "OBS_VALID_BEG",  "OBS_VALID_END",
   "FCST_VAR",       "FCST_UNITS",
   "FCST_LEV",       "OBS_VAR",
   "OBS_UNITS",      "OBS_LEV",
   "OBTYPE",         "VX_MASK",
   "INTERP_MTHD",    "INTERP_PNTS",
   "FCST_THRESH",    "OBS_THRESH",
   "COV_THRESH",     "ALPHA",
   "LINE_TYPE"
};

static const char * fho_columns [] = {
   "TOTAL",       "F_RATE",      "H_RATE",
   "O_RATE"
};

static const char * ctc_columns [] = {
   "TOTAL",       "FY_OY",       "FY_ON",
   "FN_OY",       "FN_ON"
};

static const char * cts_columns [] = {
   "TOTAL",
   "BASER",       "BASER_NCL",   "BASER_NCU",   "BASER_BCL",   "BASER_BCU",
   "FMEAN",       "FMEAN_NCL",   "FMEAN_NCU",   "FMEAN_BCL",   "FMEAN_BCU",
   "ACC",         "ACC_NCL",     "ACC_NCU",     "ACC_BCL",     "ACC_BCU",
   "FBIAS",       "FBIAS_BCL",   "FBIAS_BCU",
   "PODY",        "PODY_NCL",    "PODY_NCU",    "PODY_BCL",    "PODY_BCU",
   "PODN",        "PODN_NCL",    "PODN_NCU",    "PODN_BCL",    "PODN_BCU",
   "POFD",        "POFD_NCL",    "POFD_NCU",    "POFD_BCL",    "POFD_BCU",
   "FAR",         "FAR_NCL",     "FAR_NCU",     "FAR_BCL",     "FAR_BCU",
   "CSI",         "CSI_NCL",     "CSI_NCU",     "CSI_BCL",     "CSI_BCU",
   "GSS",         "GSS_BCL",     "GSS_BCU",
   "HK",          "HK_NCL",      "HK_NCU",      "HK_BCL",      "HK_BCU",
   "HSS",         "HSS_BCL",     "HSS_BCU",
   "ODDS",        "ODDS_NCL",    "ODDS_NCU",    "ODDS_BCL",    "ODDS_BCU",
   "LODDS",       "LODDS_NCL",   "LODDS_NCU",   "LODDS_BCL",   "LODDS_BCU",
   "ORSS",        "ORSS_NCL",    "ORSS_NCU",    "ORSS_BCL",    "ORSS_BCU",
   "EDS",         "EDS_NCL",     "EDS_NCU",     "EDS_BCL",     "EDS_BCU",
   "SEDS",        "SEDS_NCL",    "SEDS_NCU",    "SEDS_BCL",    "SEDS_BCU",
   "EDI",         "EDI_NCL",     "EDI_NCU",     "EDI_BCL",     "EDI_BCU",
   "SEDI",        "SEDI_NCL",    "SEDI_NCU",    "SEDI_BCL",    "SEDI_BCU",
   "BAGSS",       "BAGSS_BCL",   "BAGSS_BCU"
};

static const char * mctc_columns [] = {
   "TOTAL",       "N_CAT",       "Fi_Oj",       "EC_VALUE"
};

static const char * mcts_columns [] = {
   "TOTAL",       "N_CAT",
   "ACC",         "ACC_NCL",     "ACC_NCU",     "ACC_BCL",     "ACC_BCU",
   "HK",          "HK_BCL",      "HK_BCU",
   "HSS",         "HSS_BCL",     "HSS_BCU",
   "GER",         "GER_BCL",     "GER_BCU",
   "HSS_EC",      "HSS_EC_BCL",  "HSS_EC_BCU",  "EC_VALUE"
};

static const char * cnt_columns [] = {
   "TOTAL",
   "FBAR",             "FBAR_NCL",             "FBAR_NCU",          "FBAR_BCL",          "FBAR_BCU",
   "FSTDEV",           "FSTDEV_NCL",           "FSTDEV_NCU",        "FSTDEV_BCL",        "FSTDEV_BCU",
   "OBAR",             "OBAR_NCL",             "OBAR_NCU",          "OBAR_BCL",          "OBAR_BCU",
   "OSTDEV",           "OSTDEV_NCL",           "OSTDEV_NCU",        "OSTDEV_BCL",        "OSTDEV_BCU",
   "PR_CORR",          "PR_CORR_NCL",          "PR_CORR_NCU",       "PR_CORR_BCL",       "PR_CORR_BCU",
   "SP_CORR",          "KT_CORR",              "RANKS",             "FRANK_TIES",        "ORANK_TIES",
   "ME",               "ME_NCL",               "ME_NCU",            "ME_BCL",            "ME_BCU",
   "ESTDEV",           "ESTDEV_NCL",           "ESTDEV_NCU",        "ESTDEV_BCL",        "ESTDEV_BCU",
   "MBIAS",            "MBIAS_BCL",            "MBIAS_BCU",
   "MAE",              "MAE_BCL",              "MAE_BCU",
   "MSE",              "MSE_BCL",              "MSE_BCU",
   "BCMSE",            "BCMSE_BCL",            "BCMSE_BCU",
   "RMSE",             "RMSE_BCL",             "RMSE_BCU",
   "E10",              "E10_BCL",              "E10_BCU",
   "E25",              "E25_BCL",              "E25_BCU",
   "E50",              "E50_BCL",              "E50_BCU",
   "E75",              "E75_BCL",              "E75_BCU",
   "E90",              "E90_BCL",              "E90_BCU",
   "EIQR",             "EIQR_BCL",             "EIQR_BCU",
   "MAD",              "MAD_BCL",              "MAD_BCU",
   "ANOM_CORR",        "ANOM_CORR_NCL",        "ANOM_CORR_NCU",     "ANOM_CORR_BCL",     "ANOM_CORR_BCU",
   "ME2",              "ME2_BCL",              "ME2_BCU",
   "MSESS",            "MSESS_BCL",            "MSESS_BCU",
   "RMSFA",            "RMSFA_BCL",            "RMSFA_BCU",
   "RMSOA",            "RMSOA_BCL",            "RMSOA_BCU",
   "ANOM_CORR_UNCNTR", "ANOM_CORR_UNCNTR_BCL", "ANOM_CORR_UNCNTR_BCU"
};

static const char * sl1l2_columns [] = {
   "TOTAL",       "FBAR",        "OBAR",
   "FOBAR",       "FFBAR",       "OOBAR",
   "MAE"
};

static const char * sal1l2_columns [] = {
   "TOTAL",       "FABAR",       "OABAR",
   "FOABAR",      "FFABAR",      "OOABAR",
   "MAE"
};

static const char * vl1l2_columns [] = {
   "TOTAL",       "UFBAR",       "VFBAR",
   "UOBAR",       "VOBAR",       "UVFOBAR",
   "UVFFBAR",     "UVOOBAR",     "F_SPEED_BAR",
   "O_SPEED_BAR",
};

static const char * val1l2_columns [] = {
   "TOTAL",       "UFABAR",      "VFABAR",
   "UOABAR",      "VOABAR",      "UVFOABAR",
   "UVFFABAR",    "UVOOABAR"
};


static const char * vcnt_columns [] = {

   "TOTAL",

   "FBAR",         "FBAR_BCL",         "FBAR_BCU",
   "OBAR",         "OBAR_BCL",         "OBAR_BCU",
   "FS_RMS",       "FS_RMS_BCL",       "FS_RMS_BCU",
   "OS_RMS",       "OS_RMS_BCL",       "OS_RMS_BCU",
   "MSVE",         "MSVE_BCL",         "MSVE_BCU",
   "RMSVE",        "RMSVE_BCL",        "RMSVE_BCU",
   "FSTDEV",       "FSTDEV_BCL",       "FSTDEV_BCU",
   "OSTDEV",       "OSTDEV_BCL",       "OSTDEV_BCU",
   "FDIR",         "FDIR_BCL",         "FDIR_BCU",
   "ODIR",         "ODIR_BCL",         "ODIR_BCU",
   "FBAR_SPEED",   "FBAR_SPEED_BCL",   "FBAR_SPEED_BCU",
   "OBAR_SPEED",   "OBAR_SPEED_BCL",   "OBAR_SPEED_BCU",
   "VDIFF_SPEED",  "VDIFF_SPEED_BCL",  "VDIFF_SPEED_BCU",
   "VDIFF_DIR",    "VDIFF_DIR_BCL",    "VDIFF_DIR_BCU",
   "SPEED_ERR",    "SPEED_ERR_BCL",    "SPEED_ERR_BCU",
   "SPEED_ABSERR", "SPEED_ABSERR_BCL", "SPEED_ABSERR_BCU",
   "DIR_ERR",      "DIR_ERR_BCL",      "DIR_ERR_BCU",
   "DIR_ABSERR",   "DIR_ABSERR_BCL",   "DIR_ABSERR_BCU",

};

static const char * pct_columns [] = {
   "TOTAL",       "N_THRESH",    "THRESH_",
   "OY_",         "ON_"
};

static const char * pstd_columns [] = {
   "TOTAL",       "N_THRESH",    "BASER",
   "BASER_NCL",   "BASER_NCU",   "RELIABILITY",
   "RESOLUTION",  "UNCERTAINTY", "ROC_AUC",
   "BRIER",       "BRIER_NCL",   "BRIER_NCU",
   "BRIERCL",     "BRIERCL_NCL", "BRIERCL_NCU",
   "BSS",         "BSS_SMPL",    "THRESH_",
};

static const char * pjc_columns [] = {
   "TOTAL",       "N_THRESH",    "THRESH_",
   "OY_TP_",      "ON_TP_",      "CALIBRATION_",
   "REFINEMENT_", "LIKELIHOOD_", "BASER_"
};

static const char * prc_columns [] = {
   "TOTAL",       "N_THRESH",    "THRESH_",
   "PODY_",       "POFD_"
};

static const char * eclv_columns [] = {
   "TOTAL",       "BASER",    "VALUE_BASER",
   "N_PNT",       "CL_",      "VALUE_"
};

static const char * mpr_columns [] = {
   "TOTAL",       "INDEX",       "OBS_SID",
   "OBS_LAT",     "OBS_LON",     "OBS_LVL",
   "OBS_ELV",     "FCST",        "OBS",
   "OBS_QC",      "CLIMO_MEAN",  "CLIMO_STDEV",
   "CLIMO_CDF"
};

static const char * nbrctc_columns [] = {
   "TOTAL",       "FY_OY",       "FY_ON",
   "FN_OY",       "FN_ON"
};

static const char * nbrcts_columns [] = {
   "TOTAL",
   "BASER",       "BASER_NCL",   "BASER_NCU",   "BASER_BCL",   "BASER_BCU",
   "FMEAN",       "FMEAN_NCL",   "FMEAN_NCU",   "FMEAN_BCL",   "FMEAN_BCU",
   "ACC",         "ACC_NCL",     "ACC_NCU",     "ACC_BCL",     "ACC_BCU",
   "FBIAS",       "FBIAS_BCL",   "FBIAS_BCU",
   "PODY",        "PODY_NCL",    "PODY_NCU",    "PODY_BCL",    "PODY_BCU",
   "PODN",        "PODN_NCL",    "PODN_NCU",    "PODN_BCL",    "PODN_BCU",
   "POFD",        "POFD_NCL",    "POFD_NCU",    "POFD_BCL",    "POFD_BCU",
   "FAR",         "FAR_NCL",     "FAR_NCU",     "FAR_BCL",     "FAR_BCU",
   "CSI",         "CSI_NCL",     "CSI_NCU",     "CSI_BCL",     "CSI_BCU",
   "GSS",         "GSS_BCL",     "GSS_BCU",
   "HK",          "HK_NCL",      "HK_NCU",      "HK_BCL",      "HK_BCU",
   "HSS",         "HSS_BCL",     "HSS_BCU",
   "ODDS",        "ODDS_NCL",    "ODDS_NCU",    "ODDS_BCL",    "ODDS_BCU",
   "LODDS",       "LODDS_NCL",   "LODDS_NCU",   "LODDS_BCL",   "LODDS_BCU",
   "ORSS",        "ORSS_NCL",    "ORSS_NCU",    "ORSS_BCL",    "ORSS_BCU",
   "EDS",         "EDS_NCL",     "EDS_NCU",     "EDS_BCL",     "EDS_BCU",
   "SEDS",        "SEDS_NCL",    "SEDS_NCU",    "SEDS_BCL",    "SEDS_BCU",
   "EDI",         "EDI_NCL",     "EDI_NCU",     "EDI_BCL",     "EDI_BCU",
   "SEDI",        "SEDI_NCL",    "SEDI_NCU",    "SEDI_BCL",    "SEDI_BCU",
   "BAGSS",       "BAGSS_BCL",   "BAGSS_BCU"
};

static const char * nbrcnt_columns [] = {
   "TOTAL",
   "FBS",         "FBS_BCL",     "FBS_BCU",
   "FSS",         "FSS_BCL",     "FSS_BCU",
   "AFSS",        "AFSS_BCL",    "AFSS_BCU",
   "UFSS",        "UFSS_BCL",    "UFSS_BCU",
   "F_RATE",      "F_RATE_BCL",  "F_RATE_BCU",
   "O_RATE",      "O_RATE_BCL",  "O_RATE_BCU"
};

static const char * grad_columns [] = {
   "TOTAL",
   "FGBAR",       "OGBAR",       "MGBAR",
   "EGBAR",       "S1",          "S1_OG",
   "FGOG_RATIO",  "DX",          "DY"
};

static const char * dmap_columns [] = {
   "TOTAL",       "FY",          "OY",
   "FBIAS",       "BADDELEY",    "HAUSDORFF",
   "MED_FO",      "MED_OF",      "MED_MIN",     "MED_MAX",     "MED_MEAN",
   "FOM_FO",      "FOM_OF",      "FOM_MIN",     "FOM_MAX",     "FOM_MEAN",
   "ZHU_FO",      "ZHU_OF",      "ZHU_MIN",     "ZHU_MAX",     "ZHU_MEAN",
   "G",           "GBETA",       "BETA_VALUE"
};

static const char * isc_columns [] = {
   "TOTAL",
   "TILE_DIM",    "TILE_XLL",    "TILE_YLL",
   "NSCALE",      "ISCALE",      "MSE",
   "ISC",         "FENERGY2",    "OENERGY2",
   "BASER",       "FBIAS"
};

static const char * ecnt_columns [] = {
   "TOTAL",       "N_ENS",       "CRPS",
   "CRPSS",       "IGN",         "ME",
   "RMSE",        "SPREAD",      "ME_OERR",
   "RMSE_OERR",   "SPREAD_OERR", "SPREAD_PLUS_OERR",
   "CRPSCL",      "CRPS_EMP",    "CRPSCL_EMP",
   "CRPSS_EMP"
};

static const char * rps_columns [] = {
   "TOTAL",       "N_PROB",      "RPS_REL",
   "RPS_RES",     "RPS_UNC",     "RPS",
   "RPSS",        "RPSS_SMPL",   "RPS_COMP"
};

static const char * rhist_columns [] = {
   "TOTAL",       "N_RANK",      "RANK_"
};

static const char * phist_columns [] = {
   "TOTAL",       "BIN_SIZE",    "N_BIN",
   "BIN_"
};

static const char * orank_columns [] = {
   "TOTAL",       "INDEX",            "OBS_SID",
   "OBS_LAT",     "OBS_LON",          "OBS_LVL",
   "OBS_ELV",     "OBS",              "PIT",
   "RANK",        "N_ENS_VLD",        "N_ENS",
   "ENS_",        "OBS_QC",           "ENS_MEAN",
   "CLIMO_MEAN",  "SPREAD",           "ENS_MEAN_OERR",
   "SPREAD_OERR", "SPREAD_PLUS_OERR", "CLIMO_STDEV"
};

static const char * ssvar_columns [] = {
   "TOTAL",       "N_BIN",       "BIN_i",
   "BIN_N",       "VAR_MIN",     "VAR_MAX",
   "VAR_MEAN",    "FBAR",        "OBAR",
   "FOBAR",       "FFBAR",       "OOBAR",
   "FBAR_NCL",    "FBAR_NCU",
   "FSTDEV",      "FSTDEV_NCL",  "FSTDEV_NCU",
   "OBAR_NCL",    "OBAR_NCU",
   "OSTDEV",      "OSTDEV_NCL",  "OSTDEV_NCU",
   "PR_CORR",     "PR_CORR_NCL", "PR_CORR_NCU",
   "ME",          "ME_NCL",      "ME_NCU",
   "ESTDEV",      "ESTDEV_NCL",  "ESTDEV_NCU",
   "MBIAS",       "MSE",         "BCMSE",
   "RMSE"
};

static const char * genmpr_columns [] = {
   "TOTAL",       "INDEX",       "STORM_ID",
   "AGEN_INIT",   "AGEN_FHR",
   "AGEN_LAT",    "AGEN_LON",    "AGEN_DLAND",
   "BGEN_LAT",    "BGEN_LON",    "BGEN_DLAND",
   "GEN_DIST",    "GEN_TDIFF",   "INIT_TDIFF",
   "DEV_CAT",     "OPS_CAT"
};

static const char * relp_columns [] = {
   "TOTAL",       "N_ENS",       "RELP_"
};

static const char * job_summary_columns [] = {
   "TOTAL",
   "MEAN",        "MEAN_NCL",    "MEAN_NCU",     "MEAN_BCL",    "MEAN_BCU",
   "STDEV",       "STDEV_BCL",   "STDEV_BCU",
   "MIN",
   "P10",         "P25",         "P50",          "P75",         "P90",
   "MAX",         "IQR",         "RANGE",
   "WMO_TYPE",    "WMO_MEAN",    "WMO_WEIGHTED_MEAN"
};

static const char * job_go_columns [] = {
   "GO_INDEX"
};

static const char * job_ss_columns [] = {
   "SS_INDEX"
};

static const char * job_wdir_columns [] = {
   "TOTAL",
   "FBAR",        "OBAR",        "ME",           "MAE"
};

static const char * job_ramp_columns [] = {
   "TYPE",
   "FCOLUMN",    "OCOLUMN",
   "FTIME",      "OTIME",
   "FEXACT",     "OEXACT",
   "FTHRESH",    "OTHRESH",
   "WINDOW_BEG", "WINDOW_END"
};

static const char * job_ramp_mpr_columns [] = {
   "TOTAL", "INDEX",
   "INIT",  "LEAD", "VALID",
   "FPRV",  "FCUR", "FDLT", "FRAMP",
   "OPRV",  "OCUR", "ODLT", "ORAMP",
   "CATEGORY"
};

////////////////////////////////////////////////////////////////////////

static const int max_stat_col           = 100;

static const int n_header_columns       = sizeof(hdr_columns)/sizeof(*hdr_columns);

static const int n_fho_columns          = sizeof(fho_columns)/sizeof(*fho_columns);
static const int n_ctc_columns          = sizeof(ctc_columns)/sizeof(*ctc_columns);
static const int n_cts_columns          = sizeof(cts_columns)/sizeof(*cts_columns);
static const int n_mctc_columns         = sizeof(mctc_columns)/sizeof(*mctc_columns);
static const int n_mcts_columns         = sizeof(mcts_columns)/sizeof(*mcts_columns);

static const int n_cnt_columns          = sizeof(cnt_columns)/sizeof(*cnt_columns);
static const int n_sl1l2_columns        = sizeof(sl1l2_columns)/sizeof(*sl1l2_columns);
static const int n_sal1l2_columns       = sizeof(sal1l2_columns)/sizeof(*sal1l2_columns);

static const int n_vl1l2_columns        = sizeof(vl1l2_columns)/sizeof(*vl1l2_columns);
static const int n_val1l2_columns       = sizeof(val1l2_columns)/sizeof(*val1l2_columns);
static const int n_vcnt_columns         = sizeof(vcnt_columns)/sizeof(*vcnt_columns);

static const int n_pct_columns          = sizeof(pct_columns)/sizeof(*pct_columns);
static const int n_pstd_columns         = sizeof(pstd_columns)/sizeof(*pstd_columns);
static const int n_pjc_columns          = sizeof(pjc_columns)/sizeof(*pjc_columns);
static const int n_prc_columns          = sizeof(prc_columns)/sizeof(*prc_columns);
static const int n_eclv_columns         = sizeof(eclv_columns)/sizeof(*eclv_columns);

static const int n_mpr_columns          = sizeof(mpr_columns)/sizeof(*mpr_columns);

static const int n_nbrctc_columns       = sizeof(nbrctc_columns)/sizeof(*nbrctc_columns);
static const int n_nbrcts_columns       = sizeof(nbrcts_columns)/sizeof(*nbrcts_columns);
static const int n_nbrcnt_columns       = sizeof(nbrcnt_columns)/sizeof(*nbrcnt_columns);

static const int n_grad_columns         = sizeof(grad_columns)/sizeof(*grad_columns);
static const int n_dmap_columns         = sizeof(dmap_columns)/sizeof(*dmap_columns);
static const int n_isc_columns          = sizeof(isc_columns)/sizeof(*isc_columns);

static const int n_job_summary_columns  = sizeof(job_summary_columns)/sizeof(*job_summary_columns);
static const int n_job_go_columns       = sizeof(job_go_columns)/sizeof(*job_go_columns);
static const int n_job_ss_columns       = sizeof(job_ss_columns)/sizeof(*job_ss_columns);
static const int n_job_wdir_columns     = sizeof(job_wdir_columns)/sizeof(*job_wdir_columns);
static const int n_job_ramp_columns     = sizeof(job_ramp_columns)/sizeof(*job_ramp_columns);
static const int n_job_ramp_mpr_columns = sizeof(job_ramp_mpr_columns)/sizeof(*job_ramp_mpr_columns);

static const int n_ecnt_columns         = sizeof(ecnt_columns)/sizeof(*ecnt_columns);
static const int n_rps_columns          = sizeof(rps_columns)/sizeof(*rps_columns);

static const int n_rhist_columns        = sizeof(rhist_columns)/sizeof(*rhist_columns);
static const int n_phist_columns        = sizeof(phist_columns)/sizeof(*phist_columns);
static const int n_orank_columns        = sizeof(orank_columns)/sizeof(*orank_columns);
static const int n_ssvar_columns        = sizeof(ssvar_columns)/sizeof(*ssvar_columns);
static const int n_relp_columns         = sizeof(relp_columns)/sizeof(*relp_columns);

static const int n_genmpr_columns       = sizeof(genmpr_columns)/sizeof(*genmpr_columns);

////////////////////////////////////////////////////////////////////////

inline int get_n_mctc_columns  (int n) { return(3  + n*n);             } // n = N_CAT
inline int get_n_pct_columns   (int n) { return(3  + 3*(max(1, n)-1)); } // n = N_THRESH
inline int get_n_pstd_columns  (int n) { return(17 +    max(1, n)   ); } // n = N_THRESH
inline int get_n_pjc_columns   (int n) { return(3  + 7*(max(1, n)-1)); } // n = N_THRESH
inline int get_n_prc_columns   (int n) { return(3  + 3*(max(1, n)-1)); } // n = N_THRESH
inline int get_n_eclv_columns  (int n) { return(4  + 2*n);             } // n = N_PNT
inline int get_n_rhist_columns (int n) { return(2  + n);               } // n = N_RANK
inline int get_n_phist_columns (int n) { return(3  + n);               } // n = N_BINS
inline int get_n_relp_columns  (int n) { return(2  + n);               } // n = N_ENS
inline int get_n_orank_columns (int n) { return(20 + n);               } // n = N_ENS

////////////////////////////////////////////////////////////////////////

#endif   /*  __STAT_COLUMN_DEFS_H__  */

////////////////////////////////////////////////////////////////////////
