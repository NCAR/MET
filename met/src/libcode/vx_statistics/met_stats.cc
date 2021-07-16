// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

using namespace std;

#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include "met_stats.h"
#include "pair_data_point.h"
#include "compute_ci.h"
#include "grib_strings.h"
#include "vx_util.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////
//
// Code for class CIInfo
//
////////////////////////////////////////////////////////////////////////

CIInfo::CIInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

CIInfo::~CIInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

CIInfo::CIInfo(const CIInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

CIInfo & CIInfo::operator=(const CIInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void CIInfo::init_from_scratch() {

   v_ncl = (double *) 0;
   v_ncu = (double *) 0;

   v_bcl = (double *) 0;
   v_bcu = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void CIInfo::clear() {

   n   = 0;
   v   = bad_data_double;
   vif = 1.0;

   if(v_ncl) { delete [] v_ncl; v_ncl = (double *) 0; }
   if(v_ncu) { delete [] v_ncu; v_ncu = (double *) 0; }

   if(v_bcl) { delete [] v_bcl; v_bcl = (double *) 0; }
   if(v_bcu) { delete [] v_bcu; v_bcu = (double *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void CIInfo::set_bad_data() {
   int i;

   v   = bad_data_double;
   vif = 1.0;

   for(i=0; i<n; i++) {
      v_ncl[i] = v_ncu[i] = bad_data_double;
      v_bcl[i] = v_bcu[i] = bad_data_double;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void CIInfo::assign(const CIInfo &c) {
   int i;

   clear();

   allocate_n_alpha(c.n);

   v   = c.v;
   vif = c.vif;

   for(i=0; i<c.n; i++) {
      v_ncl[i] = c.v_ncl[i];
      v_ncu[i] = c.v_ncu[i];
      v_bcl[i] = c.v_bcl[i];
      v_bcu[i] = c.v_bcu[i];
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void CIInfo::allocate_n_alpha(int i) {
   int j;

   n = i;

   if(n > 0) {
      v_ncl = new double [n];
      v_ncu = new double [n];
      v_bcl = new double [n];
      v_bcu = new double [n];

      if(!v_ncl || !v_ncu || !v_bcl || !v_bcu) {
         mlog << Error << "\nCIInfo::allocate_n_alpha() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }
   }

   // Initialize the values
   for(j=0; j<n; j++) {
      v_ncl[j] = bad_data_double;
      v_ncu[j] = bad_data_double;
      v_bcl[j] = bad_data_double;
      v_bcu[j] = bad_data_double;
   }

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class CTSInfo
//
////////////////////////////////////////////////////////////////////////

CTSInfo::CTSInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

CTSInfo::~CTSInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

CTSInfo::CTSInfo(const CTSInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

CTSInfo & CTSInfo::operator=(const CTSInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::init_from_scratch() {

   alpha = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::clear() {

   n_alpha = 0;
   if(alpha) { delete [] alpha; alpha = (double *) 0; }

   cts.zero_out();
   fthresh.clear();
   othresh.clear();

   baser.clear();
   fmean.clear();
   acc.clear();
   fbias.clear();
   pody.clear();
   podn.clear();
   pofd.clear();
   far.clear();
   csi.clear();
   gss.clear();
   bagss.clear();
   hk.clear();
   hss.clear();
   odds.clear();
   lodds.clear();
   orss.clear();
   eds.clear();
   seds.clear();
   edi.clear();
   sedi.clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::assign(const CTSInfo &c) {
   int i;

   clear();

   cts = c.cts;
   fthresh = c.fthresh;
   othresh = c.othresh;

   allocate_n_alpha(c.n_alpha);
   for(i=0; i<c.n_alpha; i++) { alpha[i] = c.alpha[i]; }

   baser = c.baser;
   fmean = c.fmean;
   acc = c.acc;
   fbias = c.fbias;
   pody = c.pody;
   podn = c.podn;
   pofd = c.pofd;
   far = c.far;
   csi = c.csi;
   gss = c.gss;
   bagss = c.bagss;
   hk = c.hk;
   hss = c.hss;
   odds = c.odds;
   lodds = c.lodds;
   orss = c.orss;
   eds = c.eds;
   seds = c.seds;
   edi = c.edi;
   sedi = c.sedi;

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::allocate_n_alpha(int i) {

   n_alpha = i;

   if(n_alpha > 0) {

      alpha = new double [n_alpha];

      if(!alpha) {
         mlog << Error << "\nCTSInfo::allocate_n() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }

      baser.allocate_n_alpha(n_alpha);
      fmean.allocate_n_alpha(n_alpha);
      acc.allocate_n_alpha(n_alpha);
      fbias.allocate_n_alpha(n_alpha);
      pody.allocate_n_alpha(n_alpha);
      podn.allocate_n_alpha(n_alpha);
      pofd.allocate_n_alpha(n_alpha);
      far.allocate_n_alpha(n_alpha);
      csi.allocate_n_alpha(n_alpha);
      gss.allocate_n_alpha(n_alpha);
      bagss.allocate_n_alpha(n_alpha);
      hk.allocate_n_alpha(n_alpha);
      hss.allocate_n_alpha(n_alpha);
      odds.allocate_n_alpha(n_alpha);
      lodds.allocate_n_alpha(n_alpha);
      orss.allocate_n_alpha(n_alpha);
      eds.allocate_n_alpha(n_alpha);
      seds.allocate_n_alpha(n_alpha);
      edi.allocate_n_alpha(n_alpha);
      sedi.allocate_n_alpha(n_alpha);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::add(double f, double o) {
   add(f, o, bad_data_double, bad_data_double);
   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::add(double f, double o, double cmn, double csd) {

   if     ( fthresh.check(f, cmn, csd) &&  othresh.check(o, cmn, csd)) cts.inc_fy_oy();
   else if( fthresh.check(f, cmn, csd) && !othresh.check(o, cmn, csd)) cts.inc_fy_on();
   else if(!fthresh.check(f, cmn, csd) &&  othresh.check(o, cmn, csd)) cts.inc_fn_oy();
   else if(!fthresh.check(f, cmn, csd) && !othresh.check(o, cmn, csd)) cts.inc_fn_on();

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::compute_stats() {

   baser.v = cts.oy_tp();
   fmean.v = cts.fy_tp();
   acc.v   = cts.accuracy();
   fbias.v = cts.fbias();
   pody.v  = cts.pod_yes();
   podn.v  = cts.pod_no();
   pofd.v  = cts.pofd();
   far.v   = cts.far();
   csi.v   = cts.csi();
   gss.v   = cts.gss();
   bagss.v = cts.bagss();
   hk.v    = cts.hk();
   hss.v   = cts.hss();
   odds.v  = cts.odds();
   lodds.v = cts.lodds();
   orss.v  = cts.orss();
   eds.v   = cts.eds();
   seds.v  = cts.seds();
   edi.v   = cts.edi();
   sedi.v  = cts.sedi();

   return;
}

////////////////////////////////////////////////////////////////////////

void CTSInfo::compute_ci() {
   int i;

   //
   // Compute confidence intervals for each alpha value specified
   //
   for(i=0; i<n_alpha; i++) {

      //
      // Compute confidence intervals for the scores based on
      // proportions
      //
      compute_proportion_ci(baser.v, cts.n(), alpha[i], baser.vif,
                            baser.v_ncl[i], baser.v_ncu[i]);
      compute_proportion_ci(fmean.v, cts.n(), alpha[i], fmean.vif,
                            fmean.v_ncl[i], fmean.v_ncu[i]);
      compute_proportion_ci(acc.v, cts.n(), alpha[i], acc.vif,
                            acc.v_ncl[i], acc.v_ncu[i]);
      compute_proportion_ci(pody.v, cts.n(), alpha[i], pody.vif,
                            pody.v_ncl[i], pody.v_ncu[i]);
      compute_proportion_ci(podn.v, cts.n(), alpha[i], podn.vif,
                            podn.v_ncl[i], podn.v_ncu[i]);
      compute_proportion_ci(pofd.v, cts.n(), alpha[i], pofd.vif,
                            pofd.v_ncl[i], pofd.v_ncu[i]);
      compute_proportion_ci(far.v, cts.n(), alpha[i], far.vif,
                            far.v_ncl[i], far.v_ncu[i]);
      compute_proportion_ci(csi.v, cts.n(), alpha[i], csi.vif,
                            csi.v_ncl[i], csi.v_ncu[i]);

      //
      // Compute a confidence interval for Hanssen and Kuipers discriminant
      //
      compute_hk_ci(hk.v, alpha[i], hk.vif,
                    cts.fy_oy(), cts.fy_on(), cts.fn_oy(), cts.fn_on(),
                    hk.v_ncl[i], hk.v_ncu[i]);

      //
      // Compute a confidence interval for the odds ratio
      //
      compute_woolf_ci(odds.v, alpha[i],
                       cts.fy_oy(), cts.fy_on(), cts.fn_oy(), cts.fn_on(),
                       odds.v_ncl[i], odds.v_ncu[i]);

      cts.lodds_ci(alpha[i], lodds.v_ncl[i], lodds.v_ncu[i]);
      cts.orss_ci (alpha[i], orss.v_ncl[i],  orss.v_ncu[i]);
      cts.eds_ci  (alpha[i], eds.v_ncl[i],   eds.v_ncu[i]);
      cts.seds_ci (alpha[i], seds.v_ncl[i],  seds.v_ncu[i]);
      cts.edi_ci  (alpha[i], edi.v_ncl[i],   edi.v_ncu[i]);
      cts.sedi_ci (alpha[i], sedi.v_ncl[i],  sedi.v_ncu[i]);

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

double CTSInfo::get_stat(const char *stat_name) {
   double v = bad_data_double;

        if(strcmp(stat_name, "TOTAL") == 0) v = cts.n();
   else if(strcmp(stat_name, "BASER") == 0) v = cts.baser();
   else if(strcmp(stat_name, "FMEAN") == 0) v = cts.fmean();
   else if(strcmp(stat_name, "ACC"  ) == 0) v = cts.accuracy();
   else if(strcmp(stat_name, "FBIAS") == 0) v = cts.fbias();
   else if(strcmp(stat_name, "PODY" ) == 0) v = cts.pod_yes();
   else if(strcmp(stat_name, "PODN" ) == 0) v = cts.pod_no();
   else if(strcmp(stat_name, "POFD" ) == 0) v = cts.pofd();
   else if(strcmp(stat_name, "FAR"  ) == 0) v = cts.far();
   else if(strcmp(stat_name, "CSI"  ) == 0) v = cts.csi();
   else if(strcmp(stat_name, "GSS"  ) == 0) v = cts.gss();
   else if(strcmp(stat_name, "HK"   ) == 0) v = cts.hk();
   else if(strcmp(stat_name, "HSS"  ) == 0) v = cts.hss();
   else if(strcmp(stat_name, "ODDS" ) == 0) v = cts.odds();
   else if(strcmp(stat_name, "LODDS") == 0) v = cts.lodds();
   else if(strcmp(stat_name, "ORSS" ) == 0) v = cts.orss();
   else if(strcmp(stat_name, "EDS"  ) == 0) v = cts.eds();
   else if(strcmp(stat_name, "SEDS" ) == 0) v = cts.seds();
   else if(strcmp(stat_name, "EDI"  ) == 0) v = cts.edi();
   else if(strcmp(stat_name, "SEDI" ) == 0) v = cts.sedi();
   else if(strcmp(stat_name, "BAGSS") == 0) v = cts.bagss();
   else {
      mlog << Error << "\nCTSInfo::get_stat() -> "
           << "unknown categorical statistic name \"" << stat_name
           << "\"!\n\n";
      exit(1);
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////
//
// Code for class MCTSInfo
//
////////////////////////////////////////////////////////////////////////

MCTSInfo::MCTSInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

MCTSInfo::~MCTSInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

MCTSInfo::MCTSInfo(const MCTSInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

MCTSInfo & MCTSInfo::operator=(const MCTSInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::init_from_scratch() {

   alpha = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::clear() {

   n_alpha = 0;
   if(alpha) { delete [] alpha; alpha = (double *) 0; }

   cts.zero_out();
   fthresh.clear();
   othresh.clear();

   acc.clear();
   hk.clear();
   hss.clear();
   hss_ec.clear();
   ger.clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::assign(const MCTSInfo &c) {
   int i;

   clear();

   cts = c.cts;
   fthresh = c.fthresh;
   othresh = c.othresh;

   allocate_n_alpha(c.n_alpha);
   for(i=0; i<c.n_alpha; i++) { alpha[i] = c.alpha[i]; }

   acc = c.acc;
   hk = c.hk;
   hss = c.hss;
   hss_ec = c.hss_ec;
   ger = c.ger;

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::allocate_n_alpha(int i) {

   n_alpha = i;

   if(n_alpha > 0) {

      alpha = new double [n_alpha];

      if(!alpha) {
         mlog << Error << "\nMCTSInfo::allocate_n() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }

      acc.allocate_n_alpha(n_alpha);
      hk.allocate_n_alpha(n_alpha);
      hss.allocate_n_alpha(n_alpha);
      hss_ec.allocate_n_alpha(n_alpha);
      ger.allocate_n_alpha(n_alpha);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::set_fthresh(const ThreshArray &ta) {

   // Check that the thresholds are formatted correctly
   ta.check_bin_thresh();
   fthresh = ta;

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::set_othresh(const ThreshArray &ta) {

   // Check that the thresholds are formatted correctly
   ta.check_bin_thresh();
   othresh = ta;

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::add(double f, double o) {
   add(f, o, bad_data_double, bad_data_double);
   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::add(double f, double o, double cmn, double csd) {
   int r, c;

   // Find the row and column for the forecast and observation values.
   r = fthresh.check_bins(f, cmn, csd);
   c = othresh.check_bins(o, cmn, csd);

   // Increment the corresponding contingency table entry.
   cts.inc_entry(r, c);

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::compute_stats() {

   acc.v    = cts.gaccuracy();
   hk.v     = cts.gkuiper();
   hss.v    = cts.gheidke();
   hss_ec.v = cts.gheidke_ec(cts.ec_value());
   ger.v    = cts.gerrity();

   return;
}

////////////////////////////////////////////////////////////////////////

void MCTSInfo::compute_ci() {
   int i;

   //
   // Compute confidence intervals for each alpha value specified
   //
   for(i=0; i<n_alpha; i++) {

      //
      // Compute confidence intervals for the scores based on
      // proportions
      //
      compute_proportion_ci(acc.v, cts.total(), alpha[i], acc.vif,
                            acc.v_ncl[i], acc.v_ncu[i]);
   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class CNTInfo
//
////////////////////////////////////////////////////////////////////////

CNTInfo::CNTInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

CNTInfo::~CNTInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

CNTInfo::CNTInfo(const CNTInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

CNTInfo & CNTInfo::operator=(const CNTInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void CNTInfo::init_from_scratch() {

   alpha = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void CNTInfo::clear() {

   n = 0;
   n_alpha = 0;
   if(alpha) { delete [] alpha; alpha = (double *) 0; }

   fthresh.clear();
   othresh.clear();
   logic = SetLogic_None;

   fbar.clear();
   fstdev.clear();
   obar.clear();
   ostdev.clear();
   pr_corr.clear();
   sp_corr.clear();
   kt_corr.clear();
   anom_corr.clear();
   rmsfa.clear();
   rmsoa.clear();
   anom_corr_uncntr.clear();
   me.clear();
   me2.clear();
   estdev.clear();
   mbias.clear();
   mae.clear();
   mse.clear();
   msess.clear();
   bcmse.clear();
   rmse.clear();
   e10.clear();
   e25.clear();
   e50.clear();
   e75.clear();
   e90.clear();
   eiqr.clear();
   mad.clear();

   n_ranks = frank_ties = orank_ties = 0;

   return;
}

////////////////////////////////////////////////////////////////////////

void CNTInfo::assign(const CNTInfo &c) {
   int i;

   clear();

   fthresh = c.fthresh;
   othresh = c.othresh;
   logic   = c.logic;

   n = c.n;
   allocate_n_alpha(c.n_alpha);
   for(i=0; i<c.n_alpha; i++) { alpha[i] = c.alpha[i]; }

   fbar             = c.fbar;
   fstdev           = c.fstdev;
   obar             = c.obar;
   ostdev           = c.ostdev;
   pr_corr          = c.pr_corr;
   sp_corr          = c.sp_corr;
   kt_corr          = c.kt_corr;
   anom_corr        = c.anom_corr;
   rmsfa            = c.rmsfa;
   rmsoa            = c.rmsoa;
   anom_corr_uncntr = c.anom_corr_uncntr;
   me               = c.me;
   me2              = c.me2;
   estdev           = c.estdev;
   mbias            = c.mbias;
   mae              = c.mae;
   mse              = c.mse;
   msess            = c.msess;
   bcmse            = c.bcmse;
   rmse             = c.rmse;
   e10              = c.e10;
   e25              = c.e25;
   e50              = c.e50;
   e75              = c.e75;
   e90              = c.e90;
   eiqr             = c.eiqr;
   mad              = c.mad;
   n_ranks          = c.n_ranks;
   frank_ties       = c.frank_ties;
   orank_ties       = c.orank_ties;

   return;
}

////////////////////////////////////////////////////////////////////////

void CNTInfo::allocate_n_alpha(int i) {

   n_alpha = i;

   if(n_alpha > 0) {

      alpha = new double [n_alpha];

      if(!alpha) {
         mlog << Error << "\nCNTInfo::allocate_n_alpha() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }

      fbar.allocate_n_alpha(n_alpha);
      fstdev.allocate_n_alpha(n_alpha);
      obar.allocate_n_alpha(n_alpha);
      ostdev.allocate_n_alpha(n_alpha);
      pr_corr.allocate_n_alpha(n_alpha);
      sp_corr.allocate_n_alpha(n_alpha);
      kt_corr.allocate_n_alpha(n_alpha);
      anom_corr.allocate_n_alpha(n_alpha);
      rmsfa.allocate_n_alpha(n_alpha);
      rmsoa.allocate_n_alpha(n_alpha);
      anom_corr_uncntr.allocate_n_alpha(n_alpha);
      me.allocate_n_alpha(n_alpha);
      me2.allocate_n_alpha(n_alpha);
      estdev.allocate_n_alpha(n_alpha);
      mbias.allocate_n_alpha(n_alpha);
      mae.allocate_n_alpha(n_alpha);
      mse.allocate_n_alpha(n_alpha);
      msess.allocate_n_alpha(n_alpha);
      bcmse.allocate_n_alpha(n_alpha);
      rmse.allocate_n_alpha(n_alpha);
      e10.allocate_n_alpha(n_alpha);
      e25.allocate_n_alpha(n_alpha);
      e50.allocate_n_alpha(n_alpha);
      e75.allocate_n_alpha(n_alpha);
      e90.allocate_n_alpha(n_alpha);
      eiqr.allocate_n_alpha(n_alpha);
      mad.allocate_n_alpha(n_alpha);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void CNTInfo::compute_ci() {
   int i;
   double cv_normal_l, cv_normal_u;
   double cv_chi2_l, cv_chi2_u;
   double v, cl, cu;

   //
   // Compute the confidence interval for each alpha value specified
   // In computing the confidence intervals, the spatial correlation between
   // adjacent points is ignored, and certain assumptions of normality are
   // made.  The user must decide if the computation method is appropriate
   // for the chosen field.
   //
   for(i=0; i<n_alpha; i++) {

      //
      // Check for the degenerate case
      //
      if(n <= 1) {
         fbar.v_ncl[i]          = fbar.v_ncu[i]          = bad_data_double;
         fstdev.v_ncl[i]        = fstdev.v_ncu[i]        = bad_data_double;
         obar.v_ncl[i]          = obar.v_ncu[i]          = bad_data_double;
         ostdev.v_ncl[i]        = ostdev.v_ncu[i]        = bad_data_double;
         pr_corr.v_ncl[i]       = pr_corr.v_ncu[i]       = bad_data_double;
         anom_corr.v_ncl[i]     = anom_corr.v_ncu[i]     = bad_data_double;
         me.v_ncl[i]            = me.v_ncu[i]            = bad_data_double;
         estdev.v_ncl[i]        = estdev.v_ncu[i]        = bad_data_double;
         continue;
      }

      //
      // Compute the critical values for the Normal or Student's-T distribution
      // based on the sample size
      //
      if(n >= large_sample_threshold) {
         cv_normal_l = normal_cdf_inv(alpha[i]/2.0, 0.0, 1.0);
         cv_normal_u = normal_cdf_inv(1.0 - (alpha[i]/2.0), 0.0, 1.0);
      }
      //
      // If the number of samples is less than the large sample threshold,
      // use the T-distribution
      //
      else {
         cv_normal_l = students_t_cdf_inv(alpha[i]/2.0, n-1);
         cv_normal_u = students_t_cdf_inv(1.0 - (alpha[i]/2.0), n-1);
      }

      //
      // Compute the critical values for the Chi Squared distribution
      //
      cv_chi2_l = chi2_cdf_inv(alpha[i]/2.0, n-1);
      cv_chi2_u = chi2_cdf_inv(1.0 - (alpha[i]/2.0), n-1);

      //
      // Compute confidence interval for forecast mean using VIF
      //
      v = fbar.vif*fstdev.v*fstdev.v;
      fbar.v_ncl[i] = fbar.v + cv_normal_l*sqrt(v)/sqrt((double) n);
      fbar.v_ncu[i] = fbar.v + cv_normal_u*sqrt(v)/sqrt((double) n);

      //
      // Compute confidence interval for forecast standard deviation,
      // assuming normality of the forecast values
      //
      v = (n-1)*fstdev.v*fstdev.v/cv_chi2_u;
      if(v < 0) fstdev.v_ncl[i] = bad_data_double;
      else      fstdev.v_ncl[i] = sqrt(v);

      v = (n-1)*fstdev.v*fstdev.v/cv_chi2_l;
      if(v < 0) fstdev.v_ncu[i] = bad_data_double;
      else      fstdev.v_ncu[i] = sqrt(v);

      //
      // Compute confidence interval for observation mean using VIF
      //
      v = obar.vif*ostdev.v*ostdev.v;
      obar.v_ncl[i] = obar.v + cv_normal_l*sqrt(v)/sqrt((double) n);
      obar.v_ncu[i] = obar.v + cv_normal_u*sqrt(v)/sqrt((double) n);

      //
      // Compute confidence interval for observation standard deviation
      // assuming normality of the observation values
      //
      v = (n-1)*ostdev.v*ostdev.v/cv_chi2_u;
      if(v < 0) ostdev.v_ncl[i] = bad_data_double;
      else      ostdev.v_ncl[i] = sqrt(v);

      v = (n-1)*ostdev.v*ostdev.v/cv_chi2_l;
      if(v < 0) ostdev.v_ncu[i] = bad_data_double;
      else      ostdev.v_ncu[i] = sqrt(v);

      //
      // Compute confidence interval for the pearson correlation coefficient
      //
      if(is_bad_data(pr_corr.v) || n <= 3 ||
         is_eq(pr_corr.v, 1.0)  || is_eq(pr_corr.v, -1.0)) {
         pr_corr.v_ncl[i] = bad_data_double;
         pr_corr.v_ncu[i] = bad_data_double;
      }
      else {
         v = 0.5*log((1 + pr_corr.v)/(1 - pr_corr.v));
         cl = v + cv_normal_l/sqrt((double) (n-3));
         cu = v + cv_normal_u/sqrt((double) (n-3));
         pr_corr.v_ncl[i] = (pow(vx_math_e, 2*cl) - 1)/(pow(vx_math_e, 2*cl) + 1);
         pr_corr.v_ncu[i] = (pow(vx_math_e, 2*cu) - 1)/(pow(vx_math_e, 2*cu) + 1);
      }

      //
      // Compute confidence interval for the anomaly correlation coefficient
      //
      if(is_bad_data(anom_corr.v) || n <= 3 ||
         is_eq(anom_corr.v, 1.0)  || is_eq(anom_corr.v, -1.0)) {
         anom_corr.v_ncl[i] = bad_data_double;
         anom_corr.v_ncu[i] = bad_data_double;
      }
      else {
         v = 0.5*log((1 + anom_corr.v)/(1 - anom_corr.v));
         cl = v + cv_normal_l/sqrt((double) (n-3));
         cu = v + cv_normal_u/sqrt((double) (n-3));
         anom_corr.v_ncl[i] = (pow(vx_math_e, 2*cl) - 1)/(pow(vx_math_e, 2*cl) + 1);
         anom_corr.v_ncu[i] = (pow(vx_math_e, 2*cu) - 1)/(pow(vx_math_e, 2*cu) + 1);
      }

      //
      // Compute confidence interval for mean error using VIF
      //
      v = me.vif*estdev.v*estdev.v;
      me.v_ncl[i] = me.v + cv_normal_l*sqrt(v)/sqrt((double) n);
      me.v_ncu[i] = me.v + cv_normal_u*sqrt(v)/sqrt((double) n);

      //
      // Compute confidence interval for the error standard deviation
      //
      v = (n-1)*estdev.v*estdev.v/cv_chi2_u;
      if(v < 0) estdev.v_ncl[i] = bad_data_double;
      else      estdev.v_ncl[i] = sqrt(v);

      v = (n-1)*estdev.v*estdev.v/cv_chi2_l;
      if(v < 0) estdev.v_ncu[i] = bad_data_double;
      else      estdev.v_ncu[i] = sqrt(v);

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

double CNTInfo::get_stat(const char *stat_name) {
   double v = bad_data_double;

        if(strcmp(stat_name, "TOTAL"           ) == 0) v = n;
   else if(strcmp(stat_name, "FBAR"            ) == 0) v = fbar.v;
   else if(strcmp(stat_name, "FSTDEV"          ) == 0) v = fstdev.v;
   else if(strcmp(stat_name, "OBAR"            ) == 0) v = obar.v;
   else if(strcmp(stat_name, "OSTDEV"          ) == 0) v = ostdev.v;
   else if(strcmp(stat_name, "PR_CORR"         ) == 0) v = pr_corr.v;
   else if(strcmp(stat_name, "SP_CORR"         ) == 0) v = sp_corr.v;
   else if(strcmp(stat_name, "KT_CORR"         ) == 0) v = kt_corr.v;
   else if(strcmp(stat_name, "RANKS"           ) == 0) v = n_ranks;
   else if(strcmp(stat_name, "FRANK_TIES"      ) == 0) v = frank_ties;
   else if(strcmp(stat_name, "ORANK_TIES"      ) == 0) v = orank_ties;
   else if(strcmp(stat_name, "ME"              ) == 0) v = me.v;
   else if(strcmp(stat_name, "ESTDEV"          ) == 0) v = estdev.v;
   else if(strcmp(stat_name, "MBIAS"           ) == 0) v = mbias.v;
   else if(strcmp(stat_name, "MAE"             ) == 0) v = mae.v;
   else if(strcmp(stat_name, "MSE"             ) == 0) v = mse.v;
   else if(strcmp(stat_name, "BCMSE"           ) == 0) v = bcmse.v;
   else if(strcmp(stat_name, "RMSE"            ) == 0) v = rmse.v;
   else if(strcmp(stat_name, "E10"             ) == 0) v = e10.v;
   else if(strcmp(stat_name, "E25"             ) == 0) v = e25.v;
   else if(strcmp(stat_name, "E50"             ) == 0) v = e50.v;
   else if(strcmp(stat_name, "E75"             ) == 0) v = e75.v;
   else if(strcmp(stat_name, "E90"             ) == 0) v = e90.v;
   else if(strcmp(stat_name, "EIQR"            ) == 0) v = eiqr.v;
   else if(strcmp(stat_name, "MAD  "           ) == 0) v = mad.v;
   else if(strcmp(stat_name, "ANOM_CORR"       ) == 0) v = anom_corr.v;
   else if(strcmp(stat_name, "ME2"             ) == 0) v = me2.v;
   else if(strcmp(stat_name, "MSESS"           ) == 0) v = msess.v;
   else if(strcmp(stat_name, "RMSFA"           ) == 0) v = rmsfa.v;
   else if(strcmp(stat_name, "RMSOA"           ) == 0) v = rmsoa.v;
   else if(strcmp(stat_name, "ANOM_CORR_UNCNTR") == 0) v = anom_corr_uncntr.v;
   else {
      mlog << Error << "\nCNTInfo::get_stat() -> "
           << "unknown continuous statistic name \"" << stat_name
           << "\"!\n\n";
      exit(1);
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////
//
// Code for class SL1L2Info
//
////////////////////////////////////////////////////////////////////////

SL1L2Info::SL1L2Info() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

SL1L2Info::~SL1L2Info() {
   clear();
}

////////////////////////////////////////////////////////////////////////

SL1L2Info::SL1L2Info(const SL1L2Info &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

SL1L2Info & SL1L2Info::operator=(const SL1L2Info &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

SL1L2Info & SL1L2Info::operator+=(const SL1L2Info &c) {
   SL1L2Info s_info;

   s_info.scount  = scount + c.scount;

   if(s_info.scount > 0) {
      s_info.fbar  = (fbar*scount  + c.fbar*c.scount) /s_info.scount;
      s_info.obar  = (obar*scount  + c.obar*c.scount) /s_info.scount;
      s_info.fobar = (fobar*scount + c.fobar*c.scount)/s_info.scount;
      s_info.ffbar = (ffbar*scount + c.ffbar*c.scount)/s_info.scount;
      s_info.oobar = (oobar*scount + c.oobar*c.scount)/s_info.scount;

      if(is_bad_data(mae) || is_bad_data(c.mae)) {
         s_info.mae = bad_data_double;
      }
      else {
         s_info.mae = (mae*scount + c.mae*c.scount)/s_info.scount;
      }
   }

   s_info.sacount  = sacount + c.sacount;

   if(s_info.sacount > 0) {
      s_info.fabar  = (fabar*sacount  + c.fabar*c.sacount) /s_info.sacount;
      s_info.oabar  = (oabar*sacount  + c.oabar*c.sacount) /s_info.sacount;
      s_info.foabar = (foabar*sacount + c.foabar*c.sacount)/s_info.sacount;
      s_info.ffabar = (ffabar*sacount + c.ffabar*c.sacount)/s_info.sacount;
      s_info.ooabar = (ooabar*sacount + c.ooabar*c.sacount)/s_info.sacount;

      if(is_bad_data(mae) || is_bad_data(c.mae)) {
         s_info.mae = bad_data_double;
      }
      else {
         s_info.mae = (mae*sacount + c.mae*c.sacount)/s_info.sacount;
      }
   }

   assign(s_info);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void SL1L2Info::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void SL1L2Info::zero_out() {

   // SL1L2 Quantities
   fbar    = obar   = 0.0;
   fobar   = ffbar  = oobar  = 0.0;
   scount  = 0;

   // SAL1L2 Quantities
   fabar   = oabar  = 0.0;
   foabar  = ffabar = ooabar = 0.0;
   sacount = 0;

   mae     = 0.0;

   return;
}

////////////////////////////////////////////////////////////////////////

void SL1L2Info::clear() {

   fthresh.clear();
   othresh.clear();
   logic = SetLogic_None;

   zero_out();

   return;
}

////////////////////////////////////////////////////////////////////////

void SL1L2Info::assign(const SL1L2Info &c) {

   clear();

   fthresh = c.fthresh;
   othresh = c.othresh;
   logic   = c.logic;

   // SL1L2 Quantities
   fbar    = c.fbar;
   obar    = c.obar;
   fobar   = c.fobar;
   ffbar   = c.ffbar;
   oobar   = c.oobar;
   scount  = c.scount;

   // SAL1L2 Quantities
   fabar   = c.fabar;
   oabar   = c.oabar;
   foabar  = c.foabar;
   ffabar  = c.ffabar;
   ooabar  = c.ooabar;
   sacount = c.sacount;

   mae     = c.mae;

   return;
}

////////////////////////////////////////////////////////////////////////

void SL1L2Info::set(const PairDataPoint &pd_all) {
   int i;
   double f, o, c, wgt, wgt_sum;
   PairDataPoint pd;

   // Check for mismatch
   if(pd_all.f_na.n() != pd_all.o_na.n()) {
      mlog << Error << "\nSL1L2Info::set() -> "
           << "forecast and observation count mismatch ("
           << pd_all.f_na.n() << " != " << pd_all.o_na.n()
           << ")\n\n";
      exit(1);
   }

   // Initialize
   zero_out();

   // Apply continuous filtering thresholds to subset pairs
   pd = pd_all.subset_pairs_cnt_thresh(fthresh, othresh, logic);

   // Check for no matched pairs to process
   if(pd.n_obs == 0) return;

   // Get the sum of the weights
   wgt_sum = pd.wgt_na.sum();

   // Loop through the pair data and compute sums
   for(i=0; i<pd.n_obs; i++) {

      f   = pd.f_na[i];
      o   = pd.o_na[i];
      c   = pd.cmn_na[i];
      wgt = pd.wgt_na[i]/wgt_sum;

      // Skip bad data values in the forecast or observation fields
      if(is_bad_data(f) || is_bad_data(o)) continue;

      // SL1L2 sums
      fbar  += wgt*f;
      obar  += wgt*o;
      fobar += wgt*f*o;
      ffbar += wgt*f*f;
      oobar += wgt*o*o;
      mae   += wgt*fabs(f-o);
      scount++;

      // SAL1L2 sums
      if(!is_bad_data(c)) {
         fabar  += wgt*(f-c);
         oabar  += wgt*(o-c);
         foabar += wgt*(f-c)*(o-c);
         ffabar += wgt*(f-c)*(f-c);
         ooabar += wgt*(o-c)*(o-c);
         sacount++;
      }
   }

   if(scount == 0) {
      mlog << Error << "\nSL1L2Info::set() -> "
           << "count is zero!\n\n";
      exit(1);
   }

   // Compute the mean SAL1L2 values
   if(sacount == 0) {
      fabar  = bad_data_double;
      oabar  = bad_data_double;
      foabar = bad_data_double;
      ffabar = bad_data_double;
      ooabar = bad_data_double;
   }

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class VL1L2Info
//
////////////////////////////////////////////////////////////////////////

VL1L2Info::VL1L2Info() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

VL1L2Info::~VL1L2Info() {
   clear();
}

////////////////////////////////////////////////////////////////////////

VL1L2Info::VL1L2Info(const VL1L2Info &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

VL1L2Info & VL1L2Info::operator=(const VL1L2Info &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

VL1L2Info & VL1L2Info::operator+=(const VL1L2Info &c) {
   VL1L2Info v_info;

   v_info.vcount  = vcount + c.vcount;

   if(v_info.vcount > 0) {
      v_info.uf_bar      = (uf_bar*vcount      + c.uf_bar*c.vcount)     /v_info.vcount;
      v_info.vf_bar      = (vf_bar*vcount      + c.vf_bar*c.vcount)     /v_info.vcount;
      v_info.uo_bar      = (uo_bar*vcount      + c.uo_bar*c.vcount)     /v_info.vcount;
      v_info.vo_bar      = (vo_bar*vcount      + c.vo_bar*c.vcount)     /v_info.vcount;
      v_info.uvfo_bar    = (uvfo_bar*vcount    + c.uvfo_bar*c.vcount)   /v_info.vcount;
      v_info.uvff_bar    = (uvff_bar*vcount    + c.uvff_bar*c.vcount)   /v_info.vcount;
      v_info.uvoo_bar    = (uvoo_bar*vcount    + c.uvoo_bar*c.vcount)   /v_info.vcount;
      v_info.f_speed_bar = (f_speed_bar*vcount + c.f_speed_bar*c.vcount)/v_info.vcount;
      v_info.o_speed_bar = (o_speed_bar*vcount + c.o_speed_bar*c.vcount)/v_info.vcount;
   }

   v_info.vacount  = vacount + c.vacount;

   if(v_info.vacount > 0) {
      v_info.ufa_bar   = (ufa_bar*vacount   + c.ufa_bar*c.vacount)  /v_info.vacount;
      v_info.vfa_bar   = (vfa_bar*vacount   + c.vfa_bar*c.vacount)  /v_info.vacount;
      v_info.uoa_bar   = (uoa_bar*vacount   + c.uoa_bar*c.vacount)  /v_info.vacount;
      v_info.voa_bar   = (voa_bar*vacount   + c.voa_bar*c.vacount)  /v_info.vacount;
      v_info.uvfoa_bar = (uvfoa_bar*vacount + c.uvfoa_bar*c.vacount)/v_info.vacount;
      v_info.uvffa_bar = (uvffa_bar*vacount + c.uvffa_bar*c.vacount)/v_info.vacount;
      v_info.uvooa_bar = (uvooa_bar*vacount + c.uvooa_bar*c.vacount)/v_info.vacount;
   }

   v_info.calc_ncep_stats();

   assign(v_info);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void VL1L2Info::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////


void VL1L2Info::zero_out() {

   //
   // VL1L2 Quantities
   //

   uf_bar      = 0.0;
   vf_bar      = 0.0;
   uo_bar      = 0.0;
   vo_bar      = 0.0;
   uvfo_bar    = 0.0;
   uvff_bar    = 0.0;
   uvoo_bar    = 0.0;
   f_speed_bar = 0.0;
   o_speed_bar = 0.0;

   f_bar       = 0.0;
   o_bar       = 0.0;
   me          = 0.0;
   mse         = 0.0;
   rmse        = 0.0;
   speed_bias  = 0.0;

   FBAR        = 0.0;
   OBAR        = 0.0;

   FS_RMS      = 0.0;
   OS_RMS      = 0.0;

    MSVE       = 0.0;
   RMSVE       = 0.0;

   FSTDEV      = 0.0;
   OSTDEV      = 0.0;

   FDIR        = 0.0;
   ODIR        = 0.0;

   FBAR_SPEED  = 0.0;
   OBAR_SPEED  = 0.0;

   VDIFF_SPEED = 0.0;
   VDIFF_DIR   = 0.0;

   SPEED_ERR   = 0.0;
   SPEED_ABSERR = 0.0;

   DIR_ERR     = 0.0;
   DIR_ABSERR  = 0.0;

   vcount      = 0;

   //
   // VAL1L2 Quantities
   //

   ufa_bar     = 0.0;
   vfa_bar     = 0.0;
   uoa_bar     = 0.0;
   voa_bar     = 0.0;
   uvfoa_bar   = 0.0;
   uvffa_bar   = 0.0;
   uvooa_bar   = 0.0;

   vacount     = 0;

   return;
}

////////////////////////////////////////////////////////////////////////

void VL1L2Info::clear() {

   fthresh.clear();
   othresh.clear();
   logic = SetLogic_None;

   zero_out();

   return;
}

////////////////////////////////////////////////////////////////////////

void VL1L2Info::assign(const VL1L2Info &c) {

   clear();

   fthresh = c.fthresh;
   othresh = c.othresh;
   logic   = c.logic;

   // VL1L2 Quantities
   uf_bar      = c.uf_bar;
   vf_bar      = c.vf_bar;
   uo_bar      = c.uo_bar;
   vo_bar      = c.vo_bar;
   uvfo_bar    = c.uvfo_bar;
   uvff_bar    = c.uvff_bar;
   uvoo_bar    = c.uvoo_bar;
   f_speed_bar = c.f_speed_bar;
   o_speed_bar = c.o_speed_bar;

   f_bar     = c.f_bar;
   o_bar     = c.o_bar;

   me        = c.me;

   mse       = c.mse;

   speed_bias = c.speed_bias;

   vcount    = c.vcount;

   // VAL1L2 Quantities
   ufa_bar   = c.ufa_bar;
   vfa_bar   = c.vfa_bar;
   uoa_bar   = c.uoa_bar;
   voa_bar   = c.voa_bar;
   uvfoa_bar = c.uvfoa_bar;
   uvffa_bar = c.uvffa_bar;
   uvooa_bar = c.uvooa_bar;
   vacount   = c.vacount;

      //
      //  NCEP stats
      //

   FBAR = c.FBAR;
   OBAR = c.OBAR;

   FS_RMS = c.FS_RMS;
   OS_RMS = c.OS_RMS;

    MSVE = c.MSVE;
   RMSVE = c.RMSVE;

   FSTDEV = c.FSTDEV;
   OSTDEV = c.OSTDEV;

   FDIR = c.FDIR;
   ODIR = c.ODIR;

   FBAR_SPEED = c.FBAR_SPEED;
   OBAR_SPEED = c.OBAR_SPEED;

   VDIFF_SPEED = c.VDIFF_SPEED;
   VDIFF_DIR = c.VDIFF_DIR;

   SPEED_ERR = c.SPEED_ERR;
   SPEED_ABSERR = c.SPEED_ABSERR;

   DIR_ERR = c.DIR_ERR;
   DIR_ABSERR = c.DIR_ABSERR;

   return;
}

////////////////////////////////////////////////////////////////////////


void VL1L2Info::calc_ncep_stats() {
   double u_diff, v_diff;
   int n = vcount;

   u_diff       = uf_bar - uo_bar;
   v_diff       = vf_bar - vo_bar;

   FBAR         = f_speed_bar;
   OBAR         = o_speed_bar;

   FS_RMS       = sqrt(uvff_bar);
   OS_RMS       = sqrt(uvoo_bar);

   MSVE         = uvff_bar - 2.0*uvfo_bar + uvoo_bar;

   RMSVE        = sqrt(MSVE);

   FSTDEV       = compute_stdev(f_speed_bar*n, uvff_bar*n, n);
   OSTDEV       = compute_stdev(o_speed_bar*n, uvoo_bar*n, n);

   FDIR         = convert_u_v_to_wdir(uf_bar, vf_bar);
   ODIR         = convert_u_v_to_wdir(uo_bar, vo_bar);

   FBAR_SPEED   = convert_u_v_to_wind(uf_bar, vf_bar);
   OBAR_SPEED   = convert_u_v_to_wind(uo_bar, vo_bar);

   VDIFF_SPEED  = convert_u_v_to_wind(u_diff, v_diff);

   VDIFF_DIR    = convert_u_v_to_wdir(u_diff, v_diff);

   SPEED_ERR    = FBAR_SPEED - OBAR_SPEED;

   SPEED_ABSERR = fabs(SPEED_ERR);

   DIR_ERR      = atan2d(vf_bar*uo_bar - uf_bar*vo_bar, uf_bar*uo_bar + vf_bar*vo_bar);

   DIR_ABSERR   = fabs(DIR_ERR);

   return;
}

////////////////////////////////////////////////////////////////////////

void VL1L2Info::set(const PairDataPoint &pd_u_all,
                    const PairDataPoint &pd_v_all) {
   int i;
   double uf, vf, uo, vo, uc, vc, wgt, wgt_sum;
   double u_diff, v_diff;
   PairDataPoint pd_u, pd_v;

   // Initialize
   zero_out();

   // Check that the number of pairs are the same
   if(pd_u_all.f_na.n() != pd_u_all.o_na.n() ||
      pd_u_all.f_na.n() != pd_v_all.f_na.n() ||
      pd_v_all.f_na.n() != pd_v_all.o_na.n()) {
      mlog << Error << "\nVL1L2Info::set() -> "
           << "unequal number of UGRD and VGRD pairs ("
           << pd_u_all.f_na.n() << " != " << pd_u_all.o_na.n()
           << ")\n\n";
      exit(1);
   }

   // Apply wind speed filtering thresholds to subset pairs
   subset_wind_pairs(pd_u_all, pd_v_all, fthresh, othresh, logic, pd_u, pd_v);

   // Check for no matched pairs to process
   if(pd_u.n_obs == 0) return;

   // Get the sum of the weights
   wgt_sum = pd_u.wgt_na.sum();

   // Loop through the filtered pair data compute partial sums
   for(i=0; i<pd_u.f_na.n(); i++)  {

      // Retrieve the U,V values
      uf = pd_u.f_na[i];
      vf = pd_v.f_na[i];
      uo = pd_u.o_na[i];
      vo = pd_v.o_na[i];
      uc = pd_u.cmn_na[i];
      vc = pd_v.cmn_na[i];

      u_diff = uf - uo;
      v_diff = vf - vo;

      wgt = pd_u.wgt_na[i]/wgt_sum;

      // VL1L2 sums
      vcount     += 1;

      uf_bar     += wgt*uf;
      vf_bar     += wgt*vf;
      uo_bar     += wgt*uo;
      vo_bar     += wgt*vo;

      uvfo_bar   += wgt*(uf*uo + vf*vo);
      uvff_bar   += wgt*(uf*uf + vf*vf);
      uvoo_bar   += wgt*(uo*uo + vo*vo);


      f_bar      += wgt*sqrt(uf*uf + vf*vf);
      o_bar      += wgt*sqrt(uo*uo + vo*vo);

      me         += wgt*sqrt(u_diff*u_diff + v_diff*v_diff);

      mse        += wgt*(u_diff*u_diff + v_diff*v_diff);

      speed_bias += wgt*(sqrt(uf*uf + vf*vf) - sqrt(uo*uo + vo*vo));

         //
         //  new stuff from vector stats whitepaper
         //

      f_speed_bar += wgt*sqrt(uf*uf + vf*vf);
      o_speed_bar += wgt*sqrt(uo*uo + vo*vo);

      // VAL1L2 sums
      if(!is_bad_data(uc) && !is_bad_data(vc)) {
         vacount   += 1;
         ufa_bar   += wgt*(uf-uc);
         vfa_bar   += wgt*(vf-vc);
         uoa_bar   += wgt*(uo-uc);
         voa_bar   += wgt*(vo-vc);
         uvfoa_bar += wgt*((uf-uc)*(uo-uc) + (vf-vc)*(vo-vc));
         uvffa_bar += wgt*((uf-uc)*(uf-uc) + (vf-vc)*(vf-vc));
         uvooa_bar += wgt*((uo-uc)*(uo-uc) + (vo-vc)*(vo-vc));
      }

   }  // end for i

   if(vcount > 0) calc_ncep_stats();

   // Check for 0 points
   if(vcount == 0) {

      uf_bar        = bad_data_double;
      vf_bar        = bad_data_double;
      uo_bar        = bad_data_double;
      vo_bar        = bad_data_double;
      uvfo_bar      = bad_data_double;
      uvff_bar      = bad_data_double;
      uvoo_bar      = bad_data_double;

      me            = bad_data_double;
      mse           = bad_data_double;
      rmse          = bad_data_double;
      speed_bias    = bad_data_double;

      FBAR          = bad_data_double;
      OBAR          = bad_data_double;

      FS_RMS        = bad_data_double;
      OS_RMS        = bad_data_double;

       MSVE         = bad_data_double;
      RMSVE         = bad_data_double;

      FSTDEV        = bad_data_double;
      OSTDEV        = bad_data_double;

      FDIR          = bad_data_double;
      ODIR          = bad_data_double;

      FBAR_SPEED    = bad_data_double;
      OBAR_SPEED    = bad_data_double;

      VDIFF_SPEED   = bad_data_double;
      VDIFF_DIR     = bad_data_double;

      SPEED_ERR     = bad_data_double;
      SPEED_ABSERR  = bad_data_double;

      DIR_ERR       = bad_data_double;
      DIR_ABSERR    = bad_data_double;

   } else {
      rmse          = sqrt(mse);
   }

   if(vacount == 0) {
      ufa_bar   = bad_data_double;
      vfa_bar   = bad_data_double;
      uoa_bar   = bad_data_double;
      voa_bar   = bad_data_double;
      uvfoa_bar = bad_data_double;
      uvffa_bar = bad_data_double;
      uvooa_bar = bad_data_double;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

double VL1L2Info::get_stat(const char *stat_name) {
   double v = bad_data_double;

        if(strcmp(stat_name, "TOTAL"       ) == 0) v = vcount;
   else if(strcmp(stat_name, "FBAR"        ) == 0) v = FBAR;
   else if(strcmp(stat_name, "OBAR"        ) == 0) v = OBAR;
   else if(strcmp(stat_name, "FS_RMS"      ) == 0) v = FS_RMS;
   else if(strcmp(stat_name, "OS_RMS"      ) == 0) v = OS_RMS;
   else if(strcmp(stat_name, "MSVE"        ) == 0) v = MSVE;
   else if(strcmp(stat_name, "RMSVE"       ) == 0) v = RMSVE;
   else if(strcmp(stat_name, "FSTDEV"      ) == 0) v = FSTDEV;
   else if(strcmp(stat_name, "OSTDEV"      ) == 0) v = OSTDEV;
   else if(strcmp(stat_name, "FDIR"        ) == 0) v = FDIR;
   else if(strcmp(stat_name, "ODIR"        ) == 0) v = ODIR;
   else if(strcmp(stat_name, "FBAR_SPEED"  ) == 0) v = FBAR_SPEED;
   else if(strcmp(stat_name, "OBAR_SPEED"  ) == 0) v = OBAR_SPEED;
   else if(strcmp(stat_name, "VDIFF_SPEED" ) == 0) v = VDIFF_SPEED;
   else if(strcmp(stat_name, "VDIFF_DIR"   ) == 0) v = VDIFF_DIR;
   else if(strcmp(stat_name, "SPEED_ERR"   ) == 0) v = SPEED_ERR;
   else if(strcmp(stat_name, "SPEED_ABSERR") == 0) v = SPEED_ABSERR;
   else if(strcmp(stat_name, "DIR_ERR"     ) == 0) v = DIR_ERR;
   else if(strcmp(stat_name, "DIR_ABSERR"  ) == 0) v = DIR_ABSERR;
   else {
      mlog << Error << "\nVL1L2Info::get_stat() -> "
           << "unknown continuous statistic name \"" << stat_name
           << "\"!\n\n";
      exit(1);
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////
//
// Code for class NBRCTSInfo
//
////////////////////////////////////////////////////////////////////////

NBRCTSInfo::NBRCTSInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

NBRCTSInfo::~NBRCTSInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

NBRCTSInfo::NBRCTSInfo(const NBRCTSInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

NBRCTSInfo & NBRCTSInfo::operator=(const NBRCTSInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void NBRCTSInfo::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCTSInfo::clear() {

   fthresh.clear();
   othresh.clear();
   cthresh.clear();
   cts_info.clear();
   nbr_wdth = bad_data_int;

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCTSInfo::assign(const NBRCTSInfo &c) {

   clear();

   fthresh  = c.fthresh;
   othresh  = c.othresh;
   cthresh  = c.cthresh;
   cts_info = c.cts_info;
   nbr_wdth = c.nbr_wdth;

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCTSInfo::allocate_n_alpha(int i) {

   cts_info.allocate_n_alpha(i);

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class NBRCNTInfo
//
////////////////////////////////////////////////////////////////////////

NBRCNTInfo::NBRCNTInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

NBRCNTInfo::~NBRCNTInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

NBRCNTInfo::NBRCNTInfo(const NBRCNTInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

NBRCNTInfo & NBRCNTInfo::operator=(const NBRCNTInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

NBRCNTInfo & NBRCNTInfo::operator+=(const NBRCNTInfo &c) {
   NBRCNTInfo n_info;
   double den;

   n_info.sl1l2_info.scount = sl1l2_info.scount + c.sl1l2_info.scount;

   if(n_info.sl1l2_info.scount) {

      //
      // Aggregate FBS as a weighted average
      //
      if(is_bad_data(sl1l2_info.scount*fbs.v) ||
         is_bad_data(c.sl1l2_info.scount*c.fbs.v)) {
         n_info.fbs.v = bad_data_double;
      }
      else {
         n_info.fbs.v = (sl1l2_info.scount*fbs.v + c.sl1l2_info.scount*c.fbs.v) /
                         n_info.sl1l2_info.scount;
      }

      //
      // Aggregate the denominator for FSS as a weighted average and then
      // recompute FSS using the aggregated FBS and denominator
      //
      if(is_bad_data(fbs.v) || is_bad_data(c.fbs.v) ||
         is_bad_data(fss.v) || is_bad_data(c.fss.v) ||
         is_eq(fss.v, 1.0)  || is_eq(c.fss.v, 1.0)) {
         den = bad_data_double;
      }
      else {
         den = (  sl1l2_info.scount *  fbs.v / (1.0 - fss.v  ) +
                c.sl1l2_info.scount *c.fbs.v / (1.0 - c.fss.v)) /
                n_info.sl1l2_info.scount;
      }

      //
      // If FSS cannot be aggregated numerically, just keep its current value.
      // Or in the case of bad data, use the input FSS value.
      //
      if(is_bad_data(den) || is_eq(den, 0.0)) {
         n_info.fss.v = (is_bad_data(fss.v) ? c.fss.v : fss.v);

      }
      else {
         n_info.fss.v = 1.0 - n_info.fbs.v / den;
      }

      //
      // Aggregate F_RATE and O_RATE as weighted averages
      //
      n_info.f_rate.v = (sl1l2_info.scount*f_rate.v + c.sl1l2_info.scount*c.f_rate.v) /
                         n_info.sl1l2_info.scount;
      n_info.o_rate.v = (sl1l2_info.scount*o_rate.v + c.sl1l2_info.scount*c.o_rate.v) /
                         n_info.sl1l2_info.scount;

      //
      // Recompute AFSS and UFSS using the aggregated rates
      //
      n_info.afss.v = compute_afss(n_info.f_rate.v, n_info.o_rate.v);
      n_info.ufss.v = compute_ufss(n_info.o_rate.v);
   }

   assign(n_info);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void NBRCNTInfo::init_from_scratch() {

   alpha = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCNTInfo::clear() {

   n_alpha = 0;
   if(alpha) { delete [] alpha; alpha = (double *) 0; }

   fthresh.clear();
   othresh.clear();

   fbs.clear();
   fss.clear();
   afss.clear();
   ufss.clear();
   f_rate.clear();
   o_rate.clear();
   sl1l2_info.clear();
   nbr_wdth = bad_data_int;

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCNTInfo::assign(const NBRCNTInfo &c) {

   clear();

   fthresh    = c.fthresh;
   othresh    = c.othresh;
   fbs        = c.fbs;
   fss        = c.fss;
   afss       = c.afss;
   ufss       = c.ufss;
   f_rate     = c.f_rate;
   o_rate     = c.o_rate;
   sl1l2_info = c.sl1l2_info;
   nbr_wdth   = c.nbr_wdth;

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCNTInfo::allocate_n_alpha(int i) {

   n_alpha = i;

   if(n_alpha > 0) {

      alpha = new double [n_alpha];

      if(!alpha) {
         mlog << Error << "\nNBRCNTInfo::allocate_n_alpha() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }

      fbs.allocate_n_alpha(i);
      fss.allocate_n_alpha(i);
      afss.allocate_n_alpha(i);
      ufss.allocate_n_alpha(i);
      f_rate.allocate_n_alpha(i);
      o_rate.allocate_n_alpha(i);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void NBRCNTInfo::compute_stats() {
   double num, den;

   //
   // Compute FBS
   //
   fbs.v = sl1l2_info.ffbar + sl1l2_info.oobar - 2.0*sl1l2_info.fobar;

   //
   // Compute FSS
   //
   num = fbs.v;
   den = sl1l2_info.ffbar + sl1l2_info.oobar;

   if(is_eq(den, 0.0)) fss.v = bad_data_double;
   else                fss.v = 1.0 - (num/den);

   //
   // Compute F_RATE and O_RATE
   //
   afss.v = compute_afss(f_rate.v, o_rate.v);
   ufss.v = compute_ufss(o_rate.v);

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class ISCInfo
//
////////////////////////////////////////////////////////////////////////

ISCInfo::ISCInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

ISCInfo::~ISCInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

ISCInfo::ISCInfo(const ISCInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

ISCInfo & ISCInfo::operator=(const ISCInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::init_from_scratch() {

   mse_scale = (double *) 0;
   isc_scale = (double *) 0;
   fen_scale = (double *) 0;
   oen_scale = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::clear() {

   cts.zero_out();
   fthresh.clear();
   othresh.clear();

   mse   = isc   = bad_data_double;
   fen   = oen   = bad_data_double;
   baser = fbias = bad_data_double;

   tile_dim = bad_data_int;
   tile_xll = bad_data_int;
   tile_yll = bad_data_int;

   n_scale  = 0;
   total    = 0;

   if(mse_scale) { delete [] mse_scale; mse_scale = (double *) 0; }
   if(isc_scale) { delete [] isc_scale; isc_scale = (double *) 0; }
   if(fen_scale) { delete [] fen_scale; fen_scale = (double *) 0; }
   if(oen_scale) { delete [] oen_scale; oen_scale = (double *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::zero_out() {
   int i;

   cts.zero_out();
   total = 0;

   mse   = isc   = 0.0;
   fen   = oen   = 0.0;
   baser = fbias = 0.0;

   tile_dim = bad_data_int;
   tile_xll = bad_data_int;
   tile_yll = bad_data_int;


   if(n_scale > 0) {
      for(i=0; i<=n_scale; i++) {
         mse_scale[i] = 0.0;
         isc_scale[i] = 0.0;
         fen_scale[i] = 0.0;
         oen_scale[i] = 0.0;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::assign(const ISCInfo &c) {
   int i;

   clear();

   cts = c.cts;

   fthresh = c.fthresh;
   othresh  = c.othresh;

   mse   = c.mse;
   isc   = c.isc;
   fen   = c.fen;
   oen   = c.oen;
   baser = c.baser;
   fbias = c.fbias;
   total = c.total;

   tile_dim = c.tile_dim;
   tile_xll = c.tile_xll;
   tile_yll = c.tile_yll;

   allocate_n_scale(c.n_scale);

   for(i=0; i<n_scale; i++) {
      mse_scale[i] = c.mse_scale[i];
      isc_scale[i] = c.isc_scale[i];
      fen_scale[i] = c.fen_scale[i];
      oen_scale[i] = c.oen_scale[i];
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::allocate_n_scale(int i) {
   int j;

   if((n_scale = i) == 0) return;

   mse_scale = new double [n_scale+1];
   isc_scale = new double [n_scale+1];
   fen_scale = new double [n_scale+1];
   oen_scale = new double [n_scale+1];

   if(!mse_scale || !isc_scale || !fen_scale || !oen_scale) {
      mlog << Error << "\nISCInfo::allocate_n_scale() -> "
           << "Memory allocation error!\n\n";
      exit(1);
   }

   // Initialize the values
   for(j=0; j<=n_scale; j++) {
      mse_scale[j] = bad_data_double;
      isc_scale[j] = bad_data_double;
      fen_scale[j] = bad_data_double;
      oen_scale[j] = bad_data_double;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::compute_isc() {
   double den;
   int i;

   // Get the Total, Base Rate, and Frequency Bias
   total = cts.n();
   fbias = cts.fbias();
   baser = cts.baser();

   // Compute the denominator for ISC
   den = fbias*baser*(1.0 - baser) + baser*(1.0 - fbias*baser);

   // Compute ISC for the whole field
   if(is_bad_data(fbias) ||
      is_bad_data(baser) ||
      is_bad_data(mse)   ||
      is_eq(den, 0.0)) isc = bad_data_double;
   else                isc = 1.0 - mse/den;

   // Compute the Intensity-Scale score for each scale
   den /= (n_scale+1);
   for(i=0; i<=n_scale; i++) {

      if(is_bad_data(fbias)        ||
         is_bad_data(baser)        ||
         is_bad_data(mse_scale[i]) ||
         is_eq(den, 0.0)) isc_scale[i] = bad_data_double;
      else                isc_scale[i] = 1.0 - mse_scale[i]/den;
   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

void ISCInfo::compute_isc(int i) {
   double den;

   // Get the Total, Base Rate, and Frequency Bias
   total = cts.n();
   fbias = cts.fbias();
   baser = cts.baser();

   // Compute the denominator for ISC
   den = fbias*baser*(1.0 - baser) + baser*(1.0 - fbias*baser);

   // Compute ISC for the whole field
   if(i < 0) {
      if(is_bad_data(fbias) ||
         is_bad_data(baser) ||
         is_bad_data(mse)   ||
         is_eq(den, 0.0)) isc = bad_data_double;
      else                isc = 1.0 - mse/den;
   }

   // Compute the Intensity-Scale score for each scale
   else {
      den /= (n_scale+1);

      if(is_bad_data(fbias)        ||
         is_bad_data(baser)        ||
         is_bad_data(mse_scale[i]) ||
         is_eq(den, 0.0)) isc_scale[i] = bad_data_double;
      else                isc_scale[i] = 1.0 - mse_scale[i]/den;
   }

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class PCTInfo
//
////////////////////////////////////////////////////////////////////////

PCTInfo::PCTInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

PCTInfo::~PCTInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

PCTInfo::PCTInfo(const PCTInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

PCTInfo & PCTInfo::operator=(const PCTInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::init_from_scratch() {

   alpha = (double *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::clear() {

   n_alpha = 0;
   if(alpha) { delete [] alpha; alpha = (double *) 0; }

   pct.zero_out();
   climo_pct.zero_out();
   fthresh.clear();
   othresh.clear();

   total = bad_data_int;
   baser.clear();
   reliability = bad_data_double;
   resolution = bad_data_double;
   uncertainty = bad_data_double;
   roc_auc = bad_data_double;
   brier.clear();
   briercl.clear();
   bss = bad_data_double;
   bss_smpl = bad_data_double;

   return;
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::assign(const PCTInfo &c) {
   int i;

   clear();

   pct       = c.pct;
   climo_pct = c.climo_pct;
   fthresh   = c.fthresh;
   othresh   = c.othresh;

   allocate_n_alpha(c.n_alpha);
   for(i=0; i<c.n_alpha; i++) { alpha[i] = c.alpha[i]; }

   total       = c.total;
   baser       = c.baser;
   reliability = c.reliability;
   resolution  = c.resolution;
   uncertainty = c.uncertainty;
   roc_auc     = c.roc_auc;
   brier       = c.brier;
   briercl     = c.briercl;
   bss         = c.bss;
   bss_smpl    = c.bss_smpl;

   return;
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::allocate_n_alpha(int i) {

   n_alpha = i;

   if(n_alpha > 0) {

      alpha = new double [n_alpha];

      if(!alpha) {
         mlog << Error << "\nPCTInfo::allocate_n() -> "
              << "Memory allocation error!\n\n";
        exit(1);
      }

      baser.allocate_n_alpha(n_alpha);
      brier.allocate_n_alpha(n_alpha);
      briercl.allocate_n_alpha(n_alpha);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::compute_stats() {

   total       = pct.n();
   baser.v     = pct.baser();
   reliability = pct.reliability();
   resolution  = pct.resolution();
   uncertainty = pct.uncertainty();
   roc_auc     = pct.roc_auc();
   brier.v     = pct.brier_score();
   briercl.v   = climo_pct.brier_score();
   bss_smpl    = pct.bss_smpl();

   //
   // Compute the brier skill score
   //
   if(is_bad_data(brier.v) || is_bad_data(briercl.v) ||
      is_eq(briercl.v, 0.0)) {
      bss = bad_data_double;
   }
   else {
      bss = 1.0 - (brier.v / briercl.v);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void PCTInfo::compute_ci() {
   int i;
   double halfwidth;

   //
   // Compute confidence intervals for each alpha value specified
   //
   for(i=0; i<n_alpha; i++) {

      compute_proportion_ci(baser.v, pct.n(), alpha[i], baser.vif,
                            baser.v_ncl[i], baser.v_ncu[i]);

      // Compute brier CI using the VIF
      if(!is_bad_data(brier.v)) {
         if(!is_bad_data(halfwidth = pct.brier_ci_halfwidth(alpha[i]))) {
            halfwidth *= sqrt(brier.vif);
            brier.v_ncl[i] = brier.v - halfwidth;
            brier.v_ncu[i] = brier.v + halfwidth;
         }
      }

      // Compute climatological brier CI using the VIF
      if(!is_bad_data(briercl.v)) {
         if(!is_bad_data(halfwidth = climo_pct.brier_ci_halfwidth(alpha[i]))) {
            halfwidth *= sqrt(briercl.vif);
            briercl.v_ncl[i] = briercl.v - halfwidth;
            briercl.v_ncu[i] = briercl.v + halfwidth;
         }
      }

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class GRADInfo
//
////////////////////////////////////////////////////////////////////////

GRADInfo::GRADInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

GRADInfo::~GRADInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

GRADInfo::GRADInfo(const GRADInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

GRADInfo & GRADInfo::operator=(const GRADInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

GRADInfo & GRADInfo::operator+=(const GRADInfo &c) {
   GRADInfo g_info;

   if(dx != c.dx || dy != c.dy) {
      mlog << Error << "\nGRADInfo::operator+=() -> "
           << "the gradient DX (" << dx << " vs " << c.dx
           << ") and DY (" << dy << " vs " << c.dy
           << ") sizes must remain constant!\n\n";
      exit(1);
   }

   g_info.dx    = dx;
   g_info.dy    = dy;
   g_info.total = total + c.total;

   if(g_info.total > 0) {
      g_info.fgbar = (fgbar*total + c.fgbar*c.total) / g_info.total;
      g_info.ogbar = (ogbar*total + c.ogbar*c.total) / g_info.total;
      g_info.mgbar = (mgbar*total + c.mgbar*c.total) / g_info.total;
      g_info.egbar = (egbar*total + c.egbar*c.total) / g_info.total;
   }

   assign(g_info);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void GRADInfo::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void GRADInfo::clear() {

   dx    = dy    = 0;
   fgbar = ogbar = 0.0;
   mgbar = egbar = 0.0;
   total = 0;

   return;
}

////////////////////////////////////////////////////////////////////////

void GRADInfo::assign(const GRADInfo &c) {

   clear();

   // Gradient sizes
   dx = c.dx;
   dy = c.dy;

   // Gradient partial sums
   fgbar = c.fgbar;
   ogbar = c.ogbar;
   mgbar = c.mgbar;
   egbar = c.egbar;
   total = c.total;

   return;
}

////////////////////////////////////////////////////////////////////////

double GRADInfo::s1() const {
   double v;

   if(is_bad_data(egbar) || is_bad_data(mgbar) || is_eq(mgbar, 0.0)) {
      v = bad_data_double;
   }
   else {
      v = 100.0 * egbar / mgbar;
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////

double GRADInfo::s1_og() const {
   double v;

   if(is_bad_data(egbar) || is_bad_data(ogbar) || is_eq(ogbar, 0.0)) {
      v = bad_data_double;
   }
   else {
      v = 100.0 * egbar / ogbar;
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////

double GRADInfo::fgog_ratio() const {
   double v;

   if(is_bad_data(fgbar) || is_bad_data(ogbar) || is_eq(ogbar, 0.0)) {
      v = bad_data_double;
   }
   else {
      v = fgbar / ogbar;
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////

void GRADInfo::set(int grad_dx, int grad_dy,
                   const NumArray &fgx_na, const NumArray &fgy_na,
                   const NumArray &ogx_na, const NumArray &ogy_na,
                   const NumArray &wgt_na) {
   int i;
   double wgt, wgt_sum;

   // Check for mismatch
   if(fgx_na.n() != fgy_na.n() ||
      fgx_na.n() != ogx_na.n() ||
      fgx_na.n() != ogy_na.n() ||
      fgx_na.n() != wgt_na.n()) {
      mlog << Error << "\nGRADInfo::set() -> "
           << "count mismatch ("
           << fgx_na.n() << ", " << fgy_na.n() << ", "
           << ogx_na.n() << ", " << ogy_na.n() << ", "
           <<   wgt_na.n() << ")\n\n";
      exit(1);
   }

   // Initialize
   clear();

   // Store the gradient size
   dx = grad_dx;
   dy = grad_dy;

   // Check for no matched pairs to process
   if(fgx_na.n() == 0) return;

   // Get the sum of the weights
   wgt_sum = wgt_na.sum();

   // Loop through the pairs and compute sums
   for(i=0; i<fgx_na.n(); i++) {

      // Skip bad data
      if(is_bad_data(fgx_na[i]) || is_bad_data(fgy_na[i]) ||
         is_bad_data(ogx_na[i]) || is_bad_data(ogy_na[i])) continue;

      // Get current weight
      wgt = wgt_na[i]/wgt_sum;

      // Gradient sums
      fgbar += wgt * (fabs(fgx_na[i]) + fabs(fgy_na[i]));
      ogbar += wgt * (fabs(ogx_na[i]) + fabs(ogy_na[i]));
      mgbar += wgt * (max(fabs(fgx_na[i]), fabs(ogx_na[i])) +
                      max(fabs(fgy_na[i]), fabs(ogy_na[i])));
      egbar += wgt * (fabs(fgx_na[i] - ogx_na[i]) +
                      fabs(fgy_na[i] - ogy_na[i]));
      total++;
   }

   if(total == 0) {
      mlog << Error << "\nGRADInfo::set() -> "
           << "count is zero!\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class DMAPInfo
//
////////////////////////////////////////////////////////////////////////

DMAPInfo::DMAPInfo() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

DMAPInfo::~DMAPInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

DMAPInfo::DMAPInfo(const DMAPInfo &c) {

   init_from_scratch();

   assign(c);
}

////////////////////////////////////////////////////////////////////////

DMAPInfo & DMAPInfo::operator=(const DMAPInfo &c) {

   if(this == &c) return(*this);

   assign(c);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void DMAPInfo::init_from_scratch() {

   clear();
   reset_options();

   return;
}

////////////////////////////////////////////////////////////////////////

void DMAPInfo::clear() {

   fthresh.clear();
   othresh.clear();
   total    = fy = oy = foy = 0;
   baddeley = hausdorff = bad_data_double;
   med_fo   = med_of = med_min = med_max = med_mean = bad_data_double;
   fom_fo   = fom_of = fom_min = fom_max = fom_mean = bad_data_double;
   zhu_fo   = zhu_of = zhu_min = zhu_max = zhu_mean = bad_data_double;
   g        = gbeta  = beta_value = bad_data_double;

   return;
}

////////////////////////////////////////////////////////////////////////

void DMAPInfo::reset_options() {
   baddeley_p = 2;                      // Exponent for lp-norm
   baddeley_max_dist = bad_data_double; // Maximum distance constant
   fom_alpha = 0.1;                     // FOM Alpha
   zhu_weight = 0.5;                    // Zhu Weight
   beta_value = bad_data_double;        // G-Beta Value

   return;
}

////////////////////////////////////////////////////////////////////////

void DMAPInfo::assign(const DMAPInfo &c) {

   clear();

   fthresh = c.fthresh;
   othresh = c.othresh;

   total = c.total;
   fy    = c.fy;
   oy    = c.oy;
   foy   = c.foy;

   baddeley  = c.baddeley;
   hausdorff = c.hausdorff;

   med_fo   = c.med_fo;
   med_of   = c.med_of;
   med_min  = c.med_min;
   med_max  = c.med_max;
   med_mean = c.med_mean;

   fom_fo   = c.fom_fo;
   fom_of   = c.fom_of;
   fom_min  = c.fom_min;
   fom_max  = c.fom_max;
   fom_mean = c.fom_mean;

   zhu_fo   = c.zhu_fo;
   zhu_of   = c.zhu_of;
   zhu_min  = c.zhu_min;
   zhu_max  = c.zhu_max;
   zhu_mean = c.zhu_mean;

   g     = c.g;
   gbeta = c.gbeta;

   baddeley_p = c.baddeley_p;
   baddeley_max_dist = c.baddeley_max_dist;
   fom_alpha  = c.fom_alpha;
   zhu_weight = c.zhu_weight;
   beta_value = c.beta_value;

   return;
}

////////////////////////////////////////////////////////////////////////

double DMAPInfo::fbias() const {
   double v;

   if(oy == 0) v = bad_data_double;
   else        v = (double) fy / oy;

   return(v);
}

////////////////////////////////////////////////////////////////////////
   
void DMAPInfo::set(const SingleThresh &fthr, const SingleThresh &othr,
                   const NumArray &fdmap_na, const NumArray &odmap_na,
                   const NumArray &fthr_na,  const NumArray &othr_na) {

   // Check for mismatch
   if(fdmap_na.n() != odmap_na.n() ||
      fdmap_na.n() != fthr_na.n()  ||
      fdmap_na.n() != othr_na.n()) {
      mlog << Error << "\nDMAPInfo::set() -> "
           << "count mismatch ("
           << fdmap_na.n() << ", " << odmap_na.n() << ", "
           << fthr_na.n()  << ", " << othr_na.n() << ")\n\n";
      exit(1);
   }

   // Initialize
   clear();

   // Store the thresholds
   fthresh = fthr;
   othresh = othr;

   // Compute actual DMAP statistics here.
   int max_events;
   int non_zero_count;
   double f_distance, o_distance, zhu_common, abs_diff_distance;
   double sum_event_diff = 0.0;
   double fom_fo_sum = 0.0;
   double fom_of_sum = 0.0;
   double med_fo_sum = 0.0;
   double med_of_sum = 0.0;
   double baddeley_delta_sum = 0.0;

   non_zero_count = 0;

   mlog << Debug(4) << " DMAP.Options: baddeley_p=" << baddeley_p
        << ", baddeley_max_dist=" << baddeley_max_dist
        << ", fom_alpha=" << fom_alpha
        << ", zhu_weight=" << zhu_weight << "\n";

   for (int i=0; i<fdmap_na.n(); i++) {

      // Skip bad data
      if (is_bad_data(fdmap_na[i]) || is_bad_data(odmap_na[i]) ||
          is_bad_data(fthr_na[i])  || is_bad_data(othr_na[i])) continue;

      // Forecast event
      if (fthr_na[i] > 0) {
         fy++;
         med_of_sum += odmap_na[i];
         fom_of_sum += 1 / (1 + odmap_na[i] * odmap_na[i] * fom_alpha);
      }

      // Observation
      if (othr_na[i] > 0) {
         oy++;
         med_fo_sum += fdmap_na[i];
         fom_fo_sum += 1 / (1 + fdmap_na[i] * fdmap_na[i] * fom_alpha);
      }

      // Forecast and observation event
      if (fthr_na[i] > 0 && othr_na[i] > 0) foy++;

      sum_event_diff += (fthr_na[i] - othr_na[i]) * (fthr_na[i] - othr_na[i]);

      f_distance = (!is_bad_data(baddeley_max_dist) &&
                    fdmap_na[i] > baddeley_max_dist ?
                    baddeley_max_dist : fdmap_na[i]);
      o_distance = (!is_bad_data(baddeley_max_dist) &&
                    odmap_na[i] > baddeley_max_dist ?
                    baddeley_max_dist : odmap_na[i]);
      abs_diff_distance = abs(f_distance - o_distance);
      if (abs_diff_distance > 0.0) {
         baddeley_delta_sum += pow((double)abs_diff_distance, baddeley_p);
         non_zero_count++;
      }

      // Distance metrics
      abs_diff_distance = abs(fdmap_na[i] - odmap_na[i]);
      if (hausdorff < abs_diff_distance) hausdorff = abs_diff_distance;

      // Increment counter
      total++;
   }

   if(total == 0) {
      mlog << Error << "\nDMAPInfo::set() -> "
           << "count is zero!\n\n";
      exit(1);
   }

   max_events = max(fy, oy);

   // Mean error distance
   med_fo = (oy == 0 ? bad_data_double : med_fo_sum / oy);
   med_of = (fy == 0 ? bad_data_double : med_of_sum / fy);
   if(!is_bad_data(med_fo) && !is_bad_data(med_of)) {
      med_max  = max(med_fo, med_of);
      med_min  = min(med_fo, med_of);
      med_mean = (med_fo + med_of) / 2;
   }

   // Distance metrics
   baddeley = pow(baddeley_delta_sum/(double)total, 1.0/baddeley_p);

   // Pratt's Figure of Merit
   if(max_events > 0) {
      if(oy > 0) fom_fo = fom_fo_sum / max_events;
      if(fy > 0) fom_of = fom_of_sum / max_events;
      if(!is_bad_data(fom_fo) && !is_bad_data(fom_of)) {
         fom_max = max(fom_fo, fom_of);
         fom_min = min(fom_fo, fom_of);
         fom_mean = (fom_fo + fom_of) / 2;
      }
   }

   // Zhu Metric
   zhu_common = zhu_weight * sqrt(sum_event_diff / total);
   zhu_fo = (is_bad_data(med_fo) ? bad_data_double : zhu_common + (1-zhu_weight) * med_fo);
   zhu_of = (is_bad_data(med_of) ? bad_data_double : zhu_common + (1-zhu_weight) * med_of);
   if(!is_bad_data(zhu_fo) && !is_bad_data(zhu_of)) {
      zhu_max  = max(zhu_fo, zhu_of);
      zhu_min  = min(zhu_fo, zhu_of);
      zhu_mean = (zhu_fo + zhu_of) / 2;
   }

   // G and G-Beta
   // Reference:
   //   Gilleland, E.: Novel measures for summarizing high-resolution forecast performance,
   //     Adv. Stat. Clim. Meteorol. Oceanogr., 7, 13–34,
   //     https://doi.org/10.5194/ascmo-7-13-2021, 2021.

   // If not set by the user, default maximum distance to the number of pairs
   double max_dist = (is_bad_data(baddeley_max_dist) ?
                      (double) total : baddeley_max_dist);

   double g_med_fo = (oy == 0 ? max_dist : med_fo);
   double g_med_of = (fy == 0 ? max_dist : med_of);
   int    g_y1     = fy + oy - 2 * foy;
   double g_y2     = g_med_fo * oy + g_med_of * fy;
   double g_y      = g_y1 * g_y2;
   g               = pow(g_y, 1.0 / 3.0);

   // If not set by the user, set beta_value as the number of pairs divided by 2
   if(is_bad_data(beta_value)) beta_value = (double) total / 2.0;

   gbeta           = max(1.0 - g_y / beta_value, 0.0);

   mlog << Debug(4) << " DMAP: nf=" << fy << ", no=" << oy << ", total=" << total
        << "\tbaddeley=" << baddeley << ", hausdorff=" << hausdorff
        << "\n\tmed_fo=" << med_fo   << ", med_of="    << med_of
        << ", med_min="  << med_min  << ", med_max="   << med_max << ", med_mean="   << med_mean
        << "\n\tfom_fo=" << fom_fo   << ", fom_of="    << fom_of
        << ", fom_min="  << fom_min  << ", fom_max="   << fom_max << ", fom_mean="   << fom_mean
        << "\n\tzhu_fo=" << zhu_fo   << ", zhu_of="    << zhu_of
        << ", zhu_min="  << zhu_min  << ", zhu_max="   << zhu_max << ", zhu_mean="   << zhu_mean
        << "\n\tg="      << g        << ", gbeta="     << gbeta   << ", beta_value=" << beta_value
        << "\n";

   return;
}

////////////////////////////////////////////////////////////////////////

void DMAPInfo::set_options(const int _baddeley_p, const double _baddeley_max_dist,
                           const double _fom_alpha, const double _zhu_weight,
                           const double _beta_value) {
   baddeley_p = _baddeley_p;
   baddeley_max_dist = _baddeley_max_dist;
   fom_alpha = _fom_alpha;
   zhu_weight = _zhu_weight;
   beta_value = _beta_value;
}

////////////////////////////////////////////////////////////////////////
//
// Begin code for misc functions
//
////////////////////////////////////////////////////////////////////////

int parse_message_type(const char *msg_typ_str, char **&msg_typ_arr) {
   char tmp_str[max_str_len];
   char *c = (char *) 0;
   int n, i;

   // Compute the number of tokens in the string based on " "
   n = num_tokens(msg_typ_str, " ");

   // Check for no tokens in string
   if(n == 0) return(0);

   // Allocate space for the list of tokens
   msg_typ_arr = new char * [n];

   // Initialize the temp string for use in tokenizing
   strcpy(tmp_str, msg_typ_str);

   // Tokenize the string and store the double values
   c = strtok(tmp_str, " ");
   msg_typ_arr[0] = new char [strlen(c)+1];
   strcpy(msg_typ_arr[0], c);

   // Parse remaining tokens
   for(i=1; i<n; i++) {
      c = strtok(0, " ");
      msg_typ_arr[i] = new char [strlen(c)+1];
      strcpy(msg_typ_arr[i], c);
   }

   return(n);
}

////////////////////////////////////////////////////////////////////////

int parse_dbl_list(const char *dbl_str, double *&dbl_arr) {
   char tmp_str[max_str_len];
   char *c = (char *) 0;
   int n, i;

   // Compute the number of tokens in the string based on " "
   n = num_tokens(dbl_str, " ");

   // Check for no tokens in string
   if(n == 0) return(0);

   // Allocate space for the list of tokens
   dbl_arr = new double [n];

   // Initialize the temp string for use in tokenizing
   strcpy(tmp_str, dbl_str);

   // Tokenize the string and store the double values
   c = strtok(tmp_str, " ");
   dbl_arr[0] = atof(c);

   // Parse remaining tokens
   for(i=1; i<n; i++) dbl_arr[i] = atof(strtok(0, " "));

   return(n);
}

////////////////////////////////////////////////////////////////////////

int parse_int_list(const char *int_str, int *&int_arr) {
   char tmp_str[max_str_len];
   char *c = (char *) 0;
   int n, i;

   // Compute the number of tokens in the string based on " "
   n = num_tokens(int_str, " ");

   // Check for no tokens in string
   if(n == 0) return(0);

   // Allocate space for the list of tokens
   int_arr = new int [n];

   // Initialize the temp string for use in tokenizing
   strcpy(tmp_str, int_str);

   // Tokenize the string and store the integer values
   c = strtok(tmp_str, " ");
   int_arr[0] = nint(atof(c));

   // Parse remaining tokens
   for(i=1; i<n; i++) int_arr[i] = nint(atof(strtok(0, " ")));

   return(n);
}

////////////////////////////////////////////////////////////////////////

int max_int(const int *v_int, int n) {
   int i, v_max;

   if(n <= 0) return(0);

   v_max = v_int[0];
   for(i=1; i<n; i++) if(v_int[i] > v_max) v_max = v_int[i];

   return(v_max);
}

////////////////////////////////////////////////////////////////////////

int min_int(const int *v_int, int n) {
   int i, v_min;

   if(n <= 0) return(0);

   v_min = v_int[0];
   for(i=1; i<n; i++) if(v_int[i] < v_min) v_min = v_int[i];

   return(v_min);
}

////////////////////////////////////////////////////////////////////////
//
// Compute variance from sums of squares
//
////////////////////////////////////////////////////////////////////////

double compute_variance(double sum, double sum_sq, int n) {
   double v;

   if(n <= 1) {
      v = bad_data_double;
   }
   else {
      v = (sum_sq - sum*sum/(double) n)/((double) (n - 1));

           if(is_eq(v, 0.0)) v = 0.0;
      else if(v < 0)         v = bad_data_double;
   }

   return(v);
}

////////////////////////////////////////////////////////////////////////
//
// Compute standard deviation from sums of squares
//
////////////////////////////////////////////////////////////////////////

double compute_stdev(double sum, double sum_sq, int n) {
   return(square_root(compute_variance(sum, sum_sq, n)));
}

////////////////////////////////////////////////////////////////////////
//
// Compute correlation from sums of squares
//
////////////////////////////////////////////////////////////////////////

double compute_corr(double f, double o, double ff, double oo, double fo,
                    int n) {
   double v, c;

   v = (n*ff - f*f)*(n*oo - o*o);

   // Check for divide by zero
   if(v < 0 || is_eq(v, 0.0)) {
      c = bad_data_double;
   }
   else {
      c = ((n*fo) - (f*o))/sqrt(v);

      // Check the computed range
           if(c >  1) c =  1.0;
      else if(c < -1) c = -1.0;
   }

   return(c);
}

////////////////////////////////////////////////////////////////////////
//
// Compute the uncentered anomaly correlation from sums of squares
// without adjusting for the mean error.
//
////////////////////////////////////////////////////////////////////////

double compute_anom_corr_uncntr(double ffa, double ooa, double foa) {
   double v, c;

   v = ffa*ooa;

   // Check for square root of negative number
   if(v < 0) {
      c = bad_data_double;
   }
   else {
      c = foa/sqrt(v);

      // Check the computed range
           if(c >  1) c =  1.0;
      else if(c < -1) c = -1.0;
   }

   return(c);
}

////////////////////////////////////////////////////////////////////////

double compute_afss(double f_rate, double o_rate) {
   double num, den, afss;

   //
   // Compute Asymptotic Fractions Skill Score
   //
   num = 2.0*f_rate*o_rate;
   den = f_rate*f_rate + o_rate*o_rate;

   if(is_bad_data(f_rate) || is_bad_data(o_rate) || is_eq(den, 0.0)) {
      afss = bad_data_double;
   }
   else {
      afss = num/den;
   }

   return(afss);
}

////////////////////////////////////////////////////////////////////////

double compute_ufss(double o_rate) {
   double ufss;

   //
   // Compute Uniform Fractions Skill Score
   //
   if(is_bad_data(o_rate)) ufss = bad_data_double;
   else                    ufss = 0.5 + o_rate/2.0;

   return(ufss);
}

///////////////////////////////////////////////////////////////////////////////
//
// Convert a field of data into its corresponding field of ranks of that data.
// Return the number of valid data points that were ranked and keep track
// of the number of rank ties.
//
///////////////////////////////////////////////////////////////////////////////

int compute_rank(const DataPlane &dp, DataPlane &dp_rank, double *data_rank, int &ties) {
   int x, y, n, i;
   double *data = (double *) 0, v;
   int *data_loc = (int *) 0;

   // Arrays to store the raw data values to be ranked, their locations,
   // and their computed ranks.  The ranks are stored as doubles since
   // they can be set to 0.5 in the case of ties.
   data      = new double [dp.nx()*dp.ny()];
   data_loc  = new int    [dp.nx()*dp.ny()];

   // Search the input field for valid data and keep track of its location
   n = 0;
   for(x=0; x<dp.nx(); x++) {
      for(y=0; y<dp.ny(); y++) {

         v = dp.get(x, y);
         if(!is_bad_data(v)) {
            data[n] = v;
            data_loc[n] = dp.two_to_one(x, y);
            n++;
         }
      }
   }

   // Compute the rank of the data and store the ranks in the data_rank array
   // Keep track of the number of ties in the ranks
   ties = do_rank(data, data_rank, n);

   // Set up the dp_rank object
   dp_rank.set_size(dp.nx(), dp.ny());

   // Assign the ranks to the dp_rank field
   for(i=0; i<n; i++) {
      dp_rank.one_to_two(data_loc[i], x, y);
      dp_rank.set(data_rank[i], x, y);
   }

   // Deallocate memory
   if(data)      { delete [] data;      data = (double *) 0;  }
   if(data_loc)  { delete [] data_loc;  data_loc = (int *) 0; }

   return(n);
}

////////////////////////////////////////////////////////////////////////
