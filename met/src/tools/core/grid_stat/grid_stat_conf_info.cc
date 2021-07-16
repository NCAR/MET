// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

using namespace std;

#include <dirent.h>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <cmath>

#include "grid_stat_conf_info.h"
#include "configobjecttype_to_string.h"

#include "vx_data2d_factory.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////
//
//  Code for class GridStatConfInfo
//
////////////////////////////////////////////////////////////////////////

GridStatConfInfo::GridStatConfInfo() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

GridStatConfInfo::~GridStatConfInfo() {

   clear();
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::init_from_scratch() {

   // Initialize pointers
   vx_opt = (GridStatVxOpt *) 0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::clear() {
   int i;

   // Initialize values
   model.clear();
   obtype.clear();
   mask_map.clear();
   grid_weight_flag = GridWeightType_None;
   tmp_dir.clear();
   output_prefix.clear();
   version.clear();

   for(i=0; i<n_txt; i++) output_flag[i] = STATOutputType_None;

   nc_info.clear();

   output_ascii_flag = false;
   output_nc_flag = false;

   // Deallocate memory
   if(vx_opt) { delete [] vx_opt; vx_opt = (GridStatVxOpt *) 0; }

   // Reset counts
   n_vx = 0;

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::read_config(const char *default_file_name,
                                   const char *user_file_name) {

   // Read the config file constants
   conf.read(replace_path(config_const_filename).c_str());

   // Read the default config file
   conf.read(default_file_name);

   // Read the user-specified config file
   conf.read(user_file_name);

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::process_config(GrdFileType ftype,
                                      GrdFileType otype) {
   int i, j, n_fvx, n_ovx;
   Dictionary *fdict = (Dictionary *) 0;
   Dictionary *odict = (Dictionary *) 0;
   Dictionary i_fdict, i_odict;

   // Dump the contents of the config file
   if(mlog.verbosity_level() >= 5) conf.dump(cout);

   // Initialize
   clear();

   // Conf: version
   version = parse_conf_version(&conf);

   // Conf: model
   model = parse_conf_string(&conf, conf_key_model);

   // Conf: obtype
   obtype = parse_conf_string(&conf, conf_key_obtype);

   // Conf: grid_weight_flag
   grid_weight_flag = parse_conf_grid_weight_flag(&conf);

   // Conf: tmp_dir
   tmp_dir = parse_conf_tmp_dir(&conf);

   // Conf: output_prefix
   output_prefix = conf.lookup_string(conf_key_output_prefix);

   // Conf: fcst.field and obs.field
   fdict = conf.lookup_array(conf_key_fcst_field);
   odict = conf.lookup_array(conf_key_obs_field);

   // Determine the number of fields (name/level) to be verified
   n_fvx = parse_conf_n_vx(fdict);
   n_ovx = parse_conf_n_vx(odict);

   // Check for a valid number of verification tasks
   if(n_fvx == 0 || n_fvx != n_ovx) {
      mlog << Error << "\nGridStatConfInfo::process_config() -> "
           << "The number of verification tasks in \""
           << conf_key_obs_field << "\" (" << n_ovx
           << ") must be non-zero and match the number in \""
           << conf_key_fcst_field << "\" (" << n_fvx << ").\n\n";
      exit(1);
   }

   // Allocate memory for the verification task options
   n_vx   = n_fvx;
   vx_opt = new GridStatVxOpt [n_vx];

   // Check for consistent number of climatology fields
   check_climo_n_vx(&conf, n_vx);

   // Parse settings for each verification task
   for(i=0; i<n_vx; i++) {

      // Get the current dictionaries
      i_fdict = parse_conf_i_vx_dict(fdict, i);
      i_odict = parse_conf_i_vx_dict(odict, i);

      // Process the options for this verification task
      vx_opt[i].process_config(ftype, i_fdict, otype, i_odict);
   }

   // Summarize output flags across all verification tasks
   process_flags();

   // If VL1L2 or VAL1L2 is requested, set the uv_index
   if(output_flag[i_vl1l2]  != STATOutputType_None ||
      output_flag[i_val1l2] != STATOutputType_None) {

      for(i=0; i<n_vx; i++) {

         // Process u-wind
         if(vx_opt[i].fcst_info->is_u_wind() &&
            vx_opt[i].obs_info->is_u_wind()) {

            // Search for corresponding v-wind
            for(j=0; j<n_vx; j++) {
               if(vx_opt[j].fcst_info->is_v_wind()      &&
                  vx_opt[j].obs_info->is_v_wind()       &&
                  vx_opt[i].fcst_info->req_level_name() ==
                  vx_opt[j].fcst_info->req_level_name() &&
                  vx_opt[i].obs_info->req_level_name()  ==
                  vx_opt[j].obs_info->req_level_name()  &&
                  vx_opt[i].is_uv_match(vx_opt[j])) {

                  vx_opt[i].fcst_info->set_uv_index(j);
                  vx_opt[i].obs_info->set_uv_index(j);
               }
            }
         }
         // Process v-wind
         else if(vx_opt[i].fcst_info->is_v_wind() &&
                 vx_opt[i].obs_info->is_v_wind()) {

            // Search for corresponding u-wind
            for(j=0; j<n_vx; j++) {
               if(vx_opt[j].fcst_info->is_u_wind()      &&
                  vx_opt[j].obs_info->is_u_wind()       &&
                  vx_opt[i].fcst_info->req_level_name() ==
                  vx_opt[j].fcst_info->req_level_name() &&
                  vx_opt[i].obs_info->req_level_name()  ==
                  vx_opt[j].obs_info->req_level_name()  &&
                  vx_opt[i].is_uv_match(vx_opt[j])) {

                  vx_opt[i].fcst_info->set_uv_index(j);
                  vx_opt[i].obs_info->set_uv_index(j);
               }
            }
         }
      } // end for i
   } // end if

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::process_flags() {
   int i, j;

   // Initialize
   for(i=0; i<n_txt; i++) output_flag[i] = STATOutputType_None;
   nc_info.set_all_false();
   output_ascii_flag = false;

   // Loop over the verification tasks
   for(i=0; i<n_vx; i++) {

      // Summary of output_flag settings
      for(j=0; j<n_txt; j++) {
         if(vx_opt[i].output_flag[j] == STATOutputType_Both) {
            output_flag[j] = STATOutputType_Both;
            output_ascii_flag = true;
         }
         else if(vx_opt[i].output_flag[j] == STATOutputType_Stat &&
                           output_flag[j] == STATOutputType_None) {
            output_flag[j] = STATOutputType_Stat;
            output_ascii_flag = true;
         }
      }

      // Summary of nc_info flag settings
      if(vx_opt[i].nc_info.do_latlon)       nc_info.do_latlon       = true;
      if(vx_opt[i].nc_info.do_raw)          nc_info.do_raw          = true;
      if(vx_opt[i].nc_info.do_diff)         nc_info.do_diff         = true;
      if(vx_opt[i].nc_info.do_climo)        nc_info.do_climo        = true;
      if(vx_opt[i].nc_info.do_climo_cdp)    nc_info.do_climo_cdp    = true;
      if(vx_opt[i].nc_info.do_weight)       nc_info.do_weight       = true;
      if(vx_opt[i].nc_info.do_nbrhd)        nc_info.do_nbrhd        = true;
      if(vx_opt[i].nc_info.do_fourier)      nc_info.do_fourier      = true;
      if(vx_opt[i].nc_info.do_gradient)     nc_info.do_gradient     = true;
      if(vx_opt[i].nc_info.do_distance_map) nc_info.do_distance_map = true;
      if(vx_opt[i].nc_info.do_apply_mask)   nc_info.do_apply_mask   = true;
   }

   // Check output_ascii_flag
   if(!output_ascii_flag) {
       mlog << Debug(3)
            <<"\nNo STAT output types requested, "
            << "proceeding with ASCII output flag set to false.\n";
   }

   // Set output_nc_flag
   output_nc_flag = !nc_info.all_false();

   // Check for at least one output data type
   if(!output_ascii_flag && !output_nc_flag) {
      mlog << Error << "\nGridStatConfInfo::process_flags() -> "
           << "At least one output STAT or NetCDF type must be "
           << " requested in \"" << conf_key_output_flag << "\" or \""
           << conf_key_nc_pairs_flag << "\".\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatConfInfo::process_masks(const Grid &grid) {
   int i, j;
   MaskPlane mp;
   ConcatString name;

   mlog << Debug(2)
        << "Processing masking regions.\n";

   // Mapping of grid definition strings to mask names
   map<ConcatString,ConcatString> grid_map;
   map<ConcatString,ConcatString> poly_map;

   // Initiailize
   mask_map.clear();

   // Process the masks for each vx task
   for(i=0; i<n_vx; i++) {

      // Initialize
      vx_opt[i].mask_name.clear();

      // Parse the masking grids
      for(j=0; j<vx_opt[i].mask_grid.n_elements(); j++) {

         // Process new grid masks
         if(grid_map.count(vx_opt[i].mask_grid[j]) == 0) {
            mlog << Debug(3)
                 << "Processing grid mask: "
                 << vx_opt[i].mask_grid[j] << "\n";
            parse_grid_mask(vx_opt[i].mask_grid[j], grid, mp, name);
            grid_map[vx_opt[i].mask_grid[j]] = name;
            mask_map[name] = mp;
         }

         // Store the name for the current grid mask
         vx_opt[i].mask_name.add(grid_map[vx_opt[i].mask_grid[j]]);

      } // end for j

      // Parse the masking polylines
      for(j=0; j<vx_opt[i].mask_poly.n_elements(); j++) {

         // Process new poly mask
         if(poly_map.count(vx_opt[i].mask_poly[j]) == 0) {
            mlog << Debug(3)
                 << "Processing poly mask: "
                 << vx_opt[i].mask_poly[j] << "\n";
            parse_poly_mask(vx_opt[i].mask_poly[j], grid, mp, name);
            poly_map[vx_opt[i].mask_poly[j]] = name;
            mask_map[name] = mp;
         }

         // Store the name for the current poly mask
         vx_opt[i].mask_name.add(poly_map[vx_opt[i].mask_poly[j]]);

      } // end for j

      // Check that at least one verification masking region is provided
      if(vx_opt[i].mask_name.n_elements() == 0) {
         mlog << Error << "\nGridStatConfInfo::process_masks() -> "
              << "At least one grid or polyline verification masking "
              << "region must be provided for verification task number "
              << i+1 << ".\n\n";
         exit(1);
      }

      // Check for unique mask names
      check_mask_names(vx_opt[i].mask_name);

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::n_txt_row(int i_txt_row) const {
   int i, n;

   // Loop over the tasks and sum the line counts for this line type
   for(i=0, n=0; i<n_vx; i++) n += vx_opt[i].n_txt_row(i_txt_row);

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::n_stat_row() const {
   int i, n;

   // Loop over the line types and sum the line counts
   for(i=0, n=0; i<n_txt; i++) n += n_txt_row(i);

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_cat_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_cat_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_cnt_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_cnt_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_wind_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_wind_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_fprob_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_fprob_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_oprob_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_oprob_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_eclv_points() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_eclv_points());

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatConfInfo::get_max_n_cov_thresh() const {
   int i, n;

   for(i=0,n=0; i<n_vx; i++) n = max(n, vx_opt[i].get_n_cov_thresh());

   return(n);
}

////////////////////////////////////////////////////////////////////////
//
//  Code for class GridStatVxOpt
//
////////////////////////////////////////////////////////////////////////

GridStatVxOpt::GridStatVxOpt() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

GridStatVxOpt::~GridStatVxOpt() {

   clear();
}

////////////////////////////////////////////////////////////////////////

void GridStatVxOpt::init_from_scratch() {

   // Initialize pointers
   fcst_info   = (VarInfo *)    0;
   obs_info    = (VarInfo *)    0;

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatVxOpt::clear() {
   int i;

   // Initialize values

   desc.clear();
   var_name.clear();
   var_suffix.clear();

   mpr_sa.clear();
   mpr_ta.clear();

   fcat_ta.clear();
   ocat_ta.clear();

   fcnt_ta.clear();
   ocnt_ta.clear();
   cnt_logic = SetLogic_None;

   fwind_ta.clear();
   owind_ta.clear();
   wind_logic = SetLogic_None;

   mask_grid.clear();
   mask_poly.clear();
   mask_name.clear();

   eclv_points.clear();

   cdf_info.clear();

   ci_alpha.clear();

   boot_info.clear();
   interp_info.clear();
   nbrhd_info.clear();

   wave_1d_beg.clear();
   wave_1d_end.clear();

   grad_dx.clear();
   grad_dy.clear();

   baddeley_p = bad_data_int;
   baddeley_max_dist = bad_data_double;
   fom_alpha = bad_data_double;
   zhu_weight = bad_data_double;
   beta_value = bad_data_double;

   hss_ec_value = bad_data_double;
   rank_corr_flag = false;

   for(i=0; i<n_txt; i++) output_flag[i] = STATOutputType_None;

   nc_info.clear();

   // Deallocate memory
   if(fcst_info) { delete fcst_info; fcst_info = (VarInfo *) 0; }
   if(obs_info)  { delete obs_info;  obs_info  = (VarInfo *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatVxOpt::process_config(
        GrdFileType ftype, Dictionary &fdict,
        GrdFileType otype, Dictionary &odict) {
   int i, n;
   VarInfoFactory info_factory;
   map<STATLineType,STATOutputType>output_map;
   InterpMthd mthd;

   // Initialize
   clear();

   // Allocate new VarInfo objects
   fcst_info = info_factory.new_var_info(ftype);
   obs_info  = info_factory.new_var_info(otype);

   // Set the VarInfo objects
   fcst_info->set_dict(fdict);
   obs_info->set_dict(odict);

   // Dump the contents of the current VarInfo
   if(mlog.verbosity_level() >= 5) {
      mlog << Debug(5)
           << "Parsed forecast field:\n";
      fcst_info->dump(cout);
      mlog << Debug(5)
           << "Parsed observation field:\n";
      obs_info->dump(cout);
   }

   // No support for wind direction
   if(fcst_info->is_wind_direction() || obs_info->is_wind_direction()) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "wind direction may not be verified using grid_stat.\n\n";
      exit(1);
   }

   // Check that the observation field does not contain probabilities
   if(obs_info->is_prob()) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "the observation field cannot contain probabilities.\n\n";
      exit(1);
   }

   // Conf: desc
   desc = parse_conf_string(&odict, conf_key_desc);

   // Conf: nc_pairs_var_name
   var_name = parse_conf_string(&odict, conf_key_nc_pairs_var_name, false);

   // Check to see if the deprecated nc_pairs_var_str has been used

   // Conf: nc_pairs_var_str (deprecated)
   var_suffix = odict.lookup_string(conf_key_nc_pairs_var_str, false, false);

   // If found, print a warning.
   if(odict.last_lookup_status()) {
      mlog << Warning << "\nGridStatVxOpt::process_config() -> \""
           << conf_key_nc_pairs_var_str << "\" (" << var_suffix
           << ") is deprecated and replaced by \""
           << conf_key_nc_pairs_var_suffix << "\".\n\n";
   }
   // Otherwise, parse the new option
   else {
      // Conf: nc_pairs_var_suffix
      var_suffix = parse_conf_string(&odict, conf_key_nc_pairs_var_suffix, false);
   }

   // Conf: output_flag
   output_map = parse_conf_output_flag(&odict, txt_file_type, n_txt);

   // Populate the output_flag array with map values
   for(i=0; i<n_txt; i++) output_flag[i] = output_map[txt_file_type[i]];

   // Conf: mpr_column and mpr_thresh
   mpr_sa = odict.lookup_string_array(conf_key_mpr_column);
   mpr_ta = odict.lookup_thresh_array(conf_key_mpr_thresh);

   // Conf: cat_thresh
   fcat_ta = fdict.lookup_thresh_array(conf_key_cat_thresh);
   ocat_ta = odict.lookup_thresh_array(conf_key_cat_thresh);

   // Conf: cnt_thresh
   fcnt_ta = process_perc_thresh_bins(
                fdict.lookup_thresh_array(conf_key_cnt_thresh));
   ocnt_ta = process_perc_thresh_bins(
                odict.lookup_thresh_array(conf_key_cnt_thresh));

   // Conf: cnt_logic
   cnt_logic = check_setlogic(
      int_to_setlogic(fdict.lookup_int(conf_key_cnt_logic)),
      int_to_setlogic(odict.lookup_int(conf_key_cnt_logic)));

   // Conf: wind_thresh
   fwind_ta = process_perc_thresh_bins(
                 fdict.lookup_thresh_array(conf_key_wind_thresh));
   owind_ta = process_perc_thresh_bins(
                 odict.lookup_thresh_array(conf_key_wind_thresh));

   // Conf: wind_logic
   wind_logic = check_setlogic(
      int_to_setlogic(fdict.lookup_int(conf_key_wind_logic)),
      int_to_setlogic(odict.lookup_int(conf_key_wind_logic)));

   // Dump the contents of the current thresholds
   if(mlog.verbosity_level() >= 5) {
      mlog << Debug(5)
           << "Parsed thresholds:\n"
           << "Matched pair filter columns:     " << write_css(mpr_sa) << "\n"
           << "Matched pair filter thresholds:  " << mpr_ta.get_str() << "\n"
           << "Forecast categorical thresholds: " << fcat_ta.get_str() << "\n"
           << "Observed categorical thresholds: " << ocat_ta.get_str() << "\n"
           << "Forecast continuous thresholds:  " << fcnt_ta.get_str() << "\n"
           << "Observed continuous thresholds:  " << ocnt_ta.get_str() << "\n"
           << "Continuous threshold logic:      " << setlogic_to_string(cnt_logic) << "\n"
           << "Forecast wind speed thresholds:  " << fwind_ta.get_str() << "\n"
           << "Observed wind speed thresholds:  " << owind_ta.get_str() << "\n"
           << "Wind speed threshold logic:      " << setlogic_to_string(wind_logic) << "\n";
   }

   // Verifying a probability field
   if(fcst_info->is_prob()) {
      fcat_ta = string_to_prob_thresh(fcat_ta.get_str().c_str());
   }

   // Check for equal threshold length for non-probability fields
   if(!fcst_info->is_prob() &&
      fcat_ta.n_elements() != ocat_ta.n_elements()) {

      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The number of thresholds for each field in \"fcst."
           << conf_key_cat_thresh
           << "\" must match the number of thresholds for each "
           << "field in \"obs." << conf_key_cat_thresh << "\".\n\n";
      exit(1);
   }

   // Add default continuous thresholds until the counts match
   n = max(fcnt_ta.n_elements(), ocnt_ta.n_elements());
   while(fcnt_ta.n_elements() < n) fcnt_ta.add(na_str);
   while(ocnt_ta.n_elements() < n) ocnt_ta.add(na_str);

   // Add default wind speed thresholds until the counts match
   n = max(fwind_ta.n_elements(), owind_ta.n_elements());
   while(fwind_ta.n_elements() < n) fwind_ta.add(na_str);
   while(owind_ta.n_elements() < n) owind_ta.add(na_str);

   // Verifying with multi-category contingency tables
   if(!fcst_info->is_prob() &&
      (output_flag[i_mctc] != STATOutputType_None ||
       output_flag[i_mcts] != STATOutputType_None)) {
      check_mctc_thresh(fcat_ta);
      check_mctc_thresh(ocat_ta);
   }

   // Conf: mask_grid
   mask_grid = odict.lookup_string_array(conf_key_mask_grid);

   // Conf: mask_poly
   mask_poly = odict.lookup_string_array(conf_key_mask_poly);

   // Conf: eclv_points
   eclv_points = parse_conf_eclv_points(&odict);

   // Conf: climo_cdf
   cdf_info = parse_conf_climo_cdf(&odict);

   // Conf: ci_alpha
   ci_alpha = parse_conf_ci_alpha(&odict);

   // Conf: boot
   boot_info = parse_conf_boot(&odict);

   // Conf: interp
   interp_info = parse_conf_interp(&odict, conf_key_interp);

   // Loop through interpolation options
   for(i=0; i<interp_info.n_interp; i++) {

      mthd = string_to_interpmthd(interp_info.method[i].c_str());

      // Check for unsupported interpolation methods
      if(mthd == InterpMthd_DW_Mean ||
         mthd == InterpMthd_LS_Fit  ||
         mthd == InterpMthd_Bilin) {
         mlog << Error << "\nGridStatVxOpt::process_config() -> "
              << "Interpolation methods DW_MEAN, LS_FIT, and BILIN are "
              << "not supported in Grid-Stat.\n\n";
         exit(1);
      }

      // Check for valid interpolation width
      if(interp_info.width[i] < 1 || interp_info.width[i]%2 == 0) {
         mlog << Error << "\nGridStatVxOpt::process_config() -> "
              << "The interpolation width must be set to odd values "
              << "greater than or equal to 1 ("
              << interp_info.width[i] << ").\n\n";
         exit(1);
      }
   } // end for i

   // Conf: nbrhd
   nbrhd_info = parse_conf_nbrhd(&odict, conf_key_nbrhd);

   // Conf: fourier
   Dictionary * d = odict.lookup_dictionary(conf_key_fourier);
   wave_1d_beg = d->lookup_int_array(conf_key_wave_1d_beg);
   wave_1d_end = d->lookup_int_array(conf_key_wave_1d_end);

   // Check for the same length
   if(wave_1d_beg.n_elements() != wave_1d_end.n_elements()) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The fourier wave_1d_beg and wave_1d_end arrays must have "
           << "the same length (" << wave_1d_beg.n_elements() << " != "
           << wave_1d_end.n_elements() << ").\n\n";
      exit(1);
   }

   // Conf: gradient
   d = odict.lookup_dictionary(conf_key_gradient);
   grad_dx = d->lookup_int_array(conf_key_dx);
   grad_dy = d->lookup_int_array(conf_key_dy);

   // Check for the same length
   if(grad_dx.n_elements() != grad_dy.n_elements()) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The gradient dx and dy arrays must have "
           << "the same length (" << grad_dx.n_elements() << " != "
           << grad_dy.n_elements() << ").\n\n";
      exit(1);
   }

   // Conf: distance_map
   d = odict.lookup_dictionary(conf_key_distance_map);

   baddeley_p = d->lookup_int(conf_key_baddeley_p);
   if(baddeley_p < 1) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The \"" << conf_key_baddeley_p << "\" option ("
           << baddeley_p << ") must be set >= 1.\n\n";
      exit(1);
   }

   baddeley_max_dist = d->lookup_double(conf_key_baddeley_max_dist, false);
   if(!is_bad_data(baddeley_max_dist) && baddeley_max_dist < 0) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The \"" << conf_key_baddeley_max_dist << "\" option ("
           << baddeley_max_dist << ") must be set >= 0.\n\n";
      exit(1);
   }

   fom_alpha = d->lookup_double(conf_key_fom_alpha);
   if(fom_alpha <=0 || fom_alpha > 1) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The \"" << conf_key_fom_alpha << "\" option ("
           << fom_alpha << ") must be set > 0 and <= 1.\n\n";
      exit(1);
   }

   zhu_weight = d->lookup_double(conf_key_zhu_weight);
   if(zhu_weight <=0 || zhu_weight > 1) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The \"" << conf_key_zhu_weight << "\" option ("
           << zhu_weight << ") must be set > 0 and <= 1.\n\n";
      exit(1);
   }

   beta_value = d->lookup_double(conf_key_beta_value);
   if(!is_bad_data(beta_value) && beta_value <=0) {
      mlog << Error << "\nGridStatVxOpt::process_config() -> "
           << "The \"" << conf_key_beta_value << "\" option ("
           << beta_value << ") must be set > 0.\n\n";
      exit(1);
   }

   // Conf: hss_ec_value
   hss_ec_value = odict.lookup_double(conf_key_hss_ec_value);

   // Conf: rank_corr_flag
   rank_corr_flag = odict.lookup_bool(conf_key_rank_corr_flag);

   // Conf: nc_pairs_flag
   parse_nc_info(odict);

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatVxOpt::parse_nc_info(Dictionary &odict) {
   const DictionaryEntry * e = (const DictionaryEntry *) 0;

   e = odict.lookup(conf_key_nc_pairs_flag);

   if(!e) {
      mlog << Error
           << "\nGridStatVxOpt::parse_nc_info() -> "
           << "lookup failed for key \"" << conf_key_nc_pairs_flag
           << "\"\n\n";
      exit(1);
   }

   const ConfigObjectType type = e->type();

   if(type == BooleanType) {
      bool value = e->b_value();

      if(!value) nc_info.set_all_false();

      return;
   }

   // It should be a dictionary
   if(type != DictionaryType) {
      mlog << Error
           << "\nGridStatVxOpt::parse_nc_info() -> "
           << "bad type (" << configobjecttype_to_string(type)
           << ") for key \"" << conf_key_nc_pairs_flag << "\"\n\n";
      exit(1);
   }

   // Parse the various entries
   Dictionary * d = e->dict_value();

   nc_info.do_latlon       = d->lookup_bool(conf_key_latlon_flag);
   nc_info.do_raw          = d->lookup_bool(conf_key_raw_flag);
   nc_info.do_diff         = d->lookup_bool(conf_key_diff_flag);
   nc_info.do_climo        = d->lookup_bool(conf_key_climo_flag);
   nc_info.do_climo_cdp    = d->lookup_bool(conf_key_climo_cdp_flag);
   nc_info.do_weight       = d->lookup_bool(conf_key_weight);
   nc_info.do_nbrhd        = d->lookup_bool(conf_key_nbrhd);
   nc_info.do_fourier      = d->lookup_bool(conf_key_fourier);
   nc_info.do_gradient     = d->lookup_bool(conf_key_gradient);
   nc_info.do_distance_map = d->lookup_bool(conf_key_distance_map);
   nc_info.do_apply_mask   = d->lookup_bool(conf_key_apply_mask_flag);

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Check the settings that would impact the number of matched pairs
// when searching for U/V matches.
//
////////////////////////////////////////////////////////////////////////

bool GridStatVxOpt::is_uv_match(const GridStatVxOpt &v) const {
   bool match = true;

   //
   // The following do not impact matched pairs:
   //    desc, var_name, var_suffix,
   //    mpr_sa, mpr_ta,
   //    fcat_ta, ocat_ta,
   //    fcnt_ta, ocnt_ta, cnt_logic,
   //    fwind_ta, owind_ta, wind_logic,
   //    eclv_points, cdf_info, ci_alpha
   //    boot_info, nbrhd_info,
   //    wave_1d_beg, wave_1d_end, grad_dx, grad_dy,
   //    hss_ec_value, rank_corr_flag, output_flag, nc_info
   //

   if(!(mask_grid   == v.mask_grid  ) ||
      !(mask_poly   == v.mask_poly  ) ||
      !(mask_name   == v.mask_name  ) ||
      !(interp_info == v.interp_info)) {
      match = false;
   }

   return(match);
}

////////////////////////////////////////////////////////////////////////

void GridStatVxOpt::set_perc_thresh(const PairDataPoint &pd) {

   //
   // Compute percentiles for forecast and observation thresholds,
   // but not for wind speed or climatology CDF thresholds.
   //
   if(!fcat_ta.need_perc() && !ocat_ta.need_perc() &&
      !fcnt_ta.need_perc() && !ocnt_ta.need_perc()) return;

   //
   // Sort the input arrays
   //
   NumArray fsort = pd.f_na;
   NumArray osort = pd.o_na;
   NumArray csort = pd.cmn_na;
   fsort.sort_array();
   osort.sort_array();
   csort.sort_array();

   //
   // Compute percentiles
   //
   fcat_ta.set_perc(&fsort, &osort, &csort, &fcat_ta, &ocat_ta);
   ocat_ta.set_perc(&fsort, &osort, &csort, &fcat_ta, &ocat_ta);
   fcnt_ta.set_perc(&fsort, &osort, &csort, &fcnt_ta, &ocnt_ta);
   ocnt_ta.set_perc(&fsort, &osort, &csort, &fcnt_ta, &ocnt_ta);

   return;
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::n_txt_row(int i_txt_row) const {
   int n = 0;
   int n_bin;

   // Range check
   if(i_txt_row < 0 || i_txt_row >= n_txt) {
      mlog << Error << "\nGridStatVxOpt::n_txt_row(int) -> "
           << "range check error for " << i_txt_row << "\n\n";
      exit(1);
   }

   // Check if this output line type is requested
   if(output_flag[i_txt_row] == STATOutputType_None) return(0);

   bool prob_flag = fcst_info->is_prob();
   bool vect_flag = (fcst_info->is_u_wind() && obs_info->is_u_wind());

   // Determine row multiplier for climatology bins
   if(cdf_info.write_bins) {
      n_bin = get_n_cdf_bin();
      if(n_bin > 1) n_bin++;
   }
   else {
      n_bin = 1;
   }

   // Switch on the index of the line type
   switch(i_txt_row) {

      case(i_fho):
      case(i_ctc):
         // Number of FHO or CTC lines =
         //    Masks * Smoothing Methods * Thresholds
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_cat_thresh());
         break;

      case(i_cts):
         // Number of CTS lines =
         //    Masks * Smoothing Methods * Thresholds * Alphas
         n = (prob_flag ? 0:
              get_n_mask() * get_n_interp() * get_n_cat_thresh() *
              get_n_ci_alpha());
         break;

      case(i_mctc):
         // Number of MCTC lines =
         //    Masks * Smoothing Methods
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_interp());
         break;

      case(i_mcts):
         // Number of MCTS lines =
         //    Masks * Smoothing Methods * Alphas
         n = (prob_flag ? 0:
              get_n_mask() * get_n_interp() * get_n_ci_alpha());
         break;

      case(i_cnt):
         // Number of CNT lines =
         //    Masks * (Smoothing Methods + Fourier Waves) *
         //    Thresholds * Climo Bins * Alphas
         n = (prob_flag ? 0 :
              get_n_mask() * (get_n_interp() + get_n_wave_1d()) *
              get_n_cnt_thresh() * n_bin * get_n_ci_alpha());
         break;

      case(i_sl1l2):
      case(i_sal1l2):
         // Number of SL1L2 or SAL1L2 lines =
         //    Masks * (Smoothing Methods + Fourier Waves) *
         //    Thresholds * Climo Bins
         n = (prob_flag ? 0 :
              get_n_mask() * (get_n_interp() + get_n_wave_1d()) *
              get_n_cnt_thresh() * n_bin);
         break;

      case(i_vl1l2):
      case(i_val1l2):
      case(i_vcnt):
         // Number of VL1L2 or VAL1L2 lines =
         //    Masks * (Smoothing Methods + Fourier Waves) * Thresholds
         n = (!vect_flag ? 0 :
              get_n_mask() * (get_n_interp() + get_n_wave_1d()) *
              get_n_wind_thresh());
         break;

      case(i_nbrctc):
         // Number of NBRCTC lines =
         //    Masks * Thresholds * Neighborhoods * Frac Thresholds
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_cat_thresh() * get_n_nbrhd_wdth() *
              get_n_cov_thresh());
         break;

      case(i_nbrcts):
         // Number of NBRCTS lines =
         //    Masks * Thresholds * Neighborhoods * Frac Thresholds *
         //    Alphas
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_cat_thresh() * get_n_nbrhd_wdth() *
              get_n_cov_thresh() * get_n_ci_alpha());
         break;

      case(i_nbrcnt):
         // Number of NBRCNT lines =
         //    Masks * Thresholds * Neighborhoods * Alphas
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_cat_thresh() * get_n_nbrhd_wdth() *
              get_n_ci_alpha());
         break;

      case(i_pct):
      case(i_pjc):
      case(i_prc):
         // Number of PCT, PJC, or PRC lines =
         //    Masks * Smoothing Methods * Thresholds * Climo Bins
         n = (!prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_oprob_thresh() *
              n_bin);
         break;

      case(i_pstd):
         // Number of PSTD lines =
         //    Masks * Smoothing Methods * Thresholds * Alphas *
         //    Climo Bins
         n = (!prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_oprob_thresh() *
              get_n_ci_alpha() * n_bin);
         break;

      case(i_eclv):
         // Number of CTC -> ECLV lines =
         //    Masks * Smoothing Methods * Thresholds
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_cat_thresh());

         // Number of PCT -> ECLV lines =
         //    Probability Fields * Masks * Smoothing Methods *
         //    Max Observation Probability Thresholds *
         //    Max Forecast Probability Thresholds * Climo Bins
         n += (!prob_flag ? 0 :
               get_n_mask() * get_n_interp() * get_n_oprob_thresh() *
               get_n_fprob_thresh() * n_bin);
         break;

      case(i_grad):
         // Number of GRAD lines =
         //    Masks * Smoothing Methods * Gradient Sizes
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_grad());
         break;

      case(i_dmap):
         // Number of DMAP lines =
         //    Masks * Smoothing Methods * Thresholds
         n = (prob_flag ? 0 :
              get_n_mask() * get_n_interp() * get_n_cat_thresh());
         break;

      default:
         mlog << Error << "\nGridStatVxOpt::n_txt_row(int) -> "
              << "unexpected output type index value: " << i_txt_row
              << "\n\n";
         exit(1);
         break;
   }

   return(n);
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::get_n_cnt_thresh() const {
   return((!fcst_info || fcst_info->is_prob()) ?
          0 : fcnt_ta.n_elements());
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::get_n_cat_thresh() const {
   return((!fcst_info || fcst_info->is_prob()) ?
          0 : fcat_ta.n_elements());
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::get_n_wind_thresh() const {
   return((!fcst_info || fcst_info->is_prob()) ?
          0 : fwind_ta.n_elements());
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::get_n_fprob_thresh() const {
   return((!fcst_info || !fcst_info->is_prob()) ?
          0 : fcat_ta.n_elements());
}

////////////////////////////////////////////////////////////////////////

int GridStatVxOpt::get_n_oprob_thresh() const {
   return((!fcst_info || !fcst_info->is_prob()) ?
          0 : ocat_ta.n_elements());
}



////////////////////////////////////////////////////////////////////////
//
//  Code for struct GridStatNcOutInfo
//
////////////////////////////////////////////////////////////////////////

GridStatNcOutInfo::GridStatNcOutInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

void GridStatNcOutInfo::clear() {

   set_all_true();

   return;
}

////////////////////////////////////////////////////////////////////////

bool GridStatNcOutInfo::all_false() const {

   bool status = do_latlon       || do_raw        || do_diff     ||
                 do_climo        || do_climo_cdp  || do_weight   ||
                 do_nbrhd        || do_fourier    || do_gradient ||
                 do_distance_map || do_apply_mask;

   return(!status);
}

////////////////////////////////////////////////////////////////////////

void GridStatNcOutInfo::set_all_false() {

   do_latlon       = false;
   do_raw          = false;
   do_diff         = false;
   do_climo        = false;
   do_climo_cdp    = false;
   do_weight       = false;
   do_nbrhd        = false;
   do_fourier      = false;
   do_gradient     = false;
   do_distance_map = false;
   do_apply_mask   = false;

   return;
}

////////////////////////////////////////////////////////////////////////

void GridStatNcOutInfo::set_all_true() {

   do_latlon       = true;
   do_raw          = true;
   do_diff         = true;
   do_climo        = true;
   do_climo_cdp    = true;
   do_weight       = true;
   do_nbrhd        = true;
   do_fourier      = true;
   do_gradient     = true;
   do_distance_map = true;
   do_apply_mask   = true;

   return;
}

////////////////////////////////////////////////////////////////////////
