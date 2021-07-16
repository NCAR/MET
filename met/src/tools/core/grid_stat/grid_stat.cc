// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////
//
//   Filename:   grid_stat.cc
//
//   Description:
//      Based on user-specified parameters, this tool compares a
//      a gridded forecast fields to gridded analyses. It computes
//      many verification scores and statistics, including confidence
//      intervals, to summarize the comparison.
//
//   Mod#   Date      Name           Description
//   ----   ----      ----           -----------
//   000    04/18/07  Halley Gotway  New
//   001    10/05/07  Halley Gotway  Add checks to only compute
//                    confidence intervals for sufficiently large
//                    sample sizes.
//   002    10/19/07  Halley Gotway  Add ability to smooth the forecast
//                    field using the interp_method and interp_width
//                    parameters
//   003    01/11/08  Halley Gotway  Modify sprintf statements which
//                    use the GRIB code abbreviation string
//   004    02/06/08  Halley Gotway  Modify to read the updated NetCDF
//                    output of PCP-Combine
//   005    02/11/08  Halley Gotway  Remove BIAS from the CNT line
//                    since it's the same as the ME.
//   006    02/12/08  Halley Gotway  Fix bug in writing COP line to
//                    write out OY and ON rather than FY and FN.
//   007    02/14/08  Halley Gotway  Apply the smoothing operation to
//                    both forecast and observation fields rather than
//                    just the forecast.
//   008    02/25/08  Halley Gotway  Write the STAT lines using the
//                    routines in the vx_met_util library.
//   009    02/25/08  Halley Gotway  Add support for the NBRCTS and
//                    NBRCNT line types.
//   010    07/01/08  Halley Gotway  Add the rank_corr_flag to the
//                    config file to disable computing rank
//                    correlations.
//   011    09/23/08  Halley Gotway  Change argument sequence for the
//                    get_grib_record routine.
//   012    10/31/08  Halley Gotway  Use the get_file_type routine to
//                    determine the input file types.
//   013    11/14/08  Halley Gotway  Restructure organization of Grid-
//                    Stat and move several functions to libraries.
//   014    03/13/09  Halley Gotway  Add support for verifying
//                    probabilistic forecasts.
//   015    10/22/09  Halley Gotway  Fix output_flag cut and paste bug.
//   016    05/03/10  Halley Gotway  Don't write duplicate NetCDF
//                    matched pair variables or probabilistic
//                    difference fields.
//   017    05/06/10  Halley Gotway  Add interp_flag config parameter.
//   018    05/27/10  Halley Gotway  Add -fcst_valid, -fcst_lead,
//                    -obs_valid, and -obs_lead command line options.
//   019    06/08/10  Halley Gotway  Add support for multi-category
//                    contingency tables.
//   020    06/30/10  Halley Gotway  Enhance grid equality checks.
//   021    07/27/10  Halley Gotway  Add lat/lon variables to NetCDF.
//   022    10/28/11  Holmes         Added use of command line class to
//                    parse the command line arguments.
//   023    11/14/11  Holmes         Added code to enable reading of
//                    multiple config files.
//   024    12/09/11  Halley Gotway  When zero matched pairs are found
//                    print the status message before continuing.
//   025    02/02/12  Bullock        Changed default output directory to "."
//   026    04/30/12  Halley Gotway  Switch to using vx_config library.
//   027    04/30/12  Halley Gotway  Move -fcst_valid, -fcst_lead,
//                    -obs_valid, and -obs_lead command line options
//                    to config file.
//   028    08/21/13  Halley Gotway  Fix sizing of output tables for 12
//                    or more probabilstic thresholds.
//   029    05/20/14  Halley Gotway  Add AFSS, UFSS, F_RATE, and O_RATE
//                    to the NBRCNT line type.
//   030    11/12/14  Halley Gotway  Pass the obtype entry from the
//                    from the config file to the output files.
//   031    02/19/15  Halley Gotway  Add automated regridding.
//   032    08/04/15  Halley Gotway  Add conditional continuous verification.
//   033    09/04/15  Halley Gotway  Add climatology and SAL1L2 and VAL1L2
//                    output line types.
//   034    05/10/16  Halley Gotway  Add grid weighting.
//   035    05/20/16  Prestopnik J   Removed -version (now in command_line.cc)
//   036    02/14/17  Win            MET-621 enhancement support- additional
//                    nc_pairs_flag 'apply_mask'
//   037    05/15/17  Prestopnik P   Add shape for regrid, nbrhd and interp
//   038    06/26/17  Halley Gotway  Add ECLV line type.
//   039    08/18/17  Halley Gotway  Add fourier decomposition.
//   040    08/18/17  Halley Gotway  Add GRAD output line type.
//   041    11/01/17  Halley Gotway  Add binned climatologies.
//   042    12/15/17  Halley Gotway  Restructure config logic to make
//                    all options settable for each verification task.
//   043    02/14/17  Halley Gotway  Add nbrhd.field option to skip the
//                    computation of fractional coverage fields.
//   044    02/22/19  Halley Gotway  Make gradient dx/dy configurable.
//   045    04/08/19  Halley Gotway  Add percentile thresholds.
//   046    04/01/19  Fillmore       Add FCST and OBS units.
//   047    06/19/19  Halley Gotway  Add DMAP output line type.
//   048    10/14/19  Halley Gotway  Add support for climo distribution
//                    percentile thresholds.
//   049    11/15/19  Halley Gotway  Apply climatology bins to
//                    continuous and probabilistic statistics.
//   050    03/02/20  Halley Gotway  Add nc_pairs_var_name and rename
//                    nc_pairs_var_str to nc_pairs_var_suffix.
//   051    03/28/21  Halley Gotway  Add mpr_column and mpr_thresh
//                    filtering options.
//   052    05/28/21  Halley Gotway  Add MCTS HSS_EC output.
//
////////////////////////////////////////////////////////////////////////

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <ctype.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "grid_stat.h"

#include "vx_statistics.h"
#include "vx_nc_util.h"
#include "vx_regrid.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////

static void process_command_line(int, char **);
static void setup_first_pass    (const DataPlane &);

static void setup_txt_files(unixtime, int);
static void setup_table    (AsciiTable &);
static void setup_nc_file  (const GridStatNcOutInfo &, unixtime, int);

static void build_outfile_name(unixtime, int, const char *,
                               ConcatString &);

static void process_scores();

static void get_mask_points(const GridStatVxOpt &,
                            const MaskPlane &, const DataPlane *,
                            const DataPlane *, const DataPlane *,
                            const DataPlane *, const DataPlane *,
                            PairDataPoint &);

static void do_cts       (CTSInfo *&,   int, const PairDataPoint *);
static void do_mcts      (MCTSInfo &,   int, const PairDataPoint *);
static void do_cnt_sl1l2 (const GridStatVxOpt &, const PairDataPoint *);
static void do_vl1l2     (VL1L2Info *&, int, const PairDataPoint *, const PairDataPoint *);
static void do_pct       (const GridStatVxOpt &, const PairDataPoint *);

static void do_nbrcts(NBRCTSInfo *&, int, int, int, const PairDataPoint *);
static void do_nbrcnt(NBRCNTInfo &,  int, int, int, const PairDataPoint *, const PairDataPoint *);

static void write_nc(const ConcatString &, const DataPlane &, int,
                     const ConcatString &, int, FieldType);
static void write_nbrhd_nc(const DataPlane &, const DataPlane &, int,
                           const SingleThresh &, const SingleThresh &,
                           int, int);

static void add_var_att_local(NcVar *, const char *, const ConcatString);

static void finish_txt_files();

static void clean_up();

static void usage();
static void set_outdir(const StringArray &);
static void set_compress(const StringArray &);
static bool read_data_plane(VarInfo* info, DataPlane& dp, Met2dDataFile* mtddf,
                            const ConcatString &filename);

////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

   // Set handler to be called for memory allocation error
   set_new_handler(oom);

   // Process the command line arguments
   process_command_line(argc, argv);

   // Compute the scores and write them out
   process_scores();

   // Close the text files and deallocate memory
   clean_up();

   return(0);
}

////////////////////////////////////////////////////////////////////////

void process_command_line(int argc, char **argv) {
   CommandLine cline;
   GrdFileType ftype, otype;
   ConcatString default_config_file;
   DataPlane dp;

   // Set the default output directory
   out_dir = replace_path(default_out_dir);

   // Check for zero arguments
   if(argc == 1) usage();

   // Parse the command line into tokens
   cline.set(argc, argv);

   // Set the usage function
   cline.set_usage(usage);

   // Add the options function calls
   cline.add(set_outdir,   "-outdir",   1);
   cline.add(set_compress, "-compress", 1);

   // Parse the command line
   cline.parse();

   // Check for error. There should be three arguments left:
   // forecast, observation, and config filenames
   if(cline.n() != 3) usage();

   // Store the input file names
   fcst_file   = cline[0];
   obs_file    = cline[1];
   config_file = cline[2];

   // Create the default config file name
   default_config_file = replace_path(default_config_filename);

   // List the config files
   mlog << Debug(1)
        << "Default Config File: " << default_config_file << "\n"
        << "User Config File: "    << config_file << "\n";

   // Read the config files
   conf_info.read_config(default_config_file.c_str(), config_file.c_str());

   // Get the forecast and observation file types from config, if present
   ftype = parse_conf_file_type(conf_info.conf.lookup_dictionary(conf_key_fcst));
   otype = parse_conf_file_type(conf_info.conf.lookup_dictionary(conf_key_obs));

   // Read forecast file
   if(!(fcst_mtddf = mtddf_factory.new_met_2d_data_file(fcst_file.c_str(), ftype))) {
      mlog << Error << "\nTrouble reading forecast file \""
           << fcst_file << "\"\n\n";
      exit(1);
   }

   // Read observation file
   if(!(obs_mtddf = mtddf_factory.new_met_2d_data_file(obs_file.c_str(), otype))) {
      mlog << Error << "\nTrouble reading observation file \""
           << obs_file << "\"\n\n";
      exit(1);
   }

   // Store the input data file types
   ftype = fcst_mtddf->file_type();
   otype = obs_mtddf->file_type();

   // Process the configuration
   conf_info.process_config(ftype, otype);

   // For python types read the first field to set the grid
   if(is_python_grdfiletype(ftype)) {
      if(!fcst_mtddf->data_plane(*conf_info.vx_opt[0].fcst_info, dp)) {
         mlog << Error << "\nTrouble reading data from forecast file \""
              << fcst_file << "\"\n\n";
         exit(1);
      }
   }

   if(is_python_grdfiletype(otype)) {
      if(!obs_mtddf->data_plane(*conf_info.vx_opt[0].obs_info, dp)) {
         mlog << Error << "\nTrouble reading data from observation file \""
              << obs_file << "\"\n\n";
         exit(1);
      }
   }

   // Determine the verification grid
   grid = parse_vx_grid(conf_info.vx_opt[0].fcst_info->regrid(),
                        &(fcst_mtddf->grid()), &(obs_mtddf->grid()));

   // Compute weight for each grid point
   parse_grid_weight(grid, conf_info.grid_weight_flag, wgt_dp);

   // Set the model name
   shc.set_model(conf_info.model.c_str());

   // Set the obtype column
   shc.set_obtype(conf_info.obtype.c_str());

   // Use the first verification task to set the random number generator
   // and seed value for bootstrap confidence intervals
   rng_set(rng_ptr,
           conf_info.vx_opt[0].boot_info.rng.c_str(),
           conf_info.vx_opt[0].boot_info.seed.c_str());

   // List the input files
   mlog << Debug(1)
        << "Forecast File: "    << fcst_file << "\n"
        << "Observation File: " << obs_file  << "\n";

   return;
}

////////////////////////////////////////////////////////////////////////

void setup_first_pass(const DataPlane &dp) {

   // Unset the flag
   is_first_pass = false;

   // Process masks grids and polylines in the config file
   conf_info.process_masks(grid);

   // Create output text files as requested in the config file
   // only if the ascii_output_flag is true.
   if(conf_info.get_output_ascii_flag()){
      setup_txt_files(dp.valid(), dp.lead());
   }

   // If requested, create a NetCDF file to store the gridded
   // verification data
   if(conf_info.get_output_nc_flag()) {
      setup_nc_file(conf_info.nc_info, dp.valid(), dp.lead());
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void setup_txt_files(unixtime valid_ut, int lead_sec) {
   int  i, max_col, max_prob_col, max_mctc_col, n_prob, n_cat, n_eclv;
   ConcatString base_name;

   // Create output file names for the stat file and optional text files
   build_outfile_name(valid_ut, lead_sec, "", base_name);

   /////////////////////////////////////////////////////////////////////
   //
   // Setup the output STAT file
   //
   /////////////////////////////////////////////////////////////////////

   // Get the maximum number of data columns
   n_prob = conf_info.get_max_n_fprob_thresh();
   n_cat  = conf_info.get_max_n_cat_thresh() + 1;
   n_eclv = conf_info.get_max_n_eclv_points();

   max_prob_col = get_n_pjc_columns(n_prob);
   max_mctc_col = get_n_mctc_columns(n_cat);

   // Determine the maximum number of data columns
   max_col = ( max_prob_col > max_stat_col ? max_prob_col : max_stat_col );
   max_col = ( max_mctc_col > max_col      ? max_mctc_col : max_col );

   // Add the header columns
   max_col += n_header_columns + 1;

   // Initialize file stream
   stat_out = (ofstream *) 0;

   // Build the file name
   stat_file << base_name << stat_file_ext;

   // Create the output STAT file
   open_txt_file(stat_out, stat_file.c_str());

   // Setup the STAT AsciiTable
   stat_at.set_size(conf_info.n_stat_row() + 1, max_col);
   setup_table(stat_at);

   // Write the text header row
   write_header_row((const char **) 0, 0, 1, stat_at, 0, 0);

   // Initialize the row index to 1 to account for the header
   i_stat_row = 1;

   /////////////////////////////////////////////////////////////////////
   //
   // Setup each of the optional output text files
   //
   /////////////////////////////////////////////////////////////////////

   // Loop through output file type
   for(i=0; i<n_txt; i++) {

      // Only set it up if requested in the config file
      if(conf_info.output_flag[i] == STATOutputType_Both) {

         // Initialize file stream
         txt_out[i] = (ofstream *) 0;

         // Build the file name
         txt_file[i] << base_name << "_" << txt_file_abbr[i]
                     << txt_file_ext;

         // Create the output text file
         open_txt_file(txt_out[i], txt_file[i].c_str());

         // Get the maximum number of columns for this line type
         switch(i) {

            case(i_mctc):
               max_col = get_n_mctc_columns(n_cat) + n_header_columns + 1;
               break;

            case(i_pct):
               max_col = get_n_pct_columns(n_prob)  + n_header_columns + 1;
               break;

            case(i_pstd):
               max_col = get_n_pstd_columns(n_prob) + n_header_columns + 1;
               break;

            case(i_pjc):
               max_col = get_n_pjc_columns(n_prob)  + n_header_columns + 1;
               break;

            case(i_prc):
               max_col = get_n_prc_columns(n_prob)  + n_header_columns + 1;
               break;

            case(i_eclv):
               max_col = get_n_eclv_columns(n_eclv)  + n_header_columns + 1;
               break;

            default:
               max_col = n_txt_columns[i]  + n_header_columns + 1;
               break;
         } // end switch

         // Setup the text AsciiTable
         txt_at[i].set_size(conf_info.n_txt_row(i) + 1, max_col);
         setup_table(txt_at[i]);

         // Write the text header row
         switch(i) {

            case(i_mctc):
               write_mctc_header_row(1, n_cat, txt_at[i], 0, 0);
               break;

            case(i_pct):
               write_pct_header_row(1, n_prob, txt_at[i], 0, 0);
               break;

            case(i_pstd):
               write_pstd_header_row(1, n_prob, txt_at[i], 0, 0);
               break;

            case(i_pjc):
               write_pjc_header_row(1, n_prob, txt_at[i], 0, 0);
               break;

            case(i_prc):
               write_prc_header_row(1, n_prob, txt_at[i], 0, 0);
               break;

            case(i_eclv):
               write_eclv_header_row(1, n_eclv, txt_at[i], 0, 0);
               break;

            default:
               write_header_row(txt_columns[i], n_txt_columns[i], 1,
                  txt_at[i], 0, 0);
               break;
         } // end switch

         // Initialize the row index to 1 to account for the header
         i_txt_row[i] = 1;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void setup_table(AsciiTable &at) {

   // Justify the STAT AsciiTable objects
   justify_stat_cols(at);

   // Set the precision
   at.set_precision(conf_info.conf.output_precision());

   // Set the bad data value
   at.set_bad_data_value(bad_data_double);

   // Set the bad data string
   at.set_bad_data_str(na_str);

   // Don't write out trailing blank rows
   at.set_delete_trailing_blank_rows(1);

   return;
}

////////////////////////////////////////////////////////////////////////

void setup_nc_file(const GridStatNcOutInfo & nc_info,
                   unixtime valid_ut, int lead_sec) {

   // Create output NetCDF file name
   build_outfile_name(valid_ut, lead_sec, "_pairs.nc", out_nc_file);

   // Create a new NetCDF file and open it
   nc_out = open_ncfile(out_nc_file.c_str(), true);

   if(IS_INVALID_NC_P(nc_out)) {
      mlog << Error << "\nsetup_nc_file() -> "
           << "trouble opening output NetCDF file "
           << out_nc_file << "\n\n";
      exit(1);
   }

   // Add global attributes
   write_netcdf_global(nc_out, out_nc_file.c_str(), program_name,
                       conf_info.model.c_str(), conf_info.obtype.c_str());
   if(nc_info.do_diff) {
      add_att(nc_out, "Difference", "Forecast Value - Observation Value");
   }

   // Add the projection information
   write_netcdf_proj(nc_out, grid);

   // Define Dimensions
   lat_dim = add_dim(nc_out,"lat", (long) grid.ny());
   lon_dim = add_dim(nc_out,"lon", (long) grid.nx());

   // Add the lat/lon variables
   if(nc_info.do_latlon) {
      write_netcdf_latlon(nc_out, &lat_dim, &lon_dim, grid);
   }

   // Add grid weight variable
   if(nc_info.do_weight) {
      write_netcdf_grid_weight(nc_out, &lat_dim, &lon_dim,
                               conf_info.grid_weight_flag, wgt_dp);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void build_outfile_name(unixtime valid_ut, int lead_sec,
                        const char *suffix, ConcatString &str) {

   //
   // Create output file name
   //

   // Append the output directory and program name
   str << cs_erase << out_dir << "/" << program_name;

   // Append the output prefix, if defined
   if(conf_info.output_prefix.nonempty())
      str << "_" << conf_info.output_prefix;

   // Append the timing information
   str << "_"
       << sec_to_hhmmss(lead_sec) << "L_"
       << unix_to_yyyymmdd_hhmmss(valid_ut) << "V";

   // Append the suffix
   str << suffix;

   return;
}

////////////////////////////////////////////////////////////////////////

void process_scores() {
   int i, j, k, m, n;
   int n_cat, n_wind, n_cov;
   double dmin, dmax;
   ConcatString cs;

   // Forecast and observation fields
   MaskPlane mask_mp;
   DataPlane fcst_dp,        obs_dp;
   DataPlane fcst_dp_smooth, obs_dp_smooth;
   DataPlane fcst_dp_thresh, obs_dp_thresh;

   // Climatology mean and standard deviation
   DataPlane cmn_dp, csd_dp;
   DataPlane cmn_dp_smooth;

   // Paired forecast, observation, climatology, and weight values
   PairDataPoint pd;

   // Thresholded fractional coverage pairs
   PairDataPoint pd_thr;

   // Allocate memory in one big chunk based on grid size
   pd.extend(grid.nx()*grid.ny());

   if(conf_info.output_flag[i_nbrctc] != STATOutputType_None ||
      conf_info.output_flag[i_nbrcts] != STATOutputType_None ||
      conf_info.output_flag[i_nbrcnt] != STATOutputType_None ||
      conf_info.output_flag[i_dmap]   != STATOutputType_None) {
      pd_thr.extend(grid.nx()*grid.ny());
   }

   // Objects to handle vector winds
   DataPlane fu_dp, ou_dp;
   DataPlane fu_dp_smooth, ou_dp_smooth;
   DataPlane cmnu_dp, csdu_dp, cmnu_dp_smooth;
   PairDataPoint pd_u;

   CTSInfo    *cts_info    = (CTSInfo *) 0;
   MCTSInfo    mcts_info;
   VL1L2Info  *vl1l2_info  = (VL1L2Info *) 0;
   NBRCNTInfo  nbrcnt_info;
   NBRCTSInfo *nbrcts_info = (NBRCTSInfo *) 0;
   GRADInfo    grad_info;
   DMAPInfo    dmap_info;

   // Store the maximum number of each threshold type
   n_cat  = conf_info.get_max_n_cat_thresh();
   n_wind = conf_info.get_max_n_wind_thresh();
   n_cov  = conf_info.get_max_n_cov_thresh();

   // Allocate space for output statistics types
   cts_info    = new CTSInfo    [n_cat];
   vl1l2_info  = new VL1L2Info  [n_wind];
   nbrcts_info = new NBRCTSInfo [n_cov];

   // Compute scores for each verification task and write output_flag
   for(i=0; i<conf_info.get_n_vx(); i++) {

      // Read the gridded data from the input forecast file
      if(!read_data_plane(conf_info.vx_opt[i].fcst_info,
                          fcst_dp, fcst_mtddf, fcst_file)) continue;

      mlog << Debug(3)
           << "Reading forecast data for "
           << conf_info.vx_opt[i].fcst_info->magic_str() << ".\n";

      // Set the forecast lead time
      shc.set_fcst_lead_sec(fcst_dp.lead());

      // Set the forecast valid time
      shc.set_fcst_valid_beg(fcst_dp.valid());
      shc.set_fcst_valid_end(fcst_dp.valid());

      // Read the gridded data from the input observation file
      if(!read_data_plane(conf_info.vx_opt[i].obs_info,
                          obs_dp, obs_mtddf, obs_file)) continue;

      mlog << Debug(3)
           << "Reading observation data for "
           << conf_info.vx_opt[i].obs_info->magic_str() << ".\n";

      // Set the observation lead time
      shc.set_obs_lead_sec(obs_dp.lead());

      // Set the observation valid time
      shc.set_obs_valid_beg(obs_dp.valid());
      shc.set_obs_valid_end(obs_dp.valid());

      // Check that the valid times match
      if(fcst_dp.valid() != obs_dp.valid()) {

         mlog << Warning << "\nprocess_scores() -> "
              << "Forecast and observation valid times do not match "
              << unix_to_yyyymmdd_hhmmss(fcst_dp.valid()) << " != "
              << unix_to_yyyymmdd_hhmmss(obs_dp.valid()) << " for "
              << conf_info.vx_opt[i].fcst_info->magic_str() << " versus "
              << conf_info.vx_opt[i].obs_info->magic_str() << ".\n\n";
      }

      // Check that the accumulation intervals match
      if(conf_info.vx_opt[i].fcst_info->level().type() == LevelType_Accum &&
         conf_info.vx_opt[i].obs_info->level().type()  == LevelType_Accum &&
         fcst_dp.accum() != obs_dp.accum()) {

         mlog << Warning << "\nprocess_scores() -> "
              << "Forecast and observation accumulation times "
              << "do not match " << sec_to_hhmmss(fcst_dp.accum())
              << " != " << sec_to_hhmmss(obs_dp.accum())
              << " for " << conf_info.vx_opt[i].fcst_info->magic_str()
              << " versus " << conf_info.vx_opt[i].obs_info->magic_str()
              << ".\n\n";
      }

      // Read climatology data
      cmn_dp = read_climo_data_plane(
                  conf_info.conf.lookup_array(conf_key_climo_mean_field, false),
                  i, fcst_dp.valid(), grid);
      csd_dp = read_climo_data_plane(
                  conf_info.conf.lookup_array(conf_key_climo_stdev_field, false),
                  i, fcst_dp.valid(), grid);

      mlog << Debug(3)
           << "Found " << (cmn_dp.nx() == 0 ? 0 : 1)
           << " climatology mean and " << (csd_dp.nx() == 0 ? 0 : 1)
           << " climatology standard deviation field(s) for forecast "
           << conf_info.vx_opt[i].fcst_info->magic_str() << ".\n";

      // Apply MPR threshold filters
      if(conf_info.vx_opt[i].mpr_sa.n() > 0) {
         apply_mpr_thresh_mask(fcst_dp, obs_dp, cmn_dp, csd_dp,
                               conf_info.vx_opt[i].mpr_sa,
                               conf_info.vx_opt[i].mpr_ta);
      }

      // Setup the first pass through the data
      if(is_first_pass) setup_first_pass(fcst_dp);

      // Store the description
      shc.set_desc(conf_info.vx_opt[i].desc.c_str());

      // Store the forecast variable name
      shc.set_fcst_var(conf_info.vx_opt[i].fcst_info->name_attr());

      // Store the forecast variable units
      shc.set_fcst_units(conf_info.vx_opt[i].fcst_info->units_attr());

      // Set the forecast level name
      shc.set_fcst_lev(conf_info.vx_opt[i].fcst_info->level_attr().c_str());

      // Store the observation variable name
      shc.set_obs_var(conf_info.vx_opt[i].obs_info->name_attr());

      // Store the observation variable units
      shc.set_obs_units(conf_info.vx_opt[i].obs_info->units_attr());

      // Set the observation level name
      shc.set_obs_lev(conf_info.vx_opt[i].obs_info->level_attr().c_str());

      mlog << Debug(2) << "\n" << sep_str << "\n\n";

      // Pointer to current interpolation object
      InterpInfo * interp = &conf_info.vx_opt[i].interp_info;

      // Loop through and apply each of the smoothing operations
      for(j=0; j<conf_info.vx_opt[i].get_n_interp(); j++) {

         // Current interpolation method
         InterpMthd interp_mthd = string_to_interpmthd(interp->method[j].c_str());

         // Create grid template to find the number of points
         GridTemplateFactory gtf;
         GridTemplate* gt = gtf.buildGT(interp->shape, interp->width[j]);

         shc.set_interp_mthd(interp_mthd, interp->shape);
         int interp_pnts = gt->size();
         shc.set_interp_pnts(interp_pnts);
         delete gt;

         // If requested in the config file, smooth the forecast field
         if(interp->field == FieldType_Fcst ||
            interp->field == FieldType_Both) {
            smooth_field(fcst_dp, fcst_dp_smooth, interp_mthd,
                         interp->width[j], interp->shape, interp->vld_thresh,
                         interp->gaussian);
         }
         // Do not smooth the forecast field
         else {
            fcst_dp_smooth = fcst_dp;
         }

         // If requested in the config file, smooth the observation field
         if(interp->field == FieldType_Obs ||
            interp->field == FieldType_Both) {
            smooth_field(obs_dp, obs_dp_smooth, interp_mthd,
                         interp->width[j], interp->shape, interp->vld_thresh,
                         interp->gaussian);
         }
         // Do not smooth the observation field
         else {
            obs_dp_smooth = obs_dp;
         }

         // Loop through the masks to be applied
         for(k=0; k<conf_info.vx_opt[i].get_n_mask(); k++) {

            // Store the current mask
            mask_mp = conf_info.mask_map[conf_info.vx_opt[i].mask_name[k]];

            // Turn off the mask for missing data values
            mask_bad_data(mask_mp, fcst_dp_smooth);
            mask_bad_data(mask_mp, obs_dp_smooth);
            if(cmn_dp.nx() == fcst_dp_smooth.nx() &&
               cmn_dp.ny() == fcst_dp_smooth.ny()) {
               mask_bad_data(mask_mp, cmn_dp);
            }
            if(csd_dp.nx() == fcst_dp_smooth.nx() &&
               csd_dp.ny() == fcst_dp_smooth.ny()) {
               mask_bad_data(mask_mp, csd_dp);
            }

            // Apply the current mask to the current fields
            get_mask_points(conf_info.vx_opt[i], mask_mp,
                            &fcst_dp_smooth, &obs_dp_smooth,
                            &cmn_dp, &csd_dp, &wgt_dp, pd);

            // Set the mask name
            shc.set_mask(conf_info.vx_opt[i].mask_name[k].c_str());

            mlog << Debug(2)
                 << "Processing " << conf_info.vx_opt[i].fcst_info->magic_str()
                 << " versus " << conf_info.vx_opt[i].obs_info->magic_str()
                 << ", for smoothing method "
                 << shc.get_interp_mthd() << "("
                 << shc.get_interp_pnts_str()
                 << "), over region " << shc.get_mask()
                 << ", using " << pd.f_na.n() << " matched pairs.\n";

            // Continue if no pairs were found
            if(pd.f_na.n() == 0) continue;

            // Process percentile thresholds
            conf_info.vx_opt[i].set_perc_thresh(pd);

            // Compute CTS scores
            if(!conf_info.vx_opt[i].fcst_info->is_prob()                       &&
                conf_info.vx_opt[i].fcat_ta.n() > 0                            &&
               (conf_info.vx_opt[i].output_flag[i_fho]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_ctc]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_cts]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_eclv] != STATOutputType_None)) {

               // Initialize
               for(m=0; m<n_cat; m++) cts_info[m].clear();

               // Compute CTS
               do_cts(cts_info, i, &pd);

               // Loop through all of the thresholds
               for(m=0; m<conf_info.vx_opt[i].fcat_ta.n(); m++) {

                  // Write out FHO
                  if(conf_info.vx_opt[i].output_flag[i_fho] != STATOutputType_None &&
                     cts_info[m].cts.n() > 0) {

                     write_fho_row(shc, cts_info[m],
                        conf_info.vx_opt[i].output_flag[i_fho],
                        stat_at, i_stat_row,
                        txt_at[i_fho], i_txt_row[i_fho]);
                  }

                  // Write out CTC
                  if(conf_info.vx_opt[i].output_flag[i_ctc] != STATOutputType_None &&
                     cts_info[m].cts.n() > 0) {

                     write_ctc_row(shc, cts_info[m],
                        conf_info.vx_opt[i].output_flag[i_ctc],
                        stat_at, i_stat_row,
                        txt_at[i_ctc], i_txt_row[i_ctc]);
                  }

                  // Write out CTS
                  if(conf_info.vx_opt[i].output_flag[i_cts] != STATOutputType_None &&
                     cts_info[m].cts.n() > 0) {

                     write_cts_row(shc, cts_info[m],
                        conf_info.vx_opt[i].output_flag[i_cts],
                        stat_at, i_stat_row,
                        txt_at[i_cts], i_txt_row[i_cts]);
                  }

                  // Write out ECLV
                  if(conf_info.vx_opt[i].output_flag[i_eclv] != STATOutputType_None &&
                     cts_info[m].cts.n() > 0) {

                     write_eclv_row(shc, cts_info[m], conf_info.vx_opt[i].eclv_points,
                        conf_info.vx_opt[i].output_flag[i_eclv],
                        stat_at, i_stat_row,
                        txt_at[i_eclv], i_txt_row[i_eclv]);
                  }
               } // end for m
            } // end Compute CTS

            // Compute MCTS scores
            if(!conf_info.vx_opt[i].fcst_info->is_prob()                       &&
                conf_info.vx_opt[i].fcat_ta.n() > 1                            &&
               (conf_info.vx_opt[i].output_flag[i_mctc] != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_mcts] != STATOutputType_None)) {

               // Initialize
               mcts_info.clear();

               // Compute MCTS
               do_mcts(mcts_info, i, &pd);

               // Write out MCTC
               if(conf_info.vx_opt[i].output_flag[i_mctc] != STATOutputType_None &&
                  mcts_info.cts.total() > 0) {

                  write_mctc_row(shc, mcts_info,
                     conf_info.vx_opt[i].output_flag[i_mctc],
                     stat_at, i_stat_row,
                     txt_at[i_mctc], i_txt_row[i_mctc]);
               }

               // Write out MCTS
               if(conf_info.vx_opt[i].output_flag[i_mcts] != STATOutputType_None &&
                  mcts_info.cts.total() > 0) {

                  write_mcts_row(shc, mcts_info,
                     conf_info.vx_opt[i].output_flag[i_mcts],
                     stat_at, i_stat_row,
                     txt_at[i_mcts], i_txt_row[i_mcts]);
               }
            } // end Compute MCTS

            // Compute CNT, SL1L2, and SAL1L2 scores
            if(!conf_info.vx_opt[i].fcst_info->is_prob() &&
               (conf_info.vx_opt[i].output_flag[i_cnt]    != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_sl1l2]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_sal1l2] != STATOutputType_None)) {
               do_cnt_sl1l2(conf_info.vx_opt[i], &pd);
            }

            // Compute VL1L2 and VAL1L2 partial sums for UGRD,VGRD
            if(!conf_info.vx_opt[i].fcst_info->is_prob()                         &&
                conf_info.vx_opt[i].fcst_info->is_v_wind()                       &&
                conf_info.vx_opt[i].fcst_info->uv_index() >= 0                   &&
               (conf_info.vx_opt[i].output_flag[i_vl1l2]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_val1l2] != STATOutputType_None) ) {

               // Store the forecast variable name
               shc.set_fcst_var(ugrd_vgrd_abbr_str);

               // Store the observation variable name
               shc.set_obs_var(ugrd_vgrd_abbr_str);

               // Initialize
               for(m=0; m<n_wind; m++) vl1l2_info[m].clear();

               // Index for UGRD
               int ui = conf_info.vx_opt[i].fcst_info->uv_index();

               // Read forecast data for UGRD
               if(!read_data_plane(conf_info.vx_opt[ui].fcst_info,
                                   fu_dp, fcst_mtddf, fcst_file)) continue;

               // Read observation data for UGRD
               if(!read_data_plane(conf_info.vx_opt[ui].obs_info,
                                   ou_dp, obs_mtddf, obs_file)) continue;

               // Read climatology data for UGRD
               cmnu_dp = read_climo_data_plane(
                           conf_info.conf.lookup_array(conf_key_climo_mean_field, false),
                           ui, fcst_dp.valid(), grid);
               csdu_dp = read_climo_data_plane(
                           conf_info.conf.lookup_array(conf_key_climo_stdev_field, false),
                           ui, fcst_dp.valid(), grid);

               // If requested in the config file, smooth the forecast
               // and climatology U-wind fields
               if(interp->field == FieldType_Fcst ||
                  interp->field == FieldType_Both) {
                  smooth_field(fu_dp, fu_dp_smooth, interp_mthd,
                               interp->width[j], interp->shape, interp->vld_thresh,
                               interp->gaussian);
               }
               // Do not smooth the forecast field
               else {
                  fu_dp_smooth = fu_dp;
               }

               // If requested in the config file, smooth the observation
               // U-wind field
               if(interp->field == FieldType_Obs ||
                  interp->field == FieldType_Both) {
                  smooth_field(ou_dp, ou_dp_smooth,
                               interp_mthd, interp->width[j],
                               interp->shape, interp->vld_thresh,
                               interp->gaussian);
               }
               // Do not smooth the observation field
               else {
                  ou_dp_smooth = ou_dp;
               }

               // Apply the current mask to the U-wind fields
               get_mask_points(conf_info.vx_opt[i], mask_mp,
                               &fu_dp_smooth, &ou_dp_smooth,
                               &cmnu_dp, &csdu_dp, &wgt_dp, pd_u);

               // Compute VL1L2
               do_vl1l2(vl1l2_info, i, &pd_u, &pd);

               // Loop through all of the wind speed thresholds
               for(m=0; m<conf_info.vx_opt[i].fwind_ta.n(); m++) {

                  // Write out VL1L2
                  if(conf_info.vx_opt[i].output_flag[i_vl1l2] != STATOutputType_None &&
                     vl1l2_info[m].vcount > 0) {

                     write_vl1l2_row(shc, vl1l2_info[m],
                        conf_info.vx_opt[i].output_flag[i_vl1l2],
                        stat_at, i_stat_row,
                        txt_at[i_vl1l2], i_txt_row[i_vl1l2]);
                  }

                  // Write out VAL1L2
                  if(conf_info.vx_opt[i].output_flag[i_val1l2] != STATOutputType_None &&
                     vl1l2_info[m].vacount > 0) {

                     write_val1l2_row(shc, vl1l2_info[m],
                        conf_info.vx_opt[i].output_flag[i_val1l2],
                        stat_at, i_stat_row,
                        txt_at[i_val1l2], i_txt_row[i_val1l2]);
                  }


                  // Write out VCNT
                  if(conf_info.vx_opt[i].output_flag[i_vcnt] != STATOutputType_None &&
                     vl1l2_info[m].vcount > 0) {

                     write_vcnt_row(shc, vl1l2_info[m],
                        conf_info.vx_opt[i].output_flag[i_vcnt],
                        stat_at, i_stat_row,
                        txt_at[i_vcnt], i_txt_row[i_vcnt]);
                  }

               } // end for m

               // Reset the forecast variable name
               shc.set_fcst_var(conf_info.vx_opt[i].fcst_info->name_attr());

               // Reset the observation variable name
               shc.set_obs_var(conf_info.vx_opt[i].obs_info->name_attr());

            } // end Compute VL1L2

            // Compute PCT counts and scores
            if(conf_info.vx_opt[i].fcst_info->is_prob()                        &&
               (conf_info.vx_opt[i].output_flag[i_pct]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_pstd] != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_pjc]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_prc]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_eclv] != STATOutputType_None)) {
               do_pct(conf_info.vx_opt[i], &pd);
            }

            // Reset the verification masking region
            shc.set_mask(conf_info.vx_opt[i].mask_name[k].c_str());

         } // end for k

         // Set the data smoothing info
         ConcatString mthd = conf_info.vx_opt[i].interp_info.method[j];
         int pnts = interp_pnts;

         // Write out the data fields if requested in the config file
         if(conf_info.vx_opt[i].nc_info.do_raw) {
            write_nc((string)"FCST", fcst_dp_smooth, i, mthd, pnts,
                     conf_info.vx_opt[i].interp_info.field);
            write_nc((string)"OBS",  obs_dp_smooth, i, mthd, pnts,
                     conf_info.vx_opt[i].interp_info.field);
         }
         if(conf_info.vx_opt[i].nc_info.do_diff) {
            write_nc((string)"DIFF", subtract(fcst_dp_smooth, obs_dp_smooth),
                     i, mthd, pnts, conf_info.vx_opt[i].interp_info.field);
         }
         if(conf_info.vx_opt[i].nc_info.do_climo && !cmn_dp.is_empty()) {
            write_nc((string)"CLIMO_MEAN", cmn_dp, i, mthd, pnts,
                     conf_info.vx_opt[i].interp_info.field);
         }
         if(conf_info.vx_opt[i].nc_info.do_climo && !csd_dp.is_empty()) {
            write_nc((string)"CLIMO_STDEV", csd_dp, i, mthd, pnts,
                     conf_info.vx_opt[i].interp_info.field);
         }
         if(conf_info.vx_opt[i].nc_info.do_climo && !cmn_dp.is_empty() && !csd_dp.is_empty()) {
            write_nc((string)"CLIMO_CDF", normal_cdf(obs_dp, cmn_dp, csd_dp),
                     i, mthd, pnts, conf_info.vx_opt[i].interp_info.field);
         }

         // Write out the fields of requested climo distribution percentile threshold values
         if(conf_info.vx_opt[i].nc_info.do_climo_cdp && !cmn_dp.is_empty() && !csd_dp.is_empty()) {

            // Construct one list of all thresholds
            ThreshArray ta;
            ta.add(conf_info.vx_opt[i].fcnt_ta);
            ta.add(conf_info.vx_opt[i].ocnt_ta);
            ta.add(conf_info.vx_opt[i].fcat_ta);
            ta.add(conf_info.vx_opt[i].ocat_ta);

            // Process all the simple thresholds
            vector<Simple_Node> simp;
            ta.get_simple_nodes(simp);

            // Process all CDP thresholds except 0 and 100
            for(vector<Simple_Node>::iterator it = simp.begin();
                it != simp.end(); it++) {
               if(it->ptype() == perc_thresh_climo_dist &&
                  !is_eq(it->pvalue(), 0.0) &&
                  !is_eq(it->pvalue(), 100.0)) {
                  cs << cs_erase << "CLIMO_CDP" << nint(it->pvalue());
                  write_nc(cs, normal_cdf_inv(it->pvalue()/100.0, cmn_dp, csd_dp),
                           i, mthd, pnts, conf_info.vx_opt[i].interp_info.field);
               }
            } // end for it
         }

         // Compute gradient statistics if requested in the config file
         if(!conf_info.vx_opt[i].fcst_info->is_prob() &&
            conf_info.vx_opt[i].output_flag[i_grad] != STATOutputType_None) {

            // Allocate memory in one big chunk based on grid size
            DataPlane fgx_dp, fgy_dp, ogx_dp, ogy_dp;
            PairDataPoint pd_gx, pd_gy;
            pd_gx.extend(grid.nx()*grid.ny());
            pd_gy.extend(grid.nx()*grid.ny());

            // Loop over gradient Dx/Dy
            for(k=0; k<conf_info.vx_opt[i].get_n_grad(); k++) {

               int dx = conf_info.vx_opt[i].grad_dx[k];
               int dy = conf_info.vx_opt[i].grad_dy[k];

               // Compute the gradients
               fgx_dp = gradient(fcst_dp_smooth, 0, dx);
               fgy_dp = gradient(fcst_dp_smooth, 1, dy);
               ogx_dp = gradient(obs_dp_smooth,  0, dx);
               ogy_dp = gradient(obs_dp_smooth,  1, dy);

               // Loop through the masks to be applied
               for(m=0; m<conf_info.vx_opt[i].get_n_mask(); m++) {

                  // Store the current mask
                  mask_mp = conf_info.mask_map[conf_info.vx_opt[i].mask_name[m]];

                  // Turn off the mask for missing data values
                  mask_bad_data(mask_mp, fgx_dp);
                  mask_bad_data(mask_mp, fgy_dp);
                  mask_bad_data(mask_mp, ogx_dp);
                  mask_bad_data(mask_mp, ogy_dp);

                  // Apply the current mask to the current fields
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fgx_dp, &ogx_dp,
                                  0, 0, &wgt_dp, pd_gx);
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fgy_dp, &ogy_dp,
                                  0, 0, &wgt_dp, pd_gy);

                  // Set the mask name
                  shc.set_mask(conf_info.vx_opt[i].mask_name[m].c_str());

                  mlog << Debug(2) << "Computing Gradient DX(" << dx << ")/DY("
                       << dy << ") Statistics " << "over region " << shc.get_mask()
                       << ", using " << pd_gx.f_na.n() << " matched pairs.\n";

                  // Continue if no pairs were found
                  if(pd_gx.f_na.n() == 0) continue;

                  // Compute GRAD
                  grad_info.set(dx, dy, pd_gx.f_na, pd_gy.f_na,
                                pd_gx.o_na, pd_gy.o_na, pd_gx.wgt_na);

                  // Write out GRAD
                  if(conf_info.vx_opt[i].output_flag[i_grad] != STATOutputType_None &&
                     grad_info.total > 0) {

                     write_grad_row(shc, grad_info,
                        conf_info.vx_opt[i].output_flag[i_grad],
                        stat_at, i_stat_row,
                        txt_at[i_grad], i_txt_row[i_grad]);
                  }
               } // end for m (n_mask)

               // Write out the gradients if requested in the config file
               if(conf_info.vx_opt[i].nc_info.do_gradient) {
                  ConcatString cs;
                  cs << cs_erase << "FCST_XGRAD_" << dx;
                  write_nc(cs, fgx_dp, i, mthd, pnts,
                           conf_info.vx_opt[i].interp_info.field);
                  cs << cs_erase << "FCST_YGRAD_" << dy;
                  write_nc(cs, fgy_dp, i, mthd, pnts,
                           conf_info.vx_opt[i].interp_info.field);
                  cs << cs_erase << "OBS_XGRAD_" << dx;
                  write_nc(cs,  ogx_dp, i, mthd, pnts,
                           conf_info.vx_opt[i].interp_info.field);
                  cs << cs_erase << "OBS_YGRAD_" << dy;
                  write_nc(cs,  ogy_dp, i, mthd, pnts,
                           conf_info.vx_opt[i].interp_info.field);
               }

            } // end for k (n_grad)

         } // end if GRAD

         // Compute distance map statistics if requested in the config file
         if(!conf_info.vx_opt[i].fcst_info->is_prob() &&
             conf_info.vx_opt[i].output_flag[i_dmap] != STATOutputType_None) {

            // Allocate memory in one big chunk based on grid size
            DataPlane fcst_dp_dmap, obs_dp_dmap;
            pd.extend(grid.nx()*grid.ny());

            // Loop over the categorical thresholds
            for(k=0; k<conf_info.vx_opt[i].fcat_ta.n(); k++) {

               // Initialize
               dmap_info.clear();
               fcst_dp_dmap.clear();
               obs_dp_dmap.clear();

               // Loop through the masks to be applied
               for(m=0; m<conf_info.vx_opt[i].get_n_mask(); m++) {

                  // Set the mask name
                  shc.set_mask(conf_info.vx_opt[i].mask_name[m].c_str());

                  // Store the current mask
                  mask_mp = conf_info.mask_map[conf_info.vx_opt[i].mask_name[m]];

                  // Process percentile thresholds
                  if(conf_info.vx_opt[i].fcat_ta.need_perc() ||
                     conf_info.vx_opt[i].ocat_ta.need_perc()) {

                     // Apply the current mask
                     get_mask_points(conf_info.vx_opt[i], mask_mp,
                                     &fcst_dp, &obs_dp,
                                     &cmn_dp, 0, 0, pd);

                     // Process percentile thresholds
                     conf_info.vx_opt[i].set_perc_thresh(pd);
                  }

                  // Store the thresholds
                  dmap_info.fthresh = conf_info.vx_opt[i].fcat_ta[k];
                  dmap_info.othresh = conf_info.vx_opt[i].ocat_ta[k];

                  // Compute forecast distance map, if necessary.
                  if(fcst_dp_dmap.is_empty() ||
                     fcst_dp_thresh.is_empty() ||
                     conf_info.vx_opt[i].fcat_ta[k].need_perc()) {

                     fcst_dp_thresh = fcst_dp;
                     fcst_dp_thresh.threshold(conf_info.vx_opt[i].fcat_ta[k]);
                     fcst_dp_dmap = distance_map(fcst_dp_thresh);

                     // Write out the distance map if requested in the config file
                     if(conf_info.vx_opt[i].nc_info.do_distance_map) {
                        ConcatString cs;
                        cs << cs_erase << "FCST_DMAP_"
                           << conf_info.vx_opt[i].fcat_ta[k].get_abbr_str();
                        write_nc(cs, fcst_dp_dmap, i, mthd, pnts,
                                 conf_info.vx_opt[i].interp_info.field);
                     }
                  }

                  // Compute observation distance map, if necessary.
                  if(obs_dp_dmap.is_empty() ||
                     obs_dp_thresh.is_empty() ||
                     conf_info.vx_opt[i].ocat_ta[k].need_perc()) {

                     obs_dp_thresh = obs_dp;
                     obs_dp_thresh.threshold(conf_info.vx_opt[i].ocat_ta[k]);
                     obs_dp_dmap = distance_map(obs_dp_thresh);

                     // Write out the distance map if requested in the config file
                     if(conf_info.vx_opt[i].nc_info.do_distance_map) {
                        ConcatString cs;
                        cs << cs_erase << "OBS_DMAP_"
                           << conf_info.vx_opt[i].ocat_ta[k].get_abbr_str();
                        write_nc(cs, obs_dp_dmap, i, mthd, pnts,
                                 conf_info.vx_opt[i].interp_info.field);
                     }
                  }

                  mlog << Debug(2) << "Computing Distance Map Statistics over region "
                       << shc.get_mask() << ".\n";

                  // Apply the current mask to the distance map and
                  // thresholded fields
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fcst_dp_dmap, &obs_dp_dmap,
                                  0, 0, 0, pd);
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fcst_dp_thresh, &obs_dp_thresh,
                                  0, 0, 0, pd_thr);

                  dmap_info.set_options(
                        conf_info.vx_opt[i].baddeley_p,
                        conf_info.vx_opt[i].baddeley_max_dist,
                        conf_info.vx_opt[i].fom_alpha,
                        conf_info.vx_opt[i].zhu_weight,
                        conf_info.vx_opt[i].beta_value);

                  // Compute DMAP statistics
                  dmap_info.set(conf_info.vx_opt[i].fcat_ta[k],
                                conf_info.vx_opt[i].ocat_ta[k],
                                pd.f_na, pd.o_na, pd_thr.f_na, pd_thr.o_na);

                  // Write out DMAP
                  if(conf_info.vx_opt[i].output_flag[i_dmap] != STATOutputType_None &&
                     dmap_info.total > 0) {

                     write_dmap_row(shc, dmap_info,
                        conf_info.vx_opt[i].output_flag[i_dmap],
                        stat_at, i_stat_row,
                        txt_at[i_dmap], i_txt_row[i_dmap]);
                  }
               } // end for m (n_mask)

            } // end for k

         } // end if DMAP

      } // end for j (n_interp)

      // Loop through and apply the Neighborhood methods for each of the
      // neighborhood widths if requested in the config file
      if(!conf_info.vx_opt[i].fcst_info->is_prob() &&
         (conf_info.vx_opt[i].output_flag[i_nbrctc] != STATOutputType_None ||
          conf_info.vx_opt[i].output_flag[i_nbrcts] != STATOutputType_None ||
          conf_info.vx_opt[i].output_flag[i_nbrcnt] != STATOutputType_None)) {

         // Pointer to current interpolation object
         NbrhdInfo * nbrhd = &conf_info.vx_opt[i].nbrhd_info;

         // Initialize
         for(j=0; j<n_cov; j++) nbrcts_info[j].clear();

         // Loop through and apply each of the neighborhood widths
         for(j=0; j<nbrhd->width.n(); j++) {

            shc.set_interp_mthd(InterpMthd_Nbrhd, nbrhd->shape);
            shc.set_interp_wdth(nbrhd->width[j]);

            // Loop through and apply each of the raw threshold values
            for(k=0; k<conf_info.vx_opt[i].fcat_ta.n(); k++) {

               // Initialize
               fcst_dp_smooth.clear();
               fcst_dp_thresh.clear();
               obs_dp_smooth.clear();
               obs_dp_thresh.clear();

               // Loop through the masks to be applied
               for(m=0; m<conf_info.vx_opt[i].get_n_mask(); m++) {

                  // Set the mask name
                  shc.set_mask(conf_info.vx_opt[i].mask_name[m].c_str());

                  // Store the current mask
                  mask_mp = conf_info.mask_map[conf_info.vx_opt[i].mask_name[m]];

                  // Process percentile thresholds
                  if(conf_info.vx_opt[i].fcat_ta.need_perc() ||
                     conf_info.vx_opt[i].ocat_ta.need_perc()) {

                     // Apply the current mask
                     get_mask_points(conf_info.vx_opt[i], mask_mp,
                                     &fcst_dp, &obs_dp,
                                     &cmn_dp, 0, 0, pd);

                     // Process percentile thresholds
                     conf_info.vx_opt[i].set_perc_thresh(pd);
                  }

                  // Process the forecast data
                  if(nbrhd->field == FieldType_Fcst ||
                     nbrhd->field == FieldType_Both) {

                     // Compute fractional coverage field, if necessary.
                     if(fcst_dp_smooth.is_empty() ||
                        conf_info.vx_opt[i].fcat_ta[k].need_perc()) {

                        // Compute fractional coverage
                        fractional_coverage(fcst_dp, fcst_dp_smooth,
                                            nbrhd->width[j], nbrhd->shape,
                                            conf_info.vx_opt[i].fcat_ta[k],
                                            nbrhd->vld_thresh);

                        // Compute the binary threshold field
                        fcst_dp_thresh = fcst_dp;
                        fcst_dp_thresh.threshold(conf_info.vx_opt[i].fcat_ta[k]);
                     }
                  }
                  // Do not compute fractional coverage
                  else {
                     mlog << Debug(3) << "Skipping forecast fractional "
                          << "coverage computations since \"nbrhd.field\" = "
                          << fieldtype_to_string(nbrhd->field) << "\n";
                     fcst_dp_smooth = fcst_dp;
                     fcst_dp_thresh = fcst_dp;
                     fcst_dp_thresh.set_constant(bad_data_double);

                     // Check range of forecast values
                     fcst_dp_smooth.data_range(dmin, dmax);
                     if(dmin < 0.0 || dmax > 1.0) {
                        mlog << Warning << "\nThe range of forecast "
                             << "fractional coverage values [" << dmin
                             << ", " << dmax << "] falls outside the "
                             << "expected [0, 1] range.\n\n";
                     }
                  }

                  // Process the observation data
                  if(nbrhd->field == FieldType_Obs ||
                     nbrhd->field == FieldType_Both) {

                     // Compute fractional coverage field, if necessary.
                     if(obs_dp_smooth.is_empty() ||
                        conf_info.vx_opt[i].ocat_ta[k].need_perc()) {

                        // Compute fractional coverage
                        fractional_coverage(obs_dp, obs_dp_smooth,
                                            nbrhd->width[j], nbrhd->shape,
                                            conf_info.vx_opt[i].ocat_ta[k],
                                            nbrhd->vld_thresh);

                        // Compute the binary threshold field
                        obs_dp_thresh = obs_dp;
                        obs_dp_thresh.threshold(conf_info.vx_opt[i].ocat_ta[k]);
                     }
                  }
                  // Do not compute fractional coverage
                  else {
                     mlog << Debug(3) << "Skipping observation fractional "
                          << "coverage computations since \"nbrhd.field\" = "
                          << fieldtype_to_string(nbrhd->field) << "\n";
                     obs_dp_smooth = obs_dp;
                     obs_dp_thresh = obs_dp;
                     obs_dp_thresh.set_constant(bad_data_double);

                     // Check range of observation values
                     obs_dp_smooth.data_range(dmin, dmax);
                     if(dmin < 0.0 || dmax > 1.0) {
                        mlog << Warning << "\nThe range of observation "
                             << "fractional coverage values [" << dmin
                             << ", " << dmax << "] falls outside the "
                             << "expected [0, 1] range.\n\n";
                     }
                  }

                  // Turn off the mask for bad forecast or observation values
                  mask_bad_data(mask_mp, fcst_dp_smooth);
                  mask_bad_data(mask_mp, obs_dp_smooth);

                  if(nbrhd->field == FieldType_Fcst ||
                     nbrhd->field == FieldType_Both) {
                     mask_bad_data(mask_mp, fcst_dp_thresh);
                  }

                  if(nbrhd->field == FieldType_Obs ||
                     nbrhd->field == FieldType_Both) {
                     mask_bad_data(mask_mp, obs_dp_thresh);
                  }

                  // Apply the current mask to the fractional coverage
                  // and thresholded fields
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fcst_dp_smooth, &obs_dp_smooth,
                                  0, 0, &wgt_dp, pd);
                  get_mask_points(conf_info.vx_opt[i], mask_mp,
                                  &fcst_dp_thresh, &obs_dp_thresh,
                                  0, 0, 0, pd_thr);

                  // Store climatology values as bad data
                  pd.cmn_na.add_const(bad_data_double, pd.f_na.n());
                  pd.csd_na.add_const(bad_data_double, pd.f_na.n());

                  mlog << Debug(2) << "Processing "
                       << conf_info.vx_opt[i].fcst_info->magic_str()
                       << " versus "
                       << conf_info.vx_opt[i].obs_info->magic_str()
                       << ", for smoothing method "
                       << shc.get_interp_mthd() << "("
                       << shc.get_interp_pnts_str()
                       << "), raw thresholds of "
                       << conf_info.vx_opt[i].fcat_ta[k].get_str()
                       << " and "
                       << conf_info.vx_opt[i].ocat_ta[k].get_str()
                       << ", over region " << shc.get_mask()
                       << ", using " << pd.f_na.n()
                       << " matched pairs.\n";

                  // Continue if no pairs were found
                  if(pd.f_na.n() == 0) continue;

                  // Compute NBRCTS scores
                  if(conf_info.vx_opt[i].output_flag[i_nbrctc] != STATOutputType_None ||
                     conf_info.vx_opt[i].output_flag[i_nbrcts] != STATOutputType_None) {

                     // Initialize
                     for(n=0; n<conf_info.vx_opt[i].get_n_cov_thresh(); n++) {
                        nbrcts_info[n].clear();
                     }

                     do_nbrcts(nbrcts_info, i, j, k, &pd);

                     // Loop through all of the thresholds
                     for(n=0; n<conf_info.vx_opt[i].get_n_cov_thresh(); n++) {

                        // Only write out if n > 0
                        if(nbrcts_info[n].cts_info.cts.n() > 0) {

                           // Write out NBRCTC
                           if(conf_info.vx_opt[i].output_flag[i_nbrctc] != STATOutputType_None) {

                              write_nbrctc_row(shc, nbrcts_info[n],
                                 conf_info.vx_opt[i].output_flag[i_nbrctc],
                                 stat_at, i_stat_row,
                                 txt_at[i_nbrctc], i_txt_row[i_nbrctc]);
                           }

                           // Write out NBRCTS
                           if(conf_info.vx_opt[i].output_flag[i_nbrcts] != STATOutputType_None) {

                              write_nbrcts_row(shc, nbrcts_info[n],
                                 conf_info.vx_opt[i].output_flag[i_nbrcts],
                                 stat_at, i_stat_row,
                                 txt_at[i_nbrcts], i_txt_row[i_nbrcts]);
                           }
                        } // end if
                     } // end for n
                  } // end compute NBRCTS

                  // Compute NBRCNT scores
                  if(conf_info.vx_opt[i].output_flag[i_nbrcnt] != STATOutputType_None) {

                     // Initialize
                     nbrcnt_info.clear();

                     do_nbrcnt(nbrcnt_info, i, j, k, &pd, &pd_thr);

                     // Write out NBRCNT
                     if(nbrcnt_info.sl1l2_info.scount > 0 &&
                        conf_info.vx_opt[i].output_flag[i_nbrcnt]) {

                        write_nbrcnt_row(shc, nbrcnt_info,
                           conf_info.vx_opt[i].output_flag[i_nbrcnt],
                           stat_at, i_stat_row,
                           txt_at[i_nbrcnt], i_txt_row[i_nbrcnt]);
                     }
                  } // end compute NBRCNT

                  // Write out the neighborhood fractional coverage fields
                  // if requested in the config file
                  if(conf_info.vx_opt[i].nc_info.do_nbrhd) {
                     write_nbrhd_nc(fcst_dp_smooth, obs_dp_smooth, i,
                        conf_info.vx_opt[i].fcat_ta[k], conf_info.vx_opt[i].ocat_ta[k],
                        m, conf_info.vx_opt[i].nbrhd_info.width[j]);
                  }
               } // end for m

               // If not computing fractional coverage, all finished with thresholds
               if(nbrhd->field == FieldType_None) break;

            } // end for k

            // If not computing fractional coverage, all finished with neighborhood widths
            if(nbrhd->field == FieldType_None) break;

         } // end for j
      } // end if

      // Loop through and apply each Fourier wave
      for(j=0; j<conf_info.vx_opt[i].get_n_wave_1d(); j++) {

         // Apply Fourier decomposition
         fcst_dp_smooth = fcst_dp;
         obs_dp_smooth  = obs_dp;
         cmn_dp_smooth  = cmn_dp;

         // Reset climo spread since it does not apply to Fourier decomposition
         if(csd_dp.nx() == fcst_dp_smooth.nx() &&
            csd_dp.ny() == fcst_dp_smooth.ny()) {
            csd_dp.set_constant(bad_data_double);
         }

         if(!fcst_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                      conf_info.vx_opt[i].wave_1d_end[j]) ||
             !obs_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                      conf_info.vx_opt[i].wave_1d_end[j])) {
            mlog << Debug(2)
                 << "Skipping Fourier decomposition for waves "
                 << conf_info.vx_opt[i].wave_1d_beg[j] << " to "
                 << conf_info.vx_opt[i].wave_1d_end[j] << " due to the presence "
                 << "of bad data values.\n";
            continue;
         }

         // Decompose the climatology, if provided
         if(cmn_dp_smooth.nx() == fcst_dp_smooth.nx() &&
            cmn_dp_smooth.ny() == fcst_dp_smooth.ny()) {
            if(!cmn_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                        conf_info.vx_opt[i].wave_1d_end[j])) {
               mlog << Debug(2)
                    << "Skipping Fourier decomposition for waves "
                    << conf_info.vx_opt[i].wave_1d_beg[j] << " to "
                    << conf_info.vx_opt[i].wave_1d_end[j] << " due to the presence "
                    << "of bad climatology data values.\n";
               continue;
            }
         }

         // Build string for INTERP_MTHD column
         cs << cs_erase << "WV1_" << conf_info.vx_opt[i].wave_1d_beg[j];
         if(conf_info.vx_opt[i].wave_1d_beg[j] !=
            conf_info.vx_opt[i].wave_1d_end[j]) {
            cs << "-" << conf_info.vx_opt[i].wave_1d_end[j];
         }

         shc.set_interp_mthd(cs, GridTemplateFactory::GridTemplate_None);
         shc.set_interp_pnts(bad_data_int);

         // Loop through the masks to be applied
         for(k=0; k<conf_info.vx_opt[i].get_n_mask(); k++) {

            // Store the current mask
            mask_mp = conf_info.mask_map[conf_info.vx_opt[i].mask_name[k]];

            // Turn off the mask for missing data values
            mask_bad_data(mask_mp, fcst_dp_smooth);
            mask_bad_data(mask_mp, obs_dp_smooth);
            if(cmn_dp_smooth.nx() == fcst_dp_smooth.nx() &&
               cmn_dp_smooth.ny() == fcst_dp_smooth.ny()) {
               mask_bad_data(mask_mp, cmn_dp_smooth);
            }
            if(csd_dp.nx() == fcst_dp_smooth.nx() &&
               csd_dp.ny() == fcst_dp_smooth.ny()) {
               mask_bad_data(mask_mp, csd_dp);
            }

            // Apply the current mask to the current fields
            get_mask_points(conf_info.vx_opt[i], mask_mp,
                            &fcst_dp_smooth, &obs_dp_smooth,
                            &cmn_dp_smooth, &csd_dp, &wgt_dp, pd);

            // Set the mask name
            shc.set_mask(conf_info.vx_opt[i].mask_name[k].c_str());

            mlog << Debug(2)
                 << "Processing " << conf_info.vx_opt[i].fcst_info->magic_str()
                 << " versus " << conf_info.vx_opt[i].obs_info->magic_str()
                 << ", for Fourier wave number " << shc.get_interp_mthd()
                 << ", over region " << shc.get_mask()
                 << ", using " << pd.f_na.n() << " matched pairs.\n";

            // Continue if no pairs were found
            if(pd.f_na.n() == 0) continue;

            // Compute CNT, SL1L2, and SAL1L2 scores
            if(!conf_info.vx_opt[i].fcst_info->is_prob() &&
               (conf_info.vx_opt[i].output_flag[i_cnt]    != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_sl1l2]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_sal1l2] != STATOutputType_None)) {
               do_cnt_sl1l2(conf_info.vx_opt[i], &pd);
            }

            // Compute VL1L2 and VAL1L2 partial sums for UGRD,VGRD
            if(!conf_info.vx_opt[i].fcst_info->is_prob()                         &&
                conf_info.vx_opt[i].fcst_info->is_v_wind()                       &&
                conf_info.vx_opt[i].fcst_info->uv_index() >= 0                   &&
               (conf_info.vx_opt[i].output_flag[i_vl1l2]  != STATOutputType_None ||
                conf_info.vx_opt[i].output_flag[i_val1l2] != STATOutputType_None) ) {

               // Store the forecast variable name
               shc.set_fcst_var(ugrd_vgrd_abbr_str);

               // Store the observation variable name
               shc.set_obs_var(ugrd_vgrd_abbr_str);

               // Initialize
               for(m=0; m<n_wind; m++) vl1l2_info[m].clear();

               // Index for UGRD
               int ui = conf_info.vx_opt[i].fcst_info->uv_index();

               // Read forecast data for UGRD
               if(!read_data_plane(conf_info.vx_opt[ui].fcst_info,
                                   fu_dp, fcst_mtddf, fcst_file)) continue;

               // Read observation data for UGRD
               if(!read_data_plane(conf_info.vx_opt[ui].obs_info,
                                   ou_dp, obs_mtddf, obs_file)) continue;

               // Read climatology data for UGRD
               cmnu_dp = read_climo_data_plane(
                           conf_info.conf.lookup_array(conf_key_climo_mean_field, false),
                           ui, fcst_dp.valid(), grid);

               // Apply Fourier decomposition to the U-wind fields
               fu_dp_smooth   = fu_dp;
               ou_dp_smooth   = ou_dp;
               cmnu_dp_smooth = cmnu_dp;

               if(!fu_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                          conf_info.vx_opt[i].wave_1d_end[j]) ||
                  !ou_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                          conf_info.vx_opt[i].wave_1d_end[j])) {
                  mlog << Debug(2)
                       << "Skipping Fourier decomposition for waves "
                       << conf_info.vx_opt[i].wave_1d_beg[j] << " to "
                       << conf_info.vx_opt[i].wave_1d_end[j] << " due to the "
                       << "presence of bad data values in U-wind.\n";
                  continue;
               }

               // Decompose the U-wind climatology field, if provided
               if(cmnu_dp_smooth.nx() == fu_dp_smooth.nx() &&
                  cmnu_dp_smooth.ny() == fu_dp_smooth.ny()) {
                  if(!cmnu_dp_smooth.fitwav_1d(conf_info.vx_opt[i].wave_1d_beg[j],
                                               conf_info.vx_opt[i].wave_1d_end[j])) {
                     mlog << Debug(2)
                          << "Skipping Fourier decomposition for waves "
                          << conf_info.vx_opt[i].wave_1d_beg[j] << " to "
                          << conf_info.vx_opt[i].wave_1d_end[j] << " due to the presence "
                          << "of bad climatology data values in U-wind.\n";
                     continue;
                  }
               }

               // Apply the current mask to the U-wind fields
               get_mask_points(conf_info.vx_opt[i], mask_mp,
                               &fu_dp_smooth, &ou_dp_smooth,
                               &cmnu_dp_smooth, 0, &wgt_dp, pd_u);

               // Compute VL1L2
               do_vl1l2(vl1l2_info, i, &pd_u, &pd);

               // Loop through all of the wind speed thresholds
               for(m=0; m<conf_info.vx_opt[i].fwind_ta.n(); m++) {

                  // Write out VL1L2
                  if(conf_info.vx_opt[i].output_flag[i_vl1l2] != STATOutputType_None &&
                     vl1l2_info[m].vcount > 0) {

                     write_vl1l2_row(shc, vl1l2_info[m],
                        conf_info.vx_opt[i].output_flag[i_vl1l2],
                        stat_at, i_stat_row,
                        txt_at[i_vl1l2], i_txt_row[i_vl1l2]);
                  }

                  // Write out VAL1L2
                  if(conf_info.vx_opt[i].output_flag[i_val1l2] != STATOutputType_None &&
                     vl1l2_info[m].vacount > 0) {

                     write_val1l2_row(shc, vl1l2_info[m],
                        conf_info.vx_opt[i].output_flag[i_val1l2],
                        stat_at, i_stat_row,
                        txt_at[i_val1l2], i_txt_row[i_val1l2]);
                  }
               } // end for m

               // Reset the forecast variable name
               shc.set_fcst_var(conf_info.vx_opt[i].fcst_info->name_attr());

               // Reset the observation variable name
               shc.set_obs_var(conf_info.vx_opt[i].obs_info->name_attr());

            } // end Compute VL1L2

         } // end for k

         // Write out the data fields if requested in the config file
         if(conf_info.vx_opt[i].nc_info.do_fourier) {

            // Write out the data fields if requested in the config file
            if(conf_info.vx_opt[i].nc_info.do_raw) {
               write_nc((string)"FCST", fcst_dp_smooth, i, shc.get_interp_mthd(),
                        bad_data_int,  FieldType_Both);
               write_nc((string)"OBS",  obs_dp_smooth, i, shc.get_interp_mthd(),
                        bad_data_int, FieldType_Both);
            }
            if(conf_info.vx_opt[i].nc_info.do_diff) {
               write_nc((string)"DIFF", subtract(fcst_dp_smooth, obs_dp_smooth),
                        i, shc.get_interp_mthd(), bad_data_int,
                        FieldType_Both);
            }
            if(conf_info.vx_opt[i].nc_info.do_climo && !cmn_dp_smooth.is_empty()) {
               write_nc((string)"CLIMO_MEAN", cmn_dp_smooth, i, shc.get_interp_mthd(),
                        bad_data_int,  FieldType_Both);
            }
         } // end if

      } // end for j

   } // end for i

   // Check for no data
   if(is_first_pass) {
      mlog << Error << "\nprocess_scores() -> "
           << "no requested data found!  Exiting...\n\n";
      exit(1);
   }

   mlog << Debug(2) << "\n" << sep_str << "\n\n";

   // Deallocate memory
   if(cts_info)    { delete [] cts_info;    cts_info    = (CTSInfo *)    0; }
   if(vl1l2_info)  { delete [] vl1l2_info;  vl1l2_info  = (VL1L2Info *)  0; }
   if(nbrcts_info) { delete [] nbrcts_info; nbrcts_info = (NBRCTSInfo *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void get_mask_points(const GridStatVxOpt &vx_opt,
                     const MaskPlane &mask_mp,
                     const DataPlane *fcst_ptr, const DataPlane *obs_ptr,
                     const DataPlane *cmn_ptr,  const DataPlane *csd_ptr,
                     const DataPlane *wgt_ptr,  PairDataPoint &pd) {

   // Initialize
   pd.erase();

   // Store the climo CDF info
   pd.set_climo_cdf_info(vx_opt.cdf_info);
   
   // Apply the mask the data fields or fill with default values
   apply_mask(*fcst_ptr, mask_mp, pd.f_na);
   apply_mask(*obs_ptr,  mask_mp, pd.o_na);
   pd.n_obs = pd.o_na.n();

   if(cmn_ptr) apply_mask(*cmn_ptr, mask_mp, pd.cmn_na);
   else        pd.cmn_na.add_const(bad_data_double, pd.n_obs);
   if(csd_ptr) apply_mask(*csd_ptr, mask_mp, pd.csd_na);
   else        pd.csd_na.add_const(bad_data_double, pd.n_obs);
   if(wgt_ptr) apply_mask(*wgt_ptr, mask_mp, pd.wgt_na);
   else        pd.wgt_na.add_const(default_grid_weight, pd.n_obs);

   if(cmn_ptr && csd_ptr) pd.add_climo_cdf();

   return;
}

////////////////////////////////////////////////////////////////////////

void do_cts(CTSInfo *&cts_info, int i_vx,
            const PairDataPoint *pd_ptr) {
   int i, j, n_cts;

   //
   // Set up the CTSInfo thresholds and alpha values
   //
   n_cts = conf_info.vx_opt[i_vx].fcat_ta.n();
   for(i=0; i<n_cts; i++) {
      cts_info[i].fthresh = conf_info.vx_opt[i_vx].fcat_ta[i];
      cts_info[i].othresh = conf_info.vx_opt[i_vx].ocat_ta[i];
      cts_info[i].allocate_n_alpha(conf_info.vx_opt[i_vx].get_n_ci_alpha());

      for(j=0; j<conf_info.vx_opt[i_vx].get_n_ci_alpha(); j++) {
         cts_info[i].alpha[j] = conf_info.vx_opt[i_vx].ci_alpha[j];
      }
   }

   mlog << Debug(2) << "Computing Categorical Statistics.\n";

   //
   // Compute the counts, stats, normal confidence intervals, and
   // bootstrap confidence intervals
   //
   if(conf_info.vx_opt[i_vx].boot_info.interval == BootIntervalType_BCA) {
      compute_cts_stats_ci_bca(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         cts_info, n_cts,
         conf_info.vx_opt[i_vx].output_flag[i_cts] != STATOutputType_None,
         conf_info.vx_opt[i_vx].rank_corr_flag, conf_info.tmp_dir.c_str());
   }
   else {
      compute_cts_stats_ci_perc(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         conf_info.vx_opt[i_vx].boot_info.rep_prop,
         cts_info, n_cts,
         conf_info.vx_opt[i_vx].output_flag[i_cts] != STATOutputType_None,
         conf_info.vx_opt[i_vx].rank_corr_flag, conf_info.tmp_dir.c_str());
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void do_mcts(MCTSInfo &mcts_info, int i_vx,
             const PairDataPoint *pd_ptr) {
   int i;

   //
   // Set up the MCTSInfo size, thresholds, and alpha values
   //
   mcts_info.cts.set_size(conf_info.vx_opt[i_vx].fcat_ta.n() + 1);
   mcts_info.cts.set_ec_value(conf_info.vx_opt[i_vx].hss_ec_value);
   mcts_info.set_fthresh(conf_info.vx_opt[i_vx].fcat_ta);
   mcts_info.set_othresh(conf_info.vx_opt[i_vx].ocat_ta);
   mcts_info.allocate_n_alpha(conf_info.vx_opt[i_vx].get_n_ci_alpha());

   for(i=0; i<conf_info.vx_opt[i_vx].get_n_ci_alpha(); i++) {
      mcts_info.alpha[i] = conf_info.vx_opt[i_vx].ci_alpha[i];
   }

   mlog << Debug(2) << "Computing Multi-Category Statistics.\n";

   //
   // If there are no matched pairs to process, return
   //
   if(pd_ptr->f_na.n() == 0 || pd_ptr->o_na.n() == 0) return;

   //
   // Compute the counts, stats, normal confidence intervals, and
   // bootstrap confidence intervals
   //
   if(conf_info.vx_opt[i_vx].boot_info.interval == BootIntervalType_BCA) {
      compute_mcts_stats_ci_bca(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         mcts_info,
         conf_info.vx_opt[i_vx].output_flag[i_mcts] != STATOutputType_None,
         conf_info.vx_opt[i_vx].rank_corr_flag, conf_info.tmp_dir.c_str());
   }
   else {
      compute_mcts_stats_ci_perc(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         conf_info.vx_opt[i_vx].boot_info.rep_prop,
         mcts_info,
         conf_info.vx_opt[i_vx].output_flag[i_mcts] != STATOutputType_None,
         conf_info.vx_opt[i_vx].rank_corr_flag, conf_info.tmp_dir.c_str());
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void do_cnt_sl1l2(const GridStatVxOpt &vx_opt, const PairDataPoint *pd_ptr) {
   int i, j, k, n_bin;
   PairDataPoint pd_thr, pd;
   SL1L2Info *sl1l2_info = (SL1L2Info *) 0;
   CNTInfo   *cnt_info   = (CNTInfo *)   0;

   mlog << Debug(2)
        << "Computing Scalar Partial Sums and Continuous Statistics.\n";

   // Determine the number of climo CDF bins
   n_bin = (pd_ptr->cmn_na.n_valid() > 0 &&
            pd_ptr->csd_na.n_valid() > 0 ?
            vx_opt.get_n_cdf_bin() : 1);

   if(n_bin > 1) {
      mlog << Debug(2)
           << "Applying " << n_bin << " climatology bins.\n";
   }

   // Set flags
   bool do_sl1l2    = (vx_opt.output_flag[i_sl1l2]  != STATOutputType_None ||
                       vx_opt.output_flag[i_sal1l2] != STATOutputType_None);
   bool do_cnt      = (vx_opt.output_flag[i_cnt]    != STATOutputType_None);
   bool precip_flag = (vx_opt.fcst_info->is_precipitation() &&
                       vx_opt.obs_info->is_precipitation());

   // Allocate memory
   cnt_info   = new CNTInfo   [n_bin];
   sl1l2_info = new SL1L2Info [n_bin];

   // Process each continuous filtering threshold
   for(i=0; i<vx_opt.fcnt_ta.n(); i++) {

      // Apply continuous filtering thresholds to subset pairs
      pd_thr = pd_ptr->subset_pairs_cnt_thresh(vx_opt.fcnt_ta[i],
                                               vx_opt.ocnt_ta[i],
                                               vx_opt.cnt_logic);

      // Check for no matched pairs to process
      if(pd_thr.n_obs == 0) continue;

      // Process the climo CDF bins
      for(j=0; j<n_bin; j++) {

         // Initialize
         if(do_sl1l2) sl1l2_info[j].clear();
         if(do_cnt)   cnt_info[j].clear();

         // Apply climo CDF bins logic to subset pairs
         if(n_bin > 1) pd = subset_climo_cdf_bin(pd_thr,
                               vx_opt.cdf_info.cdf_ta, j);
         else          pd = pd_thr;

         // Check for no matched pairs to process
         if(pd.n_obs == 0) continue;

         // Compute and write SL1L2 and SAL1L2 output
         if(do_sl1l2) {

            // Store thresholds
            sl1l2_info[j].fthresh = vx_opt.fcnt_ta[i];
            sl1l2_info[j].othresh = vx_opt.ocnt_ta[i];
            sl1l2_info[j].logic   = vx_opt.cnt_logic;

            // Compute partial sums
            sl1l2_info[j].set(pd);

            // Write out SL1L2
            if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
               vx_opt.output_flag[i_sl1l2] != STATOutputType_None &&
               sl1l2_info[j].scount > 0) {

               write_sl1l2_row(shc, sl1l2_info[j],
                  vx_opt.output_flag[i_sl1l2],
                  j, n_bin, stat_at, i_stat_row,
                  txt_at[i_sl1l2], i_txt_row[i_sl1l2]);
            }

            // Write out SAL1L2
            if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
               vx_opt.output_flag[i_sal1l2] != STATOutputType_None &&
               sl1l2_info[j].sacount > 0) {

               write_sal1l2_row(shc, sl1l2_info[j],
                  vx_opt.output_flag[i_sal1l2],
                  j, n_bin, stat_at, i_stat_row,
                  txt_at[i_sal1l2], i_txt_row[i_sal1l2]);
            }
         } // end do_sl1l2

         // Compute and write CNT output
         if(do_cnt) {

            // Store thresholds
            cnt_info[j].fthresh = vx_opt.fcnt_ta[i];
            cnt_info[j].othresh = vx_opt.ocnt_ta[i];
            cnt_info[j].logic   = vx_opt.cnt_logic;

            // Setup the CNTInfo alpha values
            cnt_info[j].allocate_n_alpha(vx_opt.get_n_ci_alpha());
            for(k=0; k<vx_opt.get_n_ci_alpha(); k++) {
               cnt_info[j].alpha[k] = vx_opt.ci_alpha[k];
            }

            // Compute the stats, normal confidence intervals, and
            // bootstrap confidence intervals
            if(vx_opt.boot_info.interval == BootIntervalType_BCA) {
               compute_cnt_stats_ci_bca(rng_ptr, pd,
                  precip_flag, vx_opt.rank_corr_flag,
                  vx_opt.boot_info.n_rep,
                  cnt_info[j], conf_info.tmp_dir.c_str());
            }
            else {
               compute_cnt_stats_ci_perc(rng_ptr, pd,
                  precip_flag, vx_opt.rank_corr_flag,
                  vx_opt.boot_info.n_rep,
                  vx_opt.boot_info.rep_prop,
                  cnt_info[j], conf_info.tmp_dir.c_str());
            }

            // Write out CNT
            if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
               vx_opt.output_flag[i_cnt] != STATOutputType_None &&
               cnt_info[j].n > 0) {

               write_cnt_row(shc, cnt_info[j],
                  vx_opt.output_flag[i_cnt], j, n_bin,
                  stat_at, i_stat_row, txt_at[i_cnt], i_txt_row[i_cnt]);
            }
         } // end if do_cnt
      } // end for j (n_bin)

      // Write the mean of the climo CDF bins
      if(n_bin > 1) {

         // Compute SL1L2 climo CDF bin means
         if(vx_opt.output_flag[i_sl1l2]  != STATOutputType_None ||
            vx_opt.output_flag[i_sal1l2] != STATOutputType_None) {

            SL1L2Info sl1l2_mean;
            compute_sl1l2_mean(sl1l2_info, n_bin, sl1l2_mean);

            // Write out SL1L2
            if(vx_opt.output_flag[i_sl1l2]  != STATOutputType_None &&
               sl1l2_mean.scount > 0) {

               write_sl1l2_row(shc, sl1l2_mean,
                  vx_opt.output_flag[i_sl1l2],
                  -1, n_bin, stat_at, i_stat_row,
                  txt_at[i_sl1l2], i_txt_row[i_sl1l2]);
            }

            // Write out SAL1L2
            if(vx_opt.output_flag[i_sal1l2] != STATOutputType_None &&
               sl1l2_mean.sacount > 0) {

               write_sal1l2_row(shc, sl1l2_mean,
                  vx_opt.output_flag[i_sal1l2],
                  -1, n_bin, stat_at, i_stat_row,
                  txt_at[i_sal1l2], i_txt_row[i_sal1l2]);
            }
         }

         // Compute CNT climo CDF bin means
         if(vx_opt.output_flag[i_cnt] != STATOutputType_None) {

            CNTInfo cnt_mean;
            compute_cnt_mean(cnt_info, n_bin, cnt_mean);

            if(cnt_mean.n > 0) {

               write_cnt_row(shc, cnt_mean,
                  vx_opt.output_flag[i_cnt],
                  -1, n_bin, stat_at, i_stat_row,
                  txt_at[i_cnt], i_txt_row[i_cnt]);
            }
         }
      } // end if n_bin > 1

   } // end for i (fcnt_ta)

   // Dealloate memory
   if(sl1l2_info) { delete [] sl1l2_info; sl1l2_info = (SL1L2Info *) 0; }
   if(cnt_info)   { delete [] cnt_info;   cnt_info   = (CNTInfo *)   0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void do_vl1l2(VL1L2Info *&v_info, int i_vx,
              const PairDataPoint *pd_u_ptr,
              const PairDataPoint *pd_v_ptr) {
   int i;

   // Check that the number of pairs are the same
   if(pd_u_ptr->n_obs != pd_u_ptr->n_obs) {
      mlog << Error << "\ndo_vl1l2() -> "
           << "unequal number of UGRD and VGRD pairs ("
           << pd_u_ptr->n_obs << " != " << pd_v_ptr->n_obs
           << ")\n\n";
      exit(1);
   }

   //
   // Set all of the VL1L2Info objects
   //
   for(i=0; i<conf_info.vx_opt[i_vx].fwind_ta.n(); i++) {

      //
      // Store thresholds
      //
      v_info[i].zero_out();
      v_info[i].fthresh = conf_info.vx_opt[i_vx].fwind_ta[i];
      v_info[i].othresh = conf_info.vx_opt[i_vx].owind_ta[i];
      v_info[i].logic   = conf_info.vx_opt[i_vx].wind_logic;

      //
      // Compute partial sums
      //
      v_info[i].set(*pd_u_ptr, *pd_v_ptr);

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

void do_pct(const GridStatVxOpt &vx_opt, const PairDataPoint *pd_ptr) {
   int i, j, k, n_bin;
   PairDataPoint pd;
   PCTInfo *pct_info = (PCTInfo *) 0;

   mlog << Debug(2)
        << "Computing Probabilistic Statistics.\n";

   // Determine the number of climo CDF bins
   n_bin = (pd_ptr->cmn_na.n_valid() > 0 && pd_ptr->csd_na.n_valid() > 0 ?
            vx_opt.get_n_cdf_bin() : 1);

   if(n_bin > 1) {
      mlog << Debug(2)
           << "Applying " << n_bin << " climatology bins.\n";
   }

   // Allocate memory
   pct_info = new PCTInfo [n_bin];

   // Process each probabilistic observation threshold
   for(i=0; i<vx_opt.ocat_ta.n(); i++) {

      // Process the climo CDF bins
      for(j=0; j<n_bin; j++) {

         // Initialize
         pct_info[j].clear();

         // Apply climo CDF bins logic to subset pairs
         if(n_bin > 1) pd = subset_climo_cdf_bin(*pd_ptr,
                               vx_opt.cdf_info.cdf_ta, j);
         else          pd = *pd_ptr;

         // Store thresholds
         pct_info[j].fthresh = vx_opt.fcat_ta;
         pct_info[j].othresh = vx_opt.ocat_ta[i];
         pct_info[j].allocate_n_alpha(vx_opt.get_n_ci_alpha());

         for(k=0; k<vx_opt.get_n_ci_alpha(); k++) {
            pct_info[j].alpha[k] = vx_opt.ci_alpha[k];
         }

         // Compute the probabilistic counts and statistics
         compute_pctinfo(pd, vx_opt.output_flag[i_pstd], pct_info[j]);

         // Check for no matched pairs to process
         if(pd.n_obs == 0) continue;

         // Write out PCT
         if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
            vx_opt.output_flag[i_pct] != STATOutputType_None) {
            write_pct_row(shc, pct_info[j],
               vx_opt.output_flag[i_pct],
               j, n_bin, stat_at, i_stat_row,
               txt_at[i_pct], i_txt_row[i_pct]);
         }

         // Write out PSTD
         if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
            vx_opt.output_flag[i_pstd] != STATOutputType_None) {
            write_pstd_row(shc, pct_info[j],
               vx_opt.output_flag[i_pstd],
               j, n_bin, stat_at, i_stat_row,
               txt_at[i_pstd], i_txt_row[i_pstd]);
         }

         // Write out PJC
         if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
            vx_opt.output_flag[i_pjc] != STATOutputType_None) {
            write_pjc_row(shc, pct_info[j],
               vx_opt.output_flag[i_pjc],
               j, n_bin, stat_at, i_stat_row,
               txt_at[i_pjc], i_txt_row[i_pjc]);
         }

         // Write out PRC
         if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
            vx_opt.output_flag[i_prc] != STATOutputType_None) {
            write_prc_row(shc, pct_info[j],
               vx_opt.output_flag[i_prc],
               j, n_bin, stat_at, i_stat_row,
               txt_at[i_prc], i_txt_row[i_prc]);
         }

         // Write out ECLV
         if((n_bin == 1 || vx_opt.cdf_info.write_bins) &&
            vx_opt.output_flag[i_eclv] != STATOutputType_None) {
            write_eclv_row(shc, pct_info[j], vx_opt.eclv_points,
               vx_opt.output_flag[i_eclv],
               j, n_bin, stat_at, i_stat_row,
               txt_at[i_eclv], i_txt_row[i_eclv]);
         }
      } // end for j (n_bin)

      // Write the mean of the climo CDF bins
      if(n_bin > 1) {

         PCTInfo pct_mean;
         compute_pct_mean(pct_info, n_bin, pct_mean);

         // Write out PSTD
         if(vx_opt.output_flag[i_pstd] != STATOutputType_None) {
            write_pstd_row(shc, pct_mean,
               vx_opt.output_flag[i_pstd],
               -1, n_bin, stat_at, i_stat_row,
               txt_at[i_pstd], i_txt_row[i_pstd]);
         }
      } // end if n_bin > 1

   } // end for i (ocnt_ta)

   // Dealloate memory
   if(pct_info) { delete [] pct_info; pct_info = (PCTInfo *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void do_nbrcts(NBRCTSInfo *&nbrcts_info,
               int i_vx, int i_wdth, int i_thresh,
               const PairDataPoint *pd_ptr) {
   int i, j, n_nbrcts;
   FieldType field = conf_info.vx_opt[i_vx].nbrhd_info.field;

   //
   // Set up the NBRCTSInfo thresholds and alpha values
   //
   n_nbrcts = conf_info.vx_opt[i_vx].get_n_cov_thresh();
   for(i=0; i<n_nbrcts; i++) {

      if(field == FieldType_Fcst || field == FieldType_Both) {
         nbrcts_info[i].fthresh =
            conf_info.vx_opt[i_vx].fcat_ta[i_thresh];
      }
      else {
         nbrcts_info[i].fthresh.set_na();
      }

      if(field == FieldType_Obs || field == FieldType_Both) {
         nbrcts_info[i].othresh =
            conf_info.vx_opt[i_vx].ocat_ta[i_thresh];
      }
      else {
         nbrcts_info[i].othresh.set_na();
      }
      nbrcts_info[i].cthresh =
         conf_info.vx_opt[i_vx].nbrhd_info.cov_ta[i];

      nbrcts_info[i].cts_info.fthresh =
         conf_info.vx_opt[i_vx].nbrhd_info.cov_ta[i];
      nbrcts_info[i].cts_info.othresh =
         conf_info.vx_opt[i_vx].nbrhd_info.cov_ta[i];
      nbrcts_info[i].allocate_n_alpha(conf_info.vx_opt[i_vx].get_n_ci_alpha());

      for(j=0; j<conf_info.vx_opt[i_vx].get_n_ci_alpha(); j++) {
         nbrcts_info[i].cts_info.alpha[j] = conf_info.vx_opt[i_vx].ci_alpha[j];
      }
   }

   //
   // Keep track of the neighborhood width
   //
   for(i=0; i<n_nbrcts; i++) {
      nbrcts_info[i].nbr_wdth = conf_info.vx_opt[i_vx].nbrhd_info.width[i_wdth];
   }

   mlog << Debug(2) << "Computing Neighborhood Categorical Statistics.\n";

   //
   // Compute the stats, normal confidence intervals, and
   // bootstrap confidence intervals
   //
   if(conf_info.vx_opt[i_vx].boot_info.interval == BootIntervalType_BCA) {
      compute_nbrcts_stats_ci_bca(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         nbrcts_info, n_nbrcts,
         conf_info.vx_opt[i_vx].output_flag[i_nbrcts] != STATOutputType_None,
         conf_info.tmp_dir.c_str());
   }
   else {
      compute_nbrcts_stats_ci_perc(rng_ptr, *pd_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         conf_info.vx_opt[i_vx].boot_info.rep_prop,
         nbrcts_info, n_nbrcts,
         conf_info.vx_opt[i_vx].output_flag[i_nbrcts] != STATOutputType_None,
         conf_info.tmp_dir.c_str());
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void do_nbrcnt(NBRCNTInfo &nbrcnt_info,
               int i_vx, int i_wdth, int i_thresh,
               const PairDataPoint *pd_ptr,
               const PairDataPoint *pd_thr_ptr) {
   int i;
   FieldType field = conf_info.vx_opt[i_vx].nbrhd_info.field;

   //
   // Set up the NBRCNTInfo threshold and alpha values
   //
   if(field == FieldType_Fcst || field == FieldType_Both) {
      nbrcnt_info.fthresh = conf_info.vx_opt[i_vx].fcat_ta[i_thresh];
   }
   else {
      nbrcnt_info.fthresh.set_na();
   }

   if(field == FieldType_Obs || field == FieldType_Both) {
      nbrcnt_info.othresh = conf_info.vx_opt[i_vx].ocat_ta[i_thresh];
   }
   else {
      nbrcnt_info.othresh.set_na();
   }

   nbrcnt_info.allocate_n_alpha(conf_info.vx_opt[i_vx].get_n_ci_alpha());
   for(i=0; i<conf_info.vx_opt[i_vx].get_n_ci_alpha(); i++) {
      nbrcnt_info.alpha[i] = conf_info.vx_opt[i_vx].ci_alpha[i];
   }

   //
   // Keep track of the neighborhood width
   //
   nbrcnt_info.nbr_wdth = conf_info.vx_opt[i_vx].nbrhd_info.width[i_wdth];

   mlog << Debug(2) << "Computing Neighborhood Continuous Statistics.\n";

   //
   // Compute the stats, normal confidence intervals, and
   // bootstrap confidence intervals
   //
   if(conf_info.vx_opt[i_vx].boot_info.interval == BootIntervalType_BCA) {
      compute_nbrcnt_stats_ci_bca(rng_ptr, *pd_ptr, *pd_thr_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         nbrcnt_info,
         conf_info.vx_opt[i_vx].output_flag[i_nbrcnt] != STATOutputType_None,
         conf_info.tmp_dir.c_str());
   }
   else {
      compute_nbrcnt_stats_ci_perc(rng_ptr, *pd_ptr, *pd_thr_ptr,
         conf_info.vx_opt[i_vx].boot_info.n_rep,
         conf_info.vx_opt[i_vx].boot_info.rep_prop,
         nbrcnt_info,
         conf_info.vx_opt[i_vx].output_flag[i_nbrcnt] != STATOutputType_None,
         conf_info.tmp_dir.c_str());
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nc(const ConcatString &field_name, const DataPlane &dp,
              int i_vx, const ConcatString &interp_mthd,
              int interp_pnts, FieldType field_type) {
   int i, x, y, n, n_masks;
   ConcatString var_name, var_suffix, interp_str, mask_str;
   ConcatString fcst_name, obs_name, fcst_obs_name;
   ConcatString fcst_long_name, obs_long_name;
   ConcatString long_att, level_att, units_att;
   NcVar nc_var;
   bool apply_mask;

   // Set output forecast and observation names
   if(conf_info.vx_opt[i_vx].var_name.length() > 0) {
      fcst_name = obs_name = fcst_obs_name = conf_info.vx_opt[i_vx].var_name;
   }
   else {
      fcst_name     << conf_info.vx_opt[i_vx].fcst_info->name_attr() << "_"
                    << conf_info.vx_opt[i_vx].fcst_info->level_attr();
      obs_name      << conf_info.vx_opt[i_vx].obs_info->name_attr() << "_"
                    << conf_info.vx_opt[i_vx].obs_info->level_attr();
      fcst_obs_name << fcst_name << "_" << obs_name;
   }

   // Set output forecast and observation long names
   fcst_long_name = (conf_info.vx_opt[i_vx].fcst_info->long_name_attr().empty() ?
                     conf_info.vx_opt[i_vx].fcst_info->name_attr() :
                     conf_info.vx_opt[i_vx].fcst_info->long_name_attr());
   fcst_long_name << " at " << conf_info.vx_opt[i_vx].fcst_info->level_attr();

   obs_long_name  = (conf_info.vx_opt[i_vx].obs_info->long_name_attr().empty() ?
                     conf_info.vx_opt[i_vx].obs_info->name_attr() :
                     conf_info.vx_opt[i_vx].obs_info->long_name_attr());
   obs_long_name  << " at " << conf_info.vx_opt[i_vx].obs_info->level_attr();

   // Append nc_pairs_var_suffix config file entry
   if(conf_info.vx_opt[i_vx].var_suffix.length() > 0) {
      var_suffix << "_" << conf_info.vx_opt[i_vx].var_suffix;
   }

   // Append smoothing info for all but nearest neighbor
   if(interp_pnts > 1 ||
      interp_mthd == interpmthd_to_string(InterpMthd_Gaussian)) {
      interp_str << "_" << interp_mthd << "_" << interp_pnts;
   }
   // Append Fourier decomposition info
   else if(is_bad_data(interp_pnts)) {
      interp_str << "_" << interp_mthd;
   }

   // Determine the number of masking regions
   // MET-621: When the nc_pairs flag 'apply_mask'= FALSE in the config
   // file, generate the output NetCDF file for the full domain only.
   // The default behavior is to generate fields for each masking region.
   apply_mask = conf_info.vx_opt[i_vx].nc_info.do_apply_mask;
   n_masks    = (apply_mask ? conf_info.vx_opt[i_vx].get_n_mask() : 1);

   // Allocate memory
   float *data = (float *) 0;
   data = new float [grid.nx()*grid.ny()];

   // Set the NetCDF compression level
   int deflate_level = compress_level;
   if (deflate_level < 0) deflate_level = conf_info.get_compression_level();

   // Process each of the masking regions
   for(i=0; i<n_masks; i++) {

      // If apply_mask is true, create fields for each masking region.
      // Otherwise create only fields for the FULL domain.
      mask_str = (apply_mask ? conf_info.vx_opt[i_vx].mask_name[i] : "FULL");

      // Build the NetCDF variable name
      if(field_name == "FCST") {
         var_name  << cs_erase << field_name << "_"
                   << fcst_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Fcst ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  = fcst_long_name;
         level_att = shc.get_fcst_lev();
         units_att = conf_info.vx_opt[i_vx].fcst_info->units_attr();
      }
      else if(field_name == "OBS") {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Obs ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  = obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(field_name == "DIFF") {
         var_name  << cs_erase << field_name << "_"
                   << fcst_obs_name << var_suffix << "_" << mask_str << interp_str;
         long_att  << cs_erase
                   << fcst_long_name << " and "
                   << obs_long_name;
         level_att << cs_erase
                   << shc.get_fcst_lev() << " and "
                   << shc.get_obs_lev();
         units_att << cs_erase
                   << conf_info.vx_opt[i_vx].fcst_info->units_attr() << " and "
                   << conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(field_name == "CLIMO_MEAN") {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         // Append interpolation string for Fourier decomposition
         if(interp_str.nonempty()) {
            if(strncmp(interp_str.c_str(), "_WV", 3) == 0) var_name << interp_str;
         }
         long_att  << cs_erase
                   << "Climatology mean for "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(field_name == "CLIMO_STDEV") {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         long_att  << cs_erase
                   << "Climatology standard deviation for "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(field_name == "CLIMO_CDF") {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         long_att  << cs_erase
                   << "Climatology cumulative distribution function for "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(strncmp(field_name.c_str(), "CLIMO_CDP", 9) == 0) {
         var_name  << cs_erase
                   << field_name << "_"
                   << conf_info.vx_opt[i_vx].obs_info->name_attr() << "_"
                   << conf_info.vx_opt[i_vx].obs_info->level_attr()
                   << var_suffix << "_" << mask_str;
         long_att  << cs_erase
                   << "Climatology distribution percentile thresholds for "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(check_reg_exp("FCST_XGRAD_", field_name.c_str()) ||
              check_reg_exp("FCST_YGRAD_", field_name.c_str())) {
         var_name  << cs_erase << field_name << "_"
                   << fcst_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Fcst ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  << cs_erase
                   << "Gradient of "
                   << fcst_long_name;
         level_att = shc.get_fcst_lev();
         units_att = conf_info.vx_opt[i_vx].fcst_info->units_attr();
      }
      else if(check_reg_exp("OBS_XGRAD_", field_name.c_str()) ||
              check_reg_exp("OBS_YGRAD_", field_name.c_str())) {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Obs ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  << cs_erase
                   << "Gradient of "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else if(check_reg_exp("FCST_DMAP_", field_name.c_str())) {
         var_name  << cs_erase << field_name << "_"
                   << fcst_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Fcst ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  << cs_erase
                   << "Distance Map for "
                   << fcst_long_name;
         level_att = shc.get_fcst_lev();
         units_att = conf_info.vx_opt[i_vx].fcst_info->units_attr();
      }
      else if(check_reg_exp("OBS_DMAP_", field_name.c_str())) {
         var_name  << cs_erase << field_name << "_"
                   << obs_name << var_suffix << "_" << mask_str;
         if(field_type == FieldType_Obs ||
            field_type == FieldType_Both) {
            var_name << interp_str;
         }
         long_att  << cs_erase
                   << "Distance Map for "
                   << obs_long_name;
         level_att = shc.get_obs_lev();
         units_att = conf_info.vx_opt[i_vx].obs_info->units_attr();
      }
      else {
         mlog << Error << "\nwrite_nc() -> "
              << "unexpected field name \"" << field_name << "\"\n\n";
         exit(1);
      }

      // Skip differences for probability forecasts
      if(field_name == "DIFF" &&
         conf_info.vx_opt[i_vx].fcst_info->is_prob()) continue;

      // Skip variable names that have already been written
      if(nc_var_sa.has(var_name)) continue;

      // Otherwise, add to the list of previously defined variables
      nc_var_sa.add(var_name);

      // Define the variable
      nc_var = add_var(nc_out, (string) var_name,
                       ncFloat, lat_dim, lon_dim, deflate_level);

      // Add variable attributes
      add_var_att_local(&nc_var, "name", nc_var.getName());
      add_var_att_local(&nc_var, "long_name", long_att);
      add_var_att_local(&nc_var, "level", level_att);
      add_var_att_local(&nc_var, "units", units_att);
      write_netcdf_var_times(&nc_var, dp);
      add_att(&nc_var, "_FillValue", bad_data_float);
      add_var_att_local(&nc_var, "desc", conf_info.vx_opt[i_vx].desc);
      add_var_att_local(&nc_var, "masking_region", mask_str);
      add_var_att_local(&nc_var, "smoothing_method", interp_mthd);
      add_att(&nc_var, "smoothing_neighborhood", interp_pnts);

      // Pointer to the current mask
      MaskPlane * mask_ptr =
         &conf_info.mask_map[conf_info.vx_opt[i_vx].mask_name[i]];

      // Store the data
      for(x=0; x<grid.nx(); x++) {
         for(y=0; y<grid.ny(); y++) {

            n = DefaultTO.two_to_one(grid.nx(), grid.ny(), x, y);

            // Check apply_mask
            if(!apply_mask ||
               (apply_mask && mask_ptr->s_is_on(x, y))) {
               data[n] = dp.get(x, y);
            }
            else {
               data[n] = bad_data_float;
            }

         } // end for y
      } // end for x

      // Write out the data
      if(!put_nc_data_with_dims(&nc_var, &data[0], grid.ny(), grid.nx())) {
         mlog << Error << "\nwrite_nc() -> "
              << "error writing NetCDF variable name " << var_name
              << "\n\n";
         exit(1);
      }

   } // end for i

   // Deallocate and clean up
   if(data) { delete [] data; data = (float *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrhd_nc(const DataPlane &fcst_dp, const DataPlane &obs_dp,
                    int i_vx, const SingleThresh &fcst_st,
                    const SingleThresh &obs_st, int i_mask, int wdth) {
   int n, x, y;
   int fcst_flag, obs_flag;
   ConcatString fcst_var_name, obs_var_name, var_suffix, mask_str;
   ConcatString fcst_name, obs_name, att_str, mthd_str, nbrhd_str;
   bool apply_mask;

   // Set output forecast and observation names
   if(conf_info.vx_opt[i_vx].var_name.length() > 0) {
      fcst_name = obs_name = conf_info.vx_opt[i_vx].var_name;
   }
   else {
      fcst_name     << conf_info.vx_opt[i_vx].fcst_info->name_attr() << "_"
                    << conf_info.vx_opt[i_vx].fcst_info->level_attr();
      obs_name      << conf_info.vx_opt[i_vx].obs_info->name_attr() << "_"
                    << conf_info.vx_opt[i_vx].obs_info->level_attr();
   }

   // Append nc_pairs_var_suffix config file entry
   if(conf_info.vx_opt[i_vx].var_suffix.length() > 0) {
      var_suffix << "_" << conf_info.vx_opt[i_vx].var_suffix;
   }

   // Store the apply_mask option
   apply_mask = conf_info.vx_opt[i_vx].nc_info.do_apply_mask;

   float *fcst_data = (float *) 0;
   float *obs_data  = (float *) 0;

   NcVar fcst_var;
   NcVar obs_var;

   // Get the interpolation strings
   mthd_str = interpmthd_to_string(InterpMthd_Nbrhd);
   if(wdth > 1) nbrhd_str << "_" << mthd_str << "_" << wdth*wdth;

   int deflate_level = compress_level;
   if (deflate_level < 0) deflate_level = conf_info.get_compression_level();

   // If apply_mask is true, create fields for each masking region.
   // Otherwise create only fields for the FULL domain.
   mask_str = (apply_mask ? conf_info.vx_opt[i_vx].mask_name[i_mask] : "FULL");

   // Build the forecast variable name
   fcst_var_name << cs_erase << "FCST_"
                 << fcst_name << var_suffix << "_" << mask_str << "_"
                 << fcst_st.get_abbr_str() << nbrhd_str;

   // Build the observation variable name
   obs_var_name << cs_erase << "OBS_"
                << obs_name << var_suffix << "_" << mask_str << "_"
                << obs_st.get_abbr_str() << nbrhd_str;

   // Figure out which fields should be written
   fcst_flag = !nc_var_sa.has(fcst_var_name);
   obs_flag  = !nc_var_sa.has(obs_var_name);

   // Check for nothing to do
   if(!fcst_flag && !obs_flag) return;

   // Allocate memory for the forecast and observation fields
   fcst_data = new float [grid.nx()*grid.ny()];
   obs_data  = new float [grid.nx()*grid.ny()];

   // Add the forecast variable
   if(fcst_flag) {

      // Define the forecast variable
      fcst_var = add_var(nc_out, (string) fcst_var_name, ncFloat,
                         lat_dim, lon_dim, deflate_level);

      // Add to the list of previously defined variables
      nc_var_sa.add(fcst_var_name);

      // Add variable attributes for the forecast field
      add_var_att_local(&fcst_var, "name", shc.get_fcst_var());
      att_str << cs_erase << conf_info.vx_opt[i_vx].fcst_info->name_attr()
              << " at "
              << conf_info.vx_opt[i_vx].fcst_info->level_attr();
      add_var_att_local(&fcst_var, "long_name", att_str);
      add_var_att_local(&fcst_var, "level", shc.get_fcst_lev());
      add_var_att_local(&fcst_var, "units", conf_info.vx_opt[i_vx].fcst_info->units_attr());
      write_netcdf_var_times(&fcst_var, fcst_dp);
      add_att(&fcst_var, "_FillValue", bad_data_float);
      add_var_att_local(&fcst_var, "masking_region", mask_str);
      add_var_att_local(&fcst_var, "smoothing_method", mthd_str);
      add_var_att_local(&fcst_var, "threshold", fcst_st.get_abbr_str());
      add_att(&fcst_var, "smoothing_neighborhood", wdth*wdth);
   } // end fcst_flag

   // Add the observation variable
   if(obs_flag) {

      // Define the observation variable
      obs_var  = add_var(nc_out, (string)obs_var_name,  ncFloat,
                            lat_dim, lon_dim, deflate_level);

      // Add to the list of previously defined variables
      nc_var_sa.add(obs_var_name);

      // Add variable attributes for the observation field
      add_var_att_local(&obs_var, "name", obs_var.getName());
      att_str << cs_erase
              << conf_info.vx_opt[i_vx].obs_info->name_attr() << " at "
              << conf_info.vx_opt[i_vx].obs_info->level_attr();
      add_var_att_local(&obs_var, "long_name", att_str);
      add_var_att_local(&obs_var, "level", shc.get_obs_lev());
      add_var_att_local(&obs_var, "units", conf_info.vx_opt[i_vx].obs_info->units_attr());
      write_netcdf_var_times(&obs_var, obs_dp);
      add_att(&obs_var, "_FillValue", bad_data_float);
      add_var_att_local(&obs_var, "masking_region", mask_str);
      add_var_att_local(&obs_var, "smoothing_method", mthd_str);
      add_var_att_local(&obs_var, "threshold", obs_st.get_abbr_str());
      add_att(&obs_var, "smoothing_neighborhood", wdth*wdth);
   } // end obs_flag

   // Pointer to the current mask
   MaskPlane * mask_ptr =
      &conf_info.mask_map[conf_info.vx_opt[i_vx].mask_name[i_mask]];

   // Store the forecast and observation values
   for(x=0; x<grid.nx(); x++) {
      for(y=0; y<grid.ny(); y++) {

         n = DefaultTO.two_to_one(grid.nx(), grid.ny(), x, y);

         // Check apply_mask
         if(!apply_mask ||
            (apply_mask && mask_ptr->s_is_on(x, y))) {
             fcst_data[n] = fcst_dp.get(x, y);
              obs_data[n] =  obs_dp.get(x, y);
         }
         else {
            fcst_data[n] = bad_data_float;
             obs_data[n] = bad_data_float;
         }
      } // end for y
   } // end for x

   // Write out the forecast field
   if(fcst_flag) {
      if(!put_nc_data_with_dims(&fcst_var, &fcst_data[0], grid.ny(), grid.nx())) {
         mlog << Error << "\nwrite_nbrhd_nc() -> "
              << "error with the fcst_var->put for forecast variable "
              << fcst_var_name << "\n\n";
         exit(1);
      }
   }

   // Write out the observation field
   if(obs_flag) {
      if(!put_nc_data_with_dims(&obs_var, &obs_data[0], grid.ny(), grid.nx())) {
         mlog << Error << "\nwrite_nbrhd_nc() -> "
              << "error with the obs_var->put for observation variable "
              << obs_var_name << "\n\n";
         exit(1);
      }
   }

   // Deallocate and clean up
   if(fcst_data) { delete [] fcst_data; fcst_data = (float *) 0; }
   if(obs_data)  { delete [] obs_data;  obs_data  = (float *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

void add_var_att_local(NcVar *var, const char *att_name,
                       const ConcatString att_value) {

   if(att_value.nonempty()) add_att(var, att_name, att_value.c_str());
   else                     add_att(var, att_name, na_str);

   return;
}

////////////////////////////////////////////////////////////////////////

void finish_txt_files() {
   int i;

   // Write out the contents of the STAT AsciiTable and
   // close the STAT output files
   if(stat_out) {
      *stat_out << stat_at;
      close_txt_file(stat_out, stat_file.c_str());
   }

   // Finish up each of the optional text files
   for(i=0; i<n_txt; i++) {

      // Only write the table if requested in the config file
      if(conf_info.output_flag[i] == STATOutputType_Both) {

         // Write the AsciiTable to a file
         if(txt_out[i]) {
            *txt_out[i] << txt_at[i];
            close_txt_file(txt_out[i], txt_file[i].c_str());
         }
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void clean_up() {

   // Close the output text files that were open for writing
   finish_txt_files();

   // Close the output NetCDF file
   if(nc_out) {

      // List the NetCDF file after it is finished
      mlog << Debug(1) << "Output file: " << out_nc_file << "\n";

      //nc_out->close();
      delete nc_out;
      nc_out = (NcFile *) 0;
   }

   // Deallocate memory for data files
   if(fcst_mtddf) { delete fcst_mtddf; fcst_mtddf = (Met2dDataFile *) 0; }
   if(obs_mtddf)  { delete obs_mtddf;  obs_mtddf  = (Met2dDataFile *) 0; }

   // Deallocate memory for the random number generator
   rng_free(rng_ptr);

   return;
}

////////////////////////////////////////////////////////////////////////

void usage() {

   cout << "\n*** Model Evaluation Tools (MET" << met_version
        << ") ***\n\n"

        << "Usage: " << program_name << "\n"
        << "\tfcst_file\n"
        << "\tobs_file\n"
        << "\tconfig_file\n"
        << "\t[-outdir path]\n"
        << "\t[-log file]\n"
        << "\t[-v level]\n"
        << "\t[-compress level]\n\n"

        << "\twhere\t\"fcst_file\" is a gridded forecast file containing "
        << "the field(s) to be verified (required).\n"

        << "\t\t\"obs_file\" is a gridded observation file containing "
        << "the verifying field(s) (required).\n"

        << "\t\t\"config_file\" is a GridStatConfig file containing "
        << "the desired configuration settings (required).\n"

        << "\t\t\"-outdir path\" overrides the default output directory "
        << "(" << out_dir << ") (optional).\n"

        << "\t\t\"-log file\" outputs log messages to the specified "
        << "file (optional).\n"

        << "\t\t\"-v level\" overrides the default level of logging ("
        << mlog.verbosity_level() << ") (optional).\n"

        << "\t\t\"-compress level\" overrides the compression level of NetCDF variable ("
        << conf_info.get_compression_level() << ") (optional).\n\n" << flush;

   exit(1);
}

////////////////////////////////////////////////////////////////////////

void set_outdir(const StringArray & a) {
   out_dir = a[0];
}

////////////////////////////////////////////////////////////////////////

void set_compress(const StringArray & a) {
   compress_level = atoi(a[0].c_str());
}

////////////////////////////////////////////////////////////////////////

bool read_data_plane(VarInfo* info, DataPlane& dp, Met2dDataFile* mtddf,
                     const ConcatString &filename) {

   bool status = mtddf->data_plane(*info, dp);

   if(!status) {
      mlog << Warning << "\nread_data_plane() -> "
           << info->magic_str()
           << " not found in file: " << filename
           << "\n\n";
      return false;
   }

   // Regrid, if necessary
   if(!(mtddf->grid() == grid)) {
      mlog << Debug(1)
           << "Regridding field "
           << info->magic_str()
           << " to the verification grid.\n";
      dp = met_regrid(dp, mtddf->grid(), grid, info->regrid());
   }

   // Rescale probabilities from [0, 100] to [0, 1]
   if(info->p_flag()) rescale_probability(dp);

   return(true);
}

////////////////////////////////////////////////////////////////////////
