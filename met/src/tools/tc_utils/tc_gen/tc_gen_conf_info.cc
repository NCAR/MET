// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2020
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

#include "tc_gen_conf_info.h"

#include "vx_nc_util.h"
#include "apply_mask.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////
//
//  Code for struct TCGenNcOutInfo
//
////////////////////////////////////////////////////////////////////////

TCGenNcOutInfo::TCGenNcOutInfo() {
   clear();
}

////////////////////////////////////////////////////////////////////////

TCGenNcOutInfo & TCGenNcOutInfo::operator+=(const TCGenNcOutInfo &t) {

   if(t.do_latlon)       do_latlon       = true;
   if(t.do_fcst_genesis) do_fcst_genesis = true;
   if(t.do_fcst_fy_oy)   do_fcst_fy_oy   = true;
   if(t.do_fcst_fy_on)   do_fcst_fy_on   = true;
   if(t.do_fcst_tracks)  do_fcst_tracks  = true;
   if(t.do_best_genesis) do_best_genesis = true;
   if(t.do_best_fy_oy)   do_best_fy_oy   = true;
   if(t.do_best_fn_oy)   do_best_fn_oy   = true;
   if(t.do_best_tracks)  do_best_tracks  = true;

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void TCGenNcOutInfo::clear() {

   set_all_false();

   return;
}

////////////////////////////////////////////////////////////////////////

bool TCGenNcOutInfo::all_false() const {

   bool status = do_latlon       ||
                 do_fcst_genesis || do_fcst_fy_oy  ||
                 do_fcst_fy_on   || do_fcst_tracks ||
                 do_best_genesis || do_best_fy_oy  ||
                 do_best_fn_oy   || do_best_tracks;

   return(!status);
}

////////////////////////////////////////////////////////////////////////

void TCGenNcOutInfo::set_all_false() {

   do_latlon       = false;
   do_fcst_genesis = false;
   do_fcst_fy_oy   = false;
   do_fcst_fy_on   = false;
   do_fcst_tracks  = false;
   do_best_genesis = false;
   do_best_fy_oy   = false;
   do_best_fn_oy   = false;
   do_best_tracks  = false;

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenNcOutInfo::set_all_true() {

   do_latlon       = true;
   do_fcst_genesis = true;
   do_fcst_fy_oy   = true;
   do_fcst_fy_on   = true;
   do_fcst_tracks  = true;
   do_best_genesis = true;
   do_best_fy_oy   = true;
   do_best_fn_oy   = true;
   do_best_tracks  = true;

   return;
}

////////////////////////////////////////////////////////////////////////
//
//  Code for class TCGenVxOpt
//
////////////////////////////////////////////////////////////////////////

TCGenVxOpt::TCGenVxOpt() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

TCGenVxOpt::~TCGenVxOpt() {

   clear();
}

////////////////////////////////////////////////////////////////////////

void TCGenVxOpt::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenVxOpt::clear() {

   Desc.clear();
   Model.clear();
   StormId.clear();
   StormName.clear();
   InitBeg = InitEnd = (unixtime) 0;
   ValidBeg = ValidEnd = (unixtime) 0;
   InitHour.clear();
   Lead.clear();
   VxMaskName.clear();
   VxPolyMask.clear();
   VxGridMask.clear();
   VxAreaMask.clear();
   DLandThresh.clear();
   GenesisMatchRadius = bad_data_double;
   GenesisHitRadius = bad_data_double;
   GenesisHitBeg = GenesisHitEnd = bad_data_int;
   GenesisInitDSec = bad_data_int;
   DiscardFlag = false;
   DevFlag = OpsFlag = false;
   CIAlpha = bad_data_double;
   OutputMap.clear();
   NcInfo.clear();
   ValidGenesisDHrThresh.clear();
   BestUniqueFlag = false;

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenVxOpt::process_config(Dictionary &dict) {
   int i, beg, end;
   Dictionary *dict2 = (Dictionary *) 0;
   ConcatString file_name;
   StringArray sa;
   bool status;

   // Conf: desc
   Desc = parse_conf_string(&dict, conf_key_desc);

   // Conf: model
   Model = parse_conf_tc_model(&dict);

   // Conf: storm_id
   StormId = dict.lookup_string_array(conf_key_storm_id);

   // Conf: storm_name
   StormName = dict.lookup_string_array(conf_key_storm_name);

   // Conf: init_beg, init_end
   InitBeg = dict.lookup_unixtime(conf_key_init_beg);
   InitEnd = dict.lookup_unixtime(conf_key_init_end);

   // Conf: valid_beg, valid_end
   ValidBeg = dict.lookup_unixtime(conf_key_valid_beg);
   ValidEnd = dict.lookup_unixtime(conf_key_valid_end);

   // Conf: init_hour
   sa = dict.lookup_string_array(conf_key_init_hour);
   for(i=0; i<sa.n_elements(); i++) {
      InitHour.add(timestring_to_sec(sa[i].c_str()));
   }

   // Conf: lead
   sa = dict.lookup_string_array(conf_key_lead);
   for(i=0; i<sa.n_elements(); i++) {
      Lead.add(timestring_to_sec(sa[i].c_str()));
   }

   // Conf: vx_mask
   if(nonempty(dict.lookup_string(conf_key_vx_mask).c_str())) {
      file_name = replace_path(dict.lookup_string(conf_key_vx_mask));
      mlog << Debug(2) << "Masking File: " << file_name << "\n";
      parse_poly_mask(file_name, VxPolyMask, VxGridMask, VxAreaMask,
                      VxMaskName);
   }

   // Conf: dland_thresh
   DLandThresh = dict.lookup_thresh(conf_key_dland_thresh);

   // Conf: genesis_match_radius
   GenesisMatchRadius =
      dict.lookup_double(conf_key_genesis_match_radius);

   // Conf: genesis_hit_radius
   GenesisHitRadius =
      dict.lookup_double(conf_key_genesis_hit_radius);

   // Conf: genesis_hit_window
   dict2 = dict.lookup_dictionary(conf_key_genesis_hit_window);
   parse_conf_range_int(dict2, beg, end);
   GenesisHitBeg = beg*sec_per_hour;
   GenesisHitEnd = end*sec_per_hour;

   // Conf: genesis_minus_init_diff
   GenesisInitDSec = nint(
      dict.lookup_double(conf_key_genesis_minus_init_diff) *
      sec_per_hour);

   // Conf: discard_init_post_genesis_flag
   DiscardFlag =
      dict.lookup_bool(conf_key_discard_init_post_genesis_flag);

   // Conf: dev_method_flag and ops_method_flag
   DevFlag = dict.lookup_bool(conf_key_dev_method_flag);
   OpsFlag = dict.lookup_bool(conf_key_ops_method_flag);

   if(!DevFlag && !OpsFlag) {
      mlog << Error << "\nTCGenVxOpt::process_config() -> "
           << "at least one of " << conf_key_dev_method_flag
           << " or " << conf_key_ops_method_flag
           << " must be set to true!\n\n";
      exit(1);
   }
   
   // Conf: ci_alpha
   CIAlpha = dict.lookup_double(conf_key_ci_alpha);

   // Conf: output_flag
   OutputMap = parse_conf_output_flag(&dict, txt_file_type, n_txt);

   for(i=0, status=false; i<OutputMap.size(); i++) {
      if(OutputMap[txt_file_type[i]] != STATOutputType_None) {
         status = true;
         break;
      }
   }

   // Check for at least one output line type
   if(!status) {
      mlog << Error << "\nTCGenVxOpt::process_config() -> "
           << "at least one output line type must be requested in \""
           << conf_key_output_flag << "\"!\n\n";
      exit(1);
   }

   // Conf: nc_pairs_flag
   parse_nc_info(dict);

   // Conf: valid_minus_genesis_diff_thresh
   ValidGenesisDHrThresh =
      dict.lookup_thresh(conf_key_valid_minus_genesis_diff_thresh);

   // Conf: bset_unique_flag
   BestUniqueFlag =
      dict.lookup_bool(conf_key_best_unique_flag);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenVxOpt::parse_nc_info(Dictionary &odict) {
   const DictionaryEntry * e = (const DictionaryEntry *) 0;

   e = odict.lookup(conf_key_nc_pairs_flag);

   if(!e) {
      mlog << Error << "\nTCGenVxOpt::parse_nc_info() -> "
           << "lookup failed for key \"" << conf_key_nc_pairs_flag
           << "\"\n\n";
      exit(1);
   }

   const ConfigObjectType type = e->type();

   if(type == BooleanType) {
      bool value = e->b_value();

      if(!value) NcInfo.set_all_false();

      return;
   }

   // It should be a dictionary
   if(type != DictionaryType) {
      mlog << Error << "\nTCGenVxOpt::parse_nc_info() -> "
           << "bad type for key \"" << conf_key_nc_pairs_flag
           << "\"\n\n";
      exit(1);
   }

   // Parse the various entries
   Dictionary * d = e->dict_value();

   NcInfo.do_latlon       = d->lookup_bool(conf_key_latlon_flag);
   NcInfo.do_fcst_genesis = d->lookup_bool(conf_key_fcst_genesis);
   NcInfo.do_fcst_fy_oy   = d->lookup_bool(conf_key_fcst_fy_oy);
   NcInfo.do_fcst_fy_on   = d->lookup_bool(conf_key_fcst_fy_on);
   NcInfo.do_fcst_tracks  = d->lookup_bool(conf_key_fcst_tracks);
   NcInfo.do_best_genesis = d->lookup_bool(conf_key_best_genesis);
   NcInfo.do_best_fy_oy   = d->lookup_bool(conf_key_best_fy_oy);
   NcInfo.do_best_fn_oy   = d->lookup_bool(conf_key_best_fn_oy);
   NcInfo.do_best_tracks  = d->lookup_bool(conf_key_best_tracks);

   return;
}

////////////////////////////////////////////////////////////////////////

bool TCGenVxOpt::is_keeper(const GenesisInfo &gi) const {
   bool keep = true;

   // ATCF ID processed elsewhere

   // Only check basin, storm ID, cyclone number, and storm name for
   // BEST and operational tracks.

   if(gi.is_best_track() || gi.is_oper_track()) {

      // Check storm id
      if(StormId.n() > 0 &&
         !has_storm_id(StormId, gi.basin(), gi.cyclone(), gi.init()))
         keep = false;

      // Check storm name
      if(StormName.n() > 0 && !StormName.has(gi.storm_name()))
         keep = false;
   }

   if(!keep) return(keep);

   // Only check intialization and lead times for forecast and
   // operational tracks.

   if(!gi.is_best_track() || gi.is_oper_track()) {

      // Initialization time window
      if((InitBeg > 0 &&  InitBeg > gi.init()) ||
         (InitEnd > 0 &&  InitEnd < gi.init()))
         keep = false;

      // Initialization hours
      if(InitHour.n() > 0 && !InitHour.has(gi.init_hour()))
         keep = false;

      // Lead times
      if(Lead.n() > 0 && !Lead.has(gi.genesis_lead()))
         keep = false;
   }

   if(!keep) return(keep);

   // Valid time window
   if((ValidBeg > 0 && ValidBeg > gi.valid_min()) ||
      (ValidEnd > 0 && ValidEnd < gi.valid_max()))
      keep = false;

   // Poly masking
   if(VxPolyMask.n_points() > 0 &&
     !VxPolyMask.latlon_is_inside(gi.lat(), gi.lon()))
      keep = false;

   // Area masking
   if(!VxAreaMask.is_empty()) {
      double x, y;
      VxGridMask.latlon_to_xy(gi.lat(), -1.0*gi.lon(), x, y);
      if(x < 0 || x >= VxGridMask.nx() ||
         y < 0 || y >= VxGridMask.ny()) {
         keep = false;
      }
      else {
         keep = VxAreaMask(nint(x), nint(y));
      }
   }

   // Distance to land
   if((DLandThresh.get_type() != no_thresh_type) &&
      (is_bad_data(gi.dland()) || !DLandThresh.check(gi.dland())))
      keep = false;

   // Return the keep status
   return(keep);
}

////////////////////////////////////////////////////////////////////////

STATOutputType TCGenVxOpt::output_map(STATLineType t) const {
   return(OutputMap.at(t));
}

////////////////////////////////////////////////////////////////////////
//
//  Code for class TCGenConfInfo
//
////////////////////////////////////////////////////////////////////////

TCGenConfInfo::TCGenConfInfo() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

TCGenConfInfo::~TCGenConfInfo() {

   clear();
}

////////////////////////////////////////////////////////////////////////

void TCGenConfInfo::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenConfInfo::clear() {

   for(size_t i=0; i<VxOpt.size(); i++) VxOpt[i].clear();
   InitFreqHr = bad_data_int;
   ValidFreqHr = bad_data_int;
   FcstSecBeg = bad_data_int;
   FcstSecEnd = bad_data_int;
   MinDur = bad_data_int;
   FcstEventInfo.clear();
   BestEventInfo.clear();
   OperTechnique.clear();
   DLandFile.clear();
   DLandGrid.clear();
   DLandData.clear();
   BasinFile.clear();
   BasinGrid.clear();
   BasinData.clear();
   NcOutGrid.clear();
   OutputMap.clear();
   NcInfo.clear();
   Version.clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenConfInfo::read_config(const char *default_file_name,
                                const char *user_file_name) {

   // Read the config file constants
   Conf.read(replace_path(config_const_filename).c_str());

   // Read the default config file
   Conf.read(default_file_name);

   // Read the user-specified config file
   Conf.read(user_file_name);

   // Process the configuration file
   process_config();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenConfInfo::process_config() {
   Dictionary *dict = (Dictionary *) 0;
   TCGenVxOpt vx_opt;
   bool status;
   int i, beg, end;

   // Conf: init_freq
   InitFreqHr = Conf.lookup_int(conf_key_init_freq);

   if(InitFreqHr <= 0) {
      mlog << Error << "\nTCGenConfInfo::process_config() -> "
           << "\"" << conf_key_init_freq << "\" must be greater than "
           << "zero!\n\n";
      exit(1);
   }

   // Conf: valid_freq
   ValidFreqHr = Conf.lookup_int(conf_key_valid_freq);

   if(ValidFreqHr <= 0) {
      mlog << Error << "\nTCGenConfInfo::process_config() -> "
           << "\"" << conf_key_valid_freq << "\" must be greater than "
           << "zero!\n\n";
      exit(1);
   }

   // Conf: fcst_hr_window
   dict = Conf.lookup_dictionary(conf_key_fcst_hr_window);
   parse_conf_range_int(dict, beg, end);
   FcstSecBeg = beg*sec_per_hour;
   FcstSecEnd = end*sec_per_hour;

   // Conf: min_duration
   MinDur = Conf.lookup_int(conf_key_min_duration);

   // Conf: fcst_genesis
   dict = Conf.lookup_dictionary(conf_key_fcst_genesis);
   FcstEventInfo = parse_conf_genesis_event_info(dict);

   // Conf: best_genesis
   dict = Conf.lookup_dictionary(conf_key_best_genesis);
   BestEventInfo = parse_conf_genesis_event_info(dict);

   // Conf: oper_technique
   OperTechnique = Conf.lookup_string(conf_key_oper_technique);

   // Conf: DLandFile
   DLandFile = Conf.lookup_string(conf_key_dland_file);

   // Conf: BasinFile
   BasinFile = Conf.lookup_string(conf_key_basin_file);

   // Conf: NcOutGrid
   parse_grid_mask(Conf.lookup_string(conf_key_nc_pairs_grid),
                   NcOutGrid);

   // Conf: Version
   Version = Conf.lookup_string(conf_key_version);
   check_met_version(Version.c_str());

   // Conf: Filter
   dict = Conf.lookup_array(conf_key_filter, false);

   // If no filters are specified, use the top-level settings
   if(dict->n_entries() == 0) {
      vx_opt.process_config(Conf);
      VxOpt.push_back(vx_opt);
   }
   // Loop through the array entries
   else {
      for(i=0; i<dict->n_entries(); i++) {
         vx_opt.clear();
         vx_opt.process_config(*((*dict)[i]->dict_value()));
         VxOpt.push_back(vx_opt);
      }
   }

   // Check that each filter has a unique description
   if(VxOpt.size() > 1) {
      StringArray desc_sa;
      for(i=0; i<VxOpt.size(); i++) {
         if(desc_sa.has(VxOpt[i].Desc)) {
             mlog << Error << "\nTCGenConfInfo::process_config() -> "
                  << "a unique \"" << conf_key_desc
                  << "\" entry must be specified for each \""
                  << conf_key_filter << "\" array entry!\n\n";
             exit(1);
         }
         else {
            desc_sa.add(VxOpt[i].Desc);
         }
      } // end for i
   }

   // Loop through the filters
   for(size_t j=0; j<VxOpt.size(); j++) {

      // Update the summary OutputMap and NcInfo
      process_flags(VxOpt[j].OutputMap, VxOpt[j].NcInfo);

      // If not already set, define the valid time window relative to the
      // initialization time window.
      if(VxOpt[j].InitBeg != 0 && VxOpt[j].ValidBeg == 0) {
         VxOpt[j].ValidBeg = VxOpt[j].InitBeg + FcstSecBeg;
         mlog << Debug(3) << "For filter " << j+1 << " setting "
              << conf_key_valid_beg << " ("
              << unix_to_yyyymmdd_hhmmss(VxOpt[j].ValidBeg)
              <<  ") = " << conf_key_init_beg << " ("
              << unix_to_yyyymmdd_hhmmss(VxOpt[j].InitBeg)
              << ") + " << conf_key_fcst_hr_window << ".beg ("
              << FcstSecBeg/sec_per_hour << ").\n";
      }

      if(VxOpt[j].InitEnd != 0 && VxOpt[j].ValidEnd == 0) {
         VxOpt[j].ValidEnd = VxOpt[j].InitEnd + FcstSecEnd;
         mlog << Debug(3) << "For filter " << j+1 << " setting "
              << conf_key_valid_end << " ("
              << unix_to_yyyymmdd_hhmmss(VxOpt[j].ValidEnd)
              <<  ") = " << conf_key_init_end << " ("
              << unix_to_yyyymmdd_hhmmss(VxOpt[j].InitEnd)
              << ") + " << conf_key_fcst_hr_window << ".end ("
              << FcstSecEnd/sec_per_hour << ").\n";
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void TCGenConfInfo::process_flags(
        const map<STATLineType,STATOutputType> &m,
        const TCGenNcOutInfo &n) {
   int i;

   // Initialize output map
   if(OutputMap.empty()) OutputMap = m;

   // Update output map
   for(i=0; i<n_txt; i++) {
      if(m.at(txt_file_type[i]) == STATOutputType_Both) {
         OutputMap[txt_file_type[i]] = STATOutputType_Both;
      }
      else if(m.at(txt_file_type[i]) == STATOutputType_Stat &&
              OutputMap[txt_file_type[i]] == STATOutputType_None) {
         OutputMap[txt_file_type[i]] = STATOutputType_Stat;
      }
   }

   // Update NcInfo flags
   NcInfo += n;
   
   return;
}
////////////////////////////////////////////////////////////////////////

double TCGenConfInfo::compute_dland(double lat, double lon) {
   double x_dbl, y_dbl, dist;
   int x, y;

   // Load the distance to land data, if needed.
   if(DLandData.is_empty()) {
      load_tc_dland(DLandFile, DLandGrid, DLandData);
   }

   // Convert lat,lon to x,y
   DLandGrid.latlon_to_xy(lat, lon, x_dbl, y_dbl);

   // Round to nearest int
   x = nint(x_dbl);
   y = nint(y_dbl);

   // Check range
   if(x < 0 || x >= DLandGrid.nx() ||
      y < 0 || y >= DLandGrid.ny())   dist = bad_data_double;
   else                               dist = DLandData.get(x, y);

   return(dist);
}

////////////////////////////////////////////////////////////////////////

ConcatString TCGenConfInfo::compute_basin(double lat, double lon) {
   double x_dbl, y_dbl, dist;
   int x, y, id;
   ConcatString abbr;

   // Load the basin data, if needed.
   if(BasinData.is_empty()) {
      load_tc_basin(BasinFile, BasinGrid, BasinData);
   }

   // Convert lat,lon to x,y
   BasinGrid.latlon_to_xy(lat, lon, x_dbl, y_dbl);

   // Round to nearest int
   x = nint(x_dbl);
   y = nint(y_dbl);

   // Basin ID
   id = ((x < 0 || x >= BasinGrid.nx() ||
          y < 0 || y >= BasinGrid.ny()) ?
         bad_data_int :
         nint(BasinData.get(x, y)));

   // Convert basin ID to string
   switch(id) {
      case 0:
         abbr = "NONE";
         break;
      case 1:
         // Atlantic
         abbr = "AL";
         break;
      case 2:
         // Eastern Pacific
         abbr = "EP";
         break;
      case 3:
         // Central Pacific
         abbr = "CP";
         break;
      case 4:
         // Western Pacific
         abbr = "WP";
         break;
      case 5:
         // Northern Indian Ocean
         abbr = "NI";
         break;
      case 6:
         // Southern Indian Ocean
         abbr = "SI";
         break;
      case 7:
         // Autstralia
         abbr = "AU";
         break;
      case 8:
         // Southern Pacific
         abbr = "SP";
         break;
      default:
         abbr = "UNKNOWN";
         break;
   }

   return(abbr);
}

////////////////////////////////////////////////////////////////////////

STATOutputType TCGenConfInfo::output_map(STATLineType t) const {
   return(OutputMap.at(t));
}

////////////////////////////////////////////////////////////////////////
//
// Code for class GenCTCInfo
//
////////////////////////////////////////////////////////////////////////

GenCTCInfo::GenCTCInfo() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

GenCTCInfo::~GenCTCInfo() {

   clear();
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::init_from_scratch() {

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::clear() {

   Model.clear();
   FcstBeg = FcstEnd = (unixtime) 0;
   BestBeg = BestEnd = (unixtime) 0;

   CTSDev.clear();
   CTSOps.clear();

   VxOpt     = (const TCGenVxOpt *) 0;
   NcOutGrid = (const Grid *) 0;

   ValidGenesisDHrThresh.clear();
   BestUniqueFlag = false;

   BestDevHitMap.clear();
   BestOpsHitMap.clear();

   DpMap.clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::set_vx_opt(const TCGenVxOpt *vx_opt,
                            const Grid *nc_out_grid) {

   if(!vx_opt) return;

   // Store pointer
   VxOpt = vx_opt;

   // Store config options
   ValidGenesisDHrThresh = VxOpt->ValidGenesisDHrThresh;
   BestUniqueFlag        = VxOpt->BestUniqueFlag;

   // Setup alpha value
   if(VxOpt->DevFlag) {
      CTSDev.allocate_n_alpha(1);
      CTSDev.alpha[0] = VxOpt->CIAlpha;
   }
   if(VxOpt->OpsFlag) {
      CTSOps.allocate_n_alpha(1);
      CTSOps.alpha[0] = VxOpt->CIAlpha;
   }

   // Setup NetCDF pairs output fields
   if(!VxOpt->NcInfo.all_false()) {
      NcOutGrid = nc_out_grid;

      // Initialize data plane of all zeros
      DataPlane dp;
      dp.set_size(NcOutGrid->nx(), NcOutGrid->ny(), 0.0);

      // Add map entries for requested outputs
      if(VxOpt->NcInfo.do_fcst_genesis) DpMap[fgen_str] = dp;
      if(VxOpt->NcInfo.do_fcst_tracks)  DpMap[ftrk_str] = dp;
      if(VxOpt->NcInfo.do_best_genesis) DpMap[bgen_str] = dp;
      if(VxOpt->NcInfo.do_best_tracks)  DpMap[btrk_str] = dp;
      if(VxOpt->DevFlag) {
         if(VxOpt->NcInfo.do_fcst_fy_oy) DpMap[fdev_fyoy_str] = dp;
         if(VxOpt->NcInfo.do_fcst_fy_on) DpMap[fdev_fyon_str] = dp;
         if(VxOpt->NcInfo.do_best_fy_oy) DpMap[bdev_fyoy_str] = dp;
         if(VxOpt->NcInfo.do_best_fn_oy) DpMap[bdev_fnoy_str] = dp;
      }
      if(VxOpt->OpsFlag) {
         if(VxOpt->NcInfo.do_fcst_fy_oy) DpMap[fops_fyoy_str] = dp;
         if(VxOpt->NcInfo.do_fcst_fy_on) DpMap[fops_fyon_str] = dp;
         if(VxOpt->NcInfo.do_best_fy_oy) DpMap[bops_fyoy_str] = dp;
         if(VxOpt->NcInfo.do_best_fn_oy) DpMap[bops_fnoy_str] = dp;
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::inc_dev(bool f, bool o,
                         const GenesisInfo *fgi,
                         const GenesisInfo *bgi) {

   // Hits
   if(f && o) {
      CTSDev.cts.inc_fy_oy();
      inc_pnt(fgi->lat(), fgi->lon(), fdev_fyoy_str);
      BestDevHitMap[bgi] += 1;
      
      // Count all BEST track genesis pairs
      if(!BestUniqueFlag) {
         inc_pnt(bgi->lat(), bgi->lon(), bdev_fyoy_str);
      }
   }
   // False Alarms
   else if(f && !o) {
      CTSDev.cts.inc_fy_on();
      inc_pnt(fgi->lat(), fgi->lon(), fdev_fyon_str);
   }
   // Misses
   else if(!f && o) {
      CTSDev.cts.inc_fn_oy();

      // Count all BEST track genesis pairs
      if(!BestUniqueFlag) {
         inc_pnt(bgi->lat(), bgi->lon(), bdev_fnoy_str);
      }
   }
   // Correct Negatives (should be none)
   else {
      CTSDev.cts.inc_fn_on();
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::inc_ops(bool f, bool o,
                         const GenesisInfo *fgi,
                         const GenesisInfo *bgi) {

   // Hits
   if(f && o) {
      CTSOps.cts.inc_fy_oy();
      inc_pnt(fgi->lat(), fgi->lon(), fops_fyoy_str);
      BestOpsHitMap[bgi] += 1;

      // Count all BEST track genesis pairs
      if(!BestUniqueFlag) {
         inc_pnt(bgi->lat(), bgi->lon(), bops_fyoy_str);
      }
   }
   // False Alarms
   else if(f && !o) {
      CTSOps.cts.inc_fy_on();
      inc_pnt(fgi->lat(), fgi->lon(), fops_fyon_str);
   }
   // Misses
   else if(!f && o) {
      CTSOps.cts.inc_fn_oy();

      // Count all BEST track genesis pairs
      if(!BestUniqueFlag) {
         inc_pnt(bgi->lat(), bgi->lon(), bops_fnoy_str);
      }
   }
   // Correct Negatives (should be none)
   else {
      CTSOps.cts.inc_fn_on();
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::inc_best_unique() {

   // Only process when the flag is set
   if(!BestUniqueFlag) return;

   map<const GenesisInfo *,int>::iterator it;

   // Count the dev BEST track hits and false alarms
   for(it=BestDevHitMap.begin(); it!=BestDevHitMap.end(); it++) {

      // Zero hits is a miss
      if(it->second == 0) {
         inc_pnt(it->first->lat(), it->first->lon(), bdev_fnoy_str);
      }
      // Otherwise, it's a hit
      else {
         inc_pnt(it->first->lat(), it->first->lon(), bdev_fyoy_str);
      }
   }

   // Count the ops BEST track hits and false alarms
   for(it=BestOpsHitMap.begin(); it!=BestOpsHitMap.end(); it++) {

      // Zero hits is a miss
      if(it->second == 0) {
         inc_pnt(it->first->lat(), it->first->lon(), bops_fnoy_str);
      }
      // Otherwise, it's a hit
      else {
         inc_pnt(it->first->lat(), it->first->lon(), bops_fyoy_str);
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::add_fcst_gen(const GenesisInfo &gi) {

   // Track the range of valid times
   if(FcstBeg == 0 || FcstBeg > gi.valid_min()) FcstBeg = gi.valid_min();
   if(FcstEnd == 0 || FcstEnd < gi.valid_max()) FcstEnd = gi.valid_max();
   
   // Count the genesis and track points
   inc_pnt(gi.lat(), gi.lon(), fgen_str);
   inc_trk(gi, ftrk_str);

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::add_best_gen(const GenesisInfo &gi) {

   // Nothing to do if this genesis has already been counted
   if(BestDevHitMap.count(&gi) > 0 ||
      BestOpsHitMap.count(&gi) > 0) return;

   // Add hit counter entries for this storm
   BestDevHitMap[&gi] = 0;
   BestOpsHitMap[&gi] = 0;

   // Track the range of valid times
   if(BestBeg == 0 || BestBeg > gi.valid_min()) BestBeg = gi.valid_min();
   if(BestEnd == 0 || BestEnd < gi.valid_max()) BestEnd = gi.valid_max();

   // Count the genesis and track points
   inc_pnt(gi.lat(), gi.lon(), bgen_str);
   inc_trk(gi, btrk_str);

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::inc_pnt(double lat, double lon, const string &s) {

   // Nothing to do if there is no DataPlane map entry
   if(DpMap.count(s) == 0) return;

   int x, y;
   double x_dbl, y_dbl;

   NcOutGrid->latlon_to_xy(lat, -1.0*lon, x_dbl, y_dbl);
   x = nint(x_dbl);
   y = nint(y_dbl);

   // Only increment points on the grid
   if(x >= 0 && x < NcOutGrid->nx() &&
      y >= 0 && y < NcOutGrid->ny()) {
      DpMap[s].set((DpMap[s])(x, y) + 1, x, y);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void GenCTCInfo::inc_trk(const GenesisInfo &gi, const string &s) {

   // Nothing to do if there is no DataPlane map entry
   if(DpMap.count(s) == 0) return;

   // Loop through the track points
   for(int i=0; i<gi.n_points(); i++) {

      // Count points whose valid time is close enough to genesis time
      int dhr = (gi[i].valid() - gi.genesis_time())/sec_per_hour;
      if(ValidGenesisDHrThresh.check(dhr)) {
         inc_pnt(gi[i].lat(), gi[i].lon(), s);
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////
