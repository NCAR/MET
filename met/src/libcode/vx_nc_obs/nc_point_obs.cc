// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2020
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


////////////////////////////////////////////////////////////////////////


using namespace std;

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cmath>
#include <time.h>

#include "vx_log.h"

#include "nc_point_obs.h"
#include "write_netcdf.h"

////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////


   //
   //  Code for class MetNcPointObs
   //


////////////////////////////////////////////////////////////////////////

MetNcPointObs::MetNcPointObs() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

MetNcPointObs::~MetNcPointObs() {
   close();
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs::init_from_scratch() {
   obs_nc = (NcFile *) 0;

   nobs = 0;
   nhdr = 0;
   qty_length = 0;
   keep_nc = false;
   use_var_id = false;
   use_arr_vars = false;

   //close();
   return;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::check_nc(const char *nc_name, const char *caller) {
   bool exit_on_error = false;
   bool valid = obs_vars.is_valid(exit_on_error);
   if (!valid) {
      mlog << Error << "\n" << (0 != caller ? caller : "") << " -> "
           << "missing core data from the netCDF file: "
           << nc_name << "\n\n";
      exit(1);
   }
   return valid;
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs::close() {
   if ( !keep_nc && obs_nc ) {
      delete obs_nc;
      obs_nc = (NcFile *) 0;
   }

   //obs_vars.reset();
   obs_data.clear();
   header_data.clear();
   return;
}


////////////////////////////////////////////////////////////////////////
// Requies nc_summary
//void MetNcPointObs::get_dim_counts(int *obs_cnt, int *hdr_cnt) {
//   reset_hdr_buffer = get_hdr_obs_count(obs_nc, &out_data, obs_cnt, hdr_cnt);
//}

////////////////////////////////////////////////////////////////////////

int MetNcPointObs::get_qty_length() {
   qty_length = get_nc_string_length(&obs_vars.obs_qty_tbl_var);
   return qty_length;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::get_header(int header_offset, float hdr_arr[HDR_ARRAY_LEN],
      ConcatString &hdr_typ_str, ConcatString &hdr_sid_str, ConcatString &hdr_vld_str) {
   int hdr_idx;

   // Read the corresponding header array for this observation
   hdr_arr[0] = header_data.lat_array[header_offset];
   hdr_arr[1] = header_data.lon_array[header_offset];
   hdr_arr[2] = header_data.elv_array[header_offset];

   // Read the corresponding header type for this observation
   hdr_idx = use_arr_vars ? header_offset : header_data.typ_idx_array[header_offset];
   hdr_typ_str = header_data.typ_array[hdr_idx];

   // Read the corresponding header Station ID for this observation
   hdr_idx = use_arr_vars ? header_offset : header_data.sid_idx_array[header_offset];
   hdr_sid_str = header_data.sid_array[hdr_idx];

   // Read the corresponding valid time for this observation
   hdr_idx = use_arr_vars ? header_offset : header_data.vld_idx_array[header_offset];
   hdr_vld_str =  header_data.vld_array[hdr_idx];

   return true;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::get_header_type(int header_offset, int hdr_typ_arr[HDR_TYPE_ARR_LEN]) {
   int hdr_idx;
   // Read the header integer types
   hdr_typ_arr[0] = (header_data.prpt_typ_array.n() > header_offset ?
                     header_data.prpt_typ_array[header_offset] : bad_data_int);
   hdr_typ_arr[1] = (header_data.irpt_typ_array.n() > header_offset ?
                     header_data.irpt_typ_array[header_offset] : bad_data_int);
   hdr_typ_arr[2] = (header_data.inst_typ_array.n() > header_offset ?
                     header_data.inst_typ_array[header_offset] : bad_data_int);

   return true;
}


////////////////////////////////////////////////////////////////////////
//
//int MetNcPointObs::get_hdr_arr_len() {
//   int hdr_arr_len = IS_INVALID_NC(obs_vars.hdr_arr_dim)
//                      ? HDR_ARRAY_LEN : get_dim_size(&obs_vars.hdr_arr_dim);
//   return hdr_arr_len;
//}
//
////////////////////////////////////////////////////////////////////////////
//
//int MetNcPointObs::get_obs_arr_len() {
//   int obs_arr_len = IS_INVALID_NC(obs_vars.hdr_arr_dim)
//                     ? OBS_ARRAY_LEN : get_dim_size(&obs_vars.obs_arr_dim);
//}
//
////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::get_lats(float *hdr_lats) {
   for (int idx=0; idx<nhdr; idx++) {
      hdr_lats[idx] = header_data.lat_array[idx];
   }
   return true;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::get_lons(float *hdr_lons) {
   for (int idx=0; idx<nhdr; idx++) {
      hdr_lons[idx] = header_data.lon_array[idx];
   }
   return true;
}

////////////////////////////////////////////////////////////////////////

//StringArray MetNcPointObs::get_qty_data() {
//   const char *method_name = "MetNcPointObs::get_qty_data()";
//   if (IS_VALID_NC(obs_vars.obs_qty_tbl_var)) {
//      if (!get_nc_data_to_array(&obs_vars.obs_qty_tbl_var, &obs_qty_array)) {
//         mlog << Error << "\n" << method_name
//              << "trouble getting obs_qty\n\n";
//         exit(1);
//      }
//   }
//
//   return obs_qty_array;
//}
//
//////////////////////////////////////////////////////////////////////////
//
//StringArray MetNcPointObs::get_var_names(const char *caller) {
//   StringArray _var_names;
//   if (use_var_id) {
//      ConcatString from_caller;
//      if (0 < strlen(caller)) {
//         from_caller.add("from ");
//         from_caller.add(caller);
//      }
//      if (0 < var_data.var_names.n()) _var_names = obs_data.var_names;
//      else {
//        if (get_nc_data_to_array(&obs_vars.obs_var, &_var_names)) {
//           obs_data.var_names = _var_names;
//        }
//        else {
//           mlog << Error << "\nMetNcPointObs::read_var_names() "
//                << from_caller << " -> "
//                << "trouble getting variable names from "
//                << nc_var_obs_var << "\n\n";
//           exit(1);
//        }
//      }
//   }
//   return _var_names;
//}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::open(const char * filename) {
   return set_netcdf(open_ncfile(filename));
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::read_dim_headers() {
   bool status = false;
   if( IS_VALID_NC_P(obs_nc) ) {
      status = true;
      obs_vars.read_dims_vars(obs_nc);
      nobs = obs_vars.obs_cnt = GET_NC_SIZE(obs_vars.obs_dim);
      nhdr = obs_vars.hdr_cnt = GET_NC_SIZE(obs_vars.hdr_dim);
      //obs_vars.set_hdr_count(nhdr);
      //obs_vars.set_obs_count(nobs);
      obs_vars.use_var_id = use_var_id = IS_VALID_NC(obs_vars.obs_vid_var);
      use_arr_vars = IS_VALID_NC(obs_vars.obs_arr_var);
      header_data = obs_vars.read_header_data();
      //FIXME seg fault!!! get_nc_pb_hdr_data(obs_vars, &header_data);
   }
   return status;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::read_obs_data() {
   bool status = read_obs_data_numbers() && read_obs_data_strings();
   return status;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::read_obs_data(int buf_size, int start, float *obs_arr_block,
                                  int *obs_qty_idx_block, char *obs_qty_str_block) {
  return obs_vars.read_obs_data(buf_size, start, qty_length, obs_arr_block,
                                obs_qty_idx_block, obs_qty_str_block);
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::read_obs_data_numbers() {
   bool status = false;
   if( IS_VALID_NC_P(obs_nc) ) {
      status = obs_data.read_obs_data_numbers(obs_vars);
   }
   return status;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::read_obs_data_strings() {
   bool status = false;
   if( IS_VALID_NC_P(obs_nc) ) {
      status = obs_data.read_obs_data_strings(obs_vars);
   }
   return status;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs::set_netcdf(NcFile *nc_file, bool _keep_nc) {
   close();
   keep_nc = _keep_nc;
   obs_nc = nc_file;
   return IS_VALID_NC_P(obs_nc);
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs::set_using_var_id(bool using_var_id) {
   use_var_id = obs_vars.use_var_id = using_var_id; 
}

////////////////////////////////////////////////////////////////////////
