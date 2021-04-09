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
//#include <unistd.h>
//#include <stdlib.h>
#include <string.h>
#include <cstdio>
//#include <cmath>
//#include <time.h>

#include "vx_log.h"

#include "nc_point_obs.h"
#include "nc_summary.h"
#include "write_netcdf.h"

////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////


   //
   //  Code for class MetNcPointObs2Write
   //

////////////////////////////////////////////////////////////////////////

MetNcPointObs2Write::MetNcPointObs2Write() {
   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

MetNcPointObs2Write::~MetNcPointObs2Write() {
   close();
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::init_from_scratch() {
   MetNcPointObs::init_from_scratch();
   raw_hdr_cnt = 0;
   reset_hdr_buffer = false;
}

bool MetNcPointObs2Write::add_header(const char *hdr_typ, const char *hdr_sid,
                                     const time_t hdr_vld, const float hdr_lat,
                                     const float hdr_lon, const float hdr_elv)
{
   bool added = false;
   bool new_vld = false;
   const char *method_name = "MetNcPointObs2Write::add_header() ";

   // Can't filter duplicated one because header index was
   // assigned before checking
   int hdr_index;
   add_header_strings(hdr_typ, hdr_sid);

   if (header_data.min_vld_time == -1 || header_data.min_vld_time > hdr_vld) {
      if (header_data.min_vld_time == -1) header_data.max_vld_time = hdr_vld; 
      header_data.min_vld_time = hdr_vld;
      new_vld = true;
   }
   else if (header_data.max_vld_time < hdr_vld) {
      header_data.max_vld_time = hdr_vld;
      new_vld = true;
   }
   if (new_vld || !header_data.vld_num_array.has(hdr_vld, hdr_index, false)) {
      hdr_index = header_data.vld_array.n_elements();
      header_data.vld_array.add(unix_to_yyyymmdd_hhmmss(hdr_vld)); // Valid time
      header_data.vld_num_array.add(hdr_vld);   // Valid time
   }
   header_data.vld_idx_array.add(hdr_index);       // Index of Valid time

   header_data.lat_array.add(hdr_lat);  // Latitude
   header_data.lon_array.add(hdr_lon);  // Longitude
   header_data.elv_array.add(hdr_elv);  // Elevation
   data_buffer.cur_hdr_idx++;
   mlog << Debug(9) << method_name << "header is added (cur_index="
        << data_buffer.cur_hdr_idx << ")\n";
   added = true;
   return added;
}

///////////////////////////////////////////////////////////////////////////////

bool MetNcPointObs2Write::add_header_prepbufr(const int pb_report_type,
                                              const int in_report_type,
                                              const int instrument_type)
{
   bool added = true;
   const char *method_name = "add_header_prepbufr() ";
   // Can't filter duplicated one because header index was
   // assigned before checking
   header_data.prpt_typ_array.add(pb_report_type);
   header_data.irpt_typ_array.add(in_report_type);
   header_data.inst_typ_array.add(instrument_type);
   mlog << Debug(9) << method_name << "pb_report_type: " << pb_report_type 
        << ", in_report_type: " << in_report_type << ", instrument_type: "
        << instrument_type << "\n";
   return added;
}

///////////////////////////////////////////////////////////////////////////////

bool MetNcPointObs2Write::add_header_strings(const char *hdr_typ, const char *hdr_sid)
{
   bool added = false;
   
   // Can't filter duplicated one because header index was
   // assigned before checking
   int hdr_index;
   if (!header_data.typ_array.has(hdr_typ, hdr_index, false)) {
      hdr_index = header_data.typ_array.n_elements();
      header_data.typ_array.add(hdr_typ);          // Message type
      added = true;
   }
   header_data.typ_idx_array.add(hdr_index);       // Index of Message type
   
   if (!header_data.sid_array.has(hdr_sid, hdr_index, false)) {
      hdr_index = header_data.sid_array.n_elements();
      header_data.sid_array.add(hdr_sid);          // Station ID
      added = true;
   }
   header_data.sid_idx_array.add(hdr_index);       // Index of Station ID
 
   return added;
}

///////////////////////////////////////////////////////////////////////////////

bool MetNcPointObs2Write::add_header_vld(const char *hdr_vld)
{
   bool added = false;
   
   // Can't filter duplicated one because header index was
   // assigned before checking
   int hdr_index;
   if (!header_data.vld_array.has(hdr_vld, hdr_index, false)) {
      hdr_index = header_data.typ_array.n_elements();
      header_data.vld_array.add(hdr_vld);          // Valid time
      added = true;
   }
   return added;
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::create_pb_hdrs(int pb_hdr_count) {
   raw_hdr_cnt = pb_hdr_count;
   obs_vars.create_pb_hdrs(obs_nc, pb_hdr_count);
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::get_dim_counts(int *obs_cnt, int *hdr_cnt) {
   reset_hdr_buffer = get_hdr_obs_count(obs_nc, &out_data, obs_cnt, hdr_cnt);
}

////////////////////////////////////////////////////////////////////////
// If raw_hdr_cnt is greater than 0, skip updating header index for obs.

void MetNcPointObs2Write::init_buffer() {
   //string method_name = "init_buffer(() ";

   data_buffer.obs_data_idx    = 0;
   data_buffer.obs_data_offset = 0;
   data_buffer.hdr_data_idx    = 0;
   data_buffer.hdr_data_offset = 0;
   data_buffer.pb_hdr_data_offset = 0;

   strcpy(data_buffer.prev_hdr_typ_buf, "NotDefined");
   strcpy(data_buffer.prev_hdr_sid_buf, "NotDefined");
   strcpy(data_buffer.prev_hdr_vld_buf, "NotDefined");
   for (int index=0; index<HDR_ARRAY_LEN; index++)
      data_buffer.prev_hdr_arr_buf[index] = 0.0;

   header_data.clear();
}

///////////////////////////////////////////////////////////////////////////////
// If raw_hdr_cnt is greater than 0, skip updating header index for obs.

bool MetNcPointObs2Write::init_netcdf(int obs_count, int hdr_count,
                                     string program_name) {
   string method_name = "init_netcdf() ";

   if (reset_hdr_buffer) {
      bool hdr_cnt = get_hdr_index();
      bool reset_array = true;
      reset_header_buffer(hdr_cnt, reset_array);
      mlog << Debug(5) << method_name << "reset headers (" << hdr_cnt << ") raw data.\n";
   }

   //nobs = obs_count;
   nobs = obs_vars.obs_cnt = obs_count;
   nhdr = obs_vars.hdr_cnt = hdr_count;
   //obs_vars.set_obs_count(obs_count);
   //obs_vars.deflate_level = out_data.deflate_level;
   
   obs_vars.create_hdr_vars(obs_nc, hdr_count);
   obs_vars.create_obs_vars(obs_nc);

   //
   // Add global attributes
   //
   write_netcdf_global(obs_nc, obs_nc->getName().c_str(), program_name.c_str());

   return reset_hdr_buffer;
}

////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::init_obs_vars(bool using_var_id, int deflate_level,
                                        bool attr_agl) {
   use_var_id = using_var_id;
   obs_vars.reset(using_var_id);
   obs_vars.attr_agl = attr_agl;
//FIXME   out_data.deflate_level = obs_vars.deflate_level = deflate_level;

   header_data.clear();
}

///////////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::reset_header_buffer(int buf_size, bool reset_all) {
   for (int i=0; i<buf_size; i++) {
      for (int j=0; j<HEADER_STR_LEN; j++) {
         data_buffer.hdr_typ_str_buf[i][j] = bad_data_char;
         data_buffer.hdr_sid_str_buf[i][j] = bad_data_char;
         data_buffer.hdr_vld_str_buf[i][j] = bad_data_char;
      }
      for (int j=HEADER_STR_LEN; j<HEADER_STR_LEN2; j++) {
         data_buffer.hdr_typ_str_buf[i][j] = bad_data_char;
         data_buffer.hdr_sid_str_buf[i][j] = bad_data_char;
      }
      for (int j=0; j<HDR_ARRAY_LEN; j++) {
         data_buffer.hdr_arr_buf[i][j] = FILL_VALUE;
      }
      
      data_buffer.hdr_lat_buf[i] = FILL_VALUE;
      data_buffer.hdr_lon_buf[i] = FILL_VALUE;
      data_buffer.hdr_elv_buf[i] = FILL_VALUE;
   }
   
   if (reset_all) {
      data_buffer.cur_hdr_idx = 0;
      header_data.clear();
   }
}


// Saves the headers at NcHeaderData header_data
//
void MetNcPointObs2Write::write_arr_headers() {
   int cur_hdr_idx = data_buffer.cur_hdr_idx;
   int buf_size = (cur_hdr_idx > OBS_BUFFER_SIZE) ? OBS_BUFFER_SIZE : cur_hdr_idx;
   const string method_name = "  write_arr_headers()";
   
   mlog << Debug(5) << method_name << "  hdr_count: " << cur_hdr_idx
        << ", typ_idx_array: " << header_data.typ_idx_array.n_elements()
        << ", sid_idx_array: " << header_data.sid_idx_array.n_elements()
        << ", vld_idx_array: " << header_data.vld_idx_array.n_elements()
        << ", lat_array: " << header_data.lat_array.n_elements()
        << ", lon_array: " << header_data.lon_array.n_elements()
        << ", elv_array: " << header_data.elv_array.n_elements()
        << "\n";

   int hdr_data_idx = 0;
   bool is_pb_hdr = (0 < header_data.prpt_typ_array.n_elements())
         && !IS_INVALID_NC(obs_vars.hdr_prpt_typ_var);
   // FIXME nc_data_buffer.obs_vars = obs_vars;
   data_buffer.hdr_buf_size = buf_size;
   data_buffer.hdr_data_idx = hdr_data_idx;
   int pb_raw_hdr_cnt = raw_hdr_cnt;
   if (0 == pb_raw_hdr_cnt) pb_raw_hdr_cnt = obs_vars.raw_hdr_cnt;
   for (int index=0; index<cur_hdr_idx; index++) {
      // Message type
      data_buffer.hdr_typ_buf[hdr_data_idx] = header_data.typ_idx_array[index];
      
      // Station ID
      data_buffer.hdr_sid_buf[hdr_data_idx] = header_data.sid_idx_array[index];
      
      // Valid Time
      data_buffer.hdr_vld_buf[hdr_data_idx] = header_data.vld_idx_array[index];
      
      // Write the header array which consists of the following:
      //    LAT LON ELV
      data_buffer.hdr_lat_buf[hdr_data_idx] = (float) header_data.lat_array[index];
      data_buffer.hdr_lon_buf[hdr_data_idx] = (float) header_data.lon_array[index];
      data_buffer.hdr_elv_buf[hdr_data_idx] = (float) header_data.elv_array[index];
      
      if (is_pb_hdr && index < raw_hdr_cnt) {
         data_buffer.hdr_prpt_typ_buf[hdr_data_idx] = (float) header_data.prpt_typ_array[index];
         data_buffer.hdr_irpt_typ_buf[hdr_data_idx] = (float) header_data.irpt_typ_array[index];
         data_buffer.hdr_inst_typ_buf[hdr_data_idx] = (float) header_data.inst_typ_array[index];
      }
      
      hdr_data_idx++;
      data_buffer.hdr_data_idx = hdr_data_idx;
      
      if (hdr_data_idx >= buf_size) {
         obs_vars.write_header_to_nc(data_buffer, hdr_data_idx, is_pb_hdr);
         hdr_data_idx = data_buffer.hdr_data_idx;
      }
   }

   write_buf_headers();
}

///////////////////////////////////////////////////////////////////////////////
//FIXME to remove

void MetNcPointObs2Write::write_buf_headers () {
   if (0 < data_buffer.hdr_data_idx) {
      obs_vars.write_header_to_nc(data_buffer, data_buffer.hdr_data_idx);
   }
   obs_vars.write_table_vars(header_data, data_buffer);
}

///////////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::write_header(const char *hdr_typ, const char *hdr_sid,
                                       const time_t hdr_vld, const float hdr_lat,
                                       const float hdr_lon, const float hdr_elv)
{
   // Can't filter duplicated one because header index was
   // assigned before checking
   int hdr_index;
   bool new_vld = false;
   int hdr_data_idx = data_buffer.hdr_data_idx;
   const char *method_name = "MetNcPointObs2Write::write_header";

   // Message type
   if (!header_data.typ_array.has(hdr_typ, hdr_index, false)) {
      hdr_index = header_data.typ_array.n_elements();
      header_data.typ_array.add(hdr_typ);
   }
   data_buffer.hdr_typ_buf[hdr_data_idx] = hdr_index;
   
   // Station ID
   if (!header_data.sid_array.has(hdr_sid, hdr_index, false)) {
      hdr_index = header_data.sid_array.n_elements();
      header_data.sid_array.add(hdr_sid);
   }
   data_buffer.hdr_sid_buf[hdr_data_idx] = hdr_index;
   
   // Valid Time
   if (header_data.min_vld_time == -1 || header_data.min_vld_time > hdr_vld) {
      if (header_data.min_vld_time == -1) header_data.max_vld_time = hdr_vld; 
      header_data.min_vld_time = hdr_vld;
      new_vld = true;
   }
   else if (header_data.max_vld_time < hdr_vld) {
      header_data.max_vld_time = hdr_vld;
      new_vld = true;
   }
   
   if (new_vld || !header_data.vld_num_array.has(hdr_vld, hdr_index, false)) {
      hdr_index = header_data.vld_array.n_elements();
      header_data.vld_array.add(unix_to_yyyymmdd_hhmmss(hdr_vld));
      header_data.vld_num_array.add(hdr_vld);
   }
   data_buffer.hdr_vld_buf[hdr_data_idx] = hdr_index;
   
   // Write the header array which consists of the following:
   //    LAT LON ELV
   data_buffer.hdr_lat_buf[hdr_data_idx] = (float) hdr_lat;
   data_buffer.hdr_lon_buf[hdr_data_idx] = (float) hdr_lon;
   data_buffer.hdr_elv_buf[hdr_data_idx] = (float) hdr_elv;
   
   hdr_data_idx++;
   data_buffer.hdr_data_idx = hdr_data_idx;
   data_buffer.cur_hdr_idx++;
   mlog << Debug(9) << method_name << "header is added (cur_index="
        << data_buffer.cur_hdr_idx << ")\n";

   if (hdr_data_idx >= OBS_BUFFER_SIZE) {
      obs_vars.write_header_to_nc(data_buffer, OBS_BUFFER_SIZE);
   }
}
      
///////////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::write_observation()
{
   if (0 < data_buffer.obs_data_idx){
      //CHECKME data_buffer.obs_vars = obs_vars;
      obs_vars.write_obs_buffer(data_buffer, data_buffer.obs_data_idx);
   }
}

///////////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::write_observation(const float obs_arr[OBS_ARRAY_LEN],
                                            const char *obs_qty)
{
   int qty_index;
   int obs_data_idx = data_buffer.obs_data_idx;
   if (!data_buffer.qty_data_array.has(obs_qty, qty_index, false)) {
      qty_index = data_buffer.qty_data_array.n_elements();
      data_buffer.qty_data_array.add(obs_qty);
   }
   data_buffer.qty_idx_buf[obs_data_idx] = qty_index;
      
   //for (int idx=0; idx<OBS_ARRAY_LEN; idx++) {
   //   data_buffer.obs_data_buf[obs_data_idx][idx] = obs_arr[idx];
   //}
   data_buffer.obs_hid_buf[obs_data_idx] = obs_arr[0];
   data_buffer.obs_vid_buf[obs_data_idx] = obs_arr[1];
   data_buffer.obs_lvl_buf[obs_data_idx] = obs_arr[2];
   data_buffer.obs_hgt_buf[obs_data_idx] = obs_arr[3];
   data_buffer.obs_val_buf[obs_data_idx] = obs_arr[4];
   data_buffer.obs_data_idx++;
   data_buffer.cur_obs_idx++;
   
   if (data_buffer.obs_data_idx >= OBS_BUFFER_SIZE) {
      obs_vars.write_obs_buffer(data_buffer, OBS_BUFFER_SIZE);
   }
}

///////////////////////////////////////////////////////////////////////////////

void MetNcPointObs2Write::write_obs_data()
{
   string method_name = "write_obs_data() ";
   bool do_summary = out_data.summary_info.flag;
   bool do_save_raw_data = out_data.summary_info.raw_data;
   bool do_header = (out_data.processed_hdr_cnt == 0);
  
   mlog << Debug(5) << method_name << "do_header: " << (do_header ? "true" : "false")
        << ", do_summary: " << (do_summary ? "true" : "false")
        << ", save_raw_data: " << (do_save_raw_data ? "true" : "false")
        << "\n";
  
   if (!do_summary || (do_summary && do_save_raw_data)) {
     mlog << Debug(5) << method_name << "writing " 
          << (int)out_data.observations.size() << " raw data...\n";
     write_obs_data(out_data.observations, do_header);
   }
   if (do_summary) {
     mlog << Debug(5) << method_name << "writing summary"
          << (do_save_raw_data ? " " : " (summary only)") << "...\n";
     bool tmp_do_header = true ;
     write_obs_data(out_data.summary_obs->getSummaries(),
                    tmp_do_header);
   }
  
   int obs_buf_index = get_obs_index();
   if (obs_buf_index > 0) {
     obs_vars.write_obs_buffer(data_buffer, obs_buf_index);
   }
  
   //return true;
}

///////////////////////////////////////////////////////////////////////////////

int MetNcPointObs2Write::write_obs_data(const vector< Observation > observations,
                                        const bool do_header)
{
   int prev_hdr_idx = -1;
   string prev_header_type = "";
   string prev_station_id = "";
   ConcatString obs_qty;
   int headerOffset = data_buffer.cur_hdr_idx;
   //float obs_arr[obs_arr_len];
   const string method_name = "  write_obs_data()";

   int obs_buf_size = observations.size();
   if (obs_buf_size > OBS_BUFFER_SIZE) obs_buf_size = OBS_BUFFER_SIZE;
   
   //if (reset) {
   //   data_buffer.obs_vars = obs_vars;
   //   data_buffer.obs_buf_size = obs_buf_size;
   //   data_buffer.obs_data_idx = 0;
   //   data_buffer.obs_data_offset = 0;
   //   data_buffer.hdr_data_idx = 0;
   //   data_buffer.hdr_data_offset = 0;
   //   
   //   data_buffer.processed_count =0;
   //}
   float obs_arr[OBS_ARRAY_LEN];
   bool header_to_vector = IS_INVALID_NC(obs_vars.hdr_arr_var)
         || IS_INVALID_NC(obs_vars.hdr_lat_var);
   mlog << Debug(5) << method_name << "  obs_count: " << obs_buf_size
        << "  do_header: " << do_header
        << "  header_to_vector: " << header_to_vector << "\n";
   for (vector< Observation >::const_iterator obs = observations.begin();
        obs != observations.end(); ++obs)
   {
      data_buffer.processed_count++;
      
      if (do_header) {
         if (obs->getHeaderIndex() != prev_hdr_idx) {
            mlog << Debug(9) << method_name << "  obs->getHeaderIndex(): "
                 << obs->getHeaderIndex() << " at obs " << data_buffer.processed_count << "\n";
            prev_hdr_idx = obs->getHeaderIndex();
            if (header_to_vector) {
               add_header(
                     obs->getHeaderType().c_str(),
                     obs->getStationId().c_str(),
                     obs->getValidTime(),
                     obs->getLatitude(),
                     obs->getLongitude(),
                     obs->getElevation());
            }
            else {
               write_header(
                     obs->getHeaderType().c_str(),
                     obs->getStationId().c_str(),
                     obs->getValidTime(),
                     obs->getLatitude(),
                     obs->getLongitude(),
                     obs->getElevation());
            }
         }
      }
      
      obs_arr[0] = obs->getHeaderIndex();
      obs_arr[1] = (use_var_id ? obs->getVarCode() : obs->getGribCode());
      obs_arr[2] = obs->getPressureLevel();
      obs_arr[3] = obs->getHeight();
      obs_arr[4] = obs->getValue();
      obs_qty = (obs->getQualityFlag().length() == 0 ? na_str : obs->getQualityFlag().c_str());
      if (do_header) obs_arr[0] += headerOffset;

      write_observation(obs_arr, obs_qty.text());
      
   } // endfor - obs
   
   if (data_buffer.obs_data_idx > 0) {
      obs_vars.write_obs_buffer(data_buffer, data_buffer.obs_data_idx);
   }

   //Caller handles writing headers

   return data_buffer.processed_count;
}

////////////////////////////////////////////////////////////////////////

bool MetNcPointObs2Write::write_to_netcdf(StringArray obs_names, StringArray obs_units,
                                          StringArray obs_descs) {
   const char *method_name = "  write_to_netcdf() ";
   //if( using_obs_vars ) 
   //else 
   //write_observations(&obs_vars, &out_data);
   write_obs_data();
   
   obs_vars.create_table_vars (obs_nc, header_data, data_buffer);
   write_arr_headers();

   int var_count = obs_names.n_elements();
   if (var_count > 0) {
      int unit_count = obs_units.n();
      obs_vars.create_obs_name_vars (obs_nc, var_count, unit_count);
      obs_vars.write_obs_var_names(obs_names);
      if( unit_count > 0 ) obs_vars.write_obs_var_units(obs_units);
      if( obs_descs.n() > 0 ) obs_vars.write_obs_var_descriptions(obs_descs);
      mlog << Debug(7) << method_name << var_count
           << " variable names were saved\n";
   }
   else mlog << Debug(7) << method_name << "no var_names\n";
   //else mlog << Debug(7) << method_name << "use_var_id is false\n";

}

////////////////////////////////////////////////////////////////////////
