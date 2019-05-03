// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2015
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

#ifndef  __PB2NC_CONF_INFO_H__
#define  __PB2NC_CONF_INFO_H__

////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "vx_config.h"
#include "vx_analysis_util.h"
#include "vx_grid.h"
#include "vx_util.h"
#include "vx_cal.h"
#include "vx_math.h"

////////////////////////////////////////////////////////////////////////

class PB2NCConfInfo {

   private:

      void init_from_scratch();

   public:

      // PB2NC configuration object
      MetConfig conf;

      // Store data parsed from the PB2NC configuration object
      StringArray  message_type;        // Obseration message type
      StringArray  station_id;          // Observation location station id
      int          beg_ds;              // Time range of observations to be retained,
      int          end_ds;              // defined relative to the PrepBufr center time (seconds)      
      Grid         grid_mask;           // Grid masking region
      MaskPoly     poly_mask;           // Lat/Lon polyline masking region
      double       beg_elev;            // Range of observing location elevations to be retained
      double       end_elev;
      NumArray     pb_report_type;      // PrepBufr report type to be retained
      NumArray     in_report_type;      // Input report types to be retained
      NumArray     instrument_type;     // Instrument types to be retained
      double       beg_level;           // Range of level values to be retained
      double       end_level;
      NumArray     level_category;      // Level categories to be retained
      NumArray     obs_grib_code;       // GRIB1 code for the observation type
      int          quality_mark_thresh; // Quality marks to be retained
      bool         event_stack_flag;    // True for top, false for bottom
      ConcatString tmp_dir;             // Directory for temporary files
      ConcatString version;             // Config file version

      bool         anyair_flag;         // Flags for specific message types
      bool         anysfc_flag;
      bool         onlysf_flag;

      // More information on the PrepBufr file format:
      // http://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc

      PB2NCConfInfo();
     ~PB2NCConfInfo();

      void clear();

      void read_config(const char *, const char *);
      void process_config();
};

////////////////////////////////////////////////////////////////////////

#endif   /*  __PB2NC_CONF_INFO_H__  */

////////////////////////////////////////////////////////////////////////