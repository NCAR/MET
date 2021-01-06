// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2020
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

#ifndef  __TC_GEN_CONF_INFO_H__
#define  __TC_GEN_CONF_INFO_H__

////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include "vx_statistics.h"
#include "vx_tc_util.h"
#include "vx_config.h"
#include "vx_grid.h"
#include "vx_util.h"

////////////////////////////////////////////////////////////////////////

// Indices for the output flag types in the configuration file

static const int i_fho = 0;
static const int i_ctc = 1;
static const int i_cts = 2;

static const int n_txt = 3;

// Text file type
static const STATLineType txt_file_type[n_txt] = {
   stat_fho,        //  0
   stat_ctc,        //  1
   stat_cts         //  2
};

////////////////////////////////////////////////////////////////////////

struct TCGenNcOutInfo {

   bool do_latlon;
   bool do_fcst_genesis;
   bool do_fcst_fy_oy;
   bool do_fcst_fy_on;
   bool do_fcst_tracks;
   bool do_best_genesis;
   bool do_best_fy_oy;
   bool do_best_fn_oy;
   bool do_best_tracks;

      //////////////////////////////////////////////////////////////////

   TCGenNcOutInfo();
   TCGenNcOutInfo & operator+=(const TCGenNcOutInfo &);

   void clear();   //  sets everything to true

   bool all_false() const;

   void set_all_false();
   void set_all_true();
};

////////////////////////////////////////////////////////////////////////

class TCGenVxOpt {

   private:

      void init_from_scratch();

   public:

      TCGenVxOpt();
     ~TCGenVxOpt();

      //////////////////////////////////////////////////////////////////

      ConcatString Desc;     // Description string
      StringArray  Model;    // Forecast ATCF ID's

      // Analysis track filtering criteria
      StringArray StormId;   // List of storm ids
      StringArray StormName; // List of storm names

      // Timing information
      unixtime  InitBeg, InitEnd;
      unixtime  ValidBeg, ValidEnd;
      NumArray  InitHour;
      NumArray  Lead;

      // Polyline masking region
      ConcatString VxMaskName;
      MaskPoly     VxPolyMask;
      Grid         VxGridMask;
      MaskPlane    VxAreaMask;

      // Distance to land threshold
      SingleThresh DLandThresh;

      // Temporal and spatial matching criteria
      int GenesisSecBeg, GenesisSecEnd;
      double GenesisRadius;
      int GenesisInitDSec;

      // Scoring methods
      bool DevFlag;
      bool OpsFlag;

      // Output file options
      double CIAlpha;
      map<STATLineType,STATOutputType> OutputMap;
      TCGenNcOutInfo NcInfo;
      int NcTrackPtsBeg, NcTrackPtsEnd;

      //////////////////////////////////////////////////////////////////

      void clear();

      void process_config(Dictionary &);
      void parse_nc_info(Dictionary &);

      bool is_keeper(const GenesisInfo &) const;

      STATOutputType output_map(STATLineType) const;
};

////////////////////////////////////////////////////////////////////////

class TCGenConfInfo {

   private:

      void init_from_scratch();

   public:

      TCGenConfInfo();
     ~TCGenConfInfo();

      //////////////////////////////////////////////////////////////////

      // TCPairs configuration object
      MetConfig Conf;

      // Vector of vx task filtering options [n_vx]
      std::vector<TCGenVxOpt> VxOpt;

      // Forecast initialization frequency in hours
      int InitFreqHr;

      // Forecast lead time frequency in hours
      int ValidFreqHr;

      // Begin and end forecast hours for genesis
      int LeadSecBeg, LeadSecEnd;

      // Minimum track duration
      int MinDur;

      // Genesis event criteria
      GenesisEventInfo FcstEventInfo;
      GenesisEventInfo BestEventInfo;
      ConcatString     OperTechnique;

      // Gridded data file containing distances to land
      ConcatString DLandFile;
      Grid         DLandGrid;
      DataPlane    DLandData;

      // Gridded data file containing TC basins
      ConcatString BasinFile;
      Grid         BasinGrid;
      DataPlane    BasinData;

      // Grid for NetCDF output file
      Grid NcOutGrid;

      // Summary of output file options across all filters
      map<STATLineType,STATOutputType> OutputMap;
      TCGenNcOutInfo NcInfo;
   
      // Config file version
      ConcatString Version;

      //////////////////////////////////////////////////////////////////

      void clear();

      void read_config(const char *, const char *);

      void process_config();
      void process_flags(const map<STATLineType,STATOutputType> &,
                         const TCGenNcOutInfo &);

      double compute_dland(double lat, double lon);

      ConcatString compute_basin(double lat, double lon);

      int n_vx() const;
      int compression_level();
   
      STATOutputType output_map(STATLineType) const;
};

////////////////////////////////////////////////////////////////////////

inline int TCGenConfInfo::n_vx()         const { return(VxOpt.size());          }
inline int TCGenConfInfo::compression_level()  { return(Conf.nc_compression()); }

////////////////////////////////////////////////////////////////////////

class GenCTCInfo {

   private:

      void init_from_scratch();

   public:

      GenCTCInfo();
     ~GenCTCInfo();

      //////////////////////////////////////////////////////////////////

   ConcatString Model;
   unixtime FcstBeg, FcstEnd;
   unixtime BestBeg, BestEnd;
   CTSInfo CTSDev, CTSOps;

   const TCGenVxOpt* VxOpt;
   const Grid *NcOutGrid;

   DataPlane FcstGenesisDp, FcstTrackDp;
   DataPlane BestGenesisDp, BestTrackDp;

   DataPlane FcstDevFYOYDp, FcstDevFYONDp;
   DataPlane BestDevFYOYDp, BestDevFNOYDp;

   DataPlane FcstOpsFYOYDp, FcstOpsFYONDp;
   DataPlane BestOpsFYOYDp, BestOpsFNOYDp;

      //////////////////////////////////////////////////////////////////

   void clear();

   void set_vx_opt(const TCGenVxOpt *, const Grid *);

   void inc_dev(bool, bool,
                const GenesisInfo *, const GenesisInfo *);
   void inc_ops(bool, bool,
                const GenesisInfo *, const GenesisInfo *);

   void add_fcst_gen(const GenesisInfo &, int, int);
   void add_best_gen(const GenesisInfo &, int, int);

   void inc_pnt(double, double, DataPlane &);
   void inc_trk(const GenesisInfo &, int, int, DataPlane &);
};

////////////////////////////////////////////////////////////////////////

#endif   /*  __TC_GEN_CONF_INFO_H__  */

////////////////////////////////////////////////////////////////////////
