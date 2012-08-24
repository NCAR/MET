// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2012
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
#include <map>
#include <utility>

#include "tc_stat_job.h"

#include "met_stats.h"
#include "vx_tc_util.h"
#include "vx_log.h"
#include "vx_util.h"
#include "vx_math.h"

////////////////////////////////////////////////////////////////////////

static const char *TCStatJobType_FilterStr  = "filter";
static const char *TCStatJobType_SummaryStr = "summary";

////////////////////////////////////////////////////////////////////////

// Functions for parsing command line options
extern void add_string(const char *, StringArray &);
extern void add_unixtime(const char *, TimeArray &);
extern void add_seconds(const char *, NumArray &);

// Delimiter for separating multiple command line options
static const char *ArgsDelim = ",";

////////////////////////////////////////////////////////////////////////
//
// Code for class TCStatJobFactory
//
////////////////////////////////////////////////////////////////////////

TCStatJob *TCStatJobFactory::new_tc_stat_job_type(const char *type_str) {
   TCStatJob *job = (TCStatJob *) 0;
   TCStatJobType type = NoTCStatJobType;
   
   // Determine the TCStatJobType
   type = string_to_tcstatjobtype(type_str);

   // Switch on job type and instantiate the appropriate class.
   // The TCStatJob object is allocated and needs to be deleted by caller.
   switch(type) {

      case TCStatJobType_Filter:
         job = new TCStatJobFilter;
         break;

      case TCStatJobType_Summary:
         job = new TCStatJobSummary;
         break;         

      case NoTCStatJobType:
      default:
         mlog << Error << "\nTCStatJobFactory::new_tc_stat_job_type() -> "
              << "unsupported job type \"" << type_str << "\"\n\n";
         exit(1);
         break;
   } // end switch

   return(job);
}

///////////////////////////////////////////////////////////////////////////////

TCStatJob *TCStatJobFactory::new_tc_stat_job(const char *jobstring) {
   TCStatJob *job = (TCStatJob *) 0;
   StringArray a;
   ConcatString type_str = na_str;
   ConcatString err_str;
   int i;

   // Parse the jobstring into an array
   a.parse_wsss(jobstring);

   // Check if the job type is specified
   if(a.has("-job", i)) type_str = a[i+1];

   // Allocate a new job of the requested type
   job = new_tc_stat_job_type(type_str);

   // Parse the jobstring into the new job
   a = job->parse_job_command(jobstring);

   // Check for unused arguments
   if(a.n_elements() > 0) {

      // Build list of unknown args
      for(i=0; i<a.n_elements()-1; i++) err_str << a[i] << " ";
      err_str << a[a.n_elements()-1];
     
      mlog << Error
           << "\nTCStatJob *TCStatJobFactory::new_tc_stat_job(const char *jobstring) -> "
           << "unsupported job command options \""
           << err_str << "\".\n\n";
      exit(1);
   }
   
   return(job);
}

////////////////////////////////////////////////////////////////////////
//
// Code for class TCStatJob
//
////////////////////////////////////////////////////////////////////////

TCStatJob::TCStatJob() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

TCStatJob::~TCStatJob() {

   clear();
}

////////////////////////////////////////////////////////////////////////

TCStatJob::TCStatJob(const TCStatJob &j) {

   init_from_scratch();

   assign(j);
}

////////////////////////////////////////////////////////////////////////

TCStatJob & TCStatJob::operator=(const TCStatJob &j) {

   if(this == &j) return(*this);

   assign(j);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::init_from_scratch() {
   
   DumpOut = (ofstream *) 0;
   JobOut  = (ofstream *) 0;

   // Ignore case when performing comparisons
   AModel.set_ignore_case(1);
   BModel.set_ignore_case(1);
   StormId.set_ignore_case(1);
   Basin.set_ignore_case(1);
   Cyclone.set_ignore_case(1);
   StormName.set_ignore_case(1);
   InitMask.set_ignore_case(1);
   ValidMask.set_ignore_case(1);
   LineType.set_ignore_case(1);
   TrackWatchWarn.set_ignore_case(1);
   ColumnThreshName.set_ignore_case(1);
   ColumnStrName.set_ignore_case(1);
   InitThreshName.set_ignore_case(1);
   InitStrName.set_ignore_case(1);

   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::clear() {

   JobType = NoTCStatJobType;

   PairArray.clear();

   AModel.clear();
   BModel.clear();
   StormId.clear();
   Basin.clear();
   Cyclone.clear();
   StormName.clear();
   InitBeg  = InitEnd  = (unixtime) 0;
   InitInc.clear();
   InitExc.clear();
   InitHour.clear();
   Lead.clear();
   ValidBeg = ValidEnd = (unixtime) 0;
   ValidInc.clear();
   ValidExc.clear();
   ValidHour.clear();
   InitMask.clear();
   ValidMask.clear();
   LineType.clear();
   TrackWatchWarn.clear();
   ColumnThreshName.clear();
   ColumnThreshVal.clear();
   ColumnStrName.clear();
   ColumnStrVal.clear();
   InitThreshName.clear();
   InitThreshVal.clear();
   InitStrName.clear();
   InitStrVal.clear();
   
   DumpFile.clear();
   close_dump_file();
   JobOut = (ofstream *) 0;

   // Set to default values
   WaterOnly        = default_water_only;
   RapidInten       = default_rapid_inten;
   RapidIntenThresh = default_rapid_inten_thresh;
   Landfall         = default_landfall;
   LandfallBeg      = default_landfall_beg;
   LandfallEnd      = default_landfall_end;
   MatchPoints      = default_match_points;

   OutInitMask.clear();
   OutValidMask.clear();
   
   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::assign(const TCStatJob & j) {

   clear();
   
   JobType = j.JobType;

   AModel = j.AModel;
   BModel = j.BModel;
   StormId = j.StormId;
   Basin = j.Basin;
   Cyclone = j.Cyclone;
   StormName = j.StormName;
   InitBeg = j.InitBeg;
   InitEnd = j.InitEnd;
   InitInc = j.InitInc;
   InitExc = j.InitExc;
   InitHour = j.InitHour;
   Lead = j.Lead;
   ValidBeg = j.ValidBeg;
   ValidEnd = j.ValidEnd;
   ValidInc = j.ValidInc;
   ValidExc = j.ValidExc;
   ValidHour = j.ValidHour;
   InitMask = j.InitMask;
   ValidMask = j.ValidMask;
   WaterOnly = j.WaterOnly;
   LineType = j.LineType;
   TrackWatchWarn = j.TrackWatchWarn;
   ColumnThreshName = j.ColumnThreshName;
   ColumnThreshVal = j.ColumnThreshVal;
   ColumnStrName = j.ColumnStrName;
   ColumnStrVal = j.ColumnStrVal;
   InitThreshName = j.InitThreshName;
   InitThreshVal = j.InitThreshVal;
   InitStrName = j.InitStrName;
   InitStrVal = j.InitStrVal;

   DumpFile = j.DumpFile;
   open_dump_file();

   RapidInten = j.RapidInten;
   RapidIntenThresh = j.RapidIntenThresh;

   Landfall = j.Landfall;
   LandfallBeg = j.LandfallBeg;
   LandfallEnd = j.LandfallEnd;

   MatchPoints = j.MatchPoints;
   
   OutInitMask = j.OutInitMask;
   OutValidMask = j.OutValidMask;

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::dump(ostream & out, int depth) const {
   Indent prefix(depth);

   out << prefix << "JobType = " << tcstatjobtype_to_string(JobType)
       << "\n";

   out << prefix << "AModel ...\n";
   AModel.dump(out, depth + 1);

   out << prefix << "BModel ...\n";
   BModel.dump(out, depth + 1);

   out << prefix << "StormId ...\n";
   StormId.dump(out, depth + 1);

   out << prefix << "Basin ...\n";
   Basin.dump(out, depth + 1);

   out << prefix << "Cyclone ...\n";
   Cyclone.dump(out, depth + 1);

   out << prefix << "StormName ...\n";
   StormName.dump(out, depth + 1);
   
   out << prefix << "InitBeg = " << unix_to_yyyymmdd_hhmmss(InitBeg) << "\n";
   out << prefix << "InitEnd = " << unix_to_yyyymmdd_hhmmss(InitEnd) << "\n";

   out << prefix << "InitInc ...\n";
   InitInc.dump(out, depth + 1);
   
   out << prefix << "InitExc ...\n";
   InitExc.dump(out, depth + 1);
   
   out << prefix << "InitHour ...\n";
   InitHour.dump(out, depth + 1);
   
   out << prefix << "Lead ...\n";
   Lead.dump(out, depth + 1);   
   
   out << prefix << "ValidBeg = " << unix_to_yyyymmdd_hhmmss(ValidBeg) << "\n";
   out << prefix << "ValidEnd = " << unix_to_yyyymmdd_hhmmss(ValidEnd) << "\n";

   out << prefix << "ValidInc ...\n";
   ValidInc.dump(out, depth + 1);
   
   out << prefix << "ValidExc ...\n";
   ValidExc.dump(out, depth + 1);
   
   out << prefix << "ValidHour ...\n";
   ValidHour.dump(out, depth + 1);
   
   out << prefix << "InitMask ...\n";
   InitMask.dump(out, depth + 1);

   out << prefix << "ValidMask ...\n";
   ValidMask.dump(out, depth + 1);

   out << prefix << "LineType ...\n";
   LineType.dump(out, depth + 1);

   out << prefix << "WaterOnly = " << bool_to_string(WaterOnly) << "\n";
   
   out << prefix << "TrackWatchWarn ...\n";
   TrackWatchWarn.dump(out, depth + 1);
   
   out << prefix << "ColumnThreshName ...\n";
   ColumnThreshName.dump(out, depth + 1);

   out << prefix << "ColumnThreshVal ...\n";
   ColumnThreshVal.dump(out, depth + 1);

   out << prefix << "ColumnStrName ...\n";
   ColumnStrName.dump(out, depth + 1);

   out << prefix << "ColumnStrVal ...\n";
   ColumnStrVal.dump(out, depth + 1);

   out << prefix << "InitThreshName ...\n";
   InitThreshName.dump(out, depth + 1);

   out << prefix << "InitThreshVal ...\n";
   InitThreshVal.dump(out, depth + 1);

   out << prefix << "InitStrName ...\n";
   InitStrName.dump(out, depth + 1);

   out << prefix << "InitStrVal ...\n";
   InitStrVal.dump(out, depth + 1);

   out << prefix << "RapidInten = " << bool_to_string(RapidInten) << "\n";

   out << prefix << "RapidIntenThresh = " << RapidIntenThresh.get_str() << "\n";

   out << prefix << "Landfall = " << bool_to_string(Landfall) << "\n";

   out << prefix << "LandfallBeg = " << LandfallBeg << "\n";

   out << prefix << "LandfallEnd = " << LandfallEnd << "\n";

   out << prefix << "MatchPoints = " << bool_to_string(MatchPoints) << "\n";
   
   out << prefix << "OutInitMask = " << (OutInitMask.name() ? OutInitMask.name() : na_str) << "\n";

   out << prefix << "OutValidMask = " << (OutValidMask.name() ? OutValidMask.name() : na_str) << "\n";

   out << prefix << "DumpFile = " << (DumpFile ? DumpFile.text() : na_str) << "\n";
   
   out.flush();

   return;
}

////////////////////////////////////////////////////////////////////////

bool TCStatJob::is_keeper_track(const TrackPairInfo &tpi,
                                TCLineCounts &n) const {
   bool keep = true;
   int i, j, i_init, offset;
   double v_dbl;
   ConcatString v_str;
   StringArray sa;

   // Check TrackWatchWarn
   if(TrackWatchWarn.n_elements() > 0 &&
      !TrackWatchWarn.has(watchwarntype_to_string(tpi.track_watch_warn())))
      keep = false;
  
   // Check for the special string ALL:
   //    HUWARN, TSWARN, HUWATCH, TSWATCH
   if(TrackWatchWarn.has("ALL") &&
      (tpi.track_watch_warn() == HurricaneWarn     ||
       tpi.track_watch_warn() == TropicalStormWarn ||
       tpi.track_watch_warn() == HurricaneWatch    ||
       tpi.track_watch_warn() == TropicalStormWatch))
       keep = true;

   // Update counts
   if(!keep) n.RejTrackWatchWarn += tpi.n_points();

   // Get the index of the track initialization point
   i_init = tpi.i_init();
   
   // If ADECK initialization filter is requested, check for a valid init time
   if((InitThreshName.n_elements() > 0 ||
       InitStrName.n_elements()    > 0 ||
       OutInitMask.n_points()      > 0) &&
       is_bad_data(i_init)) {
      keep = false;
      n.RejInitThresh += tpi.n_points();
   }
   
   // Check InitThresh
   if(keep == true) {
     
      // Loop through the numeric init column thresholds
      for(i=0; i<InitThreshName.n_elements(); i++) {

         // Get the numeric init column value
         v_dbl = get_column_double(*tpi.line(i_init), InitThreshName[i]);

         // Check the column threshold
         if(is_bad_data(v_dbl) || !InitThreshVal[i].check(v_dbl)) {
           keep = false;
           n.RejInitThresh += tpi.n_points();
           break;
         }
      } // end for i
   }

   // Check InitStr
   if(keep == true) {
  
      // Loop through the column string matching
      for(i=0; i<InitStrName.n_elements(); i++) {

         // Construct list of all entries for the current column name
         sa.clear();
         for(j=0; j<InitStrName.n_elements(); j++)
            if(strcasecmp(InitStrName[i], InitStrName[j]) == 0)
               sa.add(InitStrVal[j]);

         // Determine the column offset and retrieve the value
         offset = determine_column_offset(tpi.line(i_init)->type(), InitStrName[i]);
         v_str  = tpi.line(i_init)->get_item(offset);

         // Check the string value
         if(!sa.has(v_str)) {
            keep = false;
            n.RejInitStr += tpi.n_points();
            break;
         }
      } // end for i
   }

   // Check OutInitMask
   if(keep == true) {

      // A mask has been specified
      if(OutInitMask.n_points() > 0) {

         // Check if ADECK init location falls inside the mask
         if(is_bad_data(tpi.adeck()[i_init].lat()) ||
            is_bad_data(tpi.adeck()[i_init].lon()) ||
            !OutInitMask.latlon_is_inside_dege(
               tpi.adeck()[i_init].lat(),
               tpi.adeck()[i_init].lon())) {
            keep = false;
            n.RejOutInitMask++;
         }
      }
   }

   // Update counts
   if(!keep) n.NKeep -= tpi.n_points();

   return(keep);
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::filter_track(TrackPairInfo &tpi,
                            TCLineCounts &n) const {
   int n_rej;
   
   // Check RapidInten
   if(RapidInten == true) {

      // Determine the rapid intensification points
      n_rej = tpi.check_rapid_inten(RapidIntenThresh);

      // Update counts
      n.RejRapidInten += n_rej;
      n.NKeep         -= n_rej;
   }

   // Check Landfall
   if(Landfall == true) {

      // Determine the landfall points
      n_rej = tpi.check_landfall(LandfallBeg, LandfallEnd);

      // Update counts
      n.RejLandfall += n_rej;
      n.NKeep       -= n_rej;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

bool TCStatJob::is_keeper_line(const TCStatLine &line,
                               TCLineCounts &n) const {
   bool keep = true;
   int i, j, offset;
   double v_dbl, alat, alon, blat, blon, adland, bdland;
   ConcatString v_str;
   StringArray sa;

   // Check TC-STAT header columns
   if(AModel.n_elements() > 0 &&
     !AModel.has(line.amodel()))        { keep = false; n.RejAModel++;    }
   else if(BModel.n_elements() > 0 &&
     !BModel.has(line.bmodel()))        { keep = false; n.RejBModel++;    }
   else if(StormId.n_elements() > 0 &&
     !has_storm_id(StormId, line.basin(), line.cyclone(), line.init()))
                                        { keep = false; n.RejStormId++;   }
   else if(Basin.n_elements() > 0 &&
     !Basin.has(line.basin()))          { keep = false; n.RejBasin++;     }
   else if(Cyclone.n_elements() > 0 &&
     !Cyclone.has(line.cyclone()))      { keep = false; n.RejCyclone++;   }
   else if(StormName.n_elements() > 0 &&
     !StormName.has(line.storm_name())) { keep = false; n.RejStormName++; }
   else if(InitBeg > 0 &&
      line.init() < InitBeg)            { keep = false; n.RejInit++;      }
   else if(InitEnd > 0 &&
      line.init() > InitEnd)            { keep = false; n.RejInit++;      }
   else if(InitInc.n_elements() > 0 &&
     !InitInc.has(line.init()))         { keep = false; n.RejInit++;      }
   else if(InitExc.n_elements() > 0 &&
     InitExc.has(line.init()))          { keep = false; n.RejInit++;      }
   else if(InitHour.n_elements() > 0 &&
     !InitHour.has(line.init_hour()))   { keep = false; n.RejInitHour++;  }
   else if(Lead.n_elements() > 0 &&
     !Lead.has(line.lead()))            { keep = false; n.RejLead++;      }
   else if(ValidBeg > 0 &&
      line.valid() < ValidBeg)          { keep = false; n.RejValid++;     }
   else if(ValidEnd > 0 &&
      line.valid() > ValidEnd)          { keep = false; n.RejValid++;     }
   else if(ValidInc.n_elements() > 0 &&
     !ValidInc.has(line.valid()))       { keep = false; n.RejValid++;     }
   else if(ValidExc.n_elements() > 0 &&
     ValidExc.has(line.valid()))        { keep = false; n.RejValid++;     }
   else if(ValidHour.n_elements() > 0 &&
     !ValidHour.has(line.valid_hour())) { keep = false; n.RejValidHour++; }
   else if(InitMask.n_elements() > 0 &&
     !InitMask.has(line.init_mask()))   { keep = false; n.RejInitMask++;  }
   else if(ValidMask.n_elements() > 0 &&
     !ValidMask.has(line.valid_mask())) { keep = false; n.RejValidMask++; }
   else if(LineType.n_elements() > 0 &&
     !LineType.has(line.line_type()))   { keep = false; n.RejLineType++;  }

   // Check WaterOnly
   if(keep == true && WaterOnly == true) {
      adland = atof(line.get_item("ADLAND"));
      bdland = atof(line.get_item("BDLAND"));

      if(is_bad_data(adland) || adland <= 0 ||
         is_bad_data(bdland) || bdland <= 0) {
         keep = false;
         n.RejWaterOnly++;
      }
   }

   // Check ColumnThresh
   if(keep == true) {

      // Loop through the numeric column thresholds
      for(i=0; i<ColumnThreshName.n_elements(); i++) {

         // Get the numeric column value
         v_dbl = get_column_double(line, ColumnThreshName[i]);

         // Check the column threshold
         if(is_bad_data(v_dbl) || !ColumnThreshVal[i].check(v_dbl)) {
           keep = false;
           n.RejColumnThresh++;
           break;
         }
      } // end for i
   }

   // Check ColumnStr
   if(keep == true) {

      // Loop through the column string matching
      for(i=0; i<ColumnStrName.n_elements(); i++) {

         // Construct list of all entries for the current column name
         sa.clear();
         for(j=0; j<ColumnStrName.n_elements(); j++)
            if(strcasecmp(ColumnStrName[i], ColumnStrName[j]) == 0)
               sa.add(ColumnStrVal[j]);

         // Determine the column offset and retrieve the value
         offset = determine_column_offset(line.type(), ColumnStrName[i]);
         v_str  = line.get_item(offset);

         // Check the string value
         if(!sa.has(v_str)) {
            keep = false;
            n.RejColumnStr++;
            break;
         }
      } // end for i
   }

   // Retrieve the ADECK and BDECK lat/lon values
   alat = atof(line.get_item("ALAT"));
   alon = atof(line.get_item("ALON"));
   blat = atof(line.get_item("BLAT"));
   blon = atof(line.get_item("BLON"));

   // Check MatchPoints
   if(keep == true && MatchPoints == true) {

      // Check that the ADECK and BDECK locations are both defined
      if(is_bad_data(alat) || is_bad_data(alon) ||
         is_bad_data(blat) || is_bad_data(blon)) {
         keep = false;
         n.RejMatchPoints++;
      }
   }
   
   // Check OutValidMask
   if(keep == true && OutValidMask.n_points() > 0) {

      // Check ADECK locations
      if(is_bad_data(alat) ||
         is_bad_data(alon) ||
         !OutValidMask.latlon_is_inside_dege(alat, alon)) {
         keep = false;
         n.RejOutValidMask++;
      }
   }

   return(keep);
}

////////////////////////////////////////////////////////////////////////

double TCStatJob::get_column_double(const TCStatLine &line,
                                    const ConcatString &column) const {
   StringArray sa;
   ConcatString in;
   double v, v_cur;
   bool abs_flag = false;
   int i;

   // Check for absolute value
   if(strncasecmp(column, "ABS", 3) == 0) {
      abs_flag = true;
      sa = column.split("()");
      in = sa[1];
   }
   else {
      in = column;
   }

   // Split the input column name on hyphens for differences
   sa = in.split("-");

   // Get the first value
   v = atof(line.get_item(sa[0]));

   // If multiple columns, compute the requested difference
   if(sa.n_elements() > 1) {

      // Loop through the column
      for(i=1; i<sa.n_elements(); i++) {

         // Get the current column value
         v_cur = atof(line.get_item(sa[i]));

         // Compute the difference, checking for bad data
         if(is_bad_data(v) || is_bad_data(v_cur)) v  = bad_data_double;
         else                                     v -= v_cur;
      }
   }

   // Apply absolute value, if requested
   if(abs_flag && !is_bad_data(v)) v = fabs(v);

   return(v);
}

////////////////////////////////////////////////////////////////////////

StringArray TCStatJob::parse_job_command(const char *jobstring) {
   StringArray a, b;
   const char * c = (const char *) 0;
   int i;

   // Parse the jobstring into a StringArray
   a.parse_wsss(jobstring);

   // Loop over the StringArray elements
   for(i=0; i<a.n_elements(); i++) {

      c = a[i];

      // Check for a job command option
      if(c[0] != '-') {
         b.add(a[i]);
         continue;
      }

      // Check job command options
           if(strcasecmp(c, "-job"               ) == 0) { JobType = string_to_tcstatjobtype(a[i+1]); a.shift_down(i, 1); }
      else if(strcasecmp(c, "-amodel"            ) == 0) { add_string(a[i+1], AModel);                a.shift_down(i, 1); }
      else if(strcasecmp(c, "-bmodel"            ) == 0) { add_string(a[i+1], BModel);                a.shift_down(i, 1); }
      else if(strcasecmp(c, "-storm_id"          ) == 0) { add_string(a[i+1], StormId);               a.shift_down(i, 1); }
      else if(strcasecmp(c, "-basin"             ) == 0) { add_string(a[i+1], Basin);                 a.shift_down(i, 1); }
      else if(strcasecmp(c, "-cyclone"           ) == 0) { add_string(a[i+1], Cyclone);               a.shift_down(i, 1); }
      else if(strcasecmp(c, "-storm_name"        ) == 0) { add_string(a[i+1], StormName);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_beg"          ) == 0) { InitBeg = timestring_to_unix(a[i+1]);      a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_end"          ) == 0) { InitEnd = timestring_to_unix(a[i+1]);      a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_inc"          ) == 0) { add_unixtime(a[i+1], InitInc);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_exc"          ) == 0) { add_unixtime(a[i+1], InitExc);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_hour"         ) == 0) { add_seconds(a[i+1], InitHour);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-lead"              ) == 0) { add_seconds(a[i+1], Lead);                 a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_beg"         ) == 0) { ValidBeg = timestring_to_unix(a[i+1]);     a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_end"         ) == 0) { ValidEnd = timestring_to_unix(a[i+1]);     a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_inc"         ) == 0) { add_unixtime(a[i+1], ValidInc);            a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_exc"         ) == 0) { add_unixtime(a[i+1], ValidExc);            a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_hour"        ) == 0) { add_seconds(a[i+1], ValidHour);            a.shift_down(i, 1); }
      else if(strcasecmp(c, "-init_mask"         ) == 0) { add_string(a[i+1], InitMask);              a.shift_down(i, 1); }
      else if(strcasecmp(c, "-valid_mask"        ) == 0) { add_string(a[i+1], ValidMask);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-line_type"         ) == 0) { add_string(a[i+1], LineType);              a.shift_down(i, 1); }
      else if(strcasecmp(c, "-water_only"        ) == 0) { WaterOnly = string_to_bool(a[i+1]);        a.shift_down(i, 1); }
      else if(strcasecmp(c, "-track_watch_warn"  ) == 0) { add_string(a[i+1], TrackWatchWarn);        a.shift_down(i, 1); }
      else if(strcasecmp(c, "-column_thresh"     ) == 0) { ColumnThreshName.add(a[i+1]);
                                                           ColumnThreshVal.add(a[i+2]);               a.shift_down(i, 2); }
      else if(strcasecmp(c, "-column_str"        ) == 0) { ColumnStrName.add(a[i+1]);
                                                           ColumnStrVal.add(a[i+2]);                  a.shift_down(i, 2); }
      else if(strcasecmp(c, "-init_thresh"       ) == 0) { InitThreshName.add(a[i+1]);
                                                           InitThreshVal.add(a[i+2]);                 a.shift_down(i, 2); }
      else if(strcasecmp(c, "-init_str"          ) == 0) { InitStrName.add(a[i+1]);
                                                           InitStrVal.add(a[i+2]);                    a.shift_down(i, 2); }
      else if(strcasecmp(c, "-rapid_inten"       ) == 0) { RapidInten = string_to_bool(a[i+1]);       a.shift_down(i, 1); }
      else if(strcasecmp(c, "-rapid_inten_thresh") == 0) { RapidIntenThresh.set(a[i+1]);              a.shift_down(i, 1); }
      else if(strcasecmp(c, "-landfall"          ) == 0) { Landfall = string_to_bool(a[i+1]);         a.shift_down(i, 1); }
      else if(strcasecmp(c, "-landfall_beg"      ) == 0) { LandfallBeg = atoi(a[i+1]);                a.shift_down(i, 1); }
      else if(strcasecmp(c, "-landfall_end"      ) == 0) { LandfallEnd = atoi(a[i+1]);                a.shift_down(i, 1); }
      else if(strcasecmp(c, "-match_points"      ) == 0) { MatchPoints = string_to_bool(a[i+1]);      a.shift_down(i, 1); }
      else if(strcasecmp(c, "-out_init_mask"     ) == 0) { set_mask(OutInitMask, a[i+1]);             a.shift_down(i, 1); }
      else if(strcasecmp(c, "-out_valid_mask"    ) == 0) { set_mask(OutValidMask, a[i+1]);            a.shift_down(i, 1); }
      else if(strcasecmp(c, "-dump_row"          ) == 0) { DumpFile = a[i+1]; open_dump_file();       a.shift_down(i, 1); }
      else                                               {                                            b.add(a[i]);        }
   }

   return(b);
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::set_mask(MaskPoly &p, const char *file) {
  
   ConcatString poly_file = replace_path(file);
   p.clear();
   p.load(poly_file);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::open_dump_file() {

   close_dump_file();

   if(DumpFile.empty()) return;

   DumpOut = new ofstream;
   DumpOut->open(DumpFile);

   if(!DumpOut) {
      mlog << Error << "\nTCStatJob::open_dump_file()-> "
           << "can't open the output file \"" << DumpFile
           << "\" for writing!\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::close_dump_file() {

   if(DumpOut) {

      mlog << Debug(1)
           << "Creating output dump file: " << DumpFile << "\n";

      DumpOut->close();
      delete DumpOut;
      DumpOut = (ofstream *) 0;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::write_dump_file() {
  
   if(!DumpOut) return;
   
   int i, i_row;
   TcHdrColumns tchc;
   AsciiTable out_at;
        
   // Initialize the output AsciiTable
   out_at.set_size(PairArray.n_points() + 1,
                   n_tc_header_cols + n_tc_mpr_cols);

   // Setup the output AsciiTable
   out_at.set_table_just(LeftJust);
   out_at.set_precision(default_precision);
   out_at.set_bad_data_value(bad_data_double);
   out_at.set_bad_data_str(na_str);
   out_at.set_delete_trailing_blank_rows(1);

   // Write the TCMPR header row
   write_tc_mpr_header_row(1, out_at, 0, 0);
   
   // Loop over and write the pairs
   for(i=0, i_row=1; i<PairArray.n_pairs(); i++) {

      // Check for no remaining points
      if(PairArray[i].n_points() == 0) continue;
     
      // Setup header columns
      tchc.clear();
      tchc.set_adeck_model(PairArray[i].adeck().technique());
      tchc.set_bdeck_model(PairArray[i].bdeck().technique());
      tchc.set_basin(PairArray[i].bdeck().basin());
      tchc.set_cyclone(PairArray[i].bdeck().cyclone());
      tchc.set_storm_name(PairArray[i].bdeck().storm_name());
   
      if(OutInitMask.n_points() > 0)  tchc.set_init_mask(OutInitMask.name());
      else                            tchc.set_init_mask(na_str);
      if(OutValidMask.n_points() > 0) tchc.set_valid_mask(OutValidMask.name());
      else                            tchc.set_valid_mask(na_str);

      // Write the TrackPairInfo object
      write_tc_mpr_row(tchc, PairArray[i], out_at, i_row);
   }

   // Write the AsciiTable to the file
   *DumpOut << out_at;

   // Close the dump file
   close_dump_file();
   
   return;
}

////////////////////////////////////////////////////////////////////////

ConcatString TCStatJob::serialize() const {
   ConcatString s;
   int i;

   // Initialize the jobstring
   s.clear();
   s.set_precision(default_precision);

   if(JobType != NoTCStatJobType)
      s << "-job " << tcstatjobtype_to_string(JobType) << " ";
   for(i=0; i<AModel.n_elements(); i++)
      s << "-amodel " << AModel[i] << " ";
   for(i=0; i<BModel.n_elements(); i++)
      s << "-bmodel " << BModel[i] << " ";
   for(i=0; i<StormId.n_elements(); i++)
      s << "-storm_id " << StormId[i] << " ";
   for(i=0; i<Basin.n_elements(); i++)
      s << "-basin " << Basin[i] << " ";
   for(i=0; i<Cyclone.n_elements(); i++)
      s << "-cyclone " << Cyclone[i] << " ";
   for(i=0; i<StormName.n_elements(); i++)
      s << "-storm_name " << StormName[i] << " ";
   if(InitBeg > 0)
      s << "-init_beg " << unix_to_yyyymmdd_hhmmss(InitBeg) << " ";
   if(InitEnd > 0)
      s << "-init_end " << unix_to_yyyymmdd_hhmmss(InitEnd) << " ";
   for(i=0; i<InitInc.n_elements(); i++)
      s << "-init_inc " << unix_to_yyyymmdd_hhmmss(InitInc[i]) << " ";
   for(i=0; i<InitExc.n_elements(); i++)
      s << "-init_exc " << unix_to_yyyymmdd_hhmmss(InitExc[i]) << " ";
   for(i=0; i<InitHour.n_elements(); i++)
      s << "-init_hour " << sec_to_hhmmss(InitHour[i]) << " ";
   for(i=0; i<Lead.n_elements(); i++)
      s << "-lead " << sec_to_hhmmss(Lead[i]) << " ";
   if(ValidBeg > 0)
      s << "-valid_beg " << unix_to_yyyymmdd_hhmmss(ValidBeg) << " ";
   if(ValidEnd > 0)
      s << "-valid_end " << unix_to_yyyymmdd_hhmmss(ValidEnd) << " ";
   for(i=0; i<ValidInc.n_elements(); i++)
      s << "-valid_inc " << unix_to_yyyymmdd_hhmmss(ValidInc[i]) << " ";
   for(i=0; i<ValidExc.n_elements(); i++)
      s << "-valid_exc " << unix_to_yyyymmdd_hhmmss(ValidExc[i]) << " ";
   for(i=0; i<ValidHour.n_elements(); i++)
      s << "-valid_hour " << sec_to_hhmmss(ValidHour[i]) << " ";
   for(i=0; i<InitMask.n_elements(); i++)
      s << "-init_mask " << InitMask[i] << " ";
   for(i=0; i<ValidMask.n_elements(); i++)
      s << "-valid_mask " << ValidMask[i] << " ";
   for(i=0; i<LineType.n_elements(); i++)
      s << "-line_type " << LineType[i] << " ";
   if(WaterOnly != default_water_only)
      s << "-water_only " << bool_to_string(WaterOnly) << " ";
   for(i=0; i<TrackWatchWarn.n_elements(); i++)
      s << "-track_watch_warn " << TrackWatchWarn[i] << " ";
   for(i=0; i<ColumnThreshName.n_elements(); i++)
      s << "-column_thresh " << ColumnThreshName[i] << " "
                             << ColumnThreshVal[i].get_str() << " ";
   for(i=0; i<ColumnStrName.n_elements(); i++)
      s << "-column_str " << ColumnStrName[i] << " "
                          << ColumnStrVal[i] << " ";
   for(i=0; i<InitThreshName.n_elements(); i++)
      s << "-init_thresh " << InitThreshName[i] << " "
                           << InitThreshVal[i].get_str() << " ";
   for(i=0; i<InitStrName.n_elements(); i++)
      s << "-init_str " << InitStrName[i] << " "
                        << InitStrVal[i] << " ";
   if(RapidInten != default_rapid_inten)
      s << "-rapid_inten " << bool_to_string(RapidInten) << " ";
   if(!(RapidIntenThresh == default_rapid_inten_thresh))
      s << "-rapid_inten_thresh " << RapidIntenThresh.get_str() << " ";
   if(Landfall != default_landfall)
      s << "-landfall " << bool_to_string(Landfall) << " ";
   if(LandfallBeg != default_landfall_beg ||
      LandfallEnd != default_landfall_end)
      s << "-landfall_beg " << LandfallBeg << " "
        << "-landfall_end " << LandfallEnd << " ";
   if(MatchPoints != default_match_points)
      s << "-match_points " << bool_to_string(MatchPoints) << " ";
   if(OutInitMask.n_points() > 0)
      s << "-out_init_mask " << OutInitMask.file_name() << " ";
   if(OutValidMask.n_points() > 0)
      s << "-out_valid_mask " << OutValidMask.file_name() << " ";
   if(DumpFile.length() > 0)
      s << "-dump_row " << DumpFile << " ";
   
   return(s);
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::do_job(const StringArray &file_list,
                       TCLineCounts &n) {
   int i;
   
   // Loop through the input file list
   for(i=0; i<file_list.n_elements(); i++) {
      process_tc_stat_file(file_list[i], n);
   }

   // Write the dump file
   if(DumpOut) write_dump_file();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::process_tc_stat_file(const char *path, TCLineCounts &n) {
   LineDataFile in;
   TCStatLine line;
   TrackPairInfo tpi;
   int index;

   // Open the data file
   if(!(in.open(path))) {
      mlog << Error << "\nTCStatJob::process_tc_stat_file(const char *path) -> "
           << "can't open file \"" << path << "\" for reading\n\n";
      exit(1);
   }

   mlog << Debug(3)
        << "Reading file: " << path << "\n";

   // Process each TC-STAT line
   while(in >> line) {

      // Increment counts
      n.NRead++;
      n.NKeep++;

      // Get the index value for the current line
      index = atoi(line.get_item("INDEX"));

      // Check for the first point of a track pair
      if(index == 1) {
         process_track_pair(tpi, n);
         tpi.clear();
      }

      // Add the current line to the track pair
      tpi.add(line);
   }

   // Store the last track pair
   process_track_pair(tpi, n);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJob::process_track_pair(TrackPairInfo &tpi, TCLineCounts &n) {
   int i;
   TrackPairInfo tpi_subset;

   // Check for no points
   if(tpi.n_points() == 0) return;

   // Determine if the whole track can be discarded
   if(!is_keeper_track(tpi, n)) return;

   // Filter down the track if necessary
   filter_track(tpi, n);

   // Loop over the track points to check line by line
   for(i=0; i<tpi.n_lines(); i++) {

      // Continue is this line is not be kept
      if(!tpi.keep(i)) continue;

      // Check if this line is not to be kept
      if(!is_keeper_line(*tpi.line(i), n)) {
         tpi.set_keep(i, 0);
         n.NKeep -= 1;
      }
   }

   // Retrieve the subset of track pair points marked for retention
   tpi_subset = tpi.keep_subset();
   
   // Add track pair to the PairArray
   if(tpi_subset.n_points() > 0) PairArray.add(tpi_subset);
   
   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class TCStatJobFilter
//
////////////////////////////////////////////////////////////////////////

TCStatJobFilter::TCStatJobFilter() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

TCStatJobFilter::~TCStatJobFilter() {

   clear();
}

////////////////////////////////////////////////////////////////////////

TCStatJobFilter::TCStatJobFilter(const TCStatJobFilter &j) {

   init_from_scratch();

   assign(j);
}

////////////////////////////////////////////////////////////////////////

TCStatJobFilter & TCStatJobFilter::operator=(const TCStatJobFilter &j) {

   if(this == &j) return(*this);

   assign(j);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void TCStatJobFilter::init_from_scratch() {

   TCStatJob::init_from_scratch();
  
   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobFilter::clear() {

   TCStatJob::clear();
   
   JobType = TCStatJobType_Filter;

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobFilter::assign(const TCStatJobFilter & j) {

   TCStatJob::assign(j);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobFilter::do_job(const StringArray &file_list,
                             TCLineCounts &n) {

   // Check that the -dump_row option has been supplied
   if(!DumpOut) {
      mlog << Error << "\nTCStatJobFilter::do_job() -> "
           << "this function may only be called when using the "
           << "-dump_row option in the job command line: "
           << serialize() << "\n\n";
      exit(1);
   }

   // Call the parent's do_job() to filter the data
   TCStatJob::do_job(file_list, n);
   
   // Process the filter output
   if(JobOut) do_output(*JobOut);
   else       do_output(cout);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobFilter::do_output(ostream &out) {
   ConcatString line;

   // Build a simple output line
   line << "FILTER: " << serialize() << "\n";
   out  << line << "\n";

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for class TCStatJobSummary
//
////////////////////////////////////////////////////////////////////////

TCStatJobSummary::TCStatJobSummary() {

   init_from_scratch();
}

////////////////////////////////////////////////////////////////////////

TCStatJobSummary::~TCStatJobSummary() {

   clear();
}

////////////////////////////////////////////////////////////////////////

TCStatJobSummary::TCStatJobSummary(const TCStatJobSummary &j) {

   init_from_scratch();

   assign(j);
}

////////////////////////////////////////////////////////////////////////

TCStatJobSummary & TCStatJobSummary::operator=(const TCStatJobSummary &j) {

   if(this == &j) return(*this);

   assign(j);

   return(*this);
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::init_from_scratch() {

   TCStatJob::init_from_scratch();
   
   // Ignore case when performing comparisons
   ReqColumn.set_ignore_case(1);
   Column.set_ignore_case(1);
   Case.set_ignore_case(1);
   
   clear();

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::clear() {

   TCStatJob::clear();

   JobType = TCStatJobType_Summary;

   ReqColumn.clear();
   Column.clear();
   Case.clear();
   SummaryMap.clear();

   // Set to default value
   OutAlpha  = default_tc_alpha;
   FSPThresh = default_fsp_thresh;
   
   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::assign(const TCStatJobSummary & j) {

   TCStatJob::assign(j);

   ReqColumn = j.ReqColumn;
   Column = j.Column;
   Case = j.Case;
   SummaryMap = j.SummaryMap;
   OutAlpha = j.OutAlpha;
   FSPThresh = j.FSPThresh;

   return;
}

////////////////////////////////////////////////////////////////////////

StringArray TCStatJobSummary::parse_job_command(const char *jobstring) {
   StringArray a, b;
   const char * c = (const char *) 0;
   int i;

   // Call the parent and store any unused options
   a = TCStatJob::parse_job_command(jobstring);

   // Loop over the StringArray elements
   for(i=0; i<a.n_elements(); i++) {

      // Point at the current entry
      c = a[i];
      
      // Check for a job command option
      if(c[0] != '-') {
         b.add(a[i]);
         continue;
      }

      // Check job command options
           if(strcasecmp(c, "-column"    ) == 0) { add_string(a[i+1], ReqColumn);
                                                   add_column(a[i+1]);            a.shift_down(i, 1); }
      else if(strcasecmp(c, "-by"        ) == 0) { add_string(a[i+1], Case);      a.shift_down(i, 1); }
      else if(strcasecmp(c, "-out_alpha" ) == 0) { OutAlpha = atof(a[i+1]);       a.shift_down(i, 1); }
      else if(strcasecmp(c, "-fsp_thresh") == 0) { FSPThresh.set(a[i+1]);         a.shift_down(i, 1); }
      else                                       {                                b.add(a[i]);        }
   }

   return(b);
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::add_column(const char *c) {
   int i, j;
   StringArray sa;

   // Parse entries into a StringArray
   add_string(c, sa);

   // Loop over the entries, handling special column names
   for(i=0; i<sa.n_elements(); i++) {

      // Track errors
      if(strcasecmp(sa[i], "TRACK") == 0) {
         for(j=0; j<n_tc_cols_track; j++) Column.add(tc_cols_track[j]);
      }
      // Wind errors
      else if(strcasecmp(sa[i], "WIND") == 0) {
         for(j=0; j<n_tc_cols_wind; j++) Column.add(tc_cols_wind[j]);
      }
      // Track and Intensity (TI)
      else if(strcasecmp(sa[i], "TI") == 0) {
         for(j=0; j<n_tc_cols_ti; j++) Column.add(tc_cols_ti[j]);
      }
      // Along and Cross Track (AC)
      else if(strcasecmp(sa[i], "AC") == 0) {
         for(j=0; j<n_tc_cols_ac; j++) Column.add(tc_cols_ac[j]);
      }
      // Lat/Lon Difference (XY)
      else if(strcasecmp(sa[i], "XY") == 0) {
         for(j=0; j<n_tc_cols_xy; j++) Column.add(tc_cols_xy[j]);
      }
      // Otherwise, just add the column name
      else {
         Column.add(sa[i]);
      }
   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

ConcatString TCStatJobSummary::serialize() const {
   ConcatString s;
   int i;

   // Call the parent
   s = TCStatJob::serialize();

   // Add summary job-specific options
   for(i=0; i<ReqColumn.n_elements(); i++)
      s << "-column " << ReqColumn[i] << " ";
   for(i=0; i<Case.n_elements(); i++)
      s << "-by " << Case[i] << " ";
   if(OutAlpha != default_tc_alpha)
      s << "-out_alpha " << OutAlpha << " ";
   if(!(FSPThresh == default_fsp_thresh))
      s << "-fsp_thresh " << FSPThresh.get_str();

   return(s);
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::do_job(const StringArray &file_list,
                              TCLineCounts &n) {
   int i;

   // Check that the -column option has been supplied
   if(Column.n_elements() == 0) {
      mlog << Error << "\nTCStatJobSummary::do_job() -> "
           << "this function may only be called when using the "
           << "-column option in the job command line: "
           << serialize() << "\n\n";
      exit(1);
   }

   // Loop through the input file list
   for(i=0; i<file_list.n_elements(); i++) {
      process_tc_stat_file(file_list[i], n);
   }

   // Process the PairArray entries
   process_pair_array();

   // Write the dump file
   if(DumpOut) write_dump_file();

   // Process the summary output
   if(JobOut) do_output(*JobOut);
   else       do_output(cout);

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::process_pair_array() {
   int i, j, k, l, lead;
   map<ConcatString,MapData,cs_cmp> cur_map;
   map<ConcatString,MapData,cs_cmp>::iterator it;
   ConcatString key, cur;
   MapData data;
   double val;

   // Loop over the TrackPairInfoArray entries
   for(i=0; i<PairArray.n_pairs(); i++) {

      // Initialize the map
      cur_map.clear();
     
      // Loop over TCStatLines for each TrackPairInfo entry
      // and construct a map
      for(j=0; j<PairArray[i].n_lines(); j++) {
        
         // Add summary info to the current map
         for(k=0; k<Column.n_elements(); k++) {

            // Build the key and get the current column value
            key = Column[k];
            val = get_column_double(*PairArray[i].line(j), Column[k]);

            // Add case information to the key
            for(l=0; l<Case.n_elements(); l++) {

               cur = PairArray[i].line(j)->get_item(Case[l]);

               // For bad data, use the NA string
               if(is_bad_data(atoi(cur))) cur = na_str;

               // Special handling for lead time:
               // Switch 2-digit hours to 3-digit hours so that the
               // summary job output is sorted nicely.
               if(strcasecmp(Case[l], "LEAD") == 0 && cur != na_str &&
                  abs(lead = hhmmss_to_sec(cur)) < 100*sec_per_hour) {

                  // Handle positive and negative lead times
                  key << (lead < 0 ? ":-0" : ":0")
                      << sec_to_hhmmss(abs(lead));
               }
               // Otherwise, just append the current case info
               else {
                  key << ":" << cur;
               }
               
            } // end for l

            // Key is not yet defined in the map
            if((it = cur_map.find(key)) == cur_map.end()) {

               // Add the MapData to the current map
               clear_map_data(data);
               data.Val.add(val);
               data.Hdr.add(PairArray[i].line(j)->header());
               data.AModel.add(PairArray[i].line(j)->amodel());
               data.Init.add(PairArray[i].line(j)->init());
               data.Lead.add(PairArray[i].line(j)->lead());
               data.Valid.add(PairArray[i].line(j)->valid());
               cur_map[key] = data;
            }
            // Key is defined in the map
            else {

               // Add the value for the existing key
               it->second.Val.add(val);
               it->second.Hdr.add(PairArray[i].line(j)->header());
               it->second.AModel.add(PairArray[i].line(j)->amodel());
               it->second.Init.add(PairArray[i].line(j)->init());
               it->second.Lead.add(PairArray[i].line(j)->lead());
               it->second.Valid.add(PairArray[i].line(j)->valid());
            }
            
         } // end for k
      } // end for j

      // Add the current map
      add_map(cur_map);

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::add_map(map<ConcatString,MapData,cs_cmp>&m) {
   map<ConcatString,MapData,cs_cmp>::iterator it_in;
   map<ConcatString,MapData,cs_cmp>::iterator it;

   // Loop over the input map entries
   for(it_in = m.begin(); it_in != m.end(); it_in++) {

      // Current map key is not yet defined in SummaryMap
      if((it = SummaryMap.find(it_in->first)) == SummaryMap.end()) {

         mlog << Debug(5)
              << "Summary Map Insert (" << it_in->first << ") "
              << it_in->second.Val.n_elements() << " values: "
              << it_in->second.Val.serialize() << "\n";

         // Add the pair to the map
         SummaryMap[it_in->first] = it_in->second;
      }
      // Current map key is defined in SummaryMap
      else {

         mlog << Debug(5)
              << "Summary Map Add (" << it_in->first << ") "
              << it_in->second.Val.n_elements() << " values: "
              << it_in->second.Val.serialize() << "\n";

         // Add the value for the existing key
         it->second.Val.add(it_in->second.Val);
         it->second.Hdr.add(it_in->second.Hdr);
         it->second.AModel.add(it_in->second.AModel);
         it->second.Init.add(it_in->second.Init);
         it->second.Lead.add(it_in->second.Lead);
         it->second.Valid.add(it_in->second.Valid);
      }
   } // end for it_in

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::do_output(ostream &out) {
   map<ConcatString,MapData,cs_cmp>::iterator it;
   StringArray sa;
   ConcatString case_info, line;
   AsciiTable out_at;
   NumArray v, lead, index;
   TimeArray init, valid;
   CIInfo mean_ci, stdev_ci;
   int i, r, c, dsec, tind;

   double fsp;
   NumArray fsp_total, fsp_best, fsp_ties;

   // Setup the output table
   out_at.set_size((int) SummaryMap.size() + 1, 23);
   out_at.set_table_just(LeftJust);
   out_at.set_precision(default_precision);
   out_at.set_bad_data_value(bad_data_double);
   out_at.set_bad_data_str(na_str);
   out_at.set_delete_trailing_blank_rows(1);

   // Set up CIInfo objects
   mean_ci.allocate_n_alpha(1);
   stdev_ci.allocate_n_alpha(1);
   
   // Initialize row and column indices
   r = c = 0;

   // Write the header row
   out_at.set_entry(r, c++, "COL_NAME:");
   out_at.set_entry(r, c++, "COLUMN");
   out_at.set_entry(r, c++, "CASE");
   out_at.set_entry(r, c++, "TOTAL");
   out_at.set_entry(r, c++, "VALID");
   out_at.set_entry(r, c++, "MEAN");
   out_at.set_entry(r, c++, "MEAN_NCL");
   out_at.set_entry(r, c++, "MEAN_NCU");
   out_at.set_entry(r, c++, "STDEV");
   out_at.set_entry(r, c++, "MIN");   
   out_at.set_entry(r, c++, "P10");
   out_at.set_entry(r, c++, "P25");
   out_at.set_entry(r, c++, "P50");
   out_at.set_entry(r, c++, "P75");
   out_at.set_entry(r, c++, "P90");
   out_at.set_entry(r, c++, "MAX");
   out_at.set_entry(r, c++, "SUM");
   out_at.set_entry(r, c++, "TS_INT");
   out_at.set_entry(r, c++, "TS_IND");
   out_at.set_entry(r, c++, "FSP_TOTAL");
   out_at.set_entry(r, c++, "FSP_BEST");
   out_at.set_entry(r, c++, "FSP_TIES");
   out_at.set_entry(r, c++, "FSP");

   // Compute the frequency of superior performance
   compute_fsp(fsp_total, fsp_best, fsp_ties);
   
   // Loop over the map entries and popluate the output table
   for(it=SummaryMap.begin(),r=1; it!=SummaryMap.end(); it++,r++) {

      // Initialize column index
      c = 0;

      // Split the current map key
      sa = it->first.split(":");

      // Rebuild the case info
      if(sa.n_elements() == 1) {
         case_info = "NA";
      }
      else {
         case_info = sa[1];
         for(i=2; i<sa.n_elements(); i++) case_info << ":" << sa[i];
      }

      // Get the valid subset of data
      v.clear();
      init.clear();
      lead.clear();
      valid.clear();
      for(i=0; i<it->second.Val.n_elements(); i++) {
         if(!is_bad_data(it->second.Val[i])) {
            v.add(it->second.Val[i]);
            init.add(it->second.Init[i]);
            lead.add(it->second.Lead[i]);
            valid.add(it->second.Valid[i]);
         }
      }

      // Build index array
      index.clear();
      for(i=0; i<v.n_elements(); i++) index.add(i);

      // Compute mean and standard deviation
      compute_mean_stdev(v, index, 1, OutAlpha, mean_ci, stdev_ci);

      // Compute the FSP value
      if(fsp_total[r-1] != 0) fsp = fsp_best[r-1]/fsp_total[r-1];
      else                    fsp = bad_data_double;

      // Compute time to independence for time-series data
      if(is_time_series(init, lead, valid, dsec))
         tind = compute_time_to_indep(v, dsec);
      else
         tind = dsec = bad_data_double;
                                
      // Write the table row
      out_at.set_entry(r, c++, "SUMMARY:");      
      out_at.set_entry(r, c++, sa[0]);
      out_at.set_entry(r, c++, case_info);
      out_at.set_entry(r, c++, it->second.Val.n_elements());
      out_at.set_entry(r, c++, v.n_elements());
      out_at.set_entry(r, c++, mean_ci.v);
      out_at.set_entry(r, c++, mean_ci.v_ncl[0]);
      out_at.set_entry(r, c++, mean_ci.v_ncu[0]);
      out_at.set_entry(r, c++, stdev_ci.v);
      out_at.set_entry(r, c++, v.min());      
      out_at.set_entry(r, c++, v.percentile_array(0.10));
      out_at.set_entry(r, c++, v.percentile_array(0.25));
      out_at.set_entry(r, c++, v.percentile_array(0.50));
      out_at.set_entry(r, c++, v.percentile_array(0.75));
      out_at.set_entry(r, c++, v.percentile_array(0.90));
      out_at.set_entry(r, c++, v.max());
      out_at.set_entry(r, c++, v.sum());
      if(is_bad_data(dsec)) out_at.set_entry(r, c++, dsec);
      else                  out_at.set_entry(r, c++, sec_to_hhmmss(dsec));
      if(is_bad_data(tind)) out_at.set_entry(r, c++, tind);
      else                  out_at.set_entry(r, c++, sec_to_hhmmss(tind));
      out_at.set_entry(r, c++, nint(fsp_total[r-1]));
      out_at.set_entry(r, c++, nint(fsp_best[r-1]));
      out_at.set_entry(r, c++, nint(fsp_ties[r-1]));
      out_at.set_entry(r, c++, fsp);
   }

   // Build a simple output line
   line << "JOB_LIST: " << serialize() << "\n";
   out  << line;

   // Write the table
   out << out_at << "\n" << flush;

   return;
}

////////////////////////////////////////////////////////////////////////

void TCStatJobSummary::compute_fsp(NumArray &total, NumArray &best,
                                   NumArray &ties) {
   map<ConcatString,MapData,cs_cmp>::iterator it;
   StringArray case_list;
   double v;
   double best_val;
   ConcatString best_mod;
   int i, j, k, n;

   // Initialize counts
   total.clear();
   best.clear();
   ties.clear();

   // Loop over the SummaryMap and initialize counts
   for(it=SummaryMap.begin(); it!=SummaryMap.end(); it++) {

      // Initialize counts to zero
      total.add(0);
      best.add(0);
      ties.add(0);
   } // end for it
     
   // Only compute FSP when AMODEL is in the case information
   if(!Case.has("AMODEL")) {
      mlog << Debug(4)
           << "Skipping frequency of superior performance since "
           << "the case information does not contain \"AMODEL\".\n";
      return;
   }

   // Compute FSP cases as the intersection of the cases for each entry
   for(it=SummaryMap.begin(); it!=SummaryMap.end(); it++) {

      // Initialize case list to the first entry
      if(it==SummaryMap.begin()) {
         case_list = it->second.Hdr;
      }
      // Take the intersection of the remaining entries
      else {
         case_list = intersection(case_list, it->second.Hdr);
      }
   } // end for it   

   mlog << Debug(4)
        << "Computing frequency of superior performance for "
        << Column.n_elements() << " columns and "
        << case_list.n_elements() << " cases.\n";   

   // Loop over the columns being summarized
   for(i=0; i<Column.n_elements(); i++) {

      // Check if FSP should be computed for this column
      if(strcasestr(Column[i], "-")   == NULL &&
         strcasestr(Column[i], "ERR") == NULL) {
         mlog << Debug(4)
              << "Skipping frequency of superior performance for "
              << "column \"" << Column[i] << "\" since it is not an "
              << "error or difference.\n";
         continue;
      }

      // Loop over the cases
      for(j=0; j<case_list.n_elements(); j++) {

         // Initialize top performer
         best_mod = na_str;
         best_val = bad_data_double;

         // Loop over the SummaryMap
         for(it=SummaryMap.begin(); it!=SummaryMap.end(); it++) {

            // Loop over the values for this entry
            for(k=0; k<it->second.Hdr.n_elements(); k++) {
              
               // Check if entry matches the current case
               if(strncasecmp(Column[i], it->first,
                              strlen(Column[i])) != 0 ||
                  strcmp(case_list[j], it->second.Hdr[k]) != 0)
                  continue;

               // Store the value
               v = it->second.Val[k];

               // Check for bad data
               if(is_bad_data(v)) continue;

               // Check for an meaningful improvment
               if(is_bad_data(best_val) ||
                  FSPThresh.check(abs(best_val) - abs(v))) {
                  best_val = v;
                  best_mod = it->second.AModel[k];
               }
               // Check for ties
               else if(abs(v) <= abs(best_val)) {
                  best_val = v;
                  best_mod = "TIE";
               }

            } // end for k
         } // end for it

         mlog << Debug(4)
              << "For case \"" << Column[i] << ":" << case_list[j]
              << "\" superior performance of " << best_val
              << " by \"" << best_mod << "\".\n";

         // Loop over the SummaryMap to update counts
         for(it=SummaryMap.begin(),n=0; it!=SummaryMap.end(); it++,n++) {
           
            // Check if entry has the current case
            if(strncasecmp(Column[i], it->first,
                           strlen(Column[i])) != 0 ||
               !it->second.Hdr.has(case_list[j]))
               continue;

            // Update total
            total.set(n, total[n]+1);

            // See if this entry is the top performer
            if(strcmp(best_mod, it->second.AModel[0]) == 0)
               best.set(n, best[n]+1);
            
            // See if there was a tie
            if(strcmp(best_mod, "TIE") == 0)
               ties.set(n, ties[n]+1);

         } // end for it
              
      } // end for j
   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////
//
// Code for misc functions
//
////////////////////////////////////////////////////////////////////////

TCStatJobType string_to_tcstatjobtype(const char *s) {
   TCStatJobType t;

        if(strcasecmp(s, TCStatJobType_FilterStr) == 0)  t = TCStatJobType_Filter;
   else if(strcasecmp(s, TCStatJobType_SummaryStr) == 0) t = TCStatJobType_Summary;
   else                                              t = NoTCStatJobType;

   return(t);
}

////////////////////////////////////////////////////////////////////////

ConcatString tcstatjobtype_to_string(const TCStatJobType t) {
   const char *s = (const char *) 0;

   switch(t) {
      case TCStatJobType_Filter:  s = TCStatJobType_FilterStr;  break;
      case TCStatJobType_Summary: s = TCStatJobType_SummaryStr; break;
      default:                    s = na_str;                   break;
   }

   return(ConcatString(s));
}

////////////////////////////////////////////////////////////////////////

void clear_map_data(MapData &d) {
   d.Val.clear();
   d.Hdr.clear();
   d.AModel.clear();
   d.Init.clear();
   d.Lead.clear();
   d.Valid.clear();

   return;
}

////////////////////////////////////////////////////////////////////////
//
// A time-series exists if either the init, valid, or lead time are
// fixed and the other two entries vary by a fixed amount.
//
////////////////////////////////////////////////////////////////////////

bool is_time_series(const TimeArray &init, const NumArray &lead,
                    const TimeArray &valid, int &dsec) {
   int i, dinit, dlead, dvalid;

   // Initialize
   dsec = bad_data_int;

   // The arrays should all be of the same length > 1
   if(init.n_elements() != lead.n_elements() ||
      init.n_elements() != valid.n_elements() ||
      init.n_elements() < 2) {
      mlog << Debug(4)
           << "Skipping time-series computations since the array "
           << "length differ.\n";
      return(false);
   }

   // Initialize time spacing
   dinit  = init[1]  - init[0];
   dlead  = lead[1]  - lead[0];
   dvalid = valid[1] - valid[0];
   
   // Loop over the entries to determine the time spacing
   for(i=0; i<init.n_elements()-1; i++) {

      // Make sure the time spacing remains fixed
      if(dinit != (init[i+1] - init[i])) {
         mlog << Debug(4)
              << "Skipping time-series computations since the "
              << "initialization time spacing changed: " << dinit
              << " != " << (init[i+1] - init[i]) << "\n";
         return(false);
      }
      else if(dlead != (lead[i+1] - lead[i])) {
         mlog << Debug(4)
              << "Skipping time-series computations since the "
              << "lead time spacing changed: " << dlead
              << " != " << (lead[i+1] - lead[i]) << "\n";
         return(false);
      }
      else if(dvalid != (valid[i+1] - valid[i])) {
         mlog << Debug(4)
              << "Skipping time-series computations since the "
              << "valid time spacing changed: " << dvalid
              << " != " << (valid[i+1] - valid[i]) << "\n";
         return(false);
      }
   }
    
   // Check for one fixed and the others varying by a non-zero amount
   if(dinit == 0 && dlead == dvalid && dlead > 0) {
      dsec = dlead;
      mlog << Debug(4)
           << "Computing time-series for initialization time \""
           << unix_to_yyyymmdd_hhmmss(init[0]) << "\" and spacing \""
           << sec_to_hhmmss(dsec) << "\".\n";
   }
   else if(dlead == 0 && dvalid == dinit && dvalid > 0) {
      dsec = dvalid;
      mlog << Debug(4)
           << "Computing time-series for valid time \""
           << unix_to_yyyymmdd_hhmmss(valid[0]) << "\" and spacing \""
           << sec_to_hhmmss(dsec) << "\".\n";
   }
   else if(dvalid == 0 && dinit == dlead && dinit > 0) {
      dsec = dinit;
      mlog << Debug(4)
           << "Computing time-series for lead time \""
           << sec_to_hhmmss(lead[0]) << "\" and spacing \""
           << sec_to_hhmmss(dsec) << "\".\n";
   }

   return(!is_bad_data(dsec));
}

////////////////////////////////////////////////////////////////////////

int compute_time_to_indep(const NumArray &val, int ds) {
   double mean, exp_runs, eff_size, tind;
   int i, n_abv, n_run_abv, n_bel, n_run_bel;
   bool cur_abv, prv_abv;

   // Compute the mean value of the time-series
   mean = val.mean();

   // Initialize
   n_abv = n_run_abv = 0;
   n_bel = n_run_bel = 0;
   prv_abv = false;
   
   // Count the number of values above and belwo the mean
   for(i=0; i<val.n_elements(); i++) {

      // Store the current state
      cur_abv = (val[i] >= mean);

      // Increment counters
      if(cur_abv) n_abv++;
      else        n_bel++;

      // Count runs above and below
      if(i == 1) {
         if(cur_abv) n_run_abv++;
         else        n_run_bel++;
      }
      else {
              if(cur_abv  && !prv_abv) n_run_abv++;
         else if(!cur_abv &&  prv_abv) n_run_bel++;
      }

      // Store previous state
      prv_abv = cur_abv;

   } // end for i
     
   // Calculate expected number of runs
   exp_runs = 1.0 + 2.0*(n_abv * n_bel)/(n_abv + n_bel);

   // Calculate effective sample size, time to independence
   eff_size = val.n_elements()*(n_run_abv + n_run_bel)/exp_runs;
   tind     = ds*val.n_elements()/eff_size;
   
   return(nint(tind));
}

////////////////////////////////////////////////////////////////////////

StringArray intersection(const StringArray &s1, const StringArray &s2) {
   StringArray s;
   int i;

   // Add elements common to both list
   for(i=0; i<s1.n_elements(); i++) {
      if(s2.has(s1[i])) s.add(s1[i]);
   }

   return(s);
}

////////////////////////////////////////////////////////////////////////

void add_string(const char *c, StringArray &sa) {
   ConcatString cs;

   // Parse input list into StringArray
   cs = c;
   sa.add(cs.split(ArgsDelim));
   
   return;
}

////////////////////////////////////////////////////////////////////////

void add_unixtime(const char *c, TimeArray &ta) {
   ConcatString cs;
   StringArray sa;
   int i;

   // Parse input list of timestrings into TimeArray
   cs = c;
   sa = cs.split(ArgsDelim);
   for(i=0; i<sa.n_elements(); i++) ta.add(timestring_to_unix(sa[i]));

   return;
}

////////////////////////////////////////////////////////////////////////

void add_seconds(const char *c, NumArray &na) {
   ConcatString cs;
   StringArray sa;
   int i;

   // Parse input list of timestrings into NumArray of seconds
   cs = c;
   sa = cs.split(ArgsDelim);
   for(i=0; i<sa.n_elements(); i++) na.add(timestring_to_sec(sa[i]));

   return;
}

////////////////////////////////////////////////////////////////////////
