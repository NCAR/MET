// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2019
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////
//
//   Filename:   pcp_combine.cc
//
//   Description:
//      Based on the command line options, this tool combines one or
//      more gridded data files into a single gridded data file and
//      writes the output in NetCDF format.
//
//      The tool may be run in four different modes:
//         sum, add, subtract, derive
//
//      The sum command requires the user to specify an initialization
//      time, input accumulation interval, valid time, and output
//      accumulation interval.
//
//      The add, subtract, and derive commands require the user to
//      specify a list of input files.  Each input file may be followed
//      by an accumulation interval or a config file string, describing
//      the data to be processed.  Alternatively, use the -field command
//      line option to set all config file strings to the same thing.
//
//      The subtract command requires exactly two input files while the
//      add and derive commands support one or more input files.
//
//      In all cases, the last argument is the output NetCDF file name.
//
//   Mod#   Date      Name           Description
//   ----   ----      ----           -----------
//   000    01-24-07  Halley Gotway  New
//   001    12-07-07  Halley Gotway  Change time format from
//                    YYYY-MM-DD_HH:MM:SS to YYYYMMDD_HHMMSS.
//   002    01-31-08  Halley Gotway  Add support for the -add and
//                    -subtract options.
//   003    09/23/08  Halley Gotway  Change argument sequence for the
//                    GRIB record access routines.
//   004    02/20/09  Halley Gotway  Append _HH to the variable name
//                    for non-zero accumulation times.
//   005    12/23/09  Halley Gotway  Call the library read_pds routine.
//   006    05/21/10  Halley Gotway  Enhance to search multiple
//                    -pcp_dir directory arguments.
//   007    06/25/10  Halley Gotway  Allow times to be specified in
//                    HH[MMSS] and YYYYMMDD[_HH[MMSS]] format.
//   008    06/30/10  Halley Gotway  Enhance grid equality checks.
//   009    07/27/10  Halley Gotway  Enhance to allow addition of any
//                    number of input files/accumulation intervals.
//                    Add lat/lon variables to NetCDF.
//   010    04/19/11  Halley Gotway  Bugfix for -add option.
//   011    10/20/11  Holmes         Added use of command line class to
//                    parse the command line arguments.
//   012    11/14/11  Halley Gotway  Bugfix for -add option when
//                    when handling missing data values.
//   013    12/21/11  Bullock        Ported to new repository.
//   014    03/07/12  Halley Gotway  Bugfix in get_field() function and
//                    remove unnecessary time strings.
//   015    04/12/12  Oldenburg      Support for all gridded data types.
//   016    01/23/13  Halley Gotway  Update usage statement and code
//                    cleanup.
//   017    10/17/13  Halley Gotway  Bugfix for closing file handles
//                    during pcpdir search.
//   018    04/16/14  Halley Gotway  Bugfix for the -varname option.
//   019    05/20/16  Prestopnik J   Removed -version (now in
//                    command_line.cc)
//   020    12/02/16  Halley Gotway  Change init and accumulation
//                    subtraction errors to warnings.
//   021    03/01/19  Halley Gotway  Add -derive command line option.
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <netcdf>

using namespace netCDF;

#include "vx_log.h"
#include "vx_data2d_factory.h"
#include "vx_data2d.h"
#include "vx_nc_util.h"
#include "vx_grid.h"
#include "vx_statistics.h"
#include "vx_util.h"
#include "vx_cal.h"
#include "vx_math.h"

////////////////////////////////////////////////////////////////////////

static ConcatString program_name;

static const char zero_time_str   [] = "00000000_000000";
static const char default_pcp_dir [] = ".";
static const char default_reg_exp [] = ".*";
static const char derive_options  [] =
                     "sum, min, max, range, mean, stdev, vld_count";

// Run Command enumeration
enum RunCommand { sum = 0, add = 1, sub = 2, der = 3 };

// Variables for top-level command line arguments
static RunCommand run_command = sum;
static int verbosity = 2;

// Variables common to all commands
static int          n_files;
static ConcatString field_string = "";
static ConcatString out_filename;
static StringArray  out_var_name;
static MetConfig    config;
static VarInfo *    var_info = (VarInfo *) 0;
static int          compress_level = -1;

// Variables for the sum command
static unixtime     init_time;
static int          in_accum;
static unixtime     valid_time;
static int          out_accum;
static StringArray  pcp_dir;
static ConcatString pcp_reg_exp = default_reg_exp;

// Variables for the add, subtract, and derive commands
static StringArray  file_list;
static StringArray  field_list;
static StringArray  derive_list;

// Output NetCDF file
NcFile *nc_out = (NcFile *) 0;
NcDim   lat_dim;
NcDim   lon_dim;

////////////////////////////////////////////////////////////////////////

static void process_command_line(int, char **);

static void process_sum_args(const CommandLine &);
static void process_add_sub_der_args(const CommandLine &);

static void do_sum_command();
static void do_add_command();
static void do_sub_command();
static void do_der_command();

static void sum_data_files(Grid &, DataPlane &);
static int  search_pcp_dir(const char *, const unixtime,
                           ConcatString &);

static void get_field(const char * filename, const int get_accum,
                      const unixtime get_init_ut,
                      const unixtime get_valid_ut,
                      Grid & grid, DataPlane & plane);

static void get_field(const char * filename, const char * cur_field,
                      const unixtime get_init_ut,
                      const unixtime get_valid_ut,
                      Grid & grid, DataPlane & plane);

static void open_nc(const Grid &);
static void write_nc_data(unixtime, unixtime, int, const DataPlane &,
                          const char *var_name = '\0',
                          const char *long_name_prefix = '\0');
static void close_nc();

static ConcatString parse_config_str(const char *);
static bool is_timestring(const char *);

static void usage();
static void set_sum(const StringArray &);
static void set_add(const StringArray &);
static void set_subtract(const StringArray &);
static void set_derive(const StringArray &);
static void set_logfile(const StringArray &);
static void set_verbosity(const StringArray &);
static void set_pcpdir(const StringArray &);
static void set_pcprx(const StringArray &);
static void set_field(const StringArray & a);
static void set_name(const StringArray & a);
static void set_compress(const StringArray &);

////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
   program_name = get_short_name(argv[0]);

   //
   // Set handler to be called for memory allocation error
   //
   set_new_handler(oom);

   //
   // Process the command line arguments
   //
   process_command_line(argc, argv);

   //
   // Perform the requested job command
   //
   if     (run_command == sum) do_sum_command();
   else if(run_command == add) do_add_command();
   else if(run_command == sub) do_sub_command();
   else                        do_der_command();

   return(0);
}

////////////////////////////////////////////////////////////////////////

void process_command_line(int argc, char **argv) {
   CommandLine cline;

   //
   // Check for zero arguments
   //
   if (argc == 1) usage();

   //
   // Default to running the sum command
   //
   run_command = sum;

   //
   // parse the command line into tokens
   //
   cline.set(argc, argv);

   //
   // Set the usage function
   //
   cline.set_usage(usage);

   //
   // Add the options function calls
   //
   cline.add(set_sum,        "-sum",        0);
   cline.add(set_add,        "-add",        0);
   cline.add(set_subtract,   "-subtract",   0);
   cline.add(set_derive,     "-derive",     1);
   cline.add(set_pcpdir,     "-pcpdir",     1);
   cline.add(set_pcprx,      "-pcprx",      1);
   cline.add(set_field,      "-field",      1);
   cline.add(set_name,       "-name",       1);
   cline.add(set_name,       "-varname",    1);
   cline.add(set_logfile,    "-log",        1);
   cline.add(set_verbosity,  "-v",          1);
   cline.add(set_compress,   "-compress",   1);

   //
   // Parse the command line
   //
   cline.parse();

   //
   // Set the verbosity level
   //
   mlog.set_verbosity_level(verbosity);

   //
   // Process the specific command arguments
   //
   if(run_command == sum) process_sum_args(cline);
   else                   process_add_sub_der_args(cline);

   //
   // If pcp_dir is not set, set to the default.
   //
   if(pcp_dir.n_elements() == 0) pcp_dir.add(default_pcp_dir);

   //
   // Initialize the MetConfig object
   //
   config.read(replace_path(config_const_filename));

   //
   // Check the -name option.
   //
   if(run_command != der && out_var_name.n() > 1) {
      mlog << Error << "\nprocess_command_line() -> "
           << "Unless running the \"-derive\" command, only 1 output "
           << "variable \"-name\" (" << write_css(out_var_name)
           << ") may be specified!\n\n";
      exit(1);
   }
   else if(run_command == der &&
           out_var_name.n() > 1 &&
           out_var_name.n() != derive_list.n()) {
      mlog << Error << "\nprocess_command_line() -> "
           << "When running the \"-derive\" command with the "
           << "\"-name \" option, the number of output variable names ("
           << out_var_name.n() << ") must match the number of derived "
           << "fields (" << derive_list.n() << ")!\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void process_sum_args(const CommandLine & cline) {

   //
   // Check for enough command line arguments
   //   init_time, in_accum, valid_time, out_accum, out_file
   //
   if(cline.n() != 5) {
      mlog << Error << "\nprocess_sum_args() -> "
           << "expected exactly 5 arguments but got " << cline.n()
           << "!\n\n";
      exit(1);
   }

   //
   // Init time
   //
   if(strcmp(cline[0], zero_time_str) == 0) {
      init_time = (unixtime) 0;
   }
   else {
      init_time = timestring_to_unix(cline[0]);
   }

   //
   // Input accumulation
   //
   in_accum = timestring_to_sec(cline[1]);

   //
   // Valid time
   //
   if(strcmp(cline[2], zero_time_str) == 0) {
      valid_time = (unixtime) 0;
   }
   else {
      valid_time = timestring_to_unix(cline[2]);
   }

   //
   // Output accumulation
   //
   out_accum = timestring_to_sec(cline[3]);

   //
   // Out file
   //
   out_filename = cline[4];

   //
   // Check that accumulation intervals are greater than zero
   //
   if(in_accum <= 0 || out_accum <= 0) {
      mlog << Error << "\nprocess_sum_args() -> "
           << "The input accumulation interval (" << cline[1]
           << ") and output accumulation interval (" << cline[3]
           << ") must be greater than zero.\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void process_add_sub_der_args(const CommandLine & cline) {
   int i;

   //
   // Check for enough command line arguments
   //
   if(cline.n() < 2) {
      mlog << Error << "\nprocess_add_sub_der_args() -> "
           << "expected at least 2 arguments but got " << cline.n()
           << "!\n\n";
      exit(1);
   }

   //
   // Store last entry as the output file
   //
   out_filename = cline[cline.n()-1];

   //
   // If the -field command line option was used, process remaining
   // arguments as a list of file names.
   //
   if(!field_string.empty()) {
      mlog << Debug(2)
           << "Since the \"-field\" command line option was used, "
           << "parsing the command line arguments as a list of "
           << "files.\n";

      //
      // If one input file was specified, check for an ascii file list.
      //
      if(cline.n() == 2) {
         Met2dDataFileFactory mtddf_factory;
         Met2dDataFile *mtddf = (Met2dDataFile *) 0;
         config.read_string(parse_config_str(field_string));
         GrdFileType type = parse_conf_file_type(&config);

         //
         // Attempt to read the first file as a gridded data file.
         // If the read was successful, store the file name.
         // Otherwise, process as an ascii file list.
         //
         if((mtddf = mtddf_factory.new_met_2d_data_file(cline[0],
                                                        type))) {
            file_list.add(cline[0]);
         }
         else {
            mlog << Debug(1)
                 << "Parsing input file names from ASCII file list: "
                 << cline[0] << "\n";
            file_list = parse_ascii_file_list(cline[0]);
         }

         //
         // Cleanup
         //
         if(mtddf) { delete mtddf; mtddf = (Met2dDataFile *) 0; }
      }
      //
      // Otherwise, store list of multiple input files.
      //
      else {
         for(i=0; i<(cline.n()-1); i++) {
            file_list.add(cline[i]);
         }
      }

      //
      // Store the corresponding configuration strings.
      //
      for(i=0; i<file_list.n(); i++) {
         field_list.add(field_string);
      }
   }
   //
   // If the -field command line option was not used, process remaining
   // arguments as pairs of filenames followed by configuration strings.
   //
   else {
      mlog << Debug(2)
           << "Since the \"-field\" command line option was not used, "
           << "parsing the command line arguments a list of files, "
           << "each followed by a configuration string.\n";

      //
      // Without -field, at least 3 arguments are required.
      //
      if(cline.n() < 3) {
         mlog << Error << "\nprocess_add_sub_der_args() -> "
              << "when not using the \"-field\" option, at least 3 "
              << "arguments are required!\n\n";
         exit(1);
      }

      for(i=0, n_files=0; i<(cline.n() - 1); i+=2) {
         file_list.add(cline[i]);

         //
         // Check if this is actually a file name
         //
         if(file_exists(cline[i+1])) {
            mlog << Error << "\nprocess_add_sub_der_args() -> "
                 << "file name used when config string expected ("
                 << cline[i+1]
                 << "). Did you forget the \"-field\" option?\n\n";
            exit(1);
         }
         field_list.add(cline[i+1]);
      }
   }

   //
   // Store the number of input files.
   //
   n_files = file_list.n();

   return;
}

////////////////////////////////////////////////////////////////////////

void do_sum_command() {
   DataPlane plane;
   Grid grid;
   int lead_time;
   ConcatString init_time_str;

   //
   // Compute the lead time
   //
   lead_time = valid_time - init_time;

   //
   // Build init time string
   //
   if(init_time != 0) {
      init_time_str = unix_to_yyyymmdd_hhmmss(init_time);
   }
   else {
      init_time_str = zero_time_str;
   }

   mlog << Debug(2)
        << "Performing sum command: "
        << "Init/In_Accum/Valid/Out_Accum Times = "
        << init_time_str  << "/" << sec_to_hhmmss(in_accum)  << "/"
        << unix_to_yyyymmdd_hhmmss(valid_time) << "/"
        << sec_to_hhmmss(out_accum) << "\n";

   //
   // Check that the output accumulation time is not greater than
   // the lead time, except when init_time = 0 for observations.
   //
   if(out_accum > lead_time && init_time != (unixtime) 0) {
      mlog << Error << "\ndo_sum_command() -> "
           << "The output accumulation time ("
           << sec_to_hhmmss(out_accum)
           << ") cannot be greater than the lead time ("
           << sec_to_hhmmss(lead_time) << ").\n\n";
      exit(1);
   }

   //
   // Check that the output accumulation time is divisible by the input
   // accumulation time.
   //
   if(out_accum%in_accum != 0) {
      mlog << Error << "\ndo_sum_command() -> "
           << "The output accumulation time ("
           << sec_to_hhmmss(out_accum)
           << ") must be divisible by the input accumulation "
           << "time (" << sec_to_hhmmss(in_accum) << ").\n\n";
      exit(1);
   }

   //
   // Check that the lead time is divisible by the the input
   // accumulation time except when init_time = 0 for observations.
   //
   if(lead_time%in_accum != 0 && init_time != (unixtime) 0) {
      mlog << Error << "\ndo_sum_command() -> "
           << "The lead time (" << sec_to_hhmmss(lead_time)
           << ") must be divisible by the input accumulation time ("
           << sec_to_hhmmss(in_accum) << ").\n\n";
      exit(1);
   }

   //
   // Find and sum up the matching precipitation files
   //
   sum_data_files(grid, plane);

   //
   // Write output.
   //
   open_nc(grid);
   write_nc_data(init_time, valid_time, out_accum, plane);
   close_nc();

   return;
}

////////////////////////////////////////////////////////////////////////

void sum_data_files(Grid & grid, DataPlane & plane) {
   int i, j, x, y;
   DataPlane part;
   double v_sum, v_part;
   Grid gr;
   unixtime     * pcp_times = (unixtime *) 0;
   int          * pcp_recs  = (int *) 0;
   ConcatString * pcp_files = (ConcatString *) 0;

   //
   // Compute the number of forecast precipitation files to be found,
   // and allocate memory to store their names and times
   //
   n_files   = out_accum/in_accum;
   pcp_times = new unixtime [n_files];
   pcp_recs  = new int [n_files];
   pcp_files = new ConcatString [n_files];

   mlog << Debug(2)
        << "Searching for " << n_files << " files "
        << "with accumulation times of " << sec_to_hhmmss(in_accum)
        << " to sum to a total accumulation time of "
        << sec_to_hhmmss(out_accum) << ".\n";

   //
   // Compute the valid times for the precipitation files
   // to be found.
   //
   for(i=0; i<n_files; i++) {
      pcp_times[i] = valid_time - i*in_accum;
   }

   //
   // Search for each file time.
   //
   for(i=0; i<n_files; i++) {

      //
      // Search in each directory for the current file time.
      //
      for(j=0; j<pcp_dir.n_elements(); j++) {

         pcp_recs[i] = search_pcp_dir(pcp_dir[j], pcp_times[i],
                                      pcp_files[i]);

         if(pcp_recs[i] != -1) {

            mlog << Debug(1)
                 << "[" << (i+1) << "] File " << pcp_files[i]
                 << " matches valid time of "
                 << unix_to_yyyymmdd_hhmmss(pcp_times[i]) << "\n";

            break;

         } // if

      } // end for j

      //
      // Check for no matching file found
      //
      if(pcp_recs[i] == -1) {

         mlog << Error << "\nsum_data_files() -> "
              << "Cannot find a file with a valid time of "
              << unix_to_yyyymmdd_hhmmss(pcp_times[i])
              << " and accumulation time of "
              << sec_to_hhmmss(in_accum) << " matching the regular "
              << "expression \"" << pcp_reg_exp << "\"\n\n";
         exit(1);
      }

   } // end for i

   /////////////////////////////

   //
   // Open each of the files found and parse the data.
   //
   for(i=0; i<n_files; i++) {

      mlog << Debug(1)
           << "[" << (i+1) << "] Reading input file: " << pcp_files[i]
           << "\n";

      //
      // Read data for the file.
      //
      get_field(pcp_files[i], in_accum, init_time, pcp_times[i], gr,
                part);

      //
      // For the first file processed store the grid, allocate memory
      // to store the precipitation sums, and initialize the sums
      //
      if(i == 0) {

         grid = gr;

         plane = part;

      }
      else {

         //
         // Check to make sure the grid stays the same
         //
         if(!(grid == gr)) {
            mlog << Error << "\nsum_data_files() -> "
                 << "the input fields must be on the same grid.\n"
                 << grid.serialize() << "\n" << gr.serialize()
                 << "\n\n";
            exit(1);
         }

         //
         // Increment the precipitation sums keeping track of the bad
         // data values
         //
         for(x=0; x<grid.nx(); x++) {

            for(y=0; y<grid.ny(); y++) {

               v_sum = plane(x, y);

               if(is_bad_data(v_sum)) continue;

               v_part = part(x, y);

               if(is_bad_data(v_part) ) {

                  plane.set(bad_data_float, x, y);

                  continue;

               }

               plane.set(v_sum + v_part, x, y);

            }   //  for y

         }   //  for x

      } // end else

   } // end for i

   //
   // Deallocate any memory that was allocated above
   //
   if(pcp_files) { delete [] pcp_files; pcp_files = (ConcatString *) 0; }
   if(pcp_times) { delete [] pcp_times; pcp_times = (unixtime *) 0; }
   if(pcp_recs ) { delete [] pcp_recs;  pcp_recs  = (int *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

int search_pcp_dir(const char *cur_dir, const unixtime cur_ut,
                   ConcatString & cur_file) {
   int i_rec;
   struct dirent *dirp = (struct dirent *) 0;
   DIR *dp = (DIR *) 0;

   //
   // Find the files matching the specified regular expression with
   // the correct valid and accumulation times.
   //
   if((dp = met_opendir(cur_dir)) == NULL ) {
      mlog << Error << "\nsearch_pcp_dir() -> "
           << "Cannot open precipitation directory "
           << cur_dir << "\n\n";
      exit(1);
   }

   //
   // Initialize the record index to not found.
   //
   i_rec = -1;

   //
   // Process each file contained in the directory.
   //
   while((dirp = readdir(dp)) != NULL) {

      //
      // Ignore any hidden files.
      //
      if(dirp->d_name[0] == '.') continue;

      //
      // Check the file name for a matching regular expression.
      //
      if(check_reg_exp(pcp_reg_exp, dirp->d_name) == true) {

         //
         // Check the current file for matching initialization,
         // valid, lead, and accumulation times.
         //
         cur_file << cs_erase << cur_dir << '/' << dirp->d_name;

         Met2dDataFileFactory factory;
         Met2dDataFile * mtddf = (Met2dDataFile *) 0;
         VarInfoFactory var_fac;
         VarInfo * cur_var = (VarInfo *) 0;

         //
         // Create a data file object.
         //
         mtddf = factory.new_met_2d_data_file(cur_file);
         if(!mtddf) {
            mlog << Warning << "search_pcp_dir() -> "
                 << "can't open data file \"" << cur_file << "\"\n";
            continue;
         }

         //
         // Create a VarInfo object from the data file.
         //
         cur_var = var_fac.new_var_info(mtddf->file_type());
         if(!cur_var) {
            mlog << Warning << "search_pcp_dir() -> "
                 << "unable to determine filetype of \"" << cur_file
                 << "\"\n";
            continue;
         }

         //
         // Initialize the VarInfo object with a field dictionary and
         // the requested timing information.
         //
         ConcatString accum_dict = field_string;
         if(field_string.empty()) {
            accum_dict.format("name=\"APCP\";level=\"A%s\";",
                              sec_to_hhmmss(in_accum).text());
         }
         config.read_string(accum_dict.text());
         cur_var->set_dict(config);

         cur_var->set_valid(cur_ut);
         cur_var->set_init(init_time);
         cur_var->set_lead(init_time ?
                           cur_ut - init_time : bad_data_int);

         //
         // Look for a VarInfo record match in the data file.
         //
         i_rec = mtddf->index(*cur_var);

         //
         // Cleanup.
         //
         if(mtddf)   { delete mtddf;   mtddf   = (Met2dDataFile *) 0; }
         if(cur_var) { delete cur_var; cur_var = (VarInfo *)       0; }

         //
         // Check for a valid match.
         //
         if(-1 != i_rec) break;

      } // end if

   } // end while

   if(closedir(dp) < 0) {
      mlog << Error << "\nsearch_pcp_dir() -> "
           << "Cannot close precipitation directory "
           << cur_dir << "\n\n";
      exit(1);
   }

   return(i_rec);
}

////////////////////////////////////////////////////////////////////////

void do_add_command() {
   Grid grid1, grid2;
   DataPlane total, part;
   double total_value, part_value;
   unixtime nc_init_time, nc_valid_time;
   int i, x, y, nc_accum;

   mlog << Debug(2)
        << "Performing addition command for " << n_files << " files.\n";

   //
   // Read the first field
   //
   get_field(file_list[0], field_list[0], 0, 0, grid1, total);

   // Initialize output times
   nc_init_time  = total.init();
   nc_valid_time = total.valid();
   nc_accum      = total.accum();

   //
   // Loop through the rest of the input files.
   // i starts at one here, not zero.
   //
   for(i=1; i<n_files; i++) {

      //
      // Read current field
      //
      get_field(file_list[i], field_list[i], 0, 0, grid2, part);

      //
      // Check for the same grid dimensions
      //
      if(grid1 != grid2 ) {
         mlog << Error << "\ndo_add_command() -> "
              << "the input fields must be on the same grid.\n"
              << grid1.serialize() << "\n" << grid2.serialize()
              << "\n\n";
         exit(1);
      }

      // Output init time
      if(nc_init_time != part.init()) nc_init_time = (unixtime) 0;

      // Output valid time
      if(nc_valid_time < part.valid()) nc_valid_time = part.valid();

      // Output accumulation time
      nc_accum += part.accum();

      //
      // Increment sums for each grid point
      //
      for(x=0; x<(part.nx()); x++) {
         for(y=0; y<(part.ny()); y++) {

            total_value = total(x, y);

            if(is_bad_data(total_value)) continue;

            part_value = part(x, y);

            if(is_bad_data(part_value) ) {

               total.set(bad_data_float, x, y);

               continue;

            }

            total.set(total_value + part_value, x, y);

         } // end for y
      } // end for x
   } // end for i

   //
   // Write output.
   //
   open_nc(grid1);
   write_nc_data(nc_init_time, nc_valid_time, nc_accum, total);
   close_nc();

   return;

}

////////////////////////////////////////////////////////////////////////

void do_sub_command() {
   DataPlane plus, minus, difference;
   Grid grid1, grid2;
   unixtime nc_init_time, nc_valid_time;
   int x, y, nc_accum;
   double v_plus, v_minus;

   //
   // Check for exactly two input files
   //
   if(n_files != 2) {
      mlog << Error << "\ndo_sub_command() -> "
           << "you must specify exactly two input files for "
           << "subtraction.\n\n";
      exit(1);
   }

   //
   // Read the two specified data files
   //
   get_field(file_list[0], field_list[0], 0, 0, grid1, plus);
   get_field(file_list[1], field_list[1], 0, 0, grid2, minus);

   //
   // Check for the same grid dimensions
   //
   if(grid1 != grid2 ) {
      mlog << Error << "\ndo_sub_command() -> "
           << "the input fields must be on the same grid.\n"
           << grid1.serialize() << "\n" << grid2.serialize() << "\n\n";
      exit(1);
   }

   //
   // Compute output accumulation, initialization, and valid times
   // for the subtract command.
   //
   mlog << Debug(2) << "Performing subtraction command.\n";

   //
   // Output valid time
   //
   nc_valid_time = plus.valid();

   //
   // Output initialization time
   // Warning if init_time1 != init_time2.
   //
   if(plus.init() != minus.init()) {
      mlog << Warning << "\ndo_sub_command() -> "
           << "the initialization times do not match ("
           << unix_to_yyyymmdd_hhmmss(plus.init()) <<  " != "
           << unix_to_yyyymmdd_hhmmss(minus.init())
           << ") for subtraction.  Using the first value.\n\n";
   }
   nc_init_time = plus.init();

   //
   // Output accumulation time
   // Warning if accum1 < accum2.
   //
   if(plus.accum() < minus.accum()) {
      mlog << Warning << "\ndo_sub_command() -> "
           << "the first accumulation interval is less than the "
           << "second (" << sec_to_hhmmss(plus.accum()) << " < "
           << sec_to_hhmmss(minus.accum()) << ") for subtraction.\n\n";
   }
   nc_accum = plus.accum() - minus.accum();

   //
   // Allocate space to store the differences
   //
   difference.set_size(plus.nx(), plus.ny());

   //
   // Perform the specified command for each grid point
   //
   difference.set_constant(bad_data_float);

   for(x=0; x<(plus.nx()); x++) {
      for(y=0; y<(plus.ny()); y++) {

         v_plus = plus(x, y);

         if(is_bad_data(v_plus)) continue;

         v_minus = minus(x, y);

         if(is_bad_data(v_minus)) continue;

         difference.set(v_plus - v_minus, x, y);

      } // end for y
   } // end for x

   //
   // Write output.
   //
   open_nc(grid1);
   write_nc_data(nc_init_time, nc_valid_time, nc_accum, difference);
   close_nc();

   return;

}

////////////////////////////////////////////////////////////////////////

void do_der_command() {
   Grid grid, cur_grid;
   DataPlane cur_dp, der_dp;
   DataPlane min_dp, max_dp, sum_dp, sum_sq_dp, vld_dp;
   unixtime nc_init_time, nc_valid_time;
   int i, j, nxy, nc_accum, nc_accum_sum;
   ConcatString out_name;
   double v;

   //
   // Check for at least one derivation option.
   //

   if(derive_list.n() == 0) {
      mlog << Error << "\ndo_der_command() -> "
           << "at least one derivation option must be specified!\n\n";
      exit(1);
   }

   mlog << Debug(2)
        << "Performing derivation command (" << write_css(derive_list)
        << ") for " << n_files << " files.\n";

   //
   // Loop through the input files.
   //
   for(i=0; i<n_files; i++) {

      //
      // Read the current field.
      //
      get_field(file_list[i], field_list[i], 0, 0, cur_grid, cur_dp);

      //
      // Initialize.
      //
      if(grid.nx() == 0 || grid.ny() == 0) {

         // Initialize the grid
         grid = cur_grid;

         // Initialize the timing information
         nc_init_time  = cur_dp.init();
         nc_valid_time = cur_dp.valid();
         nc_accum      = cur_dp.accum();
         nc_accum_sum  = cur_dp.accum();

         //
         // Initialize to the first field.
         //
         min_dp = cur_dp;
         max_dp = cur_dp;
         sum_dp = cur_dp;

         //
         // Initialize good data values to 0.
         //
         sum_sq_dp = cur_dp;
         for(j=0, nxy=grid.nx()*grid.ny(); j<nxy; j++) {
            if(!is_bad_data(cur_dp.data()[j])) {
               sum_dp.buf()[j] = 0.0;
               sum_sq_dp.buf()[j] = 0.0;
            }
         }

         //
         // Initialize valid data counts to 0.
         //
         vld_dp = cur_dp;
         vld_dp.set_constant(0.0);
      }
      //
      // Check for grid mismatch
      //
      else if(grid != cur_grid) {
         mlog << Error << "\ndo_der_command() -> "
              << "the input fields must be on the same grid.\n"
              << grid.serialize() << "\n" << cur_grid.serialize()
              << "\n\n";
         exit(1);
      }
      //
      // Update timing information
      //
      else {

         // Output init time:
         //    If it changes, reset to 0.
         if(nc_init_time != cur_dp.init()) {
            nc_init_time = (unixtime) 0;
         }

         // Output valid time:
         //    If it changes, keep track of the maximum.
         if(nc_valid_time < cur_dp.valid()) {
            nc_valid_time = cur_dp.valid();
         }

         // Output accumulation time:
         //    If it changes, reset to 0.
         if(nc_accum != cur_dp.accum()) {
            nc_accum = 0;
         }
         nc_accum_sum += cur_dp.accum();
      }

      //
      // Update sums and counts.
      //
      for(j=0, nxy=grid.nx()*grid.ny(); j<nxy; j++) {

         // Get current data value
         v = cur_dp.data()[j];

         // Update valid counts.
         if(!is_bad_data(v)) vld_dp.buf()[j] += 1;
         else                continue;

         // Update min/max fields which may contain bad data.
         if(!is_bad_data(min_dp.buf()[j]) && v < min_dp.buf()[j]) {
            min_dp.buf()[j] = v;
         }
         if(!is_bad_data(max_dp.buf()[j]) && v > max_dp.buf()[j]) {
            max_dp.buf()[j] = v;
         }

         // Update sum and sum of squares which does not have bad data.
         sum_dp.buf()[j]    += v;
         sum_sq_dp.buf()[j] += v*v;
      }
   }

   //
   // Compute derived fields and write the output.
   //
   open_nc(grid);

   //
   // Loop through the derived fields.
   //
   for(i=0; i<derive_list.n(); i++) {

      //
      // Build the output variable name.
      //
      if(out_var_name.n() == derive_list.n()) {
         out_name = out_var_name[i];
      }
      else {
         out_name = out_var_name[0];
         out_name << "_" << derive_list[i];
      }

      //
      // Write the current derived field.
      //
      if(strcasecmp(derive_list[i], "sum") == 0) {
         write_nc_data(nc_init_time, nc_valid_time, nc_accum_sum,
                       sum_dp, out_name, "Sum of ");
      }
      else if(strcasecmp(derive_list[i], "min") == 0) {
         write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                       min_dp, out_name, "Minimum Value of ");
      }
      else if(strcasecmp(derive_list[i], "max") == 0) {
      write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                    max_dp, out_name, "Maximum Value of ");
      }
      else if(strcasecmp(derive_list[i], "range") == 0) {
         der_dp = max_dp;
         for(i=0, nxy=grid.nx()*grid.ny(); i<nxy; i++) {
            if(is_bad_data(max_dp.data()[i]) ||
               is_bad_data(min_dp.data()[i])) {
               der_dp.buf()[i] = bad_data_double;
            }
            else {
               der_dp.buf()[i] = max_dp.data()[i] - min_dp.data()[i];
            }
         }
         write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                       der_dp, out_name, "Range of ");
      }
      else if(strcasecmp(derive_list[i], "mean") == 0) {
         der_dp = sum_dp;
         for(i=0, nxy=grid.nx()*grid.ny(); i<nxy; i++) {
            if(is_bad_data(sum_dp.data()[i]) ||
               is_bad_data(vld_dp.data()[i]) ||
               is_eq(vld_dp.data()[i], 0.0)) {
               der_dp.buf()[i] = bad_data_double;
            }
            else {
               der_dp.buf()[i] = sum_dp.data()[i]/vld_dp.data()[i];
            }
         }
         write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                       der_dp, out_name, "Mean Value of ");
      }
      else if(strcasecmp(derive_list[i], "stdev") == 0) {
         der_dp = sum_dp;
         for(i=0, nxy=grid.nx()*grid.ny(); i<nxy; i++) {
            double s  = sum_dp.data()[i];
            double sq = sum_sq_dp.data()[i];
            double n  = vld_dp.data()[i];
            if(is_bad_data(s) || is_bad_data(sq) ||
               is_bad_data(n) || n <= 1) {
               der_dp.buf()[i] = bad_data_double;
            }
            else {
               v = (sq - s*s/n)/(n-1);
               if(is_eq(v, 0.0)) v = 0.0;
               der_dp.buf()[i] = sqrt(v);
            }
         }
         write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                       der_dp, out_name, "Standard Deviation of ");
      }
      else if(strcasecmp(derive_list[i], "vld_count") == 0) {
         write_nc_data(nc_init_time, nc_valid_time, nc_accum,
                       vld_dp, out_name, "Valid Data Count of ");
      }
   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

void get_field(const char * filename, const int get_accum,
               const unixtime get_init_ut, const unixtime get_valid_ut,
               Grid & grid, DataPlane & plane) {
   get_field(filename, sec_to_hhmmss(get_accum), get_init_ut,
             get_valid_ut, grid, plane);
}

////////////////////////////////////////////////////////////////////////

void get_field(const char *filename, const char *cur_field,
               const unixtime get_init_ut, const unixtime get_valid_ut,
               Grid & grid, DataPlane & plane) {
   Met2dDataFileFactory factory;
   Met2dDataFile *mtddf = (Met2dDataFile *) 0;
   GrdFileType ftype;
   VarInfoFactory var_fac;
   VarInfo *cur_var;

   //
   // Build the field config string
   //
   ConcatString config_str = parse_config_str(cur_field);

   mlog << Debug(1)
        << "Reading data (" << config_str
        << ") from input file: " << filename << "\n";

   //
   // Parse the config string
   //
   config.read_string(config_str);

   //
   // Get the gridded file type from config string, if present
   //
   ftype = parse_conf_file_type(&config);

   //
   // Open the data file and build a VarInfo object
   //
   mtddf = factory.new_met_2d_data_file(filename, ftype);
   if(!mtddf) {
      mlog << Error << "\nget_field() -> "
           << "can't open data file \"" << filename << "\"\n\n";
      exit(1);
   }

   cur_var = var_fac.new_var_info(mtddf->file_type());
   if(!cur_var) {
      mlog << Error << "\nget_field() -> "
           << "unable to determine filetype of \"" << filename
           << "\"\n\n";
      exit(1);
   }

   //
   // Initialize the VarInfo object with a config
   //
   cur_var->set_dict(config);

   //
   // Set the VarInfo timing object
   //
   if(get_valid_ut != 0) cur_var->set_valid(get_valid_ut);
   if(get_init_ut  != 0) cur_var->set_init(get_init_ut);

   //
   // Read the record of interest into a DataPlane object
   //
   if(!mtddf->data_plane(*cur_var, plane)) {
      mlog << Error << "\nget_field() -> "
           << "can't get data plane from file \"" << filename
           << "\"\n\n";
      exit(1);
   }

   grid = mtddf->grid();

   //
   // Set the global var_info, if needed.
   //
   if(!var_info) {
      var_info = var_fac.new_var_info(mtddf->file_type());
      *var_info = *cur_var;
   }

   //
   // Set the output variable name, if needed.
   //
   if(out_var_name.n() == 0) {
      ConcatString cs;
      cs = var_info->magic_str();
      cs = str_replace_all(cs, "(", "");
      cs = str_replace_all(cs, ")", "");
      cs = str_replace_all(cs, "*", "");
      cs = str_replace_all(cs, ",", "");
      cs = str_replace_all(cs, "/", "_");
      out_var_name.add(cs);
   }

   //
   // Cleanup.
   //
   if(mtddf)   { delete mtddf;   mtddf   = (Met2dDataFile *) 0; }
   if(cur_var) { delete cur_var; cur_var = (VarInfo *)       0; }

   return;

}

////////////////////////////////////////////////////////////////////////

void open_nc(const Grid &grid) {
   ConcatString command_str;

   // List the output file
   mlog << Debug(1)
        << "Writing output file: " << out_filename << "\n";

   // Create a new NetCDF file and open it.
   nc_out = open_ncfile(out_filename, true);

   if(IS_INVALID_NC_P(nc_out)) {
      mlog << Error << "\nopen_nc() -> "
           << "trouble opening output file " << out_filename
           << "\n\n";
      delete nc_out;
      nc_out = (NcFile *) 0;
      exit(1);
   }

   // Add global attributes
   write_netcdf_global(nc_out, out_filename, program_name);

   if(run_command == sum) {
      command_str << cs_erase
                  << "Sum: " << n_files
                  << " files with accumulations of "
                  << sec_to_hhmmss(in_accum) << '.';
   } else if(run_command == add) {
      command_str << cs_erase
                  << "Addition: " << n_files << " files.";
   }
   else if(run_command == sub) {
      command_str << cs_erase
                  << "Subtraction: "
                  << file_list[0] << " minus " << file_list[1];
   }
   else { // run_command == der
      command_str << cs_erase
                  << "Derive: " << write_css(derive_list) << " of "
                  << n_files << " files.";
   }

   add_att(nc_out, "RunCommand", (const char *) command_str);

   // Add the projection information
   write_netcdf_proj(nc_out, grid);

   // Define Dimensions
   lat_dim = add_dim(nc_out, "lat", (long) grid.ny());
   lon_dim = add_dim(nc_out, "lon", (long) grid.nx());

   // Add the lat/lon variables
   write_netcdf_latlon(nc_out, &lat_dim, &lon_dim, grid);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nc_data(unixtime nc_init, unixtime nc_valid, int nc_accum,
                   const DataPlane &cur_dp, const char *out_name,
                   const char *long_name_prefix) {
   ConcatString var_str;
   ConcatString tmp_str, tmp2_str;
   NcVar nc_var;

   //
   // Choose the output variable name to be written.
   //

   //
   // Use the argument to this function, specified for -derive.
   //
   if(out_name) {
      var_str = out_name;
   }
   //
   // If the -name command line option was used or the accumulation
   // interval is zero, use the out_var_name.
   //
   else if(out_var_name.n() > 0 || nc_accum == 0) {
      var_str = out_var_name[0];
   }
   //
   // Otherwise, append the precipitation accumulation interval.
   //
   else {

      // Store up to the first underscore
      tmp_str        = out_var_name[0];
      StringArray sa = tmp_str.split("_");
      tmp_str        = sa[0];

      // For an hourly accumulation interval, append _HH
      if(nc_accum % sec_per_hour == 0) {
         var_str.set_precision(2);
         var_str << cs_erase << tmp_str << '_'
                 << HH(nc_accum/sec_per_hour);
      }

      // For any other accumulation interval, append _HHMMSS
      else {
         tmp2_str = sec_to_hhmmss(nc_accum);
         var_str << cs_erase << tmp_str << '_' << tmp2_str;
      }
   }

   int deflate_level = compress_level;
   if (deflate_level < 0) deflate_level = config.nc_compression();

   // Define Variable
   nc_var = add_var(nc_out, (const char *) var_str, ncFloat,
                    lat_dim, lon_dim, deflate_level);

   // Add variable attributes
   add_att(&nc_var, "name",  (const char *) var_str);
   if(long_name_prefix) tmp_str = long_name_prefix;
   else                 tmp_str.clear();
   tmp_str << var_info->long_name();
   add_att(&nc_var, "long_name", (const char *) tmp_str);

   // Ouput level string
   if(nc_accum != 0) {
      if(nc_accum%sec_per_hour == 0) {
         var_str << cs_erase << 'A' << (nc_accum/sec_per_hour);
      } else {
         var_str << cs_erase << 'A' << sec_to_hhmmss(nc_accum);
      }
   } else {
      var_str << cs_erase << var_info->level().name();
   }

   add_att(&nc_var, "level", (const char *) var_str);
   add_att(&nc_var, "units", (const char *) var_info->units());
   add_att(&nc_var, "_FillValue", bad_data_float);

   //
   // Add initialization, valid, and accumulation time info as
   // attributes to the nc_var
   //
   if(nc_init == (unixtime) 0) nc_init = nc_valid;

   //
   // Write out the times
   //
   write_netcdf_var_times(&nc_var, nc_init, nc_valid, nc_accum);

   //
   // Write the precip data
   //
   if(!put_nc_data_with_dims(&nc_var, cur_dp.data(),
                             cur_dp.ny(), cur_dp.nx())) {
      mlog << Error << "\nwrite_nc_data() -> "
           << "error with nc_var->put()\n\n";
      exit(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void close_nc() {

   //
   // Clean up.
   //
   if(nc_out)    { delete nc_out;   nc_out   = (NcFile *)  0; }
   if(var_info ) { delete var_info; var_info = (VarInfo *) 0; }

   return;
}

////////////////////////////////////////////////////////////////////////

ConcatString parse_config_str(const char *s) {
   ConcatString config_str;

   if(is_timestring(s)) {
      config_str.format("name=\"APCP\"; level=\"A%s\";", s);
   }
   else {
      config_str = s;
   }

   return(config_str);
}

////////////////////////////////////////////////////////////////////////

bool is_timestring(const char * text) {

   if(is_hh(text))     return(true);
   if(is_hhmmss(text)) return(true);

   return(false);
}

////////////////////////////////////////////////////////////////////////

void usage() {

   cout << "\n*** Model Evaluation Tools (MET" << met_version
        << ") ***\n\n"

        << "Usage: " << program_name << "\n"
        << "\t[-sum] sum_args |\n"
        << "\t-add input_files |\n"
        << "\t-subtract input_files |\n"
        << "\t-derive stat_list input_files\n"
        << "\tout_file\n"
        << "\t[-field string]\n"
        << "\t[-name list]\n"
        << "\t[-log file]\n"
        << "\t[-v level]\n"
        << "\t[-compress level]\n\n"

        << "\twhere\t\"-sum sum_args\" indicates that data from "
        << "multiple files should be summed up using the arguments "
        << "provided.\n"

        << "\t\t\"-add input_files\" indicates that data from one or "
        << "more input files should be added together.\n"

        << "\t\t\"-subtract input_files\" indicates that data from "
        << "exactly two input files should be subtracted.\n"

        << "\t\t\"-derive stat_list input_files\" indicates that "
        << "the comma-separated list of statistics in \"stat_list\"\n"
        << "\t\t(" << derive_options << ") should be derived from one "
        << "or more input files.\n"

        << "\t\t\"out_file\" is the name of the output NetCDF file to "
        << "be written (required).\n"

        << "\t\t\"-field string\" defines the data to be extracted "
        << "from the input files (optional).\n"

        << "\t\t\"-name list\" is a comma-separated list of output "
        << "variable names to be written to the \"out_file\" "
        << "(optional).\n"

        << "\t\t\"-log file\" write log messages to the specified file "
        << "(optional).\n"

        << "\t\t\"-v level\" overrides the default level of logging ("
        << verbosity << ") (optional).\n"

        << "\t\t\"-compress level\" overrides the compression level of "
        << "NetCDF variable ("  << config.nc_compression()
        << ") (optional).\n\n"

        << "\tSUM_ARGS:\n"
        << "\t\tinit_time\n"
        << "\t\tin_accum\n"
        << "\t\tvalid_time\n"
        << "\t\tout_accum\n"
        << "\t\t[-pcpdir path]\n"
        << "\t\t[-pcprx reg_exp]\n\n"

        << "\t\twhere\t\"init_time\" is the initialization time of the "
        << "input data files in YYYYMMDD[_HH[MMSS]] format "
        << "(required).\n"

        << "\t\t\t\"in_accum\" is the accumulation interval of the "
        << "input data files in HH[MMSS] format (required).\n"

        << "\t\t\t\"valid_time\" is the desired output valid time in "
        << "YYYYMMDD[_HH[MMSS]] format (required).\n"

        << "\t\t\t\"out_accum\" is the desired output accumulation "
        << "interval in HH[MMSS] format (required).\n"

        << "\t\t\t\"-pcpdir path\" overrides the default search "
        << "directory (" << default_pcp_dir << ") (optional).\n"

        << "\t\t\t\"-pcprx reg_exp\" overrides the default regular "
        << "expression for input file naming convention ("
        << default_reg_exp << ") (optional).\n\n"

        << "\t\tNote:\tSpecifying \"-sum\" is not required since it is "
        << "the default behavior.\n"
        << "\t\t\tNote: Set \"init_time\" to 00000000_000000 when "
        << "summing observation files.\n\n"

        << "\tINPUT_FILES:\n"
        << "\t\tfile_1 config_str_1 ... file_n config_str_n | \n"
        << "\t\tfile_1 ... file_n |\n"
        << "\t\tfile_list\n\n"

        << "\t\twhere\t\"file_i\" is the name of the i-th input "
        << "gridded data file.\n"

        << "\t\t\t\"config_str_i\" is the field to be extracted from "
        << "the i-th gridded data file.\n"

        << "\t\t\t\"file_list\" is an ASCII file containing a list of "
        << "gridded data files.\n\n"

        << "\t\tNote:\tFor \"-subtract\", exactly 2 input files must "
        << "be specified.\n"
        << "\t\tNote:\tThe \"-config\" option is required unless a "
        << "\"config_str\" is specified for each input file.\n"
        << "\t\tNote:\tThe \"config_str\" and \"-field\" strings may "
        << "be set to a timestring in HH[MMSS] format\n"
        << "\t\t\tfor accumulated precipitation or a full "
        << "configuration string.\n"

        << "\n" << flush;

   exit(1);
}

////////////////////////////////////////////////////////////////////////

void set_sum(const StringArray &) {
   run_command = sum;
}

////////////////////////////////////////////////////////////////////////

void set_add(const StringArray &) {
   run_command = add;
}

////////////////////////////////////////////////////////////////////////

void set_subtract(const StringArray &) {
   run_command = sub;
}

////////////////////////////////////////////////////////////////////////

void set_derive(const StringArray & a) {
   run_command = der;

   StringArray sa;
   sa.add_css(a[0]);

   //
   // Parse the options and set the corresponding flags
   //
   for(int i=0; i<sa.n(); i++) {
      if(!strcasestr(derive_options, sa[i])) {
         mlog << Error << "\nset_derive() -> "
           << "\"" << sa[i] << "\" is not a supported option for the "
           << "\"-derive\" command!\n\n";
         exit(1);
      }
      derive_list.add(sa[i]);
   }
}

////////////////////////////////////////////////////////////////////////

void set_logfile(const StringArray & a) {
   ConcatString filename;

   filename = a[0];

   mlog.open_log_file(filename);
}

////////////////////////////////////////////////////////////////////////

void set_verbosity(const StringArray & a) {
   verbosity = atoi(a[0]);
}

////////////////////////////////////////////////////////////////////////

void set_pcpdir(const StringArray & a) {
   pcp_dir.add(a[0]);
}

////////////////////////////////////////////////////////////////////////

void set_pcprx(const StringArray & a) {
   pcp_reg_exp = a[0];
}

////////////////////////////////////////////////////////////////////////

void set_field(const StringArray & a) {
   field_string = a[0];
}

////////////////////////////////////////////////////////////////////////

void set_name(const StringArray & a) {
   out_var_name.add_css(a[0]);
}

////////////////////////////////////////////////////////////////////////

void set_compress(const StringArray & a) {
   compress_level = atoi(a[0]);
}

////////////////////////////////////////////////////////////////////////