

////////////////////////////////////////////////////////////////////////


using namespace std;

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cmath>

#include "data2d_grib.h"
#include "data2d_grib_utils.h"
#include "grib_utils.h"
#include "data2d_utils.h"

#include "vx_math.h"
#include "vx_log.h"


////////////////////////////////////////////////////////////////////////

extern bool is_grid_relative(const GribRecord &);

extern int  get_bit_from_octet(unsigned char, int);


////////////////////////////////////////////////////////////////////////


   //
   //  Code for class MetGrib1DataFile
   //


////////////////////////////////////////////////////////////////////////


MetGrib1DataFile::MetGrib1DataFile()

{

grib1_init_from_scratch();

}


////////////////////////////////////////////////////////////////////////


MetGrib1DataFile::~MetGrib1DataFile()

{

close();

}


////////////////////////////////////////////////////////////////////////


MetGrib1DataFile::MetGrib1DataFile(const MetGrib1DataFile &)

{

mlog << Error << "\n  MetGrib1DataFile::MetGrib1DataFile(const MetGrib1DataFile &) -> "
     << "should never be called!\n\n";

exit ( 1 );

// grib1_init_from_scratch();
//
// assign(f);

}


////////////////////////////////////////////////////////////////////////


MetGrib1DataFile & MetGrib1DataFile::operator=(const MetGrib1DataFile &)

{

mlog << Error << "\n  MetGrib1DataFile::operator=(const MetGrib1DataFile &) -> "
     << "should never be called!\n\n";

exit ( 1 );

// if ( this == &f )  return ( * this );
//
// assign(f);

return ( * this );

}


////////////////////////////////////////////////////////////////////////


void MetGrib1DataFile::grib1_init_from_scratch()

{
  
GF = (GribFile *) 0;

Plane.clear();

close();

return;

}


////////////////////////////////////////////////////////////////////////


void MetGrib1DataFile::close()

{

if ( GF )  { delete GF;  GF = (GribFile *) 0; }

CurrentRecord.reset();

Plane.clear();

mtddf_clear();

return;

}


////////////////////////////////////////////////////////////////////////


bool MetGrib1DataFile::open(const char * _filename)

{

close();

GF = new GribFile;

if ( ! (GF->open(_filename)) )  {

   mlog << Error << "\n  MetGrib1DataFile::open(const char *) -> "
        << "unable to open grib1 file \"" << _filename << "\"\n\n";

   // exit ( 1 );

   close();

   return ( false );

}

Filename = _filename;

   //
   //  read first record & get grid
   //
   //   Note: we assume all the records in the grib file
   //
   //    use the same grid.
   //

GF->seek_record(0);

(*GF) >> CurrentRecord;

_Grid = new Grid;

gds_to_grid(*(CurrentRecord.gds), *(_Grid));

get_data_plane(CurrentRecord, Plane);


   //
   //  done
   //

return ( true );

}


////////////////////////////////////////////////////////////////////////


void MetGrib1DataFile::dump(ostream & out, int depth) const

{

Indent prefix(depth);

out << prefix << "File = ";

if ( Filename.empty() )  out << "(nul)\n";
else                     out << '\"' << Filename << "\"\n";

if ( _Grid )  {

   out << prefix << "Grid:\n";

   _Grid->dump(out, depth + 1);

} else {

   out << prefix << "No Grid!\n";

}

   //
   //  done
   //

out.flush();

return;

}


////////////////////////////////////////////////////////////////////////


double MetGrib1DataFile::get(int x, int y) const

{

double value = Plane.get(x, y);

return ( value );

}


////////////////////////////////////////////////////////////////////////


bool MetGrib1DataFile::data_ok(int x, int y) const

{

const double value = get(x, y);

return ( value != bad_data_double );

}


////////////////////////////////////////////////////////////////////////


void MetGrib1DataFile::data_minmax(double & data_min, double & data_max) const

{

Plane.data_range(data_min, data_max);

return;

}


////////////////////////////////////////////////////////////////////////


bool MetGrib1DataFile::read_record(const int n)

{

   //
   // check range
   //
if ( (n < 0) || (n > GF->n_records()) )  {

mlog << Error << "\n  MetGrib1DataFile::read_record() -> "
     << "range check error ... n = " << n << "\n\n";

exit ( 1 );

}

   //
   // seek to the requested record
   //
GF->seek_record(n);

   //
   // read it into rec
   //
if ( ! ( (*GF) >> CurrentRecord) )  {

mlog << Error << "\n  MetGrib1DataFile::read_record() -> "
     << "trouble reading record number " << n << "\n\n";

return (false);

}

   //
   // put the current record into the plane
   //
get_data_plane(CurrentRecord, Plane);

return (true);

}


////////////////////////////////////////////////////////////////////////


int MetGrib1DataFile::read_record(const VarInfoGrib & v)

{

if ( !GF )  {

   mlog << Error << "\n  MetGrib1DataFile::read_record(const VarInfoGrib &) -> "
        << "no grib file open!\n\n";

   // exit ( 1 );

   return ( -1 );

}

int j, j_match;
int count;


count = 0;

j_match = -1;

for (j=0; j<(GF->n_records()); ++j)  {

   if ( ! read_record(j) )  {

      mlog << Error << "\n  MetGrib1DataFile::read_record(const VarInfoGrib &) -> trouble reading record!\n\n";

      // exit ( 1 );

      return ( -1 );

   }

   //
   //  if an exact match is found, return only the current record
   //
       
   if ( is_exact_match(v, CurrentRecord) )  {

      count = 1;
      
      j_match = j;

      break;

   }

   //
   //  otherwise, continue looking for range matches
   //
   
   if ( is_range_match(v, CurrentRecord) )  {

      if ( j_match < 0 )  j_match = j;

      ++count;

   }

}

if ( j_match >= 0 )  {

   if ( ! read_record(j_match) )  {

      mlog << Error << "\n  MetGrib1DataFile::read_record(const VarInfoGrib &) -> "
           << "trouble reading record!\n\n";

      // exit ( 1 );

      return ( -1 );

   }

}

   //
   //  done
   //

return ( count );

}

////////////////////////////////////////////////////////////////////////

bool MetGrib1DataFile::data_plane(VarInfo &vinfo, DataPlane &plane) {
   bool status = false;
   int n_planes = 0;
   DataPlaneArray plane_array;

   // Call data_plane_array() to retrieve all matching records
   n_planes = data_plane_array(vinfo, plane_array);

   // Check for at least one match
   if(n_planes > 0) {

      // Check if first data_plane is an exact match
      // or this isn't pressure or vertical levels
      if( ( is_eq(plane_array.lower(0), vinfo.level().lower()) &&
            is_eq(plane_array.upper(0), vinfo.level().upper()) ) ||
          ( vinfo.level().type() != LevelType_Pres &&
            vinfo.level().type() != LevelType_Vert )
        ) {
     
         // Store the first data_plane
         plane  = plane_array[0];
         status = true;

         // Check for more than one matching data_plane
         if(n_planes > 1) {
            mlog << Warning << "\n  MetGrib1DataFile::data_plane() -> "
                 << "Found " << n_planes << " matches for VarInfo \""
                 << vinfo.magic_str() << "\" in GRIB file \"" << filename()
                 << "\".  Using the first match found.\n\n";
         }
      }
   } // end if n_planes

   // Check for bad status
   if(!status) {
      mlog << Warning << "\n  MetGrib1DataFile::data_plane() -> "
           << "No exact match found for VarInfo \""
           << vinfo.magic_str() << "\" in GRIB file \""
           << filename() << "\".\n\n";
   }

   return(status);
}

////////////////////////////////////////////////////////////////////////

int MetGrib1DataFile::data_plane_array(VarInfo &vinfo,
                                       DataPlaneArray &plane_array) {
   bool status = false;
   int i, lower, upper;
   GribRecord r;
   VarInfoGrib *vinfo_grib = (VarInfoGrib *) &vinfo;
   VarInfoGrib cur_vinfo_grib;
   LevelInfo cur_level;
   DataPlane cur_plane;
   DataPlaneArray u_plane_array, v_plane_array;

   // Initialize
   plane_array.clear();

   // Loop through the records in the GRIB file looking for matches
   for(i=0; i<GF->n_records(); i++) {

      // Read the current record
      GF->seek_record(i);
      (*GF) >> r;

      // Check for a range match
      if(is_range_match(*vinfo_grib, r)) {

         mlog << Debug(3) << "MetGrib1DataFile::data_plane_array() -> "
              << "Found range match for VarInfo \""
              << vinfo.magic_str() << "\" in GRIB record "
              << i+1 << " of GRIB file \"" << filename()
              << "\".\n";

         // Get the level information for this record
         read_pds_level(r, lower, upper);

         // Initialize the current VarInfo object
         cur_vinfo_grib = *vinfo_grib;
         cur_level = cur_vinfo_grib.level();

         // Reset the level range for pressure and vertical levels
         if(cur_level.type() == LevelType_Pres ||
            cur_level.type() == LevelType_Vert) {
            cur_level.set_range(lower, upper);
            cur_vinfo_grib.set_level_info(cur_level);
         }

         // Read data for the current VarInfo object
         status = data_plane_winds(cur_vinfo_grib, cur_plane);

         // Add current record to the data plane array
         plane_array.add(cur_plane, (double) lower, (double) upper);

         if(!status) {
            mlog << Warning << "\n  MetGrib1DataFile::data_plane_array() -> "
                 << "Can't read record number " << i+1
                 << " from GRIB file \"" << filename() << "\".\n\n";
            continue;
         }

      }
   } // end for loop

   // If no matches were found, check for wind records to be derived.
   if(plane_array.n_planes() == 0) {

      // Derive wind speed and direction
      if(vinfo_grib->code() == wdir_grib_code ||
         vinfo_grib->code() == wind_grib_code) {

         mlog << Debug(3) << "MetGrib1DataFile::data_plane_array() -> "
              << "Attempt to derive winds from U and V components.\n";
        
         // Initialize the current VarInfo object
         cur_vinfo_grib = *vinfo_grib;
        
         // Retrieve U-wind, doing a rotation if necessary
         mlog << Debug(3) << "MetGrib1DataFile::data_plane_array() -> "
              << "Reading U-wind records.\n";
         cur_vinfo_grib.set_code(ugrd_grib_code);
         data_plane_array(cur_vinfo_grib, u_plane_array);

         // Retrieve V-wind, doing a rotation if necessary
         mlog << Debug(3) << "MetGrib1DataFile::data_plane_array() -> "
              << "Reading V-wind records.\n";
         cur_vinfo_grib.set_code(vgrd_grib_code);
         data_plane_array(cur_vinfo_grib, v_plane_array);

         // Derive wind speed or direction
         if(u_plane_array.n_planes() != v_plane_array.n_planes()) {
            mlog << Warning << "\n  MetGrib1DataFile::data_plane_array() -> "
                 << "when deriving winds, the number of U-wind records ("
                 << u_plane_array.n_planes() << ") does not match the "
                 << "number of V-wind record (" << v_plane_array.n_planes()
                 << " for GRIB file \"" << filename() << "\".\n\n";
            return(0);
         }

         // Loop through each of the data planes
         for(i=0; i<u_plane_array.n_planes(); i++) {

            // Check that the current level values match
            if(u_plane_array.lower(i) != v_plane_array.lower(i) ||
               u_plane_array.upper(i) != v_plane_array.upper(i)) {
                 
               mlog << Warning << "\n  MetGrib1DataFile::data_plane_array() -> "
                    << "when deriving winds for level " << i+1
                    << ", the U-wind levels ("
                    << u_plane_array.lower(i) << ", " << u_plane_array.upper(i)
                    << ") do not match the V-wind levels ("
                    << v_plane_array.lower(i) << ", " << v_plane_array.upper(i)
                    << ") in GRIB file \"" << filename() << "\".\n\n";
               plane_array.clear();
               return(0);
            }
              
            // Derive wind direction
            if(vinfo_grib->code() == wdir_grib_code) {
               derive_wdir(u_plane_array[i], v_plane_array[i], cur_plane);
            }
            // Derive wind speed
            else {
               derive_wind(u_plane_array[i], v_plane_array[i], cur_plane);
            }

            // Add the current data plane
            plane_array.add(cur_plane, u_plane_array.lower(i), u_plane_array.upper(i));

         } // end for

      } // end if wdir or wind
   } // end if n_planes == 0

   mlog << Debug(3) << "MetGrib1DataFile::data_plane_array() -> "
        << "Found " << plane_array.n_planes()
        << " GRIB records matching VarInfo \""
        << vinfo.magic_str() << "\" in GRIB file \""
        << filename() << "\".\n";
        
   return(plane_array.n_planes());
}

////////////////////////////////////////////////////////////////////////
//
// This function retrieves a single data plane as requested in the
// VarInfo object and rotates winds if necessary.
//
////////////////////////////////////////////////////////////////////////

bool MetGrib1DataFile::data_plane_winds(VarInfoGrib &vinfo_grib, DataPlane &plane) {
   bool status = false;
   GribRecord r;
   VarInfoGrib vinfo_grib_winds = vinfo_grib;
   DataPlane u2d, v2d, u2d_rot, v2d_rot;

   // Initialize the data plane
   plane.clear();

   // Attempt to read the data plane
   status = data_plane_scalar(vinfo_grib, plane, r);

   // A match was found
   if(status) {

      // Check for grid-relative vectors
      if(is_grid_relative(r) &&
         (vinfo_grib.is_u_wind() ||
          vinfo_grib.is_v_wind() ||
          vinfo_grib.is_wind_direction())) {

         mlog << Debug(3) << "MetGrib1DataFile::data_plane_winds() -> "
              << "Attempt to rotate winds from earth-relative to grid-relative.\n";
        
         // For U-wind, retrieve the corresponding V-wind, and rotate
         if(vinfo_grib.is_u_wind()) {
           
            mlog << Debug(3) << "MetGrib1DataFile::data_plane_winds() -> "
                 << "Have U-wind record, reading V-wind record.\n";
            vinfo_grib_winds.set_code(vgrd_grib_code);
            data_plane_scalar(vinfo_grib_winds, v2d, r);
            rotate_uv_grid_to_earth(plane, v2d, grid(), u2d_rot, v2d_rot);
            plane = u2d_rot;
         }
         
         // For V-wind, retrieve the corresponding U-wind, and rotate
         else if(vinfo_grib.is_v_wind()) {

            mlog << Debug(3) << "MetGrib1DataFile::data_plane_winds() -> "
                 << "Have V-wind record, reading U-wind record.\n";
            vinfo_grib_winds.set_code(ugrd_grib_code);
            data_plane_scalar(vinfo_grib_winds, u2d, r);
            rotate_uv_grid_to_earth(u2d, plane, grid(), u2d_rot, v2d_rot);
            plane = v2d_rot;
         }
         
         // For wind direction, rotate
         else if(vinfo_grib.is_wind_direction()) {
            rotate_wdir_grid_to_earth(plane, grid(), u2d);
            plane = u2d;
         }

      }
   }

   if(!status) {
      mlog << Warning << "\n  MetGrib1DataFile::data_plane_winds() -> "
           << "No exact match found for VarInfo \""
           << vinfo_grib.magic_str() << "\" in GRIB file \""
           << filename() << "\".\n\n";   
   }
   
   return(status);
}

////////////////////////////////////////////////////////////////////////
//
// This function retrieves a single data plane as reqested in the
// VarInfo object but does not attempt to rotate or derive winds.
//
////////////////////////////////////////////////////////////////////////

bool MetGrib1DataFile::data_plane_scalar(VarInfoGrib &vinfo_grib,
                                            DataPlane &plane,
                                            GribRecord &record) {
   int i;
   bool status = false;

   // Initialize the data plane
   plane.clear();

   // Loop through the records in the GRIB file looking for a match
   for(i=0; i<GF->n_records(); i++) {

      // Read the current record.
      GF->seek_record(i);
      (*GF) >> record;

      // Check for an exact match
      if(is_exact_match(vinfo_grib, record)) {

         mlog << Debug(3) << "MetGrib1DataFile::data_plane_scalar() -> "
              << "Found exact match for VarInfo \""
              << vinfo_grib.magic_str() << "\" in GRIB record "
              << i+1 << " of GRIB file \"" << filename()
              << "\".\n";

         // Read current record
         status = get_data_plane(record, plane);
         break;
      }
   } // end for loop

   if(!status) {
      mlog << Warning << "\n  MetGrib1DataFile::data_plane_scalar() -> "
           << "No exact match found for VarInfo \""
           << vinfo_grib.magic_str() << "\" in GRIB file \""
           << filename() << "\".\n\n";
   }

   return(status);
}

///////////////////////////////////////////////////////////////////////////////
//
// Check whether or not the res_flag indicates that the vectors are defined
// grid relative rather than earth relative.
//
//////////////////////////////////////////////////////////////////////////////

bool is_grid_relative(const GribRecord &r) {
   unsigned char res_flag;

   // LatLon
   if(r.gds->type == 0) {
      res_flag = r.gds->grid_type.latlon_grid.res_flag;
   }
   // Mercator
   else if(r.gds->type == 1) {
      res_flag = r.gds->grid_type.mercator.res_flag;
   }
   // LambertConf
   else if(r.gds->type == 3) {
      res_flag = r.gds->grid_type.lambert_conf.res_flag;
   }
   // Stereographic
   else if(r.gds->type == 5) {
      res_flag = r.gds->grid_type.stereographic.res_flag;
   }
   else {
      cerr << "\n\nERROR: is_grid_relative() -> "
           << "Unsupported grid type value: " << r.gds->type
           << "\n\n" << flush;
      exit(1);
   }

   //
   // Return whether the 5th bit of the res_flag (Octet 17) is on, which
   // indicates that U and V are defined relative to the grid
   //
   return(get_bit_from_octet(res_flag, 5) == 1);
}

//////////////////////////////////////////////////////////////////////////////

int get_bit_from_octet(unsigned char u, int bit) {

   //
   // Bit numbers start at 1, not 0, and
   // the most-significant bit is number 1
   //

   if((bit < 1) || (bit > 8)) {

      cerr << "\n\nERROR: get_bit_from_octet() -> "
           << "bad bit number\n\n" << flush;
      exit(1);
   }

   unsigned char mask = (unsigned char) (1 << (8 - bit));

   if(u & mask) return(1);

   return(0);
}

////////////////////////////////////////////////////////////////////////
