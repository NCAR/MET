// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2012
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

using namespace std;

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "data_plane_util.h"
#include "interp_util.h"
#include "two_to_one.h"

#include "vx_math.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////
//
// Utility functions operating on a DataPlane
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
// Check the range of probability values and make sure it's either
// [0, 1] or [0, 100].  If it's [0, 100], rescale to [0, 1].
//
////////////////////////////////////////////////////////////////////////

void rescale_probability(DataPlane &dp) {
   double v, min_v, max_v;
   int x, y;

   //
   // Get the range of data values in the field
   //
   dp.data_range(min_v, max_v);

   //
   // Check for a valid range of probability values.
   //
   if(min_v < 0.0-loose_tol || max_v > 100.0+loose_tol) {
      mlog << Error << "\n  rescale_probability() -> "
           << "invalid range of data for a probability field: ["
           << min_v << ", " << max_v << "].\n\n";
      exit(1);
   }

   //
   // If necessary, rescale data from [0, 100] to [0, 1]
   //
   if(max_v > 1.0) {

      mlog << Debug(3)
           << "Rescaling probabilistic field from [0,100] to [0,1].\n";
     
      //
      // Divide each value by 100
      //
      for(x=0; x<dp.nx(); x++) {
         for(y=0; y<dp.ny(); y++) {

            v = dp.get(x, y);
            if(!is_bad_data(v)) dp.set(v/100.0, x, y);

         } // end for y
      } // end for x

   }

   return;

}

////////////////////////////////////////////////////////////////////////
//
// Smooth the DataPlane values using the interpolation method and
// box width specified.
//
////////////////////////////////////////////////////////////////////////

void smooth_field(const DataPlane &dp, DataPlane &smooth_dp,
                  InterpMthd mthd, int wdth, double t) {
   double v;
   int x, y, x_ll, y_ll;
   
   mlog << Debug(3)
        << "Smoothing field using the " << interpmthd_to_string(mthd)
        << "(" << wdth*wdth << ") interpolation method.\n";

   // Initialize the smoothed field to the raw field
   smooth_dp = dp;

   // Check for a width value of 1 for which no smoothing is done
   if(wdth <= 1) return;

   // Otherwise, apply smoothing to each grid point
   for(x=0; x<dp.nx(); x++) {
      for(y=0; y<dp.ny(); y++) {

         // The neighborhood width must be odd, find the lower-left
         // corner of the neighborhood
         x_ll = x - (wdth - 1)/2;
         y_ll = y - (wdth - 1)/2;

         // Compute the smoothed value based on the interpolation method
         switch(mthd) {

            case(InterpMthd_Min):     // Minimum
               v = interp_min(dp, x_ll, y_ll, wdth, t);
               break;

            case(InterpMthd_Max):     // Maximum
               v = interp_max(dp, x_ll, y_ll, wdth, t);
               break;

            case(InterpMthd_Median):  // Median
               v = interp_median(dp, x_ll, y_ll, wdth, t);
               break;

            case(InterpMthd_UW_Mean): // Unweighted Mean
               v = interp_uw_mean(dp, x_ll, y_ll, wdth, t);
               break;

            // Distance-weighted mean, least-squares fit, and bilinear
            // interpolation are omitted here since they are not
            // options for gridded data

            default:
               mlog << Error << "\n  smooth_field() -> "
                    << "unexpected interpolation method encountered: "
                    << mthd << "\n\n";
               exit(1);
               break;
         }

         // Store the smoothed value
         smooth_dp.set(v, x, y);

      } // end for y
   } // end for x

   return;
}

///////////////////////////////////////////////////////////////////////////////
//
// Convert the DataPlane field to the corresponding fractional coverage using
// the threshold critea specified.
//
///////////////////////////////////////////////////////////////////////////////

void fractional_coverage(const DataPlane &dp, DataPlane &frac_dp,
                         int wdth, SingleThresh t, double vld_t) {
   int i, j, k, n, x, y, x_ll, y_ll, y_ur, xx, yy, half_width;
   double v;
   int count_vld, count_thr;
   NumArray box_na;

   mlog << Debug(3)
        << "Computing fractional coverage field using the " << t.get_str()
        << "threshold and the " << interpmthd_to_string(InterpMthd_Nbrhd)
        << "(" << wdth*wdth << ") interpolation method.\n";
   
   // Check that width is set to 1 or greater
   if(wdth < 1) {
      mlog << Error << "\n  fractional_coverage() -> "
           << "width must be set to a value of 1 or greater.\n\n";
      exit(1);
   }

   // Initialize the fractional coverage field
   frac_dp = dp;
   frac_dp.set_constant(bad_data_double);

   // Compute the box half-width
   half_width = (wdth - 1)/2;

   // Initialize the box
   for(i=0; i<wdth*wdth; i++) box_na.add(bad_data_int);

   // Compute the fractional coverage meeting the threshold criteria
   for(x=0; x<dp.nx(); x++) {

      // Find the lower-left x-coordinate of the neighborhood
      x_ll = x - half_width;

      for(y=0; y<dp.ny(); y++) {

         // Find the lower-left y-coordinate of the neighborhood
         y_ll = y - half_width;
         y_ur = y + half_width;

         // Initialize the box for this new column
         if(y == 0) {

            // Initialize counts
            count_vld = count_thr = 0;

            for(i=0; i<wdth; i++) {

               xx = x_ll + i;

               for(j=0; j<wdth; j++) {

                  yy = y_ll + j;

                  n = DefaultTO.two_to_one(wdth, wdth, i, j);

                  // Check for being off the grid
                  if(xx < 0 || xx >= dp.nx() ||
                     yy < 0 || yy >= dp.ny()) {
                     k = bad_data_int;
                  }
                  // Check the value of v to see if it meets the threshold criteria
                  else {
                     v = dp.get(xx, yy);
                     if(is_bad_data(v))  k = bad_data_int;
                     else if(t.check(v)) k = 1;
                     else                k = 0;
                  }
                  box_na.set(n, k);

                  // Increment the counts
                  if(!is_bad_data(k)) {
                     count_vld += 1;
                     count_thr += k;
                  }

               } // end for j
            } // end for i
         } // end if

         // Otherwise, update one row of the box
         else {

            // Compute the row of the neighborhood box to be updated
            j = (y - 1) % wdth;

            for(i=0; i<wdth; i++) {

               // Index into the box
               n = DefaultTO.two_to_one(wdth, wdth, i, j);

               // Get x and y values to be checked
               xx = x_ll + i;
               yy = y_ur;

               // Decrement counts for data to be replaced
               k = box_na[n];
               if(!is_bad_data(k)) {
                  count_vld -= 1;
                  count_thr -= k;
               }

               // Check for being off the grid
               if(xx < 0 || xx >= dp.nx() ||
                  yy < 0 || yy >= dp.ny()) {
                  k = bad_data_int;
               }
               // Check the value of v to see if it meets the threshold criteria
               else {
                  v = dp.get(xx, yy);
                  if(is_bad_data(v))  k = bad_data_int;
                  else if(t.check(v)) k = 1;
                  else                k = 0;
               }
               box_na.set(n, k);

               // Increment the counts
               if(!is_bad_data(k)) {
                  count_vld += 1;
                  count_thr += k;
               }

            } // end for i
         } // end else

         // Check whether enough valid grid points were found
         if((double) count_vld/(wdth*wdth) < vld_t ||
            count_vld == 0) {
            v = bad_data_double;
         }
         // Compute the fractional coverage
         else {
            v = (double) count_thr/count_vld;
         }

         // Store the fractional coverage value
         frac_dp.set(v, x, y);

      } // end for y
   } // end for x

   return;

}

////////////////////////////////////////////////////////////////////////

void apply_mask(const DataPlane &fcst_dp, const DataPlane &obs_dp,
                const DataPlane &mask_dp,
                NumArray &f_na, NumArray &o_na) {
   int x, y;

   mlog << Debug(3)
        << "Extracting a list of valid forecast/observation matched "
        << "pairs from the input DataPlane objects.\n";
   
   // Initialize the NumArray objects
   f_na.clear();
   o_na.clear();

   if(fcst_dp.nx() != obs_dp.nx()  ||
      fcst_dp.ny() != obs_dp.ny()  ||
      fcst_dp.nx() != mask_dp.nx() ||
      fcst_dp.ny() != mask_dp.ny()) {

      mlog << Error << "\n  apply_mask() -> "
           << "data dimensions do not match\n\n";
      exit(1);
   }

   // Store the pairs in NumArray objects
   for(x=0; x<fcst_dp.nx(); x++) {
      for(y=0; y<fcst_dp.ny(); y++) {

         // Skip any grid points containing bad data values for
         // either of the raw fields or where the verification
         // masking region is turned off
         if(is_bad_data(fcst_dp.get(x, y)) ||
            is_bad_data(obs_dp.get(x, y))  ||
            !mask_dp.s_is_on(x, y)) continue;

         f_na.add(fcst_dp.get(x, y));
         o_na.add(obs_dp.get(x, y));
      } // end for y
   } // end for x

   return;
}

///////////////////////////////////////////////////////////////////////////////

void mask_bad_data(DataPlane &dp, const DataPlane &dp_mask) {
   int x, y;

   if(dp.nx() != dp_mask.nx() ||
      dp.ny() != dp_mask.ny() ) {

      cerr << "\n\nERROR: mask_bad_data() -> "
           << "grid dimensions do not match\n\n" << flush;
      exit(1);
   }

   for(x=0; x<dp.nx(); x++) {
      for(y=0; y<dp.ny(); y++) {

         if(is_bad_data(dp_mask.get(x, y)))
            dp.set(bad_data_double, x, y);
      }
   }

   return;
}

////////////////////////////////////////////////////////////////////////
