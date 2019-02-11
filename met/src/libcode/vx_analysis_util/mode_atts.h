// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2019
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*



////////////////////////////////////////////////////////////////////////


   //
   //  Warning:
   //
   //     This file is machine generated
   //
   //     Do not edit by hand
   //
   //
   //     File generated by program "do_ma_att_header"
   //
   //     on Dec 28, 2007   2:11 pm MST
   //


////////////////////////////////////////////////////////////////////////


#ifndef  __MODE_ATTRIBUTES_H__
#define  __MODE_ATTRIBUTES_H__


////////////////////////////////////////////////////////////////////////


#include <iostream>

#include "vx_util.h"
#include "vx_cal.h"

#include "int_array.h"
#include "mask_poly.h"
#include "mode_line.h"


////////////////////////////////////////////////////////////////////////


class ModeAttributes {

   private:

      void init_from_scratch();

      void assign(const ModeAttributes &);

   public:

      ModeAttributes();
     ~ModeAttributes();
      ModeAttributes(const ModeAttributes &);
      ModeAttributes & operator=(const ModeAttributes &);

      void clear();

      void dump(ostream & out, int depth = 0) const;

      int is_keeper(const ModeLine &) const;

      void parse_command_line(StringArray &);

      void augment(const ModeAttributes &);

      MaskPoly * poly;

         //
         //  toggles
         //

      int is_fcst_toggle_set;
      int is_fcst;                //  true for fcst, false for obs

      int is_single_toggle_set;
      int is_single;              //  true for single, false for pair

      int is_simple_toggle_set;
      int is_simple;              //  true for simple, false for cluster

      int is_matched_toggle_set;
      int is_matched;             //  true for matched, false for unmatched

         //
         //  string array members
         //

      StringArray model;
      StringArray desc;
      StringArray fcst_thr;
      StringArray obs_thr;
      StringArray fcst_var;
      StringArray fcst_lev;
      StringArray obs_var;
      StringArray obs_lev;

         //
         //  int array members
         //

      IntArray fcst_lead;
      IntArray fcst_valid_hour;
      IntArray fcst_init_hour;
      IntArray fcst_accum;
      IntArray obs_lead;
      IntArray obs_valid_hour;
      IntArray obs_init_hour;
      IntArray obs_accum;
      IntArray fcst_rad;
      IntArray obs_rad;

         //
         //  unixtime min/max members
         //

      int fcst_valid_min_set;
      unixtime fcst_valid_min;

      int fcst_valid_max_set;
      unixtime fcst_valid_max;

      int obs_valid_min_set;
      unixtime obs_valid_min;

      int obs_valid_max_set;
      unixtime obs_valid_max;

      int fcst_init_min_set;
      unixtime fcst_init_min;

      int fcst_init_max_set;
      unixtime fcst_init_max;

      int obs_init_min_set;
      unixtime obs_init_min;

      int obs_init_max_set;
      unixtime obs_init_max;


         //
         //  int min/max members
         //

      int area_min_set;
      int area_min;

      int area_max_set;
      int area_max;

      int area_thresh_min_set;
      int area_thresh_min;

      int area_thresh_max_set;
      int area_thresh_max;

      int intersection_area_min_set;
      int intersection_area_min;

      int intersection_area_max_set;
      int intersection_area_max;

      int union_area_min_set;
      int union_area_min;

      int union_area_max_set;
      int union_area_max;

      int symmetric_diff_min_set;
      int symmetric_diff_min;

      int symmetric_diff_max_set;
      int symmetric_diff_max;


         //
         //  double min/max members
         //

      int centroid_x_min_set;
      double centroid_x_min;

      int centroid_x_max_set;
      double centroid_x_max;

      int centroid_y_min_set;
      double centroid_y_min;

      int centroid_y_max_set;
      double centroid_y_max;

      int centroid_lat_min_set;
      double centroid_lat_min;

      int centroid_lat_max_set;
      double centroid_lat_max;

      int centroid_lon_min_set;
      double centroid_lon_min;

      int centroid_lon_max_set;
      double centroid_lon_max;

      int axis_ang_min_set;
      double axis_ang_min;

      int axis_ang_max_set;
      double axis_ang_max;

      int length_min_set;
      double length_min;

      int length_max_set;
      double length_max;

      int width_min_set;
      double width_min;

      int width_max_set;
      double width_max;

      int aspect_ratio_min_set;
      double aspect_ratio_min;

      int aspect_ratio_max_set;
      double aspect_ratio_max;

      int curvature_min_set;
      double curvature_min;

      int curvature_max_set;
      double curvature_max;

      int curvature_x_min_set;
      double curvature_x_min;

      int curvature_x_max_set;
      double curvature_x_max;

      int curvature_y_min_set;
      double curvature_y_min;

      int curvature_y_max_set;
      double curvature_y_max;

      int complexity_min_set;
      double complexity_min;

      int complexity_max_set;
      double complexity_max;

      int intensity_10_min_set;
      double intensity_10_min;

      int intensity_10_max_set;
      double intensity_10_max;

      int intensity_25_min_set;
      double intensity_25_min;

      int intensity_25_max_set;
      double intensity_25_max;

      int intensity_50_min_set;
      double intensity_50_min;

      int intensity_50_max_set;
      double intensity_50_max;

      int intensity_75_min_set;
      double intensity_75_min;

      int intensity_75_max_set;
      double intensity_75_max;

      int intensity_90_min_set;
      double intensity_90_min;

      int intensity_90_max_set;
      double intensity_90_max;

      int intensity_user_min_set;
      double intensity_user_min;

      int intensity_user_max_set;
      double intensity_user_max;

      int intensity_sum_min_set;
      double intensity_sum_min;

      int intensity_sum_max_set;
      double intensity_sum_max;

      int centroid_dist_min_set;
      double centroid_dist_min;

      int centroid_dist_max_set;
      double centroid_dist_max;

      int boundary_dist_min_set;
      double boundary_dist_min;

      int boundary_dist_max_set;
      double boundary_dist_max;

      int convex_hull_dist_min_set;
      double convex_hull_dist_min;

      int convex_hull_dist_max_set;
      double convex_hull_dist_max;

      int angle_diff_min_set;
      double angle_diff_min;

      int angle_diff_max_set;
      double angle_diff_max;

      int aspect_diff_min_set;
      double aspect_diff_min;

      int aspect_diff_max_set;
      double aspect_diff_max;

      int area_ratio_min_set;
      double area_ratio_min;

      int area_ratio_max_set;
      double area_ratio_max;

      int intersection_over_area_min_set;
      double intersection_over_area_min;

      int intersection_over_area_max_set;
      double intersection_over_area_max;

      int curvature_ratio_min_set;
      double curvature_ratio_min;

      int curvature_ratio_max_set;
      double curvature_ratio_max;

      int complexity_ratio_min_set;
      double complexity_ratio_min;

      int complexity_ratio_max_set;
      double complexity_ratio_max;

      int percentile_intensity_ratio_min_set;
      double percentile_intensity_ratio_min;

      int percentile_intensity_ratio_max_set;
      double percentile_intensity_ratio_max;

      int interest_min_set;
      double interest_min;

      int interest_max_set;
      double interest_max;



         //
         //  toggle "set" functions
         //

      void set_fcst      ();
      void set_obs       ();

      void set_single    ();
      void set_pair      ();

      void set_simple    ();
      void set_cluster   ();

      void set_matched   ();
      void set_unmatched ();


         //
         //  string array "add" functions
         //

      void add_model    (const char *);

      void add_desc     (const char *);

      void add_fcst_thr (const char *);

      void add_obs_thr  (const char *);

      void add_fcst_var (const char *);

      void add_fcst_lev (const char *);

      void add_obs_var  (const char *);

      void add_obs_lev  (const char *);



         //
         //  int array "add" functions
         //

      void add_fcst_lead       (int);

      void add_fcst_valid_hour (int);

      void add_fcst_init_hour  (int);

      void add_fcst_accum      (int);

      void add_obs_lead        (int);

      void add_obs_valid_hour  (int);

      void add_obs_init_hour   (int);

      void add_obs_accum       (int);

      void add_fcst_rad        (int);

      void add_obs_rad         (int);



         //
         //  int max/min "set" functions
         //

      void set_area_min              (int);
      void set_area_max              (int);

      void set_area_thresh_min       (int);
      void set_area_thresh_max       (int);

      void set_intersection_area_min (int);
      void set_intersection_area_max (int);

      void set_union_area_min        (int);
      void set_union_area_max        (int);

      void set_symmetric_diff_min    (int);
      void set_symmetric_diff_max    (int);



         //
         //  double max/min "set" functions
         //

      void set_centroid_x_min                 (double);
      void set_centroid_x_max                 (double);

      void set_centroid_y_min                 (double);
      void set_centroid_y_max                 (double);

      void set_centroid_lat_min               (double);
      void set_centroid_lat_max               (double);

      void set_centroid_lon_min               (double);
      void set_centroid_lon_max               (double);

      void set_axis_ang_min                   (double);
      void set_axis_ang_max                   (double);

      void set_length_min                     (double);
      void set_length_max                     (double);

      void set_width_min                      (double);
      void set_width_max                      (double);

      void set_aspect_ratio_min               (double);
      void set_aspect_ratio_max               (double);

      void set_curvature_min                  (double);
      void set_curvature_max                  (double);

      void set_curvature_x_min                (double);
      void set_curvature_x_max                (double);

      void set_curvature_y_min                (double);
      void set_curvature_y_max                (double);

      void set_complexity_min                 (double);
      void set_complexity_max                 (double);

      void set_intensity_10_min               (double);
      void set_intensity_10_max               (double);

      void set_intensity_25_min               (double);
      void set_intensity_25_max               (double);

      void set_intensity_50_min               (double);
      void set_intensity_50_max               (double);

      void set_intensity_75_min               (double);
      void set_intensity_75_max               (double);

      void set_intensity_90_min               (double);
      void set_intensity_90_max               (double);

      void set_intensity_user_min             (double);
      void set_intensity_user_max             (double);

      void set_intensity_sum_min              (double);
      void set_intensity_sum_max              (double);

      void set_centroid_dist_min              (double);
      void set_centroid_dist_max              (double);

      void set_boundary_dist_min              (double);
      void set_boundary_dist_max              (double);

      void set_convex_hull_dist_min           (double);
      void set_convex_hull_dist_max           (double);

      void set_angle_diff_min                 (double);
      void set_angle_diff_max                 (double);

      void set_aspect_diff_min                (double);
      void set_aspect_diff_max                (double);

      void set_area_ratio_min                 (double);
      void set_area_ratio_max                 (double);

      void set_intersection_over_area_min     (double);
      void set_intersection_over_area_max     (double);

      void set_curvature_ratio_min            (double);
      void set_curvature_ratio_max            (double);

      void set_complexity_ratio_min           (double);
      void set_complexity_ratio_max           (double);

      void set_percentile_intensity_ratio_min (double);
      void set_percentile_intensity_ratio_max (double);

      void set_interest_min                   (double);
      void set_interest_max                   (double);



         //
         //  unixtime max/min "set" functions
         //

      void set_fcst_valid_min (unixtime);
      void set_fcst_valid_max (unixtime);

      void set_obs_valid_min  (unixtime);
      void set_obs_valid_max  (unixtime);

      void set_fcst_init_min (unixtime);
      void set_fcst_init_max (unixtime);

      void set_obs_init_min  (unixtime);
      void set_obs_init_max  (unixtime);



         //
         //  misc "set" or "add" functions functions
         //

      void set_mask(const char * filename);

};


////////////////////////////////////////////////////////////////////////


extern void mode_atts_usage(ostream &);


////////////////////////////////////////////////////////////////////////


#endif   /*  __MODE_ATTRIBUTES_H__  */


////////////////////////////////////////////////////////////////////////


