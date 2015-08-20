

////////////////////////////////////////////////////////////////////////


#ifndef  __3D_PAIR_ATT_COLUMNS_H__
#define  __3D_PAIR_ATT_COLUMNS_H__


////////////////////////////////////////////////////////////////////////


static const char * att_3d_pair_cols [] = {

   "SPACE_CENTROID_DIST", 
   "TIME_CENTROID_DELTA", 

   "AXIS_DIFF", 

   "SPEED_DELTA", 
   "DIRECTION_DIFF", 

   "VOLUME_RATIO", 

   "START_TIME_DELTA", 
   "END_TIME_DELTA", 

   "INTERSECTION_VOLUME", 

   "INTEREST", 

};


static const int n_att_3d_pair_cols = sizeof(att_3d_pair_cols)/sizeof(*att_3d_pair_cols);


////////////////////////////////////////////////////////////////////////


#endif   /*  __3D_PAIR_ATT_COLUMNS_H__  */


////////////////////////////////////////////////////////////////////////


