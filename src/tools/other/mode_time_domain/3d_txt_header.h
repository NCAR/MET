// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2019
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*



////////////////////////////////////////////////////////////////////////


#ifndef  __3D_TXT_HEADER_H__
#define  __3D_TXT_HEADER_H__


////////////////////////////////////////////////////////////////////////


static const char * header_3d_cols [] = {


   "VERSION", 
   "MODEL", 
   "DESC", 
   "FCST_LEAD", 
   "FCST_VALID", 
   // "FCST_T_DELTA", 
   // "FCST_ACCUM", 
   "OBS_LEAD", 
   "OBS_VALID", 
   // "OBS_T_DELTA", 
   // "OBS_ACCUM", 
   "T_DELTA", 
   "FCST_RAD", 
   "FCST_THR", 
   "OBS_RAD", 
   "OBS_THR", 
   "FCST_VAR", 
   "FCST_LEV", 
   "OBS_VAR", 
   "OBS_LEV", 

   // "OBJECT_ID", 
   // "OBJECT_CAT", 

};


static const int n_header_3d_cols = sizeof(header_3d_cols)/sizeof(*header_3d_cols);


////////////////////////////////////////////////////////////////////////


#endif   /*  __3D_TXT_HEADER_H__  */


////////////////////////////////////////////////////////////////////////

