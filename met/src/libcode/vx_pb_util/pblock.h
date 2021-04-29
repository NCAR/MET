// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*



////////////////////////////////////////////////////////////////////////


#ifndef  __PBLOCK_H__
#define  __PBLOCK_H__


////////////////////////////////////////////////////////////////////////


enum PadSize {

   padsize_4,
   padsize_8,

   no_padsize

};


////////////////////////////////////////////////////////////////////////


enum Action {

   block,
   unblock,

   no_action

};

////////////////////////////////////////////////////////////////////////

extern void pblock(const char *, const char *, Action);

////////////////////////////////////////////////////////////////////////


#endif   /*  __PBLOCK_H__  */


////////////////////////////////////////////////////////////////////////


