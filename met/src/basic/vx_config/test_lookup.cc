

   // *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
   // ** Copyright UCAR (c) 1992 - 2012
   // ** University Corporation for Atmospheric Research (UCAR)
   // ** National Center for Atmospheric Research (NCAR)
   // ** Research Applications Lab (RAL)
   // ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
   // *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*




////////////////////////////////////////////////////////////////////////


static const bool debug = false;


////////////////////////////////////////////////////////////////////////


using namespace std;

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <cmath>

#include "vx_util.h"

#include "config_file.h"


////////////////////////////////////////////////////////////////////////


static ConcatString program_name;


////////////////////////////////////////////////////////////////////////


int main(int argc, char * argv [])

{

program_name = get_short_name(argv[0]);

if ( argc == 1 )  {

   cerr << "\n\n  usage: " << program_name << " config_file_list\n\n";

   exit ( 1 );

}

MetConfig config;

   //
   //  redirect stderr to stdout
   //

if ( dup2(1, 2) < 0 )  {

   cerr << "\n\n  " << program_name << ": dup2 failed\n\n";

   exit ( 1 );

}

config.set_debug(debug);

int j;
bool status = false;

for (j=1; j<argc; ++j)  {   //  j starts at one, here

   status = config.read(argv[j]);

   if ( ! status )  {

      cerr << "\n\n  " << program_name << ": failed to parse config file \"" << argv[j] << "\"\n\n";

      exit ( 1 );

   }

}

   //
   //  test
   //

const DictionaryEntry * e = (const DictionaryEntry *) 0;
const DictionaryEntry * m = (const DictionaryEntry *) 0;
const DictionaryEntry * h = (const DictionaryEntry *) 0;
Dictionary * ed = (Dictionary *) 0;
Dictionary * md = (Dictionary *) 0;


cout << "\n\n";

cout << "Getting top ...\n\n";

e = config.lookup("top");

if ( !e )  {

   cerr << "\n\n  " << program_name << ": can't find top!\n\n";

   exit ( 1 );

}

// e->dump(cout);

ed = e->dict_value();

cout << "\n\n";

   //
   //  get mid
   //

m = ed->lookup("mid");

if ( !e )  {

   cerr << "\n\n  " << program_name << ": can't find mid!\n\n";

   exit ( 1 );

}

// m->dump(cout);

md = m->dict_value();

   //
   //  try to get hello from mid
   //

h = md->lookup("hello");

if ( !h )  {

   cerr << "\n\n  " << program_name << ": can't find entry!\n\n";

   exit ( 1 );

}

h->dump(cout);












cout << "\n\n";

   //
   //  done
   //

return ( 0 );

}

////////////////////////////////////////////////////////////////////////

