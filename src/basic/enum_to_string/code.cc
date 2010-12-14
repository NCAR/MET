

////////////////////////////////////////////////////////////////////////


static const bool do_name_len_decl              = true;

static const bool echo_pound_define_after_endif = true;


////////////////////////////////////////////////////////////////////////


using namespace std;


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <cstdio>
#include <cmath>

#include "code.h"

#include "cal.h"


////////////////////////////////////////////////////////////////////////


static const char * sep = "////////////////////////////////////////////////////////////////////////\n";


////////////////////////////////////////////////////////////////////////


   //
   //  these have external linkage
   //


extern bool do_prefix;

extern const char * header_suffix;

extern bool do_angle_brackets;

extern bool do_array;

extern bool do_reverse;

extern bool verbose;

extern Unixtime generation_gmt;

extern const char * header_filename;

extern const char * program_name;


////////////////////////////////////////////////////////////////////////


static void forward(const EnumInfo & e, const char * lower, ostream & out);

static void reverse(const EnumInfo & e, const char * lower, ostream & out);


static void make_lowercase(const char * in, char * out);
static void make_uppercase(const char * in, char * out);

static void warning(ofstream &);

static void patch_name(char * len_name);

static void make_array(const EnumInfo &, ostream &);


/////////////////////////////////////////////////////////////////////////


void write_header(const EnumInfo & e)

{

if ( e.n_ids() == 0 )  return;

int j;
ofstream f;
char filename[256];
char lower[256];
char upper[256];
char pound_define[256];
char junk[256];
int len, scope_len, max_len;
char * len_name = (char *) 0;


//    if ( e.scope() )  sprintf(full_id, "%s::%s", e.scope(), e.id(j));


max_len = 0;



if ( e.scope() )  scope_len = strlen(e.scope()) + 2;   //  includes "::"
else              scope_len = 0;

for (j=0; j<(e.n_ids()); ++j)  {

   len = strlen(e.id(j));

   len += scope_len;

   if ( len > max_len )  max_len = len;

}


// cout << "Max len = " << max_len << "\n";

++max_len;   //  allow for trailing nul

// max_len += 4;   //  allow for preceeding "max_"

len_name = new char [max_len + 40];


if ( e.scope() )   sprintf(len_name, "max_enum_%s_%s_len", e.scope(), e.name());
else               sprintf(len_name, "max_enum_%s_len",   e.name());


patch_name(len_name);


   //
   //  construct output file name
   //

if ( do_prefix && (e.u_scope()) )  {

   sprintf(junk, "%s_%s", e.u_scope(), e.name());

   make_lowercase(junk, lower);

} else {

   make_lowercase(e.name(), lower);

}

if ( header_suffix )   sprintf(filename, "%s_to_string%s", lower, header_suffix);
else                   sprintf(filename, "%s_to_string",   lower);

if ( verbose )  cout << program_name << ":  Making header file \"" << filename << "\"\n" << flush;

make_uppercase(e.name(), upper);

sprintf(pound_define, "__%s_TO_STRING_H__", upper);

   //
   //  open output file
   //

f.open(filename);

if ( !f )  {

   cerr << "\n\n  write_header(const EnumInfo &) -> unable to open output file \"" << filename << "\"\n\n";

   exit ( 1 );

}

   //
   //  get to work
   //

f << "\n\n"
  << sep
  << "\n\n";

warning(f);

f << "\n\n"
  << sep
  << "\n\n";

f << "#ifndef  " << pound_define << "\n"
  << "#define  " << pound_define << "\n"
  << "\n\n";

f << sep
  << "\n\n"
  << "#include \"" << e.header() << "\"\n"
  << "\n\n"
  << sep
  << "\n\n";


f << "extern void " << lower << "_to_string(const ";

if ( e.scope() )   f << (e.scope()) << "::";

f << (e.name()) << ", char * out);\n";

if ( do_reverse )  {

   f << "\n\n";

   f << "extern bool ";

   f << "string_to_" << lower << "(const char *, ";

   if ( e.scope() )   f << (e.scope()) << "::";

   f << e.name() << " &";

   f << ");\n";

}

   f << "\n\n"
     << sep
     << "\n\n";

if ( do_name_len_decl )  {

   f << "   //\n"
     << "   //  minimum string length needed to hold output values from\n"
     << "   //\n"
     << "   //    the above function ... includes trailing nul\n"
     << "   //\n"
     << "\n\n"
     << "static const int " << len_name << " = " << max_len << ";\n"
     << "\n\n"
     << sep
     << "\n\n";

}

if ( do_array )  {

   make_array(e, f);

}

// f << "#endif   //  " << pound_define << "\n"
//   << "\n\n"
//   << sep
//   << "\n\n";

f << "#endif";

if ( echo_pound_define_after_endif )  f << "   /*  " << pound_define << "  */\n";

f << "\n\n"
  << sep
  << "\n\n";


   //
   //  done
   //

f.close();

if ( len_name )  { delete [] len_name;  len_name = (char *) 0; }

return;

}


/////////////////////////////////////////////////////////////////////////


void write_source(const EnumInfo & e)

{

ofstream out;
char L = '\"';
char R = '\"';
char filename[256];
char lower[256];
char junk[256];




if ( do_angle_brackets )  {

   L = '<';
   R = '>';

}


   //
   //  construct output file name
   //


if ( do_prefix && (e.u_scope()) )  {

   sprintf(junk, "%s_%s", e.u_scope(), e.name());

   make_lowercase(junk, lower);

} else {

   make_lowercase(e.name(), lower);

}

sprintf(filename, "%s_to_string.cc", lower);

if ( verbose )  cout << program_name << ":  Making source file \"" << filename << "\"\n" << flush;

   //
   //  open output file
   //

out.open(filename);

if ( !out )  {

   cerr << "\n\n  write_source(const EnumInfo &) -> unable to open output file \"" << filename << "\"\n\n";

   exit ( 1 );

}

   //
   //  get to work
   //

out << "\n\n"
    << sep
    << "\n\n";

warning(out);

out << "\n\n"
    << sep
    << "\n\n"
    << "using namespace std;\n"
    << "\n\n"
    << "#include <string.h>\n"
    << "\n"
    << "#include " << L << lower << "_to_string";

if ( header_suffix )    out << header_suffix;

out << R 
    << "\n"
    << "\n\n"
    << sep
    << "\n\n";


forward(e, lower, out);

out << "\n\n"
    << sep
    << "\n\n";

if ( do_reverse )  {

   reverse(e, lower, out);

   out << "\n\n"
       << sep
       << "\n\n";

}


   //
   //  done
   //

out.close();

return;

}


/////////////////////////////////////////////////////////////////////////


void forward(const EnumInfo & e, const char * lower, ostream & out)

{

if ( e.n_ids() == 0 )  return;

int j, k, m;
int max_len;
char full_id[256];




out << "void " << lower << "_to_string(const ";

if ( e.scope() )  out << (e.scope()) << "::";

out << e.name() << " t, char * out)\n"
    << "\n"
    << "{\n"
    << "\n";


out << "switch ( t )  {\n"
    << "\n";


max_len = e.max_id_length();

for (j=0; j<(e.n_ids()); ++j)  {

   if ( e.scope() )  sprintf(full_id, "%s::%s", e.scope(), e.id(j));
   else              sprintf(full_id, "%s", e.id(j));

   k = strlen(e.id(j));

   out << "   " << "case " << full_id << ":   ";

   for (m=k; m<max_len; ++m)  out.put(' ');

   out << "strcpy(out, \"" << full_id << "\");";

   for (m=k; m<max_len; ++m)  out.put(' ');

   out << "   " << "break;";

     //
     //  done with this id
     //

   out << "\n";

   if ( (j%5) == 4 )  out << "\n";

}

out << "\n"
    << "   default:\n"
    << "      strcpy(out, \"(bad value)\");\n"
    << "      break;\n";

out << "\n"
    << "}   //  switch\n"
    << "\n\n";

out << "return;\n"
    << "\n";

out << "}\n";


   //
   //  done
   //

return;

}


/////////////////////////////////////////////////////////////////////////


void reverse(const EnumInfo & e, const char * lower, ostream & out)

{

int j, k, m, n;
int max_len;


max_len = e.max_id_length();


out << "bool string_to_" << lower << "(const char * text, ";

if ( e.scope() )  out << (e.scope()) << "::";

out << e.name() << " & t)\n"
    << "\n"
    << "{\n\n";


   //
   //  do a bunch of if/then/else's
   //

n = e.n_ids();

for (j=0; j<n; ++j)  {

   k = strlen(e.id(j));

   if ( j == 0 )   out << "     ";
   else            out << "else ";

   out << "if ( strcmp(text, \"";

   if ( e.scope() )  out << (e.scope()) << "::";

   out << e.id(j) << "\"";

   for (m=k; m<max_len; ++m)  out.put(' ');

   out << ") == 0 )   { t = ";

   if ( e.scope() )  out << (e.scope()) << "::";

   out << e.id(j) << ";";

   for (m=k; m<max_len; ++m)  out.put(' ');

   out << "   " << "return ( true );";

   out << " }\n";

   if ( (j != (n - 1)) && ((j%5) == 4) )  out << "\n";

}

out << "else {\n"
    << "\n"
    << "   " << "return ( false );   //  this final \"else\" clause isn't strictly needed\n"
    << "\n"
    << "}\n\n";

out << "   //\n"
    << "   //  nope\n"
    << "   //\n"
    << "\n"
    << "return ( false );\n"
    << "\n"
    << "}\n";


   //
   //  done
   //

return;

}


/////////////////////////////////////////////////////////////////////////


   //
   //  Code for misc functions
   //


/////////////////////////////////////////////////////////////////////////


void make_lowercase(const char * in, char * out)

{

int j, k;

k = strlen(in);

for (j=0; j<k; ++j)  {

   out[j] = (char) tolower(in[j]);

}

out[j] = (char) 0;


return;

}


/////////////////////////////////////////////////////////////////////////


void make_uppercase(const char * in, char * out)

{

int j, k;

k = strlen(in);

for (j=0; j<k; ++j)  {

   out[j] = (char) toupper(in[j]);

}

out[j] = (char) 0;


return;

}


/////////////////////////////////////////////////////////////////////////


void warning(ofstream & f)

{

int j;
int month, day, year, hour, minute, second;
Unixtime t;
char junk[256];
const char * ampm = "am";
const char * zone = "MST";
const char * short_name = (const char *) 0;



   //
   //  make date/time string
   //

t = generation_gmt;

t -= 7*3600;

if ( is_dst(generation_gmt) )  {

   t += 3600;

   zone = "MDT";

}

t = 60*((t + 30)/60);   //  round to nearest minute

unix_to_mdyhms(t, month, day, year, hour, minute, second);

if ( hour >= 12 )  ampm = "pm";

hour = 1 + (hour + 11)%12;

sprintf(junk, "%s %d, %d   %d:%02d %s %s", 
              month_name[month], day, year, hour, minute, ampm, zone);

   //
   //  strip the leading path from header_filename
   //

j = strlen(header_filename) - 1;

while ( (j >= 0) && (header_filename[j] != '/') )  --j;

++j;

short_name = header_filename + j;


f << "   //\n"
  << "   //  Warning:\n"
  << "   //\n"
  << "   //     This file is machine generated.\n"
  << "   //\n"
  << "   //     Do not edit by hand.\n"
  << "   //\n"
  << "   //\n"
  << "   //     Created by enum_to_string from file \"" << short_name << "\"\n"
  << "   //\n"
  << "   //     on " << junk << "\n"
  << "   //\n";



return;

}


/////////////////////////////////////////////////////////////////////////


void patch_name(char * len_name)

{

int j, n;
int pos;
char * new_name = (char *) 0;
char c;

n = strlen(len_name);

new_name = new char [n];

pos = 0;

for (j=0; j<n; ++j)  {

   c = len_name[j];

   if ( c == ':' )  {

      new_name[pos++] = '_';

      ++j;   //  skip the other ':'

   } else if ( isupper(c) )  {

      new_name[pos++] = tolower(c);

   } else {

      new_name[pos++] = c;

   }

}   //  for j

new_name[pos++] = (char) 0;

strcpy(len_name, new_name);


   //
   //  done
   //

if ( new_name )  { delete [] new_name;  new_name = (char *) 0; }

return;

}


/////////////////////////////////////////////////////////////////////////


void make_array(const EnumInfo & e, ostream & f)

{

int j, jmax;
int k;
const int max_len = e.max_id_length();



f << "static const int n_enum_" << (e.lowercase_name()) << "s = " << (e.n_ids()) << ";\n"
  << "\n\n"
  << sep
  << "\n\n";

f << "static const " << (e.name()) << " enum_" << (e.lowercase_name()) << "_array [" << (e.n_ids()) << "] = {\n\n";

jmax = e.n_ids() - 1;

for (j=0; j<=jmax; ++j)  {

   f << "   " << e.id(j);

   if ( j != jmax )   f << ", ";
   else               f << "  ";

   f << "   ";

   for (k=strlen(e.id(j)); k<max_len; ++k)  f.put(' ');

   f << "//  # " << j;

   f << "\n";

   if ( (j != jmax) && (j%5 == 4) )  f << "\n";

}

f << "\n"
  << "};\n";


   //
   //  done
   //

f << "\n\n"
  << sep 
  << "\n\n";

return;

}


/////////////////////////////////////////////////////////////////////////


