

////////////////////////////////////////////////////////////////////////


using namespace std;

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cmath>

#include "mtd_file.h"
#include "mtd_partition.h"
#include "mtd_nc_defs.h"
#include "nc_utils.h"

#include "vx_math.h"


////////////////////////////////////////////////////////////////////////


static MtdIntFile split(const MtdIntFile & mask, int & n_shapes);

static void adjust_obj_numbers(MtdIntFile & s, int delta);

static void find_overlap(const MtdIntFile & before, const MtdIntFile & after, Mtd_Partition & p);


////////////////////////////////////////////////////////////////////////


   //
   //  Code for class MtdIntFile
   //


////////////////////////////////////////////////////////////////////////


MtdIntFile::MtdIntFile()

{

int_init_from_scratch();

}


////////////////////////////////////////////////////////////////////////


MtdIntFile::~MtdIntFile()

{

clear();

}


////////////////////////////////////////////////////////////////////////


MtdIntFile::MtdIntFile(const MtdIntFile & f)

{

int_init_from_scratch();

int_assign(f);

}


////////////////////////////////////////////////////////////////////////


MtdIntFile & MtdIntFile::operator=(const MtdIntFile & f)

{

if ( this == &f )  return ( * this );

int_assign(f);

return ( * this );

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::int_init_from_scratch()

{

Data = 0;

ObjVolume = 0;

clear();

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::clear()

{

MtdFileBase::clear();

if ( Data )  { delete [] Data;  Data = 0; }

if ( ObjVolume )  { delete [] ObjVolume;  ObjVolume = 0; }

DataMin = DataMax = 0;

Radius = -1;

Threshold = -1.0;


   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::int_assign(const MtdIntFile & f)

{

clear();

int n;

base_assign(f);

DataMin = f.DataMin;
DataMax = f.DataMax;

Radius = f.Radius;

Threshold = f.Threshold;

Nobjects = f.Nobjects;

n = Nx*Ny*Nt;

if ( f.Data )  {

   Data = new int [n];

   memcpy(Data, f.Data, n*sizeof(int));

}

n = Nobjects;

if ( f.ObjVolume )  {

   ObjVolume = new int [n];

   memcpy(ObjVolume, f.ObjVolume, n*sizeof(int));

}

   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::dump(ostream & out, int depth) const

{

Indent prefix(depth);

MtdFileBase::dump(out, depth);

out << prefix << "Radius = " << Radius << '\n';

out << prefix << "Threshold = " << Threshold << '\n';


   //
   //  done
   //

out.flush();

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_size(int _nx, int _ny, int _nt)

{


if ( Data )  { delete [] Data;  Data = 0; }

int j;
const int n3 = _nx*_ny*_nt;

Data = new int [n3];

Nx = _nx;
Ny = _ny;
Nt = _nt;

int * d = Data;

for (j=0; j<n3; ++j)  *d++ = 0;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_radius(int r)

{

if ( r < 0 )  {

   cerr << "\n\n  MtdIntFile::set_radius(int) -> bad value ... " << r << "\n\n";

   exit ( 1 );

}

Radius = r;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_threshold(double t)

{

// if ( t < 0.0 )  {
// 
//    cerr << "\n\n  MtdIntFile::set_threshold(double) -> bad value ... " << t << "\n\n";
// 
//    exit ( 1 );
// 
// }

Threshold = t;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::put(const int value, int _x, int _y, int _t)

{

const int n = mtd_three_to_one(Nx, Ny, Nt, _x, _y, _t);

Data[n] = value;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_data_minmax(double _data_min, double _data_max)

{

DataMin = _data_min;

DataMax = _data_max;

return;

}


////////////////////////////////////////////////////////////////////////


bool MtdIntFile::read(const char * _filename)

{

NcFile f(_filename);

if ( ! f.is_valid() )  return ( false );

Filename = _filename;

MtdIntFile::read(f);

return ( true );

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::read(NcFile & f)

{

NcVar * var = 0;



   //
   //  read the base class stuff
   //

MtdFileBase::read(f);

   //  DataMin, DataMax

DataMin = string_att_as_double (f, min_value_att_name);
DataMax = string_att_as_double (f, max_value_att_name);

   //  Data

set_size(Nx, Ny, Nt);

var = f.get_var(data_field_name);

if ( !(var->set_cur(0, 0, 0)) )  {

   cerr << "\n\n  MtdIntFile::read() -> trouble setting corner\n\n";

   exit ( 1 );

}

const time_t t_start = time(0);   //  for timing the data read operation

if ( ! (var->get(Data, Nt, Ny, Nx)) )  {

   cerr << "\n\n  MtdIntFile::read(const char *) -> trouble getting data\n\n";

   exit ( 1 );

}

const time_t t_stop = time(0);   //  for timing the data read operation

cout << "\n\n  MtdIntFile::read(): Time to read data = " << (t_stop - t_start) << " seconds\n\n" << flush;

   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::write(NcFile & f) const

{

NcDim * nx_dim   = 0;
NcDim * ny_dim   = 0;
NcDim * nt_dim   = 0;
NcVar * data_var = 0;
const char format [] = "%d";
char junk[256];


   //
   //  write stuff from parent class
   //

MtdFileBase::write(f);

   //
   //  get the dimensions of the data field
   //

nx_dim = f.get_dim(nx_dim_name);
ny_dim = f.get_dim(ny_dim_name);
nt_dim = f.get_dim(nt_dim_name);


   //  DataMin, DataMax

sprintf(junk, format, DataMin);

f.add_att(min_value_att_name, junk);

sprintf(junk, format, DataMax);

f.add_att(max_value_att_name, junk);

   //  Radius

f.add_att(radius_att_name, Radius);

   //  Threshold

f.add_att(threshold_att_name, Threshold);

   //  Data

f.add_var(data_field_name, ncInt, nt_dim, ny_dim, nx_dim);

data_var = f.get_var(data_field_name);

if ( !(data_var->set_cur(0, 0, 0)) )  {

   cerr << "\n\n  MtdIntFile::write() -> trouble setting corner on data field\n\n";

   exit ( 1 );

}

const time_t t_start = time(0);   //  for timing the data write operation

if ( !(data_var->put(Data, Nt, Ny, Nx)) )  {

   cerr << "\n\n  MtdIntFile::write() -> trouble with put in data field\n\n";

   exit ( 1 );

}


const time_t t_stop = time(0);   //  for timing the data write operation

cout << "\n\n  MtdIntFile::write(): Time to write data = " << (t_stop - t_start) << " seconds\n\n" << flush;

   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::write(const char * _filename) const

{

NcFile f(_filename, NcFile::Replace);

if ( ! f.is_valid() )  {

   cerr << "\n\n  MtdIntFile::write(const char *) -> unable to open netcdf output file \"" << _filename << "\"\n\n";

   // exit ( 1 );

   return;

}

MtdIntFile::write(f);

   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


MtdIntFile MtdIntFile::const_t_slice(int t) const

{

if ( (t < 0) || (t >= Nt) )  {

   cerr << "\n\n  MtdIntFile MtdIntFile::const_t_slice(int) const -> range check error\n\n";

   exit ( 1 );

}

int j, n;
const int nxy = Nx*Ny;
int bytes = nxy*sizeof(int);
int fmin, fmax;
MtdIntFile f;
int * d;
int value;

f.base_assign(*this);

f.Nt = 1;

f.DeltaT = 0;

f.StartTime = StartTime + t*DeltaT;

f.Radius = Radius;

f.Data = new int [Nx*Ny];

n = mtd_three_to_one(Nx, Ny, Nt, 0, 0, t);

memcpy(f.Data, Data + n, bytes);

d = f.Data;

fmin = fmax = f.Data[0];

for (j=0; j<nxy; ++j)  {

   value = *d++;

   if ( value < fmin )  fmin = value;
   if ( value > fmax )  fmax = value;

}

f.DataMin = fmin;
f.DataMax = fmax;

   //
   //  done
   //

return ( f );

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::fatten()

   //
   //  assumes Nt = 1 (ie, a 2D slice)
   //

{

int x, y, n;
const int nxy = Nx*Ny;
int * u = new int [nxy];
int * a = 0;

a = u;

memcpy(a, Data, nxy*sizeof(int));


   //
   //  first do the top edge.
   //
   //  if (x, y) is turned on then turn on
   //
   //  (x + 1, y)
   //

y = Ny - 1;

n = mtd_three_to_one(Nx, Ny, Nt, 0, y, 0);

a = u + n;

for (x=0; x<(Nx - 1); ++x)  {

   if ( *a++ )  Data[n + 1] = 1;

   ++n;

}

   //
   //  next, do the right edge.
   //
   //  if (x, y) is turned on then turn on
   //
   //  (x, y + 1)
   //

x = Nx - 1;

n = mtd_three_to_one(Nx, Ny, Nt, x, 0, 0);

a = u + n;

for (y=0; y<(Ny - 1); ++y)  {

   if ( *a )  Data[n + Nx] = 1;

   n += Nx;

}

   //
   //  finally, do the rest.
   //
   //  if (x, y) is on, then turn on 
   //
   //     (x + 1, y), (x, y + 1), (x + 1, y + 1)
   //

for (y = 0; y<(Ny - 2); ++y)  {

   n = mtd_three_to_one(Nx, Ny, Nt, 0, y, 0);

   a = u + n;

   for (x=0; x<(Nx - 2); ++x)  {

      if ( *a )  {

         Data[n + 1]      = 1;   //  (x + 1, y)
         Data[n + Nx]     = 1;   //  (x, y + 1)
         Data[n + Nx + 1] = 1;   //  (x + 1, y + 1)

      }

      ++a;
      ++n;

   }   //  for x

}   //  for y



   //
   //  done
   //

if ( u )  { delete [] u;  u = 0; }

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::zero_border(int n)

{

if ( !Data )  {

   cerr << "\n\n  MtdIntFile::zero_border(int) -> no data field!\n\n";

   exit ( 1 );

}

if ( 2*n >= min(Nx, Ny) )  {

   cerr << "\n\n  MtdIntFile::zero_border(int) -> border size too large!\n\n";

   exit ( 1 );

}

int x, y, t;

for (t=0; t<Nt; ++t)  {

   for (x=0; x<Nx; ++x)  {

      for (y=0; y<n; ++y)  {

         put(0, x, y,          t);
         put(0, x, Ny - 1 - y, t);

      }

   }

   for (y=0; y<Ny; ++y)  {

      for (x=0; x<n; ++x)  {

         put(0, x,          y, t);
         put(0, Nx - 1 - x, y, t);

      }

   }

}


   //
   //  done
   //

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_to_zeroes()

{

if ( !Data )  {

   cerr << "\n\n  MtdIntFile::set_to_zeroes() -> no data!\n\n";

   exit ( 1 );

}

int j;
const int n3 = Nx*Ny*Nt;
int * d = Data;

for (j=0; j<n3; ++j)  *d++ = 0;

return;

}


////////////////////////////////////////////////////////////////////////


MtdIntFile MtdIntFile::split_const_t(int & n_shapes) const

{

if ( Nt != 1 )  {

   cerr << "\n\n  split_const_t(int &) -> not const-time slice!\n\n";

   exit ( 1 );

}

int k, i, ii;
int x, y, xx, yy;
int current_shape;
int value;
bool shape_assigned = false;
MtdIntFile d;
Mtd_Partition p;
const MtdIntFile & id = *this;



d.set_size(Nx, Ny, 1);

d.set_grid(*G);

d.set_to_zeroes();

n_shapes = 0;

   //
   //  shape numbers start at ONE here!!
   //

current_shape = 0;

for (y=(Ny - 2); y>=0; --y)  {

   for (x=(Nx - 2); x>=0; --x)  {

      shape_assigned = false;

      if ( !(id(x, y, 0)) )  continue;

         //
         //  check above left
         //

      xx = x - 1;  yy = y + 1;

      if ( (xx >= 0) && (yy < Ny) && id(xx, yy, 0) )  {

         i  = d( x,  y, 0);
         ii = d(xx, yy, 0);

         if ( shape_assigned )  p.merge_values(i, ii);
         else                   d.put(ii, x, y, 0);

         shape_assigned = true;

      }

         //
         //  check above
         //

      xx = x;  yy = y + 1;

      if ( (yy < Ny) && id(xx, yy, 0) )  {

         i  = d( x,  y, 0);
         ii = d(xx, yy, 0);

         if ( shape_assigned )  p.merge_values(i, ii);
         else                   d.put(ii, x, y, 0);

         shape_assigned = true;

      }

         //
         //  check upper right
         //

      xx = x + 1;  yy = y + 1;

      if ( (xx < Nx) && (yy < Ny) && id(xx, yy, 0) )  {

         i  = d( x,  y, 0);
         ii = d(xx, yy, 0);

         if ( shape_assigned )  p.merge_values(i, ii);
         else                   d.put(ii, x, y, 0);

         shape_assigned = true;

      }

         //
         //  check to the right
         //

      xx = x + 1;  yy = y;

      if ( (xx < Nx) && id(xx, yy, 0) )  {

         i  = d( x,  y, 0);
         ii = d(xx, yy, 0);

         if ( shape_assigned )  p.merge_values(i, ii);
         else                   d.put(ii, x, y, 0);

         shape_assigned = true;

      }

         //
         //  is it a new shape?
         //

      if ( !shape_assigned )  {

         ++current_shape;

         d.put(current_shape, x, y, 0);

         p.add_no_repeat(current_shape);

      }

   }   //  for x

}   //  for y




MtdIntFile q;


q.set_size(Nx, Ny, 1);

q.set_grid(id.grid());

q.set_to_zeroes();


for (x=0; x<Nx; ++x)  {

   for (y=0; y<Ny; ++y)  {

      i = d(x, y, 0);

      value = 0;

      for (k=0; k<(p.n_elements()); ++k)  {

         if ( p.has(k, i) )  { value = k + 1;  break; }

      }

      q.put(value, x, y, 0);

   }

}


n_shapes = p.n_elements();


   //
   //  done
   //

return ( q );

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::split()

{

int j, k;
int * d = 0;
MtdIntFile old;

old = ::split(*this, Nobjects);

const int Nxyt = Nx*Ny*Nt;

if ( Nobjects == 0 )  return;

old.ObjVolume = new int [Nobjects];

d = old.Data;

for (j=0; j<Nobjects; ++j)  old.ObjVolume[j] = 0;


for (j=0; j<Nxyt; ++j, ++d)  {

   k = (int) (*d);

   if ( k == 0 )  continue;

   old.ObjVolume[k - 1] += 1;

}



   //
   //  done
   //

int_assign(old);

return;

}


////////////////////////////////////////////////////////////////////////


int MtdIntFile::volume(int k) const

{

if ( !ObjVolume )  {

   cerr << "\n\n  MtdIntFile::volume(int) -> field not split!\n\n";

   exit ( 1 );

}

if ( (k < 0) || (k >= Nobjects) )  {

   cerr << "\n\n  MtdIntFile::volume(int) -> range check error!\n\n";

   exit ( 1 );

}


return ( ObjVolume[k] );

}


////////////////////////////////////////////////////////////////////////


int MtdIntFile::total_volume() const

{

if ( !ObjVolume )  {

   cerr << "\n\n  MtdIntFile::total_volume() -> field not split!\n\n";

   exit ( 1 );

}

int j;
int sum = 0;

for (j=0; j<Nobjects; ++j)  sum += ObjVolume[j];

return ( sum );

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::toss_small_objects(int min_volume)

{

int j, n_new;
int * new_to_old = (int *) 0;

new_to_old = new int[Nobjects];   //  probably too big, but that's ok

n_new = 0;

for (j=0; j<Nobjects; ++j)  {

   if ( ObjVolume[j] >= min_volume )  new_to_old[n_new++] = j;

}


sift_objects(n_new, new_to_old);


   //
   //  done
   //

delete [] new_to_old;   new_to_old = (int *) 0;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::sift_objects(const int n_new, const int * new_to_old)

{

if ( n_new == Nobjects )  return;

int j, k;
const int n3 = Nx*Ny*Nt;
int * old_to_new = (int *) 0;
int * new_volumes = (int *) 0;
// double * new_intensities = (double *) 0;
int * d = Data;


if ( n_new == 0 )  {

   cerr << "\n\n  MtdIntFile::sift_objects() -> no objects left!\n\n";

   exit ( 1 );

}


old_to_new = new int [Nobjects];


new_volumes = new int [n_new];

// new_intensities = new double [n_new];


for (j=0; j<n_new; ++j)  {

   new_volumes[j] = ObjVolume[new_to_old[j]];

   // new_intensities[j] = MaxConvIntensity[new_to_old[j]];

}

for (j=0; j<Nobjects; ++j)  old_to_new[j] = -1;

for (j=0; j<n_new; ++j)  old_to_new[new_to_old[j]] = j;

// for (j=0; j<Nobjects; ++j)  cout << j << " ... " << old_to_new[j] << '\n';

int replace_count = 0;
int zero_count = 0;

d = Data;

for (j=0; j<n3; ++j, ++d)  {

   k = (int) (*d);

   if ( k == 0 )  continue;

   k = old_to_new[k - 1];

   if ( k < 0 )    { *d = 0;      ++zero_count; }
   else            { *d = 1 + k;  ++replace_count; }

}

// cout << "replace count = " << replace_count << '\n' << flush;
// cout << "zero count    = " << zero_count    << '\n' << flush;

   //
   //  rewire
   //

Nobjects = n_new;

DataMax  = Nobjects;


delete [] ObjVolume;  ObjVolume = (int *) 0;

ObjVolume = new_volumes;


// delete [] MaxConvIntensity;  MaxConvIntensity = (double *) 0;
// 
// MaxConvIntensity = new_intensities;

   //
   //  done
   //

delete [] old_to_new;   old_to_new = (int *) 0;

return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::calc_bbox(int & x_min, int & x_max, 
                           int & y_min, int & y_max, 
                           int & t_min, int & t_max) const

{

int x, y, t, n;

x_min = y_min = t_min = 0;

x_max = Nx;
y_max = Ny;
t_max = Nt;

for (x=0; x<Nx; ++x)  {

   for (y=0; y<Ny; ++y)  {

      for (t=0; t<Nt; ++t)  {

         n = mtd_three_to_one(Nx, Ny, Nt, x, y, t);

         if ( Data[n] )  {

            if ( x < x_min )  x_min = x;
            if ( y < y_min )  y_min = y;
            if ( t < t_min )  t_min = t;

            if ( x > x_max )  x_max = x;
            if ( y > y_max )  y_max = y;
            if ( t > t_max )  t_max = t;

         }

      }   //  for t

   }   //  for y

}   //  for x



return;

}


////////////////////////////////////////////////////////////////////////


void MtdIntFile::calc_centroid(double & xbar, double & ybar, double & tbar) const

{

int x, y, t, n;
int count;

xbar = ybar = tbar = 0.0;

count = 0;

for (x=0; x<Nx; ++x)  {

   for (y=0; y<Ny; ++y)  {

      for (t=0; t<Nt; ++t)  {

         n = mtd_three_to_one(Nx, Ny, Nt, x, y, t);

         if ( Data[n] )  {

            xbar += x;
            ybar += y;
            tbar += t;

            ++count;

         }

      }   //  for t

   }   //  for y

}   //  for x

if ( count == 0 )  {

   cerr << "\n\n  MtdIntFile::calc_centroid() const -> empty object!\n\n";

   exit ( 1 );

}

xbar /= count;
ybar /= count;
tbar /= count;

return;

}



////////////////////////////////////////////////////////////////////////


void MtdIntFile::set_volumes(int n, const int * V)

{

if ( ObjVolume )  { delete [] ObjVolume;  ObjVolume = (int *) 0; }

int j;

ObjVolume = new int [n];

for (j=0; j<n; ++j)  ObjVolume[j] = V[j];

Nobjects = n;


return;

}


////////////////////////////////////////////////////////////////////////


MtdIntFile MtdIntFile::select(int n) const

{

if ( (n < 1) || (n > Nobjects) )  {

   cerr << "\n\n  MtdIntFile::select(int) -> range check error on n ... "
        << "NObjects = " << Nobjects << " ... "
        << "n = " << n << "\n\n" << flush;

   exit ( 1 );

}

int j;
int v;
int V[2];
MtdIntFile s;
const int n3 = Nx*Ny*Nt;


s = *this;

s.set_to_zeroes();

int * in  = Data;
int * out = s.Data;

v = 0;

for (j=0; j<n3; ++j)  {

   if ( *in == n )  *out = 1;

   ++in;
   ++out;

}


V[0] = v;

s.set_volumes(1, V);

return ( s );

}


////////////////////////////////////////////////////////////////////////


   //
   //  Code for misc functions
   //


////////////////////////////////////////////////////////////////////////


MtdMoments MtdIntFile::calc_moments() const

{

int x, y, t, k, n;
MtdMoments m;

for (t=0; t<Nt; ++t)  {

   for (y=0; y<Ny; ++y)  {

      for (x=0; x<Nx; ++x)  {

         n = mtd_three_to_one(Nx, Ny, Nt, x, y, t);

         k = Data[n];

         if ( k )  {

            ++(m.N);

            m.Sx += x;
            m.Sy += y;
            m.St += t;

            m.Sxx += x*x;
            m.Syy += y*y;
            m.Stt += t*t;

            m.Sxy += x*y;
            m.Sxt += x*t;

            m.Syt += y*t;

         }   //  if k

      }   //  for x

   }   //  for y

}   //  for t

m.IsCentralized = false;

return ( m );

}


////////////////////////////////////////////////////////////////////////


MtdIntFile split(const MtdIntFile & mask, int & n_shapes)

{

int j, k;
int x, y, t;
int n_so_far;
int n_before, n_after;
int v, nc;
MtdIntFile f2;
MtdIntFile before;
MtdIntFile after;
MtdIntFile rv;
Mtd_Partition p;
const int zero_border_size = 2;
// int imin, imax;



   //
   //  find the partition
   //

n_shapes = 0;

n_so_far = 0;

f2 = mask.const_t_slice(0);

f2.zero_border(zero_border_size);

before = f2.split_const_t(n_before);

n_so_far = n_before;

for (j=0; j<n_before; ++j)  p.add_no_repeat(j + 1);

for (j=1; j<(mask.nt()); ++j)  {

   f2 = mask.const_t_slice(j);

   f2.zero_border(zero_border_size);

   after = f2.split_const_t(n_after);

   adjust_obj_numbers(after, n_so_far);

   for (k=0; k<n_after; ++k)  p.add_no_repeat(n_so_far + k + 1);

   n_so_far += n_after;

   find_overlap(before, after, p);

   before   = after;
   n_before = n_after;

}

// p.dump(cout);

n_shapes = p.n_elements();

   //
   //  use the partition
   //

rv = mask;

rv.set_to_zeroes();

n_so_far = 0;

for (t=0; t<(mask.nt()); ++t)  {

   f2 = mask.const_t_slice(t);

   f2.zero_border(zero_border_size);

   before = f2.split_const_t(n_before);

   adjust_obj_numbers(before, n_so_far);

   n_so_far += n_before;

   for (x=0; x<(mask.nx()); ++x)  {

      for (y=0; y<(mask.ny()); ++y)  {

         v = before(x, y, 0);

         if ( v == 0 ) continue;

         nc = p.which_class(v);

         if ( nc < 0 )  {

            cerr << "\n\n  split(const MtdIntFile &, int &) -> can't find cell!\n\n";

            exit ( 1 );

         }

         rv.put(nc + 1, x, y, t);

      }

   }

}

//  rv.calc_minmax(imin, imax);

   //
   //  done
   //

return ( rv );

}


////////////////////////////////////////////////////////////////////////


void adjust_obj_numbers(MtdIntFile & s, int delta)

{

if ( s.nt() != 1 )  {

   cerr << "\n\n  adjust_obj_numbers() -> not const-time slice!\n\n";

   exit ( 1 );

}

if ( delta == 0 )  return;

int x, y;
int v;
const int nx = s.nx();
const int ny = s.ny();

for (x=0; x<nx; ++x)  {

   for (y=0; y<ny; ++y)  {

      v = s(x, y, 0);

      if ( v != 0 )  s.put(v + delta, x, y, 0);

   }

}

return;

}


////////////////////////////////////////////////////////////////////////


void find_overlap(const MtdIntFile & before, const MtdIntFile & after, Mtd_Partition & p)

{

int x, y;
int a, b;
const int nx = before.nx();
const int ny = before.ny();


for (x=0; x<nx; ++x)  {

   for (y=0; y<ny; ++y)  {

      a = after(x, y, 0);

      if ( !a )  continue;

      b = before(x, y, 0);

      if ( !b )  continue;

      // cout << "merging values " << a << " and " << b << "\n";

      p.merge_values(a, b);

   }

}

return;

}


////////////////////////////////////////////////////////////////////////






