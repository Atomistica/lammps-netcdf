/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Lars Pastewka (Johns Hopkins University)
------------------------------------------------------------------------- */

#include <netcdf.h>

#include "stdlib.h"
#include "string.h"
#include "dump_nc.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "fix.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace LAMMPS_NS;

enum{INT,DOUBLE};  // same as in dump_custom.cpp

const char NC_FRAME_STR[]         = "frame";
const char NC_SPATIAL_STR[]       = "spatial";
const char NC_ATOM_STR[]          = "atom";
const char NC_CELL_SPATIAL_STR[]  = "cell_spatial";
const char NC_CELL_ANGULAR_STR[]  = "cell_angular";
const char NC_LABEL_STR[]         = "label";

const char NC_TIME_STR[]          = "time";
const char NC_CELL_LENGTHS_STR[]  = "cell_lengths";
const char NC_CELL_ANGLES_STR[]   = "cell_angles";

const char NC_UNITS_STR[]         = "units";
const char NC_SCALE_FACTOR_STR[]  = "scale_factor";

const int MAX_DIMS           = 10;
const int THIS_IS_A_FIX      = -1;
const int THIS_IS_A_COMPUTE  = -2;

/* ---------------------------------------------------------------------- */

#define NCERR(x) ncerr(x, __LINE__)

/* ---------------------------------------------------------------------- */

DumpNC::DumpNC(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  ntypes = atom->ntypes;
  typenames = NULL;

  // arrays for data rearrangement

  rbuf = NULL;

  if (multiproc)
    error->all("DumpNC: Multi-processor writes are not supported.");
  if (multifile)
    error->all("DumpNC: Multiple files are not supported.");

  ncname = new char*[nfield];
  ncdims = new int[nfield];
  nc2field = new int*[nfield];

  nnc = 0;
  for (int iarg = 5; iarg < narg; iarg++) {
    int i = iarg-5;
    int idim = 0;
    int ndims = 1;
    char mangled[1024];

    strcpy(mangled, arg[iarg]);

    // name mangling
    // in the AMBER specification
    if (!strcmp(mangled, "x") || !strcmp(mangled, "y") ||
	!strcmp(mangled, "z")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      strcpy(mangled, "coordinates");
    }
    else if (!strcmp(mangled, "vx") || !strcmp(mangled, "vy") ||
	     !strcmp(mangled, "vz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      strcpy(mangled, "velocities");
    }
    // extensions to the AMBER specification
    else if (!strcmp(mangled, "type")) {
      strcpy(mangled, "Z");
    }
    else if (!strcmp(mangled, "xs") || !strcmp(mangled, "ys") ||
	!strcmp(mangled, "zs")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      strcpy(mangled, "scaled_coordinates");
    }
    else if (!strcmp(mangled, "xu") || !strcmp(mangled, "yu") ||
	!strcmp(mangled, "zu")) {
      idim = mangled[0] - 'x';
      ndims = 3;
      strcpy(mangled, "unwrapped_coordinates");
    }
    else if (!strcmp(mangled, "fx") || !strcmp(mangled, "fy") ||
	     !strcmp(mangled, "fz")) {
      idim = mangled[1] - 'x';
      ndims = 3;
      strcpy(mangled, "forces");
    }
    else if (!strcmp(mangled, "mux") || !strcmp(mangled, "muy") ||
	     !strcmp(mangled, "muz")) {
      idim = mangled[2] - 'x';
      ndims = 3;
      strcpy(mangled, "mu");
    }
    else if (!strncmp(mangled, "c_", 2)) {
      char *ptr = strchr(mangled, '[');
      if (ptr) {
	if (mangled[strlen(mangled)-1] != ']')
	  error->all("DumpNC: Missing ']' in dump command");
	*ptr = '\0';
	idim = ptr[1] - '1';
	ndims = THIS_IS_A_COMPUTE;
      }
    }
    else if (!strncmp(mangled, "f_", 2)) {
      char *ptr = strchr(mangled, '[');
      if (ptr) {
	if (mangled[strlen(mangled)-1] != ']')
	  error->all("DumpNC: Missing ']' in dump command");
	*ptr = '\0';
	idim = ptr[1] - '1';
	ndims = THIS_IS_A_FIX;
      }
    }

    // find mangled name
    int inc = -1;
    for (int j = 0; j < nnc && inc < 0; j++) {
      if (!strcmp(ncname[j], mangled)) {
	inc = j;
      }
    }

    if (inc < 0) {
      inc = nnc;
      ncname[inc] = new char[strlen(mangled)+1];
      ncdims[inc] = ndims;
      if (ndims < 0) ndims = MAX_DIMS;
      nc2field[inc] = new int[ndims];
      for (int j = 0; j < ndims; j++) {
	nc2field[inc][j] = i;
      }
      strcpy(ncname[inc], mangled);
      nnc++;
    }

    nc2field[inc][idim] = i;
  }

  ncvar = new int[nnc];
}

/* ---------------------------------------------------------------------- */

DumpNC::~DumpNC()
{
  NCERR( nc_close(ncid) );

  if (typenames) {
    for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
    delete [] typenames;
  }

  if (rbuf) memory->destroy_2d_double_array(rbuf);

  for (int i = 0; i < nnc; i++) {
    delete [] ncname[i];
    delete [] nc2field[i];
  }
  delete [] ncvar;
  delete [] ncname;
  delete [] ncdims;
  delete [] nc2field;
}

/* ---------------------------------------------------------------------- */

void DumpNC::init_style()
{
  DumpCustom::init_style();

  // now the computes and fixes have been initialized, so we can query
  // for the size of vector quantities
  for (int i = 0; i < nnc; i++) {
    if (ncdims[i] == THIS_IS_A_COMPUTE) {
      ncdims[i] = compute[field2index[nc2field[i][0]]]->size_vector;
      if (ncdims[i] > MAX_DIMS)
	error->all("DumpNC::init_style: ncdims[i] > MAX_DIMS");
    }
    else if (ncdims[i] == THIS_IS_A_FIX) {
      ncdims[i] = fix[field2index[nc2field[i][0]]]->size_vector;
      if (ncdims[i] > MAX_DIMS)
	error->all("DumpNC::init_style: ncdims[i] > MAX_DIMS");
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpNC::openfile()
{
  int ntotal;
  int dims[NC_MAX_VAR_DIMS];
  size_t index[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  double d[1];

  if (singlefile_opened) return;
  singlefile_opened = 1;

  /*
   * Get total number of atoms
   */
  MPI_Allreduce(&atom->nlocal, &ntotal, 1, MPI_INT, MPI_SUM, world);

  NCERR( nc_create(filename, NC_64BIT_OFFSET, &ncid) );

  /*
   * Dimensions
   */
  NCERR( nc_def_dim(ncid, NC_FRAME_STR, NC_UNLIMITED, &frame_dim) );
  NCERR( nc_def_dim(ncid, NC_SPATIAL_STR, 3, &spatial_dim) );
  NCERR( nc_def_dim(ncid, NC_ATOM_STR, ntotal, &atom_dim) );
  NCERR( nc_def_dim(ncid, NC_CELL_SPATIAL_STR, 3, &cell_spatial_dim) );
  NCERR( nc_def_dim(ncid, NC_CELL_ANGULAR_STR, 3, &cell_angular_dim) );
  NCERR( nc_def_dim(ncid, NC_LABEL_STR, 10, &label_dim) );

  /*
   * Variables
   */
  dims[0] = spatial_dim;
  NCERR( nc_def_var(ncid, NC_SPATIAL_STR, NC_CHAR, 1, dims, &spatial_var) );
  NCERR( nc_def_var(ncid, NC_CELL_SPATIAL_STR, NC_CHAR, 1, dims,
		    &cell_spatial_var) );
  dims[0] = spatial_dim;
  dims[1] = label_dim;
  NCERR( nc_def_var(ncid, NC_CELL_ANGULAR_STR, NC_CHAR, 2, dims,
		    &cell_angular_var) );

  dims[0] = frame_dim;
  NCERR( nc_def_var(ncid, NC_TIME_STR, NC_DOUBLE, 1, dims, &time_var) );
  dims[0] = frame_dim;
  dims[1] = cell_spatial_dim;
  NCERR( nc_def_var(ncid, NC_CELL_LENGTHS_STR, NC_DOUBLE, 2, dims,
		    &cell_lengths_var) );
  dims[0] = frame_dim;
  dims[1] = cell_angular_dim;
  NCERR( nc_def_var(ncid, NC_CELL_ANGLES_STR, NC_DOUBLE, 2, dims,
		    &cell_angles_var) );

  /*
   * Dynamic variables
   */
  dims[0] = frame_dim;
  dims[1] = atom_dim;
  dims[2] = spatial_dim;
  for (int i = 0; i < nnc; i++) {
    nc_type xtype;

    /*
     * Type mangling
     */
    if (vtype[nc2field[i][0]] == INT) {
      xtype = NC_INT;
    }
    else {
      xtype = NC_DOUBLE;
    }

    if (ncdims[i] == 3) {
      // this is needed to store x-, y- and z-coordinates
      NCERR( nc_def_var(ncid, ncname[i], xtype, 3, dims, &ncvar[i]) );
    }
    else {
      NCERR( nc_def_var(ncid, ncname[i], xtype, 2, dims, &ncvar[i]) );
    }
  }
  /*
  dims[0] = atom_dim;
  dims[1] = frame_dim;
  NCERR( nc_def_var( ncid, p%data%name_real(i), NC_DOUBLE, 2, dims,
				  real_var(i) ) );


  dims[0] = spatial_dim;
  dims[1] = atom_dim;
  dims[2] = frame_dim;
  NCERR( nc_def_var( ncid, p%data%name_real3x3(i), NC_DOUBLE,
				  3, dims, real3x3_var(i) ) );
  */

  /*
   * Attributes
   */

  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "Conventions",
			 5, "AMBER") );
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion",
			 3, "1.0") );

  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "program",
			 6, "LAMMPS") );
  NCERR( nc_put_att_text(ncid, NC_GLOBAL, "programVersion",
			 strlen(universe->version), universe->version) );

  /*
   * Fixme! Proper units here.
   */
  NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			 10, "picosecond") );
  NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			 8, "angstrom") );
  NCERR( nc_put_att_text(ncid, cell_angles_var, NC_UNITS_STR,
			 6, "degree") );

  d[0] = 1.0;
  NCERR( nc_put_att_double(ncid, time_var, NC_SCALE_FACTOR_STR,
			   NC_DOUBLE, 1, d) );
  NCERR( nc_put_att_double(ncid, cell_lengths_var, NC_SCALE_FACTOR_STR,
			   NC_DOUBLE, 1, d) );

  /*
   * Finished with definition
   */

  NCERR( nc_enddef(ncid) );

  /*
   * Write label variables
   */

  NCERR( nc_put_var_text(ncid, spatial_var, "xyz") );
  NCERR( nc_put_var_text(ncid, cell_spatial_var, "abc") );
  index[0] = 0;
  index[1] = 0;
  count[0] = 1;
  count[1] = 5;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "alpha") );
  index[0] = 1;
  count[1] = 4;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "beta") );
  index[0] = 2;
  count[1] = 5;
  NCERR( nc_put_vara_text(ncid, cell_angular_var, index, count, "gamma") );

  framei = 0;
}

/* ---------------------------------------------------------------------- */

void DumpNC::write_header(int n)
{
  ndata = n;
  blocki = 0;
}

/* ----------------------------------------------------------------------
   write data lines to file in a block-by-block style
   write head of block (mass & element name) only if has atoms of the type
------------------------------------------------------------------------- */

void DumpNC::write_data(int n, double *mybuf)
{
  int int_data[n];
  double double_data[n];

  size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  ptrdiff_t stride[NC_MAX_VAR_DIMS];

  start[0] = framei;
  start[1] = blocki;
  start[2] = 0;

  count[0] = 1;
  count[1] = n;
  count[2] = 1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 3;

  for (int i = 0; i < nnc; i++) {
    int iaux = nc2field[i][0];

    if (vtype[iaux] == INT) {
      // integers
      if (ncdims[i] == 3) {

	for (int idim = 0; idim < 3; idim++) {
	  iaux = nc2field[i][idim];

	  for (int j = 0; j < n; j++, iaux+=size_one) {
	    int_data[j] = mybuf[iaux];
	  }
      
	  start[2] = idim;
	  NCERR( nc_put_vars_int(ncid, ncvar[i], start, count, stride,
				 int_data) );
	}
      }
      else {
	for (int j = 0; j < n; j++, iaux+=size_one) {
	  int_data[j] = mybuf[iaux];
	}

	NCERR( nc_put_vara_int(ncid, ncvar[i], start, count,
			       int_data) );
      }
    }
    else {
      // doubles
      if (ncdims[i] == 3) {

	for (int idim = 0; idim < 3; idim++) {
	  iaux = nc2field[i][idim];

	  for (int j = 0; j < n; j++, iaux+=size_one) {
	    double_data[j] = mybuf[iaux];
	  }
      
	  start[2] = idim;
	  NCERR( nc_put_vars_double(ncid, ncvar[i], start, count, stride,
				    double_data) );
	}
      }
      else {
	for (int j = 0; j < n; j++, iaux+=size_one) {
	  double_data[j] = mybuf[iaux];
	}

	NCERR( nc_put_vara_double(ncid, ncvar[i], start, count,
				   double_data) );
      }
    }
  }

  blocki += n;

  if (blocki >= ndata) {
    NCERR( nc_sync(ncid) );
    framei++;
  }
}

/* ---------------------------------------------------------------------- */

int DumpNC::modify_param2(int narg, char **arg)
{
  if (strcmp(arg[0],"element") == 0) {
    if (narg != ntypes+1)
      error->all("Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[3];
      if (strlen(arg[itype]) >= 3)
	error->all("Illegal chemical element names");
      strcpy(typenames[itype],arg[itype]);
    }
    return ntypes+1;

  } else return 0;
}

/* ---------------------------------------------------------------------- */

void DumpNC::ncerr(int err, int line)
{
  if (err != NC_NOERR) {
    char errstr[1024];
    sprintf(errstr, "NetCDF failed with error '%s' in line %i of %s.",
	    nc_strerror(err), line, __FILE__);
    error->one(errstr);
  }
}
