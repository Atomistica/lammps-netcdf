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

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"
#include "universe.h"
#include "variable.h"

#include "dump_nc.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{INT,DOUBLE};  // same as in dump_custom.cpp

const char NC_FRAME_STR[]         = "frame";
const char NC_SPATIAL_STR[]       = "spatial";
const char NC_ATOM_STR[]          = "atom";
const char NC_CELL_SPATIAL_STR[]  = "cell_spatial";
const char NC_CELL_ANGULAR_STR[]  = "cell_angular";
const char NC_LABEL_STR[]         = "label";

const char NC_TIME_STR[]          = "time";
const char NC_CELL_ORIGIN_STR[]   = "cell_origin";
const char NC_CELL_LENGTHS_STR[]  = "cell_lengths";
const char NC_CELL_ANGLES_STR[]   = "cell_angles";

const char NC_UNITS_STR[]         = "units";
const char NC_SCALE_FACTOR_STR[]  = "scale_factor";

const int MAX_DIMS           = 10;
const int THIS_IS_A_FIX      = -1;
const int THIS_IS_A_COMPUTE  = -2;
const int THIS_IS_A_VARIABLE = -3;

/* ---------------------------------------------------------------------- */

#define NCERR(x) ncerr(x, __LINE__)

/* ---------------------------------------------------------------------- */

DumpNC::DumpNC(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  // arrays for data rearrangement

  rbuf = NULL;
  sort_flag = 1;
  sortcol = 0;

  if (multiproc)
    error->all(FLERR,"DumpNC: Multi-processor writes are not (yet) supported.");
  if (multifile)
    error->all(FLERR,"DumpNC: Multiple files are not supported.");

  perat = new nc_perat_t[nfield];

  for (int i = 0; i < nfield; i++) {
    perat[i].dims = 0;
  }

  n_perat = 0;
  for (int iarg = 5; iarg < narg; iarg++) {
    int i = iarg-5;
    int idim = 0;
    int ndims = 1;
    char mangled[1024];
    bool constant = false;

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
      strcpy(mangled, "atom_types");
      constant = true;
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
	  error->all(FLERR,"DumpNC: Missing ']' in dump command");
	*ptr = '\0';
	idim = ptr[1] - '1';
	ndims = THIS_IS_A_COMPUTE;
      }
    }
    else if (!strncmp(mangled, "f_", 2)) {
      char *ptr = strchr(mangled, '[');
      if (ptr) {
	if (mangled[strlen(mangled)-1] != ']')
	  error->all(FLERR,"DumpNC: Missing ']' in dump command");
	*ptr = '\0';
	idim = ptr[1] - '1';
	ndims = THIS_IS_A_FIX;
      }
    }

    // find mangled name
    int inc = -1;
    for (int j = 0; j < n_perat && inc < 0; j++) {
      if (!strcmp(perat[j].name, mangled)) {
	inc = j;
      }
    }

    if (inc < 0) {
      inc = n_perat;
      perat[inc].dims = ndims;
      if (ndims < 0) ndims = MAX_DIMS;
      for (int j = 0; j < ndims; j++) {
	perat[inc].field[j] = i;
      }
      strcpy(perat[inc].name, mangled);
      n_perat++;
    }

    perat[inc].constant = constant;
    perat[inc].ndumped = 0;
    perat[inc].field[idim] = i;
  }

  n_perframe = 0;
  perframe = NULL;

  n_buffer = 0;
  int_buffer = NULL;
  double_buffer = NULL;

  double_precision = false;
}

/* ---------------------------------------------------------------------- */

DumpNC::~DumpNC()
{
  if (me == 0 && singlefile_opened)
    NCERR( nc_close(ncid) );

  if (rbuf) memory->destroy(rbuf);

  delete [] perat;
  if (n_perframe > 0)
    delete [] perframe;

  if (int_buffer) memory->sfree(int_buffer);
  if (double_buffer) memory->sfree(double_buffer);
}

/* ---------------------------------------------------------------------- */

void DumpNC::openfile()
{
  // now the computes and fixes have been initialized, so we can query
  // for the size of vector quantities
  for (int i = 0; i < n_perat; i++) {
    if (perat[i].dims == THIS_IS_A_COMPUTE) {
      int j = field2index[perat[i].field[0]];
      if (!fix[j]->peratom_flag)
	error->all(FLERR,"DumpNC::init_style: compute does not provide per atom "
		   "data");
      perat[i].dims = compute[j]->size_peratom_cols;
      if (perat[i].dims > MAX_DIMS)
	error->all(FLERR,"DumpNC::init_style: perat[i].dims > MAX_DIMS");
    }
    else if (perat[i].dims == THIS_IS_A_FIX) {
      int j = field2index[perat[i].field[0]];
      if (!fix[j]->peratom_flag)
	error->all(FLERR,"DumpNC::init_style: fix does not provide per atom data");
      perat[i].dims = fix[j]->size_peratom_cols;
      if (perat[i].dims > MAX_DIMS)
	error->all(FLERR,"DumpNC::init_style: perat[i].dims > MAX_DIMS");
    }
  }

  // get total number of atoms
  ntotal = group->count(igroup);

  if (me == 0) {
    int dims[NC_MAX_VAR_DIMS];
    size_t index[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
    double d[1];

    if (singlefile_opened) return;
    singlefile_opened = 1;

    NCERR( nc_create(filename, NC_64BIT_OFFSET, &ncid) );
    
    // dimensions
    NCERR( nc_def_dim(ncid, NC_FRAME_STR, NC_UNLIMITED, &frame_dim) );
    NCERR( nc_def_dim(ncid, NC_SPATIAL_STR, 3, &spatial_dim) );
    NCERR( nc_def_dim(ncid, NC_ATOM_STR, ntotal, &atom_dim) );
    NCERR( nc_def_dim(ncid, NC_CELL_SPATIAL_STR, 3, &cell_spatial_dim) );
    NCERR( nc_def_dim(ncid, NC_CELL_ANGULAR_STR, 3, &cell_angular_dim) );
    NCERR( nc_def_dim(ncid, NC_LABEL_STR, 10, &label_dim) );

    // default variables
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
    NCERR( nc_def_var(ncid, NC_CELL_ORIGIN_STR, NC_DOUBLE, 2, dims,
		      &cell_origin_var) );
    NCERR( nc_def_var(ncid, NC_CELL_LENGTHS_STR, NC_DOUBLE, 2, dims,
		      &cell_lengths_var) );
    dims[0] = frame_dim;
    dims[1] = cell_angular_dim;
    NCERR( nc_def_var(ncid, NC_CELL_ANGLES_STR, NC_DOUBLE, 2, dims,
		      &cell_angles_var) );

    // variables specified in the input file
    dims[0] = frame_dim;
    dims[1] = atom_dim;
    dims[2] = spatial_dim;

    for (int i = 0; i < n_perat; i++) {
      nc_type xtype;

      // Type mangling
      if (vtype[perat[i].field[0]] == INT) {
	xtype = NC_INT;
      }
      else {
	if (double_precision)
	  xtype = NC_DOUBLE;
	else
	  xtype = NC_FLOAT;
      }

      if (perat[i].constant) {
	// this quantity will only be written once
	if (perat[i].dims == 3)
	  // this is needed to store x-, y- and z-coordinates
	  NCERR( nc_def_var(ncid, perat[i].name, xtype, 2, dims+1,
			    &perat[i].var) );
	else
	  NCERR( nc_def_var(ncid, perat[i].name, xtype, 1, dims+1,
			    &perat[i].var) );
      }
      else {
	if (perat[i].dims == 3)
	  // this is needed to store x-, y- and z-coordinates
	  NCERR( nc_def_var(ncid, perat[i].name, xtype, 3, dims,
			    &perat[i].var) );
	else
	  NCERR( nc_def_var(ncid, perat[i].name, xtype, 2, dims,
			    &perat[i].var) );
      }
    }

    // perframe variables
    for (int i = 0; i < n_perframe; i++) {
      NCERR( nc_def_var(ncid, perframe[i].name, NC_DOUBLE, 1, dims,
			&perframe[i].var) );
    }

    // attributes
    NCERR( nc_put_att_text(ncid, NC_GLOBAL, "Conventions",
			   5, "AMBER") );
    NCERR( nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion",
			   3, "1.0") );

    NCERR( nc_put_att_text(ncid, NC_GLOBAL, "program",
			   6, "LAMMPS") );
    NCERR( nc_put_att_text(ncid, NC_GLOBAL, "programVersion",
			   strlen(universe->version), universe->version) );

    // units
    if (!strcmp(update->unit_style, "lj")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     2, "lj") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     2, "lj") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     2, "lj") );
    }
    else if (!strcmp(update->unit_style, "real")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     11, "femtosecond") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     8, "Angstrom") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     8, "Angstrom") );
    }
    else if (!strcmp(update->unit_style, "metal")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     10, "picosecond") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     8, "Angstrom") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     8, "Angstrom") );
    }
    else if (!strcmp(update->unit_style, "si")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     6, "second") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     5, "meter") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     5, "meter") );
    }
    else if (!strcmp(update->unit_style, "cgs")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     6, "second") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     10, "centimeter") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     10, "centimeter") );
    }
    else if (!strcmp(update->unit_style, "electron")) {
      NCERR( nc_put_att_text(ncid, time_var, NC_UNITS_STR,
			     11, "femtosecond") );
      NCERR( nc_put_att_text(ncid, cell_origin_var, NC_UNITS_STR,
			     4, "Bohr") );
      NCERR( nc_put_att_text(ncid, cell_lengths_var, NC_UNITS_STR,
			     4, "Bohr") );
    }
    else {
      char errstr[1024];
      sprintf(errstr, "DumpNC: Unsupported unit style '%s'",
	      update->unit_style);
      error->all(FLERR,errstr);
    }

    NCERR( nc_put_att_text(ncid, cell_angles_var, NC_UNITS_STR,
			   6, "degree") );

    d[0] = update->dt;
    NCERR( nc_put_att_double(ncid, time_var, NC_SCALE_FACTOR_STR,
			     NC_DOUBLE, 1, d) );
    d[0] = 1.0;
    NCERR( nc_put_att_double(ncid, cell_origin_var, NC_SCALE_FACTOR_STR,
			     NC_DOUBLE, 1, d) );
    d[0] = 1.0;
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
}

/* ---------------------------------------------------------------------- */

void DumpNC::write_header(bigint n)
{
  size_t start[2];

  start[0] = framei;
  start[1] = 0;

  if (me == 0) {
    size_t count[2];
    double time, cell_origin[3], cell_lengths[3], cell_angles[3];

    time = update->ntimestep;
    if (domain->triclinic == 0) {
      cell_origin[0] = domain->boxlo[0];
      cell_origin[1] = domain->boxlo[1];
      cell_origin[2] = domain->boxlo[2];

      cell_lengths[0] = domain->xprd;
      cell_lengths[1] = domain->yprd;
      cell_lengths[2] = domain->zprd;

      cell_angles[0] = 90;
      cell_angles[1] = 90;
      cell_angles[2] = 90;
    }
    else {
      double cosalpha, cosbeta, cosgamma;
      double *h = domain->h;

      cell_origin[0] = domain->boxlo[0];
      cell_origin[1] = domain->boxlo[1];
      cell_origin[2] = domain->boxlo[2];

      cell_lengths[0] = domain->xprd;
      cell_lengths[1] = sqrt(h[1]*h[1]+h[5]*h[5]);
      cell_lengths[2] = sqrt(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]);

      cosalpha = (h[5]*h[4]+h[1]*h[3])/
	sqrt((h[1]*h[1]+h[5]*h[5])*(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]));
      cosbeta = h[4]/sqrt(h[2]*h[2]+h[3]*h[3]+h[4]*h[4]);
      cosgamma = h[5]/sqrt(h[1]*h[1]+h[5]*h[5]);

      cell_angles[0] = acos(cosalpha)*180.0/MY_PI;
      cell_angles[1] = acos(cosbeta)*180.0/MY_PI;
      cell_angles[2] = acos(cosgamma)*180.0/MY_PI;
    }

    count[0] = 1;
    count[1] = 3;
    NCERR( nc_put_var1_double(ncid, time_var, start, &time) );
    NCERR( nc_put_vara_double(ncid, cell_origin_var, start, count,
			      cell_origin) );
    NCERR( nc_put_vara_double(ncid, cell_lengths_var, start, count,
			      cell_lengths) );
    NCERR( nc_put_vara_double(ncid, cell_angles_var, start, count,
			      cell_angles) );
  }

  for (int i = 0; i < n_perframe; i++) {
    double data;
    int j = perframe[i].index;
    int idim = perframe[i].dim;

    if (perframe[i].type == THIS_IS_A_COMPUTE) {
      if (idim >= 0) {
	modify->compute[j]->compute_vector();
	data = modify->compute[j]->vector[idim];
      }
      else
	data = modify->compute[j]->compute_scalar();
    }
    else if (perframe[i].type == THIS_IS_A_FIX) {
      if (idim >= 0) {
	data = modify->fix[j]->compute_vector(idim);
      }
      else
	data = modify->fix[j]->compute_scalar();
    }
    else if (perframe[i].type == THIS_IS_A_VARIABLE) {
      j = input->variable->find(perframe[i].id);
      data = input->variable->compute_equal(j);
    }

    if (me == 0)
      NCERR( nc_put_var1_double(ncid, perframe[i].var, start, &data) );
  }

  ndata = n;
  blocki = 0;
}


/* ----------------------------------------------------------------------
   write data lines to file in a block-by-block style
   write head of block (mass & element name) only if has atoms of the type
------------------------------------------------------------------------- */

void DumpNC::write_data(int n, double *mybuf)
{
  size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
  ptrdiff_t stride[NC_MAX_VAR_DIMS];

  if (!int_buffer) {
    n_buffer = n;
    int_buffer = (int *)
      memory->smalloc(n*sizeof(int), "DumpNC::int_buffer");
    double_buffer = (double *)
      memory->smalloc(n*sizeof(double), "DumpNC::double_buffer");
  }

  if (n > n_buffer) {
    n_buffer = n;
    int_buffer = (int *)
      memory->srealloc(int_buffer, n*sizeof(int), "DumpNC::int_buffer");
    double_buffer = (double *)
      memory->srealloc(double_buffer, n*sizeof(double),
		       "DumpNC::double_buffer");
  }

  start[0] = framei;
  start[1] = blocki;
  start[2] = 0;

  count[0] = 1;
  count[1] = n;
  count[2] = 1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 3;

  for (int i = 0; i < n_perat; i++) {
    int iaux = perat[i].field[0];

    if (vtype[iaux] == INT) {
      // integers
      if (perat[i].dims == 3) {

	for (int idim = 0; idim < 3; idim++) {
	  iaux = perat[i].field[idim];

	  for (int j = 0; j < n; j++, iaux+=size_one) {
	    int_buffer[j] = mybuf[iaux];
	  }
      
	  start[2] = idim;

	  if (perat[i].constant) {
	    if (perat[i].ndumped < ntotal) {
	      NCERR( nc_put_vars_int(ncid, perat[i].var,
				     start+1, count+1, stride+1,
				     int_buffer) );
	      perat[i].ndumped += n;
	    }
	  }
	  else
	    NCERR( nc_put_vars_int(ncid, perat[i].var, start, count, stride,
				   int_buffer) );
	}
      }
      else {
	for (int j = 0; j < n; j++, iaux+=size_one) {
	  int_buffer[j] = mybuf[iaux];
	}

	if (perat[i].constant) {
	  if (perat[i].ndumped < ntotal) {
	    NCERR( nc_put_vara_int(ncid, perat[i].var, start+1, count+1,
				   int_buffer) );
	    perat[i].ndumped += n;
	  }
	}
	else
	  NCERR( nc_put_vara_int(ncid, perat[i].var, start, count,
				 int_buffer) );
      }
    }
    else {
      // doubles
      if (perat[i].dims == 3) {

	for (int idim = 0; idim < 3; idim++) {
	  iaux = perat[i].field[idim];

	  for (int j = 0; j < n; j++, iaux+=size_one) {
	    double_buffer[j] = mybuf[iaux];
	  }
      
	  start[2] = idim;

	  if (perat[i].constant) {
	    if (perat[i].ndumped < ntotal) {
	      NCERR( nc_put_vars_double(ncid, perat[i].var,
					start+1, count+1, stride+1,
					double_buffer) );
	      perat[i].ndumped += n;
	    }
	  }
	  else
	    NCERR( nc_put_vars_double(ncid, perat[i].var, start, count, stride,
				      double_buffer) );
	}
      }
      else {
	for (int j = 0; j < n; j++, iaux+=size_one) {
	  double_buffer[j] = mybuf[iaux];
	}

	if (perat[i].constant) {
	  if (perat[i].ndumped < ntotal) {
	    NCERR( nc_put_vara_double(ncid, perat[i].var, start+1, count+1,
				      double_buffer) );
	    perat[i].ndumped += n;
	  }
	}
	else
	  NCERR( nc_put_vara_double(ncid, perat[i].var, start, count,
				    double_buffer) );
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

int DumpNC::modify_param(int narg, char **arg)
{
  int iarg = 0;
  if (strcmp(arg[iarg],"double_precision") == 0) {
    double_precision = true;
    iarg++;
  }
  if (strcmp(arg[iarg],"global") == 0) {
    // "perframe" quantities, i.e. not per-atom stuff

    iarg++;

    n_perframe = narg-iarg;
    perframe = new nc_perframe_t[n_perframe];

    for (int i = 0; iarg < narg; iarg++, i++) {
      int n;
      char *suffix;

      n = strlen(arg[iarg]);

      if (n > 2) {
	suffix = new char[n-1];
	strcpy(suffix, arg[iarg]+2);
      }
      else {
	char errstr[1024];
	sprintf(errstr, "DumpNC::modify_param: perframe quantity '%s' must "
		"be compute, fix or variable", arg[iarg]);
	error->all(FLERR,errstr);
      }

      if (!strncmp(arg[iarg], "c_", 2)) {
	int idim = -1;
	char *ptr = strchr(suffix, '[');

	if (ptr) {
	  if (suffix[strlen(suffix)-1] != ']')
	    error->all(FLERR,"DumpNC: Missing ']' in dump modify command");
	  *ptr = '\0';
	  idim = ptr[1] - '1';
	}

	n = modify->find_compute(suffix);
	if (n < 0)
	  error->all(FLERR,"Could not find dump modify compute ID");
	if (modify->compute[n]->peratom_flag != 0)
	  error->all(FLERR,"Dump modify compute ID computes per-atom info");
	if (idim >= 0 && modify->compute[n]->vector_flag == 0)
	  error->all(FLERR,"Dump modify compute ID does not compute vector");
	if (idim < 0 && modify->compute[n]->scalar_flag == 0)
	  error->all(FLERR,"Dump modify compute ID does not compute scalar");

	perframe[i].type = THIS_IS_A_COMPUTE;
	perframe[i].dim = idim;
	perframe[i].index = n;
	strcpy(perframe[i].name, arg[iarg]);
      }
      else if (!strncmp(arg[iarg], "f_", 2)) {
	int idim = -1;
	char *ptr = strchr(suffix, '[');

	if (ptr) {
	  if (suffix[strlen(suffix)-1] != ']')
	    error->all(FLERR,"DumpNC: Missing ']' in dump modify command");
	  *ptr = '\0';
	  idim = ptr[1] - '1';
	}

	n = modify->find_fix(suffix);
	if (n < 0)
	  error->all(FLERR,"Could not find dump modify fix ID");
	if (modify->fix[n]->peratom_flag != 0)
	  error->all(FLERR,"Dump modify fix ID computes per-atom info");
	if (idim >= 0 && modify->fix[n]->vector_flag == 0)
	  error->all(FLERR,"Dump modify fix ID does not compute vector");
	if (idim < 0 && modify->fix[n]->scalar_flag == 0)
	  error->all(FLERR,"Dump modify fix ID does not compute vector");

	perframe[i].type = THIS_IS_A_FIX;
	perframe[i].dim = idim;
	perframe[i].index = n;
	strcpy(perframe[i].name, arg[iarg]);
      }
      else if (!strncmp(arg[iarg], "v_", 2)) {
	n = input->variable->find(suffix);
	if (n < 0)
	  error->all(FLERR,"Could not find dump modify variable ID");
	if (!input->variable->equalstyle(n))
	  error->all(FLERR,"Dump modify variable must be of style equal");

	perframe[i].type = THIS_IS_A_VARIABLE;
	perframe[i].dim = 1;
	perframe[i].index = n;
	strcpy(perframe[i].name, arg[iarg]);
	strcpy(perframe[i].id, suffix);
      }
      else {
	char errstr[1024];
	sprintf(errstr, "DumpNC::modify_param: perframe quantity '%s' must "
		"be compute, fix or variable", arg[iarg]);
	error->all(FLERR,errstr);
      }

      delete [] suffix;
    }

    return narg;
  } else return 0;
}

/* ---------------------------------------------------------------------- */

void DumpNC::write_prmtop()
{
  char fn[1024];
  char tmp[81];
  FILE *f;

  strcpy(fn, filename);
  strcat(fn, ".prmtop");

  f = fopen(fn, "w");
  fprintf(f, "%%VERSION  LAMMPS\n");
  fprintf(f, "%%FLAG TITLE\n");
  fprintf(f, "%%FORMAT(20a4)\n");
  memset(tmp, ' ', 76);
  tmp[76] = '\0';
  fprintf(f, "NASN%s\n", tmp);

  fprintf(f, "%%FLAG POINTERS\n");
  fprintf(f, "%%FORMAT(10I8)\n");
  fprintf(f, "%8i", ntotal);
  for (int i = 0; i < 11; i++)
    fprintf(f, "%8i", 0);
  fprintf(f, "\n");
  for (int i = 0; i < 12; i++)
    fprintf(f, "%8i", 0);
  fprintf(f, "\n");
  for (int i = 0; i < 6; i++)
    fprintf(f, "%8i", 0);
  fprintf(f, "\n");

  fprintf(f, "%%FLAG ATOM_NAME\n");
  fprintf(f, "%%FORMAT(20a4)\n");
  for (int i = 0; i < ntotal; i++) {
    fprintf(f, "%4s", "He");
    if ((i+1) % 20 == 0)
      fprintf(f, "\n");
  }

  fprintf(f, "%%FLAG CHARGE\n");
  fprintf(f, "%%FORMAT(5E16.5)\n");
  for (int i = 0; i < ntotal; i++) {
    fprintf(f, "%16.5e", 0.0);
    if ((i+1) % 5 == 0)
      fprintf(f, "\n");
  }

  fprintf(f, "%%FLAG MASS\n");
  fprintf(f, "%%FORMAT(5E16.5)\n");
  for (int i = 0; i < ntotal; i++) {
    fprintf(f, "%16.5e", 1.0);
    if ((i+1) % 5 == 0)
      fprintf(f, "\n");
  }
  fclose(f);
}

/* ---------------------------------------------------------------------- */

void DumpNC::ncerr(int err, int line)
{
  if (err != NC_NOERR) {
    char errstr[1024];
    sprintf(errstr, "NetCDF failed with error '%s' in line %i of %s.",
	    nc_strerror(err), line, __FILE__);
    error->one(FLERR,errstr);
  }
}
