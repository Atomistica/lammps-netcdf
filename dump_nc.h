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

#ifdef DUMP_CLASS

DumpStyle(nc,DumpNC)

#else

#ifndef LMP_DUMP_NC_H
#define LMP_DUMP_NC_H

#include "dump_custom.h"

namespace LAMMPS_NS {

const int NC_FIELD_NAME_MAX = 100;

class DumpNC : public DumpCustom {
 public:
  DumpNC(class LAMMPS *, int, char **);
  ~DumpNC();

 private:
  // per-atoms quantities (positions, velocities, etc.)
  struct nc_perat_t {
    int dims;                     // number of dimensions
    int field[3];                 // field indices corresponding x,y,z-comp.
    char name[NC_FIELD_NAME_MAX]; // field name
    int var;                      // NetCDF variable

    bool constant;                // is this property per file (not per frame)
    int ndumped;                  // number of enties written for this prop.
  };

  typedef void (DumpNC::*funcptr_t)(void *);

  // per-frame quantities (variables, fixes or computes)
  struct nc_perframe_t {
    char name[NC_FIELD_NAME_MAX]; // field name
    int var;                      // NetCDF variable
    int type;                     // variable, fix, compute or callback
    int index;                    // index in fix/compute list
    funcptr_t compute;            // compute function
    int dim;                      // dimension
    char id[NC_FIELD_NAME_MAX];   // variable id
  };

  int framei;                  // current frame index
  int blocki;                  // current block index
  int ndata;                   // number of data blocks to expect

  int ntotal;                  // # of atoms

  int n_perat;                 // # of netcdf per-atom properties
  nc_perat_t *perat;           // per-atom properties

  int n_perframe;              // # of global netcdf (not per-atom) fix props
  nc_perframe_t *perframe;     // global properties

  double **rbuf;               // buf of data lines for data lines rearrangemnt

  bool double_precision;       // write everything as double precision

  int n_buffer;                // size of buffer
  int *int_buffer;             // buffer for passing data to netcdf
  double *double_buffer;       // buffer for passing data to netcdf

  int ncid;

  int frame_dim;
  int spatial_dim;
  int atom_dim;
  int cell_spatial_dim;
  int cell_angular_dim;
  int label_dim;

  int spatial_var;
  int cell_spatial_var;
  int cell_angular_var;

  int time_var;
  int cell_origin_var;
  int cell_lengths_var;
  int cell_angles_var;

  void openfile();
  void write_header(bigint);
  void write_data(int, double *);
  void write_prmtop();

  virtual int modify_param(int, char **);

  void ncerr(int, int);

  void compute_step(void *);
  void compute_elapsed(void *);
  void compute_elapsed_long(void *);
};

}

#endif
#endif
