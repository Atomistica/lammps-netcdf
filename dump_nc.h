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

class DumpNC : public DumpCustom {
 public:
  DumpNC(class LAMMPS *, int, char **);
  ~DumpNC();

 private:
  int framei;                // current frame index
  int blocki;                // current block index
  int ndata;                 // number of data blocks to expect

  int ntypes;                // # of atom types
  char **typenames;	     // array of element names for each type

  int nnc;                   // # of netcdf properties
  int *ncdims;               // # of dimensions for each netcdf property
  int **nc2field;              // NC auxiliary to field mapping
  char **ncname;             // mangled names
  int *ncvar;                // NetCDF variables for auxiliary properties

  double **rbuf;             // buf of data lines for data lines rearrangement

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
  int cell_lengths_var;
  int cell_angles_var;

  void openfile();
  void write_header(int);
  void write_data(int, double *);

  int modify_param2(int, char **);

  void ncerr(int, int);
};

}

#endif
#endif
