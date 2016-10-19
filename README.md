NOTICE
======

This repository has been merged into the main LAMMPS repository. Development here
has ceased, please head to the LAMMPS repository for a current version. `dump nc`
and `dump nc/mpiio` can be found in the `USER-NC-DUMP` directory:

https://github.com/lammps/lammps/tree/master/src/USER-NC-DUMP


LAMMPS-NETCDF
=============

The files in this directory are a user-contributed package for LAMMPS.

The person who created these files is Lars Pastewka at
Karlsruhe Institute of Technology (lars.pastewka@kit.edu).
Contact him directly if you have questions.

Lars Pastewka  
Institute for Applied Materials (IAM)  
Karlsruhe Institute of Technology (KIT)  
Kaiserstraße 12, 76131 Karlsruhe  
e-mail: lars.pastewka@kit.edu

PACKAGE DESCRIPTION
-------------------

This is a LAMMPS (http://lammps.sandia.gov/) dump style for output into a NetCDF
database. The database format follows the AMBER NetCDF trajectory convention
(http://ambermd.org/netcdf/nctraj.xhtml), but includes extensions to this
convention. These extension are:
* A variable "cell_origin" (of dimension "frame", "cell_spatial") that contains
  the bottom left corner of the simulation cell.
* Any number of additional variables corresponding to per atom scalar, vector
  or tensor quantities available within LAMMPS. Tensor quantities are written in
  Voigt notation. An additional dimension "Voigt" of length 6 is created for
  this purpose.
* Possibility to output to an HDF5 database.

NetCDF files can be directly visualized with the following tools:
* Ovito (http://www.ovito.org/). Ovito supports the AMBER convention and all of
  the above extensions.
* VMD (http://www.ks.uiuc.edu/Research/vmd/).
* AtomEye (http://www.libatoms.org/). The libAtoms version of AtomEye contains
  a NetCDF reader that is not present in the standard distribution of AtomEye.

Syntax:

> dump ID group-ID nc|nc/mpiio N file args  
> dump_modify ID keyword values

with

args = list of atom attributes  
keyword = append or double or global
* append yes: append output to existing NetCDF file
* append yes at frame: append out to existing NetCDF file and start writing to
    frame given; if negative frame is counted from the end of file.
* double yes|no: output data as double instead of single precision
* global = list of global (not per atom, but per frame) quantities
    can be variables, compute or fix data prefixed with v_, c_ and f_,
    respectively.

The list of atom attributes is identical to the 'custom' dump style. The 'nc'
dump style uses the standard NetCDF library and collects all data onto a single
processor before writing. The 'nc/mpiio' dump style uses parallel-netcdf and
MPI-IO and has better performance on a larger number of processors. Note that
'nc' outputs all atoms sorted by tag while 'nc/mpiio' outputs in order of the
MPI rank.

Example 1:

> dump 1 all nc 100 traj.nc type x y z vx vy vz  
> dump_modify 1 append yes at -1 global c_thermo_pe c_thermo_temp c_thermo_press

Example 2:

> dump 1 all nc/mpiio 1000 traj.nc id type x y z

INSTALLATION
------------

In your LAMMPS src directory type:

> git clone https://github.com/pastewka/lammps-netcdf.git USER-DUMP-NC  
> make yes-user-dump-nc

Note that LAMMPS needs to be linked to NetCDF (for 'dump nc') and Parallel
netCDF (for 'dump nc/mpiio'). This requires a modification of your favorite
makefile. Please add

> EXTRA_INC += $(shell nc-config --cflags) -I/path/to/parallel-netcdf/include  
> EXTRA_LIB += $(shell nc-config --libs) -L/path/to/parallel-netcdf/lib -lpnetcdf

to the respective EXTRA_INC, EXTRA_LIB section of the makefile.

OTHER NOTES
-----------

NetCDF: http://www.unidata.ucar.edu/software/netcdf/  
Parallel netCDF: http://trac.mcs.anl.gov/projects/parallel-netcdf/

This is package is known to work with LAMMPS 15Aug14.
