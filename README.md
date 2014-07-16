LAMMPS-NETCDF
=============

The files in this directory are a user-contributed package for LAMMPS.

The person who created these files is Lars Pastewka at
Fraunhofer IWM (lars.pastewka@iwm.fraunhofer.de).
Contact him directly if you have questions.

Lars Pastewka  
Fraunhofer-Institut für Werkstoffmechanik IWM  
Wöhlerstraße 11, 79108 Freiburg  
e-mail: lars.pastewka@iwm.fraunhofer.de  

PACKAGE DESCRIPTION
-------------------

This is a LAMMPS (http://lammps.sandia.gov/) dump style for output into a NetCDF
database. The database format follows the AMBER NetCDF trajectory convention
(http://ambermd.org/netcdf/nctraj.xhtml), but includes extensions to this
convention. These extension are:
* A variable "cell_origin" (of dimension "frame", "cell_spatial") that contains
  the bottom left corner of the simulation cell.
* A variable "atom_types" (of dimension "atom") that contains the LAMMPS atom
  type.
* Any number of additional variables corresponding to per atom scalar, vector
  or tensor quantities available within LAMMPS. Tensor quantities are written in
  Voigt notation. An additional dimension "Voigt" of length 6 is created for
  this purpose.

NetCDF files can be directly visualized with the following tools:
* Ovito (http://www.ovito.org/). Ovito supports the AMBER convention and all of
  the above extensions.
* VMD (http://www.ks.uiuc.edu/Research/vmd/).
* AtomEye (http://www.libatoms.org/). The libAtoms version of AtomEye contains
  a NetCDF reader that is not present in the standard distribution of AtomEye.

Syntax:

> dump ID group-ID nc N file args  
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

The list of atom attributes is identical to the 'custom' dump style.

Example:

> dump 1 all nc 100 traj.nc type x y z vx vy vz  
> dump_modify 1 append yes at -1 global c_thermo_pe c_thermo_temp c_thermo_press

INSTALLATION
------------

In your LAMMPS src directory type:

> git clone https://github.com/pastewka/lammps-netcdf.git USER-DUMP-NC  
> make yes-user-dump-nc

Note that LAMMPS will need to be linked to NetCDF. This will require a
modification of your favorite makefile. Please add

> EXTRA_INC += $(shell nc-config --cflags)  
> EXTRA_LIB += $(shell nc-config --libs)

to the respective EXTRA_INC, EXTRA_LIB section of the makefile.

OTHER NOTES
-----------

This is package is known to work with LAMMPS 13May14.
