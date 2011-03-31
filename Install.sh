# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  cp -p dump_nc.cpp ..

  cp -p dump_nc.h ..

elif (test $1 = 0) then

  rm ../dump_nc.cpp

  rm ../dump_nc.h

fi
