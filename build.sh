#/bin/bash

MOAB_PATH=$HOME/opt/moab5.0/
#MOAB_PATH=$HOME/research/moab
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"

echo $MOAB_LIBRARY

g++  -std=c++11 mesh2mesh.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -lboost_program_options -o test2_m2m

