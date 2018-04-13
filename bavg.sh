#/bin/bash

MOAB_PATH=$HOME/opt/moab5.0/
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"

echo $MOAB_LIBRARY

g++  -std=c++11 avg.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -o avg

