#!/bin/bash

module purge
module load cmake gcc openmpi intel-mkl # REPLACE THIS WITH YOUR EQUIVALENT MODULE CALLS

# rm -r build
# mkdir build

if [[ ! -d "build" ]]; then
  mkdir build
fi

cd build

cmake .. -G "Unix Makefiles" -DCMAKE_C_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=build -DCMAKE_CXX_COMPILER=mpic++

make

cd ..

#chgrp -R wgn1_collab build
#chmod -R 770 build
