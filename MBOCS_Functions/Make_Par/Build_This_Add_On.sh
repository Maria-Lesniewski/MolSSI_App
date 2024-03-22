#!/bin/bash

#module purge
module load cmake openmpi mkl # CHANGE FOR YOUR MODULE CALLS

# rm -r build
# mkdir build

if [[ ! -d "build" ]]; then
  mkdir build
fi

cd build

cmake ../ -G "Unix Makefiles" -DCMAKE_C_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=build -DCMAKE_CXX_COMPILER=mpic++

make

cd ../
