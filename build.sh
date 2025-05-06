#!/bin/sh
cmake -B build -D CMAKE_BUILD_TYPE=Debug -DXyce_ADMS_MODELS=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -D CMAKE_INSTALL_PREFIX=~/local/Xyce -D Trilinos_ROOT=/usr/local/XyceLibs/Parallel -G Ninja
