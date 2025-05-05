#!/bin/sh
cmake -B build -D CMAKE_INSTALL_PREFIX=/usr/local/Xyce -D Trilinos_ROOT=/usr/local/XyceLibs/Parallel -G Ninja


