#!/usr/bin/env bash

export CC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
export FC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gfortran

# debug flags:
export CFLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -pipe -I${CONDA_PREFIX}/include"
export FFLAGS="-fopenmp -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -I${CONDA_PREFIX}/include -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -pipe"

export THIS_DIR="`pwd`/`dirname "$BASH_SOURCE"`"
export PROJECT_ROOT="`pwd`/`dirname "$BASH_SOURCE"`"/..

mkdir -p $THIS_DIR/build
cd $THIS_DIR/build

cmake $PROJECT_ROOT -DCMAKE_INSTALL_PREFIX=$THIS_DIR -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_FLAGS="${FFLAGS}" -DCMAKE_C_FLAGS="${CFLAGS}"
make VERBOSE=1
make install
