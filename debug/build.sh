#!/usr/bin/env bash

export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export CC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
export FC=$CONDA_PREFIX/bin/x86_64-conda_cos6-linux-gnu-gfortran

# debug flags:
export CFLAGS="-march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-all -fno-plt -Og -g -Wall -Wextra -fvar-tracking-assignments -ffunction-sections -pipe -I${CONDA_PREFIX}/include"
export FFLAGS="-fopenmp -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -I${CONDA_PREFIX}/include -fstack-protector-all -Og -g -Wall -Wextra -fcheck=all -fbacktrace -fimplicit-none -fvar-tracking-assignments -ffunction-sections"

mkdir -p $THIS_DIR/build
cd $THIS_DIR/build

cmake $PROJECT_ROOT -DCMAKE_INSTALL_PREFIX=$THIS_DIR -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_FLAGS="${FFLAGS}" -DCMAKE_C_FLAGS="${CFLAGS}"
make VERBOSE=1
make install
