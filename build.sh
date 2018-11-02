#!/usr/bin/env bash

mkdir build
cd build

if [ "$(uname)" == "Darwin" ]; then
    echo "ERROR: Building with Mac OSX is currently unsupported" 1>&2
    exit 64
elif [ "$(uname)" == "Linux" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$GFORTRAN -DBLAS_LIBRARIES=$PREFIX/lib/libopenblas.so
fi

make
make install
