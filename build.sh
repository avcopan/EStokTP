#!/usr/bin/env bash

mkdir build
cd build

if [ "$(uname)" == "Darwin" ]; then
    echo "ERROR: Building with Mac OSX is currently unsupported" 1>&2
    exit 64
elif [ "$(uname)" == "Linux" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FC -DBLAS_LIBRARIES=$PREFIX/lib/libblas.so -DCMAKE_Fortran_FLAGS="${FFLAGS}" -DCMAKE_C_FLAGS="${CFLAGS}"
fi

make VERBOSE=1
make install
cp $RECIPE_DIR/exe/*.com $PREFIX/bin/.
