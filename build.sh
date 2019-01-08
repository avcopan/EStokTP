#!/usr/bin/env bash

mkdir build
cd build

if [ "$(uname)" == "Darwin" ]; then
    echo "ERROR: Building with Mac OSX is currently unsupported" 1>&2
    exit 64
elif [ "$(uname)" == "Linux" ]; then
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_C_COMPILER=$CC -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_Fortran_FLAGS="${DEBUG_FFLAGS}" -DCMAKE_C_FLAGS="${DEBUG_CFLAGS}"
fi

make VERBOSE=1
make install
cp $RECIPE_DIR/exe/assemble_var.com        $PREFIX/bin/.
cp $RECIPE_DIR/exe/get_files.com           $PREFIX/bin/.
cp $RECIPE_DIR/exe/run_estoktpb.com        $PREFIX/bin/.
cp $RECIPE_DIR/exe/run_estoktp.com         $PREFIX/bin/.
cp $RECIPE_DIR/exe/run_estoktp_example.com $PREFIX/bin/.
cp $RECIPE_DIR/exe/extract_data_projrot.sh $PREFIX/bin/.
