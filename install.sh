#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

echo "***************************************************************************"
echo WARNING: lent is tested with gcc version 10.3.0-2 or *older*
echo Change the example lines 
echo 'export CC=`which gcc-10`'
echo 'export CXX=`which g++-10`'
echo in install.sh to point CMake to your compiler.
echo "***************************************************************************"


export CC=`which gcc-10`
export CXX=`which g++-10`
rm -rf build && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
make -j4 && make install
