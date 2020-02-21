#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

rm -rf build && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=on ..
make && make install
