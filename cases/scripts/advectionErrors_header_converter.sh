#!/usr/bin/env bash

# Use this script to replace the headers of result files
# which have a header not matching the CSV format

STUDY=$1

for case in "$STUDY"_*;
do
    cd $case

    # This will replace the advection error header
    sed -i "1s/.*/time,volume_conservation_error,L1_advection_error,L1_advection_error_norm/" advectionErrors.dat
    # This will replace spaces with commas. 'g' ensures that all occurances in a line are replaced, not only the first one
    sed -i "s/ /,/g" advectionErrors.dat 

    cd ..
done
