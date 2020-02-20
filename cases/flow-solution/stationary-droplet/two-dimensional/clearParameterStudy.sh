#!/usr/bin/bash

# Clean up script for parameter studies. Takes one argument:
#   all:            clean EVERYTHING, meaning removal of directories and
#                   *.dat result files
#   dirs:           remove all simulation directories, but keep *.dat result files
#   "studyName":    removes directories beloging to the given parameter file
#                   but keeps *.dat result files

# Clear pyFoam files
rm -f PlyParser_FoamFileParser_parsetab.py *.database
rm -rf constant

# Clear all studies and result file
if [ "$1" == "all" ]; then
    foo=$(ls | grep "_[0-9]\{5\}_")
    rm -rf $foo *.dat
# Clear all simulation directories, but keep *.dat files
elif [ "$1" == "dirs" ]; then
    foo=$(ls | grep "_[0-9]\{5\}_")
    rm -rf $foo
# Remove directories of given parameter study, but keep *.dat files
elif [ "$1" != "" ]; then
    foo=$(ls | grep "$1" | grep "_[0-9]\{5\}_")
    rm -rf $foo
fi
