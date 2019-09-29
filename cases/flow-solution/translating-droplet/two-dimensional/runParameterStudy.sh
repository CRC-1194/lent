#!/usr/bin/bash

# Description:
# Performs the setup and the execution of a parameter study for a given
# template case and a parameter file.
# Can also be used to set up a single run capable variant of the
# given parameter study
#
# Usage:
# ./runParameterStudy arg1 arg2 (arg3)
#   arg1:   either 'poly' or 'block', specifying the mesh type
#   arg2:   name of a parameter file
#   arg3:   optional: number of the variation. 
#           Use "pyFoamRunParameterVariation.py --list-variations test_case parameter_file"
#           to list the different variations and their parameter set

case="translatingDroplet2D"

# Check that mesh type is valid
if [[ "$1" != "poly" && "$1" != "block" ]]; then
    echo "Error: $1 is not a valid mesh type. Choose block or poly."
    exit 1
fi

./setDomainAndInterfaceGeometry.py $case $2

cd $case
./makeFront

rm -rf 0 0.org
cp -r "0.org_${1}Mesh" 0.org
cp -r "0.org_${1}Mesh" 0

cd ..

# Folowing expression tests if arg3 is empty
if [ -z "$3" ]; then
    pyFoamRunParameterVariation.py --every-variant-one-case-execution --create-database --no-server-process --mesh-create-script="${1}Mesh.sh" $case $2
else
    pyFoamRunParameterVariation.py --every-variant-one-case-execution --create-database --no-server-process --no-execute-solver --mesh-create-script="${1}Mesh.sh" --single-variation=$3 $case $2
fi
