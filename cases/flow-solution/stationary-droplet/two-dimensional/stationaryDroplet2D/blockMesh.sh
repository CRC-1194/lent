#!/usr/bin/bash

# This is a work-around for the include directive in some dictionaries.
# Include only works correctly when the solver is executed inside the case.
# TODO devise an actual fix for this. Unfortunately, pyFoam seems to overwrite
# dictionaries after the caseSetup script, so every means there has no effect
mkdir ../constant
cp ./constant/caseSetup ../constant/

blockMesh
