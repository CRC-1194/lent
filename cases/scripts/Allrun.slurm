#! /usr/bin/bash

STUDY=$1
MESHTYPE=$2

# -d checks if a directory of the given name exists
if [ -z "$STUDY" ] || [ -z "$MESHTYPE" ]; 
then 
    echo Error in execution. Execution: ./Allrun.slurm studyName meshType
    echo studyName: name of the parameter file to be used
    echo meshType: meshing approach to be used. Valid options are
    echo block, cartesian, poly and hexrefined
    exit 1
fi

# Check the given mesh type
if [[ "$MESHTYPE" != "block" ]] && [[ "$MESHTYPE" != "cartesian" ]] && [[ "$MESHTYPE" != "poly" ]] && [[ "$MESHTYPE" != "hexrefined" ]];
then
    echo Error: mesh-type must be block, cartesian, poly or hexrefined.
    exit 1
fi

# TODO: solve this in a more robust/elegant way
# Set relative path to sbatch scripts
SCRIPTPATH=../../../../scripts

# Setup variation folders
lent_prepare_study_variants.py -s $STUDY -m $MESHTYPE -np

# Dispatch preprocessing steps (mesh creation, front creation, field initialization)
# as dependent jobs to the compute nodes.
# The final job is the simulation itself of course.
for case in "$STUDY"_0*"$MESHTYPE";
do 
    cd $case

    echo sbatch --parsable "$SCRIPTPATH"/00-"$MESHTYPE"Mesh.sbatch 2>&1
    MESH=$(sbatch --parsable "$SCRIPTPATH"/00-"$MESHTYPE"Mesh.sbatch 2>&1) 

    echo sbatch --parsable --dependency=afterok:"$MESH" "$SCRIPTPATH"/01-lentCreateFront.sbatch  2>&1 
    CREATEFRONT=$(sbatch --parsable --dependency=afterok:"$MESH" "$SCRIPTPATH"/01-lentCreateFront.sbatch 2>&1) 

    echo sbatch --parsable --dependency=afterok:"$CREATEFRONT" "$SCRIPTPATH"/02-lentSetFields.sbatch 2>&1
    SETFIELDS=$(sbatch --parsable --dependency=afterok:"$CREATEFRONT" "$SCRIPTPATH"/02-lentSetFields.sbatch 2>&1)

    echo sbatch --parsable --dependency=afterok:"$SETFIELDS" "$SCRIPTPATH"/03-lentFoam.sbatch 2>&1
    LENTFOAM=$(sbatch --parsable --dependency=afterok:"$SETFIELDS" "$SCRIPTPATH"/03-lentFoam.sbatch 2>&1)

    echo $case

    cd ..
done
