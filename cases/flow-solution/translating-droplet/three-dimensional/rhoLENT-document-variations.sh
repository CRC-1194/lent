#!/bin/bash

for param_file in *.parameter; 
do 
    echo Storing variations for: $param_file
    pyFoamRunParameterVariation.py \
	--list-variations templateCase_translatingDroplet3D \
	$param_file > ${param_file%%.parameter}.variations
done

