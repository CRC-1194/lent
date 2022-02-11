#!/bin/bash

for study in *Translation.parameter; 
do 
    agglomerate_pyfoam_studydata_notebook.py "$study"_00000_template_copy_block/stationaryDropletResults.csv -p $study; 
done

for notebook in *Translation.ipynb; 
do
    echo "jupyter-nbconvert --execute --to=html $notebook"
    jupyter-nbconvert --execute --to=html $notebook 
done

agglomerate_pyfoam_studydata_notebook.py \
	without-rhoEquation-densityRatioInfluence.parameter_00000_template_copy_block/stationaryDropletResults.csv \
	-p densityRatioInfluence.parameter

jupyter-nbconvert --execute --to=html without-rhoEquation-densityRatioInfluence.ipynb

agglomerate_pyfoam_studydata_notebook.py \
	with-rhoEquation-densityRatioInfluence.parameter_00000_template_copy_block/stationaryDropletResults.csv \
	-p densityRatioInfluence.parameter

jupyter-nbconvert --execute --to=html with-rhoEquation-densityRatioInfluence.ipynb

agglomerate_pyfoam_studydata_notebook.py \
	with-rhoEquation-popinet2009.parameter_00000_template_copy_block/stationaryDropletResults.csv \
	-p popinet2009.parameter

jupyter-nbconvert --execute --to=html with-rhoEquation-popinet2009.ipynb

agglomerate_pyfoam_studydata_notebook.py \
	without-rhoEquation-popinet2009.parameter_00000_template_copy_block/stationaryDropletResults.csv \
	-p popinet2009.parameter

jupyter-nbconvert --execute --to=html without-rhoEquation-popinet2009.ipynb
