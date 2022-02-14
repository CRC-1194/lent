#!/bin/bash

if [ -d rhoLENT-data ]; 
then 
    echo "rhoLENT-data folder exits, exiting." 
    exit 1
fi

mkdir rhoLENT-data
cp rhoLENT-data.md rhoLENT-data
cp --parents $(find . -name stationaryDropletResults.csv) rhoLENT-data/
cp ../../../scripts/modules/*.py rhoLENT-data/
./rhoLENT-document-variations.sh
cp *.parameter *.variations rhoLENT-data/
./rhoLENT-visualize-studies.sh
cp *.ipynb *.html *.csv *.json rhoLENT-data/
tar czf rhoLENT-data.tgz rhoLENT-data
