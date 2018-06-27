#!/bin/bash

VENV=slab2env3
ENVFILE=slab2env3.yml

echo "Environment file: $ENVFILE"
echo "Creating the $VENV virtual environment:"
source deactivate
conda env create -f $ENVFILE --force

source activate $VENV

pip install git+git://github.com/usgs/MapIO.git
pip install git+git://github.com/usgs/neicmap.git
pip install multiprocess
pip install sklearn
pip install utm
pip install geopy
pip install obspy
pip install git+git://github.com/usgs/earthquake-impact-utils.git
pip install git+git://github.com/usgs/libcomcat.git
conda install -c conda-forge basemap-data-hires

source deactivate
