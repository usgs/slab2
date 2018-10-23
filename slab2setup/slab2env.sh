#!/bin/bash

VENV=slab2env
ENVFILE=slab2env.yml

echo "Environment file: $ENVFILE"
echo "Creating the $VENV virtual environment:"

source deactivate
conda env create -f $ENVFILE --force

source activate $VENV

pip install git+git://github.com/usgs/earthquake-impact-utils.git

source deactivate
