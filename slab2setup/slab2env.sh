#!/bin/bash

VENV=slab2env
ENVFILE=slab2env.yml

echo "Environment file: $ENVFILE"
echo "Creating the $VENV virtual environment:"

source ~/.bash_profile

. $_CONDA_ROOT/etc/profile.d/conda.sh

conda deactivate
conda env create -f $ENVFILE --force

conda activate $VENV

pip install git+git://github.com/usgs/MapIO.git@0.6.2
pip install git+git://github.com/usgs/libcomcat.git@1.0
pip install git+git://github.com/usgs/earthquake-impact-utils.git@0.8.2

conda deactivate
