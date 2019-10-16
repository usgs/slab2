#!/bin/bash

VENV=slab2env
ENVFILE=slab2env.yml

echo "Environment file: $ENVFILE"
echo "Creating the $VENV virtual environment:"

conda env create --name $VENV --file=$ENVFILE
