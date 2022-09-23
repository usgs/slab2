#!/bin/bash

VENV=slab2env
ENVFILE=slab2env.yml

if [ $1 == slab ]
then
    echo "Environment file: $ENVFILE"
    echo "Creating the $VENV virtual environment:"
    conda env create --name $VENV --file=$ENVFILE
fi

if [ $1 == gmt ]
then
    bash pygmtenv.sh
fi

if [ $1 == both ]
then
    echo "Environment file: $ENVFILE"
    echo "Creating the $VENV virtual environment:"
    conda env create --name $VENV --file=$ENVFILE
    conda base
    bash pygmtenv.sh
fi
