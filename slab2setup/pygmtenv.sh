echo "Creating the pygmt virtual environment:"
conda activate base
conda create --name pygmt --channel conda-forge pygmt
conda activate pygmt
conda install -c conda-forge fiona
conda install -c conda-forge shapely
conda deactivate
echo 'pygmt virtual environment has been created'