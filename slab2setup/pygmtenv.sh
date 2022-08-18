echo "Creating the pygmt virtual environment:"
conda create --name pygmt --channel conda-forge pygmt
conda env update --file pygmt.yml
echo 'pygmt virtual environment has been created'
