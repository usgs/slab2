slab2
=======

Slab2.0 is a three-dimensional compilation of global subduction geometries, separated into regional models for each major subduction zone.

Installation and Dependencies
-----------------------------

This package (and it's notebooks) depend on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. <a href="<a href="http://www.scipy.org/scipylib/index.html">http://www.scipy.org/scipylib/index.html</a>
 * neicmap, a Python library for doing various spatial calculations (distance, angle, etc.)
 * neicio, a Python library for reading/writing various spatial data formats (GMT grids, etc.)
 * libcomcat, a Python library for retrieving earthquake data from the ANSS ComCat system

The best way to install numpy,matplotlib,and scipy is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda and Enthought distributions have been successfully tested with smtools.

Most of those distributions should include <em>pip</em>, a command line tool for installing and 
managing Python packages.  You will use pip to install the other dependencies and smtools itself.  
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install neicmap:

pip install git+git://github.com/usgs/neicmap.git

To install neicio:

pip install git+git://github.com/usgs/neicio.git

To install libcomcat:

pip install git+git://github.com/usgs/libcomcat.git

To install slab2:

pip install git+git://github.com/usgs/slab2.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall slab2

To update:

pip install -U git+git://github.com/usgs/slab2.git

Examples
---------
A notebook exploring some of the functionality of geostatsmodels is included in this repository.  

http://nbviewer.ipython.org/github/usgs/slab2/blob/master/notebooks/IntroToSlab2.0.ipynb

More of these will follow.