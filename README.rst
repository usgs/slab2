slab2
========
Slab2.0 is a three-dimensional compilation of global subduction geometries, separated into regional models for each major subduction zone.
Note, this code is still under development. Please contact Gavin Hayes for more details. 

Installation and Dependencies
--------
This package depend on:
•	Mac OSX operating system
•	Python 3
•	numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/
•	matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="http://matplotlib.org/index.html
•	scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. <a href="http://www.scipy.org/scipylib/index.html
•	mapio, a Python library for reading/writing various spatial data formats (GMT grids, etc.)
•	libcomcat, a Python library for retrieving earthquake data from the ANSS ComCat system
•	GMT5, an open source collection of tools for manipulating geographic and Cartesian data sets <a href=” http://gmt.soest.hawaii.edu/projects/gmt/wiki/Documentation 

The best way to install numpy, matplotlib,and scipy is to use one of the Python distributions described here:

http://www.scipy.org/install.html

Anaconda and Enthought distributions have been successfully tested with smtools. Anaconda was used by Slab2 programmers. 

Most of those distributions should include pip, a command line tool for installing and managing Python packages. You will use pip to install the other dependencies and smtools itself.

You may need to open a new terminal window to ensure that the newly installed versions of python and pip are in your path.

To install mapio:

See https://github.com/usgs/MapIO

To install libcomcat:

See https://github.com/usgs/libcomcat

To install slab2:

Clone the slab2 repository, and then follow instructions in the Slab2Instructions doc. 
