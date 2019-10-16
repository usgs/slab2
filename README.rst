Introduction
----------------------------------------
Slab2.0 is a three-dimensional (3D) compilation of global subduction zone geometries, separated into regional models for each major subduction zone. Refer to Hayes et al., 2018, Science for more details. This code contains functions for making 3D slab geometries and for updating the input for making such geometries. Please note that this code is still under development. Send all inquiries to Dr. Gavin P. Hayes <ghayes@usgs.gov>.

Installation and Dependencies
----------------------------------------
A slab2 environment (slab2env) will be created during installation. This will install all the packages and libraries that slab2 code depends on. You must have the slab2env activated in order to run any slab2 code. The install has been tested only on MacOSX, and not on Windows or Linux.

This package depends on:
   * Mac OSX operating system
   * Python3
   * Anaconda
   * pip
   * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/</a>
   * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="http://matplotlib.org/index.html</a>
   * scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. <a href="http://www.scipy.org/scipylib/index.html</a>
   * mapio, a Python library for reading/writing various spatial data formats (GMT grids, etc.)
   * libcomcat, a Python library for retrieving earthquake data from the ANSS ComCat system
   * GMT5, an open source collection of tools for manipulating geographic and Cartesian data sets <a href=”http://gmt.soest.hawaii.edu/projects/gmt/wiki/Documentation</a> 

Anaconda was used by Slab2 programmers. Most of those distributions should include pip, a command line tool for installing and managing Python packages.

You may need to open a new terminal window to ensure that the newly installed versions of python and pip are in your path. Prior to running slab2 code, make sure to add GMT5 to your path.

To install slab2:

Clone the slab2 repository, and then follow instructions in Slab2Instructions.pdf under Installation.