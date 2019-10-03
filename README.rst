Introduction
----------------------------------------
Slab2.0 is a three-dimensional (3D) compilation of global subduction zone geometries, separated into regional models for each major subduction zone. Refer to Hayes et al., 2018, Science for more details. This code contains functions for making 3D slab geometries and for updating the input for making such geometries. Please note that this code is still under development. Send all inquiries to Dr. Gavin P. Hayes <ghayes@usgs.gov>.

Installation and Dependencies
----------------------------------------
A slab2 environment (slab2env) will be created during installation. This will install all the packages and libraries that slab2 code depends on. You must have the slab2env activated in order to run any slab2 code. The install has been tested only on MacOSX, and not on Windows or Linux.

This package depends on:
   * Xcode and the Xcode command line tools
   * Anaconda
   * pip
   * GMT5, an open source collection of tools for manipulating geographic and Cartesian data sets <a href=”http://gmt.soest.hawaii.edu/projects/gmt/wiki/Documentation</a> 

Anaconda was used by Slab2 programmers. Most of those distributions should include pip, a command line tool for installing and managing Python packages.

You may need to open a new terminal window to ensure that the newly installed versions of python and pip are in your path. Prior to running slab2 code, make sure to add GMT5 to your path.

To install slab2:

Clone the slab2 repository, and then follow instructions in Slab2Instructions.pdf under Installation.