**This Repository Has Moved**

This repository is outdated. Instead please use https://code.usgs.gov/ghsc/esi/usgs-slab-models. 


Introduction
----------------------------------------
Slab2.0 is a three-dimensional (3D) compilation of global subduction zone geometries, separated into regional models for each major subduction zone. Refer to Hayes et al., 2018, Science for more details. This code contains functions for making 3D slab geometries and for updating the input for making such geometries. Please note that this code is still under development. To bypass installing and running the code in this repository, please visit the USGS Slab Model Cloud-based web application at https://earthquake.usgs.gov/slab2/ to generate a Slab2 model. Send all inquiries to Dr. Gavin P. Hayes <ghayes@usgs.gov> and Dr. Kirstie L. Haynie <khaynie@usgs.gov>.

Installation and Dependencies
----------------------------------------
A Slab2 environment (slab2env) will be created during installation. This will install all the packages and libraries that Slab2 code depends on. You must have the slab2env activated in order to run any Slab2 code that creates the models. The install has been tested only on MacOSX, and not on Windows or Linux.
In addition to the slab2env, an environment used for plotting will also be created (pygmt). To run any of the code that creates plots for the Slab2 models, the pygmt environment must be activated.


This package depends on:
   * Python3
   * Anaconda
   * pip
   * `numpy <http://www.numpy.org/>`_, the fundamental package for scientific computing with Python.
   * `matplotlib <http://matplotlib.org>`_, a Python 2D plotting library which produces publication quality figures. 
   * `scipy <https://scipy.org>`_, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. 
   * mapio, a Python library for reading/writing various spatial data formats (GMT grids, etc.)
   * libcomcat, a Python library for retrieving earthquake data from the ANSS ComCat system
   * `pygmt <https://www.pygmt.org/latest/>`_, a Python library which wraps GMT, allowing for easy manipulation of geographic and cartesian datasets.

Follow the steps below to run the Slab2 code.

1. **Clone the Slab2 Repository:**

   Use the following, but provide the appropriate link:

    ``git clone link/to/repository``

2. **Create the anaconda environments**

   In the terminal within the slab2/slab2code/slab2setup directory, enter:
   
    ``bash slab2env.sh both``

   Two anaconda environments, slab2env and pygmt, will be created. The slab2env should be active when using code that creates the Slab2 output data, while the pygmt environment should be active when plotting this data.

   To only create the Slab2 environment, use:

   ``bash slab2env.sh slab``

   To only create the PyGMT environment, use:

   ``bash slab2env.sh gmt``

3. **Create an input file**

   The slab2env environment must be active, and the active directory should be slab2/slab2code

   In the terminal, enter: 

    ``conda activate slab2env``

    ``python sd2.py -p [slab] -d [Slab2_location]/[input date]database -f [input file csv]``


4. **To create a Slab2 model**

   The slab2env environment must be active, and the active directory should be slab2/slab2code

   In the terminal, enter: 

    ``conda activate slab2env``

    ``python slab2.py -p [path to parameter file]``

   This will create the model with the  most recent input file

   By appending "-c 3" to the end of the command, three cores will be used.
   
   To change input parameters, please refer to number 5 of Creating a Slab Model within `Slab2Instructions.pdf <./Slab2Instructions.pdf>`_.

   Regions *mue, sul, phi,* and *ryu* depend on depths from other slabs. Please see number 3 of Creating a Slab Model within `Slab2Instructions.pdf <./Slab2Instructions.pdf>`_ for further details.

5. **To plot the Slab2 modeling results**

   Navigate to slab2/slab2code/plotting.

   The pygmt environment must be active when using code in this directory

   To make an overview map, use map.py

   In the terminal, enter:

    ``conda activate pygmt``

    ``python map.py [slab] [output date]``

   To make cross sections, use xsec.py

   In the terminal, enter: 

    ``conda activate pygmt``

    ``python map.py [slab] [output date] [input date] all``


For more complete examples and instructions, please see `Slab2Instructions.pdf <./Slab2Instructions.pdf>`_.

