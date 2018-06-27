#!/usr/bin/env python

### Created DEP.5.26.16 ###
### Edited DEP.6.3.16 to include shallim ###
### Edited DEP.7.26.16 to change error calculation ###

## tomo_slab stands for Tomography Slab2.0.  I intend for this to be the master script that will query a given tomography model for the slab center points.  These data will be input into the Slab2.0 code that will then migrate the slab based on the calibration function.

## Necessary inputs will be 1) Tomography file in a four column text file; 2) slab1.0 grid file for depth.

## Output consists of one output file which is an ascii file of with columns of the form 'lon lat depth uncertainty strike dip'

# 1) [slab]_TO_[model].csv

###########################

### Section 0: Import libraries

import neicio.gmt as gmt  ## Grid functions Mike created
import numpy as np  ## Allows for work on arrays like matlab
import slabfunctions as sf  ## Functions that I wrote for this code
import math  ## I assume this allows for simple math
from scipy.interpolate import griddata  ## Allows me to interpolate grids for plotting
import matplotlib.pyplot as plt  ## Plotting tools!
from scipy import stats  ## In order to do bin statistics for error analysis on function
import csv

############################

### Section 1: Define inputs

## Input files

tomofile = 'Tomo_files/ust_RWP14.tomo'  ## Tomography model in ascii format 'Lat Lon Dep(positive) dVp'
slabfile = 'Slab_files/cas_slab_guide_DEP.7.29.16.grd' ## Slab guide file in netcdf format

## Output file

ID = 'RWP14'  # Identifier for tomography model to be used in output files
sz = 'cas'  # Code for the subduction zone associated with this tomography model

outfile = '%s_TO_%s.csv' % (sz,ID)  ## DON'T CHANGE  
outfile2 = 'par_%s_%s.out' % (ID,sz)  ## Paramter file

# Calibration region boundaries

minlat, maxlat = 39, 50
minlon, maxlon = -123, -119

# Sample volume limits

posmag, negmag = 150, -150  # km in positive and negative directions, respectively, from slab input ***

# Node spacing for slab

hnode, vnode = 0.2, 1   # horizontal node spacing (degrees) for grid and vertical (km) sampling along vector (0.2 for regional, 0.5 for global)

# Node spacing for upsampled tomography

HRES, VRES = 0.05, 5.0  ## H in degrees, V in km

# Anomaly value threshold

thresh = 0.1 # % velocity perturbation ***

# Shallow limit for looking for slab

shallim = 100  # Depth (km) that will apply to the slab guide

# Minimum mumber of slab sample points

minnum = 20
maxnum = 200

# depth threhsolds

mindepthcut = 60
maxdepthcut = 500
###########################

### Section 2: Load input files

print "Loading input files..."

Tomo = np.loadtxt(tomofile) ## Load tomography file into a numpy array
Tomo[:,1][Tomo[:,1]>180] -= 360
Slabgrid = gmt.GMTGrid(slabfile) ## Load netcdf file of slabguide into a "grid"
if Slabgrid.geodict['xmin'] > 180:
	Slabgrid.geodict['xmin'] -= 360
if Slabgrid.geodict['xmax'] > 180:
	Slabgrid.geodict['xmax'] -= 360

Slabgrid.griddata = Slabgrid.griddata * -1  ## WILL NEED TO CHANGE THIS IF INPUT IS IN NEGATIVE DEPTH

print "Building strike and dip grids..."

Strikegrid, Dipgrid = sf.mkSDgrd(Slabgrid)  ## Strikes follow right hand rule, dips are positive in degrees

############################

### Section 3: Create grid of slab guide info

## Create grid structure

Lons = np.arange(minlon,maxlon,hnode)
Lats = np.arange(minlat,maxlat,hnode)
m = len(Lons)
n = len(Lats)

Lons, Lats = np.meshgrid(Lons,Lats)
XYgrid = np.zeros((m*n,5))  ## 5 columns for 'Lon Lat Dep Dip Str'
XYgrid[:,0] = Lons.flatten()
XYgrid[:,1] = Lats.flatten()

## Populate grid with slab guide depth, dip, and strike

for i in range(0,len(XYgrid)):
    lon = XYgrid[i,0]
    lat = XYgrid[i,1]

    XYgrid[i,2] = Slabgrid.getValue(lat,lon)  ## Populate slab depth (note depths should be positive)
    XYgrid[i,3] = Dipgrid.getValue(lat,lon)  ## Populate slab dip (note dips should be positive)
    XYgrid[i,4] = Strikegrid.getValue(lat,lon) ## Populate slab strike

############################

### Section 4: Make tomography cube array

print "Buliding tomography cube..."
Cube, ulx, uly, ulz, xdim, ydim, zdim = sf.makeCube(Tomo,HRES,VRES)  ## Lon (-180 to 180), Lat (-90 to 90), Dep pos

############################

### Section 5: Find the slab in the tomography

print "Finding the slab center..."

SlabCent = np.zeros((len(XYgrid),3))  ## This array is the size of the slab guide grid...
    #meant to have one slab center point for each slab guide node

for i in range(0,(len(XYgrid))):  ## Iterate through every point in XYgrid (each node)
    if np.isnan(XYgrid[i,4]) == False:  ## First check to make sure there is a value for the slab guide
        
        if XYgrid[i,2] > shallim:  ## Make sure the slab guide is below the threshold for what you want to look at
	
        # Calculate the locations of the points to sample at that node
        
            EPs = sf.epCalc(XYgrid[i,0],XYgrid[i,1],XYgrid[i,2],XYgrid[i,3],XYgrid[i,4],posmag,negmag,vnode)
                #EPs is an array of locations to sample, with the fourth column saying how far from the center node each point is
        
            # Sample tomography at each sampling point from EPs
        
            NodeVec = np.zeros((len(EPs),5))  ## Define a new array that will contain the location and tomography anomaly value
            NodeVec[:,0:4] = EPs  ## Fill the first four columns identically to EPs to note location
        
            for j in range(0,len(EPs)):  ## Iterate through each sampling point for the present node
            
                row, col, band = sf.getRowColBand(Cube,EPs[j,0],EPs[j,1],EPs[j,2],ulx,uly,ulz,xdim,ydim,zdim)
                    #Determine the position in the tomography cube given coordinates
            
                NodeVec[j,4] = sf.get3dvalue(Cube,row,col,band)  ## Sample the tomography model at that point, fill into NodeVec

                if np.isnan(NodeVec[j,4]) == True:
                    NodeVec[j,:] = np.nan
                if NodeVec[j,4] < thresh:
                    NodeVec[j,:] = np.nan  ## Make whole lines NaN when the tomography value is below a threshold...
                        #This hopefully retains only "slab"

	    num = len(NodeVec[np.isfinite(NodeVec[:,4])])
	    if num < minnum:
		avedep = np.nan
		avelon = np.nan
		avelat = np.nan
	    elif num > maxnum:
		avedep = np.nan
		avelon = np.nan
		avelat = np.nan
	    else:
            	avedis = np.around(np.nanmean(NodeVec[:,3]))  ## Computes the average "distance from guide" of "slab-only" points
            	ind = NodeVec[:,3] == avedis  ## Indicates the index for the location that is at the average distance
            	avedep = NodeVec[ind,2]  ## Determine the depth at the "average location"
            	avelon = NodeVec[ind,0]  ## Determine the longitude at the "average location"
            	avelat = NodeVec[ind,1]  ## Determine the latitude at the "average location"
 
            if np.isnan(avedep) == False:  ## Ignores when the avedep is NaN so that the line is left as zero.  This is taken care of in the following lines.  Currently might be inefficient
                SlabCent[i,0:3] = avelat,avelon,avedep  ## Populate error function array ***using avelon currently for testing/plotting purposes
 
tmp = SlabCent[:,0]  ## These lines are used to index which avedeps are 0 and make the whole line NaN
SlabCent[tmp == 0,:] = np.nan

output = SlabCent[np.isfinite(SlabCent[:,2])]
output = output[output[:,2] < maxdepthcut]
output = output[output[:,2] > mindepthcut]

############################

### Section 6: Print results

out = open(outfile,'w')
out.write("lat,lon,depth\n")
np.savetxt(out,output,fmt='%0.4f',delimiter=',')
out.close()

params = open(outfile2,'w')
params.write('Parameter file for tomography ID: %s\n\n' % ID)
params.write('Output datafile: %s\n\n' % outfile)
params.write('Input tomography file: tomofile = %s\n' % tomofile)
params.write('Slab guide grid file: slabfile = %s\n\n' % slabfile)
params.write('minlat = %f\n' % minlat)
params.write('maxlat = %f\n' % maxlat)
params.write('minlon = %f\n' % minlon)
params.write('maxlon = %f\n' % maxlon)
params.write('posmag = %f\n' % posmag)
params.write('negmag = %f\n' % negmag)
params.write('hnode = %f\n' % hnode)
params.write('vnode = %f\n' % vnode)
params.write('HRES = %f\n' % HRES)
params.write('VRES = %f\n' % VRES)
params.write('thresh = %f\n' % thresh)
params.write('shallim = %f\n' % shallim)
params.write('minnum = %f\n' % minnum)
params.write('maxnum = %f\n' % maxnum)
params.write('mindepthcut = %f\n' % mindepthcut)
params.write('maxdepthcut = %f\n' % maxdepthcut)
params.close()
