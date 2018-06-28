#!/usr/bin/env python

### Functions required for tomo_slab.py ###

###############################################

### Module imports

import numpy as np
import math
import mapio.gmt as gmt
from scipy import interpolate
from scipy.interpolate import griddata
import os
from matplotlib import path
from scipy import ndimage

###

def epCalc(lon, lat, dep, dip, strike, posmag, negmag, step):
    tmplength = int((posmag-negmag)/step)
    EPs = np.zeros((tmplength, 4))
    
    # Rotate from strike to direction of motion
    
    if strike > 270:
        az = (strike + 90) - 360
    else:
        az = strike + 90
    az = 360 - az    # Accounts for the fact that azimuth counts goes opposite of the positive rotation of the x-axis (points north)
    
    # Convert input angles to radians
    
    latrad = math.radians(90 - lat)
    lonrad = math.radians(lon)
    azrad = math.radians(az)
    diprad = math.radians(dip)
    
    # Define initial location in spherical coordinates
    
    crad = 6371 - dep
    ctheta = latrad
    cphi = lonrad
    
    # Convert initial location to cartesian coordinates
    
    cx = crad * math.sin(ctheta) * math.cos(cphi)
    cy = crad * math.sin(ctheta) * math.sin(cphi)
    cz = crad * math.cos(ctheta)
    
    # Define lon/lat of new coordinate system
    
    if latrad < (math.pi/2):
        x1lat = abs(latrad-(math.pi/2))
        if lonrad > 0:
            x1lon = lonrad - math.pi
        else:
            x1lon = lonrad + math.pi
    else:
        x1lon = lonrad
        x1lat = latrad - (math.pi/2)
    if lonrad < (-1 * (math.pi/2)):
        x2lon = lonrad + 3 * (math.pi/2)
    else:
        x2lon = lonrad - (math.pi/2)
    x2lat = (math.pi/2)
    x3lon = lonrad
    x3lat = latrad
    
    # Calculate transformation matrix
    
    a11 = math.sin(x1lat) * math.cos(-1 * x1lon)
    a12 = math.sin(x2lat) * math.cos(-1 * x2lon)
    a13 = math.sin(x3lat) * math.cos(-1 * x3lon)
    a21 = math.sin(x1lat) * math.cos((math.pi/2) - x1lon)
    a22 = math.sin(x2lat) * math.cos((math.pi/2) - x2lon)
    a23 = math.sin(x3lat) * math.cos((math.pi/2) - x3lon)
    a31 = math.cos(x1lat)
    a32 = math.cos(x2lat)
    a33 = math.cos(x3lat)
    
    j = 0
    
    for i in range(negmag, posmag, step):
        
        # Define translation vector in spherical coordinates
        
        trad = i
        ttheta = diprad
        tphi = azrad
        
        # Convert translation vector to cartesian coordinates
        
        tx = trad * math.sin(ttheta) * math.cos(tphi)
        ty = trad * math.sin(ttheta) * math.sin(tphi)
        tz = trad * math.cos(ttheta)
        
        # Transform translation vector into base coordinate system
        
        txnew = a11 * tx + a12 * ty + a13 * tz
        tynew = a21 * tx + a22 * ty + a23 * tz
        tznew = a31 * tx + a32 * ty + a33 * tz
        
        # Add new vector to original position vector
        
        eptx = cx + txnew
        epty = cy + tynew
        eptz = cz + tznew
        
        # Convert new sum to spherical coordinates
        
        eptrad = math.sqrt(math.pow(eptx, 2) + math.pow(epty, 2) + math.pow(eptz, 2))
        eptphirad = math.atan2(epty, eptx)
        eptthetarad = math.acos(eptz / (math.sqrt(math.pow(eptx, 2) + math.pow(epty, 2) + math.pow(eptz, 2))))
        
        # Convert into lat, lon, depth
        
        eptdep = 6371 - eptrad
        eptlat = 90 - (math.degrees(eptthetarad))
        eptlon = math.degrees(eptphirad)
        
        # Populate EPs
        
        EPs[j, 0] = eptlon
        EPs[j, 1] = eptlat
        EPs[j, 2] = eptdep
        EPs[j, 3] = i
        
        j = j + 1
    
    return EPs

def mkSDgrdTomo(Slabgrid):
    gdict = Slabgrid.getGeoDict().copy()
    dx = gdict.dx * 111.19
    dy = gdict.dy * 111.19

    Xgrad, Ygrad = np.gradient(Slabgrid.getData().copy(),dx,dy)
    Maggrid = np.sqrt((Ygrad**2)+(Xgrad**2))

    # Define a grid file that is the direction perpendicular to the max gradient
    Dirgrid = np.degrees(np.arctan2(Ygrad, Xgrad))
    Dirgrid = np.where(Dirgrid < 0, Dirgrid+360, Dirgrid)

    Strgrid = (Dirgrid.copy()-90)*-1
    #Strgrid = Dirgrid.copy()+180
    Strgrid[Strgrid<0] = Strgrid[Strgrid<0]+360

    # Assign strike and dip arrays to grids with same dimensions as depth grid
    Dipgrid = gmt.GMTGrid(np.degrees(np.arctan2(Maggrid, 1)), Slabgrid.getGeoDict().copy())
    Strikegrid = gmt.GMTGrid(Strgrid, Slabgrid.getGeoDict().copy())

    return Strikegrid, Dipgrid


def makeCube(Tomo,HRES,VRES):
    # Define input into array

    data = np.zeros((len(Tomo),4))
    data[:,0] = Tomo[:,1]
    data[:,1] = Tomo[:,0]
    data[:,2] = Tomo[:,2]
    data[:,3] = Tomo[:,3]

    deplist = sorted(set(data[:,2]))

    # Set up interpolated grid

    points = data[:,0:3]
    values = data[:,3]
    ndep = len(deplist)

    lonmin = points[:,0].min()
    lonmax = points[:,0].max()
    latmin = points[:,1].min()
    latmax = points[:,1].max()

    xpts = np.arange(np.floor(lonmin),np.ceil(lonmax)+HRES,HRES)
    ypts = np.arange(np.floor(latmin),np.ceil(latmax)+HRES,HRES)

    n = len(xpts)
    m = len(ypts)

    xpts, ypts = np.meshgrid(xpts, ypts)

    xi = np.zeros((m*n,2))
    xi[:,0] = xpts.flatten()
    xi[:,1] = ypts.flatten()

    cube = np.zeros((m,n,ndep))

    for i in range(0,ndep):
        dep = deplist[i]
        inum = points[:,2] == dep
        pnum = points[inum,0:2]
        vnum = values[inum]

        zinterp = griddata(pnum,vnum,xi,method='linear')
        zinterp.shape = (m,n)
        cube[:,:,i] = zinterp

    # Create 1D interpolated functions and fill output cube

    depmin = min(deplist)
    depmax = max(deplist)
    zpts = np.arange(np.floor(depmin),np.ceil(depmax)+VRES,VRES)
    o = len(zpts)
    output = np.zeros((m,n,o))

    for i in range(0,m-1):
        for j in range(0,n-1):
            ftmp = interpolate.interp1d(deplist,cube[i,j,:],kind='linear',fill_value=np.nan,bounds_error=False)
            output[i,j,:] = ftmp(zpts)

    ulx = np.floor(lonmin)
    uly = np.floor(latmin)
    ulz = np.floor(depmin)
    xdim = HRES
    ydim = HRES
    zdim = VRES

    return output, ulx, uly, ulz, xdim, ydim, zdim

# getRowColBand was written by Mike Hearne and edited by Daniel Portner.  It is to be used in conjunction with get3dvalue (below), providing get3dvalue with the position of a point on earth with its location in a gridded data cube

def getRowColBand(cube,lon,lat,depth,ulx,uly,ulz,xdim,ydim,zdim):
    nrows,ncols,nbands = cube.shape
    col = (lon-ulx)/xdim
    row = (lat-uly)/ydim
    band = (depth-ulz)/zdim

    return (row,col,band)

# get3dvalue was written by Mike Hearne and edited by Daniel Portner.  It takes the 3d position within a cube and does a weighted average 3d interpolation to determine the value at a given point within a 3d gridded data cube

def get3dvalue(cube,row,col,band):

    nrows,ncols,nbands = cube.shape
    if row > (nrows - 1):
        return np.nan
    elif row < 0:
        return np.nan
    elif col > (ncols - 1):
        return np.nan
    elif col < 0:
        return np.nan
    elif band > (nbands - 1):
        return np.nan
    elif band < 0:
        return np.nan
    else:

        #get the values and coordinates of the 8 points that surround the input
        #point at (xi,yi,zi)
        row0 = int(np.floor(row))
        row1 = int(np.ceil(row))
        col0 = int(np.floor(col))
        col1 = int(np.ceil(col))
        band0 = int(np.floor(band))
        band1 = int(np.ceil(band))
        points = np.array([[row0,col0,band0],
                           [row0,col1,band0],
                           [row1,col0,band0],
                           [row1,col1,band0],
                           [row0,col0,band1],
                           [row0,col1,band1],
                           [row1,col0,band1],
                           [row1,col1,band1]],dtype=np.int32)
        values = [cube[row0,col0,band0],
                  cube[row1,col1,band0],
                  cube[row1,col0,band0],
                  cube[row1,col1,band0],
                  cube[row0,col0,band1],
                  cube[row1,col1,band1],
                  cube[row1,col0,band1],
                  cube[row1,col1,band1]]
        t1 = np.power((row-points[:,0]),2)
        t2 = np.power((col-points[:,1]),2)
        t3 = np.power((band-points[:,2]),2)
        d = np.sqrt(t1+t2+t3)
        sumd = np.sum(d)
        nm1 = len(values)-1
        weights = (sumd - d)/sumd
        di = np.sum((values*weights)/nm1)

        return di

def myround(x, vec):
    tmp = np.abs(vec - x)
    indx = tmp == np.min(tmp)
    minvec = vec[indx]
    val = np.max(minvec)
    return val
