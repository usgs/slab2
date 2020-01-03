#!/usr/bin/env python

### Created DEP.7.23.15 ###
### Updated DEP.6.14.16 ###

# This is my repository of functions for use in my slab tomography code

###########################

### Import third party modules

import numpy as np
import math as m
import neicio.gmt as gmt
from scipy import interpolate
from scipy.interpolate import griddata
import os
from matplotlib import path
from scipy import ndimage

###########################

# epCalc is used to calculate an array of locations and distances from the slab surface by calculating vector endpoints within the earth
# Edited DEP.7.1.16 to accommodate changed mkSDgrd

def epCalc(lon,lat,dep,dip,strike,posmag,negmag,step):
    EPs = np.zeros(((posmag-negmag)/step,4))
    
    # Rotate from strike to direction of motion
    
    if strike > 270:
        az = (strike + 90) - 360
    else:
        az = strike + 90
    az = 360 - az    # Accounts for the fact that azimuth counts goes opposite of the positive rotation of the x-axis (points north)

    # Convert input angles to radians

    latrad = m.radians(90 - lat)
    lonrad = m.radians(lon)
    azrad = m.radians(az)
    diprad = m.radians(dip)

    # Define initial location in spherical coordinates

    crad = 6371 - dep
    ctheta = latrad
    cphi = lonrad

    # Convert initial location to cartesian coordinates

    cx = crad * m.sin(ctheta) * m.cos(cphi)
    cy = crad * m.sin(ctheta) * m.sin(cphi)
    cz = crad * m.cos(ctheta)

    # Define lon/lat of new coordinate system

    if latrad < (m.pi/2):
        x1lat = abs(latrad-(m.pi/2))
        if lonrad > 0:
            x1lon = lonrad - m.pi
        else:
            x1lon = lonrad + m.pi
    else:
        x1lon = lonrad
        x1lat = latrad - (m.pi/2)
    if lonrad < (-1 * (m.pi/2)):
	x2lon = lonrad + 3 * (m.pi/2)
    else:
	x2lon = lonrad - (m.pi/2)
    x2lat = (m.pi/2)
    x3lon = lonrad
    x3lat = latrad

    # Calculate transformation matrix

    a11 = m.sin(x1lat) * m.cos(-1 * x1lon)
    a12 = m.sin(x2lat) * m.cos(-1 * x2lon)
    a13 = m.sin(x3lat) * m.cos(-1 * x3lon)
    a21 = m.sin(x1lat) * m.cos((m.pi/2) - x1lon)
    a22 = m.sin(x2lat) * m.cos((m.pi/2) - x2lon)
    a23 = m.sin(x3lat) * m.cos((m.pi/2) - x3lon)
    a31 = m.cos(x1lat)
    a32 = m.cos(x2lat)
    a33 = m.cos(x3lat) 

    j = 0

    for i in range(negmag,posmag,step):
        
        # Define translation vector in spherical coordinates

        trad = i
        ttheta = diprad
        tphi = azrad
    
    	# Convert translation vector to cartesian coordinates
    
        tx = trad * m.sin(ttheta) * m.cos(tphi)
        ty = trad * m.sin(ttheta) * m.sin(tphi)
        tz = trad * m.cos(ttheta)

    	# Transform translation vector into base coordinate system

        txnew = a11 * tx + a12 * ty + a13 * tz
        tynew = a21 * tx + a22 * ty + a23 * tz
        tznew = a31 * tx + a32 * ty + a33 * tz

    	# Add new vector to original position vector

        eptx = cx + txnew
        epty = cy + tynew
        eptz = cz + tznew

    	# Convert new sum to spherical coordinates

        eptrad = m.sqrt(m.pow(eptx,2) + m.pow(epty,2) + m.pow(eptz,2))
        eptphirad = m.atan2(epty,eptx)
        eptthetarad = m.acos(eptz / (m.sqrt(m.pow(eptx,2) + m.pow(epty,2) + m.pow(eptz,2))))

    	# Convert into lat, lon, depth

        eptdep = 6371 - eptrad
        eptlat = 90 - (m.degrees(eptthetarad))
        eptlon = m.degrees(eptphirad)

    	# Populate EPs

        EPs[j,0] = eptlon
        EPs[j,1] = eptlat
        EPs[j,2] = eptdep
        EPs[j,3] = i

        j = j + 1

    return EPs

###########################

# makeCube is used to create a 3D cube out of an input tomography dataset

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

###########################

# getRowColBand was written by Mike Hearne and edited by Daniel Portner.  It is to be used in conjunction with get3dvalue (below), providing get3dvalue with the position of a point on earth with its location in a gridded data cube

def getRowColBand(cube,lon,lat,depth,ulx,uly,ulz,xdim,ydim,zdim):
    nrows,ncols,nbands = cube.shape
    col = (lon-ulx)/xdim
    row = (lat-uly)/ydim
    band = (depth-ulz)/zdim

    return (row,col,band)

###########################

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
	row0 = np.floor(row)
	row1 = np.ceil(row)
	col0 = np.floor(col)
	col1 = np.ceil(col)
	band0 = np.floor(band)
	band1 = np.ceil(band)
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

#########################

# mkSDgrd is used to take a gridfile (particularly one meant to represent a surface, such as Slab1.0) and create grids of similar structure that represent the strike and dip of the the original grid.
# Edited DEP.7.1.16 to fix bugs

def mkSDgrd(InGrd):

	Slabgrid = InGrd
	
	dx = Slabgrid.geodict['xdim'] * 111.19  ## Define dx as the longitude increments and convert to km
	dy = Slabgrid.geodict['ydim'] * 111.19  ## Define dy as the latitude increments and convert to km
	Xgrad, Ygrad = np.gradient(Slabgrid.griddata,dx,dy)  ## Define grids of the gradients in the x and y directions
	Maggrid = np.sqrt((Xgrad**2)+(Ygrad**2))  ## Define a grid file that is the gradient magnitude at any given node
	Dirgrid = -1*(np.degrees(np.arctan2(Ygrad,Xgrad))-90)  ## Define a grid file that is the direction perpendicular to the max gradient - in degrees clockwise from N (0 - 360)
	Dirgrid[Dirgrid < 0] += 360
	
	Dipgrid = gmt.GMTGrid()
	Strikegrid = gmt.GMTGrid()
	Dipgrid.griddata = np.degrees(np.arctan2(Maggrid,1))
	Strikegrid.griddata = Dirgrid
	Dipgrid.geodict = Slabgrid.geodict.copy()
	Strikegrid.geodict = Slabgrid.geodict.copy()

	return Strikegrid, Dipgrid

#########################

# point_in_poly was written by some guy on the internet at http://geospatialpython.com/2011/01/point-in-polygon.html.  Its purpose is to identify whether or not a point is located within a provided polygon.  I am using it to create an array of 1s and 0s based on a clipping mask.  Takes in lists for x and y.

def point_in_poly(x_all,y_all,poly):

    ans = np.zeros((len(x_all),1))
    for j in range(0,len(ans)):
	x = x_all[j]
	y = y_all[j]
	n = len(poly)
	
	p1x,p1y = poly[0]
	cnt = 0
        for i in range(n+1):
            inside = False
	    p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
		    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x,p1y = p2x,p2y
	    if inside == True:
	        cnt = cnt + 1
        if (cnt % 2) == 0:
	    ans[j] = np.nan
	else:
	    ans[j] = 1
	
    return ans

#########################

# mySurface calls gmt from the command line to take an irregularly spaced dataset and create a regularized surface out of it.  DEP.6.29.16

def mySurface(data,node,T):

	### Input:
	# data: (n x 3) numpy array of form x,y,z where x,y are coordinates and z is value
	# node: float to indicate node spacing *note that this code uses only square grids
	# T: tension factor (float from 0-1) where 0 = minimum curvature solution

	### Output:
	# newgrid: gmt-style grid after NEICIO format that is a surface of input data
	# xi: numpy array with x,y coordinates for each node	

	xmin, xmax = data[:,0].min(), data[:,0].max()
	ymin, ymax = data[:,1].min(), data[:,1].max()

	xall = np.arange(np.floor(xmin),np.ceil(xmax)+node/2,node)
	yall = np.arange(np.floor(ymin),np.ceil(ymax)+node,node)

	n = len(xall)
	m = len(yall)

	xpts, ypts = np.meshgrid(xall, yall)

	xi = np.zeros((m*n,2))
	xi[:,0] = xpts.flatten()
	xi[:,1] = ypts.flatten()
 
	np.savetxt('temp.xyz',data,fmt='%0.4f',delimiter=' ')
	rflag="-R%s/%s/%s/%s" %(np.floor(xmin),np.ceil(xmax),np.floor(ymin),np.ceil(ymax))
	iflag="-I%s/%s" %(node,node)
	gflag="-Gtemp.grd"
	tflag="-T%s" %(T)

	os.system("gmt4 blockmean temp.xyz %s %s | gmt4 surface %s %s %s %s" %(rflag,iflag,rflag,iflag,gflag,tflag))

	newgrid = gmt.GMTGrid('temp.grd')
	
	os.system("rm temp.xyz")
	os.system("rm temp.grd")

	return newgrid,xi

#########################

# mkMask creates a grid of 1 or nan based on if they are within a clipping mask or not.  DEP.6.29.16

def mkMask(clipfile,xi,shape):

	### Input:
	# clipfile: a string identifying an ascii file of x,y coordinates denoting a polygon/clipping mas
	#xi: (n x 2) numpy array of coordinates that should define a regular grid
	#shape: tuple (m,n) of shape of a grid - n*m should be equal to the length of xi

	### Output:
	#mask: (m x n) numpy array of 1s and nans where 1s are within the clipping mask

	clip = np.loadtxt(clipfile)
	lons = clip[:,0]
	lons[lons > 180] = lons - 360
	
	xy = zip(lons,clip[:,1])
	poly = path.Path(xy)
	temp = poly.contains_points(xi)
	mask = np.zeros(len(temp),) * np.nan
	mask[temp] = 1
	mask.shape = shape

	return mask

#########################

# clipGrid applies a clipping mask to a neicio gmt format grid.  DEP.6.29.16

def clipGrid(ingrid,mask):

	### Input:
	# ingrid: a gmt style grid such as created by gmt.GMTgrid() or by mySurface above
	# mask: clipping mask grid from mkMask

	### Output:
	# outgrid: same as ingrid, but clipped

	outgrid = gmt.GMTGrid()
	outgrid.griddata = np.multiply(ingrid.griddata,np.flipud(mask))
	outgrid.geodict = ingrid.geodict.copy()

	return outgrid

##########################

# gaussGridFilt takes in a square grid (can't deal with nans) and filters it using a Gaussian of width provided in degrees.  DEP.6.29.16

def gaussGridFilt(ingrid,filterwidth):

	### Input:
	# ingrid: a gmt style grid such as created by gmt.GMTgrid() or by mySurface above which has not been clipped using a clipping mask.  ASSUMES SQUARE GRIDS
	# filterwidth: width (in same units as grid spacing) of gaussian filter.  equivalent to 2*sigma of the gaussian

	### Output:
	# outgrid: same as ingrid, but filtered (smoothed)

	outgrid = gmt.GMTGrid()
	outgrid.geodict = ingrid.geodict.copy()
	sigma = (filterwidth/2) / outgrid.geodict['xdim']
	outgrid.griddata = ndimage.filters.gaussian_filter(ingrid.griddata,sigma)

	return outgrid

###########################

# distCalc3D is used to calculate the 3D straight-line distance between two points within the earth given lat,lon,dep for each point.  DEP.7.6.16

def distCalc3D(lon1,lat1,dep1,lon2,lat2,dep2):

    # Convert input angles to radians

    lat1rad = m.radians(90 - lat1)
    lon1rad = m.radians(lon1)
    lat2rad = m.radians(90 - lat2)
    lon2rad = m.radians(lon2)

    # Define initial location in spherical coordinates - point1

    rad1 = 6371 - dep1
    theta1 = lat1rad
    phi1 = lon1rad

    # Convert initial location to cartesian coordinates

    x1 = rad1 * m.sin(theta1) * m.cos(phi1)
    y1 = rad1 * m.sin(theta1) * m.sin(phi1)
    z1 = rad1 * m.cos(theta1)

    # Define initial location in spherical coordinates - point2

    rad2 = 6371 - dep2
    theta2 = lat2rad
    phi2 = lon2rad

    # Convert initial location to cartesian coordinates

    x2 = rad2 * m.sin(theta2) * m.cos(phi2)
    y2 = rad2 * m.sin(theta2) * m.sin(phi2)
    z2 = rad2 * m.cos(theta2)

    ans = (m.sqrt(m.pow(x2-x1,2) + m.pow(y2-y1,2) + m.pow(z2-z1,2)))

    return ans

###########################

# pointShift is essentially epCalc, but for a single point.  It is used to calculate the endpoint of a vector within the earth given a local lat/lon/dep, strike/dip, and distance.  DEP.7.7.16

def pointShift(lon,lat,dep,dip,strike,mag):

    # Rotate from strike to direction of motion
 
    if strike > 270:
        az = (strike + 90) - 360
    else:
        az = strike + 90
    az = 360 - az    # Accounts for the fact that azimuth counts goes opposite of the positive rotation of the x-axis (points north)

    # Convert input angles to radians

    latrad = m.radians(90 - lat)
    lonrad = m.radians(lon)
    azrad = m.radians(az)
    diprad = m.radians(dip)

    # Define initial location in spherical coordinates

    crad = 6371 - dep
    ctheta = latrad
    cphi = lonrad

    # Convert initial location to cartesian coordinates

    cx = crad * m.sin(ctheta) * m.cos(cphi)
    cy = crad * m.sin(ctheta) * m.sin(cphi)
    cz = crad * m.cos(ctheta)

    # Define lon/lat of new coordinate system

    if latrad < (m.pi/2):
        x1lat = abs(latrad-(m.pi/2))
        if lonrad > 0:
            x1lon = lonrad - m.pi
        else:
            x1lon = lonrad + m.pi
    else:
        x1lon = lonrad
        x1lat = latrad - (m.pi/2)
    if lonrad < (-1 * (m.pi/2)):
        x2lon = lonrad + 3 * (m.pi/2)
    else:
        x2lon = lonrad - (m.pi/2)
    x2lat = (m.pi/2)
    x3lon = lonrad
    x3lat = latrad

    # Calculate transformation matrix

    a11 = m.sin(x1lat) * m.cos(-1 * x1lon)
    a12 = m.sin(x2lat) * m.cos(-1 * x2lon)
    a13 = m.sin(x3lat) * m.cos(-1 * x3lon)
    a21 = m.sin(x1lat) * m.cos((m.pi/2) - x1lon)
    a22 = m.sin(x2lat) * m.cos((m.pi/2) - x2lon)
    a23 = m.sin(x3lat) * m.cos((m.pi/2) - x3lon)
    a31 = m.cos(x1lat)
    a32 = m.cos(x2lat)
    a33 = m.cos(x3lat)

    # Define translation vector in spherical coordinates

    trad = mag
    ttheta = diprad
    tphi = azrad

    # Convert translation vector to cartesian coordinates

    tx = trad * m.sin(ttheta) * m.cos(tphi)
    ty = trad * m.sin(ttheta) * m.sin(tphi)
    tz = trad * m.cos(ttheta)

    # Transform translation vector into base coordinate system

    txnew = a11 * tx + a12 * ty + a13 * tz
    tynew = a21 * tx + a22 * ty + a23 * tz
    tznew = a31 * tx + a32 * ty + a33 * tz

    # Add new vector to original position vector

    eptx = cx + txnew
    epty = cy + tynew
    eptz = cz + tznew

    # Convert new sum to spherical coordinates

    eptrad = m.sqrt(m.pow(eptx,2) + m.pow(epty,2) + m.pow(eptz,2))
    eptphirad = m.atan2(epty,eptx)
    eptthetarad = m.acos(eptz / (m.sqrt(m.pow(eptx,2) + m.pow(epty,2) + m.pow(eptz,2))))

    # Convert into lat, lon, depth

    eptdep = 6371 - eptrad
    eptlat = 90 - (m.degrees(eptthetarad))
    eptlon = m.degrees(eptphirad)

    return eptlon,eptlat,eptdep

