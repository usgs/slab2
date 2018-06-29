#!/usr/bin/env python

"""This is a compilation of all functions used by distributed Slab2.0 code (and maybe some extras), including:
    2) s2d.py
    2) slab2.py
    3) tomo_slab.py

    The functions were written variably by Ginevra Moore (GLM), Maria Furtney (MF) and Daniel Portner (DEP).
    """

###############################################

### Module imports

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mapio.gmt as gmt

import math
import os
import csv
import urllib.request, urllib.error, urllib.parse
import os.path
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib import path
from scipy import ndimage
from shapely.geometry import Polygon 
from pandas import DataFrame
from obspy.geodetics.base import gps2dist_azimuth

from scipy.interpolate import LSQBivariateSpline
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import LSQSphereBivariateSpline
from functools import partial
from multiprocess import Pool
import psutil
import gc
from sklearn import mixture
from sklearn.metrics import mean_squared_error
import warnings
import matplotlib.patches as patches
import utm
from datetime import datetime
from scipy.interpolate import Rbf
from copy import deepcopy
from pylab import arccos,argsort,cross,dot,double,eigh,pi,trace,zeros
###
# The following functions cannot be translated to lon,lat:
#   getValue (imported from external library)



###############################################

### 1 ###

###############################################

## Written GLM
## Modified DEP.8.4.16

def get_grid(slr, destdi):
    
    ''' Arguments:  slr - slab indicator (i.e. sam) - indicates which slab to get
                    desdi - 'depth', 'strike', or 'dip' - indicates which grid to get
        
        Returns:    clipped depth, strike, or dip GMT grid file from Slab 1.0 '''
    
    if destdi == 'depth':
        clip = 'clip'
    elif destdi == 'strike':
        clip = 'strclip'
    elif destdi == 'dip':
        clip = 'dipclip'
    
    url = 'http://earthquake.usgs.gov/data/slab/models/%s_slab1.0_%s.grd' % (slr, clip)
    fname = 'library/slab1grids/%s_slab1.0_%s.grd' % (slr, clip)
    
    # If the grid file is not in the specified directory, it is found from
    # a query of the Slab 1.0 website
    if not os.path.isfile(fname):
        fh = urllib.request.urlopen(url)
        data = fh.read()
        fh.close()
        f = open(fname, 'wb')
        f.write(data)
        f.close()
    
    # Load file into GMT grid
    depgrid = gmt.GMTGrid.load(fname)
    return depgrid


###############################################

### 2 ###

###############################################

## Written GLM
## Modified DEP.8.5.16

def getEventsInCircle(lon, lat, radius, eventlist):
    
    ''' Arguments:  lat - latitude of grid node that is being searched over
                    lon - longitude of grid node that is being searched over
                    radius - radius of circle to search within (km)
                    eventlist - list of events to search over. Must include lat/lon info
        
        Returns:    elist - dataframe of events that are within the specified radius of
                            the lat lon node point. Has all information that the original
                            eventlist contained for nearby data points of all types. '''
    
    # Gather latitudes and longitudes for each point in eventlist        
    lons = eventlist['lon'].values*1.0
    lats = eventlist['lat'].values*1.0
    eventlist['distance'], cosangles = npcosine(lon, lat, lons, lats)

    elist = eventlist.loc[eventlist.distance <= radius]
    
    return elist

###############################################

### 3 ###

###############################################

## Written GLM
# Modified DEP.8.5.16

def getEventsInEllipse(lon, lat, strk, aval, bval, eventlist, lon1, lat1):
    
    ''' Arguments:  lat - latitude of grid node that is being searched over
                    lon - longitude of grid node that is being searched over
                    strk - local strike of slab at this grid node
                    aval - long radius of ellipse to search within
                    bval - short radius of ellipse to search within
                    eventlist - list of events to search over. Must include lat/lon info
        
        Returns:    elist - dataframe of events that are within the specified ellipse around
                            the lat lon node point. Has all information that the original
                            eventlist contained for nearby data points of all types. '''
    
    # Gather latitudes and longitudes for each point in eventlist
    lons = eventlist['lon'].values*1.0
    lats = eventlist['lat'].values*1.0
    eventlist['distance2'], az = npcosine(lon, lat, lons, lats)
    
    mdist = []
    erta = math.sqrt(1-((math.pow(bval, 2))/(math.pow(aval, 2))))
    mdist = getEllipseRad(aval, erta, az, strk)
    
    eventlist['azimuth'] = az
    elist = eventlist.loc[eventlist.distance2 <= mdist]
    elist = elist[['lon', 'lat', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
    
    '''
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111)
    ax1.plot(eventlist['lon'].values,eventlist['lat'].values,'co',label='Prefiltered')
    ax1.plot(elist['lon'].values,elist['lat'].values,'ko',label='filtered')
    ax1.plot(lon1,lat1,'ro',label='node')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.axis('equal')
    plt.grid()
    title = 'Lat: %.1f, Lon: %.1f, aval: %i, bval: %i, strk: %.1f' % (lat1, lon1, aval, bval, strk)
    ax1.set_title(title)
    ax1.legend(loc='best')
    
    figtitle = 'Output/multitest_cam/pdf%i%i_%i_%i.png' % (int(lon1*100),int(lat1*100),aval,bval)

    #fig.savefig(figtitle)
    plt.close()
    '''
    return elist

###############################################

### 4 ###

###############################################

## Written GLM

def getEllipseRad(a, e, azi, ang):
    
    ''' Arguments:  a - long radius of the ellipse
                    e - eccentricity of the ellipse
                    azi - azimuth between two points
                    ang - local strike value
        
        Returns:    d2 - maximum radius of ellipse at this point given the azimuth and strike '''
    
    d2r = math.pi/180
    d2 = (a*(1-e**2))/(1+(e*(np.cos((azi-ang)*d2r))))
    
    return d2

###############################################

### 5 ###

###############################################

## Written GLM

def getstrike(lon, lat, strgrid):
    
    ''' Arguments:  lat - latitude of grid node that is being searched over
                    lon - longitude of grid node that is being searched over
                    strgrid - GMT surface of strike values
        
        Returns:    cstr - nearest local strike at this grid node from Slab 1.0
                            Used when there is not a local strike at the exact lat/lon grid point'''
    
    # Searches in concentric circles to find the nearest strike value to this grid node
    for radius in np.arange(0, 1, .05): #degrees
        for theta in np.arange(0, 2*np.pi, np.pi/4.): #radians
            lat2 = lat + radius*math.sin(theta)
            lon2 = lon + radius*math.cos(theta)
            cstr = strgrid.getValue(lat2, lon2)
            if cstr != np.nan:
                return cstr
    return cstr

###############################################

### 6 ###

###############################################

## Written GLM

def getdip(lat, lon, dipgrid):
    
    ''' Arguments:  lat - latitude of grid node that is being searched over
                    lon - longitude of grid node that is being searched over
                    dipgrid - GMT surface of strike values
        
        Returns:    dip - nearest local strike at this grid node from Slab 1.0
                    Used when there is not a local strike at the exact lat/lon grid point'''
    
    # Searches in concentric circles to find the nearest dip value to this grid node
    # Returns nan if a valid nearby dip cannot be found
    dip = np.NaN
    for radius in np.arange(0, 1, .05): #degrees
        for theta in np.arange(0, 2*np.pi, np.pi/4.): #radians
            lat2 = lat + radius*math.sin(theta)
            lon2 = lon + radius*math.cos(theta)
            try:
                dip = dipgrid.getValue(lat2, lon2)
            except:
                dip = np.NaN
            if np.isfinite(dip) and dip is not None:
                return dip
    return dip

###############################################

### 7 ###

###############################################

## Written GLM (Translated From G. Hayes Perl Script)

def heading(lon, lat, dist, az):
    
    ''' Arguments:  lon - longitude of known point
                    lat - latitude of known point
                    dist - distance from known point
                    az - azimuth from known point
        
        Returns:    lat - latitude of new point as projected by a certain azimuth and great circle
                            distance from known lat/lon point
                    lon - longitude of new point as projected by a certain azimuth and great circle
                            distance from known lat/lon point  '''
    
    # Creating degrees/radians conversion constants
    d2r = math.pi/180
    r2d = 180/math.pi
    
    # Ensuring that distances are positive
    if(dist<0):
        dist = dist*-1
        az = az-180
    
    # Ensuring that azimuths are between 0 and 360
    if(az<0):
        az = az+360
    elif(az>360):
        az = az-360
    
    # Finding projected latitude
    b = (90-lat)*d2r
    a = (dist/111.19)*d2r
    angC = az*d2r
    c = math.cos(a)*math.cos(b)+math.sin(a)*math.sin(b)*math.cos(angC)
    c = math.acos(c)
    cdeg = c*r2d
    lat1 = 90-cdeg
    
    # Finding projected longitude
    angA = (math.cos(a)-(math.cos(b)*math.cos(c)))/(math.sin(b)*math.sin(c))
    zeroonebool = False
    if angA > 1.0:
        angA-=0.00001
        zeroonebool = True
    if angA < -1.0:
        angA+=0.00001
        zeroonebool = True
    angA = math.acos(angA)
    adeg = angA*r2d
    if(az>0 and az<=180):
        lon1 = lon+adeg
    else:
        lon1 = lon-adeg
    
    return lon1, lat1

###############################################

### 8 ###

###############################################

## Written GLM

def datelinecross(x):
    ''' Arguments:  x - longitude value (positive or negative)
        
        Returns:    x - a positive longitude. Stays the same if the input was positive,
                        is changed to positive if the input was negative '''
    
    if x<0:
        return x+360
    else:
        return x

###############################################

### 9 ###

###############################################

## Written GLM

def meridiancross(x):
    ''' Arguments:  x - longitude value (positive or negative)
        
        Returns:    x - a longitude in the -180/180 domain '''
    
    if x>180:
        return x-360
    else:
        return x

def northcross(x):
    ''' Arguments:  x - longitude value (positive or negative)
        
        Returns:    x - a longitude in the -180/180 domain '''
    
    if x<90:
        return x+360
    else:
        return x

def unnorthcross(x):
    ''' Arguments:  x - longitude value (positive or negative)
        
        Returns:    x - a longitude in the -180/180 domain '''
    
    if x>360:
        return x-360
    else:
        return x

def zerothreesixty(data):
    data['lon']=data.apply(lambda row: datelinecross(row['lon']), axis=1)
    return data

def oneeighty(data):
    data['lon']=data.apply(lambda row: meridiancross(row['lon']), axis=1)
    return data

def northernaz(data):
    data['az']=data.apply(lambda row: northcross(row['az']), axis=1)
    return data

def notnorthanymore(data):
    data['az']=data.apply(lambda row: unnorthcross(row['az']), axis=1)
    return data

###############################################

### 10 ###

###############################################

## Written GM

def cosrule(d2r, lat1, lon1, lat2, lon2):
    
    ''' Arguments:  d2r - degree to radians conversion constant
                    lat1 - latitude point that the angle is referenced from
                    lon1 - longitude point that the angle is referenced from
                    lat2 - latitude point that the angle is going to (clockwise from 0 degrees)
                    lon2 - longitude point that the angle is going to (clockwise from 0 degrees)
        
        Returns:    dist2 - great circle distance between the two lat/lon points
                    ang - the angle between the two points (clockwise from 0 degrees from lat1/lon1 point) '''
    
    if abs(lon1-lon2) < 0.00001 or abs(lat1-lat2) < 0.00001:
        lat2 = lat2+0.0001
        lon2 = lon2+0.0001
    
    cl1 = (90-lat1) * d2r
    cl2 = (90-lat2) * d2r
    dlon = (lon2-lon1) * d2r
    dist = math.cos(cl1) * math.cos(cl2) + math.sin(cl1) * math.sin(cl2) * math.cos(dlon)
    if dist < -1:
        dist = -1.0
    if dist > 1:
        dist = 1.0
    dist2 = math.acos(dist)
    if dlon > math.pi:
        dist2 = 2 * math.pi-dist2
    if dist != 0:
        ang = (math.cos(cl2) - (dist * math.cos(cl1))) / (math.sin(dist2) * math.sin(cl1))
    else:
        ang = 1.0
    if ang < -1:
        ang = -1.0
    if ang > 1:
        ang = 1.0
    ang = math.acos(ang)
    return dist2, ang

###############################################

### 10 ###

###############################################

## Written GM

def cosine(lon1, lat1, lon2, lat2):
    
    ''' Arguments:  lon1 - latitude point that the angle is referenced from
                    lat1 - longitude point that the angle is referenced from
                    lon2 - latitude point that the angle is going to (clockwise from 0 degrees)
                    lat2 - longitude point that the angle is going to (clockwise from 0 degrees)
        
        Returns:    dist - great circle distance between the two lat/lon points
                    ang - the angle between the two points (clockwise from 0 degrees from lat1/lon1 point)
                    lon1 - same longitude as input argument, used in some applications
                    lat1 - same latitude as input argument, used in some applications  '''
    
    # Ensuring that azimuths are between 0 and 360
    if lon1 < 0:
        lon1 = lon1 + 360
    if lon2 < 0:
        lon2 = lon2 + 360

    # Creating degrees/radians conversion constants
    d2r = (math.pi/180)
    r2d = (180/math.pi)
    ddlon = lon1 - lon2

    # Getting distance and angle between the two points (in degrees)
    dist, ang = cosrule(d2r, lat1, lon1, lat2, lon2)
    if lon1 > lon2 and ddlon < 180:
        ang = 2*math.pi - ang
    dist = abs(dist*r2d)
    if dist > 180:
        dist = 360 - dist
        ang = ang + math.pi
    if ang > 2*math.pi:
        ang = 2*math.pi - ang
    dist = dist * 111.19
    ang = ang * r2d
    
    return dist, ang, lat1, lon1

###############################################

### 11 ###

###############################################

## Written GM

def round_to_nths(num, n):
    return int(num*n)/n

###############################################

### 12 ###

###############################################

## Written GM
## Modified MF.8.4.16
## Edited DEP.8.12.16 removing "yes" output

def trenchDist(slat, slon, sstrike, TR_data):
    
    ''' Arguments:  slat - latitude of the grid node
                    slon - longitude of the grid node
                    sstrike - local strike at the grid node
                    TR_data - local trench file (contains lat,lon,azimuth information)
                    minlon, maxlon - longitude bounds of slab model
                    minlon, minlat - latitude bounds of slab model
        
        Returns:    lat3 - latitude of nearest trench point
                    lon3 - longitude of nearest trench point
                    dist3 - distance between node and nearest trench point
                    yes - indicates whether or not a valid point was found. If not valid, point is not used
                    raz - azimuth between this node and the nearest trench point (clockwise from 90 deg)  '''
    
    # Ensuring that all longitude values are positive in the trench file (0-360)
    TR_data['lon']=TR_data.apply(lambda row: datelinecross(row['lon']), axis=1)
    
    # Storing longitude/latitude values of trench files in array
    lons2 = TR_data['lon'].values
    lats2 = TR_data['lat'].values
    
    # Setting minimum values and longitude bounds relative to grid node
    minlon2 = slon-1
    maxlon2 = slon+1
    mindaz = 999
    dist3 = 9999
    raz = 0
    
    # Finding local strike-90 - search along this azimuth to find nearest trench point
    az = sstrike-90
    if az<0:
        az = az+360
    
    # Looping through local trench points to find nearest one to grid node
    for i in range(len(lons2)-1):
        if lons2[i]>minlon2 and lons2[i]<maxlon2:
            ssdist, saz, lat3, lon3 = cosine(slon, slat, lons2[i], lats2[i])
            ssdist = round_to_nths(ssdist, 100)
            saz = round_to_nths(saz, 100)
            if saz < 0:
                saz = saz +360
            lon3 = lons2[i]
            lat3 = lats2[i]
            dist3 = ssdist
            raz = saz
    return dist3

###############################################

### 13 ###

###############################################

## Written GM

# GAVIN'S VERSION
def makePDF3(frame, dep_range, etypes):
    frame = frame.reset_index(drop=True)
    probs = np.zeros(np.size(dep_range))
    for l in range(np.size(dep_range)):
        psum = 0
        o = 0
        test_depth = dep_range[l]
        for type in etypes:
            frame2 = frame[frame.etype == type]
            adepths = frame2['depth'].values
            variance = frame2['unc'].values
            for i in range(len(adepths)):
                dpth, var = adepths[i], variance[i]
                stdv = math.sqrt(var)
                p=(1/(stdv*((2*np.pi)**0.5)))*np.exp(-((test_depth-dpth)**2)/(2*var))
                psum = psum+p
            probs[l] = probs[l] + psum/len(frame2)
            o = o+1
        probs[l] = probs[l]/o
    return probs

###############################################

### 14.5 ###

###############################################

def finiteUnc2(data):
    ''' Arguments:  data - dataframe
        
        Returns:    data - with finite default uncertainties
        '''
    for index, row in data.iterrows():
        unc, lat, lon, etype = row['unc'], row['lat'], row['lon'], row['etype']
        if np.isfinite(unc):
            continue
        else:
            if etype == 'AA':
                row['unc'] = 5
            if etype == 'AS':
                row['unc'] = 2.5
            if etype == 'TO':
                row['unc'] = 40
            if etype == 'EQ':
                row['unc'] = 15
            if etype == 'ER':
                row['unc'] = 10
            if etype == 'RF':
                row['unc'] = 10
            if etype == 'CP':
                row['unc'] = 10
            if etype == 'BA':
                row['unc'] = 1
    return data

def finiteUnc(x):
    ''' Arguments:  x - uncertainty value
        
        Returns:    x - a finite uncertainty value, the same as the argument if initially finite '''
    
    if np.isfinite(x):
        return x
    else:
        return 40.


###############################################

### 14 ###

###############################################

## Written GM

# GINEVRA'S VERSION
def makePDF4(frame, dep_range, etypes, testprint, dpstring):
    
    ''' Arguments:  frame - filtered list of events that will be used to find the local depth
                    dep_range - range of possible events over which to crate a probability distribution
                    etypes - list of event types that are present in this list of events
        
        Returns:    probs - a summed probability distribution matching the size of dep_range    '''
    
    # Preallocating space/resetting DF indices for future processing
    frame = frame.reset_index(drop=True)
    probs = np.zeros(len(dep_range))
    manyAS = False
    N = 1
    
    '''
    try:
        frame.loc[np.isnan(frame.unc), 'unc']=40
    except:
        print ('frame',frame)
        lontest = frame['lon'].values
        lattest = frame['lat'].values
        unctest = frame['unc'].values
        deptest = frame['depth'].values
        print ('lontest,lattest,unctest,deptest,type(unctest)',lontest,lattest,unctest,deptest,type(unctest))
        for i in range(len(frame)):
            print ('lon,lat,dep,unc',lontest[i],lattest[i],deptest[i],unctest[i])
    '''
    
    # Loop through types of events (AS, EQ, ER, BA, etc)
    for itype in etypes:
        
        # Creating arrays/trimming eventlist to only be inclusive of the current type
        pdf_type = np.zeros(len(dep_range))
        if itype == 'EQ':
            frame2 = frame[(frame.etype == 'ER')|(frame.etype == 'EQ')]
        else:
            frame2 = frame[frame.etype == itype]
        adepths = frame2[dpstring].values
        variance = frame2['unc'].values
        
        # Looping through all events of this type
        for i in range(len(adepths)):
            
            # Finding probability for range of depths for each event and adding it to summed PDF
            dpth, var = adepths[i], variance[i]
            stdv = math.sqrt(var)
            pdfi = 1/(stdv * math.sqrt(2*np.pi)) * np.exp(-((dep_range-dpth) ** 2)/(2*var))
            pdf_type = pdf_type + pdfi
        
        # Normalizing summed PDF by dividing by the area under the curve
        area = np.sum(pdf_type)
        pdf_norm = pdf_type/area
        
        # For special cases, event types with low numbers of data points will overpower other event types with
        # lower uncertainties. This addresses that problem by dividing by the degree of overpowering.
        if itype == 'AS' and len(frame2) > 2:
            manyAS = True
            N = len(frame2)
        elif itype == 'BA' and len(frame2) > 2:
            manyAS = True
            N = len(frame2)
        if itype == 'ER' and len(frame2) < 2 and manyAS:
            pdf_norm = pdf_norm/N
        elif itype == 'EQ' and len(frame2) < 2 and manyAS:
            pdf_norm = pdf_norm/N
        probs = probs + pdf_norm
    
    # Normalize the PDF by dividing by the number of event types
    probs = probs / len(etypes)
    return probs

###############################################

### 15 ###

###############################################

## Written GM
## Edited MF 7.25.16
## Edited DEP.8.3.16
## Edited MF 8.9.16
## Edited GLM 11.12.16 (middle cutoff -sdr not -3*sdr)
## Edited GLM.11.16.16 (deep cutoff -4*sdr not 5*sdr)
## Edited GLM 12.7.16 - depth search taper between narrow and broad search ranges
## Edited GLM 12.12.16 - make taper less steep at onset for a more narrow range

def depthRange(loc_depth, sdr, ddr, seismo_thick, elist, slabname, these_parameters, depthwritten):
    
    ''' Arguments:  loc_depth - local depth of Slab 1.0 where available, otherwise average of
                                events within radius of node
                    sdr - range to search around loc_depth in shallow regions (below 150 km)
                    ddr - range to search below loc_depth in deep regions
                    seismo_thick - seismogenic thickness, define where to change search range
                    elist - list of events within lon/lat radius of node
                    slabname - slab ID, i.e. sam, kur, alu etc.
        
        Returns:    elist - filtered dataset with events above or below the depth bounds removed
                    sdepth - shallow bound of events to include in finding depth of slab at this node
                    ddepth - deep bound of events to include in finding deph of slab at this node      '''
    
    dontremove = elist[(elist.etype == 'AA') | (elist.etype == 'BA') | (elist.etype == 'AS') | (elist.etype == 'RF') | (elist.etype == 'CP')]
    elist = elist[(elist.etype != 'AA') & (elist.etype != 'BA') & (elist.etype != 'AS') & (elist.etype != 'RF') & (elist.etype != 'CP')]
    # Initialize search taper arrays
    ds = []
    dd = []
    dp = []
    
    # Define depths breaks, search range, and taper density
    shallow_cutoff = seismo_thick
    if slabname == 'phiz':
        middle_cutoff = 50.0
    else:
        middle_cutoff = 300.0
    smax = float(3*sdr)
    dmax = ddr
    tape_density = 40
    
    # Create taping search radii arrays
    tape = (middle_cutoff-shallow_cutoff)/float(tape_density)
    sdt = (smax-sdr)/float(tape_density)
    ddt = (dmax-sdr)/float(tape_density)
    
    '''
    for k in range(tape_density):
        # Small increase in search range until depth>150 km
        if loc_depth < 150:
            k = k/4.0
        # Slightly larger increase in search range until depth>200 km
        elif loc_depth < 200:
            k = k/2.0
        # Full linear taper between 200 & 300 km
        j = k+1
        ds.append(sdr+sdt*j)
        dd.append(sdr+ddt*j)
        dp.append(shallow_cutoff+tape*j)
    '''
    
    k = np.arange(tape_density+1)
    if loc_depth < 150:
        k = k/4.0
    elif loc_depth < 200:
        k = k/2.0
    j = k+1
    ds = sdr+sdt*j
    dd = sdr+ddt*j
    dp = shallow_cutoff+tape*j
    
    #MF adding depth limits for shallow subduction zones
    if slabname == 'sul': 
        sdepth = loc_depth - sdr
        ddepth = 150 #max
    elif slabname == 'cot': 
        sdepth = loc_depth -sdr
        ddepth = 100 #max 
    else:    
    # Defining conditions based on depth of how far to search around local depth. SUBJECT TO CHANGE
        if loc_depth <= shallow_cutoff:
            sdepth = loc_depth - sdr
            ddepth = loc_depth + sdr
        elif loc_depth <= middle_cutoff:
            for k in range(tape_density-1):
                if loc_depth <= dp[k]:
                    sdr = ds[k]
                    ddr = dd[k]
                    break
            sdepth = loc_depth - sdr
            ddepth = loc_depth + ddr
        else:
            sdepth = loc_depth - 3*sdr
            ddepth = loc_depth + ddr
    
    elist = elist[elist.depth <= ddepth]
    elist = elist[elist.depth >= sdepth]

    elist = pd.concat([elist,dontremove])
    
    return elist, sdepth, ddepth, True

def getangle(a1, b1, c1, a2, b2, c2):

    dot12 = a1*a2 + b1*b2 + c1*c2
    dot11 = a1*a1 + b1*b1 + c1*c1
    dot22 = a2*a2 + b2*b2 + c2*c2
    try:
        lengthV1 = math.sqrt(dot11)
        lengthV2 = math.sqrt(dot22)
    except:
        lengthV1 = np.sqrt(dot11)
        lengthV2 = np.sqrt(dot22)
    try:
        inner = dot12 / (lengthV1 * lengthV2)
        if inner < -0.999:
            inner += 0.001
        elif inner > 0.999:
            inner -= 0.001
        return math.acos(inner)
    except:
        inner = dot12 / (lengthV1 * lengthV2)
        inner[inner < -0.999] += 0.001
        inner[inner > 0.999] -= 0.001

        return np.arccos(inner) # caused runtime

###############################################

### 15.25 ###

###############################################

def perpfilt(rmaxS, rmaxD, rmin, elist, shallow_cutoff):
    #ms = 90.0/(rmaxS - rmin)
    #md = 90.0/(rmaxD - rmin)

    #bs = -1.0 * ms * rmin
    #bd = -1.0 * md * rmin
    
    ms = -90.0/(rmaxS - rmin)
    md = -90.0/(rmaxD - rmin)

    bs = -1.0 * ms * rmaxS
    bd = -1.0 * md * rmaxD
    
    shallowsearch = elist[elist.angle < 90]
    deepsearch = elist[elist.angle >= 90]
    shallowsearch['maxdist'] = (shallowsearch['phi'] - bs)/ms
    deepsearch['maxdist'] = (deepsearch['phi'] - bd)/md
    
    #print 'lat,lon,loc_depth,cstr,cdip,xnode,ynode,znode,rs,rd,shallowsearch',lat,lon,loc_depth,cstr,cdip,xnode,ynode,znode,rs,rd,shallowsearch
    #print 'lat,lon,loc_depth,cstr,cdip,xnode,ynode,znode,rs,rd,deepsearch',lat,lon,loc_depth,cstr,cdip,xnode,ynode,znode,rs,rd,deepsearch

    #shallowsearch = shallowsearch[(shallowsearch['xyzdist'] < shallowsearch['maxdist']) & (shallowsearch['hdist'] >= shallowsearch['vdist'])]#12.28.9-10AM
    shallowsearch = shallowsearch[(shallowsearch['xyzdist'] < shallowsearch['maxdist'])]
    #deepsearch = deepsearch[(deepsearch['xyzdist']<rmin) | ((deepsearch['xyzdist'] < deepsearch['maxdist']) & (deepsearch['hdist'] >= deepsearch['vdist']))] #12.28.9-10AM
    deepsearch = deepsearch[(deepsearch['xyzdist'] < deepsearch['maxdist'])]

    frames = [deepsearch, shallowsearch]
    elist2 = pd.concat(frames)
    elist2 = elist2.reset_index(drop=True)
    return elist2

def dualdepthperp(loc_depth, sdr, ddr, seismo_thick, elist, slabname, cstr, lon, lat, cdip, alen, blen, these_parameters, perpwritten):

    ''' Arguments:  loc_depth - local depth of Slab 1.0 where available, otherwise average of
                                events within radius of node
                    sdr - range to search around loc_depth in shallow regions (below 150 km)
                    ddr - range to search below loc_depth in deep regions
                    seismo_thick - seismogenic thickness, define where to change search range
                    elist - list of events within lon/lat radius of node
                    slabname - slab ID, i.e. sam, kur, alu etc.
        
        Returns:    elist - filtered dataset with events above or below the depth bounds removed
                    sdepth - shallow bound of events to include in finding depth of slab at this node
                    ddepth - deep bound of events to include in finding deph of slab at this node      '''
    
    # Initialize search taper arrays
    ds = []
    dd = []
    dp = []
    
    # Define depths breaks, search range, and taper density
    shallow_cutoff = seismo_thick
    if slabname == 'phiz':
        middle_cutoff = 50.0
    else:
        middle_cutoff = 300.0
    tape1 = middle_cutoff-2*shallow_cutoff
    tape2 = middle_cutoff-3*shallow_cutoff
    ddr = ddr
    smax = float(3*sdr)
    dmax = ddr
    tape_density = 40
    
    # Create taping search radii arrays
    tape = (middle_cutoff-shallow_cutoff)/float(tape_density)
    sdt = (smax-sdr)/float(tape_density)
    ddt = (dmax-sdr)/float(tape_density)
    
    k = np.arange(tape_density+1)
    if loc_depth < 150:
        k = k/4.0
    elif loc_depth < 200:
        k = k/2.0
    j = k+1
    ds = sdr+sdt*j
    dd = sdr+ddt*j
    dp = shallow_cutoff+tape*j
    
    # Define inboard/outboard searching distances
    if loc_depth <= shallow_cutoff:
        rs = sdr
        rd = sdr
    
    elif loc_depth <= middle_cutoff:
        rs = smax
        rd = dmax
        # Radii are tapered from sdr to smax and dmax where 70 < depth < 300
        for k in range(tape_density-1):
            if loc_depth <= dp[k]:
                rs = ds[k]
                rd = dd[k]
                break
    else:
        rs = smax
        rd = ddr

    radstr = cstr * math.pi/180.0
    raddip = cdip * math.pi/180.0
    xs = math.cos(radstr)
    ys = math.sin(radstr)
    hd = math.cos(raddip)
    zd = math.sin(raddip)

    zdist = elist['depth'].values - loc_depth
    elist['zdist'] = abs(zdist)
    
    elist['cosdistance'], cosangles = npcosine(lon, lat, elist['lon'].values, elist['lat'].values)
    cosangles -= 180
    cosangles[cosangles<0]+=360
    elist['outboard'] = np.logical_not(npoutboard(cstr, cosangles)) # Will need to fix when fix outboard function
    #elist['cosdistance'][elist.outboard == True] *= -1
    
    cosangles[cosangles <= 180.0] += 360.0
    cosangles -= 180.0
    elist['anglediff'] = abs(cstr - cosangles)
    elist['phiS'] = abs(elist['anglediff']-90)
    elist['cosdistance2'] = elist['cosdistance'].values * np.sin(np.radians(elist['phiS'].values))
    elist['cosdistance'] = elist['cosdistance'].values * np.cos(np.radians(elist['phiS'].values))
    elist['cosdistance'][(elist.outboard == True) & (elist.cosdistance > 0)] *= -1
    elist['cosdistance'][(elist.outboard == False) & (elist.cosdistance < 0)] *= -1
    
    elist['alldist'] = np.sqrt(elist['zdist'].values * elist['zdist'].values + elist['cosdistance'].values * elist['cosdistance'].values)
    defkeep = elist[elist.alldist<blen]
    dangle = getangle(hd, 0.0, zd, elist['cosdistance'].values, np.zeros(len(zdist)), zdist)
    elist['dangle'] = dangle * 180/math.pi
    #elist['phiS'] = abs(elist['sangle'] - 90.0)
    elist['phiD'] = abs(elist['dangle'] - 90.0)
    
    rminD = blen+0.1
    rminS = blen-0.1
    rmaxSS = alen
    rmaxSD = rs
    rmaxDD = rd
    
    if abs(rmaxSD-rminD) < 1:
        rminD -= 2
    if abs(rmaxDD-rminD) < 1:
        rminD -= 2
    if abs(rmaxSS-rminS) < 1:
        rminS -= 2

    sdepth = loc_depth-rmaxSD
    ddepth = loc_depth+rmaxDD

    phiSmax = math.atan2(rmaxSS, rminS)*180.0/math.pi
    phiDSmax = math.atan2(rmaxSD, rminD)*180.0/math.pi-90
    phiDDmax = math.atan2(rmaxDD, rminD)*180.0/math.pi-90
    
    elist['maxdepth'] = blen/np.sin(np.radians(elist['phiD'].values))

    elist.loc[(elist.phiD < phiDDmax) & (elist.outboard == True), 'maxdepth'] = rmaxDD
    elist.loc[(elist.phiD < phiDSmax) & (elist.outboard == False), 'maxdepth'] = rmaxSD
    elist.loc[(elist.maxdepth > rmaxDD) & ((elist.outboard == True) & (elist.depth > loc_depth)), 'maxdepth'] = rmaxDD
    elist.loc[(elist.maxdepth > rmaxSD) & ((elist.outboard == False) | (elist.depth < loc_depth)), 'maxdepth'] = rmaxSD
    elist = elist[(elist['alldist'] < elist['maxdepth'])]
    
    elist2 = pd.concat([defkeep, elist])
    
    elist2 = elist2.drop_duplicates(['ID'])
    elist2 = elist2[elist2.cosdistance2<alen]
    '''
    if len(elist2)>1 and cdip > 50:

        elist2.to_csv('Output/perptest/%s_%.4f_%.4f_3.csv'%(slabname,lon,lat),header=True,index=False,float_format='%0.2f',na_rep = float('nan'),sep='\t')

        print ('elist3',lon,lat,loc_depth,cstr,cdip,rmaxDD,rmaxSD)
        elist4 = elist2[elist2.outboard == True]
        elist5 = elist2[elist2.outboard != True]
        xdist4 = (elist4['lon'].values - lon)
        ydist4 = (elist4['lat'].values - lat)
        xdist5 = (elist5['lon'].values - lon)
        ydist5 = (elist5['lat'].values - lat)
        
        fig = plt.figure()
        print ('PERPENDICULAR1',lon,lat)
        thispoint = plt.plot([0],[0],'ro',label='Node Location')
        inboardpts = plt.plot(xdist5,ydist5,'yo',label='Inboard')
        outboardpts = plt.plot(xdist4,ydist4,'bo', label='Outboard')
        strike = plt.plot([-ys,ys],[-xs,xs],'k-',label='Strike')
        plt.xlabel('longitude difference west <- -> east')
        plt.ylabel('longitude difference ^ north')
        plt.grid()
        plt.axis('equal')
        title = 'Lat: %.2f, Lon: %.2f, Strike: %.2f, Dip: %.2f' % (lat,lon,cstr,cdip)
        plt.title(title)
        lontit = lon*100
        lattit = lat*100
        plt.legend(loc='best')
        figtitle = 'Output/perptest/pdf%.2f_%.2f_p.png' % (lon,lat)
        #fig.savefig(figtitle)
        plt.close()
        
        zdist4 = (elist4['depth']-loc_depth)
        zdist5 = (elist5['depth']-loc_depth)
        fig = plt.figure()
        thispoint = plt.plot([0],[0],'ro',label='Node Location')
        outboardpts = plt.plot(elist4['cosdistance'].values,zdist4,'bo', label='Outboard')
        inboardpts = plt.plot(elist5['cosdistance'].values,zdist5,'yo',label='Inboard')
        dip = plt.plot([-hd*50,hd*100],[-zd*50,zd*100],'g-',label='Dip')
        plt.xlabel('horizontal distance outboard <- -> inboard')
        plt.ylabel('vertical distance deeper <- -> shallower')
        ax = plt.gca()
        plt.axis('equal')
        ax.invert_yaxis()
        plt.grid()
        title = 'Lat: %.2f, Lon: %.2f, Strike: %.2f, Dip: %.2f, Origin Depth: %.2f' % (lat,lon,cstr,cdip,loc_depth)
        plt.title(title)
        lontit = lon*100
        lattit = lat*100
        plt.legend(loc='best')
        figtitle = 'Output/perptest/pdf%.2f_%.2f_c.png' % (lon,lat)
        #fig.savefig(figtitle)
        plt.close()
    '''
    elist2 = elist2[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
    return elist2, rs, rd, True

###############################################

### 17 ###

###############################################

## Written GM
## Edited MF.8.2.16
## Edited DEP.8.13.16
## Edited DEP.8.16.16

def getTrenchStrike(TR_data, lat, lon, tooFardist, testprint):
    
    ''' Arguments:  TR_data - list of local trench points (includes lat, lon, strike)
                    lat - latitude of this grid node
                    lon - longitude of this grid node
                    minlat, maxlat - latitude bounds of slab model
                    minlon, maxlon - longitude bounds of slab model
                    ID - output ID of the current node (for recording purposes)
        
        Returns:    az - strike of closest trench point
                    minang - angle between nearest trench point and this node (clockwise from 0 degrees)
                    mindist - distance between node and closest trench point
                    tooFar - a boolean value that states whether or not the distance is greater than a defined value
        '''
    
    # Setting minimum values and longitude bounds relative to grid node
    mindist = 9999
    minang = 9999
    tooFar = True
    
    TR_data['cosdistance'], cosangles = npcosine(lon, lat, TR_data['lon'].values, TR_data['lat'].values)
    cosangles -= 180
    cosangles[cosangles<0]+=360
    neardist = TR_data['cosdistance'].min()
    #neardist = TR_data['hdist'].min()
    neartrench = TR_data[TR_data.cosdistance == neardist]
    
    cstrB = neartrench['az'].values[0]
    tlon = neartrench['lon'].values[0]
    tlat = neartrench['lat'].values[0]
    tdepth = neartrench['depth'].values[0]
    mindist, minang, lon1, lat1 = cosine(tlon, tlat, lon, lat)
    
    cosangles[cosangles <= 90.0] += 360.0
    TR_data['cosangles'] = cosangles
    cosangles -= 90.0
    TR_data['cosangles1'] = cosangles
    TR_data['anglediff'] = np.abs(TR_data['az'].values - cosangles)

    perptrenches = TR_data[TR_data.anglediff < 5]
    perpdist = perptrenches['cosdistance'].min()
    perptrench = perptrenches[perptrenches.cosdistance == perpdist]

    try:
        cstrP = perptrench['az'].values[0]
        tlonP = perptrench['lon'].values[0]
        tlatP = perptrench['lat'].values[0]
    except:
        cstrP = cstrB
        tlonP = tlon
        tlatP = tlat
        if testprint:
            print ('no perpendicular trench, mindist,minang,lon,tlon,lat,tlat',mindist,minang,lon,tlon,lat,tlat)
        #print 'no perpendicular trench, mindist,minang,lon,tlon,lat,tlat',mindist,minang,lon,tlon,lat,tlat

    perpdist, perpang, lon2, lat2 = cosine(tlonP, tlatP, lon, lat)

    if testprint:
        print ('trenchdata lon, lat, perpdist, perpang, mindist, minang, tlonP, tlatP, tlon, tlat',lon, lat, perpdist, perpang, mindist, minang, tlonP, tlatP, tlon, tlat)
    if mindist < tooFardist:
        tooFar = False
    if (mindist < 150 and perpdist > 350) or perpdist > 700:
        #print 'used near not perp: lon,lat,mindist,minang,perpdist,perpang,cstrB,cstrP: ',lon,lat,mindist,minang,perpdist,perpang,cstrB,cstrP
        perpdist = mindist
        perpang = minang

    return cstrP, cstrB, minang, perpdist, tooFar, tlon, tlat, tdepth

###############################################

### 18 ###

###############################################

## Written GM

def make_indices(N, P):
    ''' Arguments:  N - number of latitude rows in the slab model
                    P - number of cores used for processing
        
        Returns:    index_array - array of indeces associated with respective model rows
                                    notifies each core of the part of the model it needs to complete '''
    index_array = np.zeros(P+1)
    for i in range(P+1):
        index_array[i] = i*np.floor(N/P)
    for i in range(np.mod(N, P)):
        index_array[i+1:] += 1
    return index_array

###############################################

### 19 ###

###############################################

## Written GM

def isoutboard(az, ang):
    ''' Arguments:  az - strike of the trench
                    ang - angle between trench and point (clockwise from 0)
        
        Returns:    True - if point is outboard of the trench
                    False - if point is inboard the trench       '''
    
    # Calculating difference between the strike of the trench and the angle between trench and point
    azang = az - ang
    
    # Finding whether or not the point is outboard - conditions change for different cases
    if (az >= 180 and az <= 360):
        if (azang >=0 and azang <= 180):
            return True
        else:
            return False
    elif (az >= 0 and az < 180):
        if (azang >= 0):
            return True
        elif (azang <= -180):
            return True
        else:
            return False
    else:
        return False


def npoutboard(az, ang):

    # Written GLM 4.28.17

    ''' Arguments:  az - azimuth of reference point (float)[deg]
                    ang - array of angles, clockwise from 0, from reference
                            point to other points (arr of floats)[deg]
                            
        Returns:    out - boolean array, True where points are outboard of
                            reference, false where they are inboard. '''

    azang = az-ang
    out = np.ones(len(ang), dtype=bool)

    if az >= 180 and az <= 360:
        out[(azang >= 0)&(azang <= 180)] = True
        out[(azang < 0)|(azang > 180)] = False

    elif az >= 0 and az < 180:
        out[(azang >= 0)|(azang <= -180)] = True
        out[(azang > -180)&(azang < 0)] = False

    else:
        out[azang >= -360] = False

    return out

###############################################

### 20 ###

###############################################

# Written MF 7.19.16
# Edited DEP.8.4.16

def slabpolygon(slabname, slabfile):
    
    '''
        inputting the slabname (3 character code) and slabfile will return the polygon boundaries
    '''
    
    #load file with slab polygon boundaries
    filerows = []
    with open(slabfile) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            filerows.append(row)
    csvfile.close()

    #iterate through list to match the slabname and retrieve coordinates
    slabbounds = []
    for i in range(len(filerows)):
        if slabname == filerows[i][0]:
            slabbounds = filerows[i][1:]
            slabbounds.append(slabbounds)

    return slabbounds

###############################################

### 21 ###

###############################################

# Written MF 7.18.16
# Edited DEP.8.4.16

def determine_polygon_extrema(slabname, slabfile):
    
    '''
        inputs: slabname to be referenced against stored slab coordinates in slabfile
        outputs: the maximum and minimum latitude and longitude values for the input slab
    '''
    
    #calls slabpolygon function to get bounds for this slab region
    slabbounds = slabpolygon(slabname, slabfile)
    
    #slabbbounds come in lon1,lat1,lon2,lat2... format
    #even numbers are then longitudes while odds are latitudes
    coords = np.size(slabbounds)
    
    #simple even/odd function
    def is_odd(num):
        return num & 0x1
    
    lons = []
    lats = []
    
    for i in range(coords-1):
        val = float(slabbounds[i])
        if is_odd(i):
            lats.append(val)
        else:
            lons.append(val)

    x1 = int(min(lons))
    x2 = int(max(lons))
    y1 = int(min(lats))
    y2 = int(max(lats))
    # maybe need for alu
    #if x1<0:
    #    x1 += 360
    #if x2<0:
    #    x2 += 360
    
    return x1, x2, y1, y2



def create_grid_nodes3(grid, lonmin, lonmax, latmin, latmax):

    #define grid of searches (representative of lower left corner)
    xall = np.arange(math.floor(lonmin)-2,math.ceil(lonmax)+2,grid)
    yall = np.arange(math.floor(latmin)-2,math.ceil(latmax)+2,grid)
    lons1,lats1 = np.meshgrid(xall,yall)

    #flatten into list of lower left corners
    lllons = lons1.flatten()
    lllats = lats1.flatten()

    #combine into one array (lonmin,lonmax,latmin,latmax)
    bounds = np.zeros((len(lllons),2))
    bounds[:,0] = lllons
    bounds[:,1] = lllats

    return bounds

###############################################

### 23 ###

###############################################

def createGridInPolygon2(nodes, slabname, slabfile):
    
    #acquire slabbounds
    slabbounds = slabpolygon(slabname, slabfile)
    
    coords = np.size(slabbounds)
    
    #simple even/odd function
    def is_odd(num):
        return num & 0x1
    
    lons = []
    lats = []
    
    for i in range(coords):
        val = slabbounds[i]
        if is_odd(i):
            lats.append(val)
        else:
            lons.append(val)

    #create tuple of locations (with zip) to use in contains_points
    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(nodes[:])
    mask1 = np.zeros(len(temp),)*np.nan
    mask1[temp] = 1
    keepers = []
    for i in range(len(nodes)):
        points_in_poly = np.dot(mask1[i], nodes[i])
        if i > 0:
            keepers = np.vstack((keepers, points_in_poly))
        else:
            keepers = points_in_poly

    values = []
    for i in range(len(keepers)):
        if np.isnan(keepers[i][0]) == False:
            values.append(keepers[i])

    valid_points = np.array(values)
    return valid_points

###############################################

### 24 ###

###############################################

def getDFinMask(datadf, maskdf):

    maskdf.loc[maskdf.lon < 0, 'lon'] += 360
    datadf.loc[datadf.lon < 0, 'lon'] += 360

    lons = maskdf['lon'].values*1.0
    lats = maskdf['lat'].values*1.0

    dlons = datadf['lon'].values*1.0
    dlats = datadf['lat'].values*1.0
    nodes = np.zeros((len(dlons),2))
    nodes[:,0] = dlons
    nodes[:,1] = dlats

    #create tuple of locations (with zip) to use in contains_points
    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(nodes[:])
    mask1 = np.zeros(len(temp),)*np.nan
    mask1[temp] = 1
    keepers = []
    for i in range(len(nodes)):
        points_in_poly = np.dot(mask1[i], nodes[i])
        if i > 0:
            keepers = np.vstack((keepers, points_in_poly))
        else:
            keepers = points_in_poly
    values = []
    for i in range(len(keepers)):
        if np.isnan(keepers[i][0]) == False:
            values.append(keepers[i])

    valid_points = np.array(values)

    newdf = pd.DataFrame({'lon':valid_points[:,0],'lat':valid_points[:,1],'testcol':1})
    newdf = pd.merge(newdf, datadf, left_on = ['lon','lat'], right_on = ['lon','lat'])
    cols = datadf.columns
    newdf = newdf[newdf.testcol == 1]
    newdf = newdf[cols]

    return newdf

###############################################

### 24 ###

###############################################

def getDataInPolygon(slabname, data, slabfile):

# Written MF 7.20.16
# MF edited 8.5.16 

    ''' creates a grid of 1 or nan based on if they are within a clipping mask or not. DEP.6.29.16 '''
    ''' modified to fit this script by MAF 7/18/16 '''

    ### Summary:
    #This is very similar to nodesInPolygon except that it takes in input data without needing
    #to define a regular grid. The data are then formatted for the point in polygon search,
    #and only those data which are within th slabname (polygon) are kept.  

    ### Input:
    # slabname: a 3 digit character code identifying a slab region 
    #data: the input data which may or may not be within the polygon

    ### Output:
    #contained_data: an array of coordinate pairs (lon,lat) that reside within the polygon region

    #acquire slabbounds 
    slabbounds = slabpolygon(slabname, slabfile)

    #slabbbounds come in lon1,lat1,lon2,lat2... format
    #even numbers are then longitudes while odds are latitudes
    coords = np.size(slabbounds)

    #simple even/odd function 
    def is_odd(num):
        return num & 0x1

    lons = []
    lats = []

    for i in range(coords):
        val = slabbounds[i]
        if is_odd(i):
            lats.append(val)
        else:
            lons.append(val)

    #create tuple of locations (with zip) to use in contains_points
    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(data[:])
    mask1 = np.zeros(len(temp),)*np.nan
    mask1[temp] = 1
    keepers = []
    for i in range(len(data)):
        points_in_poly = np.dot(mask1[i], data[i])
        if i > 0:
            keepers = np.vstack((keepers, points_in_poly))
        else:
            keepers = points_in_poly

    values = []
    for i in range(len(keepers)):
        if np.isnan(keepers[i][0]) == False:
            values.append(keepers[i])

    contained_data = np.array(values)

    return contained_data 

###############################################

### 25 ###

###############################################

## Written DEP

# epCalc is used to calculate an array of locations and distances from the slab surface by calculating vector endpoints within the earth
# Edited DEP.7.1.16 to accommodate changed mkSDgrd

def epCalc(lon, lat, dep, dip, strike, posmag, negmag, step):
    EPs = np.zeros(((posmag-negmag)/step, 4))
    
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

def mkSDgrd(Slabgrid):
    
    # Get depth grid boundaries, size, and spacing
    gdict = Slabgrid.getGeoDict().copy()
    nx = gdict.nx
    ny = gdict.ny
    dx = gdict.dx
    dy = gdict.dy
    xmin = gdict.xmin
    xmax = gdict.xmax
    ymin = gdict.ymin
    ymax = gdict.ymax
    
    # Create lat/lon grid from parameters above
    dlats = np.linspace(ymax, ymin, ny)
    dlons = np.linspace(xmin, xmax, nx)
    dlons[dlons<0] += 360
    
    # the guides always have negative depths, account for this
    depthgrid = Slabgrid.getData().copy()*-1.0
    
    # initialize unfactored spacing in km between degrees
    alldy = dy * 111.19
    alldx = dx * 111.19
    Xgrad, Ygrad = [],[]
    
    # loop through list of latitudes, calculate gradient at each
    for i in range(1, ny - 1):
        
        # get row at this lat, also one row N and one row S
        thisgrid = depthgrid[i - 1:i + 2,:]
        
        # calculate longitude distance at this lat
        thisy = math.radians(abs(dlats[i]))
        thisdx = alldx * math.cos(thisy)
        
        # calculate gradient at these 3 rows
        gradyi, gradxi = np.gradient(thisgrid, alldy, thisdx)
        
        # if it is the first calculation, add first two rows of this grad
        if len(Xgrad) < 1:
            Xgrad = gradxi[0:2, :]
            Ygrad = gradyi[0:2, :]
        
        # otherwise, only add the middle row of this grad
        else:
            Xgrad = np.vstack((Xgrad, gradxi[1, :]))
            Ygrad = np.vstack((Ygrad, gradyi[1, :]))

    # add the last row of the last grad to end of grid
    Xgrad = np.vstack((Xgrad, gradxi[2, :]))
    Ygrad = np.vstack((Ygrad, gradyi[2, :]))
    
    # Get gradient magnitude
    Maggrid = np.sqrt((Ygrad**2)+(Xgrad**2))
    
    # Define a grid file that is the direction perpendicular to the max gradient
    Dirgrid = np.degrees(np.arctan2(Ygrad, Xgrad))
    Dirgrid = np.where(Dirgrid < 0, Dirgrid+360, Dirgrid)

    # Assign strike and dip arrays to grids with same dimensions as depth grid
    Dipgrid = gmt.GMTGrid(np.degrees(np.arctan2(Maggrid, 1)), Slabgrid.getGeoDict().copy())
    Strikegrid = gmt.GMTGrid(Dirgrid, Slabgrid.getGeoDict().copy())

    return Strikegrid, Dipgrid

def mkSDgrd_old(Slabgrid):
    
    # Get depth grid parameters
    gdict = Slabgrid.getGeoDict().copy()
    
    # Define dx/dy as the lon/lat increments and convert to km
    dx = gdict.dx * 111.19
    dy = gdict.dy * 111.19
    
    # Define gradient of depth grid in y and x directions
    Ygrad, Xgrad = np.gradient(Slabgrid.getData().copy()*-1.0, dx, dy)
    
    # Get gradient magnitude
    Maggrid = np.sqrt((Ygrad**2)+(Xgrad**2))
    
    # Define a grid file that is the direction perpendicular to the max gradient
    Dirgrid = np.degrees(np.arctan2(Ygrad, Xgrad))
    Dirgrid = np.where(Dirgrid < 0, Dirgrid+360, Dirgrid)

    # Assign strike and dip arrays to grids with same dimensions as depth grid
    Dipgrid = gmt.GMTGrid(np.degrees(np.arctan2(Maggrid, 1)), Slabgrid.getGeoDict().copy())
    Strikegrid = gmt.GMTGrid(Dirgrid, Slabgrid.getGeoDict().copy())
        
    return Strikegrid, Dipgrid

def mkSDgrddata(xi, zi, flipornot):
    
    # get dx, dy, and list of lats from zi coordinates (listed in xi)
    xpts, ypts = xi[:, 0], xi[:, 1]
    xpts.shape = zi.shape
    ypts.shape = zi.shape
    dlats = ypts[:, 0]
    dlons = xpts[0, :]
    ny = len(dlats)
    nx = len(dlons)
    dy = abs(dlats[1] - dlats[0])
    dx = abs(dlons[1] - dlons[0])
    
    # flip array over if needed
    if flipornot == 'flip':
        depthgrid = np.flipud(zi)
    else:
        depthgrid = np.copy(zi)

    # initialize grid spacing in km
    alldy = dy * 111.19
    alldx = dx * 111.19
    Xgrad, Ygrad = [],[]

    # loop through lats and get gradient, use different lon spacing for each lat
    for i in range(1, ny - 1):
        thisgrid = depthgrid[i - 1:i + 2,:]
        thisy = math.radians(abs(dlats[i]))
        thisdx = alldx * math.cos(thisy)
        gradyi, gradxi = np.gradient(thisgrid, alldy, thisdx)
        
        # add first two lines to gradient if first loop
        if len(Xgrad) < 1:
            Xgrad = gradxi[0:2, :]
            Ygrad = gradyi[0:2, :]
        
        # otherwise, add just this row to the gradient array
        else:
            Xgrad = np.vstack((Xgrad, gradxi[1, :]))
            Ygrad = np.vstack((Ygrad, gradyi[1, :]))

    # add the last row to the gradient array
    Xgrad = np.vstack((Xgrad, gradxi[2, :]))
    Ygrad = np.vstack((Ygrad, gradyi[2, :]))

    # Get gradient magnitude
    Maggrid = np.sqrt((Ygrad**2)+(Xgrad**2))

    # Define a grid file that is the direction perpendicular to the max gradient
    Strikegrid = np.degrees(np.arctan2(Ygrad, Xgrad))
    Strikegrid = np.where(Strikegrid < 0, Strikegrid+360, Strikegrid)

    # Assign strike and dip arrays to grids with same dimensions as depth grid
    Dipgrid = np.degrees(np.arctan2(Maggrid, 1))

    # flip grids upside down if needed
    if flipornot == 'flip':
        Strikegrid = np.flipud(Strikegrid)
        Dipgrid = np.flipud(Dipgrid)

    return Strikegrid, Dipgrid

def mkSDgrddata_old(xi, zi, flipornot):
    
    dx = abs(xi[1, 0]-xi[0, 0])*111.19
    dy = dx

    if flipornot == 'flip':

        Xgrad, Ygrad = np.gradient(np.flipud(zi), dx, dy)  ## Define grids of the gradients in the x and y directions
        Maggrid = np.sqrt((Xgrad**2)+(Ygrad**2))  ## Define a grid file that is the gradient magnitude at any given node
        Strikegrid = np.degrees(np.arctan2(Xgrad, Ygrad))  ## Define a grid file that is the direction perpendicular to the max gradient - in degrees clockwise from N (0 - 360)
        Strikegrid[(Strikegrid<0)&(np.isfinite(Strikegrid))] += 360
        #Strikegrid = np.where(Dirgrid < 0, Dirgrid+360, Dirgrid) # caused runtime

        Dipgrid = np.degrees(np.arctan2(Maggrid, 1))
    
        Strikegrid = np.flipud(Strikegrid)
        Dipgrid = np.flipud(Dipgrid)
    
    else:
    
        Xgrad, Ygrad = np.gradient(zi, dx, dy)  ## Define grids of the gradients in the x and y directions
        Maggrid = np.sqrt((Xgrad**2)+(Ygrad**2))  ## Define a grid file that is the gradient magnitude at any given node
        Strikegrid = np.degrees(np.arctan2(Xgrad, Ygrad))  ## Define a grid file that is the direction perpendicular to the max gradient - in degrees clockwise from N (0 - 360)
        Strikegrid[(Strikegrid<0)&(np.isfinite(Strikegrid))] += 360
        #Strikegrid = np.where(Dirgrid < 0, Dirgrid+360, Dirgrid) # caused runtime

        Dipgrid = np.degrees(np.arctan2(Maggrid, 1))

    return Strikegrid, Dipgrid

###############################################

### 30 ###

###############################################

## Written DEP

# point_in_poly was written by some guy on the internet at http://geospatialpython.com/2011/01/point-in-polygon.html.  Its purpose is to identify whether or not a point is located within a provided polygon.  I am using it to create an array of 1s and 0s based on a clipping mask.  Takes in lists for x and y.

def point_in_poly(x_all, y_all, poly):
    
    ans = np.zeros((len(x_all), 1))
    for j in range(len(ans)):
        x = x_all[j]
        y = y_all[j]
        n = len(poly)
        
        p1x, p1y = poly[0]
        cnt = 0
        for i in range(n+1):
            inside = False
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                            if p1x == p2x or x <= xints:
                                inside = not inside
                                p1x, p1y = p2x, p2y
        if inside == True:
            cnt = cnt + 1
        if (cnt % 2) == 0:
            ans[j] = np.nan
        else:
            ans[j] = 1
                
    return ans

###############################################

### 36 ###

###############################################

## Written  DEP.7.7.16

# pointShift is essentially epCalc, but for a single point.  It is used to calculate the endpoint of a vector within the earth given a local lat/lon/dep, strike/dip, and distance.

def pointShift(lon, lat, dep, dip, strike, mag):
    
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
    
    # Define translation vector in spherical coordinates
    
    trad = mag
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
    
    return eptlon, eptlat, eptdep

###############################################

### 37 ###

###############################################

## DEP.8.3.16
## Edited DEP.8.5.16 type
## Edited GLM 11.14.16 re-indented

# findLocDep estimates slab depth.  If Slab1.0 exists and it's deep, removes shallow earthquakes as well.

def findLocDep(slab1, tooFar, elist, seismo_thick, testprint, balist, out, slab, lon, lat):

    if np.isfinite(slab1):  # Slab1.0 exists
        loc_depth = slab1
        if loc_depth > seismo_thick + 30 and slab != 'hal' and slab != 'sol':  # If deep, remove shallow.  I removed the "toofar" part
            elist = elist[elist.depth > seismo_thick]

    elif len(balist) > 0 and out:
        loc_depth =  np.mean(balist['depth'].values)

    else: # No Slab1.0
        #if tooFar:
        #    rem_shallow = elist[elist.depth > seismo_thick]
        #    depths_in_circle = rem_shallow['depth'].values
        #else:
        depths_in_circle = elist['depth'].values
        loc_depth = np.mean(depths_in_circle)

    if loc_depth < seismo_thick:
        elist = elist[(elist.depth<seismo_thick)|((elist.etype != 'EQ')&(elist.etype != 'ER'))]

    return loc_depth, elist

###############################################
                          
### 38 ###
                          
###############################################
                          
## DEP.8.4.16

# ellipseFilt filters the data by forming an elipse around the node

def ellipseFilt(elist, lat, lon, alen, blen, cstr, mdist):
    
    if len(elist)<1:
        return elist
    else:
        if mdist > 0:
            rlon, rlat = heading(lon, lat, alen, cstr)
        else:
            rlat, rlon = lat, lon

        trimmed = getEventsInEllipse(rlon, rlat, cstr, alen, blen, elist, lon, lat)

        return trimmed

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r


def trimByTrench(trimmed, outside, AA_data, lat, lon, maxID, size, TR_data, strike, mindist, testprint, slab):

    trimmed = trimmed[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
    if outside or mindist < 2:  # If outside trench, use only bathymetry and active source
        trimmed[(trimmed.etype =='BA') | (trimmed.etype =='AS')] #GLM 11.14.16
    else:
        trimmed = trimmed[trimmed.etype != 'BA']
        if size == 0 and slab != 'helz':
            AA_data['diffdist'] = np.abs(AA_data['dist'].values - mindist)
            mindiff = AA_data['diffdist'].min()
            if mindiff < 0.1:
                thisAA = AA_data[AA_data.diffdist == mindiff]
                locAA = thisAA['depth'].values[0]
                trimmed.loc[len(trimmed)+1] = ([lat, lon, locAA, 5.0, str('AA'), (maxID+1), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.0])
                maxID += 1
    return trimmed, maxID

def trimByTrench_alu(trimmed, outside, AA_data, lat, lon, maxID, size, TR_data, strike, mindist, testprint, slab):

    trimmed = trimmed[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
    if outside or mindist < 2:  # If outside trench, use only bathymetry and active source
        if slab == 'him':
            if mindist > 2:
                #trimmed = trimmed[trimmed.etype != 'CP']
                closetrench = TR_data[TR_data.cosdistance == TR_data['cosdistance'].min()]
                wdepth = closetrench['depth'].values[0]
                trimmed.loc[len(trimmed)+1] = ([lat, lon, wdepth, 5.0, str('AA'), (maxID+1), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.0])
                maxID += 1
            else:
                trimmed[(trimmed.etype =='BA') | (trimmed.etype =='AS')] #GLM 11.14.16
        else:
            trimmed[(trimmed.etype =='BA') | (trimmed.etype =='AS')] #GLM 11.14.16
    else:
        trimmed = trimmed[trimmed.etype != 'BA']
        if size == 0: # GLM 11.18.2016 filters out in some cases where it shouldn't
            AA_data['diffdist'] = np.abs(AA_data['dist'].values - mindist)
            AA_data['sdist'] = gps2dist_azimuth(lat, lon, AA_data['avlat'], AA_data['avlon'])[0]/1000.0
            if lon > AA_data['avlon'].max() or lon < AA_data['avlon'].min():
                thisAA = AA_data[AA_data.sdist == AA_data['sdist'].min()]
                thisAA = thisAA[thisAA.diffdist < 0.2]
            else:
                nearAA = AA_data[AA_data.sdist == AA_data['sdist'].min()]
                thisAA = nearAA[nearAA.diffdist < 0.2]
                #print 'nearAA!',lon,lat,nearAA,thisAA
                if len(thisAA) < 1 and ((lon > 200 and mindist < 150) or (lon < 200 and mindist < 80)):
                    nearAA = AA_data[AA_data.diffdist < 0.2]
                    thisAA = nearAA[nearAA.sdist<500]
                    #print 'toofarbut',lon,lat,thisAA
            if len(thisAA)>0:
                wdepth = 0.0
                thisAA['weights'] = thisAA['sdist'].values/thisAA['sdist'].sum()
                thisAAw = thisAA.sort_values(by=['sdist']) # is deprecated, use sort_values(by=.....)
                thisAAi = thisAA.sort_values(by=['sdist'], ascending=False) # is deprecated, use sort_values(by=.....)
                weights = thisAAi['weights'].values
                depths = thisAAw['depth'].values
                for d in range(len(thisAA)):
                    wdepth += weights[d]*depths[d]
                #print 'lon,lat,wdepth AADISTANCES',lon,lat,wdepth,thisAAw,thisAAi
                trimmed.loc[len(trimmed)+1] = ([lat, lon, wdepth, 5.0, str('AA'), (maxID+1), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.0])
                maxID += 1
                #print 'trimmed!',trimmed
    #print '555',lon,lat,trimmed
    return trimmed, maxID

###############################################
                          
### 40 ###
                          
###############################################
                          
## DEP.8.4.16
## DEP.8.5.16 edited typo
                          
# addData adds the data used for a given node to the full used data list
                          
def addData(trimmed, lat, lon):
    IDs = trimmed['ID'].values
    LONlist = np.ones(len(IDs))
    LATlist = np.ones(len(IDs))
    LONlist = LONlist * lon
    LATlist = LATlist * lat
    used = trimmed[trimmed['ID'].isin(IDs)]

    return used, IDs, LONlist, LATlist

###############################################

### 41 ###

###############################################

def rectangleIntersectsPolygon(x1, x2, y1, y2, slabfile):

# Written MF 8.4.16
# MF edited 8.5.16
# Edited GLM 11.21.16 - poly = Polygon(poly)

    def is_odd(num):
        return num & 0x1

    #create polygon from input rectangle
    rect = Polygon([(x1, y2), (x2, y2), (x2, y1), (x1, y1)])

    #read in slab boundaries 
    filerows = []
    with open(slabfile) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            filerows.append(row)
    csvfile.close()

    #loop through the slabnames and slabboundaries by row to define each slab polygon
    #then verify whether the input rectangle overlaps any of the defined slabs
    slabbounds = []
    slabname = []
    slab = []
    for i in range(len(filerows)-1):
        lats =[]
        lons = []
        slabname = filerows[i][0]
        slabbounds = filerows[i][1:]
        slabbounds.append(slabbounds)
        for j in range(1, (len(filerows[i][:]))):
            val = float(filerows[i][j])
            if is_odd(j):
                lons.append(val)
            else:
                lats.append(val)
        poly = list(zip(lons, lats))
        poly = Polygon(poly) # GLM 11.21.16
        if rect.overlaps(poly):
            slab.append(slabname)
        else:
            continue

    #if the input rectangle does not overlap with just one slab, let the user know
    if len(slab) == 0:
        print('The input boundaries do not overlap any slabs. Please try again.')
    elif len(slab) > 1:
        response = input('You have selected multiple slabs. Which slab would you like to model?: '+str(slab)+' Please enter a string: ')
        slab = response

    return slab

def noDataNeedAA(trimmed, cstr, minang, AA_data, lat, lon, maxID, TR_data, mindist, testprint, sdr, ddr, seismo_thick, slab, these_parameters, depthwritten, perpwritten, trenchlon, trenchlat, AARF, loc_depth):
    trimmed = trimmed[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
    if len(trimmed)>0:
        elistASe = trimmed[trimmed.etype == 'AS']
        elistRFe = trimmed[trimmed.etype == 'RF']
        elistCPe = trimmed[trimmed.etype == 'CP']
        assize = elistASe.size
        rfsize = elistRFe.size
        cpsize = elistCPe.size
    else:
        assize = trimmed.size
        rfsize = trimmed.size
        cpsize = trimmed.size
    test = True
    oceanside = isoutboard(cstr, minang)  # See which side of trench node is on
    length1 = len(trimmed)
    if rfsize < 1:
        trimmed, maxID = getextraRF(trimmed, slab, cstr, mindist, trenchlon, trenchlat, AARF, maxID, lon, lat)
    length2 = len(trimmed)
    if (slab == 'alu' or slab == 'him') and length1==length2 and mindist <= AA_data['dist'].max():
        trimmed, maxID = trimByTrench_alu(trimmed, oceanside, AA_data, lat, lon, maxID, assize, TR_data, cstr, mindist, testprint, slab)
    elif (length1 == length2 and mindist <= AA_data['dist'].max()) or (trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85):
        trimmed, maxID = trimByTrench(trimmed, oceanside, AA_data, lat, lon, maxID, assize, TR_data, cstr, mindist, testprint,slab)
    if len(trimmed)>0:
        trimmed = trimmed[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
        #loc_depth = trimmed['depth'].mean()
        elist, sdepth, ddepth, depthwritten = depthRange(loc_depth, sdr, ddr, seismo_thick, trimmed, slab, these_parameters, depthwritten)
        return trimmed, test, sdepth, ddepth, cstr, maxID, loc_depth, depthwritten, perpwritten
    else:
        test=False
        return trimmed, test, np.nan, np.nan, cstr, maxID, np.nan, True, True

def makealuAA(lon, lat, AA_data):

    distdata = AA_data.drop_duplicates(['sect'])
    if lon < distdata['avlon'].min():
        AA_data = AA_data[AA_data.sect == 'west']
    elif lon > distdata['avlon'].max():
        AA_data = AA_data[AA_data.sect == 'east']
    else:
        xdist = distdata['avlon'].values - lon
        ydist = distdata['avlat'].values - lat
        distdata['hdist'] = np.sqrt(xdist*xdist + ydist*ydist)
        mindist = distdata['hdist'].min()
        meddist = distdata['hdist'].median()
        minAA = distdata[distdata.hdist == mindist]
        medAA = distdata[distdata.hdist == meddist]
        minsect = minAA['sect'].values[0]
        medsect = medAA['sect'].values[0]
        if meddist > 2*mindist:
            AA_data = AA_data[AA_data.sect == minsect]
        else:
            mindata = AA_data[AA_data.sect == minsect]
            meddata = AA_data[AA_data.sect == medsect]
            maxoutmin = mindata['dist'].max()
            maxoutmed = meddata['dist'].max()
            minmin = max(mindata['dist'].min(), meddata['dist'].min())
            maxout = mindist/(mindist+meddist)*maxoutmed + meddist/(mindist+meddist)*maxoutmin
            if maxoutmed > maxoutmin:
                meddata = meddata[meddata.dist<maxout]
                add = meddata[meddata.dist > maxoutmin]
                mindata = pd.concat([mindata, add])
                mindata = mindata.reset_index(drop=True)
            else:
                mindata = mindata[mindata.dist<maxout]
                add = mindata[mindata.dist > maxoutmed]
                meddata = pd.concat([meddata, add])
                meddata = meddata.reset_index(drop=True)
            mindata = mindata[mindata.dist > minmin]
            meddata = meddata[meddata.dist > minmin]

            ind = meddata.dist.isin(mindata.dist)
            meddata = meddata[ind]
            mindata = meddata[ind]

            meddata2 = meddata[~ind]
            #if len(meddata2>0):
            #    print('MEDDATA2', meddata2)
            mindepths = mindata['depth'].values
            meddepths = meddata['depth'].values

            weighted_depths = mindist/(mindist+meddist)*meddepths + meddist/(mindist+meddist)*mindepths

            weightedAA = pd.DataFrame({'dist':mindata['dist'].values,'depth':weighted_depths})
            return weightedAA

    AA_data = AA_data[['dist', 'depth']]
    return AA_data

def getextraRF2(lon, lat, alen, blen, eventlist, elistRF, slab, cstr, mindist, trenchlon, trenchlat, AARF, maxID):

    if slab == 'alu' or slab == 'cam':
        elist = eventlist[eventlist.etype == 'RF']
        if len(elist)>0:
            #print 'elistRF2',lon,lat,elist
            elist = getEventsInCircle(lon, lat, 100, elist)
            if len(elist) > 0:
                #print 'elistRF1',lon,lat,elist
                elist = ellipseFilt(elist, lat, lon, 100, blen/3, cstr, mdist)
                if len(elist)>0:
                    mindist = elist['distance'].min()
                    elistRF = elist[elist.distance == mindist]
                    #print 'elistRF4',lon,lat,mindist,elistRF
                    if mindist > alen:
                        addunc = ((mindist-alen)/(100-alen)/2+1) * elistRF['unc'].values[0]
                        first = (mindist-alen)/(100-alen)+1
                        #print 'first',lon,lat,first,mindist,alen
                        #print 'second',lon,lat,elistRF['unc'].values[0]
                        elistRF['unc'] = addunc
                    #print 'elistRF5',lon,lat,elistRF
                    return elistRF
                else:
                    return elistRF
            else:
                return elistRF
        else:
            return elistRF
    else:
        return elistRF

def getextraRF(trimmed, slab, cstr, mindist, trenchlon, trenchlat, AARF, maxID, lon, lat):

    if slab == 'alx' or slab == 'cax':
        if (trenchlon>258.5 and trenchlon<262.41 and trenchlat>15.46 and trenchlat<16.91) or (trenchlon>212.86 and trenchlon<217.81 and trenchlat>58.026 and trenchlat<60.45):
            AARF['diffdist'] = np.abs(AARF['dist'].values - mindist)
            mindiff = AARF['diffdist'].min()
            if mindiff < 0.1:
                trimmed = trimmed[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src', 'distance']]
                thisAA = AARF[AARF.diffdist == mindiff]
                locAA = thisAA['depth'].values[0]
                trimmed.loc[len(trimmed)+1] = ([lat, lon, locAA, 10.0, str('AA'), (maxID+1), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                maxID += 1
                #print 'trimmed!',trimmed
            return trimmed, maxID
        else:
            return trimmed, maxID
    else:
        return trimmed, maxID


###############################################

### 42 ###

###############################################

## DEP.8.5.16
## Edited GLM 11.14.16 re-indented
## Edited GLM 11.15.16 - trimbytrench take mindist from gettrenchstrike
## Edited GLM 11.17.16 - return maxID
## Edited GLM 11.23.16 - remove TO where other etypes<300 km depth exist

# allFilters is what it sounds like - it applies all of the filters to the dataset, giving a final output dataset that we can use or adapt for making PDFs for a given node

def allFilters(eventlist, lat, lon, inside, slab1, strtmp, diptmp, seismo_thick, alen, blen, clen, mdist, sdr, ddr, mindip, maxID, out, AA_data, TR_data, slab, maxdist, testprint, extended, datainfo, nodeinfo, nID):

    uprad, dorad = 0, 0
    # Removing average active source gathered from other nodes
    eventlist = eventlist[eventlist.etype != 'AA'] # Separate out bathymetry data
    depthwritten, perpwritten, cutoffwritten = True, True, True
    these_parameters = []
    
    if slab == 'kur' or slab == 'jap' and len(eventlist[eventlist.etype == 'TO']) > 1:
        alen = alen*1.5
        
    if slab == 'hel' or slab == 'cas':
        TO_cutoff = 50
    elif slab == 'ryu':
        TO_cutoff = 200
        if lon > 136:
            TO_cutoff = 200
    elif slab == 'man':
        TO_cutoff = 150
    elif slab == 'sam' and lat > -15 and lat < -5:
        TO_cutoff = 100
    else:
        TO_cutoff = 300

    if slab == 'cam':
        distAA = AA_data[np.isnan(AA_data.avlon)]['dist'].max()
        depthAA = AA_data[np.isnan(AA_data.avlon)]['depth'].max()
    else:
        distAA = AA_data['dist'].max()
        depthAA = AA_data['depth'].max()
    AARF = pd.DataFrame() # GLM 05.04.17 remove eventually

    if slab == 'sum' and lat > 22:
        AA_data = AA_data[AA_data.dist < -2000]
    if slab == 'ryu' and lon > 136:
        AA_data = AA_data[AA_data.dist < 100]

    # Creating a list of events within a search radius
    elistPD = getEventsInCircle(lon, lat, clen, eventlist)
    #if testprint:
    #    idlist = list(elistPD['ID'].values)
    #    noelist = eventlist[~eventlist['ID'].isin(idlist)]
    #    addToDataInfo(noelist, nID, 'getEventsInCircle', datainfo,'df')

    if diptmp > 85 or slab == 'sam' or slab == 'him' or (slab == 'sol' and diptmp > 60):
        elistRF01 = getEventsInCircle(lon, lat, clen, eventlist[eventlist.etype == 'RF'])
    else:
        if slab == 'cam':
            elistRF01 = getEventsInCircle(lon, lat, clen*3, eventlist[eventlist.etype == 'RF'])
        elif slab == 'him':
            elistRF011 = getEventsInCircle(lon, lat, clen*3, eventlist[(eventlist.etype == 'RF')&(eventlist.src != 'schulte')])
            elistRF012 = getEventsInCircle(lon, lat, clen*2, eventlist[(eventlist.etype == 'RF')&(eventlist.src == 'schulte')])
            elistRF01 = pd.concat([elistRF011,elistRF012])
        else:
            elistRF01 = getEventsInCircle(lon, lat, clen*2, eventlist[eventlist.etype == 'RF'])

    elistCP01 = getEventsInCircle(lon, lat, blen, eventlist[eventlist.etype == 'CP'])

    #if testprint:
    #    idlistRF = list(elistRF01['ID'].values)
    #    noelist = eventlist[~((eventlist['ID'].isin(idlistRF))&(eventlist['ID'].isin(idlistRF)))]
    #    addToDataInfo(noelist, nID, 'getEventsInCircle (reciver functions)', datainfo,'df')

    balist = elistPD[elistPD.etype == 'BA']
    if (slab != 'alu' or (lon > 205 and lon < 215)) and slab != 'ryu':
        aslist = elistPD[elistPD.etype == 'AS']
    else:
        aslist = elistPD[elistPD.etype == 'blah']

    if testprint:
        f = open(nodeinfo, 'a')
        f.write('-%i- TO_cutoff, seismo_thick, distAA, depthAA, len(elistPD) len(balist),len(aslist),len(elistRF01),len(elistCP01) %i,%i,%0.2f,%0.2f,%i,%i,%i,%i,%i \n'%(nID, TO_cutoff, seismo_thick, distAA, depthAA, len(elistPD),len(balist),len(aslist),len(elistRF01),len(elistCP01)))
        f.close()

    # Get strike of, distance to, and angle between this node and nearest trench point
    if len(TR_data)>0:
        cstrB, cstrP, minang, mindist, tooFar, trenchlon, trenchlat, trenchdepth = getTrenchStrike(TR_data, lat, lon, maxdist, testprint)
        if len(AA_data)>0 and mindist <= distAA:
            AA_data['absdists'] = np.abs(AA_data['dist'].values)
            mindist1 = AA_data['absdists'].min()
            aamin = AA_data[AA_data.absdists == mindist1]
            aamindepth = aamin['depth'].values[0]
            AAdiff = aamindepth-trenchdepth
            AA_data['depth'] = AA_data['depth'].values - AAdiff

        if len(AA_data)>0 and slab == 'ryu' and trenchlat>37:
            AA_data = AA_data[AA_data.dist<10]
    
    else:
        cstrB, cstrP, minang, mindist, tooFar, trenchlon, trenchlat = strtmp, strtmp, 360, 400, False, 126.0, 1.0
        AA_data = AA_data[AA_data.dist < -2000]

    if slab != 'sam':
        if slab == 'cam' or slab == 'himz':
            elistRF0 = ellipseFilt(elistRF01, lat, lon, clen*3, blen/2, cstrB, mdist)
        else:
            elistRF0 = ellipseFilt(elistRF01, lat, lon, alen*2, blen/2, cstrB, mdist)
    else:
        elistRF0 = ellipseFilt(elistRF01, lat, lon, alen, blen/2, cstrB, mdist)

    elistCP0 = ellipseFilt(elistCP01, lat, lon, blen, blen/2, cstrB, mdist)

    if slab == 'him' and mindist > 175:# and len(elistRF0) < 100:
        test = False
        return elistPD[elistPD.etype == 'XX'], test, uprad, dorad, strtmp, diptmp, maxID, slab1

    if testprint:
        idlistRF = list(elistRF0['ID'].values)
        noelist = elistRF01[~((elistRF01['ID'].isin(idlistRF))&(elistRF01['ID'].isin(idlistRF)))]
        addToDataInfo(noelist, nID, 'getEventsInEllipse (reciver functions)', datainfo,'df')
        idlistCP = list(elistCP0['ID'].values)
        noelist = elistCP01[~((elistCP01['ID'].isin(idlistCP))&(elistCP01['ID'].isin(idlistCP)))]
        addToDataInfo(noelist, nID, 'getEventsInEllipse (control points)', datainfo,'df')
        
    if testprint:
        print ('lon,lat,slab,distAA,depthAA,mindist, len(elistPD)',lon,lat,slab,distAA,depthAA,mindist,len(elistPD))
    if slab == 'sam' and mindist > distAA:
        if testprint:
            shallowAA = elistPD[elistPD.depth<=depthAA]
            addToDataInfo(shallowAA, nID, 'shallower than average profile', datainfo,'df')
        elistPD = elistPD[elistPD.depth>depthAA]

    #if (slab == 'alu' or slab == 'sam') and len(balist) < 1 and not inside:
    if ((slab == 'alu' and len(balist) < 1 and not out and not inside) or (slab != 'alu' and not out)) and slab != 'hal' and slab != 'him' and slab != 'pam' and slab != 'hin':
        opposite = 50
        adjacent = 400
        dipcut = math.atan2(opposite, adjacent)
        depthcut = mindist*math.tan(dipcut)
        if testprint:
            shallowAA = elistPD[elistPD.depth <= depthcut]
            addToDataInfo(shallowAA, nID, 'dipdist y=mx+b', datainfo,'df')
            
        elistPD = elistPD[elistPD.depth > depthcut]
        
    if testprint:
        print(lat, lon, 'cstrB, minang, mindist, tooFar (getTrenchStrike)', cstrB, minang, mindist, tooFar, trenchlon, trenchlat, len(elistPD))
        f = open(nodeinfo, 'a')
        f.write('-%i- cstrB, minang, mindist, tooFar, trenchlon, trenchlat, len(elistPD)      %.2f, %.2f, %.2f, %s, %.2f, %.2f, %i \n'%(nID, cstrB, minang, mindist, tooFar, trenchlon, trenchlat, len(elistPD)))
        f.close()
        
    if inside or extended: # Get strike
        cstr = strtmp
        cdip = diptmp

    else:
        cstr = cstrB
        cdip = 0.0
        if slab == 'sam':
            slab1 = elistPD['depth'].mean()

    if math.isnan(cstr) or (slab == 'alu' and lat > 57 and lat < 60 and lon > 207 and lon <215):
        cstr = cstrB
        cdip = 0.0
    if cstr < 0:
        cstr += 360
    
    if len(balist)>0 and out:
        slab1 = balist['depth'].mean()
    
    elistnotTO = elistPD[(elistPD.etype != 'TO')]
    if len(elistnotTO) > 0 and slab1 <= TO_cutoff:
        if testprint:
            shallowAA = elistPD[(elistPD.etype == 'TO')]
            addToDataInfo(shallowAA, nID, 'shallow and other data available besides TO', datainfo,'df')
            
        elistPD = elistPD[(elistPD.etype != 'TO')]

    if testprint:
        print ('111',lon,lat,elistPD,mindist,distAA,trenchlon,trenchlat,slab1,elistPD['distance'].values)
        f = open(nodeinfo, 'a')
        f.write('-%i- len(elistPD),len(elistnotTO),slab1,inside,extended,out,cstr,cdip %i,%i,%.2f,%s,%s,%s,%.2f,%.2f \n'%(nID,len(elistPD),len(elistnotTO),slab1,inside,extended,out,cstr,cdip))
        f.close()
        
    if slab == 'cam':
        #if trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85 and mindist>distAA:
        #    AA_data = AA_data[np.isfinite(AA_data.avlon)]
        #else:
        #    AA_data = AA_data[np.isnan(AA_data.avlon)]
        AA_data = AA_data[np.isnan(AA_data.avlon)]

    if (len(elistPD) < 2 and mindist <= distAA and not out and len(aslist)<1) or (len(elistPD)<2 and (trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85)):
        trimmed, test, sdepth, ddepth, cstr, maxID, loc_depth, depthwritten, perpwritten = noDataNeedAA(elistPD, cstr, minang, AA_data, lat, lon, maxID, TR_data, mindist, testprint, sdr, ddr, seismo_thick, slab, these_parameters, depthwritten, perpwritten, trenchlon, trenchlat, AARF, slab1)
        
        if testprint:
            idlist = list(trimmed['ID'].values)
            noelist = elistPD[~elistPD['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'after first nodataneedAA', datainfo,'df')
            f = open(nodeinfo, 'a')
            f.write('-%i- exited first nodataneedAA, len(trimmed), %i \n'%(nID,len(trimmed)))
            f.close()
            
        return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth

    loc_depth, elist = findLocDep(slab1, tooFar, elistPD, seismo_thick, testprint, balist, out, slab, lon, lat)
    
    if len(elist)<2 and len(elistCP0)>0:
        elist = pd.concat([elist,elistCP0])
        if len(elist)<2:
            elist = pd.concat([elistCP0,elistCP0])

    if testprint:
        idlist = list(elistPD['ID'].values)
        noelist = elist[~elist['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'findLocDep', datainfo,'df')
        f = open(nodeinfo, 'a')
        f.write('-%i- after findlocdep, loc_dep, len(elist), %.2f, %i \n'%(nID,loc_depth,len(elist)))
        f.close()

    if (len(elist) < 2 and mindist <= distAA and not out and len(aslist)<1) or (len(elist)<2 and (trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85)):
        trimmed, test, sdepth, ddepth, cstr, maxID, loc_depth, depthwritten, perpwritten = noDataNeedAA(elist, cstr, minang, AA_data, lat, lon, maxID, TR_data, mindist, testprint, sdr, ddr, seismo_thick, slab, these_parameters, depthwritten, perpwritten, trenchlon, trenchlat, AARF, loc_depth)
        
        if testprint:
            idlist = list(trimmed['ID'].values)
            noelist = elist[~elist['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'after second nodataneedAA', datainfo,'df')
            f = open(nodeinfo, 'a')
            f.write('-%i- exited second nodataneedAA, len(trimmed), %i \n'%(nID,len(trimmed)))
            f.close()
          
        return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth

    elistBA1 = elist[elist.etype == 'BA'] # Separate out bathymetry data
    elistAS1 = elist[elist.etype == 'AS'] # Separate out active source data
    elistRF1 = elist[elist.etype == 'RF']
    elistCP1 = elist[elist.etype == 'CP']
    elistBA = ellipseFilt(elistBA1, lat, lon, alen, blen, cstr, mdist)  # Filter by ellipse
    elistAS2 = ellipseFilt(elistAS1, lat, lon, alen, blen, cstr, mdist)  # Filter by ellipse
    if len(elistAS2) > 1:
        elistAS=elistAS2[elistAS2.distance == elistAS2['distance'].min()]
    else:
        elistAS = elistAS2.copy()
    elistRF = ellipseFilt(elistRF1, lat, lon, alen, blen, cstr, mdist)  # Filter by ellipse
    elistCP = ellipseFilt(elistCP1, lat, lon, alen, blen, cstr, mdist)
    if testprint:
        idlist = list(elistAS2['ID'].values)
        noelist = elistAS1[~elistAS1['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'AS ellipse filt', datainfo,'df')
        
        idlist = list(elistAS['ID'].values)
        noelist = elistAS2[~elistAS2['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'only take closest AS', datainfo,'df')
        
        idlist = list(elistBA['ID'].values)
        noelist = elistBA1[~elistBA1['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'BA ellipse filt', datainfo,'df')
        
        idlist = list(elistRF['ID'].values)
        noelist = elistRF1[~elistRF1['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'RF ellipse filt', datainfo,'df')
        
        idlist = list(elistCP['ID'].values)
        noelist = elistCP1[~elistCP1['ID'].isin(idlist)]
        addToDataInfo(noelist, nID, 'CP ellipse filt', datainfo,'df')

        f = open(nodeinfo, 'a')
        f.write('-%i- filtering special data, len(elistAS2), len(elistAS), len(elistBA1), len(elistBA), len(elistRF1), len(elistRF), len(elistCP) %i %i %i %i %i %i %i \n'%(nID,len(elistAS2), len(elistAS), len(elistBA1), len(elistBA), len(elistRF1), len(elistRF), len(elistCP)))
        f.close()
            
    if len(elist)>1 and (cdip > mindip or (len(elistBA)<1 and extended)):
        trimmed1, uprad, dorad, cutoffwritten = dualdepthperp(loc_depth, sdr, ddr, seismo_thick, elist, slab, cstr, lon, lat, cdip, alen, blen, these_parameters, cutoffwritten)
        
        if testprint:
            idlist = list(trimmed1['ID'].values)
            noelist = elist[~elist['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'OG dual depth perp', datainfo,'df')
        
        to_trimmed1 = trimmed1[trimmed1.etype == 'TO']
        if len(to_trimmed1)>0:
            if slab == 'sum' or slab == 'manz':
                to_trimmed, tosdepth, toddepth, ctwrt = dualdepthperp(loc_depth, sdr, ddr, seismo_thick, to_trimmed1, slab, cstr, lon, lat, cdip, alen, blen, these_parameters, cutoffwritten)
            elif slab == 'sam' and lat > -15 and lat < -11:
                to_trimmed, tosdepth, toddepth, ctwrt = dualdepthperp(loc_depth, sdr, ddr, seismo_thick, to_trimmed1, slab, cstr, lon, lat, cdip, 100, blen, these_parameters, cutoffwritten)
            else:
                to_trimmed, tosdepth, toddepth, ctwrt = dualdepthperp(loc_depth, sdr, ddr, seismo_thick, to_trimmed1, slab, cstr, lon, lat, cdip, blen, blen, these_parameters, cutoffwritten)
            trimmed1 = pd.concat([trimmed1[trimmed1.etype != 'TO'], to_trimmed])
        else:
            to_trimmed = to_trimmed1.copy()
        elistRF = trimmed1[trimmed1.etype == 'RF']
        elistAS = trimmed1[trimmed1.etype == 'AS']
        elistBA = trimmed1[trimmed1.etype == 'BA']
        elistCP = trimmed1[trimmed1.etype == 'CP']

        if testprint:
            idlist = list(to_trimmed['ID'].values)
            noelist = to_trimmed1[~to_trimmed1['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'reduced TO proj dual depth perp', datainfo,'df')

        if len(elistAS)<1 and len(elistRF)<1 and len(elistCP)<1:
            trimmed = trimmed1.copy()
        else:
            if len(elistAS) > 0:
                elistAS=elistAS[elistAS.distance == elistAS['distance'].min()]
                trimmed = trimmed1[trimmed1.etype != 'AS']
                trimmed = pd.concat([trimmed, elistAS])
            else:
                trimmed = trimmed1.copy()

            if len(elistRF) > 0:
                elistRF=elistRF[elistRF.distance == elistRF['distance'].min()]
                trimmed = trimmed[trimmed.etype != 'RF']
                trimmed = pd.concat([trimmed, elistRF])
                
            if len(elistCP) > 0:
                elistCP=elistCP[elistCP.distance == elistCP['distance'].min()]
                trimmed = trimmed[trimmed.etype != 'CP']
                trimmed = pd.concat([trimmed, elistCP])

        if testprint:
            idlist = list(trimmed['ID'].values)
            noelist = trimmed1[~trimmed1['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'only take nearest RF and AS and CP', datainfo,'df')
            
            f = open(nodeinfo, 'a')
            f.write('-%i- after dualdepthperp loc_depth, sdr, ddr, len(elist), len(trimmed1), len(to_trimmed1), len(to_trimmed), len(elistBA), len(elistRF), len(elistAS), len(elistCP), len(trimmed), sdepth, ddepth, alen, blen, clen, cstr, cdip %.2f, %.2f, %.2f, %i, %i, %i, %i, %i, %i, %i, %i, %i, %.2f, %.2f, %i, %i, %.2f, %.2f, %.2f  \n'%(nID,loc_depth, sdr, ddr, len(elist), len(trimmed1), len(to_trimmed1), len(to_trimmed), len(elistBA), len(elistRF), len(elistAS), len(elistCP), len(trimmed), sdepth, ddepth, alen, blen, clen, cstr, cdip))
            f.close()
            
    elif len(elist) > 1:
        elist = elist[elist.etype != 'AS'] # Make new elist without active source
        elist = elist[elist.etype != 'BA'] # Make new elist without bathymetry
        elist = elist[elist.etype != 'RF']
        elist = elist[elist.etype != 'CP']
        elist2, sdepth, ddepth, depthwritten = depthRange(loc_depth, sdr, ddr, seismo_thick, elist, slab, these_parameters, depthwritten)
        
        uprad = loc_depth-sdepth
        dorad = ddepth-loc_depth
        #print 'testing depthRange lon,lat,sdepth,ddepth,loc_depth,alen,blen,sdr,ddr',lon,lat,sdepth,ddepth,loc_depth,alen,blen,sdr,ddr
        if testprint:
            idlist = list(elist2['ID'].values)
            noelist = elist[~elist['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'depthrange', datainfo,'df')
            
        elist3 = ellipseFilt(elist2, lat, lon, alen, blen, cstr, mdist)
        
        if testprint:
            idlist = list(elist3['ID'].values)
            noelist = elist2[~elist2['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'ellipsefilt', datainfo,'df')
            
        frames = [elist3, elistAS, elistBA, elistRF, elistCP]
        trimmed = pd.concat(frames)  # Add back in active source and bathymetry data
        
        if testprint:
            f = open(nodeinfo, 'a')
            f.write('-%i- after deprange and ellipsefilt loc_depth, sdr, ddr, len(elist), len(elist2), len(elist3), len(elistAS), len(elistBA), len(elistRF), len(elistCP), sdepth, ddepth, alen, blen, clen, cstr, cdip %.2f, %.2f, %.2f, %i, %i, %i, %i, %i, %i, %i, %.2f, %.2f, %i, %i, %.2f, %.2f, %.2f  \n'%(nID,loc_depth, sdr, ddr, len(elist), len(elist2), len(elist3), len(elistAS), len(elistBA), len(elistRF), len(elistCP),sdepth, ddepth, alen, blen, clen, cstr, cdip))
            f.close()
        
    elif len(elistBA)>0 or len(elistAS)>0 or len(elistRF)>0 or len(elistCP)>0:
        #print 'only RF',lon,lat
        trimmed = pd.concat([elistAS, elistBA, elistRF, elistCP])
        sdepth, ddepth = loc_depth-sdr, loc_depth+sdr
        uprad = sdr
        dorad = sdr
    else:
        sdepth, ddepth = loc_depth-sdr, loc_depth+sdr
        uprad = sdr
        dorad = sdr
        if (mindist <= distAA and not out and len(aslist)<1) or ((trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85)):
            trimmed, test, sdepth, ddepth, cstr, maxID, loc_depth, depthwritten, perpwritten = noDataNeedAA(elist, cstr, minang, AA_data, lat, lon, maxID, TR_data, mindist, testprint, sdr, ddr, seismo_thick, slab, these_parameters, depthwritten, perpwritten, trenchlon, trenchlat, AARF, loc_depth)
            
            if testprint:
                idlist = list(trimmed['ID'].values)
                noelist = elist[~elist['ID'].isin(idlist)]
                addToDataInfo(noelist, nID, 'after third nodataneedAA', datainfo,'df')
                f = open(nodeinfo, 'a')
                f.write('-%i- exited third nodataneedAA, len(trimmed), %i \n'%(nID,len(trimmed)))
                f.close()
            
            return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth

        else:  # Skip nodes with no data
            test = False
            return elist, test, uprad, dorad, cstr, cdip, maxID, loc_depth

    elistRF0, elistRF = removematches(elistRF0,elistRF)
    if len(elistRF0)>0:
        elistRF0['unc'] = elistRF0['unc'].values*(2*((clen*3-elistRF0['distance'].values)/(clen*3)))
        elistRF0.loc[elistRF0.distance < alen, 'unc'] = 10.0
        trimmed = pd.concat([trimmed,elistRF0])

    elistCP0, elistCP = removematches(elistCP0,elistCP)
    if len(elistCP0)>0:
        #elistCP0['unc'] = elistCP0['unc'].values*(2*((clen*3-elistCP0['distance'].values)/(clen*3)))
        elistCP0.loc[elistCP0.distance < alen, 'unc'] = 10.0
        trimmed = pd.concat([trimmed,elistCP0])
        if slab == 'puy' and len(trimmed[trimmed.etype == 'CP'])>1:
            cptrimmed = trimmed[trimmed.etype == 'CP']
            ottrimmed = trimmed[trimmed.etype != 'CP']
            cptrimmed=cptrimmed[cptrimmed.distance == cptrimmed['distance'].min()]
            trimmed = pd.concat([ottrimmed,cptrimmed])

    if len(trimmed)<1 and (len(elistRF)>0 or len(elistCP)>0):
        #print 'had to add again,lon,lat,cstr,cdip,mindist,len(elistRF),elistRF,trimmed',lon,lat,cstr,cdip,mindist,len(elistRF),elistRF,trimmed
        trimmed = pd.concat([trimmed, elistRF, elistCP])

    if len(trimmed) < 2 and mindist > distAA and len(elistRF)<1 and len(elistCP)<1 and ((slab != 'ryu' and slab != 'hel') or len(trimmed[trimmed.etype == 'TO'])<1):  # Skip nodes with no data
        test = False
        return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth
    elif (len(trimmed) < 2 and mindist <= distAA and not out and len(aslist)<1) or (len(trimmed)<2 and (trenchlon>258.7 and trenchlon<260.7 and trenchlat>16.0 and trenchlat<16.85) and len(elistRF)<1):
        trimmed2, test, sdepth, ddepth, cstr, maxID, loc_depth, depthwritten, perpwritten = noDataNeedAA(trimmed, cstr, minang, AA_data, lat, lon, maxID, TR_data, mindist, testprint, sdr, ddr, seismo_thick, slab, these_parameters, depthwritten, perpwritten, trenchlon, trenchlat, AARF, loc_depth)
        
        if testprint:
            idlist = list(trimmed2['ID'].values)
            noelist = trimmed[~trimmed['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'after fourth nodataneedAA', datainfo,'df')
            f = open(nodeinfo, 'a')
            f.write('-%i- exited fourth nodataneedAA, len(trimmed), %i \n'%(nID,len(trimmed2)))
            f.close()
                
        return trimmed2, test, uprad, dorad, cstr, cdip, maxID, loc_depth

    if mindist <= distAA and not out and len(aslist)<1: # GLM 11.21.16
        if testprint:
            ogtrimmed = trimmed.copy()
            maxIDbefore = maxID
        length1 = len(trimmed)
        trimmed, maxID = getextraRF(trimmed, slab, cstr, mindist, trenchlon, trenchlat, AARF, maxID, lon, lat)
        length2 = len(trimmed)
        if (slab == 'alu' or slab == 'him') and length1 == length2:
            trimmed, maxID = trimByTrench_alu(trimmed, out, AA_data, lat, lon, maxID, elistAS.size, TR_data, cstr, mindist, testprint, slab)
        elif length1 == length2:
            trimmed, maxID = trimByTrench(trimmed, out, AA_data, lat, lon, maxID, elistAS.size, TR_data, cstr, mindist, testprint, slab)
            
        if testprint:
            idlist = list(trimmed['ID'].values)
            noelist = ogtrimmed[~ogtrimmed['ID'].isin(idlist)]
            addToDataInfo(noelist, nID, 'after trimbytrench', datainfo,'df')
            f = open(nodeinfo, 'a')
            f.write('-%i- after trimbytrench, %i, %i, %i \n'%(nID,len(trimmed),maxIDbefore,maxID))
            f.close()
            
    else:
        if out:
        
            if testprint:
                noelist = trimmed[(trimmed.etype !='BA') & (trimmed.etype !='AS')]
                addToDataInfo(noelist, nID, 'removed all non BA/AS because outboard', datainfo,'df')
                
            trimmed = trimmed[(trimmed.etype =='BA') | (trimmed.etype =='AS')]
            
        else:
        
            if testprint:
                noelist = trimmed[trimmed.etype =='BA']
                addToDataInfo(noelist, nID, 'removed all BA because inboard and far from trench', datainfo,'df')
                
            trimmed = trimmed[trimmed.etype != 'BA']

    if slab == 'alu' and len(balist) < 1 and not inside and slab != 'hal' and slab != 'him' and slab != 'pam' and slab != 'hin':
        opposite = 50
        adjacent = 400
        dipcut = math.atan2(opposite, adjacent)
        depthcut = mindist*math.tan(dipcut)
        
        if testprint:
            shallowAA = trimmed[trimmed.depth <= depthcut]
            addToDataInfo(shallowAA, nID, 'dipdist y=mx+b round2', datainfo,'df')
            
        trimmed = trimmed[trimmed.depth > depthcut]

    if len(trimmed) < 1:
        test = False
        return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth

    if testprint:
        print('lon,lat,sdepth,ddepth,cstr,maxID,loc_depth,trimmed', lon, lat, sdepth, ddepth, cstr, maxID, loc_depth, trimmed)
        
    if len(trimmed[trimmed.etype != 'CP']) > 10 and slab != 'puy' and slab != 'him':
        trimmed = trimmed[trimmed.etype != 'CP']

    test = True

    if slab == 'him' and len(trimmed[trimmed.etype == 'CP'])>0:
        trimmed = trimmed[trimmed.etype != 'AA']

    if slab == 'kur' and len(trimmed[trimmed.etype == 'TO']) > 0:
        if (len(trimmed[trimmed.etype == 'EQ']) > 0 or len(trimmed[trimmed.etype == 'ER']) > 0):
            trimmed = trimmed[trimmed.etype != 'TO']
    return trimmed, test, uprad, dorad, cstr, cdip, maxID, loc_depth

###############################################

### 43 ###

###############################################

## Written MF 8.2.16 
## DEP.8.5.16 edited
## GLM 11.17.16 edited - use dataframe column names instead of indices

# avStrDipRak finds the average strike dip and rake of the shallow planes of the CMT solutions included for a given node

def avStrDipRak(trimmed):

    EQframe = trimmed[trimmed.etype == 'EQ']

    if len(EQframe > 0) and EQframe[pd.notnull(EQframe['S1'])].size > 0:
        EQ_with_cmt = EQframe[pd.notnull(EQframe['S1'])]
        EQ_with_cmt_array = np.array(EQ_with_cmt)
        S1df = EQ_with_cmt['S1'].values
        S2df = EQ_with_cmt['S2'].values
        R1df = EQ_with_cmt['R1'].values
        R2df = EQ_with_cmt['R2'].values
        D1df = EQ_with_cmt['D1'].values
        D2df = EQ_with_cmt['D2'].values
        
    
        #get strike, dip, and rake from the shallow plane for each event
            #and average these values for each node
            #note: the column orders WILL change if new data are imported

        dip = np.ones(len(EQ_with_cmt_array))*-9999
        strike = np.ones(len(EQ_with_cmt_array))*-9999
        rake = np.ones(len(EQ_with_cmt_array))*-9999
        
        for i in range(len(EQ_with_cmt_array)):
            strike1 = S1df[i]
            strike2 = S2df[i]
            rake1 = R1df[i]
            rake2 = R2df[i]
            dip1 = D1df[i]
            dip2 = D2df[i]
            if dip1 >= dip2:
                dip[i] = dip2
                strike[i] = strike2
                rake[i] = rake2
            elif dip1 < dip2:
                dip[i] = dip1
                strike[i] = strike1
                rake[i] = rake1

        dip = dip[dip>-999]
        strike = strike[strike>-999]
        rake = rake[rake>-999]
        
        avg_dip = np.mean(dip)
        avg_strike = np.mean(strike)
        avg_rake = np.mean(rake)

    else:
        avg_dip = np.nan
        avg_strike = np.nan
        avg_rake = np.nan

    return avg_strike, avg_dip, avg_rake

###############################################

### 43.5 ###

###############################################

def makeMultiDF(multi_peak2, multipeaks, lon, lat, nID):

    mdepths = multi_peak2['depths'].values
    mnodes = np.ones(len(mdepths)) * nID
    mlats = np.ones(len(mdepths)) * lat
    mlons = np.ones(len(mdepths)) * lon
    multiDF = pd.DataFrame({'lon':mlons,'lat':mlats,'depth':mdepths,'nID':mnodes})
    if len(multipeaks)>0:
        frames = [multipeaks, multiDF] # GLM 11.23.16
        multipeaks = pd.concat(frames)
        multipeaks = multipeaks.reset_index(drop=True)
        multipeaks = multipeaks[['lon', 'lat', 'depth', 'nID']]
    else:
        multipeaks = multipeaks.append(multiDF)
    return multipeaks

def makeMultiDFP(multi_peak2, multipeaks, lon, lat, nID, strike, dip, loc_depth):

    depphi = 90-abs(dip-90)
    mdepths = multi_peak2['depths'].values
    
    mdepthsOUT = mdepths[mdepths>0]
    mdepthsIN = mdepths[mdepths<=0]
    
    azOUT = az_perp(strike)
    azIN = az_other_perp(strike)
    
    perpdistOUT = mdepthsOUT*math.sin(math.radians(depphi))
    perpdistIN = mdepthsIN*math.sin(math.radians(depphi))
    
    lonsout, latsout = np.zeros(len(perpdistOUT)), np.zeros(len(perpdistOUT))
    lonsin, latsin = np.zeros(len(perpdistIN)), np.zeros(len(perpdistIN))
    
    for i in range(len(perpdistOUT)):
        if abs(perpdistOUT[i]) > 0.001:
            lonsout[i], latsout[i] = heading(lon, lat, abs(perpdistOUT[i]), azOUT)
        else:
            lonsout[i], latsout[i] = lon, lat
    for i in range(len(perpdistIN)):
        if abs(perpdistIN[i]) > 0.001:
            lonsin[i], latsin[i] = heading(lon, lat, abs(perpdistIN[i]), azIN)
        else:
            lonsin[i], latsin[i] = lon, lat

    perpdepthsOUT = mdepthsOUT*math.cos(math.radians(depphi))
    perpdepthsIN = mdepthsIN*math.cos(math.radians(depphi))
    
    mlons = np.concatenate((lonsout, lonsin))
    mlats = np.concatenate((latsout, latsin))
    perpdepths = np.concatenate((perpdepthsOUT, perpdepthsIN))
    mdepths = perpdepths+loc_depth
    
    mnodes = np.ones(len(mdepths)) * nID

    multiDF = pd.DataFrame({'lon':mlons,'lat':mlats,'depth':mdepths,'nID':mnodes})

    if len(multipeaks)>0:
        frames = [multipeaks, multiDF] # GLM 11.23.16
        multipeaks = pd.concat(frames)
        multipeaks = multipeaks.reset_index(drop=True)
        multipeaks = multipeaks[['lon', 'lat', 'depth', 'nID']]
    else:
        multipeaks = multipeaks.append(multiDF)
    return multipeaks

###############################################

### 43.75 ###

###############################################

def getLocalMax(multi_peak):

    sums = multi_peak['Summed_Values'].values
    depths = multi_peak['depths'].values
    n = len(sums)
    peaks = []
    
    if len(multi_peak)<2:
        return multi_peak
    for i in range(n):
        depthC = depths[i]
        sumC = sums[i]
        # if on the first point, must treat differently
        if i == 0:
            depthD = depths[i+1]
            sumD = sums[i+1]
            # if there are no other points around it, append peak
            if abs(depthC-depthD) > 1:
                peaks.append(depthC)
                continue
            # or, if this depth has a higher probability than the next, this must be the local peak
            elif sumC > sumD:
                peaks.append(depthC)
                continue
            # else, go to next point, additional calculations not necessary
            else:
                continue
        # if on the last point, must treat differently
        elif i == n-1:
            depthS = depths[i-1]
            sumS = sums[i-1]
            # if there are no other points around it, append peak
            if abs(depthC-depthS) > 1:
                peaks.append(depthC)
                continue
            # or, if this depth has a higher probability than the last, this must be the local peak
            elif sumC > sumS:
                peaks.append(depthC)
                continue
            # else, go to next point, additional calculations not necessary
            else:
                continue
        else:
            depthD = depths[i+1]
            depthS = depths[i-1]
            sumD = sums[i+1]
            sumS = sums[i-1]

        # if there are other points around this one but this probability is greater than both of them, append peak
        if abs(depthC-depthS) < 2 and abs(depthC-depthD) < 2:
            if sumC > sumS and sumC >sumD:
                peaks.append(depthC)
            else:
                continue

        # if there are only nearby points that are deeper, but this probability is higher, append peak
        elif abs(depthC-depthS) > 2 and abs(depthC-depthD) < 2:
            if sumC > sumD:
                peaks.append(depthC)
            else:
                continue
        # if there are only nearby points that are deeper, but this probability is higher, append peak
        if abs(depthC-depthS) < 2 and abs(depthC-depthD) > 2:
            if sumC > sumS:
                peaks.append(depthC)
            else:
                continue

        # if there are no other nearby points, this must be the local peak
        elif abs(depthC-depthS) > 2 and abs(depthC-depthD) > 2:
            peaks.append(depthC)
        
    peaks = np.array(peaks)
    multi_peak2 = multi_peak.loc[multi_peak['depths'].isin(peaks)]

    return multi_peak2

###############################################

### 44 ###

###############################################

## DEP.8.8.16
## Edited GLM 11.14.16 re-indented

# fullPDFcalc takes in a dataset of lat, lon, dep, and unc, calculates a summed pdf of the data, and then determines the peak depth and standard deviation of the dataset.

def fullPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, whichpdf, slab, cstr, cdip):

    multipeaks = pd.DataFrame()
    
    elistASe = trimmed[trimmed.etype == 'AS' ]
    elistBAe = trimmed[trimmed.etype == 'BA' ]
    elistAAe = trimmed[trimmed.etype == 'AA' ]
    elistRFe = trimmed[trimmed.etype == 'RF' ]
    elistCPe = trimmed[trimmed.etype == 'CP' ]
    elistTOe = trimmed[trimmed.etype == 'TO' ]
    
    if len(elistAAe)>0 and len(trimmed) <4:
        if abs(elistAAe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from AA',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]

    if len(elistASe)>0 and len(trimmed) <5:
        if abs(elistASe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from AS',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]

    if len(elistBAe)>0 and len(trimmed) <5:
        if abs(elistBAe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from BA',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]
    
    nantest = trimmed['depth'].values
    nantest = nantest[np.isnan(nantest)]
    
    if len(nantest) > 0 or np.isnan(sdepth) or np.isnan(ddepth) or np.isnan(loc_depth):
        print ('NAN problem?? lon,lat,nID,sdepth,ddepth,loc_depth,trimmed',lon,lat,nID,sdepth,ddepth,loc_depth,trimmed)
        peak_depth = np.nan
        stdv = np.nan
        test = False
        n = 0
        return peak_depth, stdv, test, n, multipeaks, stdv
    
    multi = False
    n = 0
    
    if len(trimmed)>1:

        # Distinguishing between different data types
        ASframe = trimmed[trimmed.etype == 'AS']
        AAframe = trimmed[trimmed.etype == 'AA']
        EQframe = trimmed[trimmed.etype == 'EQ']
        BAframe = trimmed[trimmed.etype == 'BA']
        ERframe = trimmed[trimmed.etype == 'ER']
        TOframe = trimmed[trimmed.etype == 'TO']
        RFframe = trimmed[trimmed.etype == 'RF']
        CPframe = trimmed[trimmed.etype == 'CP']

        # Adding present event types to list of event types
        #and calculate average rake, strike, and dip for output file if CMT info available
        etypes = []
        AA = False
        AS = False
        BA = False
        RF = False
        TO = False
        ER = False
        EQ = False
        CP = False
        if len(ASframe) > 0:
            etypes.append('AS')
            AS = True
        if len(AAframe) > 0:
            etypes.append('AA')
            AA = True
        if len(EQframe) > 0 or len(ERframe)>0:
            etypes.append('EQ')
            if len(EQframe) > 0:
                EQ = True
            if len(ERframe) > 0:
                ER = True
        if len(BAframe) > 0:
            etypes.append('BA')
            BA = True
        #if len(ERframe > 0):
        #    etypes.append('ER')
        #    ER = True
        if len(TOframe) > 0:
            etypes.append('TO')
            TO = True
        if len(RFframe) > 0:
            etypes.append('RF')
            RF = True
        if len(CPframe) > 0:
            etypes.append('CP')
            CP = True

        # Make PDF
        #changed the values from 15 to 50, tbd whether this is a good idea or not!
        if ddepth > 1000:
            ddepth = np.max(trimmed['depth'].values)+10
        dep_range = np.arange(sdepth - 50, ddepth + 50, 1)
        PDF = makePDF4(trimmed, dep_range, etypes, testprint, 'depth')
        PDF_df1 = DataFrame(dep_range, columns=['depths'])
        PDF_df1['Summed_Values'] = PDF
        
        # Eliminates values less than 0.001 and finds min, max, and peak depth in PDF
        
        if len(PDF_df1) > 0:
            PDF_df = PDF_df1.loc[PDF_df1.Summed_Values >= 0.001]
            if len(PDF_df) < 1:
                PDF_df = PDF_df1.loc[PDF_df1.Summed_Values >= 0.0001]
                if len(PDF_df) < 1:
                    print ('noPDF? lon,lat,nID,sdepth,ddepth,loc_depth,trimmed',lon,lat,nID,sdepth,ddepth,loc_depth,trimmed)
                    peak_depth = np.nan
                    stdv = np.nan
                    test = False
                    n = 0
                    return peak_depth, stdv, test, n, multipeaks, stdv
        
        else:
            #print 'noPDF? lon,lat,nID,sdepth,ddepth,loc_depth,trimmed',lon,lat,nID,sdepth,ddepth,loc_depth,trimmed
            peak_depth = np.nan
            stdv = np.nan
            test = False
            n = 0
            return peak_depth, stdv, test, n, multipeaks, stdv
        
        if AA or AS or BA or RF or TO or CP or (ER and EQ and slab != 'kur'):
        #if AA or AS or BA or RF or TO or ER:
            peak = PDF_df['Summed_Values'].max()
            peakbuffer = 0.1*peak
            depthbuffer = 10
            d_min = PDF_df['depths'].min()
            d_max = PDF_df['depths'].max()

            # Finding the depth associated with the peak PDF value
            peak_df = PDF_df[PDF_df.Summed_Values == peak]
            peak_depth = peak_df['depths'].values[0]
            meandepth = False

        else:
            peak_depth = PDF_df['depths'].mean()
            peak_df = PDF_df[PDF_df.depths == peak_depth]
            peak = 1
            peakbuffer = 0.01
            meandepth = True

        # GLM 11.22.16 - adding bimodal distribution condition
        PDF_df['buffer'] = PDF_df['Summed_Values'].values + peakbuffer
        multi_peak = PDF_df[PDF_df.buffer >= peak]
        #multi_peak = PDF_df[(PDF_df.buffer >= peak) & ((PDF_df.depths < peakmin) | (PDF_df.depths > peakmax))] # GLM 11.25.16
        multi_peak2 = getLocalMax(multi_peak)
        
        if len(multi_peak2)>1 and not meandepth:
            multipeaks = makeMultiDF(multi_peak2, multipeaks, lon, lat, nID)
            multi = True
            test = True
            n = len(multi_peak2)
        else:
            try:
                peak_depth = peak_depth # gets the value out of the array
                test = True
                multi = False
                n = 1
            except:
                #print 'multidepth PDF Exception: lon,lat,nID: ',lon,lat,nID
                test = False
                stdv = np.nan # GLM 11.14.16 investigate this exception if missing PDFs
                peak_depth = np.nan
                return peak_depth, stdv, test, 0, multipeaks, stdv

        # Finding standard deviation of PDF
        thissum = 0
        for d in PDF_df['depths'].values:
            residual = peak_depth - d
            thissum += residual * residual
        stdv = math.sqrt(1.0/len(PDF_df)*thissum)
        
        minonperp = PDF_df['depths'].min()
        centsurf = abs(peak_depth-minonperp)
        
        # For testing PDFs of specific points - change lat-lon ranges to use
        if testprint:
            
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(121)
            
            thispoint = ax1.plot([0], [0], 'ro', label='Node Location')
            trimmed['lonplot'] = trimmed['lon'].values-lon
            trimmed['latplot'] = trimmed['lat'].values-lat
            
            if len(BAframe)>0:
                BA2 = trimmed[trimmed.etype == 'BA']
                bap = ax1.plot(BA2['lonplot'].values, BA2['latplot'].values, 'r.', label='BA')
            if len(EQframe)>0:
                EQ2 = trimmed[trimmed.etype == 'EQ']
                eqp = ax1.plot(EQ2['lonplot'].values, EQ2['latplot'].values, 'c.', label='EQ')
            if len(ERframe)>0:
                ER2 = trimmed[trimmed.etype == 'ER']
                erp = ax1.plot(ER2['lonplot'].values, ER2['latplot'].values, 'y.', label='ER')
            if len(AAframe)>0:
                AA2 = trimmed[trimmed.etype == 'AA']
                aap = ax1.plot(AA2['lonplot'].values, AA2['latplot'].values, 'k.', label='AA')
            if len(ASframe)>0:
                AS2 = trimmed[trimmed.etype == 'AS']
                asp = ax1.plot(AS2['lonplot'].values, AS2['latplot'].values, 'm.', label='AS')
            if len(TOframe)>0:
                TO2 = trimmed[trimmed.etype == 'TO']
                top = ax1.plot(TO2['lonplot'].values, TO2['latplot'].values, 'g.', label='TO')
            if len(RFframe)>0:
                RF2 = trimmed[trimmed.etype == 'RF']
                rfp = ax1.plot(RF2['lonplot'].values, RF2['latplot'].values, 'b.', label='RF')
            if len(CPframe)>0:
                CP2 = trimmed[trimmed.etype == 'CP']
                CPp = ax1.plot(CP2['lonplot'].values, CP2['latplot'].values, color='orange',marker='.', label='CP')
        
            ax1.set_xlabel('Longitude Difference From Node Coordinate')
            ax1.set_ylabel('Latitude Difference From Node Coordinate')
            ax1.axis('equal')
            plt.grid()
            title = 'Lat: %.2f, Lon: %.2f, Strike: %.2f, Dip: %.2f, Origin Depth: %.2f' % (lat, lon, cstr, cdip, loc_depth)
            ax1.set_title(title)
            lontit = lon*100
            lattit = lat*100
            ax1.legend(loc='best')

            a2 = (lat-trimmed['lat'])*(lat-trimmed['lat'])
            b2 = (lon-trimmed['lon'])*(lon-trimmed['lon'])
            c = np.sqrt(a2+b2)/2
            
            ax2 = fig.add_subplot(122)
            if len(BAframe)>0:
                BAa2 = (lat-BAframe['lat'])*(lat-BAframe['lat'])
                BAb2 = (lon-BAframe['lon'])*(lon-BAframe['lon'])
                BAc = np.sqrt(BAa2+BAb2)/2
                bap = ax2.plot(BAc, BAframe['depth'].values, 'r.', label='BA')
            if len(EQframe)>0:
                EQa2 = (lat-EQframe['lat'])*(lat-EQframe['lat'])
                EQb2 = (lon-EQframe['lon'])*(lon-EQframe['lon'])
                EQc = np.sqrt(EQa2+EQb2)/2
                eqp = ax2.plot(EQc, EQframe['depth'].values, 'c.', label='EQ')
            if len(ERframe)>0:
                ERa2 = (lat-ERframe['lat'])*(lat-ERframe['lat'])
                ERb2 = (lon-ERframe['lon'])*(lon-ERframe['lon'])
                ERc = np.sqrt(ERa2+ERb2)/2
                erp = ax2.plot(ERc, ERframe['depth'].values, 'y.', label='ER')
            if len(AAframe)>0:
                AAframe.loc[AAframe.lon < 0, 'lon']+=360
                AAa2 = (lat-AAframe['lat'])*(lat-AAframe['lat'])
                AAb2 = (lon-AAframe['lon'])*(lon-AAframe['lon'])
                AAc = np.sqrt(AAa2+AAb2)/2
                aap = ax2.plot(AAc, AAframe['depth'].values, 'k.', label='AA')
            if len(ASframe)>0:
                ASa2 = (lat-ASframe['lat'])*(lat-ASframe['lat'])
                ASb2 = (lon-ASframe['lon'])*(lon-ASframe['lon'])
                ASc = np.sqrt(ASa2+ASb2)/2
                asp = ax2.plot(ASc, ASframe['depth'].values, 'm.', label='AS')
            if len(TOframe)>0:
                TOa2 = (lat-TOframe['lat'])*(lat-TOframe['lat'])
                TOb2 = (lon-TOframe['lon'])*(lon-TOframe['lon'])
                TOc = np.sqrt(TOa2+TOb2)/2
                top = ax2.plot(TOc, TOframe['depth'].values, 'g.', label='TO')
            if len(RFframe)>0:
                RFa2 = (lat-RFframe['lat'])*(lat-RFframe['lat'])
                RFb2 = (lon-RFframe['lon'])*(lon-RFframe['lon'])
                RFc = np.sqrt(RFa2+RFb2)/2
                rfp = ax2.plot(RFc, RFframe['depth'].values, 'b.', label='RF')
            if len(CPframe)>0:
                CPa2 = (lat-CPframe['lat'])*(lat-CPframe['lat'])
                CPb2 = (lon-CPframe['lon'])*(lon-CPframe['lon'])
                CPc = np.sqrt(CPa2+CPb2)/2
                CPp = ax2.plot(CPc, CPframe['depth'].values, color='orange',marker='.', label='CP')
        
            if sdepth<0:
                sdepth *= -1
            ax2.plot((0.1, 0.1), (loc_depth-sdepth, ddepth+loc_depth), 'b-')
            ax2.plot((0, 0.2), (loc_depth-sdepth, loc_depth-sdepth), 'b-')
            rangep = ax2.plot((0, 0.2), (ddepth+loc_depth, ddepth+loc_depth), 'b-', label='depthrange')
            locp = ax2.plot((0, np.max(c)), (loc_depth, loc_depth), 'g-', label='Slab1')
            pdfp = ax2.plot(PDF_df['Summed_Values'].values, PDF_df['depths'].values, linewidth=2, color='k', label='PDF')
            pkp = ax2.plot([peak, peak], [loc_depth-sdepth, ddepth+loc_depth], 'r--')
            pkp = ax2.plot([0, 0.5], [peak_depth, peak_depth], 'r--', label='Peak Depth')
            x1, x2, y1, y2 = ax2.axis()
            xmax = max(np.max(c), peak)
            ax2.axis((0, xmax, y1, y2))
            ax2.invert_yaxis()
            ax2.set_xlabel('Probability (PDF) Degree Distance from Node/2 (data)')
            ax2.set_ylabel('Depth')
            title = 'Lat: %.4f, Lon: %.4f, NID: %.4f' % (lat, lon, nID)
            ax2.set_title(title)
            ax2.grid()
            plt.legend(loc='best')
            lontit = lon*100
            lattit = lat*100
            figtitle = 'Output/PDF%s/%spdf%i.png' % (slab, whichpdf, nID)
            #fig.savefig(figtitle)
            plt.close()
            filetitle = 'Output/PDF%s/%sused%i.csv' % (slab, whichpdf, nID)
            trimmed.to_csv(filetitle, header=True, index=False, float_format='%0.2f', na_rep = float('nan'))


    # If there is only one event, we do not solve for the depth at that point unless it is AA, BA, or AS
    elif len(elistBAe) > 0 or len(elistASe) > 0 or len(elistAAe) > 0 or len(elistRFe) > 0 or len(elistTOe) > 0 or len(elistCPe) > 0:
        frames = [elistBAe, elistASe, elistAAe, elistRFe, elistTOe, elistCPe]
        trimmed_once = pd.concat(frames)
        all_depths = trimmed_once['depth'].values
        variance1 = trimmed_once['unc'].values
        peak_depth = np.mean(all_depths)
        stdv = np.mean(variance1)
        test = True
        n = 1
    else:
        peak_depth = np.nan
        stdv = np.nan
        test = False
        n = 0
    # GLM 11.23.16
    try:
        return peak_depth, stdv, test, n, multipeaks, centsurf
    except:
        return peak_depth, stdv, test, n, multipeaks, stdv

def slabShift_noGMT(tmp_res, node, T, trenches, taper_depth, taper_width, ages, ages_error, filterwidth, slab, maxthickness, spacing, lonname, latname, depthname, fracS, nCores, meanBA, testprint, kdeg, knot_no, rbfs, use_box):
    
    tmpdata = np.zeros((len(tmp_res), 4))
    tmpdata[:, 0], tmpdata[:, 1] = tmp_res[lonname].values, tmp_res[latname].values
    try:
        tmpdata[:, 2], tmpdata[:, 3] = tmp_res['depth'].values, tmp_res['stdv'].values
    except:
        tmpdata[:, 2], tmpdata[:, 3] = tmp_res['depth'].values, tmp_res['unc'].values

    print ('        generating shifting surface ...')
    if slab == 'sum':
        Surfgrid, xi, dl = chunksurface(tmpdata, node, T, slab, spacing, 'depth', 'time', 'test.txt', filterwidth, pd.DataFrame(), nCores, trenches, meanBA, kdeg, knot_no, rbfs, tmp_res, 'shift', 'og','lon',100,110,105)
        flipornot = 'flip'
    elif slab == 'jap':
        Surfgrid, xi, dl = chunksurface(tmpdata, node, T, slab, spacing, 'depth', 'time', 'test.txt', filterwidth, pd.DataFrame(), nCores, trenches, meanBA, kdeg, knot_no, rbfs, tmp_res, 'shift', 'og','lat',30,40,35)
        flipornot = 'flip'
    else:
        Surfgrid, xi, dl = pySurface3(tmpdata, node, T, slab, spacing, 'depth', 'time', 'test.txt', filterwidth, pd.DataFrame(), nCores, trenches, meanBA, kdeg, knot_no, rbfs, tmp_res, 'shift', 'og')
        flipornot = 'dontflip'

    sigma = (3.0/2.0) / spacing
    if slab == 'mue':
        sigma = (1.0/2.0) / spacing
    filtshifted = ndimage.filters.gaussian_filter(Surfgrid, sigma)

    strgrid3, dipgrid3 = mkSDgrddata(xi, filtshifted, flipornot)
    resdata = np.zeros((len(xi),5))
    resdata[:,0] = xi[:,0]
    resdata[:,1] = xi[:,1]
    resdata[:,2] = filtshifted.flatten()
    resdata[:,3] = strgrid3.flatten()
    resdata[:,4] = dipgrid3.flatten()
    newres = mkContourClip(tmp_res, trenches, node, resdata, False,slab)
    if len(trenches)>0:
        clip = clippingmask(newres,trenches,node,False, slab, 'first')
    else:
        clip = noTrenchPolygon(newres, node, False, slab)

    mask = maskdatag(clip, xi)
    mask.shape = Surfgrid.shape
    filtshifted = (filtshifted*mask)
    strgrid, dipgrid = mkSDgrddata(xi, filtshifted, flipornot)

    Filtgrid = pd.DataFrame({'lon':xi[:, 0],'lat':xi[:, 1],'depth':filtshifted.flatten(),'strike':strgrid.flatten(),'dip':dipgrid.flatten()})

    Filtgrid = Filtgrid[(np.isfinite(Filtgrid.depth))&(np.isfinite(Filtgrid.strike))&(np.isfinite(Filtgrid.dip))]
    #Filtgrid.to_csv('shiftingsurface.csv',header=True,index=False,na_rep=np.nan)

    # Determine age of plate at trench
    if len(trenches)>0 and slab != 'helz':
        trench_age = np.zeros((len(trenches['lon'].values), 4))
        trench_age[:, 0], trench_age[:, 1] = trenches['lon'].values, trenches['lat'].values
        trench_age[:, 0][trench_age[:, 0]>180] -= 360

        for i in range(len(trench_age)):
            if trench_age[i, 0] > 179.8 and trench_age[i, 0] < 180.1: #GLM 11.30.16 eventually fix ages file
                trench_age[i, 0] = 179.8
            trench_age[i, 2] = ages.getValue(trench_age[i, 1], trench_age[i, 0])/100 # GLM 11.16.16 [i,2] instead of [:,2]
            trench_age[i, 3] = ages_error.getValue(trench_age[i, 1], trench_age[i, 0])/100

        trench_age[:,0][trench_age[:,0]<0] += 360

        trench_age[trench_age==327.] = np.nan  # Hardwired 327 because of default value of age grid file
        ta0, ta1, ta2, ta3 = trench_age[:, 0], trench_age[:, 1], trench_age[:, 2], trench_age[:, 3]
        ta0, ta1, ta2, ta3 = ta0[np.isfinite(ta3)], ta1[np.isfinite(ta3)], ta2[np.isfinite(ta3)], ta3[np.isfinite(ta3)]

        trench_age = np.zeros((len(ta0), 4))
        trench_age[:, 0], trench_age[:, 1], trench_age[:, 2], trench_age[:, 3] = ta0, ta1, ta2, ta3

    # Determine strike, dip, nearest trench for each input datum
    all_pts = np.zeros((len(tmp_res), 10))
    
    # Fill in lat,lon,dep from original data
    all_pts[:, 0], all_pts[:, 1], all_pts[:, 2] = tmp_res['bzlon'].values, tmp_res['bzlat'].values, tmp_res['depth'].values

    surfarr = np.zeros((len(Filtgrid), 4))
    #surfarr = np.zeros((len(Filtgrid),2))
    surfarr[:, 0] = Filtgrid['lon'].values
    surfarr[:, 1] = Filtgrid['lat'].values
    surfarr[:, 2] = Filtgrid['strike'].values
    surfarr[:, 3] = Filtgrid['dip'].values

    # Fill in plate age and error
    if len(trenches) > 0 and slab != 'helz':
        all_pts[:, 5] = griddata(trench_age[:, 0:2], trench_age[:, 2], all_pts[:, 0:2], method='nearest')
        all_pts[:, 8] = griddata(trench_age[:, 0:2], trench_age[:, 3], all_pts[:, 0:2], method='nearest')
    else:
        all_pts[:, 5] = 75
        all_pts[:, 8] = 20
    
    # Fill in strike and dip from original slab center surface
    all_pts[:, 3] = griddata(surfarr[:, 0:2], surfarr[:, 2], all_pts[:, 0:2], method='nearest')
    all_pts[:, 4] = griddata(surfarr[:, 0:2], surfarr[:, 3], all_pts[:, 0:2], method='nearest')
    
    # Calculating crustal thickness
    """ References
        thermal conductivity (lithosphere): k = 3.138 W/m C (Stein and Stein, 1996)
        specific heat (lithosphere): Cp = 1.171 kJ/kg C (Stein and Stein, 1996)
        density (mantle): rhom = 3330 kg/m^3 (Stein and Stein, 1996)
        thermal diffusivity (lithosphere): kappa = k/(Cp*rhom) = ~0.8E-6 m^2/s
        thickness (lithosphere): h = 2.32 * sqrt(kappa*age) where age is in seconds and h is in meters (Turcotte and Schubert)
            **base of lithosphere defined when (T-T1)/(T0-T1) = 0.1 (see T&S eqs. 4.93 and 4.115)
    """

    k = 3.138
    Cp = 1171.
    pm = 3330.
    kappa = k / (Cp*pm)  # For explanation see above

    new_pts = np.zeros((len(all_pts), 9))

    #centsurf = tmp_res['centsurf'].values
    for i in range(len(all_pts)):

        age_sec = all_pts[i, 5] * 1000000 * 365.25 * 24 * 60 * 60  # Convert age in Myr to age in seconds
        all_pts[i, 6] = 2.32 * math.sqrt(kappa * age_sec) / 1000  # Divide by 1000 converts from meters to kilometers - thickness
        if slab == 'hal' or slab == 'pam' or slab == 'hin' or slab == 'him':
            all_pts[i,6] = 100

    max_thickness = 5000

    if taper_width == max_thickness:
        taper_width = np.max(all_pts[:, 6])
    else:
        taper_width = taper_width
    
    maxthickness = np.max(all_pts[:, 6])
    all_pts[:,9] = tmp_res['onlyto'].values

    for i in range(len(all_pts)):
        error_sec = all_pts[i, 8] * 1000000 * 365.25 * 24 * 60 * 60  # Convert error to seconds

        if all_pts[i, 2] <= taper_depth:
            new_pts[i, 0:3] = all_pts[i, 0:3]
            all_pts[i, 7] = 0
            new_pts[i, 3] = 0  # Thickness error
            continue

        elif all_pts[i, 2] > taper_depth and all_pts[i, 2] < taper_depth+taper_width:
        
            x = taper_width/math.sin(np.radians(all_pts[i, 4]))
            dzs = abs(all_pts[i,2] - taper_depth)
            dxs = dzs/math.sin(np.radians(all_pts[i, 4]))
            taper = dxs/x
            if testprint:
                print (all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper)
            taper = dzs/(2*taper_width)
            if testprint:
                print ('all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper')
                print (all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper)

        else:
            taper = 1.0

        if all_pts[i,4] > 60 and slab != 'alu':
            all_pts[i,4] = 90
        if slab == 'man' and all_pts[i, 2] > 200:
            all_pts[i,4] = 90
        if all_pts[i, 9] == 1 and (slab == 'man' or slab == 'sam'):
            all_pts[i, 7] = (all_pts[i, 6]*fracS) * taper * 1.5
        else:
            all_pts[i, 7] = (all_pts[i, 6]*fracS) * taper

        if slab == 'muez':
            all_pts[i, 4] *= 1.5
            
        new_pts[i, 0], new_pts[i, 1], new_pts[i, 2] = pointShift(all_pts[i, 0], all_pts[i, 1], all_pts[i, 2], all_pts[i, 4], all_pts[i, 3], all_pts[i, 7])

        age_sec = all_pts[i, 5] * 1000000 * 365.25 * 24 * 60 * 60
        new_pts[i, 3] = math.sqrt(math.pow((2.32 * k * taper / (math.sqrt(kappa * age_sec) * Cp * pm * 1000. * 2. )), 2)*math.pow((error_sec/10.), 2) +
                              math.pow((2.32 * age_sec * taper / (math.sqrt(kappa * age_sec) * Cp * pm * 1000. * 2. )), 2)*math.pow((k/10.), 2) +
                              math.pow((-1. * 2.32 * k * age_sec * taper / (math.pow((kappa * age_sec), (3./2.)) * Cp * 1000. * 2. )), 2)*math.pow((pm/10.), 2) +
                              math.pow((-1. * 2.32 * k * age_sec * taper / (math.pow((kappa * age_sec), (3./2.)) * pm * 1000. * 2. )), 2)*math.pow((Cp/10.), 2))
        if testprint:
            print ('new_pts[i, 0], new_pts[i, 1], new_pts[i, 2]', new_pts[i, 0], new_pts[i, 1], new_pts[i, 2])
            print ('lon,lat,depth,strike,dip,thickness,taper,taper-depth,taper-width',all_pts[i,0],all_pts[i,1],all_pts[i,2],all_pts[i,3],all_pts[i,4],all_pts[i,7],taper,taper_depth,taper_width)

    new_pts[:, 4] = all_pts[:, 7]
    try:
        new_pts[:, 5] = tmp_res['nID'].values
    except:
        new_pts[:, 5] = tmp_res['ID'].values
    new_pts[:, 6] = all_pts[:, 2]
    new_pts[:, 7] = all_pts[:, 3]
    new_pts[:, 8] = all_pts[:, 4]
    
    new_pts[:, 0][new_pts[:, 0]<0] += 360

    shift_out = pd.DataFrame({'lon':new_pts[:, 0],'lat':new_pts[:, 1],'depth':new_pts[:, 2],'shiftstd':new_pts[:, 3],'smag':new_pts[:, 4],'nID':new_pts[:, 5].astype(int),'sstr':new_pts[:, 7],'sdip':new_pts[:, 8],'thickness':all_pts[:,6]})

    return shift_out, maxthickness

def newSlabShift(tmp_res, node, T, trenches, taper_depth, taper_width, ages, ages_error, filterwidth, slab, maxthickness, spacing, lonname, latname, depthname, fracS, fill_dat, meanBA, testprint, kdeg, knot_no, rbfs, use_box):
    
    tmp_res2 = pd.concat([tmp_res, fill_dat])
    tmpdata = np.zeros((len(tmp_res2), 4))
    tmpdata[:, 0], tmpdata[:, 1] = tmp_res2[lonname].values, tmp_res2[latname].values
    try:
        tmpdata[:, 2], tmpdata[:, 3] = tmp_res2['depth'].values, tmp_res2['stdv'].values
    except:
        tmpdata[:, 2], tmpdata[:, 3] = tmp_res2['depth'].values, tmp_res2['unc'].values

    print ('        generating shifting surface ...')

    Surfgrid, xi, dl = pySurface4(tmpdata, node, T, slab, spacing, 'depth', 'time', 'test.txt', filterwidth, pd.DataFrame(), 1, trenches, meanBA, kdeg, knot_no, rbfs, tmp_res, 'shift', 'og')
    flipornot = 'dontflip'

    sigma = (3.0/2.0) / spacing
    if slab == 'mue':
        sigma = (1.0/2.0) / spacing
    filtshifted = ndimage.filters.gaussian_filter(Surfgrid, sigma)

    strgrid, dipgrid = mkSDgrddata(xi, filtshifted, flipornot)

    Filtgrid = pd.DataFrame({'lon':xi[:, 0],'lat':xi[:, 1],'depth':filtshifted.flatten(),'strike':strgrid.flatten(),'dip':dipgrid.flatten()})

    Filtgrid = Filtgrid[(np.isfinite(Filtgrid.depth))&(np.isfinite(Filtgrid.strike))&(np.isfinite(Filtgrid.dip))]
    #Filtgrid.to_csv('shiftingsurface.csv',header=True,index=False,na_rep=np.nan)

    # Determine age of plate at trench
    if len(trenches)>0 and slab != 'helz':
        trench_age = np.zeros((len(trenches['lon'].values), 4))
        trench_age[:, 0], trench_age[:, 1] = trenches['lon'].values, trenches['lat'].values
        trench_age[:, 0][trench_age[:, 0]>180] -= 360

        for i in range(len(trench_age)):
            if trench_age[i, 0] > 179.8 and trench_age[i, 0] < 180.1: #GLM 11.30.16 eventually fix ages file
                trench_age[i, 0] = 179.8
            trench_age[i, 2] = ages.getValue(trench_age[i, 1], trench_age[i, 0])/100 # GLM 11.16.16 [i,2] instead of [:,2]
            trench_age[i, 3] = ages_error.getValue(trench_age[i, 1], trench_age[i, 0])/100

        trench_age[:,0][trench_age[:,0]<0] += 360

        trench_age[trench_age==327.] = np.nan  # Hardwired 327 because of default value of age grid file
        ta0, ta1, ta2, ta3 = trench_age[:, 0], trench_age[:, 1], trench_age[:, 2], trench_age[:, 3]
        ta0, ta1, ta2, ta3 = ta0[np.isfinite(ta3)], ta1[np.isfinite(ta3)], ta2[np.isfinite(ta3)], ta3[np.isfinite(ta3)]

        trench_age = np.zeros((len(ta0), 4))
        trench_age[:, 0], trench_age[:, 1], trench_age[:, 2], trench_age[:, 3] = ta0, ta1, ta2, ta3

    # Determine strike, dip, nearest trench for each input datum
    all_pts = np.zeros((len(tmp_res), 10))
    
    # Fill in lat,lon,dep from original data
    all_pts[:, 0], all_pts[:, 1], all_pts[:, 2] = tmp_res['bzlon'].values, tmp_res['bzlat'].values, tmp_res['depth'].values

    surfarr = np.zeros((len(Filtgrid), 4))
    #surfarr = np.zeros((len(Filtgrid),2))
    surfarr[:, 0] = Filtgrid['lon'].values
    surfarr[:, 1] = Filtgrid['lat'].values
    surfarr[:, 2] = Filtgrid['strike'].values
    surfarr[:, 3] = Filtgrid['dip'].values

    # Fill in plate age and error
    if len(trenches) > 0 and slab != 'helz':
        all_pts[:, 5] = griddata(trench_age[:, 0:2], trench_age[:, 2], all_pts[:, 0:2], method='nearest')
        all_pts[:, 8] = griddata(trench_age[:, 0:2], trench_age[:, 3], all_pts[:, 0:2], method='nearest')
    else:
        all_pts[:, 5] = 75
        all_pts[:, 8] = 20
    
    # Fill in strike and dip from original slab center surface
    all_pts[:, 3] = griddata(surfarr[:, 0:2], surfarr[:, 2], all_pts[:, 0:2], method='nearest')
    all_pts[:, 4] = griddata(surfarr[:, 0:2], surfarr[:, 3], all_pts[:, 0:2], method='nearest')
    
    # Calculating crustal thickness
    """ References
        thermal conductivity (lithosphere): k = 3.138 W/m C (Stein and Stein, 1996)
        specific heat (lithosphere): Cp = 1.171 kJ/kg C (Stein and Stein, 1996)
        density (mantle): rhom = 3330 kg/m^3 (Stein and Stein, 1996)
        thermal diffusivity (lithosphere): kappa = k/(Cp*rhom) = ~0.8E-6 m^2/s
        thickness (lithosphere): h = 2.32 * sqrt(kappa*age) where age is in seconds and h is in meters (Turcotte and Schubert)
            **base of lithosphere defined when (T-T1)/(T0-T1) = 0.1 (see T&S eqs. 4.93 and 4.115)
    """

    k = 3.138
    Cp = 1171.
    pm = 3330.
    kappa = k / (Cp*pm)  # For explanation see above

    new_pts = np.zeros((len(all_pts), 9))

    #centsurf = tmp_res['centsurf'].values
    for i in range(len(all_pts)):

        age_sec = all_pts[i, 5] * 1000000 * 365.25 * 24 * 60 * 60  # Convert age in Myr to age in seconds
        all_pts[i, 6] = 2.32 * math.sqrt(kappa * age_sec) / 1000  # Divide by 1000 converts from meters to kilometers - thickness
        if slab == 'hal' or slab == 'pam' or slab == 'hin' or slab == 'him':
            all_pts[i,6] = 100

    max_thickness = 5000

    if taper_width == max_thickness:
        taper_width = np.max(all_pts[:, 6])
    else:
        taper_width = taper_width

    maxthickness = np.max(all_pts[:, 6])
    all_pts[:,9] = tmp_res['onlyto'].values

    for i in range(len(all_pts)):
        error_sec = all_pts[i, 8] * 1000000 * 365.25 * 24 * 60 * 60  # Convert error to seconds

        if all_pts[i, 2] <= taper_depth:
            new_pts[i, 0:3] = all_pts[i, 0:3]
            all_pts[i, 7] = 0
            new_pts[i, 3] = 0  # Thickness error
            continue

        elif all_pts[i, 2] > taper_depth and all_pts[i, 2] < taper_depth+taper_width:
        
            x = taper_width/math.sin(np.radians(all_pts[i, 4]))
            dzs = abs(all_pts[i,2] - taper_depth)
            dxs = dzs/math.sin(np.radians(all_pts[i, 4]))
            taper = dxs/x
            if testprint:
                print (all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper)
            taper = dzs/(2*taper_width)
            if testprint:
                print ('all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper')
                print (all_pts[i,2],all_pts[i,4],x,dxs,dzs,taper_depth-taper_width,taper_depth+taper_width,taper*2,taper_depth,taper)

        else:
            taper = 1.0

        if all_pts[i,4] > 60 and slab != 'alu':
            all_pts[i,4] = 90
        if slab == 'man' and all_pts[i, 2] > 200:
            all_pts[i,4] = 90
        if all_pts[i, 9] == 1 and (slab == 'man' or slab == 'sam'):
            all_pts[i, 7] = (all_pts[i, 6]*fracS) * taper * 1.5
        else:
            all_pts[i, 7] = (all_pts[i, 6]*fracS) * taper

        if slab == 'muez':
            all_pts[i, 4] *= 1.5
        
        new_pts[i, 0], new_pts[i, 1], new_pts[i, 2] = pointShift(all_pts[i, 0], all_pts[i, 1], all_pts[i, 2], all_pts[i, 4], all_pts[i, 3], all_pts[i, 7])

        age_sec = all_pts[i, 5] * 1000000 * 365.25 * 24 * 60 * 60
        new_pts[i, 3] = math.sqrt(math.pow((2.32 * k * taper / (math.sqrt(kappa * age_sec) * Cp * pm * 1000. * 2. )), 2)*math.pow((error_sec/10.), 2) +
                              math.pow((2.32 * age_sec * taper / (math.sqrt(kappa * age_sec) * Cp * pm * 1000. * 2. )), 2)*math.pow((k/10.), 2) +
                              math.pow((-1. * 2.32 * k * age_sec * taper / (math.pow((kappa * age_sec), (3./2.)) * Cp * 1000. * 2. )), 2)*math.pow((pm/10.), 2) +
                              math.pow((-1. * 2.32 * k * age_sec * taper / (math.pow((kappa * age_sec), (3./2.)) * pm * 1000. * 2. )), 2)*math.pow((Cp/10.), 2))
        if testprint:
            print ('new_pts[i, 0], new_pts[i, 1], new_pts[i, 2]', new_pts[i, 0], new_pts[i, 1], new_pts[i, 2])
            print ('lon,lat,depth,strike,dip,thickness,taper,taper-depth,taper-width',all_pts[i,0],all_pts[i,1],all_pts[i,2],all_pts[i,3],all_pts[i,4],all_pts[i,7],taper,taper_depth,taper_width)

    new_pts[:, 4] = all_pts[:, 7]
    try:
        new_pts[:, 5] = tmp_res['nID'].values
    except:
        new_pts[:, 5] = tmp_res['ID'].values
    new_pts[:, 6] = all_pts[:, 2]
    new_pts[:, 7] = all_pts[:, 3]
    new_pts[:, 8] = all_pts[:, 4]

    new_pts[:, 0][new_pts[:, 0]<0] += 360

    shift_out = pd.DataFrame({'lon':new_pts[:, 0],'lat':new_pts[:, 1],'depth':new_pts[:, 2],'shiftstd':new_pts[:, 3],'smag':new_pts[:, 4],'nID':new_pts[:, 5].astype(int),'sstr':new_pts[:, 7],'sdip':new_pts[:, 8],'thickness':all_pts[:,6]})

    return shift_out, maxthickness


## Written GLM 11.21.2016
def az_perp(x):
    ''' Arguments:  x - azimuth
        
        Returns:    x - the input azimuth - 90 degrees (azimuth oriented outboard the trench) '''
    if x<=90:
        return x+270
    else:
        return x-90

def az_other_perp(x):
    if x >= 270:
        return x - 270
    else:
        return x + 90

def npaz_perp(x):
    ''' Arguments:  x - azimuth
        
        Returns:    x - the input azimuth - 90 degrees (azimuth oriented outboard the trench) '''
    x += 270
    x[x > 360] -= 360
    return x

def npaz_other_perp(x):

    x -= 270
    x[x  < 0] += 360
    return x

###############################################

### 48 ###

###############################################

## Written GLM 11.21.2016
def uncraise(x):
    ''' Arguments:  x - uncertainty
        
        Returns:    x - raised to the minimum uncertainty '''
    minunc = 15.0
    if x<minunc:
        return minunc
    else:
        return x


def distchecker(maskdata, distcheck):

    for i in range(1, len(maskdata)-1):
        lastdist = distcheck[i-1]
        thisdist = distcheck[i]
        nextdist = distcheck[i+1]

        if thisdist < lastdist and thisdist < nextdist:
            maskdata = maskdata[maskdata.dists != thisdist]
    
    return maskdata


def movingav2(x,y,testprint,filtdist):
    x2 = np.copy(x)
    y2 = np.copy(y)
    n=0
    for i in range(1,len(x)-1):
        thisx = x[i]
        lastx = x[i-1]
        nextx = x[i+1]
        thisy = y[i]
        lasty = y[i-1]
        nexty = y[i+1]
        lastdiff = (thisx-lastx)*(thisx-lastx)+(thisy-lasty)*(thisy-lasty)
        nextdiff = (thisx-nextx)*(thisx-nextx)+(thisy-nexty)*(thisy-lasty)
        outdiff = (nextx-lastx)*(nextx-lastx)+(nexty-lasty)*(nexty-lasty)
        if outdiff<lastdiff*filtdist or outdiff<nextdiff*filtdist:
            x2[i] = (nextx+lastx)/2.0
            y2[i] = (nexty+lasty)/2.0
            if testprint:
                print ('dropped outdiff,lastdiff,nextdiff,lastx,thisx,nextx,lasty,thisy,nexty',outdiff,lastdiff,nextdiff,lastx,thisx,nextx,lasty,thisy,nexty)
            n+=1
        else:
            x2[i] = thisx
            y2[i] = thisy
    x2[0] = x[0]
    y2[0] = y[0]
    x2[-1] = x[-1]
    y2[-1] = y[-1]
    return x2,y2,n

def mkClipPolygon(indataOG, trench, spacing, testprint):
    
    indatadat = np.zeros((len(indataOG),3)).astype(np.float64)
    indatadat[:,0] = np.round(indataOG['lon'].values, 2)
    indatadat[:,1] = np.round(indataOG['lat'].values, 2)
    indatadat[:,2] = np.round(indataOG['depth'].values, 2)
    
    #print (indatadat)
    
    slab = trench['slab'].values[0]
    
    if slab == 'sulz' or slab == 'cotz':
        dw = 0.3
    elif slab != 'phi':
        dw = spacing*1.5
    else:
        dw = spacing
    
    dl = 20.0
    if slab == 'sol' or slab == 'sul':
        dist = 150.0
    elif slab == 'hel':
        dist = 25
    else:
        dist = 75.0
    idnot = (np.ones(len(indataOG))*-9999).astype(int)

    if slab == 'sum':
        minempty = 0
        mindepth = 100.0
        distthresh = 800
    elif slab == 'hel':
        minempty = 0
        mindepth = 70.0
        distthresh = 400
    elif slab == 'kur':
        minempty = 0
        mindepth = 100.0
        distthresh = 800
    elif slab == 'alu':
        mindepth = 10.0
        minempty = 0
        distthresh = 800
    elif slab == 'sol' or slab == 'png':
        mindepth = 10.0
        distthresh = 500
        minempty = 0
    elif slab == 'sul' or slab == 'cot':
        mindepth = 10.0
        distthresh = 300
        minempty = 0
    elif slab == 'phi' or slab == 'man':
        mindepth = 30.0
        minempty = 0
        distthresh = 800
    elif slab == 'cam':
        mindepth = 40
        minempty = 0
        distthresh = 800
    elif slab == 'cam':
        mindepth = 40
        minempty = 0
        distthresh = 1100
    elif slab == 'cas':
        mindepth = 40
        minempty = 0
        distthresh = 800
    else:
        minempty = 0
        mindepth = 25.0
        distthresh = 800
    # make trench-side of clipping mask
    trench.loc[trench.lon<0,'lon'] += 360
    trench['az90'] = npaz_perp(trench['az'].values*1.0)
    dists = np.ones(len(trench))*dist
    tlons = trench['lon'].values*1.0
    tlats = trench['lat'].values*1.0
    #lon90, lat90 = npheading(tlons,tlats,az90,dists)
    lon90, lat90=zip(*trench.apply(lambda row: heading(row['lon'], row['lat'], dist, row['az90']), axis=1))
    
    masktrench = pd.DataFrame({'lon':lon90,'lat':lat90})
    idlist = []
    for i in range (len(indatadat)):
        nodelon = indatadat[i,0]
        nodelat = indatadat[i,1]
        nodedepth = indatadat[i,2]

        if nodedepth < mindepth:
            idnot[i] = i
        elif slab == 'sam' and nodelon < 286 and nodelat > -30:
            idnot[i] = i
        elif slab == 'sam' and nodelon < 287.5 and nodelon > 285 and nodelat < -40:
            idnot[i] = i
        elif slab == 'sam' and nodelon < 287.5 and nodelon > 285 and nodelat < 10 and nodelat > 5:
            idnot[i] = i
        elif slab == 'sam' and nodelat > -40 and nodelat < 2 and nodedepth < 300:
            idnot[i] = i
        elif slab == 'png' and nodelon > 135 and nodelon < 140 and nodedepth<100:
            idnot[i] = i
        elif slab == 'hal' and nodelat > 8 and nodedepth<500:
            idnot[i] = i
        elif slab == 'sol' and nodelat > -8 and nodelat<-6.5 and nodelon>147 and nodelon<150:
            idnot[i] = i
        elif slab == 'ker' and nodelat > -25 and nodelat < -23 and nodelon <178:
            idnot[i] = i
        elif slab == 'kur' and nodelat > 45 and nodelat < 50 and nodelon <142 and nodelon > 140:
            idnot[i] = i
        elif slab == 'kur' and nodelat > 42 and nodelat < 45 and nodelon <132:
            idnot[i] = i
        elif slab == 'kur' and nodelat > 54 and nodedepth < 500 and nodelon <156:
            idnot[i] = i
        elif slab == 'cam' and nodelon > 258 and nodelon <261 and nodelat < 19.8:
            idnot[i] = i
        elif slab == 'cam' and nodelon > 262 and nodedepth < 300 and nodelon <265:
            idnot[i] = i
        elif slab == 'ita' and nodelon > 13 and nodelat < 40.5 and nodelat>40  and nodelon <14:
            idnot[i] = i
        elif slab == 'png' and nodelon > 132 and nodelat> -2  and nodelon <135:
            idnot[i] = i
        elif slab == 'puy' and nodelon > 165 and nodedepth< 150  and nodelon <166:
            idnot[i] = i
        elif slab == 'mak' and nodelon < 65 and nodelat <29.4:
            idnot[i] = i
        elif slab == 'mak' and nodelat <28 and nodelat > 65:
            idnot[i] = i
        elif slab == 'hel' and nodelon > 30 and nodelat>36 and nodedepth < 400:
            idnot[i] = i
        elif slab == 'hel' and nodelon < 25 and nodelat<39:
            idnot[i] = i
        elif slab == 'hel' and nodelon < 26 and nodelat<38:
            idnot[i] = i
        elif slab == 'sul' and nodelon < 122.5 and nodelon > 120 and nodedepth<200:
            idnot[i] = i
        elif slab == 'cas' and nodelon > 235 and nodelat < 50 and nodelat > 45 and nodedepth < 200:
            idnot[i] = i
        elif slab == 'alu' and nodelon > 199 and nodelon < 201 and nodedepth < 220:
            idnot[i] = i
        elif slab == 'alu' and nodelon > 168 and nodelon < 172:
            idnot[i] = i
        elif slab == 'man' and nodelat > 15 and nodelat < 20 and nodedepth < 150:
            idnot[i] = i
        
    idnot = idnot[idnot>-999]
    notbytrench = np.delete(indatadat, idnot, 0)
    
    lons = np.ones(len(notbytrench))*-9999
    lats = np.ones(len(notbytrench))*-9999
    northlist = np.ones(len(notbytrench))*-9999
    eastlist = np.ones(len(notbytrench))*-9999
    southlist = np.ones(len(notbytrench))*-9999
    westlist = np.ones(len(notbytrench))*-9999

    lonEmin = 999
    lonEmax = -999
    latEmin = 999
    latEmax = -999

    for i in range(len(notbytrench)):
        dw1 = dw
        if slab == 'sum' and nodelat > 12:
            dw1 = 0.1
        if slab == 'cam' and nodelon > 270:
            dw1 = 0.15
        if slab == 'png' and nodelon > 136:
            dw1 = 0.1
        if slab == 'png' and nodelon > 137.5:
            dw1 = 0.3
        if slab == 'hel' and nodelon < 22.5:
            dw1 = 0.1
        if slab == 'mak' and nodelat < 29:
            dw1 = 0.3
            
        nodelon, nodelat = notbytrench[i,0], notbytrench[i,1]
        NS = indatadat[(indatadat[:,0] < nodelon+dw1) & (indatadat[:,0] > nodelon-dw1)]
        EW = indatadat[(indatadat[:,1] < nodelat+dw1) & (indatadat[:,1] > nodelat-dw1)]
        north = NS[(NS[:,1] > nodelat) & (NS[:,1] < nodelat+dl)]
        south = NS[(NS[:,1] < nodelat) & (NS[:,1] > nodelat-dl)]
        east = EW[(EW[:,0] > nodelon) & (EW[:,0] < nodelon+dl)]
        west = EW[(EW[:,0] < nodelon) & (EW[:,0] > nodelon-dl)]
        
        n = 0
        if len(north) < 1:
            n += 1
            northlist[i] = 1
        else:
            northlist[i] = 0
        if len(south) < 1 and slab != 'mak' and slab != 'him':
            n += 1
            southlist[i] = 1
        else:
            southlist[i] = 0
        if len(east) < 1 and slab != 'mak':
            n += 1
            eastlist[i] = 1
        else:
            eastlist[i] = 0
        if len(west) < 1 and slab !='him':
            n += 1
            westlist[i] = 1
        else:
            westlist[i] = 0

        if n > minempty:
            lons[i] = nodelon
            lats[i] = nodelat
        elif slab == 'sam' and n>0 and nodelon > 289 and nodelat > -14:
            lons[i] = nodelon
            lats[i] = nodelat
        elif slab == 'kur' and n>0 and nodelat > 55:
            lons[i] = nodelon
            lats[i] = nodelat

    lonbool = lons > -999
    maskN = southlist == 0
    maskE = westlist == 0
    maskS = northlist == 0
    maskW = eastlist == 0

    northlist = northlist[lonbool]
    eastlist = eastlist[lonbool]
    southlist = southlist[lonbool]
    westlist = westlist[lonbool]
    
    lons = lons[lons>-999]
    lats = lats[lats>-999]
    
    trenchtest = masktrench[(masktrench.lat<=np.max(lats))&(masktrench.lat>=np.min(lats))&(masktrench.lat<=np.max(lons))&(masktrench.lat<=np.max(lats))]
    addfirst = masktrench.iloc[[0]]
    lastpoint = masktrench.iloc[[-1]]
    lastlon = lastpoint['lon'].values[0]
    lastlat = lastpoint['lat'].values[0]
    lastN,lastE,lastS,lastW = 1,1,1,1
    sortedlons = np.ones(len(lons))*-9999
    sortedlats = np.ones(len(lats))*-9999
    sortedangs = np.ones(len(lats))*-9999
    gotOne = True
    alons = np.array(lons)
    alats = np.array(lats)
    
    awest = np.array(westlist)
    aeast = np.array(eastlist)
    anorth = np.array(northlist)
    asouth = np.array(southlist)

    presort = pd.DataFrame({'lon':lons, 'lat':lats, 'depth':1})
    #presort.to_csv('cliptest.csv',header=True,index=False,na_rep=np.nan)

    if slab == 'car':
        lastlon = 290
        lastlat = 20
    
    n = 0
    while gotOne == True:
        dists, angs = npcosine(lastlon, lastlat, alons, alats)
        
        if n>1:
            if lastN == 1:
                maskN = asouth == 0
            else:
                maskN = np.ones(len(dists), dtype=bool)
            if lastE == 1:
                maskE = awest == 0
            else:
                maskE = np.ones(len(dists), dtype=bool)
            if lastS == 1:
                maskS = anorth == 0
            else:
                maskS = np.ones(len(dists), dtype=bool)
            if lastW == 1:
                maskW = aeast == 0
            else:
                maskW = np.ones(len(dists), dtype=bool)

            distsT = dists[maskN & maskE & maskS & maskW]
        
        if len(dists)>0:
            #print ('len(dists)',len(dists))
            #print ('dists,angs,alons,alats,lastlon,lastlat',dists,angs,alons,alats,lastlon,lastlat)
            
            if n>1 and len(distsT)>0 and slab != 'cot':
                minT = np.min(distsT)
                imindista = np.where(dists == minT)
                imindist = imindista[0][0]
            else:
                imindist = np.argmin(dists)

            if dists[imindist] < distthresh or n == 0:
            
                lastE, lastW = aeast[imindist], awest[imindist]
                lastN, lastS = anorth[imindist], asouth[imindist]
                
                lastlon, lastlat = alons[imindist], alats[imindist]
                lastang = angs[imindist]
                
                sortedlons[n] = lastlon
                sortedlats[n] = lastlat
                sortedangs[n] = lastang
                
                alons = np.delete(alons, imindist)
                alats = np.delete(alats, imindist)

                anorth = np.delete(anorth, imindist)
                aeast = np.delete(aeast, imindist)
                asouth = np.delete(asouth, imindist)
                awest = np.delete(awest, imindist)
                
                n+=1
            else:
                gotOne = False
        else:
            gotOne = False
    sortedlons = sortedlons[sortedlons>-999]
    sortedlats = sortedlats[sortedlats>-999]
    sortedangs = sortedlats[sortedlats>-999]
    #print ('sortedlons,sortedlats',sortedlons,sortedlats)

    if slab == 'mak':
        maskdata = pd.DataFrame({'lon':lons,'lat':lats})
        maskdata = maskdata.sort_values(by=['lon'],ascending = False)
    elif slab == 'cas':
        maskdata = pd.DataFrame({'lon':lons,'lat':lats})
        maskdata = maskdata.sort_values(by=['lat'],ascending = True)
    else:
        maskdata = pd.DataFrame({'lon':sortedlons,'lat':sortedlats})

    filtno = 10
    filtnum = 0
    n2 = 1

    while n2>0:
        maskdata['lon'], maskdata['lat'], n2 = movingav2(maskdata['lon'].values, maskdata['lat'].values,testprint,1)
        filtnum += 1
        print (filtnum)

    maskdata = maskdata[['lon', 'lat']]
    maskdata = maskdata.reset_index(drop=True)

    if slab == 'cot':
        clip = pd.concat([masktrench[:-3], maskdata, addfirst])
    else:
        clip = pd.concat([masktrench, maskdata, addfirst])

    return clip


def noTrenchPolygon(indataOG, spacing, testprint, slab):

    if slab == 'hal':
        indataOG = indataOG.sort_values(by=['lat'],ascending=False)
        toplat = indataOG['lat'].max()
        toplon = indataOG[indataOG.lat == toplat]['lon'].values[0]
    
    indatadat = np.zeros((len(indataOG),3)).astype(np.float64)
    indatadat[:,0] = np.round(indataOG['lon'].values,decimals=2)
    indatadat[:,1] = np.round(indataOG['lat'].values,decimals=2)
    indatadat[:,2] = np.round(indataOG['depth'].values,decimals=2)
    
    dw = spacing*1.5
    dl = 20.0
    dist = 75.0

    idnot = (np.ones(len(indataOG))*-9999).astype(int)
    minempty = 0
    if slab == 'hal':
        mindepth = 50.0
        distthresh = 8000
    else:
        mindepth = 3
        distthresh = 8000
    idlist = []
    for i in range (len(indatadat)):
        nodelon = indatadat[i,0]
        nodelat = indatadat[i,1]
        nodedepth = indatadat[i,2]
        if nodedepth < mindepth:
            idnot[i] = i

    idnot = idnot[idnot>-999]
    notbytrench = np.delete(indatadat, idnot, 0)
    
    lons = np.ones(len(notbytrench))*-9999
    lats = np.ones(len(notbytrench))*-9999
    northlist = np.ones(len(notbytrench))*-9999
    eastlist = np.ones(len(notbytrench))*-9999
    southlist = np.ones(len(notbytrench))*-9999
    westlist = np.ones(len(notbytrench))*-9999

    lonEmin = 999
    lonEmax = -999
    latEmin = 999
    latEmax = -999

    for i in range(len(notbytrench)):
        dw1 = dw
        dl1 = dl
        dw2 = spacing
        dl2 = 1.0
            
        nodelon, nodelat = notbytrench[i,0], notbytrench[i,1]
        NS = indatadat[(indatadat[:,0] < nodelon+dw1) & (indatadat[:,0] > nodelon-dw1)]
        EW = indatadat[(indatadat[:,1] < nodelat+dw1) & (indatadat[:,1] > nodelat-dw1)]
        north = NS[(NS[:,1] > nodelat) & (NS[:,1] < nodelat+dl1)]
        south = NS[(NS[:,1] < nodelat) & (NS[:,1] > nodelat-dl1)]
        east = EW[(EW[:,0] > nodelon) & (EW[:,0] < nodelon+dl1)]
        west = EW[(EW[:,0] < nodelon) & (EW[:,0] > nodelon-dl1)]
        
        n = 0
        if len(north) < 1:
            NS = indatadat[(indatadat[:,0] < nodelon+dw2) & (indatadat[:,0] > nodelon-dw2)]
            north = NS[(NS[:,1] > nodelat+dl2) & (NS[:,1] < nodelat+dl1)]
            if len(north) < 1:
                n += 1
                northlist[i] = 1
        else:
            northlist[i] = 0

        if len(south) < 1:
            NS = indatadat[(indatadat[:,0] < nodelon+dw2) & (indatadat[:,0] > nodelon-dw2)]
            south = NS[(NS[:,1] < nodelat-dl2) & (NS[:,1] > nodelat-dl1)]
            if len(south) < 1:
                n += 1
                southlist[i] = 1
        else:
            southlist[i] = 0

        if len(east) < 1:
            EW = indatadat[(indatadat[:,1] < nodelat+dw2) & (indatadat[:,1] > nodelat-dw2)]
            east = EW[(EW[:,0] > nodelon+dl2) & (EW[:,0] < nodelon+dl1)]
            if len(east) < 1:
                n += 1
                eastlist[i] = 1
        else:
            eastlist[i] = 0

        if len(west) < 1:
            EW = indatadat[(indatadat[:,1] < nodelat+dw2) & (indatadat[:,1] > nodelat-dw2)]
            west = EW[(EW[:,0] < nodelon-dl2) & (EW[:,0] > nodelon-dl1)]
            if len(west) < 1:
                n += 1
                westlist[i] = 1
        else:
            westlist[i] = 0

        if n > minempty:
            lons[i] = nodelon
            lats[i] = nodelat

            if slab == 'hin':
                northlist[i] = 1
                southlist[i] = 1
                eastlist[i] = 1
                westlist[i] = 1

    lonbool = lons > -999
    maskN = southlist == 0
    maskE = westlist == 0
    maskS = northlist == 0
    maskW = eastlist == 0

    northlist = northlist[lonbool]
    eastlist = eastlist[lonbool]
    southlist = southlist[lonbool]
    westlist = westlist[lonbool]
    
    lons = lons[lonbool]
    lats = lats[lonbool]

    lastlon = lons[0]
    lastlat = lats[0]
    firstlon = lons[0]
    firstlat = lats[0]
    lastN,lastE,lastS,lastW = 1,1,1,1
    sortedlons = np.ones(len(lons))*-9999
    sortedlats = np.ones(len(lats))*-9999
    sortedangs = np.ones(len(lats))*-9999
    gotOne = True
    alons = np.array(lons)
    alats = np.array(lats)
    
    awest = np.array(westlist)
    aeast = np.array(eastlist)
    anorth = np.array(northlist)
    asouth = np.array(southlist)

    n = 0
    while gotOne == True:
        dists, angs = npcosine(lastlon, lastlat, alons, alats)
        distf,angf,lonf,latf = cosine(lastlon,lastlat,firstlon,firstlat)
        
        if n>1:
            if lastN == 1:
                maskN = asouth == 0
            else:
                maskN = np.ones(len(dists), dtype=bool)
            if lastE == 1:
                maskE = awest == 0
            else:
                maskE = np.ones(len(dists), dtype=bool)
            if lastS == 1:
                maskS = anorth == 0
            else:
                maskS = np.ones(len(dists), dtype=bool)
            if lastW == 1:
                maskW = aeast == 0
            else:
                maskW = np.ones(len(dists), dtype=bool)

            distsT = dists[maskN & maskE & maskS & maskW]
        
        if len(dists)>0:
            #print (lastlon,lastlat,firstlon,firstlat,distf,np.min(dists))
            if np.min(dists) > distf*0.75 and n > 20:
                gotOne = False
                break
            
            if n>1 and len(distsT)>0:
                minT = np.min(distsT)
                imindista = np.where(dists == minT)
                imindist = imindista[0][0]
            else:
                imindist = np.argmin(dists)
            
            if dists[imindist] < distthresh or n == 0:
            
                lastE, lastW = aeast[imindist], awest[imindist]
                lastN, lastS = anorth[imindist], asouth[imindist]
                            
                lastlon, lastlat = alons[imindist], alats[imindist]
                lastang = angs[imindist]
                
                sortedlons[n] = lastlon
                sortedlats[n] = lastlat
                sortedangs[n] = lastang
                
                alons = np.delete(alons, imindist)
                alats = np.delete(alats, imindist)

                anorth = np.delete(anorth, imindist)
                aeast = np.delete(aeast, imindist)
                asouth = np.delete(asouth, imindist)
                awest = np.delete(awest, imindist)
                
                n+=1
            else:
                gotOne = False
        else:
            gotOne = False

    sortedlons = sortedlons[sortedlons>-999]
    sortedlats = sortedlats[sortedlats>-999]
    sortedangs = sortedlats[sortedlats>-999]
    #print ('sortedlons,sortedlats',sortedlons,sortedlats)

    maskdata = pd.DataFrame({'lon':sortedlons,'lat':sortedlats})

    filtno = 10
    filtnum = 0
    n2 = 1

    while n2>0 and slab != 'hin':
        maskdata['lon'], maskdata['lat'], n2 = movingav2(maskdata['lon'].values, maskdata['lat'].values,testprint,3)
        filtnum += 1

    maskdata = maskdata[['lon', 'lat']]
    maskdata = maskdata.reset_index(drop=True)

    maskdata.loc[len(maskdata)+1] = ([lons[0],lats[0]])

    #if slab == 'hal':
    #    maskdata = maskdata[(maskdata.lat < 5)|(maskdata.lon < 126)]
    #maskdata.to_csv('testingmask.csv',header=True,index=False)
    
    return maskdata
    
###############################################

### 50 ###

###############################################

## Written MF

def assignNodes(nodes, num_ranks, your_rank):

    len_data = len(nodes)
    q, r = divmod(len_data, num_ranks)
    my_nodes = (q*your_rank+min(your_rank, r), q*(your_rank+1)+min(your_rank+1, r))

    return my_nodes

###############################################

### 51 ###

###############################################

## Written GLM 11.22.16

def raiseEQunc(eventlistALL):
    ## GLM 11.21.16 raising earthquake uncertainties to a minimum of 10 km
    eventlist1 = eventlistALL[eventlistALL.etype != 'EQ']
    eventlist2 = eventlistALL[eventlistALL.etype == 'EQ']
    eventlist2['unc']=eventlist2.apply(lambda row: uncraise(row['unc']), axis=1)
    frames = [eventlist1, eventlist2]
    eventlistALL = pd.concat(frames)
    eventlist = eventlistALL.loc[eventlistALL.etype == 'RF', 'unc'] = 5
    eventlistALL = eventlistALL.reset_index(drop=True)
    eventlist = eventlistALL[['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', 'S1', 'D1', 'R1', 'S2', 'D2', 'R2', 'src']]
    return eventlist

###############################################

### 52 ###

###############################################

def getNearbyNodes(lat, lon, radius, eventlist):
    
    ''' Arguments:  lat - latitude of grid node that is being searched over
                    lon - longitude of grid node that is being searched over
                    radius - radius of circle to search within (km)
                    eventlist - list of events to search over. Must include lat/lon info
        
        Returns:    elist - dataframe of events that are within the specified radius of
                            the lat lon node point. Has all information that the original
                            eventlist contained for nearby data points of all types. '''
    
    # Gather latitudes and longitudes for each point in eventlist
    elat = eventlist['lat'].as_matrix()
    elon = eventlist['lon'].as_matrix()
    # Make array of length eventlist with distances from each lat/lon point from node
    a = abs(lat-eventlist['lat'])
    a2 = (lat-eventlist['lat'])*(lat-eventlist['lat'])
    b2 = (lon-eventlist['lon'])*(lon-eventlist['lon'])
    c = np.sqrt(a2+b2)
    pd.options.mode.chained_assignment = None
    # Filter points from eventlist that are outside of the radius
    eventlist['distance'] = list(c)
    elist = eventlist.loc[eventlist.distance <= radius]
    elist = elist[elist.depth > 0]
    return elist

def getNodesInEllipse(lat, lon, stk, radius, eventlist, alen):

    # Gather latitudes and longitudes for each point in eventlist
    elon = np.copy(eventlist['lon'].values)
    elat = np.copy(eventlist['lat'].values)

    aval = radius
    bval = radius/6

    rlon, rlat = heading(lon, lat, alen, stk)

    # Make array of length eventlist with azimuth from the node to each lat/lon point
    distances, az = npcosine(rlon, rlat, elon, elat)
    eventlist['distance'] = distances/111.19
    
    # Calculate maximum search distance for the associated strike and azimuth at this node
    mdist = []
    erta = math.sqrt(1-((math.pow(bval, 2))/(math.pow(aval, 2))))
    mdist = getEllipseRad(aval, erta, az, stk)
    
    # Filter points from eventlist that are outside of the ellipse
    eventlist['azimuth'] = az
    #print ('lon,lat,eventlist,mdist',lon,lat,eventlist,mdist)
    elist = eventlist.loc[eventlist.distance <= mdist]

    return elist

###############################################

### 53 ###

###############################################

## Written GLM 11.23.16
def findMultiDepth(lon, lat, nID, nonbinodes, spacing, multidepths, stk, slab, dep, alen, testprint):

    local = multidepths[multidepths.nID == nID]
    nonbinodes = nonbinodes[(nonbinodes.lon < lon+1)&(nonbinodes.lon > lon-1)]
    nonbinodes = nonbinodes[(nonbinodes.lat < lat+1)&(nonbinodes.lat > lat-1)]
    nearby = getNodesInEllipse(lat, lon, stk, spacing*2, nonbinodes, alen)
    nearby = nearby[nearby.distance > 0.02]
    #print (lon,lat,nonbinodes,local,nearby)
    

    if len(nearby)>0 and len(local)>0:
        loc_depth = np.mean(nearby['depth'].values)
        local['depdiff'] = np.abs(local['depth'].values - loc_depth)
        diffmin = local['depdiff'].min()
        peak_df = local[local.depdiff == diffmin]
        peak_depth = peak_df['depth'].values[0]
        if testprint:
            print ('ellipse yes lon,lat,loc_depth,peak_depth,nID',lon,lat,loc_depth,peak_depth,nID)
    
    elif len(local) > 0:
        local['nearS1'] = np.abs(local['depth'].values-dep)
        diffmin = local['nearS1'].min()
        peak_df = local[local.nearS1 == diffmin]
        peak_depth = peak_df['depth'].values[0]
        if testprint:
            print ('none in ellipse so Slab1  lon,lat,dep,peak_depth,nID',lon,lat,dep,peak_depth,nID)

    else:
        if testprint:
            print ('didnt go through multidepths: lon,lat,nID,local,nearby',lon,lat,nID,local,nearby)
        peak_depth = dep

    if testprint:
        fig = plt.figure()
        peaks = np.ones(len(local))*spacing/3
        plt.plot(nearby['distance'].values, nearby['depth'].values, 'bo')
        plt.plot(peaks, local['depth'].values, 'yo')
        plt.plot(spacing/2, peak_depth, 'ro')
        plt.xlabel('node distance')
        plt.ylabel('Depth')
        plt.grid()
        title = 'Lat: %.4f, Lon: %.4f, NID: %.4f' % (lat, lon, nID)
        plt.title(title)
        lontit = lon*100
        lattit = lat*100
        figtitle = 'Output/multitest_%s/pdf%i_2.png' % (slab, nID)
        #fig.savefig(figtitle)
        plt.close()

    return peak_depth

def findMultiDepthP(lon, lat, nID, nonbinodes, spacing, multidepths, stk, slab, dep, dip, alen, testprint):

    local = multidepths[multidepths.nID == nID]
    nonbinodes = nonbinodes[(nonbinodes.lon < lon+1)&(nonbinodes.lon > lon-1)]
    nonbinodes = nonbinodes[(nonbinodes.lat < lat+1)&(nonbinodes.lat > lat-1)]
    nearby = getNodesInEllipse(lat, lon, stk, spacing*2, nonbinodes, alen)
    nearby = nearby[nearby.distance > 0.02]
    #print (lon,lat,nonbinodes,local,nearby)
    
    if len(nearby)>0:
        loc_depth = np.mean(nearby['depth'].values)
        local['depdiff'] = np.abs(local['depth'].values - loc_depth)
        diffmin = local['depdiff'].min()
        peak_df = local[local.depdiff == diffmin]
        peak_depth = peak_df['depth'].values[0]
        peak_lon = peak_df['lon'].values[0]
        peak_lat = peak_df['lat'].values[0]
        if testprint:
            print ('ellipse yes lon,lat,loc_depth,peak_depth,nID,peak_lon,peak_lat',lon,lat,dep,peak_depth,nID,peak_lon,peak_lat)
    
    else:
        local['nearS1'] = np.abs(local['depth'].values-dep)
        diffmin = local['nearS1'].min()
        peak_df = local[local.nearS1 == diffmin]
        peak_depth = peak_df['depth'].values[0]
        peak_lon = peak_df['lon'].values[0]
        peak_lat = peak_df['lat'].values[0]
        if testprint:
            print ('none in ellipse so Slab1  lon,lat,dep,peak_depth,nID,peak_lon,peak_lat',lon,lat,dep,peak_depth,nID,peak_lon,peak_lat)

    if testprint:
        fig = plt.figure()
        peaks = np.ones(len(local))*spacing/3
        plt.plot(nearby['distance'].values, nearby['depth'].values, 'bo')
        plt.plot(peaks, local['depth'].values, 'yo')
        plt.plot(spacing/2, peak_depth, 'ro')
        plt.xlabel('node distance')
        plt.ylabel('Depth')
        plt.grid()
        title = 'Lat: %.4f, Lon: %.4f, NID: %.4f' % (lat, lon, nID)
        plt.title(title)
        lontit = lon*100
        lattit = lat*100
        figtitle = 'Output/multitest_%s/pdf%i_2.png' % (slab, nID)
        #fig.savefig(figtitle)
        plt.close()

    return peak_lon, peak_lat, peak_depth

###############################################

### 54 ###

###############################################

## Written GLM 12.01.16

def removePoints(donotuse, eventlist, lonmin, lonmax, latmin, latmax, printtest, datainfo, getfixed, slab):
    
    if len(donotuse) > 0:
        polyclip = makepolymask(slab,'library/misc/slab_polygons.txt')
        donotuse.loc[donotuse.lon < 0, 'lon']+=360
        eventlist.loc[eventlist.lon < 0, 'lon']+=360
        polyclip.loc[polyclip.lon < 0, 'lon']+=360
        
        pts = np.zeros((len(donotuse),2))
        pts[:, 0] = donotuse['lon'].values
        pts[:, 1] = donotuse['lat'].values
        
        mask = maskdatag(polyclip, pts)
        
        if slab == 'mue' or slab == 'car':
            donotuse['depth'] = donotuse['depth'].values*mask
            donotuse = donotuse[np.isfinite(donotuse.depth)]
            donotuse = donotuse.reset_index(drop=True)
            #donotuse.to_csv('%s_donotuse.csv'%slab,header=True,index=False)

    if getfixed:
        # Removing fixed events
        
        #fixedEQdepths = np.array([120.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 650.0])
        fixedEQdepths = np.array([10.0, 15.0, 20.0, 25.0, 33.0, 35.0, 47.0, 50.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 650.0])
        for i in fixedEQdepths:
            eventlist = eventlist[((eventlist.etype != 'EQ') & (eventlist.etype != 'ER')) | \
                    ((eventlist.etype == 'EQ') & (eventlist.depth != i)) | \
                    ((eventlist.etype == 'ER') & (eventlist.depth != i))]
    if len(donotuse) > 0:
        if printtest:
            print ('The following points were manually removed from the dataset: (list length = %i)'%len(donotuse))
        for index, r in donotuse.iterrows():
            lat, lon, depth, etype = r['lat'], r['lon'], r['depth'], r['etype']
            if slab == 'car' or slab == 'mue':
                near = eventlist[(eventlist.lon < lon+0.2) & (eventlist.lon > lon-0.2) & (eventlist.lat < lat+0.2) & (eventlist.lat > lat-0.2)]
            else:
                near = eventlist[(eventlist.lon < lon+0.5) & (eventlist.lon > lon-0.5) & (eventlist.lat < lat+0.5) & (eventlist.lat > lat-0.5)]
            for i, row in near.iterrows():
                
                latB, lonB, depthB, typeB = row['lat'], row['lon'], row['depth'], row['etype']
                d1 = abs(lat-latB)
                d2 = abs(lon-lonB)
                d3 = abs(depth-depthB)
                if d1 < 0.1 and d2 < 0.1 and d3 < 15 and etype == typeB:
                    if printtest:
                        data = row['ID']
                        addToDataInfo(data, 0, 'removepoints', datainfo, 'indiv')
                    if printtest:
                        print('lon,lat,depth,event-type', index, lonB, latB, depthB, typeB)
                    eventlist.drop(i, inplace=True)
            eventlist = eventlist.reset_index(drop=True)
    return eventlist

def doublePoints(doubleuse, eventlist, maxID):
    
    newIDs = list(range(maxID+1, maxID+1+len(doubleuse)))
    doubleuse['ID'] = newIDs
    eventlist = pd.concat([eventlist, doubleuse])
    
    maxID = eventlist['ID'].max()
    maxID += 1
    return eventlist, maxID


def mkSlabData(depgrid, strgrid, dipgrid, testprint):

    # get grid parameters
    gdict = depgrid.getGeoDict().copy()
    nx = gdict.nx
    ny = gdict.ny
    xmin = gdict.xmin
    xmax = gdict.xmax
    ymin = gdict.ymin
    ymax = gdict.ymax
    
    # print things if necessary
    if testprint:
        print('xmin,xmax,ymin,ymax', xmin, xmax, ymin, ymax)
    
    # make meshgrid in 0-360 degree longitudes
    if xmin < 0:
        xmin += 360
    if xmax < 0:
        xmax += 360
    xall = np.linspace(xmin, xmax, nx)
    yall = np.linspace(ymin, ymax, ny)
    n = len(xall)
    m = len(yall)
    if testprint:
        print('xmin,xmax,ymin,ymax', xmin, xmax, ymin, ymax)
    xpts, ypts = np.meshgrid(xall, yall)

    # move grids into flattened array
    slab1lons = xpts.flatten()
    slab1lats = ypts.flatten()
    slab1deps = np.flipud(depgrid.getData().copy()).flatten()
    slab1strs = np.flipud(strgrid.getData().copy()).flatten()
    slab1dips = np.flipud(dipgrid.getData().copy()).flatten()

    # eliminate grid coordinates with non-finite information
    slab1lons = slab1lons[np.isfinite(slab1dips)]
    slab1lats = slab1lats[np.isfinite(slab1dips)]
    slab1deps = slab1deps[np.isfinite(slab1dips)]
    slab1strs = slab1strs[np.isfinite(slab1dips)]
    slab1dips = slab1dips[np.isfinite(slab1dips)]
    slab1lons[slab1lons<0]+=360

    # store array in dataframe
    slab1data = pd.DataFrame({'lon':slab1lons,'lat':slab1lats,'depth':slab1deps,'strike':slab1strs,'dip':slab1dips})
    slab1data = slab1data[['lon', 'lat', 'depth', 'strike', 'dip']]

    return slab1data

def movingav(x):
    x2 = np.copy(x)
    for i in range(1, len(x)-1):
        thisaz = x[i]
        lastaz = x[i-1]
        nextaz = x[i+1]
        lastdiff = abs(thisaz-lastaz)
        nextdiff = abs(thisaz-nextaz)
        if thisaz < lastaz and thisaz < nextaz and (nextdiff>50 and lastdiff >50):
            x2[i] = (nextaz+lastaz)/2.0
            #print 'i,nextaz,lastaz,thisaz,x2[i]',i,nextaz,lastaz,thisaz,x2[i]
        elif thisaz > lastaz and thisaz > nextaz and (nextdiff>50 and lastdiff >50):
            x2[i] = (nextaz+lastaz)/2.0
        else:
            x2[i] = thisaz
    if abs(x[0]-x2[1]) > 90:
        x2[0] = x2[1]
    else:
        x2[0] = x[0]
    if abs(x[-1]-x2[-2]) > 90:
        x2[-1] = x2[-2]
    else:
        x2[-1] = x[-1]
    return x2

def myfilter(data):
    datalon = data.drop_duplicates(['lon'])
    datalat = data.drop_duplicates(['lat'])
    lons = datalon['lon'].values
    lats = datalat['lat'].values
    newdf = pd.DataFrame(columns = ['lon', 'lat', 'depth', 'unc'])
    for lon in lons:
        #print 'lon',lon
        thislon = data[data.lon == lon]
        thislon.sort(columns = 'lat')
        x = thislon['depth'].values
        lonsmooth = movingav(x)
        thislon['depth'] = lonsmooth
        frames = [newdf, thislon]
        newdf = pd.concat(frames)
        newdf = newdf.reset_index(drop=True)
        newdf = newdf[['lon', 'lat', 'depth', 'unc']]
    
        newdf2 = pd.DataFrame(columns = ['lon', 'lat', 'depth', 'unc'])
    for lat in lats:
        #print 'lat',lat
        thislat = newdf[newdf.lat == lat]
        thislat.sort(columns = 'lon')
        x = thislat['depth'].values
        latsmooth = movingav(x)
        thislat['depth'] = latsmooth
        frames = [newdf2, thislat]
        newdf2 = pd.concat(frames)
        newdf2 = newdf2.reset_index(drop=True)
        newdf2 = newdf2[['lon', 'lat', 'depth', 'unc']]

    return newdf2

def maskdatag(clip2, xi):

    clip = clip2.copy()
    clip.loc[clip.lon < 0, 'lon']+=360
    lons = clip['lon'].values
    lats = clip['lat'].values
    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(xi)
    mask1 = (np.zeros(len(temp),) * np.nan)
    mask1[temp] = 1
        
    return mask1


def makeErrorgrid(Surfgrid,xi,errordata):

    xpts = xi[:,0]
    ypts = xi[:,1]
    xpts.shape = Surfgrid.shape
    ypts.shape = Surfgrid.shape

    x = errordata[:,0]
    y = errordata[:,1]
    z = errordata[:,2]
    
    try:
        zi = griddata((x, y), z, (xpts, ypts), method='nearest')
    except:
        addx = np.random.rand(len(x))/1000
        x = x+addx
        y = y+addx
        z = z+addx
        zi = griddata((x, y), z, (xpts, ypts), method='nearest')

    zi.shape = Surfgrid.shape
    return zi

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def extendslightly(newdat,clip,data,dist,slab,shiftorfin,TRdata):

    pts = np.zeros((len(newdat),2))
    pts[:, 0] = newdat[:, 0]
    pts[:, 1] = newdat[:, 1]
    
    mask2 = maskdatag(clip, pts)
    maskdepths = np.multiply(newdat[:, 2], mask2)
    newdat[:, 2] = maskdepths
    newdat = newdat[~np.isnan(newdat).any(axis=1)]

    lons, lats = [], []
    for index,row in clip.iterrows():
        lon, lat = row['lon'], row['lat']
        if len(TRdata)>0 and (slab != 'sol' or lon > 150):
            loc_tr = TRdata[(TRdata.lon > lon-3) & (TRdata.lon < lon+3) & (TRdata.lat > lat-3) & (TRdata.lat < lat+3)]
            if len(loc_tr)>0:
                #loc_tr['dist'] = gps2dist_azimuth(lat, lon, loc_tr['lat'], loc_tr['lon'])[0]/1000.0
                loc_tr['dist'], tempangles = npcosine(lon, lat, loc_tr['lon'].values, loc_tr['lat'].values)
                mindist = loc_tr['dist'].min()
                loc_tr = loc_tr[loc_tr.dist == mindist]
                lonT = loc_tr['lon'].values[0]
                latT = loc_tr['lat'].values[0]
                azT = loc_tr['az'].values[0]
                thisdist, thisang, latB, lonB = cosine(lonT, latT, lon, lat)
                out = isoutboard(azT, thisang)
                if out:
                    continue
                else:
                    lons.append(lon)
                    lats.append(lat)
            else:
                lons.append(lon)
                lats.append(lat)
        else:
            lons.append(lon)
            lats.append(lat)
                    
    cliparr = np.zeros((len(lons),6))
    cliparr[:,0] = lons
    cliparr[:,1] = lats
    
    strclip = []
    for i in range(len(cliparr)-1):
        d, a, la1, lo1 = cosine(cliparr[i,0], cliparr[i,1], cliparr[i+1,0], cliparr[i+1,1])
        strclip.append(a)
    strclip.append(strclip[0])

    cliparr[:,2] = strclip

    if slab == 'sum':
        cliparr = cliparr[(cliparr[:,1] < 25)]
    if slab == 'kur' or slab == 'izu' or slab == 'jap':
        cliparr = cliparr[(cliparr[:,1] < 35)|(cliparr[:,1] > 41)]
        cliparr = cliparr[(cliparr[:,1] < 45)|(cliparr[:,1] > 50)]
    if slab == 'kerz':
        cliparr = cliparr[(cliparr[:,1] > -30)|((cliparr[:,1] > -38.5)&(cliparr[:,1] < -37.5))]
    if slab == 'izu':
        cliparr = cliparr[((cliparr[:,1] < 15)|(cliparr[:,1] > 27))|((cliparr[:,1] > 21.5)&(cliparr[:,1] < 23))]
    if slab == 'manz':
        cliparr = cliparr[(cliparr[:,1] > 1)&(cliparr[:,1] < 16)]
    if slab == 'sol':
        cliparr = cliparr[(cliparr[:,0] > 155)|(cliparr[:,0] < 152)]
    '''
    if slab == 'sam':
        cliparr = cliparr[(cliparr[:,1] < -42)|(cliparr[:,1] > -39)]
    if slab == 'sum':
        cliparr = cliparr[(cliparr[:,1] < 25)]
        cliparr = cliparr[(cliparr[:,0] < 120)]
    if slab == 'hin':
        cliparr = cliparr[(cliparr[:,1] > 36)]
    if slab == 'sol':
        cliparr = cliparr[(cliparr[:,0] < 160)]
        cliparr = cliparr[(cliparr[:,0] > 148)]
    if slab == 'alu':
        cliparr = cliparr[(cliparr[:,0] < 212)]
    '''
        
    clipstr = griddata(newdat[:, 0:2], newdat[:, 3], cliparr[:, 0:2], method='nearest')
    cliparr[:,3] = clipstr
    azstrdiff = abs(cliparr[:,2] - cliparr[:,3])
    cliparr[:,4] = azstrdiff

    #np.savetxt('cliparr0.csv', cliparr, header='lon,lat,az,strike,diff,other',fmt='%.2f', delimiter=',',comments='')
    cliparr = cliparr[((cliparr[:,4] > 160)&(cliparr[:,4] < 200))|((cliparr[:,4] < 20)|(cliparr[:,4] > 340))]

    clipstr = cliparr[:,3]
    
    #np.savetxt('cliparr1.csv', cliparr, header='lon,lat,az,strike,diff,other',fmt='%.2f', delimiter=',',comments='')
    clip180 = clipstr-90
    clip180[clip180<0]+=360
    lmax = 2.5
    wmax = 0.2
    clipdip = []
    clipstr = []
    clipdep = []
    clipdist = []
    ogdips = []
    dipadds = []
    idnot = []
    for i in range(len(cliparr)):
        plon,plat,pstr = cliparr[i,0], cliparr[i,1], clip180[i]
        testdat = projLW(lmax, wmax, plon, plat, pstr, newdat, ['lon','lat','depth','strike','dip'])
        
        if len(testdat) < 1:
            idnot.append(i)
            continue
        
        if testdat['strike'].min() < 45 and testdat['strike'].max() > 315:
            testdat.loc[testdat.strike < 90, 'strike'] += 360
        maxdip = testdat['dip'].max()
        meanstr = testdat['strike'].mean()
        maxdep = testdat['depth'].max()
        
        depthdist = testdat[testdat.depth == maxdep]
        distdepth = depthdist['dist'].values[0]
        diptest = testdat[testdat.dip == testdat['dip'].max()]
        disttest = diptest['dist'].values[0]
        
        gradval = 25
        diptest2 = testdat[(testdat.dist < disttest+gradval)&(testdat.dist > disttest)]
        if len(diptest2) < 1:
            idnot.append(i)
            continue
        
        maxdistance = diptest2['dist'].max()
        maxdist2 = diptest2[diptest2.dist == maxdistance]
        fardip = maxdist2['dip'].values[0]
        dipadd = maxdip - fardip
        gradfrac = dist/gradval
        dipadd *= gradfrac
        
        if disttest < 100:
            maxdip += dipadd
        distdepth = disttest

        cliparr[i,0], cliparr[i,1] = heading(cliparr[i,0],cliparr[i,1],distdepth,clip180[i])
        if maxdip < 90:
            clipdip.append(maxdip)
            clipstr.append(meanstr)
            clipdep.append(maxdep)
            clipdist.append(distdepth)
            ogdips.append(maxdip)
            dipadds.append(dipadd)
        else:
            idnot.append(i)

    cliparr = np.delete(cliparr,idnot,axis=0)
    cliptesting = pd.DataFrame({'lon':cliparr[:,0],'lat':cliparr[:,1],'depth':clipdep,'strike':clipstr,'dip':clipdip,'dist':clipdist,'dipadd':dipadds,'ogdip':ogdips})
    #cliptesting.to_csv('%s_projtesting.csv'%slab,header=True,index=False)
    clipdip = np.array(clipdip)
    clipstr = np.array(clipstr)
    clipdep = np.array(clipdep)
    
    #np.savetxt('cliparr2.csv', cliparr, header='lon,lat,az,strike,diff,other',fmt='%.2f', delimiter=',',comments='')

    if slab == 'phi' or slab == 'sul' or slab == 'cot' or slab == 'mue' or slab == 'cot':
        depthcutoff = 30
    else:
        depthcutoff = 90

    cliparr = cliparr[clipdep>depthcutoff]
    clipstr = clipstr[clipdep>depthcutoff]
    clipdip = clipdip[clipdep>depthcutoff]
    clipdep = clipdep[clipdep>depthcutoff]
    
    clipazs = clipstr+90.0
    clipazs[clipazs>360]-=360
    dists = np.ones(len(cliparr))*dist

    extlons = []
    extlats = []
    for i in range(len(cliparr)):
        #print (cliparr[i,0],cliparr[i,1],dists[i],clipstr[i],clipdip[i])
        extlon, extlat = heading(cliparr[i,0],cliparr[i,1],dists[i],clipazs[i])
        extlons.append(extlon)
        extlats.append(extlat)

    extdeps = clipdep + dists/np.tan(np.radians(90-clipdip))

    addarr = np.zeros((len(extlons),5))
    addarr[:,0] = extlons
    addarr[:,1] = extlats
    addarr[:,2] = extdeps
    addarr[:,3] = clipstr
    addarr[:,4] = clipdip

    addarr = addarr[addarr[:,2] < (2*np.max(data[:,2]))]
    notlist = []
    testr = 1.0
    if slab == 'van' or slab == 'sco' or slab == 'puy' or slab == 'man':
        testr = 0.5
    
    for i in range (len(addarr)):
        x,y,z,s = addarr[i,0], addarr[i,1], addarr[i,2], addarr[i,3]
        if s > 360:
            s -= 360
        nearnodes = data[(data[:,0] > x-testr)&(data[:,0] < x+testr) & \
                          (data[:,1] > y-testr)&(data[:,1] < y+testr)]
        if len(nearnodes) < 1:
            continue
        else:
            dists, angs = npcosine(x, y, nearnodes[:,0], nearnodes[:,1])
            angs -= 180
            angs[angs<0]+=360
            noutb = np.logical_not(npoutboard(s, angs))
            if len(noutb[noutb == False])>0:
                notlist.append(i)
            else:
                continue

    addarr = np.delete(addarr, notlist, axis = 0)
    
    ''' remove sam bad nodes here!!!!!!!!'''

    return addarr

def rbffill(data, sigma, lonname, latname, depthname, filter, slab, smoother, gridsp):
    
    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    x[x < 0] += 360

    xi = np.arange(np.floor(np.min(x))-3, np.ceil(np.max(x))+3, gridsp)
    yi = np.arange(np.floor(np.min(y))-3, np.ceil(np.max(y))+3, gridsp)
    xpts, ypts = np.meshgrid(xi, yi)

    interp = Rbf(x, y, z, function='linear', smooth=smoother)
    zi = interp(xpts, ypts)

    xyzip = np.zeros((len(xpts.flatten()),2))
    xyzip[:, 0] = xpts.flatten()
    xyzip[:, 1] = ypts.flatten()
    zif = ndimage.filters.gaussian_filter(zi, sigma/2)
    
    strikegrid, dipgrid = mkSDgrddata(xyzip, zif, 'flip')
    
    newdat = np.zeros((len(zif.flatten()), 5))
    newdat[:, 0], newdat[:, 1], newdat[:, 2] = xpts.flatten(), ypts.flatten(), zif.flatten()
    newdat[:, 3], newdat[:, 4] = strikegrid.flatten(), dipgrid.flatten()

    pts = np.zeros((len(newdat),2))
    pts[:, 0] = newdat[:, 0]
    pts[:, 1] = newdat[:, 1]

    return newdat

def linfill(data, sigma, lonname, latname, depthname, filter, slab, node):
    
    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    x[x < 0] += 360

    #np.savetxt('griddatatest1.csv', data, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
    # define grid.
    gridsp = node
    xi = np.arange(np.floor(np.min(x))-3, np.ceil(np.max(x))+3, gridsp)
    yi = np.arange(np.floor(np.min(y))-3, np.ceil(np.max(y))+3, gridsp)

    xpts, ypts = np.meshgrid(xi, yi)
    # grid the data.
    try:
        zi = griddata((x, y), z, (xpts, ypts), method='linear')
    except:
        addx = np.random.rand(len(x))/1000
        x = x+addx
        y = y+addx
        z = z+addx
        zi = griddata((x, y), z, (xpts, ypts), method='linear')
    
    xyzip = np.zeros((len(xpts.flatten()),2))
    xyzip[:, 0] = xpts.flatten()
    xyzip[:, 1] = ypts.flatten()
    zif = ndimage.filters.gaussian_filter(zi, sigma/2)
    
    strikegrid, dipgrid = mkSDgrddata(xyzip, zif, 'flip')
    
    newdat = np.zeros((len(zif.flatten()), 5))
    newdat[:, 0], newdat[:, 1], newdat[:, 2] = xpts.flatten(), ypts.flatten(), zif.flatten()
    newdat[:, 3], newdat[:, 4] = strikegrid.flatten(), dipgrid.flatten()

    pts = np.zeros((len(newdat),2))
    pts[:, 0] = newdat[:, 0]
    pts[:, 1] = newdat[:, 1]

    newdat = newdat[~np.isnan(newdat).any(axis=1)]
    
    return newdat

def pySurface3(data, node, T, slab, spacing, deptherr, time, these_parameters, filter, filldat, nCores, TRdata, meanBA, kdeg, knot_no, rbfs, shift_out, shiftorfin, extra):

    dataxmin = np.min(data[:,0])
    dataxmax = np.max(data[:,0])
    dataymin = np.min(data[:,1])
    dataymax = np.max(data[:,1])
    
    if len(filldat)<1:
        if slab == 'camz':
            sigma = 0.1
        else:
            sigma = 1
        rfbdata = rbffill(data, 0.01, 'lon', 'lat', 'depth', filter, slab, 10, spacing)
        if node < 0.05:
            filldat = linfill(rfbdata, 0.01, 'lon', 'lat', 'depth', 100, slab, 0.05)
        else:
            filldat = linfill(rfbdata, 0.01, 'lon', 'lat', 'depth', 100, slab, node)
        #may need to not do this for split surface
        newres = mkContourClip(shift_out, TRdata, spacing, filldat, False, slab)
        if len(TRdata)>0:
            clip2 = clippingmask(newres,TRdata,node,False, slab, 'first')
        else:
            clip2 = noTrenchPolygon(newres, node, False, slab)
    
        if extra != 'manz':
            dataadd1 = extendslightly(filldat,clip2,data,50,slab,shiftorfin,TRdata)
            dataadd2 = extendslightly(filldat,clip2,data,70,slab,shiftorfin,TRdata)
            dataadd3 = extendslightly(filldat,clip2,data,90,slab,shiftorfin,TRdata)
            dataadd4 = extendslightly(filldat,clip2,data,110,slab,shiftorfin,TRdata)

            dataadd1 = dataadd1[dataadd1[:,2] > 30]
            dataadd2 = dataadd2[dataadd2[:,2] > 30]
            dataadd3 = dataadd3[dataadd3[:,2] > 30]
            dataadd4 = dataadd4[dataadd4[:,2] > 30]
            
            extdata1 = np.vstack((dataadd1[:,:5],dataadd2[:,:5],dataadd3[:,:5],dataadd4[:,:5]))
            extdata = np.vstack((data[:,:3],dataadd1[:,:3],dataadd2[:,:3],dataadd3[:,:3],dataadd4[:,:3]))
            rfbdata = rbffill(extdata, 0.0001, 'lon', 'lat', 'depth', filter, slab, rbfs, spacing)
            #np.savetxt('%s_griddatatest21.csv'%slab, extdata1, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
            #np.savetxt('%s_griddatatest2.csv'%slab, extdata, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
        else:
            rfbdata = rbffill(data[:,:3], 0.0001, 'lon', 'lat', 'depth', filter, slab, rbfs, spacing)
        if node < 0.05:
            filldat = linfill(rfbdata, 0.0001, 'lon', 'lat', 'depth', 100, slab, 0.05)
        else:
            filldat = linfill(rfbdata, 0.0001, 'lon', 'lat', 'depth', 100, slab, node)
        
        #np.savetxt('%s_griddatatest3.csv'%slab, rfbdata, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
        #np.savetxt('%s_griddatatest4.csv'%slab, filldat, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
        
        filldat = filldat[~np.isnan(filldat).any(axis=1)]
        if slab == 'kur':
            filldat = filldat[(filldat[:,0]<dataxmax+1) & (filldat[:,0]>dataxmin-4) \
                            & (filldat[:,1]<dataymax+1) & (filldat[:,1]>dataymin-1)]
        else:
            filldat = filldat[(filldat[:,0]<dataxmax+1) & (filldat[:,0]>dataxmin-1) \
                            & (filldat[:,1]<dataymax+1) & (filldat[:,1]>dataymin-1)]
        filldat[:,3] = 100
        filldat = filldat[:, 0:4]

    data[:, 0][data[:, 0]<0] += 360
    filldat[:, 0][filldat[:, 0]<0] += 360
    xmin, xmax = np.min(data[:, 0]), np.max(data[:, 0])
    ymin, ymax = np.min(data[:, 1]), np.max(data[:, 1])
    
    deglats = (data[:, 1] - 90)*-1.0
    radlats = np.radians(deglats)
    radlons = np.radians(data[:, 0])

    rxn,rxx,ryn,ryx = np.min(radlons),np.max(radlons),np.min(radlats),np.max(radlats)

    rnode = np.radians(node)
    rbuff = np.radians(3.0)
    xall = np.arange(rxn-rbuff, rxx+rbuff, rnode)
    if slab == 'kur':
        xall = np.arange(rxn-(rbuff*2.5), rxx+rbuff, rnode)
    yall = np.arange(ryn-rbuff, ryx+rbuff, rnode)
    dl = False

    n = len(xall)
    m = len(yall)

    xpts, ypts = np.meshgrid(xall, yall)
    xi = np.zeros((m*n, 2))

    xi[:, 0] = np.degrees(xpts.flatten())
    xi[:, 1] = 90.0 - np.degrees(ypts.flatten())

    data = np.vstack((data[:,:4], filldat[:,:4]))
    data[:, 3][np.isnan(data[:, 3])] = 40

    x = data[:, 0]
    y = (data[:, 1]-90)*-1
    z = data[:, 2]
    w = 1.0/data[:, 3]

    xrad = np.radians(x)
    yrad = np.radians(y)
    yrad[yrad<0] = math.pi/2.0+np.abs(yrad[yrad<0])
    zrad = z

    ntx = int(abs(np.floor(xmin)-np.ceil(xmax))*knot_no)
    nty = int(abs(np.floor(ymin)-np.ceil(ymax))*knot_no)

    tx = np.linspace(xall.min(), xall.max(), ntx)
    ty = np.linspace(yall.min(), yall.max(), nty)
    
    if deptherr == 'depth':
        f = open(these_parameters, 'a')
        f.write('knot_no: %s \n' % str(knot_no))
        f.write('kdeg: %s \n' % str(kdeg))
        f.close()
    print ('               interpolating ....')
    lut = LSQSphereBivariateSpline(yrad, xrad, zrad, ty[1:-1], tx[1:-1], w=w)
    print ('               interpolated')
    interpdepths = lut.ev(ypts.flatten(), xpts.flatten())
    interpdepths.shape = xpts.shape

    return interpdepths, xi, dl

def pySurface4(data, node, T, slab, spacing, deptherr, time, these_parameters, filter, filldat, nCores, TRdata, meanBA, kdeg, knot_no, rbfs, shift_out, shiftorfin, extra):

    knot_no = 2
    data[:, 0][data[:, 0]<0] += 360
    if len(filldat) > 0:
        filldat[:, 0][filldat[:, 0]<0] += 360
    xmin, xmax = np.min(data[:, 0]), np.max(data[:, 0])
    ymin, ymax = np.min(data[:, 1]), np.max(data[:, 1])
    
    deglats = (data[:, 1] - 90)*-1.0
    radlats = np.radians(deglats)
    radlons = np.radians(data[:, 0])

    rxn,rxx,ryn,ryx = np.min(radlons),np.max(radlons),np.min(radlats),np.max(radlats)

    rnode = np.radians(node)
    rbuff = np.radians(3.0)
    xall = np.arange(rxn-rbuff, rxx+rbuff, rnode)
    if slab == 'kur':
        xall = np.arange(rxn-(rbuff*2.5), rxx+rbuff, rnode)
    yall = np.arange(ryn-rbuff, ryx+rbuff, rnode)
    dl = False

    n = len(xall)
    m = len(yall)

    xpts, ypts = np.meshgrid(xall, yall)
    xi = np.zeros((m*n, 2))

    xi[:, 0] = np.degrees(xpts.flatten())
    xi[:, 1] = 90.0 - np.degrees(ypts.flatten())

    if len(filldat) > 0:
        data = np.vstack((data[:,:4], filldat[:,:4]))
    data[:, 3][np.isnan(data[:, 3])] = 40

    x = data[:, 0]
    y = (data[:, 1]-90)*-1
    z = data[:, 2]
    w = 1.0/data[:, 3]

    xrad = np.radians(x)
    yrad = np.radians(y)
    yrad[yrad<0] = math.pi/2.0+np.abs(yrad[yrad<0])
    zrad = z

    ntx = int(abs(np.floor(xmin)-np.ceil(xmax))*knot_no)
    nty = int(abs(np.floor(ymin)-np.ceil(ymax))*knot_no)

    tx = np.linspace(xall.min(), xall.max(), ntx)
    ty = np.linspace(yall.min(), yall.max(), nty)
    
    if deptherr == 'depth':
        f = open(these_parameters, 'a')
        f.write('knot_no: %s \n' % str(knot_no))
        f.write('kdeg: %s \n' % str(kdeg))
        f.close()
    print ('               interpolating ....')
    lut = LSQSphereBivariateSpline(yrad, xrad, zrad, ty[1:-1], tx[1:-1], w=w)
    print ('               interpolated')
    interpdepths = lut.ev(ypts.flatten(), xpts.flatten())
    interpdepths.shape = xpts.shape

    return interpdepths, xi, dl


def gblockmean(xpts,ypts,interpdepths,tx,ty,bwidth,w,xrad,yrad):

    pbarr = np.zeros((len(xpts.flatten()),4))
    pbarr[:,0] = xpts.flatten()
    pbarr[:,1] = ypts.flatten()
    pbarr[:,2] = interpdepths.flatten()
    xyrad = np.zeros((len(xrad),2))
    xyrad[:,0] = xrad
    xyrad[:,1] = yrad
    pbarr[:,3] = griddata(xyrad, w, pbarr[:, 0:2], method='nearest')

    txx, tyy = np.meshgrid(tx,ty)
    blockx = txx.flatten()
    blocky = tyy.flatten()

    rwidth = np.radians(bwidth/2)
    blockd, blockw = [],[]
    for i in range(len(blockx)):
        rlon, rlat = blockx[i], blocky[i]
        thisarr = pbarr[(pbarr[:,0] > rlon-rwidth) & (pbarr[:,0] < rlon+rwidth) & \
                        (pbarr[:,1] > rlat-rwidth) & (pbarr[:,1] < rlat+rwidth)]
        newdep = np.average(thisarr[:,2], weights = thisarr[:,3])
        newwht = np.average(thisarr[:,3], weights = thisarr[:,3])
        blockd.append(newdep)
        blockw.append(newwht)

    blockdata = np.zeros((len(blockx),4))
    blockdata[:,0] = blockx
    blockdata[:,1] = blocky
    blockdata[:,2] = blockd
    blockdata[:,3] = blockw

    return blockdata

def perpPDFdepths(elist, cstr, cdip, lon, lat, loc_depth, maxthickness):

    hd2 = math.cos(math.radians(cdip))
    zd2 = math.sin(math.radians(cdip))
    cdip -= 90
    radstr = cstr * math.pi/180.0
    raddip = cdip * math.pi/180.0
    xs = math.cos(radstr)
    ys = math.sin(radstr)
    hd = math.cos(raddip)
    zd = math.sin(raddip)

    zdist = elist['depth'].values - loc_depth # (- -1*loc_depth)

    elist['zdist'] = abs(zdist)
    
    elist['cosdistance'], cosangles = npcosine(lon, lat, elist['lon'].values, elist['lat'].values)
    cosangles -= 180
    cosangles[cosangles<0]+=360
    elist['outboard'] = np.logical_not(npoutboard(cstr, cosangles)) # Will need to fix when fix outboard function
    
    cosangles[cosangles <= 180.0] += 360.0
    cosangles -= 180.0
    elist['anglediff'] = abs(cstr - cosangles)
    elist['phiS'] = abs(elist['anglediff']-90)
    
    elist['cosdistance'] = (elist['cosdistance'].values * np.cos(np.radians(elist['phiS'].values)))

    elist['cosdistance'][(elist.outboard == True) & (elist.cosdistance > 0)] *= -1
    elist['cosdistance'][(elist.outboard == False) & (elist.cosdistance < 0)] *= -1
    elist.loc[elist.etype == 'CP', 'cosdistance'] = 0.0
    elist.loc[elist.etype == 'RF', 'cosdistance'] = 0.0
    elist.loc[elist.etype == 'AS', 'cosdistance'] = 0.0
    elist.loc[elist.etype == 'BA', 'cosdistance'] = 0.0
    elist['alldist'] = np.sqrt(elist['zdist'].values * elist['zdist'].values + elist['cosdistance'].values * elist['cosdistance'].values)
    dangle = getangle(hd, 0.0, zd, elist['cosdistance'].values, np.zeros(len(zdist)), zdist)
    elist['dangle'] = (dangle * 180/math.pi)
    phiD = abs(elist['dangle'] - 90.0)
    
    elist['perpdistance'] = (elist['alldist'].values*np.cos(dangle))
    elist['perpdistance'][(elist.outboard == True) & (elist.cosdistance < 0)] *= -1
    elist['perpdistance'][(elist.outboard == False) & (elist.cosdistance > 0)] *= -1
    elist['perpdistance'][elist.etype == 'RF'] *= -1
    elist['perpdistance'][elist.etype == 'AS'] *= -1
    elist['perpdistance'][elist.etype == 'BA'] *= -1
    elist['perpdistance'][elist.etype == 'CP'] *= -1
    
    maxperp, minperp = elist['perpdistance'].max(), elist['perpdistance'].min()
    diffmax = maxperp-minperp
    meandist = elist['perpdistance'].mean()
    removelist = pd.DataFrame()
    
    while diffmax > maxthickness and len(elist) > 1:
        #print 'too wide!! lon,lat,loc_depth,diffmax,maxperp,minperp,meandist,maxthickness',lon,lat,loc_depth,diffmax,maxperp,minperp,meandist,maxthickness
        if abs(maxperp-meandist) > abs(minperp-meandist):
            removelist = pd.concat([removelist, elist[elist.perpdistance == maxperp]])
            elist = elist[elist.perpdistance != maxperp]
            maxperp = elist['perpdistance'].max()
        else:
            removelist = pd.concat([removelist, elist[elist.perpdistance == minperp]])
            elist = elist[elist.perpdistance != minperp]
            minperp = elist['perpdistance'].min()

        meandist = elist['perpdistance'].mean()
        diffmax = maxperp - minperp
    
    elist = elist[['ID', 'perpdistance', 'outboard', 'cosdistance']]

    #if len(removelist) > 0:
    #    print ('removelist!!',lon,lat,removelist)
    return elist, hd2, zd2, removelist


def perpPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, whichpdf, slab, strike, dip, maxthickness):

    multipeaks = pd.DataFrame()
    
    cstr, cdip = strike, dip

    elistASe = trimmed[trimmed.etype == 'AS' ]
    elistBAe = trimmed[trimmed.etype == 'BA' ]
    elistAAe = trimmed[trimmed.etype == 'AA' ]
    elistRFe = trimmed[trimmed.etype == 'RF' ]
    elistTOe = trimmed[trimmed.etype == 'TO' ]
    elistCPe = trimmed[trimmed.etype == 'CP' ]
    
    if len(elistAAe)>0 and len(trimmed) <4:
        if abs(elistAAe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from AA',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]

    if len(elistASe)>0 and len(trimmed) <5:
        if abs(elistASe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from AS',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]

    if len(elistBAe)>0 and len(trimmed) <5:
        if abs(elistBAe['depth'].mean() - trimmed['depth'].mean()) > 50:
            #print 'EQ too different from BA',lon,lat,trimmed
            trimmed = trimmed[(trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'BA')]
    
    nantest = trimmed['depth'].values
    nantest = nantest[np.isnan(nantest)]
    
    if len(nantest) > 0:
        #print 'NAN problem?? lon,lat,nID,sdepth,ddepth,loc_depth,trimmed',lon,lat,nID,sdepth,ddepth,loc_depth,trimmed
        peak_depth = np.nan
        stdv = np.nan
        test = False
        n = 0
        return lon, lat, peak_depth, stdv, test, 0, multipeaks, stdv, pd.DataFrame()
    
    multi = False
    n = 0
    
    if len(trimmed)>1:

        # Distinguishing between different data types
        ASframe = trimmed[trimmed.etype == 'AS']
        AAframe = trimmed[trimmed.etype == 'AA']
        EQframe = trimmed[trimmed.etype == 'EQ']
        BAframe = trimmed[trimmed.etype == 'BA']
        ERframe = trimmed[trimmed.etype == 'ER']
        TOframe = trimmed[trimmed.etype == 'TO']
        RFframe = trimmed[trimmed.etype == 'RF']
        CPframe = trimmed[trimmed.etype == 'CP']

        # Adding present event types to list of event types
        #and calculate average rake, strike, and dip for output file if CMT info available
        etypes = []
        AA = False
        AS = False
        BA = False
        RF = False
        TO = False
        ER = False
        EQ = False
        CP = False
        if len(ASframe) > 0:
            etypes.append('AS')
            AS = True
        if len(AAframe) > 0:
            etypes.append('AA')
            AA = True
        if len(EQframe) > 0 or len(ERframe)>0:
            etypes.append('EQ')
            if len(EQframe) > 0:
                EQ = True
            if len(ERframe) > 0:
                ER = True
        if len(BAframe) > 0:
            etypes.append('BA')
            BA = True
        #if len(ERframe > 0):
        #    etypes.append('ER')
        #    ER = True
        if len(TOframe) > 0:
            etypes.append('TO')
            TO = True
        if len(RFframe) > 0:
            etypes.append('RF')
            RF = True
        if len(CPframe) > 0:
            etypes.append('CP')
            CP = True
            
        # Make perpendicular PDF
        dip90 = dip-90
        elist, hd, zd, removelist = perpPDFdepths(trimmed, strike, dip, lon, lat, loc_depth, maxthickness)
        
        if len(elist) < 2:
            #print 'data too dispersed and sparse to resolve a depth',lon,lat,nID
            test = False
            stdv = np.nan # GLM 11.14.16 investigate this exception if missing PDFs
            peak_depth = np.nan
            return lon, lat, peak_depth, stdv, test, 0, multipeaks, stdv, removelist
        
        trimmed.merge(elist, on='ID')
        trimmed = trimmed[np.isfinite(trimmed['perpdistance'].values)]
        sPdepth = sdepth/math.sin(math.radians(dip90))
        dPdepth = ddepth/math.sin(math.radians(dip90))*-1
        if testprint:
            print ('spdepth,dpdepth',lon,lat,loc_depth,sdepth,sPdepth,ddepth,dPdepth,dip90)
        try:
            if dip<45:
                dep_range = np.arange(sPdepth-15, dPdepth+15, 1)
            else:
                dep_range = np.arange(-250, 250, 1)
        except:
            dep_range = np.arange(-250, 250, 1)

        PDF = makePDF4(trimmed, dep_range, etypes, testprint, 'perpdistance')
        PDF_df1 = pd.DataFrame(dep_range, columns=['depths'])
        PDF_df1['Summed_Values'] = PDF

        # Eliminates values less than 0.001 and finds min, max, and peak depth in PDF
        
        if len(PDF_df1)>0:
            PDF_df = PDF_df1.loc[PDF_df1.Summed_Values >= 0.001]

            if len(PDF_df) < 1:
                PDF_df = PDF_df1.loc[PDF_df1.Summed_Values >= 0.0001]
                if len(PDF_df) < 1:
                    #print ('noperp PDF? lon,lat,nID,sdepth,ddepth,loc_depth,trimmed',lon,lat,nID,sdepth,ddepth,loc_depth,trimmed)
                    peak_depth = np.nan
                    peak_lat = np.nan
                    peak_lon = np.nan
                    multipeaks = []
                    centsurf = 1
                    removelist = []
                    stdv = np.nan
                    test = False
                    n = 0
                    return peak_lon, peak_lat, peak_depth, stdv, test, n, multipeaks, centsurf, removelist
                
        else:
            #print 'noPDFperp',lon,lat,nID
            test = False
            stdv = np.nan # GLM 11.14.16 investigate this exception if missing PDFs
            peak_depth = np.nan
            return lon, lat, peak_depth, stdv, test, 0, multipeaks, stdv, removelist
        
        if len(PDF_df) < 1:
            #print 'PDF too poorly dispersed to resolve depth',lon,lat,nID
            test = False
            stdv = np.nan # GLM 11.14.16 investigate this exception if missing PDFs
            peak_depth = np.nan
            return lon, lat, peak_depth, stdv, test, 0, multipeaks, stdv, removelist
        
        if AA or AS or BA or RF or TO or CP or (ER and EQ and slab != 'kur'):
        #if AA or AS or BA or RF or TO or ER:
            peak = PDF_df['Summed_Values'].max()
            peakbuffer = 0.1*peak
            depthbuffer = 10
            d_min = PDF_df['depths'].min()
            d_max = PDF_df['depths'].max()

            # Finding the depth associated with the peak PDF value
            peak_df = PDF_df[PDF_df.Summed_Values == peak]
            peak_depth = peak_df['depths'].values[0]
            meandepth = False

        else:
            meandepth = PDF_df['depths'].mean()
            PDF_df['meandiff'] = np.abs(PDF_df['depths'].values - meandepth)
            meanmin = PDF_df['meandiff'].min()
            peak_df = PDF_df[PDF_df.meandiff == meanmin]
            peak_depth = peak_df['depths'].values[0]
            peak = peak_df['Summed_Values'].values[0]
            peakbuffer = 0.01
            meandepth = True

        # GLM 11.22.16 - adding bimodal distribution condition
        PDF_df['buffer'] = PDF_df['Summed_Values'].values + peakbuffer
        multi_peak = PDF_df[PDF_df.buffer >= peak]
        #multi_peak = PDF_df[(PDF_df.buffer >= peak) & ((PDF_df.depths < peakmin) | (PDF_df.depths > peakmax))] # GLM 11.25.16
        multi_peak2 = getLocalMax(multi_peak)
        
        if len(multi_peak2)>1 and not meandepth:
            multipeaks = makeMultiDFP(multi_peak2, multipeaks, lon, lat, nID, strike, dip, loc_depth)
            multi = True
            test = True
            n = len(multi_peak2)
        else:
            try:
                test = True
                multi = False
                n = 1
            except:
                #print 'multidepth PDF Exception: lon,lat,nID: ',lon,lat,nID
                test = False
                stdv = np.nan # GLM 11.14.16 investigate this exception if missing PDFs
                peak_depth = np.nan
                return peak_depth, stdv, test, 0, multipeaks, stdv, removelist

        # Finding standard deviation of PDF
        thissum = 0
        for d in PDF_df['depths'].values:
            residual = peak_depth - d
            thissum += residual * residual
        stdv = math.sqrt(1.0/len(PDF_df)*thissum)

        minonperp = PDF_df['depths'].min()
        maxonperp = PDF_df['depths'].max()
        centsurf = abs(peak_depth-minonperp)

        # Plotting PDF along axis perpendicular to dip
        probs = PDF_df['Summed_Values'].values*100
        ndist = PDF_df['depths'].values
        totdist = np.sqrt(probs*probs + ndist*ndist)
        alpha = np.arccos(ndist/totdist)
        phi = math.radians(dip90)+alpha
        xplot = totdist*np.cos(phi)*-1
        yplot = totdist*np.sin(phi)*-1

        peakprob = peak*100
        peakndist = peak_depth
        peaktotdist = math.sqrt(peakprob*peakprob + peakndist*peakndist)
        peakalpha = math.acos(peakndist/peaktotdist)
        peakphi = math.radians(dip90)+peakalpha
        peakxplot = peaktotdist*math.cos(peakphi)*-1
        peakyplot = peaktotdist*math.sin(peakphi)*-1

        depphi = 90-abs(dip90)
        peak_depthz = peak_depth*math.cos(math.radians(depphi))
        peak_depthx = peak_depth*math.sin(math.radians(depphi))
        
        if peak_depthz > 0:
            shiftaz = az_perp(strike)
        else:
            shiftaz = az_other_perp(strike)
        if abs(peak_depthx) > 0.001:
            peaklon, peaklat = heading(lon, lat, abs(peak_depthx), shiftaz)
        else:
            peaklon, peaklat = lon, lat


        perpdepths = PDF_df['depths'].values*math.cos(math.radians(depphi))
        PDF_df['depths'] = perpdepths+loc_depth
        peak_depth = loc_depth + peak_depthz

        # For testing PDFs of specific points - change lat-lon ranges to use
        if testprint:
        
            plotsx = math.cos(math.radians(dip90))*(minonperp-10)
            plotsz = math.sin(math.radians(dip90))*(minonperp-10)
            if sdepth>0:
                plotsz *= -1
                plotsx *= -1
            plotdx = math.cos(math.radians(dip90))*(maxonperp+10)*-1
            plotdz = math.sin(math.radians(dip90))*(maxonperp+10)*-1
            
            perppdfx = math.cos(math.radians(dip90))*PDF_df['depths'].values
            perppdfx = math.cos(math.radians(dip90))*PDF_df['depths'].values
            
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(121)
            
            thispoint = ax1.plot([0], [0], 'ro', label='Node Location')

            dip = ax1.plot([-hd*50, hd*100], [-zd*50, zd*100], 'g-', label='Dip')
            ax1.plot([plotdx, plotsx], [plotdz, plotsz], 'b-', label='PDF Axis')
            ax1.plot([plotdx, plotsx], [plotdz, plotsz], 'bs')
            
            trimmed['lonplot'] = trimmed['lon'].values-lon
            trimmed['latplot'] = trimmed['lat'].values-lat
            trimmed['depplot'] = trimmed['depth'].values-loc_depth
            
            if len(BAframe)>0:
                BA2 = trimmed[trimmed.etype == 'BA']
                bap = ax1.plot(BA2['cosdistance'].values, BA2['depplot'].values, 'r.', label='BA')
            if len(EQframe)>0:
                EQ2 = trimmed[trimmed.etype == 'EQ']
                eqp = ax1.plot(EQ2['cosdistance'].values, EQ2['depplot'].values, 'c.', label='EQ')
            if len(ERframe)>0:
                ER2 = trimmed[trimmed.etype == 'ER']
                erp = ax1.plot(ER2['cosdistance'].values, ER2['depplot'].values, 'y.', label='ER')
            if len(AAframe)>0:
                AA2 = trimmed[trimmed.etype == 'AA']
                aap = ax1.plot(AA2['cosdistance'].values, AA2['depplot'].values, 'k.', label='AA')
            if len(ASframe)>0:
                AS2 = trimmed[trimmed.etype == 'AS']
                asp = ax1.plot(AS2['cosdistance'].values, AS2['depplot'].values, 'm.', label='AS')
            if len(TOframe)>0:
                TO2 = trimmed[trimmed.etype == 'TO']
                top = ax1.plot(TO2['cosdistance'].values, TO2['depplot'].values, 'g.', label='TO')
            if len(RFframe)>0:
                RF2 = trimmed[trimmed.etype == 'RF']
                rfp = ax1.plot(RF2['cosdistance'].values, RF2['depplot'].values, 'b.', label='RF')
            if len(CPframe)>0:
                CP2 = trimmed[trimmed.etype == 'CP']
                CPp = ax1.plot(CP2['cosdistance'].values, CP2['depplot'].values, color='orange', marker='.', label='CP')
            
            ax1.plot(xplot, yplot, 'k-', label='PDFx100', linewidth=2)
            ax1.plot([-60, 60], [peak_depthz, peak_depthz], 'r--', label='Peak Depth')
            ax1.plot([-peak_depthx, -peak_depthx], [-100, 100], 'r--')
            ax1.plot([peakxplot, -peak_depthx], [peakyplot, peak_depthz], 'r-')
            ax1.set_xlabel('horizontal distance outboard <- -> inboard')
            ax1.set_ylabel('Depth deeper <- -> shallower')
            ax1.axis('equal')
            ax1.invert_yaxis()
            plt.grid()
            title = 'Lat: %.2f, Lon: %.2f, Strike: %.2f, Dip: %.2f, Origin Depth: %.2f' % (lat, lon, cstr, cdip, loc_depth)
            ax1.set_title(title)
            lontit = lon*100
            lattit = lat*100
            ax1.legend(loc='best')

            a2 = (lat-trimmed['lat'])*(lat-trimmed['lat'])
            b2 = (lon-trimmed['lon'])*(lon-trimmed['lon'])
            c = np.sqrt(a2+b2)/2
            
            ax2 = fig.add_subplot(122)
            if len(BAframe)>0:
                BAa2 = (lat-BAframe['lat'])*(lat-BAframe['lat'])
                BAb2 = (lon-BAframe['lon'])*(lon-BAframe['lon'])
                BAc = np.sqrt(BAa2+BAb2)/2
                bap = ax2.plot(BAc, BAframe['depth'].values, 'r.', label='BA')
            if len(EQframe)>0:
                EQa2 = (lat-EQframe['lat'])*(lat-EQframe['lat'])
                EQb2 = (lon-EQframe['lon'])*(lon-EQframe['lon'])
                EQc = np.sqrt(EQa2+EQb2)/2
                eqp = ax2.plot(EQc, EQframe['depth'].values, 'c.', label='EQ')
            if len(ERframe)>0:
                ERa2 = (lat-ERframe['lat'])*(lat-ERframe['lat'])
                ERb2 = (lon-ERframe['lon'])*(lon-ERframe['lon'])
                ERc = np.sqrt(ERa2+ERb2)/2
                erp = ax2.plot(ERc, ERframe['depth'].values, 'y.', label='ER')
            if len(AAframe)>0:
                AAframe.loc[AAframe.lon < 0, 'lon']+=360
                AAa2 = (lat-AAframe['lat'])*(lat-AAframe['lat'])
                AAb2 = (lon-AAframe['lon'])*(lon-AAframe['lon'])
                AAc = np.sqrt(AAa2+AAb2)/2
                aap = ax2.plot(AAc, AAframe['depth'].values, 'k.', label='AA')
            if len(ASframe)>0:
                ASa2 = (lat-ASframe['lat'])*(lat-ASframe['lat'])
                ASb2 = (lon-ASframe['lon'])*(lon-ASframe['lon'])
                ASc = np.sqrt(ASa2+ASb2)/2
                asp = ax2.plot(ASc, ASframe['depth'].values, 'm.', label='AS')
            if len(TOframe)>0:
                TOa2 = (lat-TOframe['lat'])*(lat-TOframe['lat'])
                TOb2 = (lon-TOframe['lon'])*(lon-TOframe['lon'])
                TOc = np.sqrt(TOa2+TOb2)/2
                top = ax2.plot(TOc, TOframe['depth'].values, 'g.', label='TO')
            if len(RFframe)>0:
                RFa2 = (lat-RFframe['lat'])*(lat-RFframe['lat'])
                RFb2 = (lon-RFframe['lon'])*(lon-RFframe['lon'])
                RFc = np.sqrt(RFa2+RFb2)/2
                rfp = ax2.plot(RFc, RFframe['depth'].values, 'b.', label='RF')
            if len(CPframe)>0:
                CPa2 = (lat-CPframe['lat'])*(lat-CPframe['lat'])
                CPb2 = (lon-CPframe['lon'])*(lon-CPframe['lon'])
                CPc = np.sqrt(CPa2+CPb2)/2
                CPp = ax2.plot(CPc, CPframe['depth'].values, color='orange', marker='.', label='CP')
        
            if sdepth<0:
                sdepth *= -1
            ax2.plot((0.1, 0.1), (loc_depth-sdepth, ddepth+loc_depth), 'b-')
            ax2.plot((0, 0.2), (loc_depth-sdepth, loc_depth-sdepth), 'b-')
            rangep = ax2.plot((0, 0.2), (ddepth+loc_depth, ddepth+loc_depth), 'b-', label='depthrange')
            locp = ax2.plot((0, np.max(c)), (loc_depth, loc_depth), 'g-', label='Slab1')
            pdfp = ax2.plot(PDF_df['Summed_Values'].values, PDF_df['depths'].values, linewidth=2, color='k', label='PDF')
            pkp = ax2.plot([peak, peak], [plotsz, plotdz], 'r--')
            pkp = ax2.plot([0, 0.5], [peak_depth, peak_depth], 'r--', label='Peak Depth')
            x1, x2, y1, y2 = ax2.axis()
            xmax = max(np.max(c), peak)
            ax2.axis((0, xmax, y1, y2))
            ax2.invert_yaxis()
            ax2.set_xlabel('Probability (PDF) Degree Distance from Node/2 (data)')
            ax2.set_ylabel('Depth')
            title = 'Lat: %.4f, Lon: %.4f, NID: %.4f' % (lat, lon, nID)
            ax2.set_title(title)
            ax2.grid()
            plt.legend(loc='best')
            lontit = lon*100
            lattit = lat*100
            figtitle = 'Output/PDF%s/%spdf%i.png' % (slab, whichpdf, nID)
            fig.savefig(figtitle)
            plt.close()
            filetitle = 'Output/PDF%s/%sused%i.csv' % (slab, whichpdf, nID)
            trimmed.to_csv(filetitle, header=True, index=False, float_format='%0.2f', na_rep = float('nan'))

    # If there is only one event, we do not solve for the depth at that point unless it is AA, BA, or AS
    elif len(elistBAe) > 0 or len(elistASe) > 0 or len(elistAAe) > 0 or len(elistRFe) > 0 or len(elistTOe) or len(elistCPe)> 0:
        frames = [elistBAe, elistASe, elistAAe, elistRFe, elistTOe, elistCPe]
        trimmed_once = pd.concat(frames)
        all_depths = trimmed_once['depth'].values
        variance1 = trimmed_once['unc'].values
        peak_depth = np.mean(all_depths)
        peaklat = lat
        peaklon = lon
        stdv = np.mean(variance1)
        test = True
        n = 1
    else:
        peak_depth = np.nan
        peaklat = lat
        peaklon = lon
        stdv = np.nan
        test = False
        n = 0

    try:
        return peaklon, peaklat, peak_depth, stdv, test, n, multipeaks, centsurf, removelist
    except:
        return peaklon, peaklat, peak_depth, stdv, test, n, multipeaks, stdv, pd.DataFrame()

def doublecheckEREQ(elist, lon, lat):
    if len(elist) == 2:
        if len(elist[elist.etype == 'EQ'])>0:
            if len(elist[elist.etype == 'ER'])>0:
                EQevent = elist[elist.etype == 'EQ']
                ERevent = elist[elist.etype == 'ER']
                
                EQlon = EQevent['lon'].values[0]
                EQlat = EQevent['lat'].values[0]
                EQdep = EQevent['depth'].values[0]
                
                ERlon = ERevent['lon'].values[0]
                ERlat = ERevent['lat'].values[0]
                ERdep = ERevent['depth'].values[0]
                
                if abs(EQlon-ERlon)<0.1 and abs(EQlat-ERlat)<0.1 and abs(EQdep-ERdep)<1:
                    #print ('removed EQ bc ER and EQ same event',lon, lat, elist)
                    return ERevent, False
                    
                else:
                    return elist, True
            else:
                return elist, True
        else:
            return elist, True
    else:
        return elist, True
        
def refilter4(locstr,stkthresh,lon,lat,elist,alen,blen,slab1guide1,slab,testprint):

    if len(elist[elist.etype == 'BA'])>1:
        return alen
    
    elistRFAS = elist[(elist.etype == 'RF')|(elist.etype == 'AS')|(elist.etype == 'CP')]
    elist = elist[(elist.etype != 'RF')&(elist.etype != 'AS')&(elist.etype != 'CP')]
    if len(elist)<1:
        return alen
    
    ellons = np.arange(lon-1.3,lon+1.2,0.01)
    ellats = np.arange(lat-1.3,lat+1.2,0.01)
    
    elon, elat = np.meshgrid(ellons,ellats)
    
    elon2 = elon.flatten()
    elat2 = elat.flatten()
    
    rlon, rlat = heading(lon, lat, alen, locstr)
    distance, az = npcosine(rlon, rlat, elon2, elat2)

    erta = math.sqrt(1-((math.pow(blen, 2))/(math.pow(alen, 2))))
    mdist = getEllipseRad(alen, erta, az, locstr)

    #elon2 = elon2[np.isfinite(distance)]
    #elat2 = elat2[np.isfinite(distance)]
    #distance = distance[np.isfinite(distance)]
    #distance = distance[np.isfinite(distance)]
    
    elon22 = elon2[distance <= mdist] # caused runtime
    elat22 = elat2[distance <= mdist] # caused runtime
    dist2 = distance[distance <= mdist] # caused runtime
    cdist,az2 = npcosine(lon,lat,elon22,elat22)
    
    evarr = np.zeros((len(elon22),5))
    evarr[:,0] = elon22
    evarr[:,1] = elat22
    evarr[:,4] = cdist
    
    try:
        slab1guide = slab1guide1[(slab1guide1[:,0] < lon+3)&(slab1guide1[:,0] > lon-3)&(slab1guide1[:,1] < lat+3)&(slab1guide1[:,1] > lat-3)]
        evarr[:,3] = griddata(slab1guide[:, 0:2], slab1guide[:, 3], evarr[:, 0:2], method='nearest')
        locstr = evarr[:,3][((evarr[:,0] <= lon+0.01)&(evarr[:,0] >= lon-0.01))&((evarr[:,1] <= lat+0.01)&(evarr[:,1] >= lat-0.01))][0]
    except:
        try:
            slab1guide = slab1guide1[(slab1guide1[:,0] < lon+10)&(slab1guide1[:,0] > lon-10)&(slab1guide1[:,1] < lat+10)&(slab1guide1[:,1] > lat-10)]
            evarr[:,3] = griddata(slab1guide[:, 0:2], slab1guide[:, 3], evarr[:, 0:2], method='nearest')
            locstr = evarr[:,3][((evarr[:,0] <= lon+0.01)&(evarr[:,0] >= lon-0.01))&((evarr[:,1] <= lat+0.01)&(evarr[:,1] >= lat-0.01))][0]
        except:
            if slab == 'sam' and lat > 5:
                if testprint:
                    print ('so far from slab, but kept alen',lon,lat)
                return alen
            else:
                if testprint:
                    print ('so far from slab reduced to blen',lon,lat)
                return blen

    maxstrike = np.max(evarr[:,3])
    minstrike = np.min(evarr[:,3])
    if (maxstrike > 270 and minstrike < 90) or locstr<45:
        evarr[:,3][evarr[:,3]<90]+=360
        if locstr < 90:
            locstr += 360
        maxstrike = np.max(evarr[:,3])
        minstrike = np.min(evarr[:,3])

    lontest = evarr[:,0]*1.0
    lattest = evarr[:,1]*1.0
    strtest = evarr[:,3]*1.0
    if slab == 'sam' or (slab == 'sum' and lat > 15) or slab == 'man' or slab == 'van':
        evarr[:,3][evarr[:,3]<90]+=360
    
    evstd = np.std(evarr[:,3])
    #print (lon, lat, locstr, np.min(evarr[:,3]), np.max(evarr[:,3]),np.mean(evarr[:,3]),evstd)
    distcut = 0
    while evstd > stkthresh and (alen-distcut)>blen:
        distcut += 1
        evarr = evarr[evarr[:,4]<(alen-distcut)]
        evstd = np.std(evarr[:,3])
    
    if testprint:
        print ('lon,lat,alen,blen,evstd,locstr,stkthresh',lon,lat,alen-distcut,blen,evstd,locstr,stkthresh)
    
    return alen-distcut


# Eliminates events in dfo that are found in dfm
def removematches(dfo, dfm):
    ind = (dfo.lon.isin(dfm.lon)) & (dfo.depth.isin(dfm.depth)) & (dfo.lon.isin(dfm.depth))
    dfo0 = dfo[~ind]
    dfo1 = dfo[ind]
    return dfo0, dfo1

def getTrenchInPolygon(slabname, trench, polygonFile):

    #####################################
    #written by Maria Furtney, 7/20/2016#
    #####################################

    ''' creates a grid of 1 or nan based on if they are within a clipping mask or not. DEP.6.29.16 '''
    ''' modified to fit this script by MAF 7/18/16 '''

    ### Input:
    # slabname: a 3 digit character code identifying a slab region 
    #data: the input data which may or may not be within the polygon

    ### Output:
    #contained_data: an array of coordinate pairs (lon,lat) that reside within the polygon region 
    #check if slabbounds are already defined. If not, acquire them
    slabbounds = slabpolygon(slabname, polygonFile)
       
    #slabbbounds come in lon1,lat1,lon2,lat2... format
    #even numbers are then longitudes while odds are latitudes
    coords = np.size(slabbounds)

    #simple even/odd function 
    def is_odd(num):
        return num & 0x1

    lons = []
    lats = []
    
    for i in range(coords):
        val = slabbounds[i][1:]
        if is_odd(i):
            lats.append(val)
        else: 
            lons.append(val)

    trench1 = zerothreesixty(trench)
    data = list(zip(trench1['lon'].values, trench1['lat'].values))

    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(data[:])
    mask1 = np.zeros(len(temp),)*np.nan 
    mask1[temp] = 1
    keepers = []
    for i in range(len(data)):
        points_in_poly = np.dot(mask1[i], data[i])
        if i > 0:
            keepers = np.vstack((keepers, points_in_poly))
        else: 
            keepers = points_in_poly
     
    rows_to_drop = []
    for i in range(len(keepers)):
        if np.isnan(keepers[i][0]) == True:
            rows_to_drop.append(i)

    data_to_keep = trench.drop(trench.index[rows_to_drop])
    return data_to_keep

def npcosine(lon1, lat1, lon2, lat2):

    # Written GLM 4.25.17

    ''' Arguments:  lon1 - longitude point that the angle is referenced from
                            (float)[deg]
                    lat1 - latitude point that the angle is referenced from
                            (float)[deg]
                    lon2 - array of longitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                    lat2 - array of latitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                        
        Returns:    dist - array of great circle distance between the two
                            lat/lon points (arr of floats)[km]
                    ang - array of angles between the two points (clockwise from
                            0 degrees from lat1/lon1 point) 
                            (arr of floats)[deg] '''

    # Creating degrees/radians conversion constants
    d2r = (math.pi/180.0)
    r2d = (180.0/math.pi)
    ddlon = lon1-lon2
    
    # Ensuring that azimuths are between 0 and 360
    if lon1 < 0.0:
        lon1 += 360.0
    lon2[lon2<0.0] += 360.0

    # Getting distance and angle between the two points (in degrees)
    dist, ang = npcosrule(d2r, lon1, lat1, lon2, lat2)
    ang[(lon1>lon2)&(ddlon<180.0)] = 2*math.pi-ang[(lon1>lon2)&(ddlon<180.0)] # causes runtime
    dist = np.abs(dist * r2d)

    dist[dist > 180.0] = 360-dist[dist > 180] # causes runtime
    ang[dist > 180.0] += math.pi # causes runtime
    ang[ang > 2.0*math.pi] = 2.0*math.pi - ang[ang > 2.0*math.pi] # causes runtime
    dist *= 111.19
    ang *= r2d
    
    lon2[lon2<0]+=360

    return dist, ang


def npcosrule(d2r, lon1, lat1, lon2, lat2):

    # Written GLM 4.25.17

    ''' Arguments:  d2r - degree to radians conversion constant (float)
                    lon1 - longitude point that the angle is referenced from
                            (float)[deg]
                    lat1 - latitude point that the angle is referenced from
                            (float)[deg]
                    lon2 - array of longitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                    lat2 - array of latitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
        
        Returns:    dist2 - array of great circle distance between the two 
                            lat/lon points (arr of floats)[km]
                    ang - array of angles between the two points (clockwise from
                            0 degrees from lat1/lon1 point) 
                            (arr of floats)[deg] '''

    # breaks when lat1==lat2 or lon1==lon2. Add small difference where needed
    londiff = np.abs(lon2-lon1)
    latdiff = np.abs(lat2-lat1)
    lon2[londiff<0.0001] += 0.0001 # causes runtime
    lat2[latdiff<0.0001] += 0.0001 # causes runtime

    cl1 = (90.0-lat1)*d2r
    cl2 = (90.0-lat2)*d2r
    dlon = (lon2-lon1)*d2r

    coscl2 = np.cos(cl2)
    sincl2 = np.sin(cl2)
    cosdlon = np.cos(dlon)
    
    coscl1 = math.cos(cl1)
    sincl1 = math.sin(cl1)

    dist = (coscl1 * coscl2) + (sincl1 * sincl2 * cosdlon)

    dist[dist < -1] = -1.0 # causes runtime
    dist[dist > 1] = 1.0 # causes runtime

    dist2 = np.arccos(dist)
    dist2[dlon > math.pi] = 2*math.pi - dist2[dlon > math.pi] # causes runtime
    
    ang = np.zeros(len(dist))
    num = np.zeros(len(dist))
    den = np.zeros(len(dist))

    num[dist != 0] = (coscl2[dist != 0] - (dist[dist != 0] * coscl1))
    den[dist != 0] = (np.sin(dist2[dist != 0]) * sincl1)

    ang[dist != 0] = num[dist != 0] / den[dist != 0]
    ang[dist == 0] = 1.0
    ang[ang < -1] = -1.0 # causes runtime
    ang[ang > 1] = 1.0 # causes runtime
    ang2 = np.arccos(ang)
    
    return dist2, ang2

def getlatloncutoffs(lons,lats,eventlist, testprint):
    
    lonmean = eventlist['lon'].mean()
    latmean = eventlist['lat'].mean()
    
    NWlons = lons[(lons<lonmean)&(lats>=latmean)]
    NWlats = lats[(lons<lonmean)&(lats>=latmean)]
    NWelist = eventlist[(eventlist.lon<lonmean+2)&(eventlist.lat>latmean-2)]
    
    SWlons = lons[(lons<lonmean)&(lats<latmean)]
    SWlats = lats[(lons<lonmean)&(lats<latmean)]
    SWelist = eventlist[(eventlist.lon<lonmean+2)&(eventlist.lat<latmean+2)]
    
    SElons = lons[(lons>=lonmean)&(lats<latmean)]
    SElats = lats[(lons>=lonmean)&(lats<latmean)]
    SEelist = eventlist[(eventlist.lon>lonmean-2)&(eventlist.lat<latmean+2)]
    
    NElons = lons[(lons>=lonmean)&(lats>=latmean)]
    NElats = lats[(lons>=lonmean)&(lats>=latmean)]
    NEelist = eventlist[(eventlist.lon>lonmean-2)&(eventlist.lat>latmean-2)]
    
    NWelist = NWelist.reset_index(drop=True)
    SWelist = SWelist.reset_index(drop=True)
    SEelist = SEelist.reset_index(drop=True)
    NEelist = NEelist.reset_index(drop=True)
    
    listoflons = [NWlons,SWlons,SElons,NElons]
    listoflats = [NWlats,SWlats,SElats,NElats]
    listofelists = [NWelist,SWelist,SEelist,NEelist]
    
    if testprint:
        print ('lenghts of arrays',len(NWlons),len(SWlons),len(SElons),len(NElons))

    return listoflons,listoflats,listofelists
    #return [lons],[lats],[eventlist]

def getlatloncutoffs2(lons,lats,eventlist):

    listoflons,listoflats,listofelists = getlatloncutoffs(lons,lats,eventlist)
    secondlonlist,secondlatlist,secondelist = [],[],[]
    for cut in range(len(listoflons)):
        theselons = listoflons[cut]
        theselats = listoflats[cut]
        theseevents = listofelists[cut]
        thislonlist,thislatlist,thiselist = getlatloncutoffs(theselons,theselats,theseevents)
        secondlonlist.extend(thislonlist)
        secondlatlist.extend(thislatlist)
        secondelist.extend(thiselist)

    return secondlonlist,secondlatlist,secondelist

def makeReference(slab1data,lons,lats,grid,testprint,slab):

    slab1query = np.zeros((len(slab1data),5))
    try:
        slab1query[:,0] = slab1data['lon'].values
        slab1query[:,1] = slab1data['lat'].values
        slab1query[:,2] = slab1data['depth'].values
        slab1query[:,3] = slab1data['strike'].values
        slab1query[:,4] = slab1data['dip'].values
    except:
        slab1query[:,0] = slab1data[:,0]
        slab1query[:,1] = slab1data[:,1]
        slab1query[:,2] = slab1data[:,2]
        slab1query[:,3] = slab1data[:,3]
        slab1query[:,4] = slab1data[:,4]

    #np.savetxt('%s_diptest1.csv'%slab, slab1query, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
    slab1guide = np.zeros((len(lons),6))
    slab1guide[:,0] = lons
    slab1guide[:,1] = lats
    if slab == 'izuz' or slab == 'japz' or slab == 'puyz' or slab == 'solz':
        slab1guide[:,2] = griddata(slab1query[:, 0:2], slab1query[:, 2], slab1guide[:, 0:2], method='linear')
        slab1guide[:,3] = griddata(slab1query[:, 0:2], slab1query[:, 3], slab1guide[:, 0:2], method='nearest')
        slab1guide[:,4] = griddata(slab1query[:, 0:2], slab1query[:, 4], slab1guide[:, 0:2], method='nearest')
        
        #np.savetxt('%s_diptest2.csv'%slab, slab1guide, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
        noguide = slab1guide[np.isnan(slab1guide[:,2])]
        yesguide = slab1guide[np.isfinite(slab1guide[:,2])]
        noguide[:,2] = griddata(slab1query[:, 0:2], slab1query[:, 2], noguide[:, 0:2], method='nearest')
        slab1guide = np.vstack((yesguide,noguide))
        #np.savetxt('%s_diptest3.csv'%slab, slab1guide, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
    else:
        slab1guide[:,2] = griddata(slab1query[:, 0:2], slab1query[:, 2], slab1guide[:, 0:2], method='nearest')
        slab1guide[:,3] = griddata(slab1query[:, 0:2], slab1query[:, 3], slab1guide[:, 0:2], method='nearest')
        slab1guide[:,4] = griddata(slab1query[:, 0:2], slab1query[:, 4], slab1guide[:, 0:2], method='nearest')
        #np.savetxt('%s_diptest3.csv'%slab, slab1guide, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
    slab1guide[:,2] *= -1.0
    slab1guide = np.round(slab1guide,decimals=1)
    if testprint:
        print ('slab1guide',slab1guide)
    
    return slab1guide,slab1query
    
def getslab12(slab1guide,slab1query,lon,lat,grid,depgrid,strgrid,dipgrid,testprint,TRdata,meanBA,slab):

    # If Slab1.0 exists at this lat,lon, collect the local depth, strike,
    # and dip. Takes different longitude formats into account for specific slabs

    if slab == 'phi' or slab == 'sol' or slab == 'man' or slab == 'him':
        depthresh = 30
    else:
        depthresh = 70

    if slab == 'man' or slab == 'phi' or slab == 'png' or slab == 'sul' or slab == 'cot' or slab == 'car' or slab == 'hel' or slab == 'ita' or slab == 'puy' or slab == 'mak' or slab == 'cam' or slab == 'pan' or slab == 'mue' or slab == 'sco' or slab == 'ryu' or slab == 'him':
        slab1, strtmp, diptmp, inside, extended, out, extlon, extlat = extendinginterp(slab1guide,lon,lat,slab1query,grid,TRdata,meanBA,testprint,False,depthresh,slab)
        return slab1, strtmp, diptmp, inside, extended, out, extlon, extlat
    
    try:
        try:
            slab1 = depgrid.getValue(lat, lon) * -1.0
            strtmp = strgrid.getValue(lat, lon)
            diptmp = dipgrid.getValue(lat, lon)
        except:
            slab1 = depgrid.getValue(lat, lon-360) * -1.0
            strtmp = strgrid.getValue(lat, lon-360)
            diptmp = dipgrid.getValue(lat, lon-360)
        slab1 = np.max(slab1) # gets the value out of the array
        strtmp = np.max(strtmp) # gets the value out of the array
        diptmp = np.max(diptmp) # gets the value out of the array

        if np.isfinite(slab1) and np.isfinite(strtmp) and np.isfinite(diptmp):
            if testprint:
                print ('0',lon,lat,slab1,strtmp,diptmp)
            return slab1, strtmp, diptmp, True, False, False, lon, lat
        else:
            slab1, strtmp, diptmp, inside, extended, out, extlon, extlat = extendinginterp(slab1guide,lon,lat,slab1query,grid,TRdata,meanBA,testprint,False,depthresh,slab)
            return slab1, strtmp, diptmp, inside, extended, out, extlon, extlat
            
    except:
        slab1, strtmp, diptmp, inside, extended, out, extlon, extlat = extendinginterp(slab1guide,lon,lat,slab1query,grid,TRdata,meanBA,testprint,False,depthresh,slab)
        return slab1, strtmp, diptmp, inside, extended, out, extlon, extlat
    
    
def extendinginterp(slab1guide,lon,lat,slab1query,grid,TRdata,meanBA,testprint,interp,depthresh,slab):
    
    if len(TRdata)>0 and (slab != 'sol' or lon > 150):
        loc_tr = TRdata[(TRdata.lon > lon-3) & (TRdata.lon < lon+3) & (TRdata.lat > lat-3) & (TRdata.lat < lat+3)]
        if len(loc_tr)>0:
            #loc_tr['dist'] = gps2dist_azimuth(lat, lon, loc_tr['lat'], loc_tr['lon'])[0]/1000.0
            loc_tr['dist'], tempangles = npcosine(lon, lat, loc_tr['lon'].values, loc_tr['lat'].values)
            mindist = loc_tr['dist'].min()
            loc_tr = loc_tr[loc_tr.dist == mindist]
            lonT = loc_tr['lon'].values[0]
            latT = loc_tr['lat'].values[0]
            azT = loc_tr['az'].values[0]
            thisdist, thisang, latB, lonB = cosine(lonT, latT, lon, lat)
            out = isoutboard(azT, thisang)
            if out:
                if testprint:
                    print ('outboard trench lon,lat,lonT,latT,azT,thisdist',lon,lat,lonT,latT,azT,thisdist)
                return meanBA, azT, 0.0, False, True, True, lon, lat

    if testprint:
        print ('lon,lat',lon,lat)
    thisguide = slab1guide[(slab1guide[:,0] == lon)&(slab1guide[:,1] == lat)]
    thisdepth = thisguide[0,2]
    thisstrike = thisguide[0,3]
    thisdip = thisguide[0,4]
    buffval = 1.5
    if slab == 'sol' or slab == 'ker' or slab == 'izu' or slab == 'hin' or slab == 'pam' or slab == 'man':
        buffval = 0.5
    thisquery1 = slab1query[slab1query[:,0]<lon+grid*buffval]
    thisquery1 = thisquery1[thisquery1[:,0]>lon-grid*buffval]
    thisquery1 = thisquery1[thisquery1[:,1]<lat+grid*buffval]
    thisquery1 = thisquery1[thisquery1[:,1]>lat-grid*buffval]
    
    if len(thisquery1)>0:
        if testprint:
            print ('1',lon,lat,thisdepth,thisstrike,thisdip,thisguide[0,0],thisguide[0,1])
        return thisdepth, thisstrike, thisdip, False, True, False, lon, lat

    if interp:
        thisquery = slab1query[slab1query[:,0]<lon+10]
        thisquery = thisquery[thisquery[:,0]>lon-10]
        thisquery = thisquery[thisquery[:,1]<lat+10]
        thisquery = thisquery[thisquery[:,1]>lat-10]
    else:
        thisquery = slab1query[slab1query[:,0]<lon+2]
        thisquery = thisquery[thisquery[:,0]>lon-2]
        thisquery = thisquery[thisquery[:,1]<lat+2]
        thisquery = thisquery[thisquery[:,1]>lat-2]
        if slab == 'izu':
            if len(thisquery) < 1 and lat > 23 and lat < 27:
                thisquery = slab1query[slab1query[:,0]<lon+5]
                thisquery = thisquery[thisquery[:,0]>lon-5]
                thisquery = thisquery[thisquery[:,1]<lat+5]
                thisquery = thisquery[thisquery[:,1]>lat-5]
    
    if len(thisquery)<1:
        if testprint:
            print ('2',lon,lat,thisdepth,thisstrike,thisdip)
        return np.nan, np.nan, np.nan, False, False, False, lon, lat

    if slab == 'izuz':
        if lat < 27 and lat > 23:
            thisdip *= 1.5
            if thisdip > 85:
                thisdip = 85
                
    distances,cosangles = npcosine(lon,lat,thisquery[:,0],thisquery[:,1])
    mindist = np.min(distances)
    mincosangle = cosangles[distances == mindist][0]
    minlon = thisquery[:,0][distances == mindist][0]
    minlat = thisquery[:,1][distances == mindist][0]
    cosangle = mincosangle+90
    disttest = 111.19*grid
    
    anglediff = abs(mincosangle-thisstrike)

    if anglediff>55 and anglediff<125 and thisdepth > 100:
        depthadd = mindist * math.tan(thisdip*math.pi/180.0)
        thisdepth += depthadd
        if testprint:
            print ('3',lon,lat,thisdepth,thisstrike,thisdip, minlon, minlat, mindist)
        return thisdepth, thisstrike, thisdip, False, True, False, minlon, minlat
    
    elif anglediff<=15:
        if testprint:
            print ('4',lon,lat,thisdepth,thisstrike,thisdip)
        return np.nan, np.nan, np.nan, False, False, False, lon, lat

    elif mindist<8*disttest and thisdepth > depthresh:
        depthadd = mindist * math.tan(thisdip*math.pi/180.0)
        thisdepth += depthadd
        if testprint:
            print ('5',lon,lat,thisdepth,thisstrike,thisdip)
        return thisdepth, thisstrike, thisdip, False, True, False, minlon, minlat
    else:
        if testprint:
            print ('6',lon,lat,thisdepth,thisstrike,thisdip)
        return thisdepth, thisstrike, thisdip, False, False, True, lon, lat
        #return np.nan, np.nan, np.nan, False, False, False

def npheading(lon, lat, dist, az):
    
    ''' Arguments:  lon - longitude of known point (array)
                    lat - latitude of known point (array)
                    dist - distance from known point (array)
                    az - azimuth from known point (array)
        
        Returns:    lat - latitude of new point as projected by a certain azimuth and great circle
                            distance from known lat/lon point
                    lon - longitude of new point as projected by a certain azimuth and great circle
                            distance from known lat/lon point  '''
    
    # Creating degrees/radians conversion constants
    d2r = math.pi/180
    r2d = 180/math.pi
    
    # Ensuring that distances are positive
    az[dist < 0] -= 180.0
    dist[dist < 0] *= -1.0
    
    # Ensuring that azimuths are between 0 and 360
    az[az < 0] += 360
    az[az > 0] -= 360
    
    # Finding projected latitude
    b = (90 - lat) * d2r
    a = (dist / 111.19) * d2r
    angC = az * d2r
    c = np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * np.cos(angC)
    c = np.arccos(c)
    cdeg = c * r2d
    lat1 = 90 - cdeg
    
    # Finding projected longitude
    angA = (np.cos(a) - (np.cos(b) * np.cos(c))) / (np.sin(b) * np.sin(c))
    angA[angA > 1.0] -= 0.00001
    angA[angA < -1.0] += 0.00001

    angA = np.arccos(angA)
    adeg = angA * r2d

    lon1 = np.copy(lon)
    lon1[(az > 0) & (az <= 180)] += adeg[(az > 0) & (az <= 180)]
    lon1[(az <= 0) | (az > 180)] -= adeg[(az <= 0) | (az > 180)]
    
    return lon1, lat1

def cmtfilter(data,seismo_thick,printtest,datainfo,slab):

    ''' Arguments:  data - data with all shallow/nonshallow and thrust/nonthrust earthquake
    
        Returns:    filtered - fitered dataframe which DEPENDS ON WHAT YOU DO/DONT COMMENT OUT
        
                                (1) filters only shallow earthquakes that have MT criteria which are non thrust
                                all other shallow earthquakes WITHOUT MT info are NOT filtered
                                
                                OR
                                
                                (2) filters ALL shallow earthquakes UNLESS they have MT info and that
                                MT info has the criteria of a thrust event. '''
    
    if slab == 'aluz' or slab == 'casz':
        datanotEQ = data[data.etype != 'EQ']
        data = data[data.etype == 'EQ']
    else:
        datanotEQ = data[(data.etype != 'EQ')&(data.etype != 'ER')]
        data = data[(data.etype == 'EQ')|(data.etype == 'ER')]
    # Removes non-thrust events from depths shallower than seismogenic zone
    deep_data = data[data.depth >= seismo_thick]
    
    # Includes shallow data without MT info (1) - comment out next two lines for (2)
    #dfn = data[np.isnan(data['Paz'])]
    #dfn = dfn[data.depth < seismo_thick]
    
    if printtest:
        filtered = data[np.isnan(data['S1'])]
        addToDataInfo(filtered, 0, 'shallow with no MT info', datainfo,'df')
    
    data = data[np.isfinite(data['S1'])]
    shallow_data = data[data.depth < seismo_thick]
    
    if printtest:
        #if 'Ndip' in shallow_data.columns:
        #    filtered = shallow_data[(shallow_data.Tpl<=50) | (shallow_data.Ndip>30)]
        #else:
        #    filtered = shallow_data[(shallow_data.R1<=30) | (shallow_data.R2<=30)
        #               | (shallow_data.R1>=150) | (shallow_data.R2>=150)]
        filtered = shallow_data[shallow_data.kagan < 35]
        
        addToDataInfo(filtered, 0, 'non thrust and shallow', datainfo,'df')

    # Depending on which MT info are provided, filters non-thrust, shallow events
    #if 'Ndip' in shallow_data.columns:
    #    shallow_data = shallow_data[(shallow_data.Tpl>50) & (shallow_data.Ndip<=30)]
    #else:
    #    shallow_data = shallow_data[(shallow_data.R1>30) & (shallow_data.R2>30)
    #                   & (shallow_data.R1<150) & (shallow_data.R2<150)]
    shallow_data = shallow_data[shallow_data.kagan < 35]
    
    # Includes shallow data without MT info (1) - comment out next line for (2)
    # filtered = pd.concat([deep_data, shallow_data, dfn, datanotEQ])

    # Only includes shallow thrust events (2) - uncomment line below for (2) and comment necessary lines above
    filtered = pd.concat([deep_data, shallow_data, datanotEQ])
    
    return filtered


def getSZthickness(data,folder,slab,maxdep,maxdepdiff,origorcentl,origorcentd,slaborev,savedir,lengthlim):

    (slab,sl2,date) = folder.split('_')

    # if solomon islands, don't use events on western end for calculation
    if slab == 'sol':
        data = data[data.lon > 148]

    # assign minimum number of data points to use after filters
    if slab == 'hel' or slab == 'car' or slab == 'mak':
        mindata = 40
    else:
        mindata = 20

    # get all events with moment tensors and filter by mag and depth difference
    alldata = data[np.isfinite(data.S1)]
    alldata = alldata[(alldata.mag < 10.0)&(alldata.depdiff < maxdepdiff)&(alldata.depdiff > -1*maxdepdiff)] # was (alldata.mag < 7.5)
    alldata = alldata[alldata.mdep<=maxdep]

    # filter events to only include 30 < rake < 150
    alldata = alldata[((alldata.R1>45) & (alldata.R2>45)
                   & (alldata.R1<135) & (alldata.R2<135))]

    # filter dataset by kagans angle, if not enough data points exist, broaden kagans angle filter
    maxadd = 60-35
    for i in range(maxadd):
        dat = alldata[alldata.kagan < 35+i]
        if len(dat) > mindata:
            maxkagan = 35+i
            break
    try:
        alldata = alldata[alldata.kagan < maxkagan]
    except:
        print ('not enough events within surface filter')
        print ('broadening maxkagan to max value of 60 degrees')
        # maxkagan = 100
        maxkagan = 60
        alldata = alldata[alldata.kagan < maxkagan]

    # save filtered dataset to file
    alldata = alldata[['lat','lon','depth','unc','etype','ID','mag','time','S1','D1','R1','S2','D2','R2','src','slab2str','slab2dip','slab2rak','mrr','mtt','mpp','mrt','mrp','mtp','mrrS','mttS','mppS','mrtS','mrpS','mtpS','kagan','depdiff','slab2dep','mlon','mlat','mdep']]
    alldata.to_csv('%s/%s_slab2_szt_%s.csv' % (savedir,slab,date),header=True,index=False,na_rep=np.nan)

    # make depth array for gaussian fitting
    if origorcentd == 'c':
        depths = alldata.mdep.as_matrix()
    else:
        depths = alldata.depth.as_matrix()

    if slaborev == 's':
        depths = alldata.slab2dep.as_matrix()*-1.0

    N = len(depths)
    depths.shape = (N,1)

    if N < 2:
        print ('seismogenic zone thickness was not calculated, not enough data within filters')
        szt = 40
        lsz = 10
        depthlist = []
        norm_pdfdT = 100
        lendata = 0
        if savedir != 'Output/%s_slab2_%s'%(slab,date):
            return szt, lsz, depthlist, norm_pdfdT, lendata
        else:
            return szt, lsz

    # initialize empty list for smooth depths and array of 0s for histogram rep.
    dc = [0 for i in range(maxdep)]
    depth2 = []
    i = 0

    # loop through range of depths
    for i in range(maxdep):

        # make window for smoothing
        mind = i-2
        if mind < 0:
            mind = 0
        maxd = i+2
        if maxd > maxdep:
            maxd = maxdep

        # loop through depth array
        for d in depths:
        
            # if within window, append depth to smooth depths, incrament hist. value
            if d >= mind and d <= maxd:
                dc[i] = dc[i]+1
                depth2.append(i)

    # normalize histogram value by dividing by amount of depths
    i2 = float(len(depth2))
    dc2 = [x / i2 for x in dc]

    # make smooth depth array for gaussian fitting
    N2 = len(depth2)
    depths2 = np.array(depth2)
    depths2.shape = (N2,1)

    # get mean and stdv.
    m0 = np.mean(depths2)
    sd0 = np.std(depths2)

    # calculate normal distribution from mean and stdv.
    depthlist = list(range(maxdep))
    norm_pdf = matplotlib.mlab.normpdf(depthlist,m0,sd0)

    # find bimodal distribution using 3 Gaussians
    clf3 = mixture.GaussianMixture(n_components=3, covariance_type='full')
    clf3.fit(depths2)

    # find bimodal distribution using 2 Gaussians
    clf2 = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf2.fit(depths2)

    # acquire weights, means, and covariances from contributing distributions
    m1, m2, m3 = clf3.means_
    w1, w2, w3 = clf3.weights_
    c1, c2, c3 = clf3.covariances_

    md1, md2 = clf2.means_
    wd1, wd2= clf2.weights_
    cd1, cd2 = clf2.covariances_

    # calculate standard deviations for triple and double distributions
    sd1=np.sqrt(c1)[0]
    sd2=np.sqrt(c2)[0]
    sd3=np.sqrt(c3)[0]
    
    sdd1=np.sqrt(cd1)[0]
    sdd2=np.sqrt(cd2)[0]

    # create summed PDF for triple and double distributions
    norm_pdf1 = matplotlib.mlab.normpdf(depthlist,m1,sd1)*w1
    norm_pdf2 = matplotlib.mlab.normpdf(depthlist,m2,sd2)*w2
    norm_pdf3 = matplotlib.mlab.normpdf(depthlist,m3,sd3)*w3
    norm_pdfT = norm_pdf1 + norm_pdf2 + norm_pdf3

    norm_pdfd1 = matplotlib.mlab.normpdf(depthlist,md1,sdd1)*wd1
    norm_pdfd2 = matplotlib.mlab.normpdf(depthlist,md2,sdd2)*wd2
    norm_pdfdT = norm_pdfd1 + norm_pdfd2
    
    # calculate rms for all distributions
    rms1a = math.sqrt(mean_squared_error(dc2, norm_pdf))
    rms3a = math.sqrt(mean_squared_error(dc2, norm_pdfT))
    rms2a = math.sqrt(mean_squared_error(dc2, norm_pdfdT))

    # make baseline indices and sum
    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999

    # assign percentiles for different dataset lengths
    lendata = len(alldata)
    if lendata < lengthlim:
        lobound = 0.05
        hibound = 0.9
    else:
        lobound = 0.05
        hibound = 0.95

    # loop through triple gaussian PDF and identify depth indices of important percentiles
    for i in range(len(norm_pdfT)):
        sumtest += norm_pdfT[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break

    # if the upper percentile isn't reached, define the upper bound as the deepest depth
    if sumindex == -999:
        sumindex = len(depthlist) - 1

    # get triple gaussian depths from indices
    szT_triple = depthlist[sumindex]
    tap_triple = depthlist[tapindex]
    low_triple = depthlist[lowindex]

    # reset indices and sum
    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999

    # loop through double gaussian PDF and identify depth indices of important percentiles
    for i in range(len(norm_pdfdT)):
        sumtest += norm_pdfdT[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break

    # if the upper percentile isn't reached, define the upper bound as the deepest depth
    if sumindex == -999:
        sumindex = len(depthlist) - 1

    # get double gaussian depths from indices
    szT_double = depthlist[sumindex]
    tap_double = depthlist[tapindex]
    low_double = depthlist[lowindex]

    # reset indices and sum
    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999

    # loop through normal gaussian PDF and identify depth indices of important percentiles
    for i in range(len(norm_pdf)):
        sumtest += norm_pdf[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break

    # if the upper percentile isn't reached, define the upper bound as the deepest depth
    if sumindex == -999:
        sumindex = len(depthlist) - 1

    # get normal gaussian depths from indices
    szT_single = depthlist[sumindex]
    tap_single = depthlist[tapindex]
    low_single = depthlist[lowindex]

    # plot to show fit compared to data
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(111)

    # plot data histogram
    ax1.plot([depthlist[0],depthlist[0]],[0,dc2[0]],linewidth=10,c='k',label='Data')
    for i in range(1,len(dc2)):
        ax1.plot([depthlist[i],depthlist[i]],[0,dc2[i]],linewidth=10,c='k')

    # plot two normal distributions to be summed for double normal dist.
    ax1.plot(depthlist,norm_pdfd1,linewidth=2,c='springgreen',label='m=%.2f, s=%.2f, w=%.2f'%(m1,sd1,w1))
    ax1.plot(depthlist,norm_pdfd2,linewidth=2,c='springgreen',label='m=%.2f, s=%.2f, w=%.2f'%(m2,sd2,w2))

    # plot normal gaussian distribution and double normal distribution
    ax1.plot(depthlist,norm_pdf,label='RMS1: %.4f'%rms1a,linewidth=2,c='y')
    ax1.plot(depthlist,norm_pdfdT,label='RMS2: %.4f'%rms2a,linewidth=2,c='g')

    # plot and label 10th/5th and 90th/95th percentile depths for either distribution
    if lendata < lengthlim:
        ax1.plot([szT_single,szT_single],[0,np.max(dc2)],'r:',label='90th_single (%.4f)'%szT_single,linewidth=2)
        ax1.plot([szT_double,szT_double],[0,np.max(dc2)],'r-.',label='90th_double (%.4f)'%szT_double,linewidth=2)
        ax1.plot([low_single,low_single],[0,np.max(dc2)],'c:',label='10th_single (%.4f)'%low_single,linewidth=2)
        ax1.plot([low_double,low_double],[0,np.max(dc2)],'c-.',label='10th_double (%.4f)'%low_double,linewidth=2)
    else:
        ax1.plot([szT_single,szT_single],[0,np.max(dc2)],'r:',label='95th_single (%.4f)'%szT_single,linewidth=2)
        ax1.plot([szT_double,szT_double],[0,np.max(dc2)],'r-.',label='95th_double (%.4f)'%szT_double,linewidth=2)
        ax1.plot([low_single,low_single],[0,np.max(dc2)],'c:',label='5th_single (%.4f)'%low_single,linewidth=2)
        ax1.plot([low_double,low_double],[0,np.max(dc2)],'c-.',label='5th_double (%.4f)'%low_double,linewidth=2)

    # assign seismogenic zone thickness depth using distribution with best rms fit
    if rms1a < rms2a:
        ax1.plot(szT_single,np.max(dc2),'go',linewidth=5)
        szt, tap, lsz = szT_single, tap_single, low_single
    else:
        ax1.plot(szT_double,np.max(dc2),'go',linewidth=5)
        szt, tap, lsz = szT_double, tap_double, low_double

    if len(depthlist) < 50:
        print ('less than 50 constraining events, setting default values ...')
        szt, tap, lsz = 40, 10, 10

    # normalize plot axes. assign labels, title, and legend
    ax1.set_xlim([0,65])
    ax1.legend(loc='best')
    ax1.grid()
    ax1.set_xlabel('Depths')
    ax1.set_ylabel('P')
    ax1.set_title('Depth distribution (%s) %i EQs (surfacefilt = %i km, kaganfilt = %i deg, orig = %s, depth = %s, hist= %s)'%(slab,len(alldata),maxdepdiff,maxkagan,origorcentl,origorcentd,slaborev))
    sztdata = alldata[(alldata.depth > lsz)&(alldata.depth < szt)]
    meanslabdip = sztdata['slab2dip'].mean()
    sztdata['dipdiff1'] = np.abs(sztdata['slab2dip'].values-sztdata['D1'].values)
    sztdata['dipdiff2'] = np.abs(sztdata['slab2dip'].values-sztdata['D2'].values)
    sztdata1 = sztdata[sztdata.dipdiff1 < sztdata.dipdiff2]
    sztdata2 = sztdata[sztdata.dipdiff1 >= sztdata.dipdiff2]
    sztdata1['closedip'] = sztdata1['D1'].values*1.0
    sztdata2['closedip'] = sztdata2['D2'].values*1.0
    sztdata = pd.concat([sztdata1,sztdata2])
    meanevendip = sztdata['closedip'].mean()

    # save figure
    figtitle = '%s/%s_slab2_szt_%s.png' % (savedir,slab,date)
    fig.savefig(figtitle)
    plt.close()

    if savedir != 'Output/%s_slab2_%s'%(slab,date):
        return szt, lsz, depthlist, norm_pdfdT, lendata
    else:
        return szt, lsz

def orgEQs(interface,eventlist,maxdepdiff, seismo_thick, slab, maxdep):

    eventlist = eventlist[eventlist.etype == 'EQ']
    print ('length of original eventlist',len(eventlist))
    # get IDs for interface events
    eventlist.loc[eventlist.slab2dep < 0, 'slab2dep'] *= -1.0
    if len(interface) > 0:
        interface = interface[interface.mdep < seismo_thick]
        sztIDs = list(interface['ID'].values)
    else:
        sztIDs = []
    print ('length of interface events',len(interface))
    
    # initialize bins for catalog and quality
    interA = eventlist[eventlist['ID'].isin(sztIDs)]
    upperA = pd.DataFrame()
    intraA = pd.DataFrame()
    interB = pd.DataFrame()
    upperB = pd.DataFrame()
    intraB = pd.DataFrame()
    
    # remove interface events from other catalogs
    eventlist = eventlist[~eventlist['ID'].isin(sztIDs)]
    
    # differentiate between events with/without cmt locations
    mtevents = eventlist[np.isfinite(eventlist.mdep)]
    otevents = eventlist[np.isnan(eventlist.mdep)]
    donotuse = pd.DataFrame()
    lonmin = otevents['lon'].min()
    lonmax = otevents['lon'].max()
    latmin = otevents['lat'].min()
    latmax = otevents['lat'].max()
    datainfo = 'test.txt'
    getfixed = True
    print ('length mt and ot events',len(mtevents), len(otevents))
    
    otevents = removePoints(donotuse, otevents, lonmin, lonmax, latmin, latmax, False, datainfo, getfixed, slab)
    
    print ('length mt and ot events',len(mtevents), len(otevents))
    
    # sort events with CMTs based on depth
    sztmt = mtevents[mtevents.mdep <= seismo_thick]
    depmt = mtevents[mtevents.mdep > seismo_thick]
    
    # deeper than seismogenic zone depths, split non-interface events to above/below
    if len(depmt) > 0:
        upperA1 = depmt[depmt.mdep < depmt.slab2dep-maxdepdiff]
        intraA1 = depmt[depmt.mdep >= depmt.slab2dep-maxdepdiff]
        upperA = pd.concat([upperA, upperA1])
        intraA = pd.concat([intraA, intraA1])

    # within seismogenic zone depths, split non-interface events to above/below
    if len(sztmt) > 0:
        upperA2 = sztmt[sztmt.mdep < sztmt.slab2dep]
        intraA2 = sztmt[sztmt.mdep >= sztmt.slab2dep]
        upperA = pd.concat([upperA, upperA2])
        intraA = pd.concat([intraA, intraA2])
        
    # sort events without CMTs based on depth
    sztot = otevents[otevents.depth <= seismo_thick]
    depot = otevents[otevents.depth > seismo_thick]

    # sort events without mts, but deeper than sz by slab depth
    if len(depot) > 0:
        upperB1 = depot[depot.depth < depot.slab2dep-maxdepdiff]
        intraB1 = depot[depot.depth >= depot.slab2dep-maxdepdiff]
        upperB = pd.concat([upperB, upperB1])
        intraB = pd.concat([intraB, intraB1])

    # get non MT interface events using buffer
    if len(sztot) > 0:
        interB1 = sztot[(sztot.depth <= sztot.slab2dep + maxdepdiff) & (sztot.depth >= sztot.slab2dep-maxdepdiff)]
        sztot = sztot[(sztot.depth > sztot.slab2dep + maxdepdiff) | (sztot.depth < sztot.slab2dep-maxdepdiff)]
        interB = pd.concat([interB, interB1])
        
        # split remaining events above/below slab
        if len(sztot) > 0:
            upperB2 = sztot[sztot.depth < sztot.slab2dep]
            intraB2 = sztot[sztot.depth >= sztot.slab2dep]
            upperB = pd.concat([upperB, upperB2])
            intraB = pd.concat([intraB, intraB2])

    interA['qual'] = 'A'
    upperA['qual'] = 'A'
    intraA['qual'] = 'A'

    interB['qual'] = 'B'
    upperB['qual'] = 'B'
    intraB['qual'] = 'B'

    if len(upperA) > 0:
        moveintraA = upperA[upperA.mdep > maxdep]
        upperA = upperA[upperA.mdep <= maxdep]
    if len(upperB) > 0:
        moveintraB = upperB[upperB.depth > maxdep]
        upperB = upperB[upperB.depth <= maxdep]
    
    inter = pd.concat([interA, interB])
    upper = pd.concat([upperA, upperB])
    intra = pd.concat([intraA, moveintraA, intraB, moveintraB])

    inter = inter[['lat','lon','depth','unc','etype','ID','mag','time','S1','D1','R1','S2','D2','R2','src','slab2str','slab2dip','slab2rak','mrr','mtt','mpp','mrt','mrp','mtp','mrrS','mttS','mppS','mrtS','mrpS','mtpS','kagan','depdiff','slab2dep','mlon','mlat','mdep','qual']]
    upper = upper[['lat','lon','depth','unc','etype','ID','mag','time','S1','D1','R1','S2','D2','R2','src','slab2str','slab2dip','slab2rak','mrr','mtt','mpp','mrt','mrp','mtp','mrrS','mttS','mppS','mrtS','mrpS','mtpS','kagan','depdiff','slab2dep','mlon','mlat','mdep','qual']]
    intra = intra[['lat','lon','depth','unc','etype','ID','mag','time','S1','D1','R1','S2','D2','R2','src','slab2str','slab2dip','slab2rak','mrr','mtt','mpp','mrt','mrp','mtp','mrrS','mttS','mppS','mrtS','mrpS','mtpS','kagan','depdiff','slab2dep','mlon','mlat','mdep','qual']]

    print ('length of sorted event catalogs',len(inter), len(upper), len(intra))
    
    return inter, upper, intra

def sortoverlap(slab1, slab2, date1, date2, savedir):

    ''' assumes slab1 is the underriding plate '''
    # read in currently stored files
    halupper = pd.read_csv('%s/%s_slab2_upper_%s.csv'%(savedir,slab1,date1))
    halintra = pd.read_csv('%s/%s_slab2_intra_%s.csv'%(savedir,slab1,date1))
    halinter = pd.read_csv('%s/%s_slab2_inter_%s.csv'%(savedir,slab1,date1))
    sulupper = pd.read_csv('%s/%s_slab2_upper_%s.csv'%(savedir,slab2,date2))
    sulintra = pd.read_csv('%s/%s_slab2_intra_%s.csv'%(savedir,slab2,date2))
    sulinter = pd.read_csv('%s/%s_slab2_inter_%s.csv'%(savedir,slab2,date2))
    
    # get lists of IDs
    sulupperids = list(sulupper['ID'].values)
    sulinterids = list(sulinter['ID'].values)
    halinterids = list(halinter['ID'].values)
    sulintraids = list(sulintra['ID'].values)
    
    # remove overriding data from upper file for lower data
    halupper = halupper[~halupper['ID'].isin(sulupperids)]
    halupper = halupper[~halupper['ID'].isin(sulinterids)]
    halupper = halupper[~halupper['ID'].isin(sulintraids)]
    
    # remove underriding interface info from upper intra info
    sulintra = sulintra[~sulintra['ID'].isin(halinterids)]
    
    # get intra events in lower that exist in intra list for upper, split A and B
    sulintraids = list(sulintra['ID'].values)
    halintraYover = halintra[halintra['ID'].isin(sulintraids)]
    halintraYoverA = halintraYover[halintraYover.qual == 'A']
    halintraYoverB = halintraYover[halintraYover.qual == 'B']
    
    # save non-overlapping events for later, split A and B
    halintraNover = halintra[~halintra['ID'].isin(sulintraids)]
    halintraNoverA = halintraNover[halintraNover.qual == 'A']
    halintraNoverB = halintraNover[halintraNover.qual == 'B']
    
    # get lower intra events - based on depth > lower slab depth
    halintraunderA = halintraYoverA[halintraYoverA.mdep > halintraYoverA.slab2dep]
    halintraunderB = halintraYoverB[halintraYoverB.depth > halintraYoverB.slab2dep]
    halintraunderAids = list(halintraunderA['ID'].values)
    halintraunderBids = list(halintraunderB['ID'].values)
    
    # remove events defined as intra for lower plate from overlapping segment of upper plate
    sulintra = sulintra[~sulintra['ID'].isin(halintraunderAids)]
    sulintra = sulintra[~sulintra['ID'].isin(halintraunderBids)]
    
    # make new under intraplate file with overlapping and non-overlapping events
    halintra = pd.concat([halintraNoverA, halintraunderA, halintraNoverB, halintraunderB])
    
    # rewrite upper and intra file for lower slab, intra file for upper slab
    halintra.to_csv('%s/%s_slab2_intra_%s.csv'%(savedir,slab1,date1),header=True,index=False,na_rep=np.nan)
    halupper.to_csv('%s/%s_slab2_upper_%s.csv'%(savedir,slab1,date1),header=True,index=False,na_rep=np.nan)
    sulintra.to_csv('%s/%s_slab2_intra_%s.csv'%(savedir,slab2,date2),header=True,index=False,na_rep=np.nan)


def getSZthickness1(data,folder,slab):

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    pd.options.mode.chained_assignment = None

    if slab == 'camz':
        maxdep=50
    elif slab == 'pan':
        maxdep=60
    else:
        maxdep = 65

    maxdep = 65
    if slab == 'sol':
        data = data[data.lon > 148]

    if slab == 'hel' or slab == 'car' or slab == 'mak':
        mindata = 40
    else:
        mindata = 20
        
    maxdepdiff = 20

    #testfile = 'gausstest_full_10000_48_50_2_2_0.5_0.5'
    #data = pd.read_csv('%s.csv'%testfile)
    alldata = data[np.isfinite(data.S1)]
    alldata = alldata[(alldata.mag < 7.5)&(alldata.depdiff < maxdepdiff)&(alldata.depdiff > -1*maxdepdiff)]
    #alldata = alldata[(alldata.R1>30) & (alldata.R2>30)
    #                       & (alldata.R1<150) & (alldata.R2<150)]

    alldata = alldata[alldata.mdep<=maxdep]
    
    maxadd = 120-35
    for i in range(maxadd):
        dat = alldata[alldata.kagan < 35+i]
        if len(dat) > mindata:
            maxkagan = 35+i
            break

    try:
        alldata = alldata[alldata.kagan < maxkagan]
    except:
        print ('not enough events within surface filter')
        maxkagan = 100
        alldata = alldata[alldata.kagan < maxkagan]
    
    alldata = alldata[['lat','lon','depth','unc','etype','ID','mag','time','S1','D1','R1','S2','D2','R2','src','slab2str','slab2dip','slab2rak','mrr','mtt','mpp','mrt','mrp','mtp','mrrS','mttS','mppS','mrtS','mrpS','mtpS','kagan','depdiff','slab2dep','mlon','mlat','mdep']]

    try:
        (slab,slab2k,date) = folder.split('_')
        alldata.to_csv('Output/%s/%s_slab2_szt_%s.csv' % (folder, slab, date),header=True,index=False,na_rep=np.nan)
    except:
        #print ('requested directory does not exist:','Output/%s/%s_slab2_szt_%s.csv' % (folder, slab, date))
        #print ('saving file in current directory as:', 'SZdist_%s.png' % (folder))
        alldata.to_csv('SZdist_%s.csv' % (folder),header=True,index=False,na_rep=np.nan)

    lendata = len(alldata)

    # make depth array for gaussian fitting
    depths = alldata.mdep.as_matrix()
    N = len(depths)
    depths.shape = (N,1)

    # initialize empty list for smooth depths and array of 0s for histogram rep.
    dc = [0 for i in range(maxdep)]
    depth2 = []
    i = 0

    # loop through range of depths
    for i in range(maxdep):

        # make window for smoothing
        mind = i-2
        if mind < 0:
            mind = 0
        maxd = i+2
        if maxd > maxdep:
            maxd = maxdep

        # loop through depth array
        for d in depths:
        
            # if within window, append depth to smooth depths, incrament hist. value
            if d >= mind and d <= maxd:
                dc[i] = dc[i]+1
                depth2.append(i)

    # normalize histogram value by dividing by amount of depths
    i2 = float(len(depth2))
    dc2 = [x / i2 for x in dc]

    # make smooth depth array for gaussian fitting
    N2 = len(depth2)
    depths2 = np.array(depth2)
    depths2.shape = (N2,1)

    # get mean and stdv.
    m0 = np.mean(depths2)
    sd0 = np.std(depths2)

    # calculate normal distribution from mean and stdv.
    depthlist = list(range(maxdep))
    norm_pdf = matplotlib.mlab.normpdf(depthlist,m0,sd0)

    # find bimodal distribution using 3 Gaussians
    clf3 = mixture.GaussianMixture(n_components=3, covariance_type='full')
    clf3.fit(depths2)

    clf2 = mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf2.fit(depths2)
    # acquire weights, means, and covariances from contributing distributions

    m1, m2, m3 = clf3.means_
    w1, w2, w3 = clf3.weights_
    c1, c2, c3 = clf3.covariances_

    md1, md2 = clf2.means_
    wd1, wd2= clf2.weights_
    cd1, cd2 = clf2.covariances_

    # calculate standard deviations
    sd1=np.sqrt(c1)[0]
    sd2=np.sqrt(c2)[0]
    sd3=np.sqrt(c3)[0]
    
    sdd1=np.sqrt(cd1)[0]
    sdd2=np.sqrt(cd2)[0]

    norm_pdf1 = matplotlib.mlab.normpdf(depthlist,m1,sd1)*w1
    norm_pdf2 = matplotlib.mlab.normpdf(depthlist,m2,sd2)*w2
    norm_pdf3 = matplotlib.mlab.normpdf(depthlist,m3,sd3)*w3

    norm_pdfT = norm_pdf1 + norm_pdf2 + norm_pdf3

    norm_pdfd1 = matplotlib.mlab.normpdf(depthlist,md1,sdd1)*wd1
    norm_pdfd2 = matplotlib.mlab.normpdf(depthlist,md2,sdd2)*wd2

    norm_pdfdT = norm_pdfd1 + norm_pdfd2
    
    # calculate rms for all distributions
    rms1a = math.sqrt(mean_squared_error(dc2, norm_pdf))
    rms3a = math.sqrt(mean_squared_error(dc2, norm_pdfT))
    rms2a = math.sqrt(mean_squared_error(dc2, norm_pdfdT))

    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999

    if lendata < -50:
        lobound = 0.1
        hibound = 0.9
    else:
        lobound = 0.05
        hibound = 0.95

    for i in range(len(norm_pdfT)):
        sumtest += norm_pdfT[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break

    if sumindex == -999:
        sumindex = len(depthlist) - 1
    szT_triple = depthlist[sumindex]
    tap_triple = depthlist[tapindex]
    low_triple = depthlist[lowindex]

    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999
    for i in range(len(norm_pdfdT)):
        sumtest += norm_pdfdT[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break

    if sumindex == -999:
        sumindex = len(depthlist) - 1
    szT_double = depthlist[sumindex]
    tap_double = depthlist[tapindex]
    low_double = depthlist[lowindex]

    sumtest = 0.0
    tapindex = -999
    sumindex = -999
    lowindex = -999
    for i in range(len(norm_pdf)):
        sumtest += norm_pdf[i]
        if sumtest < lobound:
            lowindex = i
        if sumtest >= 0.65 and tapindex < 0:
            tapindex = i
        if sumtest >= hibound:
            sumindex = i
            break
    if sumindex == -999:
        sumindex = len(depthlist) - 1

    szT_single = depthlist[sumindex]
    tap_single = depthlist[tapindex]
    low_single = depthlist[lowindex]

    # plot to show fit compared to data
    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(111)
    ax1.plot([depthlist[0],depthlist[0]],[0,dc2[0]],linewidth=10,c='k',label='Data')
    for i in range(1,len(dc2)):
        ax1.plot([depthlist[i],depthlist[i]],[0,dc2[i]],linewidth=10,c='k')
    #ax1.plot(depthlist,norm_pdf1,linewidth=2,c='skyblue',label='m=%.2f, s=%.2f, w=%.2f'%(m1,sd1,w1))
    #ax1.plot(depthlist,norm_pdf2,linewidth=2,c='skyblue',label='m=%.2f, s=%.2f, w=%.2f'%(m2,sd2,w2))
    #ax1.plot(depthlist,norm_pdf3,linewidth=2,c='skyblue',label='m=%.2f, s=%.2f, w=%.2f'%(m3,sd3,w3))
    ax1.plot(depthlist,norm_pdfd1,linewidth=2,c='springgreen',label='m=%.2f, s=%.2f, w=%.2f'%(m1,sd1,w1))
    ax1.plot(depthlist,norm_pdfd2,linewidth=2,c='springgreen',label='m=%.2f, s=%.2f, w=%.2f'%(m2,sd2,w2))
    ax1.plot(depthlist,norm_pdf,label='RMS1: %.4f'%rms1a,linewidth=2,c='y')
    ax1.plot(depthlist,norm_pdfdT,label='RMS2: %.4f'%rms2a,linewidth=2,c='g')
    #ax1.plot(depthlist,norm_pdfT,label='RMS3: %.4f'%rms3a,linewidth=2,c='b')
    if lendata < -50:
        ax1.plot([szT_single,szT_single],[0,np.max(dc2)],'r:',label='90th_single (%.4f)'%szT_single,linewidth=2)
        ax1.plot([szT_double,szT_double],[0,np.max(dc2)],'r-.',label='90th_double (%.4f)'%szT_double,linewidth=2)
        ax1.plot([low_single,low_single],[0,np.max(dc2)],'c:',label='10th_single (%.4f)'%low_single,linewidth=2)
        ax1.plot([low_double,low_double],[0,np.max(dc2)],'c-.',label='10th_double (%.4f)'%low_double,linewidth=2)
    else:
        ax1.plot([szT_single,szT_single],[0,np.max(dc2)],'r:',label='95th_single (%.4f)'%szT_single,linewidth=2)
        ax1.plot([szT_double,szT_double],[0,np.max(dc2)],'r-.',label='95th_double (%.4f)'%szT_double,linewidth=2)
        ax1.plot([low_single,low_single],[0,np.max(dc2)],'c:',label='5th_single (%.4f)'%low_single,linewidth=2)
        ax1.plot([low_double,low_double],[0,np.max(dc2)],'c-.',label='5th_double (%.4f)'%low_double,linewidth=2)
    #ax1.plot([szT_triple,szT_triple],[0,np.max(dc2)],'r--',label='szT_triple (%.4f)'%szT_triple,linewidth=2)
    if rms1a < rms3a:
        if rms1a < rms2a:
            ax1.plot(szT_single,np.max(dc2),'go',linewidth=5)
            szt, tap = szT_single, tap_single
        else:
            ax1.plot(szT_double,np.max(dc2),'go',linewidth=5)
            szt, tap = szT_double, tap_double
        #ax1.plot([tap_single,tap_single],[0,np.max(dc2)],'p-.',label='taper depth (%.4f)'%tap_single,linewidth=2)
    else:
        if rms3a < rms2a:
            ax1.plot(szT_triple,np.max(dc2),'go',linewidth=5)
            szt, tap = szT_triple, tap_triple
        else:
            ax1.plot(szT_double,np.max(dc2),'go',linewidth=5)
            szt, tap = szT_double, tap_double
        #ax1.plot([tap_triple,tap_triple],[0,np.max(dc2)],'p--',label='taper depth (%.4f)'%tap_triple,linewidth=2)

    ax1.set_xlim([0,65])
    ax1.legend(loc='best')
    ax1.grid()
    ax1.set_xlabel('Depths')
    ax1.set_ylabel('P')
    ax1.set_title('Depth distribution (%s) %i EQs (surfacefilt = %i km, kaganfilt = %i deg)'%(slab,len(alldata),maxdepdiff,maxkagan))
    
    try:
        figtitle = 'Output/%s/%s_slab2_szt_%s.png' % (folder,slab,date)
        #figtitle = '%s_notsmoothed.png'%testfile
        fig.savefig(figtitle)
        plt.close()
    except:
        #print ('requested directory does not exist:','Output/%s/%s_slab2_szt_%s.png' % (folder,slab,date))
        #print ('saving file in current directory as:', 'SZdist_%s.png' % (folder))
        figtitle = 'SZdist_%s.png' % (folder)
        fig.savefig(figtitle)
        plt.close()

    return szt, tap

def nodesift(nodes, spacing):

    uppernodes = pd.DataFrame()
    lowernodes = pd.DataFrame()
    deepnodes = nodes[nodes.depth > 50]
    deepnodes = deepnodes.reset_index(drop=True)
    for i, row in deepnodes.iterrows():
        
        lon, lat = row['lon'], row['lat']
        try:
            depth, strike = row['depth'], row['sstr']
        except:
            depth = row['depth']

        nearnodes = nodes[(nodes.lat < lat+2*spacing) & \
                                (nodes.lat > lat-2*spacing) & \
                                (nodes.lon < lon+2*spacing) & \
                                (nodes.lon > lon-2*spacing)]
        

        nearnodes['zdist'] = nearnodes['depth'].values - depth
        
        nearnodes = nearnodes[nearnodes.zdist < -25]
        
        try:
            nearnodes['cosdistance'], cosangles = npcosine(lon, lat, nearnodes['lon'].values, nearnodes['lat'].values)
            cosangles -= 180
            cosangles[cosangles<0]+=360
            nearnodes['outboard'] = np.logical_not(npoutboard(strike, cosangles))

            nearnodes1 = nearnodes[(nearnodes.outboard == True)]
        except:
            nearnodes1 = nearnodes[nearnodes.depth >0]

        nearnodes2 = nearnodes[(nearnodes.lat < lat+spacing) & \
                                (nearnodes.lat > lat-spacing) & \
                                (nearnodes.lon < lon+spacing) & \
                                (nearnodes.lon > lon-spacing)]

        #print (i, len(deepnodes))
        #print ('lon, lat, depth, strike, len(nearnodes1), len(nearnodes2)',lon, lat, depth, strike, len(nearnodes1), len(nearnodes2))
        nodedf = nodes[(nodes.lat == lat)&(nodes.lon == lon)&(nodes.depth == depth)]
        if len(nearnodes1)>0 or len(nearnodes2)>0:
            lowernodes = pd.concat([lowernodes, nodedf])
        else:
            uppernodes = pd.concat([uppernodes, nodedf])

    uppernodes = pd.concat([nodes[nodes.depth <= 50], uppernodes])
    return lowernodes, uppernodes


def projLW(lmax, wmax, plon, plat, pstr, data, namelist):

    lons = data[:,0]
    lats = data[:,1]
    dist, phis = npcosine(plon, plat, lons, lats)
    alph = np.abs(phis-pstr)
    dws = np.abs(dist*np.sin(np.radians(alph)))
    dls = np.abs(dist*np.cos(np.radians(alph)))

    newdat = pd.DataFrame({'lon':lons, 'lat':lats, 'dl':dls, 'dw':dws, 'alph':alph, 'dist':dist, 'phis':phis, 'strike':pstr})
    
    n = 0
    for name in namelist:
        if name != 'lon' and name != 'lat':
            newdat[name] = data[:,n]
        n += 1
        
    lonlatstr = '%s_%s_%s'%(int(plon*100),int(plat*100),int(pstr*100))
    #newdat.to_csv('projtest/%s.csv'%lonlatstr, header=True, index=False, na_rep=np.nan)
    
    newdat = newdat[(newdat.dw < wmax*111.19) & (newdat.dl < lmax*111.19) & (newdat.dl > 0.1*111.19) & (alph < 90)]
    #if len(newdat)>2:
    #    newdat.to_csv('projtest/%s_2.csv'%lonlatstr, header=True, index=False, na_rep=np.nan)
    
    namelist.append('dist')
    newdat = newdat[namelist]
    return newdat
    
def extendEdges(nodes,spacing,slab):
    
    nodedat = np.zeros((len(nodes),2))
    nodedat[:,0] = nodes['pslon'].values
    nodedat[:,1] = nodes['pslat'].values
    
    projlons = []
    projlats = []
    projdeps = []
    oglons = []
    oglats = []
    ogstrs = []
    ogstrs180 = []
    for index,row in nodes.iterrows():
    
        lon, lat, strike = row['pslon'], row['pslat'], row['sstr']
        strike90 = strike+90
        strike180 = strike+180
        strike270 = strike+270

        if strike90>360:
            strike90 -= 360
        if strike180>360:
            strike180 -= 360
        if strike270>360:
            strike270 -= 360

        projdat0 = projLW(20, spacing, lon, lat, strike, nodedat, spacing)
        projdat1 = projLW(20, spacing, lon, lat, strike180, nodedat, spacing)
        projdat2 = projLW(20, spacing, lon, lat, strike90, nodedat, spacing)
        projdat3 = projLW(20, spacing, lon, lat, strike270, nodedat, spacing)
        
        if (len(projdat0)<1 or len(projdat1)<1) and (len(projdat2)>0 and len(projdat3)>0) and (slab != 'sol' or lon > 160) and (slab != 'sum' or lat > 25):
            pplon = row['lon']
            pplat = row['lat']
            ppdep = row['depth']
            
            if len(projdat0<1):
                lonadd1, latadd1 = heading(pplon, pplat, 20, strike)
                lonadd2, latadd2 = heading(pplon, pplat, 40, strike)
                lonadd3, latadd3 = heading(pplon, pplat, 60, strike)
            else:
                lonadd1, latadd1 = heading(pplon, pplat, 20, strike180)
                lonadd2, latadd2 = heading(pplon, pplat, 40, strike180)
                lonadd3, latadd3 = heading(pplon, pplat, 60, strike180)
            
            if len(projdat0<1):
                lonadd1, latadd1 = heading(pplon, pplat, 10, strike)
                lonadd2, latadd2 = heading(pplon, pplat, 30, strike)
                lonadd3, latadd3 = heading(pplon, pplat, 50, strike)
            else:
                lonadd1, latadd1 = heading(pplon, pplat, 10, strike180)
                lonadd2, latadd2 = heading(pplon, pplat, 30, strike180)
                lonadd3, latadd3 = heading(pplon, pplat, 50, strike180)
            
            projlons.append(lonadd1)
            projlats.append(latadd1)
            projlons.append(lonadd2)
            projlats.append(latadd2)
            projlons.append(lonadd3)
            projlats.append(latadd3)
            
            projdeps.append(ppdep)
            projdeps.append(ppdep)
            projdeps.append(ppdep)

            oglons.append(pplon)
            oglons.append(pplon)
            oglons.append(pplon)

            oglats.append(pplat)
            oglats.append(pplat)
            oglats.append(pplat)
            
            ogstrs.append(strike)
            ogstrs.append(strike)
            ogstrs.append(strike)
            
            ogstrs180.append(strike180)
            ogstrs180.append(strike180)
            ogstrs180.append(strike180)
            

    tempcoords = pd.DataFrame({'lon':projlons, 'lat':projlats, 'depth':projdeps, 'stdv':100, 'oglon':oglons, 'oglat':oglats, 'strike':ogstrs, 'strike180':ogstrs180})
    #tempcoords = pd.DataFrame({'lon':oglons, 'lat':oglats, 'depth':projdeps, 'stdv':100, 'oglon':oglons, 'oglat':oglats})

    nodes = pd.concat([tempcoords, nodes])

    return nodes, tempcoords
        
def addToDataInfo(data, node, infostr, datainfo, dfORindiv):

    print ('node, infostr, data',node, infostr, data)
    

    if dfORindiv == 'df':
        data = data.reset_index(drop=True)
        f = open(datainfo, 'a')
        for index,row in data.iterrows():
            dataID = row['ID']
            f.write('%i,%i,%s \n' % (dataID,node,infostr))
        f.close()
    elif dfORindiv == 'indiv':
        f = open(datainfo, 'a')
        f.write('%i,%i,%s \n' % (data,node,infostr))
        f.close()
    else:
        print ('need to specify dfORindiv, not written',node, data)


    
def removeSZnodes(nodes, fracS, percshift, SZT):

    names = list(nodes)
    nodes['maxshift'] = nodes['thickness'].values * fracS
    nodes['pershift'] = nodes['smag']/nodes['maxshift'].values
    nodesSZ = nodes[(nodes.smag != 0) & (nodes.pershift<=percshift)]
    nodes = nodes[(nodes.smag == 0) | (nodes.pershift>percshift)]

    nodes = nodes[(nodes.smag != 0) | (nodes.psdepth < SZT - 0.05*SZT)]
    nodes = nodes[names]

    return nodes, nodesSZ


def nodeMisfit(nodes, results, clip):

    results = results[~np.isnan(results[:,3])]
    lons = nodes['lon'].values
    lats = nodes['lat'].values
    deps = nodes['depth'].values
    
    try:
        stdv = nodes['stdv'].values
    except:
        stdv = nodes['unc'].values

    whts = 1.0/(stdv*stdv)
    whts = whts/np.sum(whts)

    xi = np.zeros((len(lons),2))
    xi[:,0] = lons
    xi[:,1] = lats
    mask = maskdatag(clip, xi)
    deps = deps*mask

    whts = whts[np.isfinite(deps)]
    lons = lons[np.isfinite(deps)]
    lats = lats[np.isfinite(deps)]
    deps = deps[np.isfinite(deps)]
    
    datafit = 0
    distances = []
    rlons = []
    rlats = []
    for i in range(len(lons)):
        lon = lons[i]
        lat = lats[i]
        dep = deps[i]
        wht = whts[i]

        locr = results[(results[:,0] < lon+1) & (results[:,0] > lon-1) & \
                            (results[:,1] < lat+1) & (results[:,1] > lat-1)]

        if len(locr) < 1:
            locr = results[(results[:,0] < lon+2) & (results[:,0] > lon-2) & \
                            (results[:,1] < lat+2) & (results[:,1] > lat-2)]

        r1 = 6371 - dep
        r2 = 6371 - locr[:,3]

        p1 = np.radians(lon)
        p2 = np.radians(locr[:,0])

        t1 = np.radians(np.abs(lat - 90.0))
        t2 = np.radians(np.abs(locr[:,1] - 90.0))
        
        dist = r1*r1 + r2*r2 - 2*r1*r2*(math.sin(t1)*np.sin(t2)*np.cos(p1-p2) + math.cos(t1)*np.cos(t2))
        
        try:
            mindist = math.sqrt(np.min(dist))
            datafit += mindist*wht
            distances.append(mindist)
            rlon = locr[:,0][np.argmin(dist)]
            rlat = locr[:,1][np.argmin(dist)]
            rlons.append(rlon)
            rlats.append(rlat)
            if mindist > 100:
                #print ('mindist big',mindist)
                r2test = r2[np.argmin(dist)]
                #print ('r1, r2test',r1,r2test)
                #print ('lon, rlon',lon, rlon)
                #print ('lat, rlat',lat, rlat)
                t2test = t2[np.argmin(dist)]
                p2test = p2[np.argmin(dist)]
                #print ('math.degrees(t1), math.degrees(t2test)',math.degrees(t1), math.degrees(t2test))
                #print ('math.degrees(p1), math.degrees(p2test)',math.degrees(p1), math.degrees(p2test))

        except:
            #print ('np.min(dist)', lon, lat, rlon, rlat,dist)
            distances.append(0)
            rlons.append(lon)
            rlats.append(lat)

    datafit /= math.sqrt(len(lons))
    misfit = pd.DataFrame({'lon':lons, 'lat':lats, 'depth':deps, 'misfit':distances, 'rlon':rlons, 'rlat':rlats, 'diff':distances})
    return misfit, datafit

def depthMisfit(nodes, results, clip):

    results = results[~np.isnan(results[:,3])]
    lons = nodes['lon'].values
    lats = nodes['lat'].values
    deps = nodes['depth'].values
    try:
        stdv = nodes['stdv'].values
    except:
        stdv = nodes['unc'].values

    xi = np.zeros((len(lons),2))
    xi[:,0] = lons
    xi[:,1] = lats
    mask = maskdatag(clip, xi)
    deps = deps*mask
    
    lons = lons[np.isfinite(deps)]
    lats = lats[np.isfinite(deps)]
    stdv = stdv[np.isfinite(deps)]
    deps = deps[np.isfinite(deps)]
    
    xy = np.zeros((len(lons),2))
    xy[:,0] = lons
    xy[:,1] = lats
    
    rdep = griddata(results[:, 0:2], results[:, 3], xy, method='nearest')

    datafit = mean_squared_error(deps, rdep, 1/(stdv*stdv))/1
    depdiff = deps-rdep

    misfit = pd.DataFrame({'lon':lons, 'lat':lats, 'depth':deps, 'rdeps':rdep, 'diff':depdiff})
    return misfit, datafit

def plotLcurve(misfitdf, figsave):
    
    misfitdf['dfit'] = misfitdf['dfit'].values
    
    # take 'first derivative' of L curve
    f = misfitdf['filt'].values
    d = misfitdf['dfit'].values
    dd = np.ones(len(f))
    for i in range(len(f)-1):
        dd[i] = (d[i+1] - d[i])/d[i+1]
    dd[-1] = dd[-2]
    misfitdf['dd'] = dd

    
    # take 'second derivative' of L curve
    d2d = np.ones(len(f))
    for i in range(len(f)-1):
        d2d[i] = (dd[i+1] - dd[i])/(f[i+1] - f[i])
    d2d[-1] = d2d[-2]
    misfitdf['d2d'] = d2d
    misfitdf['actualfilter'] = 1.0/misfitdf['filt'].values
    
    # identify 'curviest part' of lcurve - minima of second derivative
    mind2d = misfitdf['d2d'].max()
    mindf = misfitdf[misfitdf.d2d == mind2d]
    bestfilt = mindf['filt'].values[0]
    bestdfit = mindf['dfit'].values[0]
    mindd = mindf['dd'].values[0]

    # plot lcurve
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(121)
    ax1.plot(misfitdf.filt,misfitdf.dfit,'.k')
    ax1.plot(bestfilt,bestdfit,'.r')
    ax1.set_ylabel('RMS')
    ax1.set_xlabel('1/filter')
    ax1.grid()
    ymin = misfitdf['dfit'].min()-(misfitdf['dfit'].min()/10)
    ax1.set_ylim([ymin,misfitdf['dfit'].max()])
    ax1.set_xlim([0,misfitdf['filt'].max()])
    title = 'L-Curve, best fit: 1/%.4f = %.4f' % (bestfilt,1/bestfilt)
    ax1.set_title(title)

    # plot lcurve slopes (first derivative)
    ax3 = fig.add_subplot(222)
    ax3.plot(misfitdf.filt,misfitdf.dd,'.k')
    ax3.plot(bestfilt,mindd,'.r')
    ax3.set_ylabel('Slope of L-Curve (dd/df)')
    ax3.grid()
    ax3.set_ylim([misfitdf['dd'].min(),misfitdf['dd'].max()])
    ax3.set_xlim([0,misfitdf['filt'].max()])
    title = 'L-Curve minimization: 1/%.4f = %.4f' % (bestfilt,1/bestfilt)
    ax3.set_title(title)

    # plot second derivative of lcurve
    ax4 = fig.add_subplot(224)
    ax4.plot(misfitdf.filt,misfitdf.d2d,'.k')
    ax4.plot(bestfilt,mind2d,'.r')
    ax4.set_ylabel('Slope of L-Curve Slopes (d2d/df2)')
    ax4.set_xlabel('1/filter')
    ax4.grid()
    ax4.set_ylim([misfitdf['d2d'].min(),misfitdf['d2d'].max()])
    ax4.set_xlim([0,misfitdf['filt'].max()])

    fig.savefig(figsave)
    plt.close()

    return 1.0/bestfilt, misfitdf



def histofdiffs(mindist, knot_no, rbfs, filt, kdeg, slab, fullfolder, date):

    absdist = np.abs(mindist['diff'].values)
    absdistmax = int(np.max(absdist))
    diffs = mindist['diff'].values
    
    # initialize empty list for smooth depths and array of 0s for histogram rep.
    dc = [0 for i in range(absdistmax)]
    diffs2 = []
    i = 0

    # loop through range of depths
    for i in range(absdistmax):

        # make window for smoothing
        mind = i-2
        if mind < 0:
            mind = 0
        maxd = i+2
        if maxd > absdistmax:
            maxd = absdistmax

        # loop through depth array
        for d in diffs:
        
            # if within window, append depth to smooth depths, incrament hist. value
            if d >= mind and d <= maxd:
                dc[i] = dc[i]+1
                diffs2.append(i)

    # normalize histogram value by dividing by amount of depths
    i2 = float(len(diffs2))
    dc2 = [x / i2 for x in dc]

    # make smooth diff array for gaussian fitting
    N2 = len(diffs2)
    diffs2 = np.array(diffs2)
    diffs2.shape = (N2,1)

    # get mean and stdv.
    m0 = np.mean(diffs2)
    sd0 = np.std(diffs2)

    # calculate normal distribution from mean and stdv.
    difflist = list(range(absdistmax))


    fig = plt.figure(figsize=(15, 10))
    ax1 = fig.add_subplot(111)
    ax1.plot([difflist[0],difflist[0]],[0,dc2[0]],linewidth=8,c='k',label='Data')
    for i in range(1,len(dc2)):
        ax1.plot([difflist[i],difflist[i]],[0,dc2[i]],linewidth=10,c='k')

    foldertitle = '%s/%s_slab2_diffhists_%s'%(fullfolder, slab, date)
    ax1.grid()
    ax1.set_xlabel('diffs')
    ax1.set_ylabel('P')
    ax1.set_title('diff distribution %s %s %s %s' % (knot_no, rbfs, filt, kdeg))
    figtitle = '%s/diffhist_%s_%s_%s_%s.png'% (foldertitle, knot_no, rbfs, filt, kdeg)
    try:
        fig.savefig(figtitle)
    except:
        os.system('mkdir %s' % (foldertitle))
        fig.savefig(figtitle)
    plt.close()

def makesudoguide(inFile):

    eventlistALL = pd.read_table('%s' % inFile, sep=',', dtype={
            'lon': np.float64, 'lat': np.float64,'depth': np.float64,
            'unc': np.float64, 'etype': str, 'ID': np.int, 'mag': np.float64,
            'S1': np.float64, 'D1': np.float64, 'R1': np.float64,
            'S2': np.float64, 'D2': np.float64, 'R2': np.float64,
            'src': str, 'time': str})

    minlon = eventlistALL['lon'].min()
    maxlon = eventlistAll['lon'].max()
    minlat = eventlistALL['lat'].min()
    maxlat = eventlistAll['lat'].max()

    meandep = eventlistAll['depth'].mean()

    xpts = np.arange(np.floor(minlon), np.ceil(maxlon), 0.2)
    ypts = np.arange(np.floor(minlat), np.ceil(maxlat), 0.2)
    
    xpts, ypts = np.meshgrid(xpts, ypts)

    zpts = np.ones(xpts.shape)

def SDRtoMT(data, strikename, dipname, rakename, mrrn, mttn, mppn, mrtn, mrpn, mtpn):

    # define degree-radian conversions
    d2r = math.pi/180.0
    r2d = 180.0/math.pi
    
    # get strike dip and rake according to df column names
    str = data[strikename].values
    dip = data[dipname].values
    rak = data[rakename].values
    
    # get exponent
    magpow = data['mag'].values * 1.5 + 16.1
    mom = np.power(np.ones(len(data))*10, magpow)

    # get tensor components
    mrr=mom*np.sin(2*dip*d2r)*np.sin(rak*d2r)
    mtt=-mom*((np.sin(dip*d2r)*np.cos(rak*d2r)*np.sin(2*str*d2r))+(np.sin(2*dip*d2r)*np.sin(rak*d2r)*(np.sin(str*d2r)*np.sin(str*d2r))))
    mpp=mom*((np.sin(dip*d2r)*np.cos(rak*d2r)*np.sin(2*str*d2r))-(np.sin(2*dip*d2r)*np.sin(rak*d2r)*(np.cos(str*d2r)*np.cos(str*d2r))))
    mrt=-mom*((np.cos(dip*d2r)*np.cos(rak*d2r)*np.cos(str*d2r))+(np.cos(2*dip*d2r)*np.sin(rak*d2r)*np.sin(str*d2r)))
    mrp=mom*((np.cos(dip*d2r)*np.cos(rak*d2r)*np.sin(str*d2r))-(np.cos(2*dip*d2r)*np.sin(rak*d2r)*np.cos(str*d2r)))
    mtp=-mom*((np.sin(dip*d2r)*np.cos(rak*d2r)*np.cos(2*str*d2r))+(0.5*np.sin(2*dip*d2r)*np.sin(rak*d2r)*np.sin(2*str*d2r)))
    
    # add components to dataframe and return
    data[mrrn] = mrr
    data[mttn] = mtt
    data[mppn] = mpp
    data[mrtn] = mrt
    data[mrpn] = mrp
    data[mtpn] = mtp

    return data

# output format: lon(0) lat(1) raw_dep(2) smooth_dep(3) strike(4) dip(5) raw_unc(6) shift_unc(7) raw-smooth(8) thickness(9)
def getKagan(inFile, used_data, seismo_thick, output, slab, folder):

    # identify interface events from list of filtered data
    idata = used_data[(used_data.depth <= seismo_thick)&(np.isfinite(used_data.S1))]
    idata = idata.reset_index(drop=True)
    
    # identify all other events in dataset
    tdata = pd.read_csv(inFile)
    tdata.loc[tdata.lon < 0, 'lon']+=360
    ndata = removePoints(idata, tdata, tdata['lon'].min(), tdata['lon'].max(),
                            tdata['lat'].min(), tdata['lat'].max(), False,
                            'removingpoints.csv', False, slab)
    
    ndata = ndata[np.isfinite(ndata.S1)]

    ndata.to_csv('ndata.csv',header=True,index=False,na_rep=np.nan)
    idata.to_csv('idata.csv',header=True,index=False,na_rep=np.nan)
    
    
    # get strike and dip of slab at all interface locations (lon, lat)
    iarr = np.zeros((len(idata),2))
    iarr[:,0] = idata['lon'].values
    iarr[:,1] = idata['lat'].values
    idata['slab2str'] = griddata(output[:, 0:2], output[:, 4], iarr[:, 0:2], method='nearest')
    idata['slab2dip'] = griddata(output[:, 0:2], output[:, 5], iarr[:, 0:2], method='nearest')
    idata['slab2rak'] = 90.0

    # get strike and dip of slab at all NON interface locations (lon, lat)
    narr = np.zeros((len(ndata),2))
    narr[:,0] = ndata['lon'].values
    narr[:,1] = ndata['lat'].values
    ndata['slab2str'] = griddata(output[:, 0:2], output[:, 4], narr[:, 0:2], method='nearest')
    ndata['slab2dip'] = griddata(output[:, 0:2], output[:, 5], narr[:, 0:2], method='nearest')
    ndata['slab2rak'] = 90.0

    # get moment tensors from strike, dip, and rake (data)
    idMT = SDRtoMT(idata, 'S1', 'D1', 'R1', 'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp')
    ndMT = SDRtoMT(ndata, 'S1', 'D1', 'R1', 'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp')
    
    # get moment tensors from strike, dip, and rake (slab2)
    isMT = SDRtoMT(idata, 'slab2str', 'slab2dip', 'slab2rak', 'mrrS', 'mttS', 'mppS', 'mrtS', 'mrpS', 'mtpS')
    nsMT = SDRtoMT(ndata, 'slab2str', 'slab2dip', 'slab2rak', 'mrrS', 'mttS', 'mppS', 'mrtS', 'mrpS', 'mtpS')
    
    idMT = idMT[np.isfinite(idMT.mrr)]
    idMT = idMT[np.isfinite(idMT.mrrS)]
    
    ndMT = ndMT[np.isfinite(ndMT.mrr)]
    ndMT = ndMT[np.isfinite(ndMT.mrrS)]
    
    # initialize kagan angle arrays
    interface_kagans = []
    allother_kagans = []
    
    # loop through interface dataset
    for index,row in idMT.iterrows():
    
        # make moment tensor from interface EVENTS
        mrrD = row['mrr']
        mttD = row['mtt']
        mppD = row['mpp']
        mrtD = row['mrt']
        mrpD = row['mrp']
        mtpD = row['mtp']
        vm1 = [mrrD, mttD, mppD, mrtD, mrpD, mtpD]
    
        # make moment tensor from local SLAB
        mrrS = row['mrrS']
        mttS = row['mttS']
        mppS = row['mppS']
        mrtS = row['mrtS']
        mrpS = row['mrpS']
        mtpS = row['mtpS']
        vm2 = [mrrS, mttS, mppS, mrtS, mrpS, mtpS]
        
        # calculate kagan angle between event and slab
        kagan = calc_theta(vm1,vm2)
        interface_kagans.append(kagan)

    # loop through interface dataset
    for index,row in ndMT.iterrows():
    
        # make moment tensor from interface EVENTS
        mrrD = row['mrr']
        mttD = row['mtt']
        mppD = row['mpp']
        mrtD = row['mrt']
        mrpD = row['mrp']
        mtpD = row['mtp']
        vm1 = [mrrD, mttD, mppD, mrtD, mrpD, mtpD]
    
        # make moment tensor from local SLAB
        mrrS = row['mrrS']
        mttS = row['mttS']
        mppS = row['mppS']
        mrtS = row['mrtS']
        mrpS = row['mrpS']
        mtpS = row['mtpS']
        vm2 = [mrrS, mttS, mppS, mrtS, mrpS, mtpS]
        
        # calculate kagan angle between event and slab
        kagan = calc_theta(vm1,vm2)
        allother_kagans.append(kagan)

    idMT['kagan'] = interface_kagans
    ndMT['kagan'] = allother_kagans

    # writing interface/non interface kagan info to file
    idMT.to_csv('%s_interfaceevents.csv'%folder,header=True,index=False,na_rep=np.nan)
    ndMT.to_csv('%s_noninterfaceevents.csv'%folder,header=True,index=False,na_rep=np.nan)
    
    # ensuring that both histograms have the same amount of bins
    ihist, ibins = np.histogram(interface_kagans, bins='auto')
    autobins = len(ibins)
    ihist, ibins = np.histogram(interface_kagans, bins=autobins)
    nhist, nbins = np.histogram(allother_kagans, bins=autobins)

    # calculating mean and stdv of each dataset
    imean = idMT['kagan'].mean()
    istd = idMT['kagan'].std()
    nmean = ndMT['kagan'].mean()
    nstd = ndMT['kagan'].std()

    # plotting histograms and statistics for this slab
    fig = plt.figure(figsize=(10, 20))

    ax1 = fig.add_subplot(211)
    ax1.hist(interface_kagans, bins=autobins, color='skyblue')
    ax1.plot([imean, imean], [0, np.max(ihist)], '-r', linewidth=4, label='mean: %.2f'%imean)
    ax1.plot([imean-istd, imean+istd], [np.max(ihist), np.max(ihist)], '-b', linewidth=4)
    ax1.plot([imean-istd, imean-istd], [np.max(ihist), np.max(ihist)-2], '-b', linewidth=4, label='stdv: %.2f'%istd)
    ax1.plot([imean+istd, imean+istd], [np.max(ihist), np.max(ihist)-2], '-b', linewidth=4)
    ax1.set_ylabel('no. Events')
    ax1.set_xlabel('Kagan Angle')
    ax1.set_xlim([0, 120])
    ax1.grid()
    title = 'Interface Events vs Slab'
    ax1.set_title(title)
    ax1.legend(loc='best')

    ax2 = fig.add_subplot(212)
    ax2.hist(allother_kagans, bins=autobins, color='skyblue')
    ax2.plot([nmean, nmean], [0, np.max(nhist)], '-r', linewidth=4, label='mean: %.2f'%nmean)
    ax2.plot([nmean-nstd, nmean+nstd], [np.max(nhist), np.max(nhist)], '-b', linewidth=4)
    ax2.plot([nmean-nstd, nmean-nstd], [np.max(nhist), np.max(nhist)-2], '-b', linewidth=4, label='stdv: %.2f'%nstd)
    ax2.plot([nmean+nstd, nmean+nstd], [np.max(nhist), np.max(nhist)-2], '-b', linewidth=4)
    ax2.set_ylabel('no. Events')
    ax2.set_xlabel('Kagan Angle')
    ax2.set_xlim([0, 120])
    ax2.grid()
    title = 'Non-Interface Events vs Slab'
    ax2.set_title(title)
    ax2.legend(loc='best')

    figtitle = '%s_kaganinfo.png'%folder
    fig.savefig(figtitle)
    plt.close()


def set_mt(vm):
    TM      = zeros([3,3],dtype='double')
    TM[0,0] = vm[0]
    TM[0,1] = vm[3]
    TM[0,2] = vm[4]
    TM[1,0] = vm[3]
    TM[1,1] = vm[1]
    TM[1,2] = vm[5]
    TM[2,0] = vm[4]
    TM[2,1] = vm[5]
    TM[2,2] = vm[2]
    return TM

def calc_eigenvec(TM):
    V,S    = eigh(TM)
    inds   = argsort(V)
    S      = S[:,inds]
    S[:,2] = cross(S[:,0],S[:,1])
    return S

def ang_from_R1R2(R1,R2):
    return arccos((trace(dot(R1,R2.transpose()))-1.)/2.)

def calc_theta(vm1,vm2):    
    V1 = calc_eigenvec(set_mt(vm1))
    V2 = calc_eigenvec(set_mt(vm2))
    th = ang_from_R1R2(V1,V2)
    for j in range(3):
        k       = (j+1)%3
        V3      = deepcopy(V2)
        V3[:,j] = -V3[:,j]
        V3[:,k] = -V3[:,k]
        x       = ang_from_R1R2(V1,V3)
        if x < th:
            th = x
    return th*180./pi


def addGuidePoints(tmp_res, slab):

    # Control points for HEL slab guide
    #if slab == 'hel':
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.242,30.750,40,100,62.788,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.332,30.780,40,100,95.744,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.306,30.770,40,100,166.586,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.348,30.780,40,100,224.819,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.551,30.830,40,100,280.623,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,37.982,30.940,40,100,354.285,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,38.414,31.050,40,100,389.988,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,38.801,31.150,40,100,420.602,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,39.371,31.290,40,100,464.208,7.51,126.1,1,65,180,False,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([100,100,39.667,31.370,40,100,479.355,7.51,126.1,1,65,180,False,20])

    # Control points for PNG slab guide
    #if slab == 'png':
    #    tmp_res.loc[len(tmp_res)+1] = ([-4.421,138.383,40,150.700,7.51,126.1,1,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([-4.230,137.393,40,163.500,7.51,126.1,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([-4.227,139.043,40,98.690,7.51,126.1,3,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([-3.950,136.913,40,123.700,7.51,126.1,4,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([-4.085,136.761,40,194.700,7.51,126.1,5,65,180,20])

    # Control Points for south south PHI slab guide
    #if slab == 'phi':
    #    tmp_res.loc[len(tmp_res)+1] = ([7.51,126.1,40,250,7.51,126.1,1,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([7.01,126.2,40,250,7.51,126.1,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([6.51,126.4,40,250,7.51,126.1,3,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([6.01,126.5,40,250,7.51,126.1,4,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([5.51,126.8,40,250,7.51,126.1,5,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([5.01,127,40,250,7.51,126.1,6,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([4.52,127.51,40,250,7.51,126.1,7,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([4.01,127.72,40,250,7.51,126.1,8,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([3.52,127.82,40,250,7.51,126.1,9,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([3.01,128.01,40,250,7.51,126.1,10,65,180,20])
    
    # Control Points for COT slab guide
    #if slab == 'cot':
    #    tmp_res.loc[len(tmp_res)+1] = ([5.01,126.01,40,100,5,127,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([4.02,126.02,40,100,4,128,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([3.01,126.03,40,100,3,128,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([4.52,126.04,40,100,5,127,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([3.51,126.05,40,100,4,128,2,65,180,20])
    
    # Control Points for SUL slab guide
    #if slab == 'sul':
    #    tmp_res.loc[len(tmp_res)+1] = ([0.01,123.01,40,150,0,123,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.02,122.51,40,150,0,122.5,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.03,122.02,40,150,0,122,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.04,121.53,40,150,0,121.5,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.05,121.04,40,150,0,121,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.04,123.21,40,150,0,123,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.03,122.75,40,150,0,122.5,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.02,122.22,40,150,0,122,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.01,121.76,40,150,0,121.5,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([0.04,121.23,40,150,0,121,2,65,180,20])
    
    # Control Points for HAL slab guide
    #if slab == 'hal':
    #    tmp_res.loc[len(tmp_res)+1] = ([7.84,125.55,40,169,7.84,125.55,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([8.84,125.75,40,169,7.84,125.55,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([9.964,124.035,40,644.5,9.964,124.035,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([8.964,123.935,40,644.5,9.964,124.035,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([8.44,124.66,40,392.8,8.44,124.66,2,65,180,20])
    
    # Control Points for north PHI slab guide
    #if slab == 'phi':
    #    tmp_res.loc[len(tmp_res)+1] = ([17.5,122.0,40,100.0,17.5,122.0,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([17.0,122.0,40,100.0,17.0,122.0,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([16.7,121.5,40,150.0,16.7,121.5,2,65,180,20])
    
    # Control Points for north RYU slab guide
    #if slab == 'ryu': #manually assigned control points

    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,33.196,130.438,20,20,312.811,33.196,130.438,1,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,33.348,130.881,20,20,279.325,33.196,130.438,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,30.738,129.406,20,20,323.973,33.196,130.438,4,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,31.957,130.124,20,20,290.487,33.196,130.438,5,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,29.167,127.247,20,20,398.473,33.196,130.438,6,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,27.286,124.337,20,20,402.107,33.196,130.438,7,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,27.462,124.745,20,20,413.269,33.196,130.438,8,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,33.840,131.574,20,20,301.649,33.196,130.438,9,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,35.007,134.758,20,20,212.353,33.196,130.438,10,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,35.507,134.958,20,20,212.353,33.196,130.438,11,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,35.741,133.551,20,20,134.220,33.196,130.438,12,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,35.341,133.051,20,20,134.220,33.196,130.438,13,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,35.041,132.551,20,20,134.220,33.196,130.438,14,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,29.041,127.051,20,20,400.220,33.196,130.438,16,65,180,20])

    #if slab == 'ryu':
    #    # adding johnny wu's model as control points in north
    #    for index, row in ryutodata.iterrows():
    #        ryulon, ryulat, ryudep = row['lon'], row['lat'], row['depth']
    #        tmp_res.loc[len(tmp_res)+1] = ([radius1,radius2,ryulat,ryulon,20,20,ryudep,33.196,130.438,index,65,180,20])
            
    # Control Points for west SUM slab guide
    #if slab == 'sum':
    #    tmp_res.loc[len(tmp_res)+1] = ([26.128,93.193,40,40.500,17.5,122.0,2,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([26.046,94.248,40,42.600,17.5,122.0,3,65,180,20])
    #    tmp_res.loc[len(tmp_res)+1] = ([24.886,94.002,40,73.800,17.5,122.0,4,65,180,20])
    return tmp_res


def getReferenceKagan(slab1data, eventlist, origorcentl, origorcentd):

    # ensure all longitude and CMT longitudes are in 0-360
    slab1data.loc[slab1data.lon < 0, 'lon'] += 360
    eventlist.loc[eventlist.lon < 0, 'lon'] += 360
    eventlist.loc[eventlist.mlon < 0, 'mlon'] += 360
    
    # identify data with and without moment tensor information
    odata = eventlist[np.isnan(eventlist.S1)]
    idata = eventlist[np.isfinite(eventlist.mlon)]
    
    # move slab geometry dataframe to array
    output = np.zeros((len(slab1data),5))
    output[:,0] = slab1data['lon'].values*1.0
    output[:,1] = slab1data['lat'].values*1.0
    output[:,2] = slab1data['strike'].values*1.0
    output[:,3] = slab1data['dip'].values*1.0
    output[:,4] = slab1data['depth'].values*1.0
    
    # get strike and dip of slab at all interface locations (lon, lat)
    iarr = np.zeros((len(idata),2))
    
    if origorcentl == 'o':
        # assign slab2 fault plane geometry at PDE origin location
        iarr[:,0] = idata['lon'].values
        iarr[:,1] = idata['lat'].values
    elif origorcentl == 'c':
        # assign slab2 fault plane geometry at CMT location
        iarr[:,0] = idata['mlon'].values
        iarr[:,1] = idata['mlat'].values
    else:
        print ('origorcentl must be o or c, enter -l o or -l c at input')
        print ('exiting ... ')
        sys.exit()

    # get slab2 fault plane geometry at each event location
    idata['slab2str'] = griddata(output[:, 0:2], output[:, 2], iarr[:, 0:2], method='nearest')
    idata['slab2dip'] = griddata(output[:, 0:2], output[:, 3], iarr[:, 0:2], method='nearest')
    idata['slab2dep'] = griddata(output[:, 0:2], output[:, 4], iarr[:, 0:2], method='nearest')
    idata['slab2rak'] = 90.0

    # get moment tensors from strike, dip, and rake (data)
    idMT = SDRtoMT(idata, 'S1', 'D1', 'R1', 'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp')
    
    # get moment tensors from strike, dip, and rake (slab2)
    isMT = SDRtoMT(idata, 'slab2str', 'slab2dip', 'slab2rak', 'mrrS', 'mttS', 'mppS', 'mrtS', 'mrpS', 'mtpS')
    
    idMT = idMT[np.isfinite(idMT.mrr)]
    idMT = idMT[np.isfinite(idMT.mrrS)]
    
    # initialize kagan angle arrays
    interface_kagans = []
    
    # loop through interface dataset
    for index,row in idMT.iterrows():
    
        # make moment tensor from interface EVENTS
        mrrD = row['mrr']
        mttD = row['mtt']
        mppD = row['mpp']
        mrtD = row['mrt']
        mrpD = row['mrp']
        mtpD = row['mtp']
        vm1 = [mrrD, mttD, mppD, mrtD, mrpD, mtpD]
    
        # make moment tensor from local SLAB
        mrrS = row['mrrS']
        mttS = row['mttS']
        mppS = row['mppS']
        mrtS = row['mrtS']
        mrpS = row['mrpS']
        mtpS = row['mtpS']
        vm2 = [mrrS, mttS, mppS, mrtS, mrpS, mtpS]
        
        # calculate kagan angle between event and slab
        kagan = calc_theta(vm1,vm2)
        interface_kagans.append(kagan)

    # add kagan angle and depth difference values to input dataset
    idMT['kagan'] = interface_kagans
    
    if origorcentd == 'o':
        idMT['depdiff'] = idMT['depth'].values - (-1*idMT['slab2dep'].values)
    elif origorcentd == 'c':
        idMT['depdiff'] = idMT['mdep'].values - (-1*idMT['slab2dep'].values)
    else:
        print ('origorcentd must be o or c, enter -d o or -d c at input')
        print ('exiting ... ')
        sys.exit()

    oarr = np.zeros((len(odata),2))
    oarr[:,0] = odata['lon'].values
    oarr[:,1] = odata['lat'].values
    odata['slab2str'] = griddata(output[:, 0:2], output[:, 2], oarr[:, 0:2], method='nearest')
    odata['slab2dip'] = griddata(output[:, 0:2], output[:, 3], oarr[:, 0:2], method='nearest')
    odata['slab2dep'] = griddata(output[:, 0:2], output[:, 4], oarr[:, 0:2], method='nearest')
    odata['slab2rak'] = 90.0

    newlist = pd.concat([idMT, odata])
    
    return newlist


def getReferenceKagan1(slab1data, eventlist):

    slab1data.loc[slab1data.lon < 0, 'lon'] += 360
    eventlist.loc[eventlist.lon < 0, 'lon'] += 360
    eventlist.loc[eventlist.mlon < 0, 'mlon'] += 360
    
    ogcolumns = eventlist.columns
    
    odata = eventlist[np.isnan(eventlist.S1)]
    idata = eventlist[np.isfinite(eventlist.mlon)]
    
    output = np.zeros((len(slab1data),5))
    output[:,0] = slab1data['lon'].values*1.0
    output[:,1] = slab1data['lat'].values*1.0
    output[:,2] = slab1data['strike'].values*1.0
    output[:,3] = slab1data['dip'].values*1.0
    output[:,4] = slab1data['depth'].values*1.0
    
    # get strike and dip of slab at all interface locations (lon, lat)
    iarr = np.zeros((len(idata),2))
    iarr[:,0] = idata['mlon'].values
    iarr[:,1] = idata['mlat'].values
    #print (iarr)
    #print (output)
    idata['slab2str'] = griddata(output[:, 0:2], output[:, 2], iarr[:, 0:2], method='nearest')
    idata['slab2dip'] = griddata(output[:, 0:2], output[:, 3], iarr[:, 0:2], method='nearest')
    idata['slab2dep'] = griddata(output[:, 0:2], output[:, 4], iarr[:, 0:2], method='nearest')
    idata['slab2rak'] = 90.0

    # get moment tensors from strike, dip, and rake (data)
    idMT = SDRtoMT(idata, 'S1', 'D1', 'R1', 'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp')
    
    # get moment tensors from strike, dip, and rake (slab2)
    isMT = SDRtoMT(idata, 'slab2str', 'slab2dip', 'slab2rak', 'mrrS', 'mttS', 'mppS', 'mrtS', 'mrpS', 'mtpS')
    
    idMT = idMT[np.isfinite(idMT.mrr)]
    idMT = idMT[np.isfinite(idMT.mrrS)]
    
    # initialize kagan angle arrays
    interface_kagans = []
    
    # loop through interface dataset
    for index,row in idMT.iterrows():
    
        # make moment tensor from interface EVENTS
        mrrD = row['mrr']
        mttD = row['mtt']
        mppD = row['mpp']
        mrtD = row['mrt']
        mrpD = row['mrp']
        mtpD = row['mtp']
        vm1 = [mrrD, mttD, mppD, mrtD, mrpD, mtpD]
    
        # make moment tensor from local SLAB
        mrrS = row['mrrS']
        mttS = row['mttS']
        mppS = row['mppS']
        mrtS = row['mrtS']
        mrpS = row['mrpS']
        mtpS = row['mtpS']
        vm2 = [mrrS, mttS, mppS, mrtS, mrpS, mtpS]
        
        # calculate kagan angle between event and slab
        kagan = calc_theta(vm1,vm2)
        interface_kagans.append(kagan)

    idMT['kagan'] = interface_kagans
    idMT['depdiff'] = idMT['mdep'].values - (-1*idMT['slab2dep'].values)
    newlist = pd.concat([idMT, odata])
    return newlist


def makenewclip(clip,shift_out,griddf,TRdata,node,slab):

    clip = clip.reset_index(drop=True)

    newlons, newlats = [], []
    for index,row in clip.iterrows():
        lon, lat = row['lon'], row['lat']
        
        #print (index,lon,lat)
        if len(TRdata)>0 and (slab != 'sol' or lon > 150):
            loc_tr = TRdata[(TRdata.lon > lon-3) & (TRdata.lon < lon+3) & (TRdata.lat > lat-3) & (TRdata.lat < lat+3)]
            if len(loc_tr)>0:
                #loc_tr['dist'] = gps2dist_azimuth(lat, lon, loc_tr['lat'], loc_tr['lon'])[0]/1000.0
                loc_tr['dist'], tempangles = npcosine(lon, lat, loc_tr['lon'].values, loc_tr['lat'].values)
                mindist = loc_tr['dist'].min()
                loc_tr = loc_tr[loc_tr.dist == mindist]
                lonT = loc_tr['lon'].values[0]
                latT = loc_tr['lat'].values[0]
                azT = loc_tr['az'].values[0]
                thisdist, thisang, latB, lonB = cosine(lonT, latT, lon, lat)
                out = isoutboard(azT, thisang)
                if out:
                    #print ('trench side of mask', lon, lat)
                    newlons.append(lon)
                    newlats.append(lat)
                    continue
        locnodes = shift_out[(shift_out.lon > lon - node*4) & \
                                (shift_out.lon < lon + node*4) & \
                                (shift_out.lat > lat - node*4) & \
                                (shift_out.lat < lat + node*4)]
        maxdep = locnodes['depth'].max()
        locgrd = griddf[(griddf.depth < maxdep + 5) & \
                        (griddf.depth > maxdep - 1) & \
                        (griddf.lon > lon - 0.5) & \
                        (griddf.lon < lon + 0.5) & \
                        (griddf.lat > lat - 0.5) & \
                        (griddf.lat < lat + 0.5)]
        if len(locgrd) < 1:
            #print ('no grid at this depth within 111 km',lon,lat,maxdep)
            newlons.append(lon)
            newlats.append(lat)
            continue
            
        locgrd['dist'] = gps2dist_azimuth(lat, lon, locgrd['lat'], locgrd['lon'])[0]
        neargrd = locgrd[locgrd.dist == locgrd['dist'].min()]
        newlon = neargrd['lon'].values[0]
        newlat = neargrd['lat'].values[0]
        
        if slab == 'man' and newlat < 18 and newlat > 17.5:
            newlons.append(lon)
            newlats.append(lat)
            continue
        
        #print ('new lon lat',lon,lat,newlon,newlat,maxdep)

        newlons.append(newlon)
        newlats.append(newlat)

    newmask = pd.DataFrame({'lon':newlons, 'lat':newlats})

    return newmask

def maketiltclip(clip,shift_out,griddf,trenches,node,slab):

    clip = clip.reset_index(drop=True)

    newlons, newlats = [], []
    for index,row in clip.iterrows():
        lon, lat = row['lon'], row['lat']
        
        #print (index,lon,lat)
        if len(trenches)>0 and (slab != 'sol' or lon > 150):
            loc_tr = trenches[(trenches.lon > lon-3) & (trenches.lon < lon+3) & (trenches.lat > lat-3) & (trenches.lat < lat+3)]
            if len(loc_tr)>0:
                #loc_tr['dist'] = gps2dist_azimuth(lat, lon, loc_tr['lat'], loc_tr['lon'])[0]/1000.0
                loc_tr['dist'], tempangles = npcosine(lon, lat, loc_tr['lon'].values, loc_tr['lat'].values)
                mindist = loc_tr['dist'].min()
                loc_tr = loc_tr[loc_tr.dist == mindist]
                lonT = loc_tr['lon'].values[0]
                latT = loc_tr['lat'].values[0]
                azT = loc_tr['az'].values[0]
                thisdist, thisang, latB, lonB = cosine(lonT, latT, lon, lat)
                out = isoutboard(azT, thisang)
                if out:
                    #print ('trench side of mask', lon, lat)
                    newlons.append(lon)
                    newlats.append(lat)
                    continue

    upclip = pd.DataFrame({'lon':newlons, 'lat':newlats})
    x00, y00 = trenches['lon'].values[-1], trenches['lat'].values[-1]
    if x00<0:
        x00 += 360

    # if trench azimuth crosses 360-0 azimuth, flip to opposing azimuth
    if trenches['az'].min() < 90 and trenches['az'].max() > 270:
        newaz = trenches['az'].values
        newaz-=180
        newaz[newaz<0]+=360
        trenches['az'] = newaz

    # find mean azimuth, ensure oriented in correct direction depending on slab
    if slab == 'izu' or slab == 'jap':
        meanstk = 170.0
    elif slab == 'man':
        meanstk = 180
    else:
        meanstk = trenches['az'].mean()
    if slab == 'phi':
        meanstk += 180
        if meanstk > 360:
            meanstk -= 360

    # find perpendicular azimuth to mean strike of trench
    if meanstk < 270.0:
        perpstk = meanstk + 90.0
    else:
        perpstk = meanstk - 270

    # offset reference point (behaves better closer or further from slab in different regions)
    if trenches['az'].min() < 90 and trenches['az'].max() > 270 or slab == 'phi' or slab == 'solz':
        x0, y0 = heading(x00, y00, 100, meanstk)

    elif trenches['az'].min() < 90 and trenches['az'].max() > 270 or slab == 'sul' or slab == 'sol':
        x0, y0 = heading(x00, y00, 800, meanstk)

    # most regions need to be flipped 180 degrees
    else:
        meanstkshift = meanstk+180
        if meanstkshift>360:
            meanstkshift -= 360
        if slab == 'izu' or slab == 'jap':
            if meanstkshift > 270:
                meanstkshift -= 270
            else:
                meanstkshift += 90
            x0, y0 = heading(x00, y00, 1800, meanstkshift)
        else:
            x0, y0 = heading(x00, y00, 800, meanstkshift)

    lons = griddf['lon'].values
    lats = griddf['lat'].values
    uncs = np.ones(len(griddf))
    deps = griddf['depth'].values
    tiltdata2,dSs,dPs = newrefframe(x0,y0,meanstk,lons,lats,deps,uncs,slab)
    tiltmask = tiltedmask(tiltdata2, slab, 10)
    tiltmask = tiltmask[tiltmask.lat > 0]

    tiltarr
    xP, yP, depi = sdptoxyz(supplement, x0, y0, meanstk)


def splitsurface(nodeFile,outFile,clipFile,trenches,node,filt,grid,slab, knot_no, kdeg, rbfs, folder):

    # import necessary data
    if slab == 'sol':
        TR_data = pd.read_csv(trenches)
        TR_data = TR_data[TR_data.slab == slab]
        TR_data = TR_data[TR_data.lon>149]
        trenchFile = 'library/misc/trenches_downsample.csv'
        trenches = pd.read_csv(trenchFile)
        trenches = trenches[trenches.slab == slab]
    elif slab == 'hin' or slab == 'pam':
        lonts = [71.34, 71.6]
        latts = [33.845, 33.89]
        azts = [270, 270]
        trenches = pd.DataFrame({'lon':lonts,'lat':latts,'az':270,'bound':'IN\EU','slab':slab, 'depth':0})
        TR_data = trenches.copy()
    else:
        trenchFile = trenches
        trenches = pd.read_csv(trenchFile)
        if slab == 'izu' or slab == 'kur':
            trenches = trenches[trenches.slab == 'jap']
        else:
            trenches = trenches[trenches.slab == slab]
        TR_data = trenches.copy()
    results = pd.read_csv(outFile)
    nodesOG = pd.read_csv(nodeFile)
    clip = pd.read_csv(clipFile)

    # define depth bounds based on slab2 dip in OG reference frame (r.f.)
    u45 = 45+5
    l45 = 45-5
    dip45 = results[(results.dip_shift_smooth < u45)&(results.dip_shift_smooth > l45)]
    dep45 = results['dep_shift_smooth'].mean()
    shdep = dep45 - 10
    dedep = dep45 + 10
    

    # define different parts of different slabs that may turn over (to minimize computation time)
    nodes = nodesOG[nodesOG.depth > dep45]
    if slab == 'sol':
        #results = results.iloc[::4, :]
        resultsout1 = results[(results.dep_shift_smooth > shdep)|(results.lon<146)|(results.lon>158)]
        results1 = results[(results.lon<146)&(results.lon>145.8)]
        results2 = results[(results.lon>158)&(results.lon<158.2)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lon>=146)&(results.lon<=158)]
        results = pd.concat([results,results1,results2])
        results = results[np.isfinite(results.dep_shift_smooth)]
        trenches = trenches[(trenches.lon>=146)&(trenches.lon<=158)]
        resultsout1['depth'] = resultsout1['dep_shift_smooth'].values*1.0
        #resultsout1.to_csv('resultsout1.csv',header=True, index=False, na_rep=np.nan)
    elif slab == 'izu' or slab == 'jap':
        #results = results.iloc[::2, :]
        resultsout = results[(results.lat>27)|(results.lat<15)]
        results2 = results[(results.lat<15)&(results.lat>14.8)]
        results1 = results[(results.lat<27.2)&(results.lat>27)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lat<=27)&(results.lat>=15)]
        results = pd.concat([results,results1,results2])
        results = results[np.isfinite(results.dep_shift_smooth)]
        results['depth'] = results['dep_shift_smooth'].values*1.0
        #results.to_csv('results.csv',header=True,index=False, na_rep=np.nan)
        trenches = trenches[(trenches.lat<=27)&(trenches.lat>=15)]
        resultsout['depth'] = resultsout['dep_shift_smooth'].values*1.0
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        nodes = nodes[(nodes.lat<=27)&(nodes.lat>=15)]
    elif slab == 'sum':
        #results1 = results[results.lat<20].iloc[::8, :]
        results2 = results[results.lat >= 20]
        results = pd.concat([results1, results2])
        resultsout = results[(results.lon<100)|(results.lon>122)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lon>=100)&(results.lon<=122)]
        results['depth'] = results['dep_shift_smooth'].values*1.0
        #results.to_csv('results.csv',header=True,index=False, na_rep=np.nan)
        trenches = trenches[(trenches.lon>=100)&(trenches.lon<=122)]
        resultsout['depth'] = resultsout['dep_shift_smooth'].values*1.0
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        nodes = nodes[(nodes.lon>=100) & (nodes.lon<=122)]
    elif slab == 'ker':
        #results = results.iloc[::8, :]
        resultsout = results[(results.dep_shift_smooth < dep45)|(results.lat>-30)]
        results1 = results[(results.lat > -30)&(results.lat < -30)&(results.dep_shift_smooth<dep45)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lat<=-30)]
        results = pd.concat([results,results1])
        results = results[np.isfinite(results.dep_shift_smooth)]
        resultsout['depth'] = resultsout['dep_shift_smooth'].values*1.0
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        nodes = nodesOG[nodesOG.lat<=-30]
    elif slab == 'manz':
        #results = results.iloc[::8, :]
        resultsout = results[(results.dep_shift_smooth < dep45)|(results.lat>21)]
        results1 = results[(results.lat>21)&(results.lat<22)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lat<=21)]
        results = pd.concat([results,results1])
        results = results[np.isfinite(results.dep_shift_smooth)]
        resultsout['depth'] = resultsout['dep_shift_smooth'].values*1.0
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        nodes = nodes[nodes.lat<=21]
    elif slab == 'puyz':
        #results = results.iloc[::8, :]
        resultsout = results[(results.dep_shift_smooth < dep45)|(results.lat<-47)]
        results = results[(results.dep_shift_smooth <= dep45)&(results.lat>=-47)]
        resultsout['depth'] = resultsout['dep_shift_smooth'].values*1.0
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        nodes = nodes[nodes.lat>=-47]
    else:
        results = results[results.dep_shift_smooth<dep45]

    # reset indices for looping through dataframes
    nodes = nodes.reset_index(drop=True)
    results = results.reset_index(drop=True)

    # if sol, only consider slab outboard trench (not opposing arm on W end)
    if slab == 'sol':


        nodes = nodes[(nodes.lon>=146) & (nodes.lon<=158)]
        nodesW = nodes[nodes.lon <= 153]
        nodesE = nodes[nodes.lon > 153]
        nodesW = nodesW.reset_index(drop=True)
        nodesW, nodesout = getoutboard(nodesW, trenches[trenches.lon < 152], slab)
        nodes = pd.concat([nodesE, nodesW])
        #nodes.to_csv('nodetest.csv',header=True,index=False, na_rep=np.nan)
        resultsW = results[results.lon <= 153]
        resultsE = results[results.lon > 153]
        resultsW = resultsW.reset_index(drop=True)
        resultsW, resultsout2 = getoutboard(resultsW, trenches[trenches.lon < 152], slab)
        resultsW['dep_shift_smooth'] = resultsW['depth'].values
        results = pd.concat([resultsE,resultsW])
        results = results[np.isfinite(results.dep_shift_smooth)]
        resultsout2['dep_shift_smooth'] = resultsout2['depth'].values*1.0
        #results.to_csv('sol_results.csv',header=True,index=False, na_rep=np.nan)
        resultsout = pd.concat([resultsout1, resultsout2])
        #resultsout2.to_csv('resultsout2.csv',header=True, index=False, na_rep=np.nan)
        #resultsout.to_csv('resultsout.csv',header=True, index=False, na_rep=np.nan)
        '''
        nodes = pd.read_csv('nodetest1.csv')
        results = pd.read_csv('results1.csv')
        resultsout = pd.read_csv('resultsout1.csv')
        resultsinout = results[results.lon<146]
        results = results[results.lon >= 146]
        resultsout = pd.concat([resultsinout,resultsout])
        nodes = nodes[nodes.lon>146]
        '''

    # define xyz coordinates of nodes
    lons = nodes['lon'].values
    lats = nodes['lat'].values
    deps = nodes['depth'].values
    uncs = nodes['stdv'].values
    sunc = nodes['shiftstd'].values
    thck = nodes['thickness'].values

    # add lon lat depth of original surface above dep45 to node list
    lons = np.hstack([lons, results['lon'].values])
    lons[lons<0]+=360
    lats = np.hstack([lats, results['lat'].values])
    uncs = np.hstack([uncs, results['dz1'].values])
    sunc = np.hstack([sunc, results['dz2'].values])
    thck = np.hstack([thck, results['thickness'].values])
    try:
        deps = np.hstack([deps, results['dep_shift_smooth'].values])
    except:
        deps = np.hstack([deps, results['depth'].values])

    # find coordinate of end of trench
    x00, y00 = trenches['lon'].values[-1], trenches['lat'].values[-1]
    if x00<0:
        x00 += 360

    # if trench azimuth crosses 360-0 azimuth, flip to opposing azimuth
    if trenches['az'].min() < 90 and trenches['az'].max() > 270:
        newaz = trenches['az'].values
        newaz-=180
        newaz[newaz<0]+=360
        trenches['az'] = newaz

    # find mean azimuth, ensure oriented in correct direction depending on slab
    if slab == 'izu' or slab == 'jap':
        meanstk = 170.0
    elif slab == 'hin':
        meanstk = 91
    elif slab == 'pam':
        meanstk = 250
    else:
        meanstk = trenches['az'].mean()
    if slab == 'phi':
        meanstk += 180
        if meanstk > 360:
            meanstk -= 360

    # find perpendicular azimuth to mean strike of trench
    if meanstk < 270.0:
        perpstk = meanstk + 90.0
    else:
        perpstk = meanstk - 270

    # offset reference point (behaves better closer or further from slab in different regions)
    if trenches['az'].min() < 90 and trenches['az'].max() > 270 or slab == 'phi' or slab == 'solz':
        x0, y0 = heading(x00, y00, 100, meanstk)

    elif trenches['az'].min() < 90 and trenches['az'].max() > 270 or slab == 'sul' or slab == 'sol' or slab == 'hin':
        x0, y0 = heading(x00, y00, 800, meanstk)

    # most regions need to be flipped 180 degrees
    else:
        meanstkshift = meanstk+180
        if meanstkshift>360:
            meanstkshift -= 360
        if slab == 'izu' or slab == 'jap':
            if meanstkshift > 270:
                meanstkshift -= 270
            else:
                meanstkshift += 90
            x0, y0 = heading(x00, y00, 1800, meanstkshift)
        else:
            x0, y0 = heading(x00, y00, 800, meanstkshift)

    tiltdata2,dSs,dPs = newrefframe(x0,y0,meanstk,lons,lats,deps,uncs,slab)
    if slab == 'sol':
        lonx = np.array([152.0,156.0])
        latx = np.array([-3.0,-6.0])
        depx = np.array([0,0])
        uncx = np.array([1,1])
        tiltdatax,dSsx,dPsx = newrefframe(x0,y0,meanstk,lonx,latx,depx,uncx,slab)
        #print (dSsx,dPsx)

    #tiltdata2.to_csv('%s_tiltdata.csv'%slab, header=True, index=False)

    if slab == 'sol':
        distcut = dSsx[0]
        distcut1 = dSsx[1]
    else:
        distcut = 0
        distcut1 = 0

    # make new labeled dataframe in new r.f. and make tilted clipping mask
    tiltdata = pd.DataFrame({'lon':dSs,'lat':deps,'depth':dPs,'unc':uncs})
    tiltmask = tiltedmask(tiltdata, slab, 10, distcut, distcut1)

    # plot data and mask for checking
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    con = ax1.scatter(tiltdata['lon'].values,tiltdata['lat'].values,c=tiltdata['depth'].values,s=10,edgecolors='none',cmap='plasma')
    prefilt = ax1.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')
    figtitle = 'tiltdata.png'
    #fig.savefig(figtitle)
    plt.close()

    # make array for making tilted surface
    data = np.zeros((len(tiltdata),4))
    data[:,0] = tiltdata['lon'].values
    data[:,1] = tiltdata['lat'].values
    data[:,2] = tiltdata['depth'].values
    data[:,3] = tiltdata['unc'].values

    # identify constants associated with OG r.f. surface
    sigma = 0.3
    rbfs2 = 10
    if slab == 'pam' or slab == 'hin':
        rbfs2 = 0.001
    spacing = grid
    node2 = 0.02
    spacing = 0.1

    newdat11 = gridthedata5(data, sigma, slab, spacing, rbfs2, dep45, meanstk)
    #np.savetxt('newdat11.csv', newdat11, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')

    # calculate spline of data and rbf filler control points
    filt2 = filt
    
    newdat, mindat, maskdat, nodemask, supplement, supp2, geomask, geosupp, dsddat = gridthedata6(data, newdat11, filt2, tiltmask, slab, kdeg, knot_no, node2, dep45, distcut, 'first', meanstk, np.zeros((1,3)),dep45,thck,uncs,sunc,tiltdata)
    minarr = np.zeros((len(mindat),3))
    minarr[:,0] = mindat['dist'].values
    minarr[:,1] = mindat['depth'].values
    minarr[:,2] = mindat['perp'].values
    maskdat2 = maskdat[(maskdat.p < 99999)&(maskdat.p > -99999)]

    supdf = pd.DataFrame({'dist':supplement[:,0], 'depth':supplement[:,1], 'perp':supplement[:,2]})
    dsddf = pd.DataFrame({'dist':dsddat[:,0], 'depth':dsddat[:,1], 'perp':dsddat[:,2], 'strike':dsddat[:,3], 'dip':dsddat[:,4], 'dz1':dsddat[:,5], 'dz2':dsddat[:,6],'dz3':dsddat[:,8],'thickness':dsddat[:,7]})
    merge1 = pd.merge(supdf, dsddf, left_on = ['dist','depth','perp'], right_on = ['dist','depth','perp'])
    #merge1.to_csv('%s_mergetest.csv'%slab,header=True,index=False,na_rep=np.nan)
    
    # move back to OG r.f. ((s,d,p) -> (x,y,z))
    xP, yP, depi = sdptoxyz(supplement, x0, y0, meanstk)
    xPM, yPM, depiM = sdptoxyz(minarr, x0, y0, meanstk)
    xPA, yPA, depiA = sdptoxyz(supp2, x0, y0, meanstk)
    # save to dataframe
    finaldat = pd.DataFrame({'lon':xP,'lat':yP,'depth':depi,'strike':merge1['strike'].values,'dip':merge1['dip'].values})
    if slab == 'ker':
        finaldat = finaldat[finaldat.lat <= -30]
    if slab == 'izu':
        finaldat = finaldat[(finaldat.lat <= 27)&(finaldat.lat >= 15)]
    if slab == 'sol':
        finaldat = finaldat[(finaldat.lon >= 146)&(finaldat.lon <= 158)]
    mindat['lon'] = xPM
    mindat['lat'] = yPM
    mindat = mindat[['lon','lat','depth','strike','dip','dist','perp','grd']]
    masknodes = pd.DataFrame({'lon':xPA, 'lat':yPA, 'depth':depiA, 'stdv':30, \
                'smag':0,'shiftstd':0,'avstr':np.nan,'avdip':np.nan, \
                'avrke':np.nan ,'psdepth':depiA ,'sstr':nodesOG['sstr'].mean(), \
                'sdip':90 ,'nID':np.ones(len(depiA))*-1 ,'pslon':xPA, \
                'pslat':yPA ,'bzlon':xPA ,'bzlat':yPA, \
                'centsurf':nodesOG['centsurf'].mean(), \
                'thickness':nodesOG['thickness'].mean(), \
                'alen':nodesOG['alen'].mean() ,'blen':nodesOG['blen'].mean(), \
                'clen':nodesOG['clen'].mean() ,'ogstr':nodesOG['ogstr'].mean(), \
                'ogdip':nodesOG['ogdip'].mean() ,'hstdv':1 ,'vstdv':1 })
    masknodes = masknodes.iloc[::4, :]

    if slab == 'ker':
        masknodes = masknodes[masknodes.lat <= -30]
    if slab == 'izu':
        masknodes = masknodes[(masknodes.lat <= 27)&(masknodes.lat >= 15)]
    if slab == 'sol':
        masknodes = masknodes[(masknodes.lon <= 158)&(masknodes.lon >= 146)]
    masknodes = masknodes.iloc[::4,:]
    masknodes1 = masknodes[['lon','lat','depth']]
    #masknodes1.to_csv('%s_apex.csv'%slab,header=True,index=False,na_rep=np.nan)
    #finaldat.to_csv('%s_pre-finalinterp.csv'%slab,header=True,index=False,na_rep=np.nan)

    # plot final interpolation with mask vs og data
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(311)
    con = ax1.scatter(tiltdata['lon'].values,tiltdata['lat'].values,c=tiltdata['depth'].values,s=10,edgecolors='none',cmap='plasma')
    ax1.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax1.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax1.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax1.plot(geomask['s'].values,geomask['d'].values,'g.',label='newmask')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    #ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax2 = fig.add_subplot(312)
    con = ax2.scatter(newdat[:,0],newdat[:,1],c=newdat[:,2],s=10,edgecolors='none',cmap='plasma')
    ax2.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax2.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax2.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax2.set_xlabel('Longitude (distance along strike)')
    ax2.set_ylabel('Latitude (Depth)')
    ax2.axis('equal')
    ax2.invert_yaxis()
    plt.grid()
    #ax2.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax3 = fig.add_subplot(313)
    con = ax3.scatter(supplement[:,0],supplement[:,1],c=supplement[:,2],s=10,edgecolors='none',cmap='plasma')
    ax3.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax3.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax3.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax3.set_xlabel('Longitude (distance along strike)')
    ax3.set_ylabel('Latitude (Depth)')
    ax3.axis('equal')
    ax3.invert_yaxis()
    plt.grid()
    #ax3.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    # save figure
    figtitle = '%s_tiltdata.png'%slab
    #fig.savefig(figtitle)
    plt.close()

    # clip original surface and nodes by mask in tilted reference frame

    (slab,slab2k,date) = folder.split('_')

    results = pd.read_csv('Output/%s/%s_slab2_res_%s.csv' % (folder, slab, date))
    nodes = pd.read_csv('Output/%s/%s_slab2_nod_%s.csv' % (folder, slab, date))

    newclip = pd.DataFrame({'lon':maskdat['s'].values, 'lat':maskdat['d'].values})
    nodeclip = pd.DataFrame({'lon':geomask['s'].values, 'lat':geomask['d'].values})

    nlons = nodes['lon'].values
    nlats = nodes['lat'].values
    ndeps = nodes['depth'].values
    nuncs = nodes['stdv'].values
    tiltnodes,dSs,dPs = newrefframe(x0,y0,meanstk,nlons,nlats,ndeps,nuncs,slab)

    nodes2, deepnodes2, tiltnodes = cliptilt(tiltnodes,newclip,nodes,finaldat,slab,'first')
    shift_out2 = pd.concat([nodes2, masknodes])
    nodes, deepnodes, tiltnodes = cliptilt(tiltnodes,nodeclip,nodes,finaldat,slab,'first')
    shift_out = pd.concat([nodes, masknodes])
    
    #shift_out.to_csv('%s_shiftout.csv'%slab,header=True,index=False)
    #shift_out2.to_csv('%s_shiftout2.csv'%slab,header=True,index=False)

    T = 0.0
    date = datetime.today().strftime('%m.%d.%y')
    now = datetime.now()
    time = '%s.%s' % (now.hour, now.minute)
    if slab == 'manz':
        node = node2

    npass = 1
    meanBA = 5.0

    if slab == 'hin' or slab == 'pam':
        TR_data = TR_data[TR_data.slab == 'xx']
    
    surfdata = np.zeros((len(shift_out), 4))
    surfdata[:, 0], surfdata[:, 1], surfdata[:, 2], surfdata[:, 3] = shift_out['lon'].values, shift_out['lat'].values, shift_out['depth'].values, shift_out['stdv'].values
    
    errordata = np.zeros((len(nodes), 4))
    errordata[:, 0], errordata[:, 1], errordata[:, 2], errordata[:, 3] = nodes['lon'].values, nodes['lat'].values, nodes['stdv'].values, np.ones(len(nodes))
    
    errordataB = np.zeros((len(nodes), 4))
    errordataB[:, 0], errordataB[:, 1], errordataB[:, 2], errordataB[:, 3] = nodes['lon'].values, nodes['lat'].values, nodes['shiftstd'].values, np.ones(len(nodes))
    
    thickdata = np.zeros((len(nodes),4))
    thickdata[:, 0], thickdata[:, 1], thickdata[:, 2], thickdata[:, 3] = nodes['lon'].values, nodes['lat'].values, nodes['thickness'].values, np.ones(len(nodes))

    if slab == 'jap':
        Surfgrid, xi, dl = chunksurface(surfdata, node, T, slab, grid, 'depth', time, 'test.txt', filt, pd.DataFrame(), npass, TR_data, meanBA, kdeg, knot_no, rbfs, shift_out,'fin','extra','lat',30,40,35)
        flipornot = 'flip'
    else:
        Surfgrid, xi, dl = pySurface3(surfdata, node, T, slab, grid, 'depth', time, 'test.txt', filt, pd.DataFrame(), npass, TR_data, meanBA, kdeg, knot_no, rbfs, shift_out,'fin','extra')
        flipornot = 'dontflip'

    Surfgrid_unmask = np.copy(Surfgrid)
    sigma = (filt/2.0) / node

    rlons = xi[:,0]
    rlats = xi[:,1]
    rdeps = Surfgrid.flatten()
    runcs = np.ones(len(rdeps))
    tiltresults,dSs,dPs = newrefframe(x0,y0,meanstk,rlons,rlats,rdeps,runcs,slab)
    results, deepresults, tiltresults = cliptilt(tiltresults,nodeclip,tiltresults,finaldat,slab,'first')
    geoarray = np.zeros((len(results),4))
    perparray = np.zeros((len(supp2),4))
    geoarray[:,0] = results['lon'].values
    geoarray[:,1] = results['lat'].values
    geoarray[:,2] = results['newlat'].values
    perparray[:,0] = xPA
    perparray[:,1] = yPA
    perparray[:,2] = depiA

    geoarray = np.vstack((geoarray,perparray))
    Surfgrid2 = griddata(geoarray[:, 0:2], geoarray[:,2], xi[:, 0:2], method = 'nearest')
    Surfgrid2.shape = Surfgrid.shape

    Errorgrid = makeErrorgrid(Surfgrid, xi, errordata)
    Errorgrid2 = makeErrorgrid(Surfgrid, xi, errordataB)
    thickgrid = makeErrorgrid(Surfgrid, xi, thickdata)

    if slab == 'hin' or slab == 'pamz':
        Surfgrid = np.copy(Surfgrid2)

    if slab == 'izu':
        filt2 = 1.5
        Filtgrid = specializufilt(Surfgrid,xi,filt,filt2,node)
        Errorgrid = specializufilt(Errorgrid,xi,filt,filt2,node)
        Errorgrid2 = specializufilt(Errorgrid2,xi,filt,filt2,node)
        thickgrid = specializufilt(thickgrid,xi,filt,filt2,node)
    else:
        Filtgrid = ndimage.filters.gaussian_filter(Surfgrid, sigma, mode='reflect')
        Errorgrid = ndimage.filters.gaussian_filter(Errorgrid, sigma, mode='reflect')
        Errorgrid2 = ndimage.filters.gaussian_filter(Errorgrid2, sigma, mode='reflect')
        thickgrid = ndimage.filters.gaussian_filter(thickgrid, sigma, mode='reflect')

    strgrid3, dipgrid3 = mkSDgrddata(xi, Filtgrid, flipornot)

    resdata = np.zeros((len(xi),5))
    resdata[:,0] = xi[:,0]
    resdata[:,1] = xi[:,1]

    resdata[:,2] = Filtgrid.flatten()
    resdata[:,3] = strgrid3.flatten()
    resdata[:,4] = dipgrid3.flatten()

    #np.savetxt('resdata.csv', resdata, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')
    #shift_out.to_csv('shiftout.csv',header=True,index=False,na_rep=np.nan)
    nxy = np.zeros((len(nodesOG),2))
    nxy[:,0] = nodesOG['lon'].values*1.0
    nxy[:,1] = nodesOG['lat'].values*1.0
    sxy = np.zeros((len(masknodes),2))
    sxy[:,0] = masknodes['lon'].values*1.0
    sxy[:,1] = masknodes['lat'].values*1.0
    suppstrs = griddata(nxy[:, 0:2], nodesOG['sstr'].values, sxy[:, 0:2], method = 'nearest')
    suppdips = griddata(nxy[:, 0:2], nodesOG['sdip'].values, sxy[:, 0:2], method = 'nearest')
    masknodes['sstr'] = suppstrs
    masknodes['sdip'] = suppdips

    clipnodes = pd.concat([nodesOG,masknodes])

    newres = mkContourClip(clipnodes, TR_data, node, resdata, False,slab)

    if len(TR_data)>0:
        clip = clippingmask(newres,TR_data,node,False, slab,'first')
    else:
        clip = noTrenchPolygon(newres, node, False, slab)

    mask = maskdatag(clip, xi)
    mask.shape = Surfgrid.shape

    strgrid3 = (strgrid3*mask)
    dipgrid3 = (dipgrid3*mask)
    Filtgrid = (Filtgrid*mask)
    Surfgrid = (Surfgrid*mask)
    Errorgrid = (Errorgrid*mask)
    Errorgrid2 = (Errorgrid2*mask)
    thickgrid = (thickgrid*mask)
    smooth_dif = Surfgrid.flatten()-Filtgrid.flatten()

    results = pd.DataFrame({'lon':xi[:, 0], 'lat':xi[:, 1], 'raw_dep':Surfgrid.flatten(), 'dep_shift_smooth':Filtgrid.flatten(), 'str_shift_smooth':strgrid3.flatten(), 'dip_shift_smooth':dipgrid3.flatten(), 'dz1':Errorgrid.flatten(), 'dz2':Errorgrid2.flatten(), 'dz3':smooth_dif.flatten(),'thickness':thickgrid.flatten()})

    rlons = results['lon'].values
    rlats = results['lat'].values
    rdeps = results['dep_shift_smooth'].values
    runcs = results['dz1'].values
    tiltresults,dSs,dPs = newrefframe(x0,y0,meanstk,rlons,rlats,rdeps,runcs,slab)
    
    geodata = np.zeros((len(tiltresults),3))
    geodata[:,0] = tiltresults['newlon'].values
    geodata[:,1] = tiltresults['newlat'].values
    geodata[:,2] = tiltresults['depth'].values
    
    geodata2 = np.zeros((len(tiltresults),3))
    geodata2[:,0] = tiltresults['lon'].values
    geodata2[:,1] = tiltresults['lat'].values
    geodata2[:,2] = tiltresults['newlat'].values

    results, deepresults, tiltresults = cliptilt(tiltresults,newclip,results,finaldat,slab,'first')
    if slab == 'izu':
        results = results[(results.lat > 15.2)|(results.dep_shift_smooth < 350)]

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(311)
    con = ax1.scatter(tiltresults['newlon'].values,tiltresults['newlat'].values,c=tiltresults['inorout'].values,s=10,edgecolors='none',cmap='plasma')
    ax1.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax1.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax1.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax2 = fig.add_subplot(312)
    con = ax2.scatter(tiltnodes['newlon'].values,tiltnodes['newlat'].values,c=tiltnodes['inorout'].values,s=10,edgecolors='none',cmap='plasma')
    ax2.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax2.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax2.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax2.set_xlabel('Longitude (distance along strike)')
    ax2.set_ylabel('Latitude (Depth)')
    ax2.axis('equal')
    ax2.invert_yaxis()
    plt.grid()
    ax2.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax3 = fig.add_subplot(313)
    con = ax3.scatter(tiltnodes['newlon'].values,tiltnodes['newlat'].values,c=tiltnodes['depth'].values,s=10,edgecolors='none',cmap='plasma')
    ax3.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax3.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax3.set_xlabel('Longitude (distance along strike)')
    ax3.set_ylabel('Latitude (Depth)')
    ax3.axis('equal')
    ax3.invert_yaxis()
    plt.grid()
    ax3.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    # save figure
    figtitle = '%s_tiltdata2.png'%slab
    #fig.savefig(figtitle)
    plt.close()
    
    ''' re making supplemental surface with geo surface'''

    # identify constants associated with OG r.f. surface
    sigma = 0.3
    rbfs2 = 10
    spacing = grid
    node2 = 0.02
    spacing = 0.1

    nlons = deepnodes2['lon'].values
    nlats = deepnodes2['lat'].values
    ndeps = deepnodes2['depth'].values
    nuncs = deepnodes2['stdv'].values
    tiltnodes,dSs,dPs = newrefframe(x0,y0,meanstk,nlons,nlats,ndeps,nuncs,slab)
    sminx = finaldat['lon'].min()
    smaxx = finaldat['lon'].max()
    sminy = finaldat['lat'].min()
    smaxy = finaldat['lat'].max()
    sminz = finaldat['depth'].min()
    smaxz = finaldat['depth'].max()

    results = results[np.isfinite(results.dep_shift_smooth)]
    rlons = results['lon'].values
    rlats = results['lat'].values
    rdeps = results['dep_shift_smooth'].values
    runcs = results['dz1'].values
    tiltresults,dSs,dPs = newrefframe(x0,y0,meanstk,rlons,rlats,rdeps,runcs,slab)
    
    if slab == 'ker':
        tiltresults = tiltresults[tiltresults.lat < smaxy + 5]
    if slab == 'izu':
        tiltresults = tiltresults[(tiltresults.lat < smaxy + 5)&(tiltresults.lat > sminy - 2)]
    if slab == 'sol':
        tiltresults = tiltresults[(tiltresults.lon < smaxx + 2)&(tiltresults.lon > sminx - 2)]

    #tiltresults.to_csv('%s_newsupptest.csv'%slab,header=True,index=False,na_rep=np.nan)

    #tiltresults = tiltresults.iloc[::4,:]
    newdata = pd.concat([tiltnodes,tiltresults])
    data = np.zeros((len(newdata),4))
    data[:,0] = newdata['newlon'].values*1.0
    data[:,1] = newdata['newlat'].values*1.0
    data[:,2] = newdata['depth'].values*1.0
    data[:,3] = newdata['unc'].values*1.0

    newdat11 = gridthedata5(data, sigma, slab, spacing, rbfs2, sminz, meanstk)
    #np.savetxt('newdat11.csv', newdat11, header='lon,lat,depth,strike,dip',fmt='%.2f', delimiter=',',comments='')

    # calculate spline of data and rbf filler control points
    filt2 = filt

    newdat, mindat, maskdat, nodemask, supplement, supp2, geomask, geosupp, dsddat = gridthedata6(data, newdat11, filt2, tiltmask, slab, kdeg, knot_no, node2, sminz, distcut, 'second', meanstk, geodata, sminz,thck,uncs,sunc,tiltdata)
    
    supdf = pd.DataFrame({'dist':supplement[:,0], 'depth':supplement[:,1], 'perp':supplement[:,2]})
    dsddf = pd.DataFrame({'dist':dsddat[:,0], 'depth':dsddat[:,1], 'perp':dsddat[:,2], 'strike':dsddat[:,3], 'dip':dsddat[:,4], 'dz1':dsddat[:,5], 'dz2':dsddat[:,6],'dz3':dsddat[:,8],'thickness':dsddat[:,7]})

    merge1 = pd.merge(supdf, dsddf, left_on = ['dist','depth','perp'], right_on = ['dist','depth','perp'])
    #merge1.to_csv('%s_mergetest.csv'%slab,header=True,index=False,na_rep=np.nan)

    # move back to OG r.f. ((s,d,p) -> (x,y,z))
    xP, yP, depi = sdptoxyz(supplement, x0, y0, meanstk)
    # save to dataframe
    finaldat = pd.DataFrame({'lon':xP,'lat':yP,'newlat':depi,'newlon':supplement[:,0],'depth':supplement[:,2], 'unc':10, 'strike':merge1['strike'].values,'dip':merge1['dip'].values,'dz1':merge1['dip'].values, 'dz1':merge1['dz1'].values, 'dz2':merge1['dz2'].values,'dz3':merge1['dz3'].values,'thickness':merge1['thickness'].values})
    shallowdat, finaldat, alldat = cliptilt(finaldat,newclip,finaldat,finaldat,slab,'first')
    
    if slab == 'ker':
        finaldat = finaldat[finaldat.lat <= -30]
    if slab == 'izu':
        finaldat = finaldat[(finaldat.lat <= 27)&(finaldat.lat >= 15)]
    if slab == 'sol':
        finaldat = finaldat[(finaldat.lon >= 146)&(finaldat.lon <= 158)]
    if slab == 'man':
        finaldat = finaldat[(finaldat.lat >= clip['lat'].min()) & (finaldat.lat <= clip['lat'].max())]
    
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(311)
    con = ax1.scatter(tiltresults['newlon'].values,tiltresults['newlat'].values,c=tiltresults['depth'].values,s=10,edgecolors='none',cmap='plasma')
    ax1.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax1.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax1.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax2 = fig.add_subplot(312)
    con = ax2.scatter(data[:,0],data[:,1],c=data[:,2],s=10,edgecolors='none',cmap='plasma')
    ax2.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax2.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax2.plot(nodemask['s'].values,nodemask['d'].values,'b.',label='newmask')
    ax2.set_xlabel('Longitude (distance along strike)')
    ax2.set_ylabel('Latitude (Depth)')
    ax2.axis('equal')
    ax2.invert_yaxis()
    plt.grid()
    ax2.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax3 = fig.add_subplot(313)
    con = ax3.scatter(finaldat['newlon'].values,finaldat['newlat'].values,c=finaldat['depth'].values,s=10,edgecolors='none',cmap='plasma')
    #con = ax3.scatter(results['newlon'].values,results['newlat'].values,c=results['depth'].values,s=10,edgecolors='none',cmap='plasma')
    ax3.plot(tiltmask['lon'].values,tiltmask['lat'].values,'k.',label='Mask')
    ax3.plot(maskdat['s'].values,maskdat['d'].values,'r.',label='newmask')
    ax3.set_xlabel('Longitude (distance along strike)')
    ax3.set_ylabel('Latitude (Depth)')
    ax3.axis('equal')
    ax3.invert_yaxis()
    plt.grid()
    ax3.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')
    
    # save figure
    figtitle = '%s_tiltdata3.png'%slab
    #fig.savefig(figtitle)
    plt.close()

    if slab == 'hin' or slab == 'pam':
        polyclip = makepolymask(slab,'library/misc/slab_polygons.txt')
        finaldat.loc[finaldat.lon < 0, 'lon']+=360
        polyclip.loc[polyclip.lon < 0, 'lon']+=360
        
        pts = np.zeros((len(finaldat),2))
        pts[:, 0] = finaldat['lon'].values
        pts[:, 1] = finaldat['lat'].values
        
        mask = maskdatag(polyclip, pts)
    
        finaldat['depth'] = finaldat['depth'].values*mask
        finaldat = finaldat[np.isfinite(finaldat.depth)]
        finaldat = finaldat[finaldat.lon < 76]
        finaldat = finaldat.reset_index(drop=True)
        #finaldat.to_csv('%s_finaldat.csv'%slab,header=True,index=False)

    ''' ... done with regeneration of supplement ... '''

    # Create output array
    print('%s_%s_%s_%s' % (slab, date, time, str(grid)))
    print("    Populating output array...")

    output = (np.zeros([len(results), 10]) * np.nan)

    output[:, 0] = results['lon'].values  # lon Longitude at node (not shifted)
    output[:, 1] = results['lat'].values  # lat Latitude at node
    output[:, 2] = results['raw_dep'].values  # dep_shift Post-shift surface depth before smoothing
    output[:, 3] = results['dep_shift_smooth'].values  # dep_shift_smooth Post-shift surface depth after smoothing
    output[:, 4] = results['str_shift_smooth'].values # str_shift_smooth Post-shift surface strike after smoothing (strike was not smoothed - only depth was smoothed)
    output[:, 5] = results['dip_shift_smooth'].values  # dip_shift_smooth Post-shift surface dip after smoothing
    output[:, 6] = results['dz1'].values # dz1 Interpolated, but unsmoothed uncertainty from raw data
    output[:, 7] = results['dz2'].values  #dz2 Interpolated, unsmoothed uncertainty from shift
    output[:, 8] = results['dz3'].values # dz3 error induced by smoothing (taken as the standard deviation of smoothed-unsmoothed)
    output[:, 9] = results['thickness'].values  #dz2 Interpolated, unsmoothed thickness
    output[:, 0][output[:, 0]<0]+=360

    finaldat['depth'] = finaldat['newlat'].values
    finaldat = finaldat[['lon','lat','depth','strike','dip','dz1','dz2','dz3','thickness']]

    if node > node2:
        node2 = 0.05
    xiold = np.copy(xi)
    xall = np.arange(np.floor(np.min(resdata[:,0])), np.ceil(np.max(resdata[:,0])), node2)
    yall = np.arange(np.floor(np.min(resdata[:,1])), np.ceil(np.max(resdata[:,1])), node2)
    xpts, ypts = np.meshgrid(xall, yall)
    xi = np.zeros((len(xpts.flatten()),2))
    xi[:,0] = xpts.flatten()
    xi[:,1] = ypts.flatten()
    interpdepths = griddata(xiold,Surfgrid_unmask.flatten(),xi,method='nearest')
    mask = maskdatag(clip, xi)
    Surfgrid = interpdepths*mask
    results = pd.DataFrame({'lon':xi[:, 0], 'lat':xi[:, 1], 'depth':Surfgrid, 'dz1':1})
    rlons = results['lon'].values
    rlats = results['lat'].values
    rdeps = results['depth'].values
    runcs = results['dz1'].values
    tiltresults,dSs,dPs = newrefframe(x0,y0,meanstk,rlons,rlats,rdeps,runcs,slab)
    results, deepresults, tiltresults = cliptilt(tiltresults,newclip,results,finaldat,slab,'second')
    outresults = results[np.isnan(results.inorout)]
    outresults = outresults[outresults.lon >= finaldat['lon'].min()]
    outresults = outresults[outresults.lon <= finaldat['lon'].max()]
    outresults = outresults[outresults.lat >= finaldat['lat'].min()]
    outresults = outresults[outresults.lat <= finaldat['lat'].max()]
    outresults = outresults.reset_index(drop=True)
    clfilons = []
    clfilats = []
    clfideps = []
    for index,row in outresults.iterrows():
        rlon,rlat = row['lon'],row['lat']
        finnear = finaldat[(finaldat.lon > rlon-node2)&(finaldat.lon < rlon+node2)&(finaldat.lat > rlat-node2)&(finaldat.lat < rlat+node2)]
        if len(finnear)>0:
            clfilons.append(rlon)
            clfilats.append(rlat)
            clfideps.append(finnear['depth'].values[0])
    if slab == 'izu':
        results = results[(results.lat > 15.2)|(results.depth < 350)]

    results = results[np.isfinite(results.inorout)]
    results = results[np.isfinite(results.depth)]
    #print ('results?',results)
    resnew = results[['lon','lat','depth']]
    finnew = pd.DataFrame({'lon':clfilons,'lat':clfilats,'depth':clfideps})

    newres = pd.concat([resnew,finnew])
    #newres.to_csv('%s_extramask.csv'%slab,header=True,index=False)

    if len(TR_data)>0:
        clip = clippingmask(newres,TR_data,node2,False, slab, 'second')
    else:
        clip = noTrenchPolygon(newres, node, False, slab)

    clip.loc[clip.lon < 0, 'lon']+=360

    return clip, output, finaldat, nodes, deepnodes2

def gridthedata5(data, sigma, slab, spacing, rbfs, dedep, meanstk):
    
    # get coordinates for finding extrema of dataset
    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    
    # define new grid for interpolating based on spacing and extent of OG r.f.
    gridsp = spacing * 111.19
    xi = np.arange(np.floor(np.min(x))-2, np.ceil(np.max(x))+2, gridsp)
    yi = np.arange(np.floor(np.min(y))-2, np.ceil(np.max(y))+2, gridsp)
    xpts, ypts = np.meshgrid(xi, yi)
    xyzip = np.zeros((len(xpts.flatten()),2))
    xyzip[:, 0] = xpts.flatten()
    xyzip[:, 1] = ypts.flatten()
    
    # separate upper and lower datasets (surface and nodes)
    dataup = data[data[:,1] <= dedep]
    datado = data[data[:,1] > dedep]

    # resample upper part of dataset to this grid size (too big)
    resdataup = np.zeros((len(xyzip),4))
    resdataup[:,0] = xpts.flatten()
    resdataup[:,1] = ypts.flatten()
    resdataup[:,2] = griddata(dataup[:, 0:2], dataup[:,2], resdataup[:, 0:2], method = 'nearest')
    resdataup[:,3] = griddata(dataup[:, 0:2], dataup[:,3], resdataup[:, 0:2], method = 'nearest')
    resdataup = resdataup[resdataup[:,1]<=dedep]

    # combine datasets back together and reassign x,y,z
    data = np.vstack((resdataup,datado))
    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    con = ax1.scatter(x,y,c=z,s=12,edgecolors='none',cmap='plasma')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')
    figtitle = 'resampled.png'
    #fig.savefig(figtitle)
    plt.close()
    
    smoother = rbfs
    
    # make radial basis function of dataset
    try:
        interp = Rbf(x, y, z, function='linear',smooth=smoother)
    
    # function fails if data is too regular - add small random noise if fails
    except:
        addx = np.random.rand(len(x))/1000
        x = x+addx
        y = y+addx
        z = z+addx
        interp = Rbf(x, y, z, function='linear', smooth=smoother)

    # evaluate the radial basis function at defined grid coordinates and filter
    
    zi = interp(xpts, ypts)

    zif = ndimage.filters.gaussian_filter(zi, sigma/2)

    # calculate strike and dip of grid in new reference frame
    Ygrad, Xgrad = np.gradient(zif, gridsp, gridsp, edge_order=2)
    strikegrid = np.degrees(np.arctan(Xgrad)) # do this for strike!!!!
    
    quadgrid = strikegrid.flatten()
    tau = meanstk
    for i in range(len(quadgrid)):
        lam = quadgrid[i]
        beta = tau-lam
        if beta < 90 and beta > 0:
            quadgrid[i] = beta
        elif beta < 0 and beta > -90:
            quadgrid[i] = abs(beta)
        else:
            beta = lam-tau
            if beta < 90 and beta > 0:
                quadgrid[i] = beta
            else:
                beta = abs(beta)
                if beta > 90:
                    a = 180-beta
                    quadgrid[i] = a

    quadgrid.shape = strikegrid.shape
    beta = np.radians(quadgrid)
    delta = np.arctan(np.abs(Ygrad))

    dipgrid = np.degrees(np.arctan(np.tan(delta)*np.cos(beta)))
    
    dipgrid = 90-dipgrid
    strikegrid = meanstk - strikegrid
    strikegrid[strikegrid < 0] += 360
    if slab == 'pam' or slab == 'man' or slab == 'hin':
        strikegrid -= 180
        strikegrid[strikegrid < 0] += 360
    
    # save information to dataframe
    newdat = np.zeros((len(zif.flatten()), 7))
    newdat[:, 0], newdat[:, 1], newdat[:, 2] = xpts.flatten(), ypts.flatten(), zif.flatten()
    newdat[:, 3], newdat[:, 4] = strikegrid.flatten(), dipgrid.flatten()
    newdat[:, 5], newdat[:, 6] = Xgrad.flatten(), Ygrad.flatten()

    # plot depth strike and dip of grid
    fig = plt.figure(figsize=(20, 25))
    ax1 = fig.add_subplot(311)
    con = ax1.scatter(newdat[:,0],newdat[:,1],c=newdat[:,2],s=50,edgecolors='none',cmap='plasma')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax2 = fig.add_subplot(312)
    con = ax2.scatter(newdat[:,0],newdat[:,1],c=newdat[:,3],s=50,edgecolors='none',cmap='plasma')
    ax2.set_ylabel('Latitude (Depth)')
    ax2.axis('equal')
    ax2.invert_yaxis()
    plt.grid()
    cbar = fig.colorbar(con)
    cbar.set_label('Strike')
    
    ax3 = fig.add_subplot(313)
    con = ax3.scatter(newdat[:,0],newdat[:,1],c=newdat[:,4],s=50,edgecolors='none',cmap='plasma')
    ax3.set_ylabel('Latitude (Depth)')
    ax3.axis('equal')
    ax3.invert_yaxis()
    plt.grid()
    cbar = fig.colorbar(con)
    cbar.set_label('Dip')

    # save figure and return rbf data
    figtitle = '%s_depthstrikedip.png'%slab
    #fig.savefig(figtitle)
    plt.close()

    return newdat

def newrefframe(x0,y0,meanstk,lons,lats,deps,uncs,slab):
    # shift from OG r.f. to new r.f ((x,y,z) -> (s,d,p))
    print ('x0,y0,meanstk',x0,y0,meanstk)
    dTs,thetas = npcosine(x0, y0, lons, lats)
    alphas = np.abs(meanstk-thetas)
    dSs = dTs*np.cos(np.radians(alphas))
    dPs = dTs*np.sin(np.radians(alphas))
    stk180 = meanstk - 180
    if stk180 < 0:
        stk180 += 360

    # not sure why I do this (need to investigate)
    if slab != 'phi' and slab != 'sol':
        dPs[(thetas>meanstk) | (thetas<stk180)] *= -1

    # add all components to dataframe and save
    tiltdata2 = pd.DataFrame({'newlon':dSs,'newlat':deps,'depth':dPs,'unc':uncs,'lon':lons,'lat':lats,'thetas':thetas,'alphas':alphas, 'dTs':dTs})
    return tiltdata2,dSs,dPs

def cliptilt(tiltresults,newclip,results,suppdat,slab,fors):
    xyzip = np.zeros((len(tiltresults),2))
    xyzip[:,0] = tiltresults['newlon'].values
    xyzip[:,1] = tiltresults['newlat'].values
    rmask = maskdataT(newclip,xyzip)
    tiltresults['inorout'] = rmask
    minlon, maxlon = suppdat['lon'].min(0), suppdat['lon'].max(0)
    minlat, maxlat = suppdat['lat'].min(0), suppdat['lat'].max(0)
    if slab == 'sol':
        tiltresults.loc[tiltresults.lon < minlon, 'inorout'] = 1
        tiltresults.loc[tiltresults.lon > maxlon, 'inorout'] = 1
    if slab == 'izu' or slab == 'jap' or slab == 'manz' or slab == 'ker' or slab == 'puyz':
        tiltresults.loc[tiltresults.lat < minlat, 'inorout'] = 1
        tiltresults.loc[tiltresults.lat > maxlat, 'inorout'] = 1
    results['inorout'] = tiltresults['inorout'].values
    deepresults = results[np.isnan(results.inorout)]
    if fors == 'first':
        results = results[np.isfinite(results.inorout)]
    return results, deepresults, tiltresults

def cliptiltopp(tiltresults,newclip,results,suppdat,slab):
    xyzip = np.zeros((len(tiltresults),2))
    xyzip[:,0] = tiltresults['newlon'].values
    xyzip[:,1] = tiltresults['newlat'].values
    rmask = maskdataT(newclip,xyzip)
    tiltresults['inorout'] = rmask
    minlon, maxlon = suppdat['lon'].min(0), suppdat['lon'].max(0)
    minlat, maxlat = suppdat['lat'].min(0), suppdat['lat'].max(0)
    if slab == 'sol':
        tiltresults.loc[tiltresults.lon < minlon, 'inorout'] = np.nan
        tiltresults.loc[tiltresults.lon > maxlon, 'inorout'] = np.nan
    if slab == 'izu' or slab == 'jap' or slab == 'manz' or slab == 'ker' or slab == 'puyz':
        tiltresults.loc[tiltresults.lat < minlat, 'inorout'] = np.nan
        tiltresults.loc[tiltresults.lat > maxlat, 'inorout'] = np.nan
    results['inorout'] = tiltresults['inorout'].values
    deepresults = results[np.isnan(results.inorout)]
    results = results[np.isfinite(results.inorout)]
    return results, deepresults, tiltresults

def nodecliptilt(tiltresults,newclip,results,suppdat):

    dxy = np.zeros((len(tiltresults),4))
    dxy[:,0] = tiltresults['newlon'].values
    dxy[:,1] = tiltresults['newlat'].values
    newclip = newclip[newclip.lat > 10]
    dxy[:,2] = griddata(newclip['lon'].values, newclip['lat'].values, dxy[:,1], method = 'nearest')
    for i in range(len(dxy)):
        if dxy[i,1] < dxy[i,2]:
            dxy[i,3] = 1
        else:
            #print (dxy[i,0],dxy[i,1], dxy[i,2])
            dxy[i,3] = 0

    tiltresults['inorout'] = dxy[:,3]
    #print (tiltresults[tiltresults.inorout == 0])
    minlon, maxlon = suppdat['lon'].min(0), suppdat['lon'].max(0)
    minlat, maxlat = suppdat['lat'].min(0), suppdat['lat'].max(0)
    tiltresults.loc[tiltresults.lon < minlon, 'inorout'] = 1
    tiltresults.loc[tiltresults.lon > maxlon, 'inorout'] = 1
    tiltresults.loc[tiltresults.lat < minlat, 'inorout'] = 1
    tiltresults.loc[tiltresults.lat > maxlat, 'inorout'] = 1
    #print (tiltresults[tiltresults.inorout == 0])
    results['inorout'] = tiltresults['inorout'].values
    deepresults = results[results.inorout == 0]
    results = results[results.inorout == 1]

    return results, deepresults, tiltresults
    
def gridthedata6(data, filldat1, filt, clip, slab, kdeg, knot_no, node, dedep, distcut, fors, meanstk, geodata, sminz,thck,uncs,sunc,tiltdata):
    
    # get coordinates for finding extrema of dataset
    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    
    errordata1 = np.zeros((len(tiltdata),3))
    errordata1[:,0] = tiltdata['lon'].values
    errordata1[:,1] = tiltdata['lat'].values
    errordata1[:,2] = uncs
    
    errordata2 = np.zeros((len(tiltdata),3))
    errordata2[:,0] = tiltdata['lon'].values
    errordata2[:,1] = tiltdata['lat'].values
    errordata2[:,2] = sunc
    
    thickdata = np.zeros((len(tiltdata),3))
    thickdata[:,0] = tiltdata['lon'].values
    thickdata[:,1] = tiltdata['lat'].values
    thickdata[:,2] = thck
    
    # define new grid for interpolating based on spacing and extent of OG r.f.
    gridsp = node * 111.19
    xi = np.arange(np.floor(np.min(x))-2, np.ceil(np.max(x))+2, gridsp)
    yi = np.arange(np.floor(np.min(y))-2, np.ceil(np.max(y))+2, gridsp)
    xpts, ypts = np.meshgrid(xi, yi)
    xyzip = np.zeros((len(xpts.flatten()),2))
    xyzip[:, 0] = xpts.flatten()
    xyzip[:, 1] = ypts.flatten()
    
    # separate upper and lower datasets (surface and nodes)
    dataup = data[data[:,1] <= dedep]
    datado = data[data[:,1] > dedep]

    # resample upper part of dataset to this grid size (too big)
    resdataup = np.zeros((len(xyzip),4))
    resdataup[:,0] = xpts.flatten()
    resdataup[:,1] = ypts.flatten()
    resdataup[:,2] = griddata(dataup[:, 0:2], dataup[:,2], resdataup[:, 0:2], method = 'nearest')
    resdataup[:,3] = griddata(dataup[:, 0:2], dataup[:,3], resdataup[:, 0:2], method = 'nearest')
    resdataup = resdataup[resdataup[:,1]<=dedep]

    # combine datasets back together and reassign x,y,z
    data = np.vstack((resdataup,datado))

    # combine filler dataset (from radial basis function) with resampled data and nodes
    filldat = np.ones((len(filldat1),4))
    filldat[:,0] = filldat1[:,0]
    filldat[:,1] = filldat1[:,1]
    filldat[:,2] = filldat1[:,2]
    filldat[:,3] = np.ones(len(filldat))*80
    data = np.vstack((data, filldat))
    data[:, 3][np.isnan(data[:, 3])] = 40

    x = data[:, 0]*1.0
    y = data[:, 1]*1.0
    z = data[:, 2]*1.0
    
    # set weights and knot numbers
    w = 1/data[:, 3]
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    ntx = int(abs(np.floor(xmin)-np.ceil(xmax))*knot_no/111.19)
    nty = int(abs(np.floor(ymin)-np.ceil(ymax))*knot_no/111.19)
    tx = np.linspace(xi.min(), xi.max(), ntx)
    ty = np.linspace(yi.min(), yi.max(), nty)

    # make least squares spline - use smoothbivaraitespline if memory errors occur
    lut = LSQBivariateSpline(x, y, z, tx[1:-1], ty[1:-1], w=w, kx=int(kdeg), ky=int(kdeg))
    #lut = SmoothBivariateSpline(x, y, z, w=w, kx=kdeg, ky=kdeg)

    # evaluate spline at established coordinates
    interpdepths2 = lut.ev(np.ravel(xpts), np.ravel(ypts), dx=0, dy=0)
    
    # filter by actual filter - one above might be different for testing
    sigma = (filt/2.0)/gridsp * 111.19
    interpdepths2.shape = xpts.shape
    
    errorgrid1 = makeErrorgrid(xpts,xyzip,errordata1)
    errorgrid2 = makeErrorgrid(xpts,xyzip,errordata2)
    thickgrid = makeErrorgrid(xpts,xyzip,thickdata)
    
    interpdepths1 = ndimage.filters.gaussian_filter(interpdepths2, sigma)
    errorgrid1 = ndimage.filters.gaussian_filter(errorgrid1, sigma)
    errorgrid2 = ndimage.filters.gaussian_filter(errorgrid2, sigma)
    thickgrid = ndimage.filters.gaussian_filter(thickgrid, sigma)
    
    thickerrorarr = np.zeros((len(thickgrid.flatten()),4))
    thickerrorarr[:,0] = errorgrid1.flatten()
    thickerrorarr[:,1] = errorgrid2.flatten()
    thickerrorarr[:,2] = thickgrid.flatten()
    thickerrorarr[:,3] = (interpdepths2-interpdepths1).flatten()

    mindat, dsddat = getzero(interpdepths1, filt, gridsp, xpts, ypts, xi, slab, meanstk)

    dsddat = np.hstack((dsddat,thickerrorarr))
    interpdepths = interpdepths1.ravel()
    
    # put (s,d,p) coordinates into one array
    newdat = np.zeros((len(interpdepths),3))
    newdat[:,0] = xpts.flatten()
    newdat[:,1] = ypts.flatten()
    newdat[:,2] = interpdepths

    # clip array based on mask - remove if going to make more dense
    pts = np.zeros((len(newdat),2))
    pts[:, 0] = newdat[:, 0]
    pts[:, 1] = newdat[:, 1]
    mask2 = maskdataT(clip, pts)
    #print ('mask2',mask2)
    maskdepths = np.multiply(newdat[:, 2], mask2)
    newdat[:, 2] = maskdepths
    newdat = newdat[~np.isnan(newdat).any(axis=1)]
    
    lomidbound = 10
    upmidbound = mindat['depth'].min()/2
    if slab == 'ker':
        upmidbound = 50
    if slab == 'pam':
        upmidbound = 50
        lomidbound = 200
    if slab == 'hin':
        upmidbound = 200
        lomidbound = 200

    maskdatN, nodemaskN, supplementN, midsupN, extraclip = makevertmask(clip, mindat, xi, newdat, upmidbound, lomidbound, slab)

    if fors == 'second':
        perpdat = np.zeros((len(xpts.flatten()),3))
        perpdat[:,0] = xpts.flatten()
        perpdat[:,1] = ypts.flatten()
        perpdat[:,2] = interpdepths2.flatten()
        mincliplon = maskdatN['s'].min()
        maxcliplon = maskdatN['s'].max()
        
        for index,row in extraclip.iterrows():
            clon, clat = row['s'], row['d']
            if clat == 0:
                continue
            if slab == 'hin' or slab == 'pam' or slab == 'man' or (clon > mincliplon+100 and clon < maxcliplon-100):
                if slab == 'sol':
                    geodata[:,2][(geodata[:,0] > clon - 10)&(geodata[:,0] < clon + 10)&(geodata[:,1] > clat-25)] = np.nan
                    perpdat[:,2][(perpdat[:,0] > clon - 10)&(perpdat[:,0] < clon + 10)&(perpdat[:,1] <= clat-25)] = np.nan
                elif slab == 'hin':
                    geodata[:,2][(geodata[:,0] > clon - 10)&(geodata[:,0] < clon + 10)&(geodata[:,1] > sminz+50)] = np.nan
                    perpdat[:,2][(perpdat[:,0] > clon - 10)&(perpdat[:,0] < clon + 10)&(perpdat[:,1] <= sminz+50)] = np.nan
                else:
                    geodata[:,2][(geodata[:,0] > clon - 10)&(geodata[:,0] < clon + 10)&(geodata[:,1] > sminz-50)] = np.nan
                    perpdat[:,2][(perpdat[:,0] > clon - 10)&(perpdat[:,0] < clon + 10)&(perpdat[:,1] <= sminz-50)] = np.nan
            elif slab == 'izu' or slab == 'sol':
                if clon < mincliplon + 100 or clon > maxcliplon - 100:
                    perpdat[:,2][(perpdat[:,0] > clon - 10)&(perpdat[:,0] < clon + 10)] = np.nan
            elif slab == 'ker':
                if clon < mincliplon+100:
                    perpdat[:,2][(perpdat[:,0] > clon - 10)&(perpdat[:,0] < clon + 10)] = np.nan
        
        geodata1 = geodata[np.isfinite(geodata[:,2])]
        perpdat1 = perpdat[np.isfinite(perpdat[:,2])]
        alldat = np.vstack((geodata1,perpdat1))
        interpdepths3 = griddata(alldat[:, 0:2], alldat[:, 2], perpdat[:, 0:2], method='nearest')

        interpdepths3.shape = xpts.shape
        
        interpdepths1 = ndimage.filters.gaussian_filter(interpdepths3, sigma)

        mindat2, dsddat = getzero(interpdepths1, filt, gridsp, xpts, ypts, xi, slab, meanstk)
        interpdepths = interpdepths1.ravel()

        thickerrorarr[:,3] = (interpdepths2-interpdepths1).flatten()
        dsddat = np.hstack((dsddat,thickerrorarr))

        # put (s,d,p) coordinates into one array
        newdat = np.zeros((len(interpdepths),3))
        newdat[:,0] = xpts.flatten()
        newdat[:,1] = ypts.flatten()
        newdat[:,2] = interpdepths

        # clip array based on mask - remove if going to make more dense
        pts = np.zeros((len(newdat),2))
        pts[:, 0] = newdat[:, 0]
        pts[:, 1] = newdat[:, 1]
        mask2 = maskdataT(clip, pts)
        #print ('mask2',mask2)
        maskdepths = np.multiply(newdat[:, 2], mask2)
        newdat[:, 2] = maskdepths
        newdat = newdat[~np.isnan(newdat).any(axis=1)]
        
        lomidbound = 10
        upmidbound = mindat['depth'].min()/2
        if slab == 'ker':
            upmidbound = 50
        if slab == 'hin' or slab == 'pam':
            upmidbound = 200
            lomidbound = 500

        maskdatN, nodemaskN, supplementN, midsupN, extraclip = makevertmask(clip, mindat, xi, newdat, upmidbound, lomidbound, slab)

        fig = plt.figure(figsize=(20, 25))
        ax1 = fig.add_subplot(311)
        con = ax1.scatter(geodata[:,0],geodata[:,1],c=geodata[:,2],s=50,edgecolors='none',cmap='plasma')
        ax1.scatter(maskdatN['s'].values,maskdatN['d'].values,c='k',s=50,edgecolors='none')
        ax1.set_ylabel('Latitude (Depth)')
        ax1.axis('equal')
        ax1.invert_yaxis()
        plt.grid()
        title = 'Perp Distance from strike-depth Plane'
        ax1.set_title(title)
        cbar = fig.colorbar(con)
        cbar.set_label('Distance from S-D plane')

        ax2 = fig.add_subplot(312)
        con = ax2.scatter(perpdat[:,0],perpdat[:,1],c=perpdat[:,2],s=50,edgecolors='none',cmap='plasma')
        ax2.scatter(maskdatN['s'].values,maskdatN['d'].values,c='k',s=50,edgecolors='none')
        ax2.set_ylabel('Latitude (Depth)')
        ax2.axis('equal')
        ax2.invert_yaxis()
        plt.grid()
        cbar = fig.colorbar(con)
        cbar.set_label('Strike')
        
        ax3 = fig.add_subplot(313)
        con = ax3.scatter(alldat[:,0],alldat[:,1],c=alldat[:,2],s=50,edgecolors='none',cmap='plasma')
        ax3.scatter(maskdatN['s'].values,maskdatN['d'].values,c='k',s=50,edgecolors='none')
        ax3.set_ylabel('Latitude (Depth)')
        ax3.axis('equal')
        ax3.invert_yaxis()
        plt.grid()
        cbar = fig.colorbar(con)
        cbar.set_label('Dip')

        figtitle = '%s_mergetest.png'%slab
        #fig.savefig(figtitle)
        plt.close()

        return newdat, mindat, maskdatN, nodemaskN, newdat, midsupN, maskdatN, midsupN, dsddat
    else:
        if slab == 'sol':
            maskdat, nodemask, supplement, midsup, extraclip = makevertmask(clip, mindat, xi, newdat, upmidbound, lomidbound, slab)
        else:
            maskdat, nodemask, supplement, midsup = makevertmask2(clip, mindat, xi, newdat, upmidbound, lomidbound,slab, distcut)
        return newdat, mindat, maskdat, nodemask, supplement, midsup, maskdatN, midsupN, dsddat

def getzero(interpdepths1, filt, gridsp, xpts, ypts, xi, slab, meanstk):

    interpdepths1.shape = xpts.shape
    
    # calculate strike and dip of grid in new reference frame
    Ygrad, Xgrad = np.gradient(interpdepths1, gridsp, gridsp, edge_order=2)
    strikegrid = np.degrees(np.arctan(Xgrad)) # do this for strike!!!!
    
    quadgrid = strikegrid.flatten()
    tau = meanstk
    for i in range(len(quadgrid)):
        lam = quadgrid[i]
        beta = tau-lam
        if beta < 90 and beta > 0:
            quadgrid[i] = beta
        elif beta < 0 and beta > -90:
            quadgrid[i] = abs(beta)
        else:
            beta = lam-tau
            if beta < 90 and beta > 0:
                quadgrid[i] = beta
            else:
                beta = abs(beta)
                if beta > 90:
                    a = 180-beta
                    quadgrid[i] = a

    quadgrid.shape = strikegrid.shape
    beta = np.radians(quadgrid)
    delta = np.arctan(np.abs(Ygrad))

    dipgrid = np.degrees(np.arctan(np.tan(delta)*np.cos(beta)))
    
    if slab == 'hin' or slab == 'pam':
        dipgrid = np.degrees(np.arctan(np.tan(delta)*np.sin(beta)))

    if slab == 'ker' or slab == 'izu' or slab == 'pam':
        dipgrid[Ygrad>0] *= -1
    else:
        dipgrid[Ygrad<0] *= -1

    dipgrid = 90-dipgrid

    strikegrid = meanstk - strikegrid
    strikegrid[strikegrid < 0] += 360
    if slab == 'pam' or slab == 'man' or slab == 'hin':
        strikegrid -= 180
        strikegrid[strikegrid < 0] += 360
    
    # save information to dataframe
    newdat = np.zeros((len(interpdepths1.ravel()), 5))
    newdat[:, 0], newdat[:, 1], newdat[:, 2] = xpts.flatten(), ypts.flatten(), interpdepths1.ravel()
    newdat[:, 3], newdat[:, 4] = strikegrid.flatten(), dipgrid.flatten()
    
    # initialize arrays for min gradient
    mindips = []
    minstrs = []
    mindeps = []
    mindist = []
    minperp = []
    mingrds = []
    
    # loop through distance along strike columns
    for x in xi:
    
        # avoid choosing points that are at the sh or de extremes of the dataset
        if slab == 'izu' or slab == 'jap':
            these = newdat[(newdat[:,0] == x)&(newdat[:,1] > 200)]
        if slab == 'sol':
            these = newdat[(newdat[:,0] == x)&(newdat[:,1] > 100)]
        else:
            these = newdat[(newdat[:,0] == x)&(newdat[:,1] > 50)]
        
        # get the perpendicular distance from sd plane in this column
        perps = these[:,2]
        grads = np.ones(len(perps))
        
        # loop through column to calculate gradient
        for i in range(2,len(perps)-2):
            p0 = perps[i-2]
            p1 = perps[i-1]
            p2 = perps[i]
            p3 = perps[i+1]
            p4 = perps[i+2]
            pb1 = p2-p1
            pf1 = p3-p2
            pb2 = p1-p0
            pf2 = p4-p3
            
            # calculate average gradient
            gr = abs((pb1+pf1+pb2+pf2)/4)

            grads[i] = gr

        # ensure that the extremes aren't the minimum gradient
        grads[0] = grads[2]
        grads[1] = grads[2]
        grads[-1] = grads[-3]
        grads[-2] = grads[-3]

        # find minimum gradient in dataset, extract associated values
        zerog = these[grads<0.03]
        mindep = zerog[:,1]
        minstr = zerog[:,3]
        minper = zerog[:,2]
        mindip = zerog[:,4]
        mindips.extend(mindip)
        minstrs.extend(minstr)
        mindeps.extend(mindep)
        minperp.extend(minper)
        mindist.extend(np.ones(len(zerog))*x)
        mingrds.extend(grads[grads<0.03])
    
    # save to array and # only take points that are going from + to - or - to +
    mindat = pd.DataFrame({'dist':mindist, 'depth':mindeps, 'perp':minperp, \
                            'strike':minstrs, 'dip':mindips, 'grd':mingrds})
    
    # filter apex depths to create a continuous line (initially jumpy)
    mindepsa = mindat['depth'].values
    deps = ndimage.filters.gaussian_filter(mindepsa, sigma=5)

    # plot depth strike and dip of grid
    fig = plt.figure(figsize=(20, 25))
    ax1 = fig.add_subplot(311)
    con = ax1.scatter(newdat[:,0],newdat[:,1],c=newdat[:,2],s=50,edgecolors='none',cmap='plasma')
    ax1.scatter(mindat['dist'].values,mindat['depth'].values,c='k',s=50,edgecolors='none')
    ax1.set_ylabel('Latitude (Depth)')
    ax1.axis('equal')
    ax1.invert_yaxis()
    plt.grid()
    title = 'Perp Distance from strike-depth Plane'
    ax1.set_title(title)
    cbar = fig.colorbar(con)
    cbar.set_label('Distance from S-D plane')

    ax2 = fig.add_subplot(312)
    con = ax2.scatter(newdat[:,0],newdat[:,1],c=newdat[:,3],s=50,edgecolors='none',cmap='plasma')
    ax2.scatter(mindat['dist'].values,mindat['depth'].values,c='k',s=50,edgecolors='none')
    ax2.set_ylabel('Latitude (Depth)')
    ax2.axis('equal')
    ax2.invert_yaxis()
    plt.grid()
    cbar = fig.colorbar(con)
    cbar.set_label('Strike')
    
    ax3 = fig.add_subplot(313)
    con = ax3.scatter(newdat[:,0],newdat[:,1],c=newdat[:,4],s=50,edgecolors='none',cmap='plasma')
    ax3.scatter(mindat['dist'].values,mindat['depth'].values,c='k',s=50,edgecolors='none')
    ax3.set_ylabel('Latitude (Depth)')
    ax3.axis('equal')
    ax3.invert_yaxis()
    plt.grid()
    cbar = fig.colorbar(con)
    cbar.set_label('Dip')

    figtitle = '%s_depthstrikedipLSQ.png'%slab
    #fig.savefig(figtitle)
    plt.close()
    
    return mindat, newdat

def makevertmask(clip, mindat, xi, newdat, upmidbound, lomidbound, slab):

    minxs = []
    minds = []
    minps = []
    xspac = xi[1] - xi[0]
    tspac = xspac*5
    if slab == 'sol':
        tspac = xspac*15
    for x in xi:
        dat = mindat[(mindat.dist < x + tspac) & (mindat.dist > x - tspac)]
        if len(dat)>0:
            mind = dat['depth'].min()
            mdat = dat[dat.depth == mind]
            minxs.append(x)
            minds.append(mdat['depth'].values[0])
            minps.append(mdat['perp'].values[0])
        else:
            minxs.append(x)
            minds.append(mindat['depth'].max())
            minps.append(-99999)

    fullmindat = pd.DataFrame({'s':minxs, 'd':minds, 'p':minps})
    darr = np.array(fullmindat['d'].values)
    if slab == 'sol':
        sigma = 2
        fullmindat['d'] = ndimage.filters.gaussian_filter(darr, sigma, mode='nearest')
    clip = clip[clip.lat > 20]

    fullclipma = pd.DataFrame({'s':minxs, 'p': 99999})
    fullclipma['d'] = griddata(clip['lon'].values, clip['lat'].values, minxs, method='nearest')

    newclip = pd.DataFrame({'s':xi[::-1], 'd':0, 'p':99999})
    newclip2 = pd.DataFrame({'s':xi[::-1], 'd':0, 'p':99999})
    supp = np.zeros((1,3))*np.nan
    supp2 = np.zeros((1,3))*np.nan
    supp3 = np.zeros((1,3))*np.nan
    for x in xi:
        gmindat = fullmindat[fullmindat.s == x]
        clipdat = fullclipma[fullclipma.s == x]
        dm = gmindat['d'].values[0]
        dc = clipdat['d'].values[0]
        pm = gmindat['p'].values[0]
        if pm != -99999 and dm < dc:
            newclip = pd.concat([newclip,gmindat])
            gmindat2 = gmindat.copy()
            gmindat2['d'] = dm-upmidbound
            newclip2 = pd.concat([newclip2,gmindat2])
            pdat = newdat[newdat[:,0] == x]
            pdatdown = pdat[pdat[:,1]>dm-50]
            pdatmid = pdat[(pdat[:,1]>dm-upmidbound)&(pdat[:,1]<dm+lomidbound)]
            
            supp = np.vstack((supp,pdatdown))
            supp2 = np.vstack((supp2,pdatmid))
        else:
            newclip = pd.concat([newclip,clipdat])
            newclip2 = pd.concat([newclip2,clipdat])

    enddf = newclip.iloc[[0]]
    newclipF = pd.concat([newclip,enddf])
    enddf2 = newclip2.iloc[[0]]
    newclip2F = pd.concat([newclip2,enddf2])
    supp = supp[~np.isnan(supp).any(axis=1)]
    supp2 = supp2[~np.isnan(supp2).any(axis=1)]
    return newclipF, newclip2F, supp, supp2, newclip

def makevertmask2(clip, mindat, xi, newdat, upmidbound, lomidbound, slab, distcut):

    mindist = np.min(newdat[:,0])
    maxdist = np.max(newdat[:,0])
    mindepth = mindat['depth'].min()
    maxdepth = np.max(newdat[:,1])
    if slab == 'man':
        mindat2 = mindat[mindat.depth > 180]
        mindepth = mindat2['depth'].min()
    
    if slab == 'izu':
        mindat2 = mindat[mindat.depth > 300]
        mindepth = mindat2['depth'].min()
    
    if slab == 'sol':
        east = mindat[(mindat.dist<distcut)&(mindat.depth > 175)]
        west = mindat[(mindat.dist>=distcut)&(mindat.depth > 175)]
        mindepthW = west['depth'].min()
        mindepthE = east['depth'].min()
    
    tapdist = 600
    upbuff = 50.0
    if slab == 'sol':
        upbuff = 25.0
    upbuff = 25.0
    slope = (maxdepth-mindepth)/tapdist
    minxs = []
    minds = []
    minps = []
    xspac = xi[1] - xi[0]
    tspac = xspac*5
    for x in xi:
        dat = newdat[(newdat[:,0] < x + tspac) & (newdat[:,0] > x - tspac)]
        if len(dat)>0:
            if slab == 'sol':
                if x < distcut:
                    mdat = dat[dat[:,1] == mindepthE]
                else:
                    mdat = dat[dat[:,1] == mindepthW]
            elif slab == 'kerz':
                if x < np.min(xi)+tapdist:
                    tdepth = -1*slope*(x-np.min(xi))+maxdepth
                    mdat = dat[(dat[:,1] < tdepth+2)&(dat[:,1] > tdepth-2)]
                    #print ('tdepth',tdepth)
                else:
                    mdat = dat[dat[:,1] == mindepth]
            elif slab == 'izuz':
                if x < np.min(xi)+tapdist:
                    tdepth = -1*slope*(x-np.min(xi))+maxdepth
                    mdat = dat[(dat[:,1] < tdepth+0.1)&(dat[:,1] > tdepth-0.1)]
                elif x > np.max(xi)-tapdist:
                    tdepth = -1*slope*(np.max(xi)-x)+maxdepth
                    mdat = dat[(dat[:,1] < tdepth+0.1)&(dat[:,1] > tdepth-0.1)]
                else:
                    mdat = dat[dat[:,1] == mindepth]
            else:
                mdat = dat[dat[:,1] == mindepth]
            minxs.append(x)
            minds.append(mdat[0,1]-upbuff)
            minps.append(mdat[0,2])
        else:
            minxs.append(x)
            minds.append(mindat['depth'].max())
            minps.append(-99999)

    fullmindat = pd.DataFrame({'s':minxs, 'd':minds, 'p':minps})
    clip = clip[clip.lat > 20]

    fullclipma = pd.DataFrame({'s':minxs, 'p': 99999})
    fullclipma['d'] = griddata(clip['lon'].values, clip['lat'].values, minxs, method='nearest')

    lomidbound = lomidbound + upbuff
    newclip = pd.DataFrame({'s':xi[::-1], 'd':0, 'p':99999})
    newclip2 = pd.DataFrame({'s':xi[::-1], 'd':0, 'p':99999})
    supp = np.zeros((1,3))*np.nan
    supp2 = np.zeros((1,3))*np.nan
    supp3 = np.zeros((1,3))*np.nan
    for x in xi:
        gmindat = fullmindat[fullmindat.s == x]
        clipdat = fullclipma[fullclipma.s == x]
        dm = gmindat['d'].values[0]
        dc = clipdat['d'].values[0]
        pm = gmindat['p'].values[0]
        if pm != -99999 and dm < dc:
            newclip = pd.concat([newclip,gmindat])
            gmindat2 = gmindat.copy()
            gmindat2['d'] = dm-upmidbound
            newclip2 = pd.concat([newclip2,gmindat2])
            pdat = newdat[newdat[:,0] == x]
            pdatmid = pdat[(pdat[:,1]>dm-upmidbound)&(pdat[:,1]<dm+lomidbound)]
            pdatdown = pdat[pdat[:,1]>dm-50]
            supp = np.vstack((supp,pdatdown))
            supp2 = np.vstack((supp2,pdatmid))
        else:
            newclip = pd.concat([newclip,clipdat])
            newclip2 = pd.concat([newclip2,clipdat])

    enddf = newclip.iloc[[0]]
    newclip = pd.concat([newclip,enddf])
    enddf2 = newclip2.iloc[[0]]
    newclip2 = pd.concat([newclip2,enddf2])
    supp = supp[~np.isnan(supp).any(axis=1)]
    supp2 = supp2[~np.isnan(supp2).any(axis=1)]
    return newclip, newclip2, supp, supp2

def sdptoxyz(data, x0, y0, meanstk):
    
    dSis = data[:,0]
    depi = data[:,1]
    dPis = data[:,2]
    
    xS = np.zeros(len(dSis))
    yS = np.zeros(len(dSis))
    xP = np.zeros(len(dSis))
    yP = np.zeros(len(dSis))
    
    for i in range(0,len(dPis)):
        xS[i],yS[i] = heading(x0,y0,dSis[i],meanstk)
        xP[i],yP[i] = heading(xS[i],yS[i],dPis[i],meanstk-90)
    
    return xP, yP, depi

def getoutboard(nodes, TRdata, slab):
    lonlist, latlist, depthlist = [],[],[]
    unclist, outlist = [],[]
    print ('sifting through all points and determining which are inboard and outboard of the trench')
    for index,row in nodes.iterrows():
        #if index%100 == 0:
        #    print ('testing index %i out of %i'%(index, len(nodes)))
        try:
            lon,lat,depth,unc = row['lon'], row['lat'], row['depth'], row['stdv']
        except:
            try:
                lon,lat,depth,unc = row['# lon'], row['lat'], row['dep_shift_smooth'], row['dz1']
            except:
                lon,lat = row['lon'], row['lat']
                depth,unc = row['dep_shift_smooth'], row['dz1']
        loc_tr = TRdata[(TRdata.lon > lon-3) & (TRdata.lon < lon+3) & \
                        (TRdata.lat > lat-3) & (TRdata.lat < lat+3)]
        if len(loc_tr)>0:
            #loc_tr['dist'] = gps2dist_azimuth(lat, lon, loc_tr['lat'], loc_tr['lon'])[0]/1000.0
            loc_tr['dist'], tempangles = npcosine(lon, lat, loc_tr['lon'].values, loc_tr['lat'].values)
            mindist = loc_tr['dist'].min()
            loc_tr = loc_tr[loc_tr.dist == mindist]
            lonT = loc_tr['lon'].values[0]
            latT = loc_tr['lat'].values[0]
            azT = loc_tr['az'].values[0]
            thisdist, thisang, latB, lonB = cosine(lonT, latT, lon, lat)
            out = isoutboard(azT, thisang)
        else:
            out = False
        lonlist.append(lon)
        latlist.append(lat)
        depthlist.append(depth)
        outlist.append(out)
        unclist.append(unc)

    shallow = pd.DataFrame({'lon':lonlist, 'lat':latlist, \
                            'depth':depthlist, 'out':outlist, 'stdv':unclist})
    shallowin = shallow[shallow.out == False]
    shallowout = shallow[shallow.out == True]

    return shallowin, shallowout

def tiltedmask(data,slab,spacing, distcut, distcut1):

    xmin = data['lon'].min()
    xmax = data['lon'].max()
    
    ymax = data['depth'].min()

    toparray = np.arange(xmin, xmax, spacing)
    yarr = []
    xarr = []
    for x in toparray:
        datai = data[(data.lon > x-6*spacing) & (data.lon < x+6*spacing)]
        if len(datai)>0:
            ymax = datai['lat'].max()
            yarr.append(ymax)
            xarr.append(x)
        else:
            continue

    maxy = data['lat'].max()
    tiltmask = pd.DataFrame({'lon':xarr, 'lat':maxy})

    if slab == 'sol':
        east = data[data.lon<distcut1]
        west = data[data.lon>distcut]
        cent = data[data.lon>=distcut1]
        cent = cent[cent.lon<=distcut]
        #print ('distcut1,distcut,west,east,cent',distcut1,distcut,west,east,cent)
        maxyE = east['lat'].max()
        maxyW = west['lat'].max()
        maxyC = cent['lat'].max()
        #print ('maxyE,maxyW,maxyC',maxyE,maxyW,maxyC)
        tiltmask.loc[tiltmask.lon < distcut1, 'lat'] = maxyE
        tiltmask.loc[tiltmask.lon > distcut, 'lat'] = maxyW
        tiltmask.loc[(tiltmask.lon >= distcut1)&(tiltmask.lon <= distcut), 'lat'] = maxyC
        maskarr = np.array(tiltmask['lat'].values)
        sigma = 5
        tiltmask['lat'] = ndimage.filters.gaussian_filter(maskarr, sigma, mode='nearest')

    toparr = toparray[::-1]
    tilttop = pd.DataFrame({'lon':toparr, 'lat':np.zeros(len(toparr))})
    tiltmask = pd.concat([tiltmask, tilttop])
    return tiltmask

def maskdataT(clip2, xi):

    clip = clip2.copy()
    #clip.loc[clip.lon < 0, 'lon']+=360
    lons = clip['lon'].values
    lats = clip['lat'].values
    xy = list(zip(lons, lats))
    poly = path.Path(xy)
    temp = poly.contains_points(xi)
    mask1 = (np.zeros(len(temp),) * np.nan)
    mask1[temp] = 1
        
    return mask1

def chunksurface(surfdata, node, T, slab, grid, depname, time, testname, filt, filldat, npass, TR_data, meanBA, kdeg, knot_no, rbfs, shift_out,finorshift,extra,latorlon,mincut,maxcut,maincut):

    if latorlon == 'lat':
        surfdata1 = surfdata[surfdata[:,1]<maxcut]
        surfdata2 = surfdata[surfdata[:,1]>mincut]
    else:
        surfdata1 = surfdata[surfdata[:,0]<maxcut]
        surfdata2 = surfdata[surfdata[:,0]>mincut]

    Surfgrid1, xi1, dl = pySurface3(surfdata1, node, T, slab, grid, depname, time, testname, filt, filldat, npass, TR_data, meanBA, kdeg, knot_no, rbfs, shift_out,finorshift,extra)
    print ('first chunk done')
    Surfgrid2, xi2, dl = pySurface3(surfdata2, node, T, slab, grid, depname, time, testname, filt, filldat, npass, TR_data, meanBA, kdeg, knot_no, rbfs, shift_out,finorshift,extra)
    print ('second chunk done')

    if slab == 'jap' and finorshift == 'fin':
        extrasig = filt/2
        sigma = (extrasig/2.0) / node
        Surfgrid2 = ndimage.filters.gaussian_filter(Surfgrid2, sigma, mode='reflect')

    griddf1 = np.zeros((len(xi1),3))
    griddf2 = np.zeros((len(xi2),3))
    griddf1[:,0] = xi1[:,0]
    griddf1[:,1] = xi1[:,1]
    griddf1[:,2] = Surfgrid1.flatten()
    griddf2[:,0] = xi2[:,0]
    griddf2[:,1] = xi2[:,1]
    griddf2[:,2] = Surfgrid2.flatten()

    if latorlon == 'lat':
        griddf1 = griddf1[griddf1[:,1]<maincut]
        griddf2 = griddf2[griddf2[:,1]>=maincut]
    else:
        griddf1 = griddf1[griddf1[:,0]<maincut]
        griddf2 = griddf2[griddf2[:,0]>=maincut]

    griddf = np.vstack((griddf1,griddf2))
    xmin, xmax = np.min(griddf[:, 0]), np.max(griddf[:, 0])
    ymin, ymax = np.min(griddf[:, 1]), np.max(griddf[:, 1])
    xall = np.arange(np.floor(xmin), np.ceil(xmax)+node, node)
    yall = np.arange(np.floor(ymin), np.ceil(ymax)+node, node)
    n = len(xall)
    m = len(yall)
    xpts, ypts = np.meshgrid(xall, yall)
    xi = np.zeros((m*n, 2))
    xi[:, 0] = xpts.flatten()
    xi[:, 1] = ypts.flatten()

    Surfgrid = griddata(griddf[:, 0:2], griddf[:, 2], xi, method='nearest')
    Surfgrid.shape = xpts.shape
    
    fig = plt.figure(figsize=(25, 20))
    ax1 = fig.add_subplot(131)
    con = ax1.scatter(griddf1[:,0],griddf1[:,1],c=griddf1[:,2],s=50,edgecolors='none',cmap='plasma')
    ax1.scatter(surfdata1[:,0],surfdata1[:,1],c='k',s=5,edgecolors='none')
    ax1.set_ylabel('Latitude')
    ax1.axis('equal')
    plt.grid()
    title = 'Longitude'
    ax1.set_title(title)
    cbar = fig.colorbar(con)
    cbar.set_label('Depth')

    ax2 = fig.add_subplot(132)
    con = ax2.scatter(xi[:,0],xi[:,1],c=Surfgrid.flatten(),s=50,edgecolors='none',cmap='plasma')
    ax2.scatter(surfdata[:,0],surfdata[:,1],c='k',s=5,edgecolors='none')
    ax2.set_ylabel('Latitude')
    ax2.axis('equal')
    plt.grid()
    title = 'Longitude'
    ax2.set_title(title)
    cbar = fig.colorbar(con)
    cbar.set_label('Depth')

    ax3 = fig.add_subplot(133)
    con = ax3.scatter(griddf2[:,0],griddf2[:,1],c=griddf2[:,2],s=50,edgecolors='none',cmap='plasma')
    ax3.scatter(surfdata2[:,0],surfdata2[:,1],c='k',s=5,edgecolors='none')
    ax3.set_ylabel('Latitude')
    ax3.axis('equal')
    plt.grid()
    title = 'Longitude'
    ax3.set_title(title)
    cbar = fig.colorbar(con)
    cbar.set_label('Depth')

    figtitle = 'cuttest.png'
    #fig.savefig(figtitle)
    plt.close()

    return Surfgrid, xi, False

def mkContourClip(nodes, trench, spacing, results, testprint,slab):

    while spacing < 0.05:
        results = results[::2]
        spacing *= 2
        
    indatadat = np.zeros((len(nodes),5)).astype(np.float64)
    indatadat[:,0] = nodes['lon'].values
    indatadat[:,1] = nodes['lat'].values
    indatadat[:,2] = nodes['depth'].values
    try:
        indatadat[:,3] = nodes['sstr'].values
        indatadat[:,4] = nodes['sdip'].values
    except:
        #print('using reference geometry')
        indatadat[:,3] = nodes['ogstr'].values
        indatadat[:,4] = nodes['ogdip'].values

    if slab == 'sol':
        indatadat[:,3] = nodes['ogstr'].values
        indatadat[:,4] = nodes['ogdip'].values
    
    dd = 10
    if slab == 'hel' or slab == 'sul' or slab == 'puy' or slab == 'cot' or slab == 'kurz' or slab == 'mak' or slab == 'hin':
        ds = 60
    elif slab == 'hinz' or slab == 'pamz':
        ds = 15
    else:
        ds = 30

    if slab == 'mak':
        dd = 50
    dxy = 5
    
    rounds = list(range(1,50))
    newres = pd.DataFrame()
    chunklen = int(np.ceil(len(indatadat)/len(rounds)))
    for round in rounds:
        rlons,rlats,rdeps = [],[],[]
        rstrs,rdips = [],[]
        beglen = chunklen*(round-1)
        endlen = chunklen*round
        if endlen > len(indatadat)-1:
            endlen = len(indatadat)-1
        for i in range(beglen,endlen):
            x,y,z,s,d = indatadat[i,0], indatadat[i,1], indatadat[i,2], indatadat[i,3], indatadat[i,4]
            temp = results[(results[:,2] < z+dd) & (results[:,2] > z-dd)]
            if len(temp) > 0:
                temp = temp[(temp[:,0] < x+dxy) & (temp[:,0] > x-dxy) & \
                                 (temp[:,1] < y+dxy) & (temp[:,1] > y-dxy)]
                if len(temp) > 0:
                    temp = temp[(temp[:,3] < s+ds) & (temp[:,3] > s-ds)]
                    if len(temp) > 0:
                        rlons.extend(temp[:,0])
                        rlats.extend(temp[:,1])
                        rdeps.extend(temp[:,2])
                        rstrs.extend(temp[:,3])
                        rdips.extend(temp[:,4])
                    else:
                        continue
                else:
                    continue
            else:
                continue

        thisres = pd.DataFrame({'lon':rlons,'lat':rlats,'depth':rdeps,'strike':rstrs,'dip':rdips})
        #print ('len(thisres),len(newres)',len(thisres),len(newres))
        thisres = thisres.drop_duplicates(['lon','lat','depth'])
        newres = pd.concat([newres,thisres])
        #print ('after duplicates',len(thisres),len(newres))
        del thisres

    newres = newres.drop_duplicates(['lon','lat','depth'])
    #newres.to_csv('%s_cliptesting.csv'%slab,header=True,index=False)

    results = np.zeros((len(newres),5))
    results[:,0] = newres['lon'].values
    results[:,1] = newres['lat'].values
    results[:,2] = newres['depth'].values
    results[:,3] = newres['strike'].values
    results[:,4] = newres['dip'].values

    locstr = griddata(results[:, 0:2], results[:, 3], indatadat[:, 0:2], method='nearest')
    locdip = griddata(results[:, 0:2], results[:, 4], indatadat[:, 0:2], method='nearest')
    dd = 300
    ds = 90
    dxy = 0.5
    dxy2 = 0.2
    
    if slab == 'cas' or slab == 'puy' or slab == 'sol' or slab == 'scoz' or slab == 'man' or slab == 'himz':
        dxy = 0.2
    
    rounds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    newres = pd.DataFrame()
    chunklen = int(np.ceil(len(indatadat)/len(rounds)))
    for round in rounds:
        rlons,rlats,rdeps = [],[],[]
        rstrs,rdips = [],[]
        beglen = chunklen*(round-1)
        endlen = chunklen*round
        if endlen > len(indatadat)-1:
            endlen = len(indatadat)-1
        for i in range(beglen,endlen):
            x,y,z = indatadat[i,0], indatadat[i,1], indatadat[i,2]
            s,d = locstr[i], locdip[i]
            if slab == 'hal' and x > 126:
                tnodes = indatadat[(indatadat[:,0] < x+dxy2) & (indatadat[:,0] > x-dxy2) & \
                             (indatadat[:,1] < y+dxy2) & (indatadat[:,1] > y-dxy2)]
            else:
                tnodes = indatadat[(indatadat[:,0] < x+dxy) & (indatadat[:,0] > x-dxy) & \
                             (indatadat[:,1] < y+dxy) & (indatadat[:,1] > y-dxy)]
            maxdep = np.max(tnodes[:,2])
            temp = results[(results[:,2] < maxdep)&(results[:,4] >= d*0.75)]
            if len(temp) > 0:
                temp1 = temp[(temp[:,0] < x+dxy) & (temp[:,0] > x-dxy) & \
                                 (temp[:,1] < y+dxy) & (temp[:,1] > y-dxy)]
                if len(temp1) > 0:
                    rlons.extend(temp1[:,0])
                    rlats.extend(temp1[:,1])
                    rdeps.extend(temp1[:,2])
                    rstrs.extend(temp1[:,3])
                    rdips.extend(temp1[:,4])
                else:
                    continue
            else:
                continue

        thisres = pd.DataFrame({'lon':rlons,'lat':rlats,'depth':rdeps,'strike':rstrs,'dip':rdips})
        thisres = thisres.drop_duplicates(['lon','lat','depth'])

        newres = pd.concat([newres,thisres])
        del thisres

    newres = newres.drop_duplicates(['lon','lat','depth'])

    if slab == 'ryu':
        newres = newres[(newres.dip > 20) | (newres.depth < 300)]
        newres = newres[(newres.lon > 123)|(newres.lat < 26)]

    if slab == 'himz':
        newres = newres[newres.lon < 90]
    
    if slab == 'hin' or slab == 'pam':
        polyclip = makepolymask(slab,'library/misc/slab_polygons.txt')
        newres.loc[newres.lon < 0, 'lon']+=360
        polyclip.loc[polyclip.lon < 0, 'lon']+=360
        
        pts = np.zeros((len(newres),2))
        pts[:, 0] = newres['lon'].values
        pts[:, 1] = newres['lat'].values
        
        mask = maskdatag(polyclip, pts)
    
        newres['depth'] = newres['depth'].values*mask
        newres = newres[np.isfinite(newres.depth)]
        newres = newres[newres.lon < 76]
        newres = newres.reset_index(drop=True)
        #newres.to_csv('%s_newres.csv'%slab,header=True,index=False)
    
    #newres.to_csv('%s_cliptest.csv'%slab,header=True,index=False)

    newres.loc[newres.lon<0,'lon'] += 360
    return newres


def clippingmask(indataOG, trench, spacing, testprint, slab, fors):
    
    indatadat = np.zeros((len(indataOG),3)).astype(np.float64)
    indatadat[:,0] = np.round(indataOG['lon'].values, 2)
    indatadat[:,1] = np.round(indataOG['lat'].values, 2)
    indatadat[:,2] = np.round(indataOG['depth'].values, 2)
    
    #print (indatadat)
    
    dw = spacing*1.5
    
    dl = 20.0
    if slab == 'solz' or slab == 'sulz':
        dist = 150.0
    elif slab == 'hel' or slab == 'ryu':
        dist = 25
    elif slab == 'him':
        dist = 0.01
    else:
        dist = 75.0

    dist = 0.1

    if slab == 'samz' or slab == 'phiz' or slab == 'cotz' or slab == 'sulz':
        dist = 0.01
    idnot = (np.ones(len(indataOG))*-9999).astype(int)

    if slab == 'sum':
        mindepth = 100
    elif slab == 'hel':
        mindepth = 70
    elif slab == 'kur':
        mindepth = 100
    elif slab == 'ker':
        mindepth = 100
    elif slab == 'alu':
        mindepth = 10.0
    elif slab == 'sol' or slab == 'png' or slab == 'pan':
        mindepth = 10.0
    elif slab == 'sul':
        mindepth = 10.0
    elif slab == 'phi' or slab == 'man':
        mindepth = 30.0
    elif slab == 'cam':
        mindepth = 40
    elif slab == 'cam':
        mindepth = 40
    elif slab == 'cas' or slab == 'sco':
        mindepth = 200
    elif slab == 'cot':
        mindepth = 40
    else:
        mindepth = 25.0

    minempty = 0
    distthresh = 1000
    
    # make trench-side of clipping mask
    trench.loc[trench.lon<0,'lon'] += 360
    trench['az90'] = npaz_perp(trench['az'].values*1.0)
    dists = np.ones(len(trench))*dist
    tlons = trench['lon'].values*1.0
    tlats = trench['lat'].values*1.0
    #lon90, lat90 = npheading(tlons,tlats,az90,dists)
    lon90, lat90=zip(*trench.apply(lambda row: heading(row['lon'], row['lat'], dist, row['az90']), axis=1))
    
    masktrench = pd.DataFrame({'lon':lon90,'lat':lat90})
    
    idlist = []
    for i in range (len(indatadat)):
        nodelon = indatadat[i,0]
        nodelat = indatadat[i,1]
        nodedepth = indatadat[i,2]

        if nodedepth < mindepth:
            idnot[i] = i
        elif slab == 'sam' and nodelon < 287.5 and nodelon > 280 and nodelat < 6 and nodelat > -10:
            idnot[i] = i
        elif slab == 'makz' and nodelat < 29.5:
            idnot[i] = i
        elif slab == 'sol' and nodelat < -7 and nodelon > 148 and nodelon < 150:
            idnot[i] = i
        elif slab == 'ryuz' and nodelon < 125 and nodedepth < 200:
            idnot[i] = i
        
    idnot = idnot[idnot>-999]
    notbytrench = np.delete(indatadat, idnot, 0)
    
    lons = np.ones(len(notbytrench))*-9999
    lats = np.ones(len(notbytrench))*-9999
    northlist = np.ones(len(notbytrench))*-9999
    eastlist = np.ones(len(notbytrench))*-9999
    southlist = np.ones(len(notbytrench))*-9999
    westlist = np.ones(len(notbytrench))*-9999

    lonEmin = 999
    lonEmax = -999
    latEmin = 999
    latEmax = -999

    for i in range(len(notbytrench)):
        dw1 = dw
        dl1 = dl
        dw2 = 1.0
        dl2 = 1.0
        if slab == 'sum':
            dl1 = 5.0
        if slab == 'van' or slab == 'phi':
            dw2 = 0.2
            dl1 = 0.5
        if slab == 'sco':
            dw2 = 0.2
            dl1 = 1.0
        if (slab == 'kur' and notbytrench[i,1] > 35 and notbytrench[i,1] < 50) or (slab == 'jap' and notbytrench[i,1] > 35 and notbytrench[i,1] < 50) or (slab == 'izu' and notbytrench[i,1] > 35):
            dw2 = 3.0
            
        nodelon, nodelat = notbytrench[i,0], notbytrench[i,1]
        NS = indatadat[(indatadat[:,0] < nodelon+dw1) & (indatadat[:,0] > nodelon-dw1)]
        EW = indatadat[(indatadat[:,1] < nodelat+dw1) & (indatadat[:,1] > nodelat-dw1)]
        north = NS[(NS[:,1] > nodelat) & (NS[:,1] < nodelat+dl1)]
        south = NS[(NS[:,1] < nodelat) & (NS[:,1] > nodelat-dl1)]
        east = EW[(EW[:,0] > nodelon) & (EW[:,0] < nodelon+dl1)]
        west = EW[(EW[:,0] < nodelon) & (EW[:,0] > nodelon-dl1)]
        
        n = 0
        if len(north) < 1:
            NS = indatadat[(indatadat[:,0] < nodelon+dw2) & (indatadat[:,0] > nodelon-dw2)]
            north = NS[(NS[:,1] > nodelat+dl2) & (NS[:,1] < nodelat+dl1)]
            if len(north) < 1:
                n += 1
                northlist[i] = 1
        else:
            northlist[i] = 0

        if len(south) < 1:
            NS = indatadat[(indatadat[:,0] < nodelon+dw2) & (indatadat[:,0] > nodelon-dw2)]
            south = NS[(NS[:,1] < nodelat-dl2) & (NS[:,1] > nodelat-dl1)]
            if len(south) < 1:
                n += 1
                southlist[i] = 1
        else:
            southlist[i] = 0

        if len(east) < 1:
            EW = indatadat[(indatadat[:,1] < nodelat+dw2) & (indatadat[:,1] > nodelat-dw2)]
            east = EW[(EW[:,0] > nodelon+dl2) & (EW[:,0] < nodelon+dl1)]
            if len(east) < 1:
                n += 1
                eastlist[i] = 1
        else:
            eastlist[i] = 0

        if len(west) < 1:
            EW = indatadat[(indatadat[:,1] < nodelat+dw2) & (indatadat[:,1] > nodelat-dw2)]
            west = EW[(EW[:,0] < nodelon-dl2) & (EW[:,0] > nodelon-dl1)]
            if len(west) < 1:
                n += 1
                westlist[i] = 1
        else:
            westlist[i] = 0

        if n > minempty:
            lons[i] = nodelon
            lats[i] = nodelat

    lonbool = lons > -999
    maskN = southlist == 0
    maskE = westlist == 0
    maskS = northlist == 0
    maskW = eastlist == 0

    northlist = northlist[lonbool]
    eastlist = eastlist[lonbool]
    southlist = southlist[lonbool]
    westlist = westlist[lonbool]
    
    lons = lons[lons>-999]
    lats = lats[lats>-999]
    
    #print (lons,lats)
    trenchtest = masktrench[(masktrench.lat<=np.max(lats))&(masktrench.lat>=np.min(lats))&(masktrench.lat<=np.max(lons))&(masktrench.lat<=np.max(lats))]
    addfirst = masktrench.iloc[[0]]
    lastpoint = masktrench.iloc[[-1]]
    lastlon = lastpoint['lon'].values[0]
    lastlat = lastpoint['lat'].values[0]
 
    firstlat = addfirst['lat'].values[0]
    firstlon = addfirst['lon'].values[0]

    lastN,lastE,lastS,lastW = 1,1,1,1
    sortedlons = np.ones(len(lons))*-9999
    sortedlats = np.ones(len(lats))*-9999
    sortedangs = np.ones(len(lats))*-9999
    gotOne = True
    alons = np.array(lons)
    alats = np.array(lats)
    
    awest = np.array(westlist)
    aeast = np.array(eastlist)
    anorth = np.array(northlist)
    asouth = np.array(southlist)

    presort = pd.DataFrame({'lon':lons, 'lat':lats, 'depth':1})
    #presort.to_csv('%s_presorted.csv'%slab,header=True,index=False,na_rep=np.nan)
    
    n = 0
    while gotOne == True and slab != 'cas' and slab != 'puy' and slab != 'mak':
        dists, angs = npcosine(lastlon, lastlat, alons, alats)
        distf,angf,lonf,latf = cosine(lastlon,lastlat,firstlon,firstlat)
        
        
        if n>1:
            if lastN == 1:
                maskN = asouth == 0
            else:
                maskN = np.ones(len(dists), dtype=bool)
            if lastE == 1:
                maskE = awest == 0
            else:
                maskE = np.ones(len(dists), dtype=bool)
            if lastS == 1:
                maskS = anorth == 0
            else:
                maskS = np.ones(len(dists), dtype=bool)
            if lastW == 1:
                maskW = aeast == 0
            else:
                maskW = np.ones(len(dists), dtype=bool)

            distsT = dists[maskN & maskE & maskS & maskW]
        
        if len(dists)>0:
            #print (lastlon,lastlat,firstlon,firstlat,distf,np.min(dists))
            if np.min(dists) > distf*0.75:
                gotOne = False
                break
            
            if n>1 and len(distsT)>0:
                minT = np.min(distsT)
                imindista = np.where(dists == minT)
                imindist = imindista[0][0]
            else:
                imindist = np.argmin(dists)
            
            if dists[imindist] < distthresh or n == 0:
            
                lastE, lastW = aeast[imindist], awest[imindist]
                lastN, lastS = anorth[imindist], asouth[imindist]
                            
                lastlon, lastlat = alons[imindist], alats[imindist]
                lastang = angs[imindist]
                
                sortedlons[n] = lastlon
                sortedlats[n] = lastlat
                sortedangs[n] = lastang
                
                alons = np.delete(alons, imindist)
                alats = np.delete(alats, imindist)

                anorth = np.delete(anorth, imindist)
                aeast = np.delete(aeast, imindist)
                asouth = np.delete(asouth, imindist)
                awest = np.delete(awest, imindist)
                
                n+=1
            else:
                gotOne = False
        else:
            gotOne = False
    sortedlons = sortedlons[sortedlons>-999]
    sortedlats = sortedlats[sortedlats>-999]
    sortedangs = sortedlats[sortedlats>-999]
    
    if slab != 'cas' and slab != 'puy' and slab != 'mak' and slab != 'him':
        maskdata = pd.DataFrame({'lon':sortedlons,'lat':sortedlats})
    else:
        maskdata = pd.DataFrame({'lon':lons,'lat':lats})
        if slab == 'cas':
            maskdata = maskdata[maskdata.lat > 38.5]
        if slab == 'manz':
            maskdata = maskdata[maskdata.lat > 12.5]
        if slab == 'puy':
            maskdata = maskdata[maskdata.lat > -50]
        if slab == 'mak':
            lons = np.arange(masktrench['lon'].min()-1.5,masktrench['lon'].max()+1.5)
            lats = np.ones(len(lons))*30
            maskdata = pd.DataFrame({'lon':lons,'lat':lats})
            maskdata = maskdata.sort_values(by=['lon'], ascending=False)
        if slab == 'him':
            trench['az270'] = trench['az90'].values - 180
            trench.loc[trench.az270 < 0, 'az270'] += 360
            lon270, lat270=zip(*trench.apply(lambda row: heading(row['lon'], row['lat'], 175, row['az270']), axis=1))
            
            maskdata = pd.DataFrame({'lon':lon270,'lat':lat270})
            maskdata = maskdata.sort_values(by=['lon'], ascending=False)
        
        
        else:
            maskdata = maskdata.sort_values(by=['lat'], ascending=True)

    #maskdata.to_csv('%s_prefiltered.csv'%slab,header=True,index=False)
    filtno = 10
    filtnum = 0
    n2 = 1

    filtmult = 2
    if slab == 'phi':
        filtmult = 1

    while n2>0:
        maskdata['lon'], maskdata['lat'], n2 = movingav2(maskdata['lon'].values, maskdata['lat'].values,testprint,filtmult)
        filtnum += 1
        #print (filtnum)

    maskdata = maskdata[['lon', 'lat']]
    maskdata = maskdata.reset_index(drop=True)

    clip = pd.concat([masktrench, maskdata, addfirst])
    if slab == 'car':
        clip = clip[clip.lon >= 289]
        clip = clip[clip.lat >= 10]
    if slab == 'mue':
        clip = clip[clip.lon >= 289]
        clip = clip[clip.lon <= 296]
    if slab == 'him':
        clip = clip[clip.lon <= 92]

    #clip.to_csv('%s_postfiltered.csv'%slab,header=True,index=False)
    cliparr = np.zeros((len(clip),2))
    cliparr[:,0] = clip['lon'].values
    cliparr[:,1] = clip['lat'].values
    
    if slab != 'alu' and slab != 'ker':
        cliparr[:,0][cliparr[:,0] > 180] -= 360
    inpolygon = createGridInPolygon2(cliparr, slab, 'library/misc/slab_polygons.txt')
    if slab == 'sum':
        inpolygon  = inpolygon[inpolygon[:,1] < 26.8]
    if slab != 'alu' and slab != 'ker':
        inpolygon[:,0][inpolygon[:,0] < 0] += 360
    
    if inpolygon[0,0] != inpolygon[0,-1] or inpolygon[1,0] != inpolygon[1,-1]:
        inpolygon = np.vstack((inpolygon,inpolygon[0, :]))
        
    inpolydf = pd.DataFrame({'lon':inpolygon[:,0], 'lat':inpolygon[:,1]})

    #inpolydf.to_csv('%s_endclip.csv'%slab,header=True,index=False)

    if slab == 'mue' or slab == 'car':
        if clip['lon'].values[0] != clip['lon'].values[-1] or clip['lat'].values[0] != clip['lat'].values[-1]:
            addfirst = clip.iloc[[0]]
            clip = pd.concat([clip,addfirst])
        return clip
    else:
        return inpolydf

def underclip(output,halgrid):

    printtest = False
    depgrid = gmt.GMTGrid.load(halgrid)
    strgrid, dipgrid = mkSDgrd(depgrid)
    halres = mkSlabData(depgrid, strgrid, dipgrid, printtest)
    halres['depth'] = halres['depth'].values*-1.0
    
    halres = halres[(halres.lon <= np.max(output[:,0])) & (halres.lon >= np.min(output[:,0])) & \
                    (halres.lat <= np.max(output[:,1])) & (halres.lat >= np.min(output[:,1]))]
    
    nanoutput = output[np.isnan(output[:,3])]
    output = output[np.isfinite(output[:,3])]
    for i in range(len(output)):
        #print (i,len(output))
        x,y,z = output[i,0],output[i,1],output[i,3]
        halhere = halres[(halres.lon < x+0.05)&(halres.lon > x-0.05) & \
                         (halres.lat < y+0.05)&(halres.lat > y-0.05)]
        if len(halhere) > 0:
            halhere = halhere[(halhere.depth < z+10)]

            if len(halhere) > 0:
                output[i,2] = np.nan
                output[i,3] = np.nan
                output[i,4] = np.nan
                output[i,5] = np.nan
                output[i,6] = np.nan
                output[i,7] = np.nan
                output[i,8] = np.nan
                output[i,9] = np.nan

    output = np.vstack((output,nanoutput))
    return output

def specialpuyfilt(Surfgrid,xi,filt1,filt2,node):

    xpts = xi[:,0]
    ypts = xi[:,1]
    xpts.shape = Surfgrid.shape
    ypts.shape = Surfgrid.shape
    
    buff = 2
    (rows,cols) = Surfgrid.shape
    upsurf = Surfgrid[ypts>-46-buff]
    dosurf = Surfgrid[ypts<=-46+buff]
    upsurfx = xpts[ypts>-46-buff]
    dosurfx = xpts[ypts<=-46+buff]
    upsurfy = ypts[ypts>-46-buff]
    dosurfy = ypts[ypts<=-46+buff]
    
    uprows = int(len(upsurf)/cols)
    dorows = int(len(dosurf)/cols)
    
    upsurf.shape = (uprows,cols)
    dosurf.shape = (dorows,cols)
    upsurfx.shape = (uprows,cols)
    dosurfx.shape = (dorows,cols)
    upsurfy.shape = (uprows,cols)
    dosurfy.shape = (dorows,cols)

    n = int(filt2/filt1) - 1
    
    sigma = (filt2/2.0) / node / math.sqrt(n+1)
    for i in range(0,n):
        dosurf = ndimage.filters.gaussian_filter(dosurf, sigma, mode='reflect')
        #print ('sigma/node*n',sigma*node*i*2)
    
    upsurf = upsurf[upsurfy>-46]
    dosurf = dosurf[dosurfy<=-46]
    upsurfx = upsurfx[upsurfy>-46]
    dosurfx = dosurfx[dosurfy<=-46]
    upsurfy = upsurfy[upsurfy>-46]
    dosurfy = dosurfy[dosurfy<=-46]
    
    uprows = int(len(upsurf)/cols)
    dorows = int(len(dosurf)/cols)
    
    upsurf.shape = (uprows,cols)
    dosurf.shape = (dorows,cols)

    Surfgrid = np.vstack((upsurf,dosurf))

    Filtgrid = ndimage.filters.gaussian_filter(Surfgrid, sigma, mode='reflect')

    return Filtgrid

def specialkurfilt(Surfgrid,xi,filt1,filt2,node):

    xpts = xi[:,0]
    ypts = xi[:,1]
    xpts.shape = Surfgrid.shape
    ypts.shape = Surfgrid.shape
    
    buff = 2
    (rows,cols) = Surfgrid.shape
    upsurf = Surfgrid[ypts>50-buff]
    dosurf = Surfgrid[ypts<=36+buff]
    misurf = Surfgrid[(ypts>36-buff)&(ypts<50+buff)]
    upsurfy = ypts[ypts>50-buff]
    dosurfy = ypts[ypts<=36+buff]
    misurfy = ypts[(ypts>36-buff)&(ypts<50+buff)]
    misurfx = xpts[(ypts>36-buff)&(ypts<50+buff)]
    
    uprows = int(len(upsurf)/cols)
    dorows = int(len(dosurf)/cols)
    mirows = int(len(misurf)/cols)
    
    upsurf.shape = (uprows,cols)
    dosurf.shape = (dorows,cols)
    misurf.shape = (mirows,cols)
    upsurfy.shape = (uprows,cols)
    dosurfy.shape = (dorows,cols)
    misurfy.shape = (mirows,cols)
    misurfx.shape = (mirows,cols)
    
    misurfL = misurf[misurfx<144+buff]
    misurfR = misurf[misurfx>144-buff]
    misurfLx = misurfx[misurfx<144+buff]
    misurfRx = misurfx[misurfx>144-buff]
    
    Lcols = int(len(misurfL)/mirows)
    Rcols = int(len(misurfR)/mirows)
    
    misurfL.shape = (mirows,Lcols)
    misurfR.shape = (mirows,Rcols)
    misurfLx.shape = (mirows,Lcols)
    misurfRx.shape = (mirows,Rcols)

    sigma = (filt2/2.0) / node
    misurfL = ndimage.filters.gaussian_filter(misurfL, sigma, mode='reflect')
    misurfR = ndimage.filters.gaussian_filter(misurfR, sigma/3, mode='reflect')
    upsurf = ndimage.filters.gaussian_filter(upsurf, sigma/3, mode='reflect')
    
    misurfL = misurfL[misurfLx<144]
    misurfR = misurfR[misurfRx>=144]
    
    Lcols = int(len(misurfL)/mirows)
    Rcols = int(len(misurfR)/mirows)
    
    misurfL.shape = (mirows,Lcols)
    misurfR.shape = (mirows,Rcols)
    upsurf.shape = (uprows,cols)
    
    #print ('misurf shape 2',misurf.shape)
    
    misurf = np.hstack((misurfL,misurfR))
    
    #print ('misurf shape 1',misurf.shape)
    #print ('misurfy shape 1',misurfy.shape)
    #print ('upsurf shape 1',upsurf.shape)
    
    upsurf = upsurf[upsurfy>50]
    dosurf = dosurf[dosurfy<=36]
    misurf = misurf[(misurfy>36)&(misurfy<=50)]
    upsurfy = upsurfy[upsurfy>50]
    dosurfy = dosurfy[dosurfy<=36]
    misurfy = misurfy[(misurfy>36)&(misurfy<=50)]
    
    uprows = int(len(upsurf)/cols)
    dorows = int(len(dosurf)/cols)
    mirows = int(len(misurf)/cols)
    
    upsurf.shape = (uprows,cols)
    dosurf.shape = (dorows,cols)
    misurf.shape = (mirows,cols)

    Surfgrid = np.vstack((upsurf,misurf,dosurf))
    sigma = (filt1/2.0) / node
    Filtgrid = ndimage.filters.gaussian_filter(Surfgrid, sigma, mode='reflect')

    return Filtgrid

def specializufilt(Surfgrid,xi,filt1,filt2,node):

    xpts = xi[:,0]
    ypts = xi[:,1]
    xpts.shape = Surfgrid.shape
    ypts.shape = Surfgrid.shape
    
    buff = 2
    (rows,cols) = Surfgrid.shape
    dosurf = Surfgrid[ypts<=36+buff]
    misurf = Surfgrid[ypts>36-buff]
    dosurfy = ypts[ypts<=36+buff]
    misurfy = ypts[ypts>36-buff]
    misurfx = xpts[ypts>36-buff]
    
    dorows = int(len(dosurf)/cols)
    mirows = int(len(misurf)/cols)
    
    dosurf.shape = (dorows,cols)
    misurf.shape = (mirows,cols)
    dosurfy.shape = (dorows,cols)
    misurfy.shape = (mirows,cols)
    misurfx.shape = (mirows,cols)
    
    misurfL = misurf[misurfx<144+buff]
    misurfR = misurf[misurfx>144-buff]
    misurfLx = misurfx[misurfx<144+buff]
    misurfRx = misurfx[misurfx>144-buff]
    
    Lcols = int(len(misurfL)/mirows)
    Rcols = int(len(misurfR)/mirows)
    
    misurfL.shape = (mirows,Lcols)
    misurfR.shape = (mirows,Rcols)
    misurfLx.shape = (mirows,Lcols)
    misurfRx.shape = (mirows,Rcols)

    sigma = (filt2/2.0) / node
    misurfL = ndimage.filters.gaussian_filter(misurfL, sigma, mode='reflect')
    misurfR = ndimage.filters.gaussian_filter(misurfR, sigma/3, mode='reflect')
    
    misurfL = misurfL[misurfLx<144]
    misurfR = misurfR[misurfRx>=144]
    
    Lcols = int(len(misurfL)/mirows)
    Rcols = int(len(misurfR)/mirows)
    
    misurfL.shape = (mirows,Lcols)
    misurfR.shape = (mirows,Rcols)
    
    misurf = np.hstack((misurfL,misurfR))
    
    dosurf = dosurf[dosurfy<=36]
    misurf = misurf[misurfy>36]
    dosurfy = dosurfy[dosurfy<=36]
    misurfy = misurfy[misurfy>36]
    
    dorows = int(len(dosurf)/cols)
    mirows = int(len(misurf)/cols)
    
    dosurf.shape = (dorows,cols)
    misurf.shape = (mirows,cols)

    Surfgrid = np.vstack((misurf,dosurf))
    sigma = (filt1/2.0) / node
    Filtgrid = ndimage.filters.gaussian_filter(Surfgrid, sigma, mode='reflect')

    return Filtgrid

def makepolymask(slabname,slabfile):

    filerows = []
    slabbounds = []
    with open(slabfile) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[0] == slabname:
                row.pop(0)
                #print ('row',row)
                for i in range(len(row)):
                    slabbounds.append(float(row[i]))

    #print (row)
    
    coords = np.size(slabbounds)
    
    #simple even/odd function
    def is_odd(num):
        return num & 0x1
    
    lons = []
    lats = []
    
    for i in range(coords):
        val = slabbounds[i]
        if is_odd(i):
            lats.append(val)
        else:
            lons.append(val)

    #print ('slabname, coords, lons, lats')
    #print (slabname)
    #print (coords)
    #print (lons)
    #print (lats)
    #print (slabbounds)
    polyclip = pd.DataFrame({'lon':lons, 'lat':lats})

    return polyclip

def preshiftfill(tmp_res, emptynodes, refdeps, mindip, dipthresh):

    # merge reference depths for nodes with data with filled node data frame
    tmp_res = pd.merge(tmp_res, refdeps)
    
    # initialize arrays for interpolating over
    fulln = np.zeros((len(tmp_res),15))
    emptn = np.zeros((len(emptynodes),15))
    
    # set original search locations and vectors (nodes with data)
    fulln[:,0] = tmp_res['lon'].values*1.0
    fulln[:,1] = tmp_res['lat'].values*1.0
    fulln[:,2] = tmp_res['ogdep'].values*1.0
    fulln[:,3] = tmp_res['ogstr'].values*1.0
    fulln[:,4] = tmp_res['ogdip'].values*1.0
    
    # set original search locations and vectors (nodes w/o data)
    emptn[:,0] = emptynodes['lon'].values*1.0
    emptn[:,1] = emptynodes['lat'].values*1.0
    emptn[:,2] = emptynodes['ogdep'].values*1.0
    emptn[:,3] = emptynodes['ogstr'].values*1.0
    emptn[:,4] = emptynodes['ogdip'].values*1.0
    
    # modify search vectors according to perp search dip bounds
    emptn[:,4][emptn[:,4] > dipthresh] = 90.0
    emptn[:,4][emptn[:,4] < mindip] = 0.0
    
    # add other info to nodes with data array
    fulln[:,5] = tmp_res['bzlon'].values*1.0
    fulln[:,6] = tmp_res['bzlat'].values*1.0
    fulln[:,7] = tmp_res['depth'].values*1.0
    fulln[:,8] = tmp_res['stdv'].values*1.0
    fulln[:,9] = tmp_res['centsurf'].values*1.0
    fulln[:,12] = tmp_res['onlyto'].values*1
    
    # get r, phi, theta values for search points and peak points (nodes w data)
    r1 = 6371 - tmp_res['ogdep'].values
    r2 = 6371 - tmp_res['depth'].values

    p1 = np.radians(tmp_res['lon'].values)
    p2 = np.radians(tmp_res['bzlon'].values)

    t1 = np.radians(np.abs(tmp_res['lat'].values - 90.0))
    t2 = np.radians(np.abs(tmp_res['bzlat'].values - 90.0))
    
    # find distance between reference point and center of benioff zone
    dist = r1*r1 + r2*r2 - 2*r1*r2*(np.sin(t1)*np.sin(t2)*np.cos(p1-p2) + np.cos(t1)*np.cos(t2))
    
    # determine shift direction (inboard or outboard from point)
    inorout = np.ones(len(dist))*-1.0
    for i in range(len(fulln)):
        reflon, reflat = fulln[i,0], fulln[i,1]
        bzlon, bzlat = fulln[i,5], fulln[i,6]
        xydist, ang, x1, y1 = cosine(reflon, reflat, bzlon, bzlat)
        outhere = isoutboard(fulln[i,3], ang)
        if outhere:
            inorout[i] = 1.0
    
    # add outboard and distance values to array
    fulln[:,10] = dist
    fulln[:,11] = inorout
    
    # interpolate values to nodes without data
    emptn[:, 8] = griddata(fulln[:, 0:2], fulln[:, 8], emptn[:, 0:2], method='nearest')
    emptn[:, 9] = griddata(fulln[:, 0:2], fulln[:, 9], emptn[:, 0:2], method='nearest')
    emptn[:, 10] = griddata(fulln[:, 0:2], fulln[:, 10], emptn[:, 0:2], method='nearest')
    emptn[:, 11] = griddata(fulln[:, 0:2], fulln[:, 11], emptn[:, 0:2], method='nearest')
    emptn[:, 12] = griddata(fulln[:, 0:2], fulln[:, 12], emptn[:, 0:2], method='nearest')
    #emptn[:, 10] *= emptn[:, 11]
    emptn[:, 10] /= 1000
    emptn[:, 12][emptn[:, 12] < 0.5] = 0
    emptn[:, 12][emptn[:, 12] >= 0.5] = 1

    # loop through empty nodes and find "interpolated" center of benioff zone
    for i in range(len(emptn)):
        ipslon, ipslat, ipsdep = emptn[i, 0], emptn[i, 1], emptn[i, 2]
        istrike, idip, idist = emptn[i, 3], emptn[i, 4], emptn[i, 10]
        ibzlon, ibzlat, ibzdep = pointShift(ipslon, ipslat, ipsdep, idip, istrike, idist)
        emptn[i, 5] = ibzlon
        emptn[i, 6] = ibzlat
        emptn[i, 7] = ibzdep

    emptynodes['bzlon'] = emptn[:, 5]
    emptynodes['bzlat'] = emptn[:, 6]
    emptynodes['depth'] = emptn[:, 7]
    emptynodes['stdv'] = emptn[:, 8] * 10
    emptynodes['centsurf'] = emptn[:, 9]
    emptynodes['onlyto'] = emptn[:, 12]
    emptynodes['smag1'] = emptn[:, 10]
    emptynodes['inorout'] = emptn[:, 11]

    emptynodes = emptynodes[np.isfinite(emptynodes.bzlon)]

    return emptynodes

