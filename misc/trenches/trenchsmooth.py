#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import argparse

def cosrule(d2r,lat1,lon1,lat2,lon2):
    
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

def cosine(lon1,lat1,lon2,lat2):
    
    ''' Arguments:  lon1 - latitude point that the angle is referenced from
                    lat1 - longitude point that the angle is referenced from
                    lon2 - latitude point that the angle is going to (clockwise from 0 degrees)
                    lat2 - longitude point that the angle is going to (clockwise from 0 degrees)
        
        Returns:    dist - great circle distance between the two lat/lon points
                    ang - the angle between the two points (clockwise from 0 degrees from lat1/lon1 point)
                    lon1 - same longitude as input argument, used in some applications
                    lat1 - same latitude as input argument, used in some applications  '''
    
    # Ensuring that azimuths are between 0 and 360
    if lon1 > 180:
        lon1 = lon1 - 360
    if lon2 > 180:
        lon2 = lon2 - 360

    # Creating degrees/radians conversion constants
    d2r = (math.pi/180)
    r2d = (180/math.pi)
    ddlon = lon1 - lon2

    # Getting distance and angle between the two points (in degrees)
    dist,ang = cosrule(d2r,lat1,lon1,lat2,lon2)
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

def datelinecross(x):
    ''' Arguments:  x - longitude value (positive or negative)
        
        Returns:    x - a positive longitude. Stays the same if the input was positive,
                        is changed to positive if the input was negative '''
    
    if x<0:
        return x+360
    else:
        return x

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
    data['lon']=data.apply(lambda row: datelinecross(row['lon']),axis=1)
    return data

def oneeighty(data):
    data['lon']=data.apply(lambda row: meridiancross(row['lon']),axis=1)
    return data

def northernaz(data):
    data['az']=data.apply(lambda row: northcross(row['az']),axis=1)
    return data

def notnorthanymore(data):
    data['az']=data.apply(lambda row: unnorthcross(row['az']),axis=1)
    return data

def movingav(x):
    n = 0
    x2 = np.copy(x)
    for i in range(1,len(x)-1):
        thisaz = x[i]
        lastaz = x[i-1]
        nextaz = x[i+1]
        lastdiff = abs(thisaz-lastaz)
        nextdiff = abs(thisaz-nextaz)
        if thisaz < lastaz and thisaz < nextaz and (nextdiff>20 and lastdiff >10):
            x2[i] = (nextaz+lastaz)/2.0
            print('changed this az 1 lastaz,thisaz,nextaz',lastdiff,lastaz,thisaz,nextaz,nextdiff)
            n+= 1
        elif thisaz > lastaz and thisaz > nextaz and (nextdiff>20 and lastdiff >10):
            x2[i] = (nextaz+lastaz)/2.0
            print('changed this az 2 lastaz,thisaz,nextaz',lastdiff,lastaz,thisaz,nextaz,nextdiff)
            n+= 1
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
    return x2, n

def filldf(data,ns_ew):
    n = 0
    newtrench = pd.DataFrame(columns = ['lon','lat','az','bound','slab'])
    lons = data['lon'].values
    lats = data['lat'].values
    azzs = data['az'].values
    bound = data['bound'].values[0]
    slab = data['slab'].values[0]
    newlons = []
    newlats = []
    newazzs = []
    for i in range(0,len(data)-1):
        thislon = lons[i]
        nextlon = lons[i+1]
        thislat = lats[i]
        nextlat = lats[i+1]
        thisazz = azzs[i]
        nextazz = azzs[i+1]
    
        dist,ang,lon1,lat1 = cosine(thislon,thislat,nextlon,nextlat)
        if dist > 20:
            n+= 1
            print('dist,ang,lon1,lat1',dist,ang,lon1,lat1)
            meanlon = (thislon+nextlon)/2.0
            meanlat = (thislat+nextlat)/2.0
            meanazz = (thisazz+nextazz)/2.0
    
            newlons.append(thislon)
            newlons.append(meanlon)
            newlats.append(thislat)
            newlats.append(meanlat)
            newazzs.append(thisazz)
            newazzs.append(meanazz)
        else:
            newlons.append(thislon)
            newlats.append(thislat)
            newazzs.append(thisazz)
    
    newlons.append(lons[-1])
    newlats.append(lats[-1])
    newazzs.append(azzs[-1])
    newdata = pd.DataFrame({'lon':newlons,'lat':newlats,'az':newazzs,'bound':bound,'slab':slab})
    return newdata, n

def main(args):

    sorf = args.smoothorfill
    newtrenchfile = args.newtrenchlist
    trenchfile = args.trenchfile
    trenches = pd.read_csv(newtrenchfile)
    slabs = trenches['slab'].values
    slablist = mylist = list(set(list(slabs)))
    newtrenches = pd.DataFrame()

    for slab in slablist:
    
        thisdf = trenches[trenches.slab == slab]
        if slab == 'alu' or slab =='ker':
            thisdf = zerothreesixty(thisdf)
        strike = np.mean(thisdf['az'].values)
        if np.min(thisdf['az'].values) < 45 and np.max(thisdf['az'].values) > 315:
            ns_ew = 'N'
            thisdf = northernaz(thisdf)
        elif (strike < 135. and strike > 45.) or (strike < 315. and strike > 225.):
            ns_ew = 'EW'
        else:
            ns_ew = 'S'

        if sorf == 's':
            n = 10
            while n > 0:
                ogstrikes = thisdf['az'].values
                avstrikes, n = movingav(ogstrikes)
                thisdf['az'] = avstrikes
        elif sorf == 'f':
            n = 10
            while n > 0:
                thisdf, n = filldf(thisdf,ns_ew)
        else:
            print('sorf must be s or f')
            print('enter -s s or -s f at input')
            print('Exiting ..... ')
            exit()
        
        thisdf = oneeighty(thisdf)
        if ns_ew == 'N':
            thisdf = notnorthanymore(thisdf)
        frames = [newtrenches,thisdf]
        newtrenches = pd.concat(frames,sort=True)

    newtrenches = newtrenches[['lon','lat','az','bound','slab']]
    newtrenches.to_csv(trenchfile,header=True,index=False,na_rep=np.nan,float_format='%.4f')


if __name__=='__main__':
    desc = '''
        This script is used to take a sorted trench segment (or many in the same file) and smooth out any outlying azimuth values calculated from sortpoints.py. Once the segments are all sufficiently smoothed, this is also intended to be used as a tool for filling in the trench locations. This is done by running the script until each point is within 20 km of the next point. The density of points is necessary for the slab2 modeling algorithm.
        
        The process requires manual additions/edits to the script itsself in 
        order to document the way that the trench points should be sorted.
        
        Required arguments include:
            -n newtrenchlist: the file containing the original unsorted list 
                              of trench coordinates
            -t trenchfile: the output file to write the new sorted list of 
                           coordinates to
            -s smoothorfill: s to smooth newtrenchlist
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-n', '--newtrenchlist', dest='newtrenchlist', type=str,
                        required=True, help='new unsorted trench coord list')
                        
    parser.add_argument('-t', '--trenchfile', dest='trenchfile', type=str,
                        required=True, help='file to save output to')
                        
    parser.add_argument('-s', '--smoothorfill', dest='smoothorfill', type=str,
                        required=True, help='s to smooth newtrenchlist, f to fill newtrenchlist')
                        
    pargs = parser.parse_args()

    main(pargs)


