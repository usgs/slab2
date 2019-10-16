#!/usr/bin/env python

import math
import os.path
import numpy as np
import pandas as pd
from neicmap.distance import sdist
import argparse

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


def cosrule(d2r,lat1,lon1,lat2,lon2):
    # Logic by Gavin Hayes
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

def cosine(lon1,lat1,lon2,lat2):
    # Logic by Gavin Hayes 
    if lon1 > 180:
        lon1 = lon1 - 360
    if lon2 > 180:
        lon2 = lon2 - 360
    
    if abs(lon1-lon2) < 0.001:
        lon2+=0.01
    if abs(lat1-lat2) < 0.001:
        lat2+=0.01

    d2r = (math.pi/180)
    r2d = (180/math.pi)
    ddlon = lon1 - lon2
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

def newaz(df):
    lons = df['lon'].values
    lats = df['lat'].values
    az = np.zeros(len(df))
    for i in range(len(lons)-1):
        lon1, lat1 = lons[i], lats[i]
        lon2, lat2 = lons[i+1], lats[i+1]
        dist, ang, lat1, lon1 = cosine(lon1,lat1,lon2,lat2)
        az[i] = ang

    az[-1] = az[-2]
    az[0] = az[1]
    df['az'] = az
    return df

def main(args):

    oldtrenchfile = args.newtrenchlist
    oldtrench = pd.read_csv(oldtrenchfile)
    oldtrench = zerothreesixty(oldtrench)
    slabs = oldtrench['slab'].values
    slablist = mylist = list(set(list(slabs)))
    print ('generating trenches for this list of slab models:', slablist)
    allslabs = pd.DataFrame()

    # sort points for each slab in slablist, sort direction changes for each region.
    for slab in slablist:
        thistrench = oldtrench[oldtrench.slab == slab]
        
        print ('making trench for:',slab)
        
        bound = thistrench['bound'].values[0]
        firstlon = thistrench['lon'].values[-1]
        firstlat = thistrench['lat'].values[-1]
        slons,slats = [],[]
        slons.append(firstlon)
        slats.append(firstlat)
        thistrench = thistrench[(thistrench.lon != firstlon)|(thistrench.lat != firstlat)]
        while (len(thistrench))>0:
            thistrench['dist'] = sdist(firstlat, firstlon, thistrench['lat'], thistrench['lon'])/1000.0
            closest = thistrench[thistrench.dist == thistrench['dist'].min()]
            firstlon = closest['lon'].values[0]
            firstlat = closest['lat'].values[0]
            thistrench = thistrench[(thistrench.lon != firstlon)|(thistrench.lat != firstlat)]
            slons.append(firstlon)
            slats.append(firstlat)
            
        thistrench = pd.DataFrame({'lon':slons,'lat':slats,'slab':slab,'bound':bound})
        
        thisnew = newaz(thistrench)
        thisnew['number'] = range(len(thisnew))
        thisnew = thisnew.sort_values(by=['number'],ascending=False)
        
        allslabs = pd.concat([allslabs,thisnew],sort=True)

    allslabs = oneeighty(allslabs)
    allslabs = allslabs[['lon','lat','az','bound','slab','number']]
    allslabs.to_csv(args.trenchfile,header=True,index=False,na_rep=np.nan)

if __name__=='__main__':
    desc = '''
        This script is used to sort the trench coordinates associated with a 
        given slab according to the right hand rule, then calculate the azimuth 
        values between consecutive points.
        
        The process requires manual additions/edits to the script itsself in 
        order to document the way that the trench points should be sorted.
        
        Required arguments include:
            -n newtrenchlist: the file containing the original unsorted list 
                              of trench coordinates
            -t trenchfile: the output file to write the new sorted list of 
                           coordinates to
            
        The script assumes that the last point listed for each slab is the start
        of the trench according to the right hand rule.
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-n', '--newtrenchlist', dest='newtrenchlist', type=str,
                        required=True, help='new unsorted trench coord list')
                        
    parser.add_argument('-t', '--trenchfile', dest='trenchfile', type=str,
                        required=True, help='file to save output to')
                        
    pargs = parser.parse_args()

    main(pargs)










