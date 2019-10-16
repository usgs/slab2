import pandas as pd
import numpy as np
import os.path
import math
from obspy.geodetics.base import gps2dist_azimuth

# Eliminates events in dfo that are found in dfm
def removematches(dfo, dfm):
    ind = dfo.lat.isin(dfm.lat) & dfo.lon.isin(dfm.lon)
    dfo = dfo[~ind]
    return dfo

# Finds events in dfo that are in dfm
def findNotMatches(dfo,dfm):
    ind = dfo.lat.isin(dfm.lat) & dfo.lon.isin(dfm.lon)
    dfo = dfo[~ind]
    return dfo

def findNotMatches2(dfo,dfm):
    ind = dfo.lat.isin(dfm.lonlat)
    dfo = dfo[~ind]
    return dfo

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
        try:
            ang = (math.cos(cl2) - (dist * math.cos(cl1))) / (math.sin(dist2) * math.sin(cl1))
        except:
            cl1 += 0.0001
            cl2 += 0.0001
            dist2 += 0.0001
            ang = (math.cos(cl2) - (dist * math.cos(cl1))) / (math.sin(dist2) * math.sin(cl1))
    else:
        ang = 1.0
    if ang < -1:
        ang = -1.0
    if ang > 1:
        ang = 1.0
    ang = math.acos(ang)
    return dist2, ang

def npcosrule(d2r, lon1, lat1, lon2, lat2):
    
    # Written GLM 4.25.17
    # Added here KLH 09/25/2019
    
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


def npcosine(lon1, lat1, lon2, lat2):
    
    # Written GLM 4.25.17
    # Added here KLH 09/25/2019
    
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

def cosine(lon1,lat1,lon2,lat2):
    # Logic by Gavin Hayes
    if lon1 > 180:
        lon1 = lon1 - 360
    if lon2 > 180:
        lon2 = lon2 - 360
    
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

def outboard(az,ang):
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

def refilter(OBdata,TRdata):
    for index,row in OBdata.iterrows():
        lonB, latB = row['lon'], row['lat']
        
        loc_tr = TRdata[(TRdata.lon > lonB-2) & (TRdata.lon < lonB+2) & (TRdata.lat > latB-2) & (TRdata.lat < latB+2)]
        #loc_tr['dist'] = gps2dist_azimuth(latB,lonB,loc_tr['lat'],loc_tr['lon'])[0]/1000.0
        loc_tr['dist'], tempangles = npcosine(lonB, latB, loc_tr['lon'].values, loc_tr['lat'].values)
        mindist = loc_tr['dist'].min()
        loc_tr = loc_tr[loc_tr.dist == mindist]
        lonT, latT, azT = loc_tr['lon'].values[0], loc_tr['lat'].values[0], loc_tr['az'].values[0]
        thisdist, thisang, lat, lon = cosine(lonT,latT,lonB,latB)

        out = outboard(azT,thisang)
        if not out:
            print('removed point',lonB,latB,lonT,latT)
            OBdata.drop(index, inplace=True)
            
    return OBdata


vcosine = np.vectorize(cosine)
voutboard = np.vectorize(outboard)

''' edit next two lines to create bathymetry files '''
# list of slabs to create bathymetry files for
slablist = ['hal','phi','sul','cot','png','sum','sol','man','ryu','kur','izu','van','ker','puy','alu','cas','cam','car','mue','sam','hel','cal','mak','sco']

# trench file to generate search around
trenches = pd.read_csv('trenches_usgs_2017_depths.csv')
''' edit previous two lines to create bathymetry files '''

bathymetry = pd.read_csv('totalBA.csv')

trenches.loc[trenches.lon < 0, 'lon']+=360
bathymetry.loc[bathymetry.lon < 0, 'lon']+=360

for slab in slablist:
    print('slab',slab)
    loc_tr = trenches[trenches.slab == slab]
    
    print('loc_tr',loc_tr)
    lonmin = loc_tr['lon'].min() - 2
    lonmax = loc_tr['lon'].max() + 2
    latmin = loc_tr['lat'].min() - 2
    latmax = loc_tr['lat'].max() + 2
    
    loc_ba = bathymetry[(bathymetry.lon < lonmax) & (bathymetry.lon > lonmin) & (bathymetry.lat < latmax) & (bathymetry.lat > latmin)]
    
    loc_ba['baID'] = np.arange(len(loc_ba))
    
    slabBA = pd.DataFrame(columns = ['lon','lat','elev','thickness','depth','baID'])
    
    print('slab, loc_ba',slab, loc_ba)

    print('len(loc_tr)',len(loc_tr))
    n = 0
    for index,row in loc_tr.iterrows():
        lonT, latT, azT = row['lon'], row['lat'], row['az']
        
        loc_ba1 = loc_ba[(loc_ba.lon > lonT-2) & (loc_ba.lon < lonT+2) & (loc_ba.lat > latT-2) & (loc_ba.lat < latT+2)]
        #loc_ba1['dist'] = gps2dist_azimuth(latT,lonT,loc_ba1['lat'],loc_ba1['lon'])[0]/1000.0
        loc_ba1['dist'], tempangles = npcosine(lonT, latT, loc_ba1['lon'].values, loc_ba1['lat'].values) # commented out KLH 09/25/2019
        slab_ba = loc_ba1[loc_ba1.dist < 200]
        lonsBA = slab_ba['lon'].values
        latsBA = slab_ba['lat'].values

        if len(slab_ba)>0:
            slab_ba['dist'], slab_ba['ang'], lats, lons = vcosine(lonT,latT,lonsBA,latsBA)
            slab_ba['outboard'] = voutboard(azT,slab_ba['ang'].values)
            slab_ba = slab_ba[slab_ba.outboard == True]
        
        #slab_ba = removematches(slab_ba,slabBA)
        BAframes = [slabBA, slab_ba]
        slabBA = pd.concat(BAframes,sort=True)
        slabBA = slabBA.reset_index(drop=True)
    
        n+=1
        print(n)


    slabBA2 = slabBA.drop_duplicates(['baID'])

    slabBA2.to_csv('%s_BAth2.csv' % slab,header=True,index=False,float_format='%0.4f')
    itterBA = slabBA2['baID'].values
    BAids = slabBA['baID'].values


    slabBA = slabBA.drop_duplicates(['baID'])

    print('len(slabBA)',len(slabBA))

    slabBA['unc'] = 5.0
    slabBA = slabBA[['lon','lat','depth','unc','elev','thickness']]

    slabBA1 = refilter(slabBA,loc_tr)
    slabBA.to_csv('%s_BA_bath_all.csv' % slab,header=True,index=False,float_format='%0.4f')
    slabBA1.to_csv('%s_BA_bath.csv' % slab,header=True,index=False,float_format='%0.4f')




