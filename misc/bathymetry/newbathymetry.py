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
            print 'removed point',lonB,latB,lonT,latT
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
    print 'slab',slab
    loc_tr = trenches[trenches.slab == slab]
    
    print 'loc_tr',loc_tr
    lonmin = loc_tr['lon'].min() - 2
    lonmax = loc_tr['lon'].max() + 2
    latmin = loc_tr['lat'].min() - 2
    latmax = loc_tr['lat'].max() + 2
    
    loc_ba = bathymetry[(bathymetry.lon < lonmax) & (bathymetry.lon > lonmin) & (bathymetry.lat < latmax) & (bathymetry.lat > latmin)]
    
    loc_ba['baID'] = np.arange(len(loc_ba))
    
    slabBA = pd.DataFrame(columns = ['lon','lat','elev','thickness','depth','baID'])
    
    print 'slab, loc_ba',slab, loc_ba

    print 'len(loc_tr)',len(loc_tr)
    n = 0
    for index,row in loc_tr.iterrows():
        lonT, latT, azT = row['lon'], row['lat'], row['az']
        
        loc_ba1 = loc_ba[(loc_ba.lon > lonT-2) & (loc_ba.lon < lonT+2) & (loc_ba.lat > latT-2) & (loc_ba.lat < latT+2)]
        #loc_ba1['dist'] = gps2dist_azimuth(latT,lonT,loc_ba1['lat'],loc_ba1['lon'])[0]/1000.0
        loc_ba1['dist'], tempangles = npcosine(lonT, latT, loc_ba1['lon'].values, loc_ba1['lat'].values)
        slab_ba = loc_ba1[loc_ba1.dist < 200]
        lonsBA = slab_ba['lon'].values
        latsBA = slab_ba['lat'].values

        if len(slab_ba)>0:
            slab_ba['dist'], slab_ba['ang'], lats, lons = vcosine(lonT,latT,lonsBA,latsBA)
            slab_ba['outboard'] = voutboard(azT,slab_ba['ang'].values)
            slab_ba = slab_ba[slab_ba.outboard == True]
        
        #slab_ba = removematches(slab_ba,slabBA)
        BAframes = [slabBA, slab_ba]
        slabBA = pd.concat(BAframes)
        slabBA = slabBA.reset_index(drop=True)
    
        n+=1
        print n


    slabBA2 = slabBA.drop_duplicates(['baID'])

    slabBA2.to_csv('%s_BAth2.csv' % slab,header=True,index=False,float_format='%0.4f')
    itterBA = slabBA2['baID'].values
    BAids = slabBA['baID'].values


    slabBA = slabBA.drop_duplicates(['baID'])

    print 'len(slabBA)',len(slabBA)

    slabBA['unc'] = 5.0
    slabBA = slabBA[['lon','lat','depth','unc','elev','thickness']]

    slabBA1 = refilter(slabBA,loc_tr)
    slabBA.to_csv('%s_BA_bath_all.csv' % slab,header=True,index=False,float_format='%0.4f')
    slabBA1.to_csv('%s_BA_bath.csv' % slab,header=True,index=False,float_format='%0.4f')




