import math
import os.path
import numpy as np
import pandas as pd


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
    print 'lendfbefore',len(df)
    for i in range(len(lons)-1):
        lon1, lat1 = lons[i], lats[i]
        lon2, lat2 = lons[i+1], lats[i+1]
        dist, ang, lat1, lon1 = cosine(lon1,lat1,lon2,lat2)
        az[i] = ang
    az[-1] = az[-2]
    az[0] = az[1]
    df['az'] = az
    print 'lendfafter',len(df)
    print df
    return df

oldtrenchfile = 'ryutrench.csv'
#oldtrenchfile = 'newtrenches.csv'
oldtrench = pd.read_csv(oldtrenchfile)
oldtrench = zerothreesixty(oldtrench)
slabs = oldtrench['slab'].values
slablist = mylist = list(set(list(slabs)))
print slablist
allslabs = pd.DataFrame()

for slab in slablist:
    thistrench = oldtrench[oldtrench.slab == slab]
    if slab == 'man':
        thistrench = thistrench.sort(['lat'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'sco':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'ryu':
        '''
        lattrench = thistrench[thistrench.lat > 25]
        lontrench1 = thistrench[thistrench.lat <= 25]
        lontrench2 = thistrench[thistrench.lat >= 25]
        lattrench = lattrench.sort(['lat'],ascending=False)
        lontrench = lontrench.sort(['lon'],ascending=False)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        thisnew = pd.concat([thisnewlat,thisnewlon])
        '''
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'van':
        lattrench = thistrench[thistrench.lat > -21]
        lontrench = thistrench[thistrench.lat <= -21]
        lattrench = lattrench.sort(['lat'],ascending=True)
        lontrench = lontrench.sort(['lon'],ascending=False)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        thisnew = pd.concat([thisnewlon,thisnewlat])
    elif slab == 'sam':
        thistrench = thistrench.sort(['lat'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'cas':
        thistrench = thistrench.sort(['lat'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'sum':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
        lattrench1 = thistrench[thistrench.lon < 105]
        lontrench1 = thistrench[(thistrench.lon >= 105) & (thistrench.lon <= 120)]
        lattrench2 = thistrench[(thistrench.lon > 120)]
        lattrench1 = lattrench1.sort(['lat'],ascending=True)
        lontrench1 = lontrench1.sort(['lon'],ascending=False)
        lattrench2 = lattrench2.sort(['lat'],ascending=False)
        thisnewlat1 = newaz(lattrench1)
        thisnewlon1 = newaz(lontrench1)
        thisnewlat2 = newaz(lattrench2)
        thisnew = pd.concat([thisnewlat2,thisnewlon1,thisnewlat1])
    elif slab == 'kur':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'him':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'sol':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'phi':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'hal':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'izu':
        lattrench = thistrench[thistrench.lat > 13]
        lontrench = thistrench[thistrench.lat <= 13]
        lattrench = lattrench.sort(['lat'],ascending=False)
        lontrench = lontrench.sort(['lon'],ascending=False)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        thisnew = pd.concat([thisnewlat,thisnewlon])
    elif slab == 'car':
        lattrench = thistrench[thistrench.lon > 300]
        lontrench = thistrench[thistrench.lon <= 300]
        lattrench = lattrench.sort(['lat'],ascending=False)
        lontrench = lontrench.sort(['lon'],ascending=True)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        thisnew = pd.concat([thisnewlon,thisnewlat])
    elif slab == 'cam':
        lattrench = thistrench[thistrench.lat > 18]
        lontrench = thistrench[thistrench.lat <= 18]
        lattrench = lattrench.sort(['lat'],ascending=True)
        lontrench = lontrench.sort(['lon'],ascending=False)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        thisnew = pd.concat([thisnewlon,thisnewlat])
    elif slab == 'hel':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'ita':
        lattrench1 = thistrench[thistrench.lat > 40]
        lattrench2 = thistrench[(thistrench.lat <= 40) & (thistrench.lon > 17)]
        lontrench = thistrench[(thistrench.lat <= 40) & (thistrench.lon <= 17)]
        lattrench = pd.concat([lattrench1,lattrench2])
        lattrench = lattrench.sort(['lat'],ascending=False)
        lontrench = lontrench.sort(['lon'],ascending=False)
        thisnewlat = newaz(lattrench)
        thisnewlon = newaz(lontrench)
        #thisnewlat = thisnewlat.reset_index(drop=True)
        #thisnewlon = thisnewlon.reset_index(drop=True)
        #mergelon1,mergelat1 = thisnewlat['lon'].values[-1],thisnewlat['lat'].values[-1]
        #mergelon2,mergelat2 = thisnewlon['lon'].values[-1],thisnewlon['lat'].values[-1]
        #dist, ang, lat1, lon1 = cosine(mergelon1,mergelat1,mergelon2,mergelat2)
        #thisnewlon.set_value(len(thisnewlon),'az',ang)
        thisnew = pd.concat([thisnewlat,thisnewlon])
    elif slab == 'mak':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'mue':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'cot':
        thistrench = thistrench.sort(['lat'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'sul':
        thistrench = thistrench.sort(['lon'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'pan':
        thistrench = thistrench.sort(['lon'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'him':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'hin':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'alu':
        thistrench = thistrench.sort(['lon'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'ker':
        thistrench = thistrench.sort(['lat'],ascending=False)
        thisnew = newaz(thistrench)
    elif slab == 'png':
        thistrench = thistrench.sort(['lon'],ascending=True)
        thisnew = newaz(thistrench)
    elif slab == 'puy':
        thistrench = thistrench.sort(['lat'],ascending=True)
        thisnew = newaz(thistrench)
    else:
        thisnew = newaz(thistrench)
    thisnew['number'] = range(len(thisnew))
    thisnew = thisnew.sort(['number'],ascending=False)
    allslabs = pd.concat([allslabs,thisnew])

allslabs2 = pd.DataFrame()

for slab in slablist:
    thistrench = allslabs[oldtrench.slab == slab]
    thisnew = newaz(thistrench)
    thisnew['number'] = range(len(thisnew))
    allslabs2 = pd.concat([allslabs2,thisnew])

allslabs2 = oneeighty(allslabs)
allslabs2 = allslabs2[['lon','lat','az','bound','slab','number']]
allslabs2.to_csv('ryutrench2.csv',header=True,index=False,na_rep=np.nan)











