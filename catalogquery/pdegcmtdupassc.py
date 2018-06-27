#import standard libraries
import obspy.imaging.beachball
import datetime
import os
import csv
import pandas as pd
import numpy as np
import fnmatch
from geopy.distance import vincenty
from math import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import path
import argparse

class Earthquake:

    '''Creates an earthquake object from which event information can be extracted'''
    
    def __init__(self,time,coords,depth,lat,lon,mag,catalog,id):
        self.time = time
        self.coords = coords
        self.depth = depth
        self.lat = lat
        self.lon = lon
        self.mag = mag
        self.catalog = catalog
        self.id = id

def getvals(row):

    '''Gathers time, lat, lon, depth, mag, information from row in dataframe.'''

    time = row['time']
    lat = row['lat']
    lon = row['lon']
    depth = row['depth']
    mag = row['mag']
    ep = (lat,lon)
    try:
        ID = row['CID']
    except:
        ID = row['GID']
    return time,ep,depth,lat,lon,mag,ID

def earthquake_string(eqo):

    ''' Puts earthquake information into a string to be written or printed
    
        Arguments:  eqo - earthquake object
        
        Returns:    eqos - a string of information stored in earthquake object input argument '''
    
    eqos = (str(eqo.lat) + ',' + str(eqo.lon) + ',' + str(eqo.depth) + ','
            + str(eqo.mag) + ',' + str(eqo.time) + ',' + eqo.catalog)
    return eqos

def find_closest(eqo, eqm1, eqm2):

    '''Determines which of two potential matches in one catalog is closer to an event in another.
        
        Arguments:  eqo - earthquake event in first catalog that matches two events in the second
                    eqm1 - the first event in the second catalog that matches eqo
                    eqm2 - the second event in the second catalog that matches eqo
                    
        Returns:    closest - the closest event weighting time first, then distance, then magnitude '''
    
    # Prints information to console to make user aware of more than one match
    print ('-------------------------------------- lat %s lon %s depth %s mag %s time \
        %s catlog' % (',',',',',',',',','))
    print ('There is more than one match for event: %s' % earthquake_string(eqo))
    print ('event1: %s' % earthquake_string(eqm1))
    print ('event2: %s' % earthquake_string(eqm2))
    
    # Gets distance between either event and the common match eqo
    darc1 = vincenty(eqo.coords, eqm1.coords).meters/1000
    darc2 = vincenty(eqo.coords, eqm2.coords).meters/1000
    dh1 = abs(eqo.depth - eqm1.depth)
    dh2 = abs(eqo.depth - eqm2.depth)
    dist1 = sqrt(darc1*darc1 + dh1*dh1)
    dist2 = sqrt(darc2*darc2 + dh2*dh2)
    
    # Gets magnitude and time differences between each event and the common match
    dtime1 = abs(eqo.time - eqm1.time)
    dtime2 = abs(eqo.time - eqm2.time)
    dmag1 = abs(eqo.mag - eqm1.mag)
    dmag2 = abs(eqo.mag - eqm2.mag)
    
    # Finds the closest match to eqo by checking time first, then distance, then magnitude
    if dtime1 < dtime2:
        closest = eqm1
    elif dtime2 < dtime1:
        closest = eqm2
    elif dtime1 == dtime2 and dist1 < dist2:
        closest = eqm1
    elif dtime1 == dtime2 and dist2 < dist1:
        closest = eqm1
    elif dmag1 == dmag2 and dist1 == dist2 and dmag1 < dmag2:
        closest = eqm1
    elif dmag1 == dmag2 and dist1 == dist2 and dmag2 < dmag1:
        closest  = eqm2

    # If all things are equal, the first event is chosen as a match by default
    else:
        print ('The two events are equidistant to the match in time, space, and magnitude.\
            The second event was therefore determined independent.')
        closest = eqm1
        return closest
    print ('>>>>closest event: %s' % earthquake_string(closest))
    return closest

def removematches(dfo, dfm):

    '''Eliminates events in dfo (dataframe) that are found in dfm (dataframe) '''
    
    ind = (dfo.time.isin(dfm.time) & dfo.lat.isin(dfm.lat) & dfo.lon.isin(dfm.lon)
           & dfo.mag.isin(dfm.mag) & dfo.depth.isin(dfm.depth))
    dfo = dfo[~ind]
    return dfo

def removematchesM(dfo, dfm):

    '''Eliminates events in dfo (dataframe) that are found in dfm (dataframe) '''
    
    ind = (dfo.mrr.isin(dfm.mrr) & dfo.mtt.isin(dfm.mtt) & dfo.mpp.isin(dfm.mpp)
           & dfo.mrt.isin(dfm.mrt) & dfo.mrp.isin(dfm.mrp) & dfo.mtp.isin(dfm.mtp))
    dfo = dfo[~ind]
    return dfo

def rid_matches(cat1c, cat2c, name1, name2, whichcat):

    ''' Compares two catalogs, identifies and removes matching events from cat2.
        
        Arguments:  cat1 - the first catalog (dataframe), no events are removed from this catalog
                    cat2 - the second catalog (dataframe), events in this catalog that are close 
                            in space, time, and magnitude to those in cat1 are filtered out
                    
        Returns:    df - a filtered version of cat2 without events that match those in cat1 '''
    
    # Setting constants that define matching criteria
    tdelta = 30
    distdelta = 100
    magdelta = 0.5
    tdeltacomp = datetime.timedelta(seconds = tdelta)
    
    # Ensuring that all times are in datetime object format & trimming catalogs to only extend
    # accross the bounds and time constraints of the other
    cat1c['time'] = pd.to_datetime(cat1c['time'])
    cat2c['time'] = pd.to_datetime(cat2c['time'])

    matches = pd.DataFrame(columns = ['clat','clon','cdepth','cmag','ctime','glat','glon','gdepth','gmag','gtime','id_no','type'])
    matches1 = pd.DataFrame(columns = ['lon','lat','depth','mag','time'])

    count = 0

    # Compares each event in cat2 to each event in cat1
    for index,row in cat1c.iterrows():
        n = 0
        
        print (index, len(cat1c))
        # Getting earthquake info from event and storing it in an Earthquake object (cat1)
        time1, ep1, depth1, lat1, lon1, mag1, ID1 = getvals(row)
        eq1 = Earthquake(time1, ep1, depth1, lat1, lon1, mag1, name1, ID1)
        cat2c2 = cat2c[(cat2c.lon < lon1 + 1)&(cat2c.lon > lon1 - 1)& \
                           (cat2c.lat < lat1 + 1)&(cat2c.lon > lat1 - 1)]
        cat2c2 = cat2c2[(cat2c2.mag < mag1+magdelta)&(cat2c2.mag > mag1-magdelta)]
        cat2c2['timediff'] = np.abs((cat2c2.time-time1).astype('timedelta64[s]'))
        cat2c2 = cat2c2[cat2c2.timediff < tdelta]
        for i, r in cat2c2.iterrows():
        
            # Getting earthquake info from event and storing it in an Earthquake object (cat1)
            time2, ep2, depth2, lat2, lon2, mag2, ID2 = getvals(r)
            eq2 = Earthquake(time2, ep2, depth2, lat2, lon2, mag2, name2, ID2)

            if vincenty(ep1,ep2).meters/1000 <= distdelta:
                
                # If there is already a match for this event, find the closest
                # The closest is stored and compared to third, fourth matches etc if they exist
                if n >= 1:
                    match = find_closest(eq1, match, eq2)
                    n += 1
                if n == 0:
                    match = eq2
                    n += 1
            else:
                pass
    
        
        # Add matching events to match dataframe
        if n > 0:
            if name2 == 'C':
                thisevent = cat2c[cat2c.time == match.time]
                id_no = thisevent['id_no'].values[0]
                mttype = thisevent['type'].values[0]
                matches.loc[len(matches)+1] = [match.lat, match.lon, match.depth, match.mag, match.time,lat1, lon1, depth1, mag1, time1,id_no,mttype]

            elif name2 == 'G':
                thisevent = cat1c[cat1c.time == time1]
                id_no = thisevent['id_no'].values[0]
                mttype = thisevent['type'].values[0]
                matches.loc[len(matches)+1] = [lat1, lon1, depth1, mag1, time1,match.lat, match.lon, match.depth, match.mag, match.time,id_no,mttype]
                
            else:
                if n > 1:
                    matches.loc[len(matches)+1] = [lat1, lon1, depth1, mag1, time1, eq1.lat, eq1.lon, eq1.depth, eq1.mag, eq1.time,count,count]
                    matches.loc[len(matches)+1] = [lat1, lon1, depth1, mag1, time1, eq2.lat, eq2.lon, eq2.depth, eq2.mag, eq2.time,count,count]
            
            matches1.loc[len(matches1)+1] = [match.lon,match.lat,match.depth,match.mag,match.time]
            count += 1
        
    # Write matches to matching file
    matchfile = 'matchfiles/cg_matches_%s-%s.csv'%(name1,name2)
    matches.to_csv(matchfile,header=True,index=False,na_rep=np.nan)
    
    if whichcat == 2:
        df = removematches(cat2c, matches1)
    else:
        df = removematches(cat1c, matches1)

    # Return filtered catalog to be written to output file (df), and list of matches from second catalog (matches)
    return df, matches

def savedupes(cat11, cat22, matches1, filename):

    ''' Arguments:  cat1 - catalog in initial format (i.e. lon, lat, depth etc)
                    cat2 - catalog in initial format (i.e. lon, lat, depth etc)
                    matches - df listing matches (i.e. clon, clat ... glon, glat)
                        assumes cat1 = clon, clat ... cat2 = glon, glat ...
                    filename - cat1name_cat2name.csv
                
        Returns:    writes dataframe with all info from both catalogs '''
    
    
    # make copies of initial dataframes so originals are not altered
    cat1 = cat11.copy()
    cat2 = cat22.copy()
    matches = matches1.copy()
    
    # rename all columns for merge (first catalog)
    cat1['lon1'] = cat1['lon'].values*1.0
    cat1['lat1'] = cat1['lat'].values*1.0
    cat1['depth1'] = cat1['depth'].values*1.0
    cat1['mag1'] = cat1['mag'].values*1.0
    cat1['time1'] = cat1['time'].values
    cat1['mrr1'] = cat1['mrr'].values*1.0
    cat1['mtt1'] = cat1['mtt'].values*1.0
    cat1['mpp1'] = cat1['mpp'].values*1.0
    cat1['mrt1'] = cat1['mrt'].values*1.0
    cat1['mrp1'] = cat1['mrp'].values*1.0
    cat1['mtp1'] = cat1['mtp'].values*1.0
    cat1['type1'] = cat1['type'].values
    
    # rename all columns for merge (second catalog)
    cat2['lon2'] = cat2['lon'].values*1.0
    cat2['lat2'] = cat2['lat'].values*1.0
    cat2['depth2'] = cat2['depth'].values*1.0
    cat2['mag2'] = cat2['mag'].values*1.0
    cat2['time2'] = cat2['time'].values
    cat2['mrr2'] = cat2['mrr'].values*1.0
    cat2['mtt2'] = cat2['mtt'].values*1.0
    cat2['mpp2'] = cat2['mpp'].values*1.0
    cat2['mrt2'] = cat2['mrt'].values*1.0
    cat2['mrp2'] = cat2['mrp'].values*1.0
    cat2['mtp2'] = cat2['mtp'].values*1.0
    cat2['type2'] = cat2['type'].values

    # merge match file with first and second catalogs
    merge1 = pd.merge(matches, cat1, left_on = ['clat','clon','cmag','cdepth','ctime'], right_on = ['lat1','lon1','mag1','depth1','time1'])
    merge2 = pd.merge(matches, cat2, left_on = ['glat','glon','gmag','gdepth','gtime'], right_on = ['lat2','lon2','mag2','depth2','time2'])
    
    # merge first and second catalogs based on columns in match file
    merge3 = pd.merge(merge1, merge2, left_on = ['glat','glon','gmag','gdepth','gtime'], right_on = ['glat','glon','gmag','gdepth','gtime'])

    # trim merged file to list matched information
    merge3 = merge3[['lat1','lon1','mag1','depth1','time1','mrr1','mtt1','mpp1','mrt1','mrp1','mtp1','type1','lat2','lon2','mag2','depth2','time2','mrr2','mtt2','mpp2','mrt2','mrp2','mtp2','type2']]

    merge3.to_csv(filename, header=True,index=False, na_rep=np.nan)

def main(args):
    # Import gCMT and the latest query of the PDE
    gcmtfile = args.gcmtfile
    comcatfile = args.ccatfile
    associatedfile = args.assofile

    gcmt1 = pd.read_csv(gcmtfile,delim_whitespace=True)
    comcat1 = pd.read_csv(comcatfile)

    comcat1['mlon'] = comcat1['moment_lon'].values*1.0
    comcat1['mlat'] = comcat1['moment_lat'].values*1.0
    comcat1['mdep'] = comcat1['moment_depth'].values*1.0

    # Removing any duplicate IDs from the comcat database (some search overlaps may exist)
    ogindex = list(comcat1.index)
    comcat = comcat1.drop_duplicates(subset=['id_no','time','lat','lon','depth'], keep="first")
    newindex = list(comcat.index)
    duprows = list(set(ogindex) - set(newindex))
    gdupe = comcat1.iloc[duprows]

    # for identifying duplicate events in the gcmt catalog
    ogindex = list(gcmt1.index)
    gcmt = gcmt1.drop_duplicates(subset=['year','month','day','hour','minute','second','lat','lon','depth'], keep="first")
    newindex = list(gcmt.index)
    duprows = list(set(ogindex) - set(newindex))
    gdupe = gcmt1.iloc[duprows]

    # Convert gCMT format to match that returned from getcsv
    gcmt['time'] = pd.to_datetime(gcmt[['year', 'month', 'day', 'hour', 'minute', 'second']])
    comcat['time'] = pd.to_datetime(comcat['time'])

    # Convert units of gCMT tensors to match with PDE (Nm)
    exponents = np.power(10,gcmt['expo'].values-7)
    mrr = gcmt['mrr'].values.astype(float)
    mtt = gcmt['mtt'].values.astype(float)
    mpp = gcmt['mpp'].values.astype(float)
    mrt = gcmt['mrt'].values.astype(float)
    mrp = gcmt['mrp'].values.astype(float)
    mtp = gcmt['mtp'].values.astype(float)
    gcmt['mrr'] = mrr*exponents
    gcmt['mtt'] = mtt*exponents
    gcmt['mpp'] = mpp*exponents
    gcmt['mrt'] = mrt*exponents
    gcmt['mrp'] = mrp*exponents
    gcmt['mtp'] = mtp*exponents
    gcmt['type'] = 'gcmt'

    # Make sure all events follow the same longitudinal coordinate scheme
    gcmt['lon'][gcmt.lon<0]+=360
    comcat['lon'][comcat.lon<0]+=360

    # Find 1990 double event (special case that doesn't get associated)
    first1990g = gcmt[(gcmt.year == 1990)&(gcmt.month == 4)&(gcmt.day == 26) & \
                    (gcmt.hour == 9)&(gcmt.minute == 37)&(gcmt.second == 14.30) & \
                    (gcmt.lat == 36.010)&(gcmt.lon == 100.270) & \
                    (gcmt.depth == 33.0)&(gcmt.mag == 6.25)]

    second1990g = gcmt[(gcmt.year == 1990)&(gcmt.month == 4)&(gcmt.day == 26) & \
                    (gcmt.hour == 9)&(gcmt.minute == 37)&(gcmt.second == 18.50) & \
                    (gcmt.lat == 35.960)&(gcmt.lon == 100.230) & \
                    (gcmt.depth == 33.0)&(gcmt.mag == 6.43)]

    # Find double event in comcat
    first1990c = comcat[comcat.id_no == 'usp00048gt']
    second1990c = comcat[comcat.id_no == 'usp00048gu']

    # Reset indices
    first1990g = first1990g.reset_index(drop=True)
    second1990g = second1990g.reset_index(drop=True)
    first1990c = first1990c.reset_index(drop=True)
    second1990c = second1990c.reset_index(drop=True)

    # Merge events in correct order
    first1990g = first1990g[['mrr','mtt','mpp','mrp','mtp','mrt','gmlon','gmlat','gmdep']]
    second1990g = second1990g[['mrr','mtt','mpp','mrp','mtp','mrt','gmlon','gmlat','gmdep']]
    first1990c = first1990c[['id_no','time','lat','lon','depth','mag']]
    second1990c = second1990c[['id_no','time','lat','lon','depth','mag']]
    first1990 = pd.concat([first1990c, first1990g], axis=1)
    second1990 = pd.concat([second1990c, second1990g], axis=1)

    # Save double event association to be added to final dataset
    wonkychina = pd.concat([first1990,second1990])
    wonkychina['type'] = 'gcmt'
    wonkychina['mlon'] = wonkychina['gmlon'].values.astype(float)
    wonkychina['mlat'] = wonkychina['gmlat'].values.astype(float)
    wonkychina['mdep'] = wonkychina['gmdep'].values.astype(float)

    # Remove events from gCMT
    gcmt = gcmt[(gcmt.year != 1990)|(gcmt.month != 4)|(gcmt.day != 26) | \
                    (gcmt.hour != 9)|(gcmt.minute != 37)|(gcmt.second != 14.30) | \
                    (gcmt.lat != 36.010)|(gcmt.lon != 100.270) | \
                    (gcmt.depth != 33.0)|(gcmt.mag != 6.25)]

    gcmt = gcmt[(gcmt.year != 1990)|(gcmt.month != 4)|(gcmt.day != 26) | \
                    (gcmt.hour != 9)|(gcmt.minute != 37)|(gcmt.second != 18.50) | \
                    (gcmt.lat != 35.960)|(gcmt.lon != 100.230) | \
                    (gcmt.depth != 33.0)|(gcmt.mag != 6.43)]

    # Remove events from ComCat
    comcat = comcat[comcat.id_no != 'usp00048gt']
    comcat = comcat[comcat.id_no != 'usp00048gu']

    # trim for testing
    '''comcat = comcat[(comcat.lon < 121)&(comcat.lon>118)&(comcat.lat < -6)&(comcat.lat>-9)]
    gcmt = gcmt[(gcmt.lon < 121)&(gcmt.lon>118)&(gcmt.lat < -6)&(gcmt.lat>-9)]'''

    # Reset indices for looping, drop duplicate ids incase there are doubles
    comcat = comcat.reset_index(drop=True)
    gcmt = gcmt.reset_index(drop=True)
    comcat['CID'] = range(len(comcat))
    comcat['GID'] = range(len(comcat))
    gcmt['CID'] = range(len(gcmt))
    gcmt['GID'] = range(len(gcmt))

    # Identify events in PDE where gCMT moment tensor information will be used
    comcat1 = comcat[(np.isnan(comcat.mrr))|((comcat.type != 'usMww')&(comcat.type != 'duputelMww'))]

    # Identify events in PDE where existing MT info will be used
    comcat0 = comcat[(np.isfinite(comcat.mrr))&((comcat.type == 'usMww')|(comcat.type == 'duputelMww'))]
    comcat1 = comcat1.reset_index(drop=True)
    comcat0 = comcat0.reset_index(drop=True)

    # Associate gCMT and PDE events with MT type usMww or duputelMww (keep usMww and duputelMww events)
    Gcmt_notusMww, usMww_gcmt = rid_matches(comcat0,gcmt,'C','G',2)
    savedupes(comcat0, gcmt, usMww_gcmt, 'matchfiles/usMww_gcmt.csv')

    # Associate all other PDE events with gCMT (use PDE origin, gCMT MTs)
    comcat_useGcmt, allPDE_gcmt = rid_matches(Gcmt_notusMww,comcat1,'G','C',2)
    savedupes(comcat1, gcmt, allPDE_gcmt, 'matchfiles/allPDE_gcmt.csv')

    # Merge gCMT moment tensors with PDE origins/id_nos
    gcmt_time = Gcmt_notusMww[['mrr','mtt','mpp','mrt','mrp','mtp','time','lat','lon','mag','depth','gmlon','gmlat','gmdep']]
    use_gcmt_MT = pd.merge(allPDE_gcmt, gcmt_time, left_on = ['glat','glon','gmag','gdepth','gtime'], right_on = ['lat','lon','mag','depth','time'])
    use_gcmt_MT = use_gcmt_MT[['id_no','ctime','clat','clon','cdepth','cmag','mrr','mtt','mpp','mrt','mrp','mtp','type','gmlon','gmlat','gmdep']]
    renamed_gcmt = pd.DataFrame({'lon':use_gcmt_MT.clon, 'lat':use_gcmt_MT.clat, \
                            'depth':use_gcmt_MT.cdepth, 'mag':use_gcmt_MT.cmag, \
                            'time':use_gcmt_MT.ctime, 'mrr':use_gcmt_MT.mrr, \
                            'mtt':use_gcmt_MT.mtt, 'mpp':use_gcmt_MT.mpp, \
                            'mrt':use_gcmt_MT.mrt, 'mrp':use_gcmt_MT.mrp, \
                            'mtp':use_gcmt_MT.mtp, 'id_no':use_gcmt_MT.id_no, \
                            'mlon':use_gcmt_MT.gmlon, 'mlat':use_gcmt_MT.gmlat, \
                            'mdep':use_gcmt_MT.gmdep, 'type':'gcmt'})


    # Write all gcmt events not matched with comcat. likely matched with usMww event
    missing_gcmt = removematchesM(Gcmt_notusMww, renamed_gcmt)
    missing_gcmt.to_csv('matchfiles/gcmtNotinPDE.csv',header=True,index=False,na_rep=np.nan)


    # Merge all files and write final associated catalog to file
    final = pd.concat([comcat0,comcat_useGcmt,renamed_gcmt,wonkychina])
    final = final[['id_no','time','lat','lon','depth','mag','mrr','mtt','mpp','mrt','mrp','mtp','type','mlon','mlat','mdep']]
    final.to_csv(associatedfile,header=True,index=False,na_rep=np.nan)

    print ('Associated file was written to',associatedfile)

# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        This is used to associate a gcmt catalog file with a PDE query
        
        It retains all PDE earthquake origins
        
        Where the PDE listing does not include moment tensor information, the 
        moment tensor listing from GCMT is used
        
        Where the PDE listing has a moment tensor without type usMww or duputelMww,
        the moment tensor listing from GCMT is used
        
        Events listed in gcmt that are not listed in the PDE are not associated
        or written to the file
        
        Events that are listed in the PDE but not listed in gcmt are written to 
        the associated file without any moment tensor information
                
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-g', '--gcmtfile', dest='gcmtfile', type=str,
                        required=True,help='file name for file containing latest GCMT catalog query')
    parser.add_argument('-c', '--ccatfile', dest='ccatfile', type=str,
                        required=True,help='file name for file containing latest PDE catalog query')
    parser.add_argument('-f', '--assofile', dest='assofile', type=str,
                        required=True,help='name of file to save association to')
    
    pargs = parser.parse_args()
    
    main(pargs)
