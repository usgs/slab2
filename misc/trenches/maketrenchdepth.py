import pandas as pd
import numpy as np
from neicmap.distance import sdist, getAzimuth
import os

trenches = pd.read_csv('trenches_usgs_2017.csv')
trenches.loc[trenches.lon < 0, 'lon'] += 360
masterdir = '../../MasterDB/bathymetry'

BAdataset = pd.DataFrame()
slablist = []
for SGfile in os.listdir(masterdir):
    try:
        slab = SGfile[0:3]
        badata = pd.read_csv('%s/%s'%(masterdir,SGfile))
        badata['slab'] = slab
        BAdataset = pd.concat([BAdataset,badata])
        slablist.append(slab)
    except:
        print ('couldnt read this file',SGfile)

newtrdata = pd.DataFrame()
for slab in slablist:

    BAdata = BAdataset[BAdataset.slab == slab]
    TRdata = trenches[trenches.slab == slab]
    
    depthlist = []
    for index,row in TRdata.iterrows():
    
        tlon,tlat = row['lon'],row['lat']
        locba = BAdata[(BAdata.lon<tlon+0.5)&(BAdata.lon>tlon-0.5) &\
                        (BAdata.lat<tlat+0.5)&(BAdata.lat>tlat-0.5)]
        if len(locba)>0:
            locba['dist'] = sdist(tlat, tlon, locba['lat'], locba['lon'])
            mindist = locba['dist'].min()
            thisba = locba[locba.dist == mindist]
            if slab == 'him':
                depthlist.append(thisba['elev'].values[0]/1000)
            else:
                depthlist.append(thisba['depth'].values[0])
        else:
            depthlist.append(np.nan)

        print (slab, tlon, tlat)
    TRdata['depth'] = depthlist


    newtrdata = pd.concat([newtrdata,TRdata])

newtrdata.to_csv('trenches_usgs_2017_depths.csv',header=True,index=False,na_rep=np.nan)

