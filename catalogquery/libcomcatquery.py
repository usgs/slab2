#!/usr/bin/env python

#imports
import os
import csv
from datetime import datetime
import numpy as np
import pandas as pd
from libcomcat.search import search
from libcomcat.utils import get_summary_data_frame
from libcomcat.utils import get_detail_data_frame
import argparse

def main(args):
    
    intervalno = args.intervalno
    directory = args.directory
    if intervalno > 11:
        if args.previous is not None:
            previous = args.previous
        else:
            print ('previous file required at input if interval no > 11.')
            print ('enter previous file with flag -p previousFile at input')
            print ('Exiting ...')
            exit()

    os.system('mkdir %s'%directory)

    if intervalno == 1:
        start = datetime(1900, 1, 1)
        finish = datetime(1970, 1, 1)
    if intervalno == 2:
        start = datetime(1970, 1, 1)
        finish = datetime(1980, 1, 1)
    if intervalno == 3:
        start = datetime(1980, 1, 1)
        finish = datetime(1990, 1, 1)
    if intervalno == 4:
        start = datetime(1990, 1, 1)
        finish = datetime(1995, 1, 1)
    if intervalno == 5:
        start = datetime(1995, 1, 1)
        finish = datetime(2000, 1, 1)
    if intervalno == 6:
        start = datetime(2000, 1, 1)
        finish = datetime(2005, 1, 1)
    if intervalno == 7:
        start = datetime(2005, 1, 1)
        finish = datetime(2010, 1, 1)
    if intervalno == 8:
        start = datetime(2010, 1, 1)
        finish = datetime(2012, 1, 1)
    if intervalno == 9:
        start = datetime(2012, 1, 1)
        finish = datetime(2014, 1, 1)
    if intervalno == 10:
        start = datetime(2014, 1, 1)
        finish = datetime(2016, 1, 1)
    if intervalno == 11:
        start = datetime(2016, 1, 1)
        finish = datetime.utcnow()
    if intervalno > 11:
        predata = pd.read_csv(previous)
        predata['time'] = pd.to_datetime(predata['time'])
        start = predata['time'].max()
        finish = datetime.utcnow()

    print (start)
    print (finish)

    #define magnitude range of earthquakes to search over for shallow earthquakes
    #(these are whatever magnitude range is defined in the catalogues)
    min_mag = 3.0
    max_mag = 9.9
    magrange = (min_mag, max_mag)

    #define depth to search over
    min_sh = 0
    max_sh = 900
    depthrange_sh = (min_sh, max_sh)

    #define grid size and extent of search (representative of lower left corner)
    if intervalno > 11:
        grid = 50.0
    else:
        grid = 10.0
    lonmin, lonmax = -180, 180
    latmin, latmax = -75, 75

    #define grid of searches (representative of lower left corner)
    xall = np.arange(lonmin,lonmax,grid)
    yall = np.arange(latmin,latmax,grid)
    lons1,lats1 = np.meshgrid(xall,yall)

    #flatten into list of lower left corners
    lllons = lons1.flatten()
    lllats = lats1.flatten()

    #define lists representing upper right corners
    urlons = lllons+grid
    urlats = lllats+grid

    #combine into one array (lonmin,lonmax,latmin,latmax)
    bounds = np.zeros((len(lllons),4))
    bounds[:,0] = lllons
    bounds[:,1] = urlons
    bounds[:,2] = lllats
    bounds[:,3] = urlats
    iterations = len(bounds)

    lllatfail = []
    lllonfail = []
    urlatfail = []
    urlonfail = []

    totdf = pd.DataFrame()
    num = 0
    for i, line in enumerate(bounds):

        bounds = line

        #to follow along with the progress of data querying
        #since there are a lot of iterations, only print every so often
        k = 100
        if i % k == 0:
            
            print ('Now querying grid %s of %s' % (i, iterations))
        
        searchlist = search(starttime=start, endtime=finish, minlatitude=bounds[2], maxlatitude=bounds[3], minlongitude=bounds[0], maxlongitude=bounds[1],  minmagnitude=3.0)

        if len(searchlist) > 0:
            detaildf = get_detail_data_frame(searchlist,get_tensors='preferred',get_moment_supplement=True)

            totdf = pd.concat([totdf,detaildf])
            print (bounds,len(detaildf),len(totdf))

            if len(totdf) > 5000:
                totdf.to_csv('%s/%s_%i.%i.csv'%(directory,directory,intervalno,num),header=True,index=False,na_rep=np.nan)
                num += 1
                totdf = pd.DataFrame()


    totdf.to_csv('%s/%s_%i.csv'%(directory,directory,intervalno),header=True,index=False,na_rep=np.nan)

# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        This is used to query the PDE either:
        
            in its entirety with 11 simultaneous runs or
            
            since the most recent event from a previous query
            
        To query the catalog in its entirety:
        
            open 11 terminal windows
            
            run each window with a different -i flag indicating different 
            date ranges:
                -i 1    >   1/1/1900 to 1/1/1970
                -i 2    >   1/1/1970 to 1/1/1980
                -i 3    >   1/1/1980 to 1/1/1990
                -i 4    >   1/1/1990 to 1/1/1995
                -i 5    >   1/1/1995 to 1/1/2000
                -i 6    >   1/1/2000 to 1/1/2005
                -i 7    >   1/1/2005 to 1/1/2010
                -i 8    >   1/1/2010 to 1/1/2012
                -i 9    >   1/1/2012 to 1/1/2014
                -i 10   >   1/1/2014 to 1/1/2016
                -i 11   >   1/1/2016 to present
            the script will take about 10 hours to complete for all date ranges
            
            files will be written to the -d directory input with different number 
            extensions
            
        To query the catalog since the most recent event in a previous query:
        
            run the script with -i flag > 11
            
            include -p previous file in argument list when running
            
            a file containing events since most recent event in the listed file 
            will be written to the -d directory
                
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-p', '--previous', dest='previous', type=str,
                        help='previous query file, required if -i flag is greater than 11')
    parser.add_argument('-d', '--directory', dest='directory', type=str,
                        required=True,help='directory to save file to')
    parser.add_argument('-i', '--intervalno', dest='intervalno', type=int,
                        required=True,help='number associated with a date range to search over')
    
    pargs = parser.parse_args()
    
    main(pargs)
