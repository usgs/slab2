#!/usr/bin/env python

#stdlib
import argparse
from datetime import datetime,timedelta
import os
import sys
import re
import pandas as pd
import numpy as np

#local library
from s2d_fnctns import *
#from slab2functions import * 
#MF 8.4.16 need to import s2f as well? or copy to s2d?

def main(args):

    '''What is executed upon running the program. 
        
        Runs through the provided input files and calls writetofile (above) for each. 
        The result is a new comma delimited file with information pertinent to Slab 2.0.  '''
    
    slab_ = args.outFile
    os.system("rm %s"%args.outFile)
    default_uncertainty = 15 #km
    seismo_thick = args.seis
    
    filelist = []
    
    # If the database is called, compiles all files within bounds into a list to be written to
    # output file.
    if args.database is not None:
        if args.bounds is not None:
            check = 1
            lonmin = args.bounds[0]
            lonmax = args.bounds[1]
            latmin = args.bounds[2]
            latmax = args.bounds[3]
            slab = rectangleIntersectsPolygon(lonmin,lonmax,latmin,latmax)
        else:
            check = 0
            slab = args.slab_polygon
        for filename in os.listdir(args.database):
            if filename.endswith('.csv'):
                slabname,etype,name = filename.split('_')
                if slabname == 'ALL':
                    slabname = args.slab_polygon
                if etype == 'SR' or etype == 'ST':
                    etype = 'AS'
                name = name[:-4]
                if slabname == slab:
                    # Changes the default uncertainty based on event type.
                    if etype == 'AS':
                        default_uncertainty = 2.5
                    elif etype == 'RF':
                        default_uncertainty = 10
                    elif etype == 'ER':
                        default_uncertainty = 10
                    elif etype == 'GP':
                        default_uncertainty = 10
                    elif etype == 'BA':
                        default_uncertainty = 1
                    elif etype == 'TO':
                        default_uncertainty = 40
                    elif etype == 'CP':
                        default_uncertainty = 10
                    f = NewFile(args.database+'/'+filename, default_uncertainty, etype, name)
                    filelist.append(f)
                else:
                    pass
            else:
                print ('The file %s was not written to the dataset - FILES MUST BE IN .CSV FORMAT' % filename)

    #if the database is called and regional preference is specified
    if args.database_r is not None:
        files = os.listdir(args.database_r)
        regional_pref = 1
        for filename in files:
            if filename.endswith('.csv'):
                slabname,etype,name = filename.split('_')
                name = name[:-4]
                # Changes the default uncertainty based on event type.
                if etype == 'AS':
                    default_uncertainty = 2.5
                elif etype == 'RF':
                    default_uncertainty = 10
                elif etype == 'GP':
                    default_uncertainty = 10
                elif etype == 'ER':
                    default_uncertainty = 10
                elif etype == 'BA':
                    default_uncertainty = 1
                elif etype == 'TO':
                    default_uncertainty = 40
                elif etype == 'CP':
                    default_uncertainty = 10
                # If bounds are provided, only adds files labeled with the associated slab region to file list.
                if args.bounds is not None:
                    check = 1
                    lonmin = args.bounds[0]
                    lonmax = args.bounds[1]
                    latmin = args.bounds[2]
                    latmax = args.bounds[3]
                    slabname = rectangleIntersectsPolygon(lonmin,lonmax,latmin,latmax)
                    if slabname in slablist:
                        # Creates file object to be written to output file
                        f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                        filelist.append(f)
                    #MF 8.4.16 this won't be an option
                        #elif slabname == 'ALL':
                        # Creates file object to be written to output file
                        #f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                        #filelist.append(f)
                    else:
                        pass
                else:
                    check = 0
                    slabname = args.slab_polygon
                    # Creates file object to be written to output file
                    f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                    filelist.append(f)
            else:
                print ('The file %s was not written to the dataset - FILES MUST BE IN .CSV FORMAT' % filename)

    #if the database is called and global preference is specified
    if args.database_g is not None:
        files = os.listdir(args.database_g)
        regional_pref = 0
        for filename in files:
            print (filename)
            if filename.endswith('.csv'):
                slabname,etype,name = filename.split('_')
                name = name[:-4]
                # Changes the default uncertainty based on event type.
                if etype == 'AS':
                    default_uncertainty = 2.5
                elif etype == 'RF':
                    default_uncertainty = 10
                elif etype == 'GP':
                    default_uncertainty = 10
                elif etype == 'ER':
                    default_uncertainty = 10
                elif etype == 'BA':
                    default_uncertainty = 1
                elif etype == 'TO':
                    default_uncertainty = 40
                elif etype == 'CP':
                    default_uncertainty = 10
                # If bounds are provided, only adds files labeled with the associated slab region to file list.
                if args.bounds is not None:
                    check = 1
                    lonmin = args.bounds[0]
                    lonmax = args.bounds[1]
                    latmin = args.bounds[2]
                    latmax = args.bounds[3]
                    slabname = rectangleIntersectsPolygon(lonmin,lonmax,latmin,latmax)
                    if slabname in slablist:
                        # Creates file object to be written to output file
                        f = NewFile(args.database+'/'+filename, default_uncertainty, etype, name)
                        filelist.append(f)
                    #MF 8.4.16 this won't be an option
                        #elif slabname == 'ALL':
                        # Creates file object to be written to output file
                        #f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                        #filelist.append(f)
                    else:
                        pass
                else:
                    check = 0
                    slabname = args.slab_polygon
                    # Creates file object to be written to output file
                    f = NewFile(args.database+'/'+filename, default_uncertainty, etype, name)
                    filelist.append(f)
                print ('The file %s was not written to the dataset - FILES MUST BE IN .CSV FORMAT' % filename)

    # Adds additional input files to file list if they are provided in command line.
    if args.input is not None:
        for file in args.input:
            try:
                f=NewFile(file[0], file[1], file[2])
                filelist.append(f)
            except:
                try:
                    f = NewFile(file[0], default_uncertainty, file[1])
                    filelist.append(f)
                except:
                    print ('Input file information must be given as infile,unc,etype \
                        or as infile,etype')

    # Notifies user and exits if no files are written to the output file.
    if len(filelist) == 0:
        print ('No input files were provided. Include argument: -i infile1,unc1,etype1, -i \
        infile2,unc2,etype2 or include argument: -d slab2database_location.')
        sys.exit(0)

    catalogs = []
    file_no = 1
    #Writes each file in filelist (events in bounds within each file) to output file
    for file in filelist:
        writetofile(file.filename, slab_, file.event_type, file.unc, args, catalogs, file_no, seismo_thick, slab, file.name)
        file_no = file_no+2

    # Makes rough plot of data for the region
    # slabplotter(args)

#Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        Writes a file in csv format with the information pertinent to Slab 2.0. Columns
            specified as:
        
        (lat,lon,depth,uncertainty,event-type,mag,time,P-azimuth,P-plunge,T-azimuth,
            T-plunge,strike1,dip1,rake1,strike2,dip2,rake2,ID)
        
        Fields are represented by 'nan' where information is not available. Minimum input information
            includes columns labeled as lat,lon,depth. If the type of event is an earthquake, time 
            and mag are also required data columns.
        Where applicable, CMT information should be represented in tensorial form with 
            columns labeled: mrr,mtt,mpp,mrt,mrp and mtp.
        
        Expected event type input:
            Earthquake: EQ
            Receiver Function: RF
            Tomography: TO
            Active Source: AS
            Bathymetry: BA
            Earthquake Relocation: ER
            Centroid Moment Tensor: MT
            GPS: GP
            
        EXAMPLE: to compile slab and moment tensor information around Southern Peru in 2013 from
            original data files:
        
        s2d.py -b -85 -60 -25 -15 -s 2013-01-01 -e 2014-01-01 -f speru13slab.csv 
            -i isc-gem.csv,15,EQ -i csn_cat2.txt,15,EQ -i so_peru_rf.txt,10,RF 
            -i s_peru_to.tsv,10,TO
        
        EXAMPLE: to compile slab and moment tensor information around Southern Peru from all 
            available files in slab2database:
        
        s2d.py -b -85 -60 -25 -15 -d /some/directory/slab2database -f speru_slab.csv
        
        The database and original files can be added in the same call of s2d.py.
        
        Note that output file and lat/lon bounds are required arguments. If mag and time bounds 
            are not specified, the values are set to all magnitudes and dates from 1900 to present.
        
        Note that when specifying a search box that crosses the -180/180 meridian, specify 
            longitudes as you would if you were not crossing that meridian (i.e., lonmin=179, 
            lonmax=-179).  The program will resolve the discrepancy.
        
        If more than one earthquake dataset is provided, the files will be compared and matching 
            events will not be written to the slab file. The program prioritizes catalogs in 
            the order in which they are entered, writing the matching events from the first entry 
            and disregarding the associated match(es) in later entries.
        If more than one match is found for a single event, the program selects the closest match 
            and determines the remaining potential matches as independent.
        
        Where CMT information is provided, non-thrust earthquakes at depths shallower than 60 km 
            are filtered out.
        
        The file is saved, and more information can be appended to it by calling to s2d.py.
        
        A local copy of s2d_fnctns.py is required to run this program.
        
        '''
    
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group2 = parser.add_mutually_exclusive_group(required=True)
    #required arguments
    parser.add_argument('-f','--output-file', dest='outFile', type=str, required=True,
                        help='Designated output file where data is written.')
    #mutually exclusive input data arguments
    group.add_argument('-d','--database', dest='database', type=str, 
                       help = 'Directory where all current slab 2.0 files are located. If this option is chosen, and if multiple earthquake catalogues exist, they will NOT be compared for matching earthquakes.')
    group.add_argument('-i','--inputFiles', dest='input', type=infile, action='append',
                        help = 'List of input files with their associated uncertainty and event-type. \
                        -i input1,unc1,etype1 -i input2,unc2,etype2')
    group.add_argument('-dg','--database_g', dest='database_g', type=str, 
                       help = 'Directory where all current slab 2.0 files are located. If this option is chosen, and if multiple earthquake catalogues exist, they will be compared and preference will be given to global catalogues over regional catalogues if there exist matching earthquakes.')
    group.add_argument('-dr','--database_r', dest='database_r', type=str, 
                       help = 'Directory where all current slab 2.0 files are located. If this option is chosen, and if multiple earthquake catalogues exist, they will be compared and preference will be given to regional catalogues over global catalogues if there exist matching earthquakes.')
    #mutually exclusive input slab boundary arguments
    group2.add_argument('-b','--bounds', metavar=('lonmin','lonmax','latmin','latmax'),
                        dest='bounds', type=float, nargs=4,
                        help='Rectangular bounds to constrain event search [lonmin lonmax latmin latmax]')
    #MF 8.3.16 edited input arguments
    group2.add_argument('-p','--slab_polygon', dest='slab_polygon', type=str,
                        help = 'Name of slab region to search. If this option is selected, default polygon boundaries are entered for the slab.')
    parser.add_argument('-t','--seis', dest='seis', type=float,
                        help = 'Depth at which earthquakes start occur within plate')
    #optional arguments
    parser.add_argument('-s','--start-time', dest='startTime', type=maketime,
                        help='Start time for search (defaults to 01/01/1900).  YYYY-mm-dd or \
                        YYYY-mm-ddTHH:MM:SS')
    parser.add_argument('-e','--end-time', dest='endTime', type=maketime,
                        help='End time for search (defaults to current date/time).  YYYY-mm-dd or \
                        YYYY-mm-ddTHH:MM:SS')
    parser.add_argument('-m','--mag-range', metavar=('minmag','maxmag'),dest='magRange',
                        type=float,nargs=2,
                        help='Min/max (authoritative) magnitude to restrict search.')

                        
    pargs = parser.parse_args()
                        
                        
    main(pargs)
