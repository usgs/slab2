import os
import csv
from datetime import datetime
import numpy as np
import pandas as pd
from libcomcat.search import search
from libcomcat.dataframes import get_summary_data_frame # KLH 09/24/2019
from libcomcat.dataframes import get_detail_data_frame # KLH 09/24/2019
from libcomcat.search import get_event_by_id
import argparse

def main(args):
    
    failed = args.failed
    queried = args.queried

    searchlist = []
    idlist = []
    badlist = []
    for filename in os.listdir(failed):
        print ('searching through terminal output %s'%filename)
        with open('%s/%s'%(failed,filename)) as fp:
            for cnt, line in enumerate(fp):
                if line[0:6] == 'Failed':
                    f,t,g,d,v,o,e,id_no = line.split()
                    try:
                        searchobj = get_event_by_id(id_no)
                        searchlist.append(searchobj)
                        idlist.append(id_no)
                    except:
                        badlist.append(id_no)

    print ('getting details for %i events: '%(len(idlist)),idlist)
    print ('these events %i still failed: '%(len(badlist)),badlist)
    detaildf = get_summary_data_frame(searchlist)
    detaildf.to_csv('%s/%sfailed.csv'%(queried,queried),header=True,index=False,na_rep=np.nan)

if __name__=='__main__':
    desc = '''
        This is used to go through terminal outputs and find failed PDE queries,
        then re-query the catalog for those IDs.
        
        Run this by typing:
            python getfailedqueries.py -f faileddirectory -q querydirectory
            
        where faileddirectory is where the terminal windows were saved and 
        querydirectory is the directory where the rest of the query files were
        saved.
                
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-f', '--failed', dest='failed', type=str,
                        required=True,help='directory where terminal outputs from original query are stored')
    parser.add_argument('-q', '--queried', dest='queried', type=str,
                        required=True,help='directory where all of the query output files were stored')
    
    pargs = parser.parse_args()
    
    main(pargs)
