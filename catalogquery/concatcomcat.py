import pandas as pd
import numpy as np
from datetime import datetime
import os

pd.options.mode.chained_assignment = None

manyortwo = input('If merging several files listed in one folder from a PDE query enter 1. If converting original comcat output to condensed slab2 format enter 1. If you are adding recent events to an existing catalog (both in condensed slab2 format) enter 2. : ')

if manyortwo == '1':
    #for combining an entire query
    database = input('enter database with several PDE query files to be merged:')

    dftot = pd.DataFrame()
    for filename in os.listdir(database):
        if filename.endswith('.csv'):
            df = pd.read_csv('%s/%s'%(database,filename))
            dftot = pd.concat([dftot,df],sort=True)
            print (len(df), len(dftot))

    places, mtypes = [],[]
    for name in dftot.columns:
        try:
            (place,mtype,der,latlond) = name.split('_')
            if latlond == 'depth':
                places.append(place)
                mtypes.append(mtype)
                colstr = '%s_%s'
        except:
            continue

    print ('places',places)
    print ('mtypes',mtypes)

    dftot['lon'] = dftot['longitude'].values*1.0
    dftot['lat'] = dftot['latitude'].values*1.0
    dftot['mag'] = dftot['magnitude'].values*1.0
    dftot['mag_type'] = dftot['magtype'].values
    dftot['id_no'] = dftot['id'].values

    newdf = pd.DataFrame()
    for i in range(len(places)):
        mrrstr = '%s_%s_mrr'%(places[i],mtypes[i])
        mttstr = '%s_%s_mtt'%(places[i],mtypes[i])
        mppstr = '%s_%s_mpp'%(places[i],mtypes[i])
        mrtstr = '%s_%s_mrt'%(places[i],mtypes[i])
        mrpstr = '%s_%s_mrp'%(places[i],mtypes[i])
        mtpstr = '%s_%s_mtp'%(places[i],mtypes[i])
        lonstr = '%s_%s_derived_longitude'%(places[i],mtypes[i])
        latstr = '%s_%s_derived_latitude'%(places[i],mtypes[i])
        depstr = '%s_%s_derived_depth'%(places[i],mtypes[i])
        if i == 0:
            thisone = dftot[np.isfinite(dftot[depstr].values)]
            others = dftot[np.isnan(dftot[depstr].values)]
        else:
            thisone = others[np.isfinite(others[depstr].values)]
            others = others[np.isnan(others[depstr].values)]

        thisone['mrr'] = thisone[mrrstr].values*1.0
        thisone['mtt'] = thisone[mttstr].values*1.0
        thisone['mpp'] = thisone[mppstr].values*1.0
        thisone['mrt'] = thisone[mrtstr].values*1.0
        thisone['mrp'] = thisone[mrpstr].values*1.0
        thisone['mtp'] = thisone[mtpstr].values*1.0
        thisone['moment_lon'] = thisone[lonstr].values*1.0
        thisone['moment_lat'] = thisone[latstr].values*1.0
        thisone['moment_depth'] = thisone[depstr].values*1.0
        thisone['type'] = '%s%s'%(places[i],mtypes[i])

        thisone = thisone[['id_no','time','lat','lon','depth','mag','mag_type','moment_lon','moment_lat','moment_depth','mrr','mtt','mpp','mrt','mrp','mtp','type']]
        
        newdf = pd.concat([newdf,thisone],sort=True)
    
    dftot = pd.concat([newdf,others],sort=True)
    dftot = dftot[['id_no','time','lat','lon','depth','mag','mag_type','moment_lon','moment_lat','moment_depth','mrr','mtt','mpp','mrt','mrp','mtp','type']]
    dftot.to_csv('%s.csv'%database,header=True,index=False,na_rep=np.nan)
    print ('All files were merged and written to %s.csv'%database)

elif manyortwo == '2':
    # for adding latest query to older file
    oldfile = input('enter name of existing catalog: ')
    newfile = input('enter name of catalog to add to existing: ')
    newcat = input('enter name of file to write combined catalogs to: ')

    new = pd.read_csv(newfile)
    old = pd.read_csv(oldfile)
    
    oldcolumns = old.columns

    newcatdf = pd.concat([old,new],sort=True)
    newcatdf = newcatdf.reset_index(drop=True)
    newcatdf = newcatdf[oldcolumns]
    
    newcatdf.to_csv(newcat,header=True,index=False,na_rep=np.nan)

else:
    print ('nothing happened, re-run and enter 1 or 2 as instructed above.')
