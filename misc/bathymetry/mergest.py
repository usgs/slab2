import pandas as pd
import numpy as np
from scipy.interpolate import griddata

def findNotMatches2(dfo,dfm):
    ind = dfo.lonlat.isin(dfm.lonlat)
    dfo = dfo[~ind]
    return dfo

bathymetry = pd.read_csv('gebco5x5.csv')
sedthick1 = pd.read_csv('sedmap.csv')
sedthick2 = pd.read_csv('sedthick2-2.csv')

bathymetry.loc[bathymetry.lon<0,'lon'] += 360
sedthick1.loc[sedthick1.lon<0,'lon'] += 360
sedthick2.loc[sedthick2.lon<0,'lon'] += 360

sedthick2 = sedthick2[np.isfinite(sedthick2.thickness)]

print len(sedthick2)
print len(bathymetry)

sea_slab = pd.merge(sedthick2, bathymetry, how='inner', on=['lat','lon'])

print len(sea_slab)

sea_slab['lonlat'] = sea_slab.lon.astype(str).str.cat(sea_slab.lat.astype(str), sep='_')
bathymetry['lonlat'] = bathymetry.lon.astype(str).str.cat(bathymetry.lat.astype(str), sep='_')

missingmerge = findNotMatches2(bathymetry,sea_slab)
missingmerge = missingmerge.reset_index(drop=True)

print len(missingmerge)
print len(bathymetry)-len(missingmerge)

dataB = np.zeros((len(missingmerge),4))
dataS = np.zeros((len(sedthick1),3))

dataB[:,0] = missingmerge['lon'].values
dataB[:,1] = missingmerge['lat'].values
dataB[:,2] = missingmerge['elev'].values

dataS[:,0] = sedthick1['lon'].values
dataS[:,1] = sedthick1['lat'].values
dataS[:,2] = sedthick1['thickness'].values

dataB[:, 3] = griddata(dataS[:, 0:2], dataS[:, 2], dataB[:, 0:2], method='nearest')

missingmerge['thickness'] = dataB[:,3]*1000

allBA = pd.concat([sea_slab, missingmerge])
allBA['depth'] = (allBA['elev'].values - allBA['thickness'].values) * -0.001

allBA = allBA[['lon','lat','depth','elev','thickness']]
allBA.to_csv('totalBA.csv',header=True,index=False,na_rep=np.nan)


