import pandas as pd
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

nc_f = 'DNA13_percent.nc'  # Your filename
nc_fid = Dataset(nc_f, 'r')

lons = nc_fid.variables['longitude'][:]
lats = nc_fid.variables['latitude'][:]
deps = nc_fid.variables['depth'][:]
vps = nc_fid.variables['vp'][:]

y,z,x = np.meshgrid(lats,deps,lons)
print '2',x.shape, y.shape, z.shape, vps.shape

all = np.zeros((len(x.flatten()),4))

all[:,0] = x.flatten()
all[:,1] = y.flatten()
all[:,2] = z.flatten()
all[:,3] = vps.flatten()

all = all[(all[:,0] < 260)&(all[:,1]>30)]

df = pd.DataFrame({'lon':all[:,0]-360,'lat':all[:,1],'depth':all[:,2],'dvs':all[:,3]})
df = df[['lat','lon','depth','dvs']]
df.to_csv('ust_RWP14.tomo',header=False,index=False,sep=' ')

data = pd.read_csv('cas_TO_RWP14.csv')
data.loc[data.lon < 0, 'lon'] += 360

secondata = pd.read_csv('/Users/ginevramoore/Documents/Slab2/Ginevra2017/dataStructure/0904database/cas_TO_BSL14.csv')
secondata.loc[secondata.lon < 0, 'lon'] += 360
lats = np.linspace(np.min(all[:,1]),np.max(all[:,1]),50)

norm = MidpointNormalize(midpoint=0)

for i in range(len(lats)):
    all1 = all[(all[:,1] < lats[i]+1)&(all[:,1] > lats[i]-1)]
    data1 = data[(data.lat < lats[i]+1)&(data.lat > lats[i]-1)]
    data2 = secondata[(secondata.lat < lats[i]+1)&(secondata.lat > lats[i]-1)]

    fig = plt.figure(figsize=(15, 10))

    ax1 = fig.add_subplot(111)
    con = ax1.scatter(all1[:,0]*111.19,all1[:,2],c=all1[:,3],s=200,edgecolors='none',cmap='PRGn', norm=norm)
    ax1.scatter(data2.lon*111.19, data2.depth, c='k', s=15,edgecolors='none',label='BSL14')
    ax1.scatter(data1.lon*111.19, data1.depth, c='r', s=10,edgecolors='none',label='RWP14')
    ax1.set_ylabel('Depth')
    ax1.set_xlabel('Longitude * 111.19')
    ax1.axis('equal')
    ax1.grid()
    ax1.invert_yaxis()
    title = 'Latitude = %f'%lats[i]
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('vp')

    figtitle = 'tomotest/tiltdata_%f.png'%lats[i]
    fig.savefig(figtitle)
    plt.close()

'''
depths = np.linspace(np.min(all[:,2]),np.max(all[:,2]),50)
    
for i in range(len(depths)):
    all1 = all[(all[:,2] < depths[i]+10)&(all[:,2] > depths[i]-10)]

    fig = plt.figure(figsize=(10, 10))

    ax1 = fig.add_subplot(111)
    con = ax1.scatter(all1[:,0],all1[:,1],c=all1[:,3],s=100,edgecolors='none',cmap='plasma')
    ax1.set_ylabel('Latitude ')
    ax1.set_xlabel('Longitude')
    ax1.axis('equal')
    ax1.grid()
    title = 'depth = %f'%depths[i]
    ax1.set_title(title)
    ax1.legend(loc='best')
    cbar = fig.colorbar(con)
    cbar.set_label('vp')

    figtitle = 'tomotest/updata_%f.png'%depths[i]
    fig.savefig(figtitle)
    plt.close()
'''
