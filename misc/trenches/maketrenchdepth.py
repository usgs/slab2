import pandas as pd
import numpy as np
import os

trenches = pd.read_csv("trenches_usgs_2017.csv")
trenches.loc[trenches.lon < 0, "lon"] += 360
masterdir = "../../MasterDB/bathymetry"

BAdataset = pd.DataFrame()
slablist = []
for SGfile in os.listdir(masterdir):
    try:
        slab = SGfile[0:3]
        badata = pd.read_csv("%s/%s" % (masterdir, SGfile))
        badata["slab"] = slab
        BAdataset = pd.concat([BAdataset, badata])
        slablist.append(slab)
    except:
        print("couldnt read this file", SGfile)

newtrdata = pd.DataFrame()
for slab in slablist:

    BAdata = BAdataset[BAdataset.slab == slab]
    TRdata = trenches[trenches.slab == slab]

    depthlist = []
    for index, row in TRdata.iterrows():

        tlon, tlat = row["lon"], row["lat"]
        locba = BAdata[
            (BAdata.lon < tlon + 0.5)
            & (BAdata.lon > tlon - 0.5)
            & (BAdata.lat < tlat + 0.5)
            & (BAdata.lat > tlat - 0.5)
        ]
        if len(locba) > 0:
            locba["dist"] = sdist(tlat, tlon, locba["lat"], locba["lon"])
            mindist = locba["dist"].min()
            thisba = locba[locba.dist == mindist]
            if slab == "him":
                depthlist.append(thisba["elev"].values[0] / 1000)
            else:
                depthlist.append(thisba["depth"].values[0])
        else:
            depthlist.append(np.nan)

        print(slab, tlon, tlat)
    TRdata["depth"] = depthlist

    newtrdata = pd.concat([newtrdata, TRdata])

newtrdata.to_csv(
    "trenches_usgs_2017_depths.csv", header=True, index=False, na_rep=np.nan
)

# taken from https://github.com/usgs/neicmap (pip install no longer works)
def sdist(lat1, lon1, lat2, lon2):
    """
    Approximate great circle distance (meters) assuming spherical Earth (6367 km radius).
    @param lat1: Latitude(s) of first point(s).
    @param lon1: Longitude(s) of first point(s).
    @param lat2: Latitude(s) of second point(s).
    @param lon2: Longitude(s) of second point(s).
    @return: Vector of great circle distances, same length as longer of two input arrays of points.
    """
    R = 6367 * 1e3  # radius of the earth in meters, assuming spheroid
    dlon = lon1 - lon2
    t1 = pow((cosd(lat2) * sind(dlon)), 2)
    t2 = pow((cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(dlon)), 2)
    t3 = sind(lat1) * sind(lat2) + cosd(lat1) * cosd(lat2) * cosd(dlon)

    dsig = numpy.arctan2(numpy.sqrt(t1 + t2), t3)

    gcdist = R * dsig
    return gcdist
