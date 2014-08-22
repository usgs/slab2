#!/usr/bin/env python

#stdlib imports

#third party imports
import numpy as np
import math

#usgs imports
from neicmap.distance import sdist,getAzimuth

def getEventsInCircle(lat,lon,radius,eventlist):
    """
    Select earthquakes from a list of events using a lat/lon center and a radius in km.
    @param lat: Latitude at center of circle.
    @param lon: Longitude at center of circle.
    @param radius: Radius of search in km.
    @param eventlist: List of event dictionaries, each containing (at least) keys lat,lon,depth,mag.
    @return: Filtered list of earthquake dictionaries.
    """
    elat = np.array([e['lat'] for e in eventlist]) #extract the latitudes from each dictionary in the list
    elon = np.array([e['lon'] for e in eventlist]) #extract the longitudes from each dictionary in the list
    d = sdist(lat,lon,elat,elon)/1000.0 #do vectorized distance calculation from center point to all lat/lon pairs
    az = []
    for i in range(0,len(elat)):
        az.append(getAzimuth(lat,lon,elat[i],elon[i]))  #vectorize the azimuth function someday

    for i in range(0,len(eventlist)):
        eventlist[i]['distance'] = d[i] #add distance field to each dictionary in the list
        eventlist[i]['azimuth'] = az[i] #add azimuth field to each dictionary in the list

    idx = (d <= radius).nonzero()[0] #array of indices where distance is less than input threshold
    newlist = (np.array(eventlist)[idx]).tolist() #return a shortened list
    
    return newlist

def getEventsInEllipse(lat,lon,str,aval,bval,eventlist):
    """
    Select earthquakes from a list of events using dist/az from a predefined lat/lon, and axes of an ellipse
    getEventsInCircle must already have been used
    @param str: Strike of slab at center of circle.
    @param aval: Long axis of ellipse in km.
    @param bval: Short axis of ellipse in km.
    @param eventlist: List of event dictionaries, each containing (at least) keys lat,lon,depth,mag.
    @return: Filtered list of earthquake dictionaries.
    """
    elat = np.array([e['lat'] for e in eventlist]) #extract the latitudes from each dictionary in the list
    elon = np.array([e['lon'] for e in eventlist]) #extract the longitudes from each dictionary in the list
    d = sdist(lat,lon,elat,elon)/1000.0 #do vectorized distance calculation from center point to all lat/lon pairs
    az = []
    for i in range(0,len(elat)):
        az.append(getAzimuth(lat,lon,elat[i],elon[i]))  #vectorize the azimuth function someday

    for i in range(0,len(eventlist)):
        eventlist[i]['distance'] = d[i] #add distance field to each dictionary in the list
        eventlist[i]['azimuth'] = az[i] #add azimuth field to each dictionary in the list

#   print 'Eventlist[1] = %.4f, %.4f, %.4f, %.4f' % (elat[1],elon[1],d[1],az[1]) 

    mdist = []
    erta = math.sqrt(1-((math.pow(bval,2))/(math.pow(aval,2))))

#    print 'ERTA = %.4f' % (erta)

    for i in range(0,len(elat)):
        mdist.append(getEllipseRad(aval,erta,az[i],str))

#    print 'Lon:', (elon[1:-1])
#    print 'Lat:', (elat[1:-1])
#    print 'Dist:', (d[1:-1])
#    print 'Azi:', (az[1:-1])
#    print 'Edist:', (mdist[1:-1])

    idx = (d <= mdist).nonzero()[0]  #array of indices where distance is less than input threshold

#    print 'Indexes:'
#    print (idx)[1:-1]

    newerlist = (np.array(eventlist)[idx]).tolist() #return a shortened list

    return newerlist

def getEllipseRad(a,e,azi,ang):
    d2r = math.pi/180
    d2 = (a*(1-e**2))/(1+(e*(math.cos((azi-ang)*d2r))))

    return d2

def heading(lat,lon,dist,az):
    """
    Project a location along a great circle by a distance dist, towards azimuth az
    """

    d2r = math.pi/180
    r2d = 180/math.pi
    if(dist<0):
        dist = dist*-1
        az = az-180

    if(az<0):
        az = az+360

    if(az>360):
        az = az-360

    b = (90-lat)*d2r
    a = (dist/111.19)*d2r
    angC = az*d2r
    c = math.cos(a)*math.cos(b)+math.sin(a)*math.sin(b)*math.cos(angC)
    c = math.acos(c)
    cdeg = c*r2d
    lat1 = 90-cdeg

    angA = (math.cos(a)-(math.cos(b)*math.cos(c)))/(math.sin(b)*math.sin(c))
    angA = math.acos(angA)
    adeg = angA*r2d

    if(az>0 and az<=180):
        lon1 = lon+adeg
    else:
        lon1 = lon-adeg

    return lat1,lon1

def main():
    eventlist = [{'lat':33.1,'lon':-118.1,'depth':15.1,'mag':5.4},
                 {'lat':33.2,'lon':-118.2,'depth':15.2,'mag':5.2},
                 {'lat':33.3,'lon':-118.3,'depth':15.3,'mag':5.1},
                 ]
    lat = 31.5
    lon = -119.5
    radius = 225.0
    shortlist = getEventsInCircle(lat,lon,radius,eventlist)
    print 'All events within radius of %.1f km from %.4f,%.4f' % (radius,lat,lon)
    for event in shortlist:
        print event

if __name__ == '__main__':
    main() #do all the work in a main function so sub functions aren't polluted by variables in global namespace
    
                  
