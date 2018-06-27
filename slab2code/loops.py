#!/usr/bin/env python

import numpy as np
import pandas as pd
import slab2functions as s2f
import psutil

def loop1(lons, lats, testarea, slab, depgrid, strgrid, dipgrid,
        slab1query, eventlist, seismo_thick,
        alen, blen, mdist, sdr, ddr, mindip, maxID2, AA_data,
        TR_data, maxdist, maxthickness, minstk,
        tomo_sets, meanBA, slab1guide, spacing, slab1data, dipthresh, datainfo, nodeinfo, j):
    
    ''' Arguments:  lons - list of longitudes (list of floats)[deg]
                    lats - list of latitudes (list of floats)[deg]
                    nodeIDs - list of nodeIDs (list of ints)
                    testarea - list defined as [lonmin,lonmax,latmin,latmax]
                            nodes within this box will print troubleshooting
                            statements and PDFs and data used to create them 
                            are saved (list of floats)[deg]
                    slab - three letter code to indicate slab of interest (str)
                    depgrid - clipped depth GMT grid file from Slab 1.0 (grd)
                    strgrid - grid of strike values calculated from 
                            Slab 1.0 (grd)
                    dipgrid - grid of dip values calculated from Slab 1.0 (grd)
                    slab1data - dataframe with columns [lon,lat,depth,strike,
                            dip] representing Slab1.0 (dataframe)
                    eventlist - dataframe with columns [lat,lon,depth,unc,etype,
                            ID,mag,S1,D1,R1,S2,D2,R2,src] representing input
                            dataset (dataframe)
                    seismo_thick - depth defining where earthquakes at and 
                            below are intra-slab as opposed to on the 
                            slab surface. (float)[km]
                    alen - long radius of ellipsoid (float)[km]
                    blen - short radius of search ellipsoid (float)[km]
                    mdist - alen*ecentricity of ellipse (float)[km]
                    sdr - shallow search range of ellipsoid (float)[km]
                    ddr - deep search range of ellipsoid (float)[km]
                    mindip - upper limit of dip to search vertically. Nodes with
                            dips greater than mindip have a perpendicular search
                            ellipsoid (float)[deg]
                    maxID2 - maxID of original input dataset (int)
                    AA_data - dataframe with columns [dist,depth].
                            if slab == 'alu' or slab == 'him', columns = 
                            [dist,depth,avlat,avlon] (dataframe)
                    TR_data - dataframe with columns [lon,lat,az,
                                bound,slab] (dataframe)
                    maxdist - useless variable, need to remove
                    maxthickness - maximum thickness of this slab as calculated
                            by a function relating thickness to age (float)[km]
                    minstk - maximum strike difference between node and farthest
                            data point constraining that node (float)[deg]
                    tomo_sets - list of tomography datasets for this slab.
                            (list of strings)
                    OGmulti - dataframe with columns [lon,lat,depth,nID] used to
                            record which nodes have a multi-modal distribution
                            and the depths to each peak (dataframe)
                    i - index representing this node. (int)
        
        Returns:    lon - lons[i], search longitude of this node (float)[deg]
                    lat - lats[i], search latitude of this node (float)[deg]
                    locdep - search depth of this node (float)[km]
                    locstr - search strike of this node (float)[deg]
                    locdip - search dip of this node (float)[deg]
                    nID - integer value representing this node (int)
                    farpoint - length of long axis of ellipse after 
                            refiltering based on strike difference from 
                            node center (float)[km]
                    used_TO - list of tomography data points used to constrain
                            nodes on this slab. 2D numpy array with columns
                            [mean(|tomo_depths-EQdepths|),index representing
                            TO source] (array of floats)
                    used_tmp - list of data points used to constrain the depth
                            at this node. 2D numpy array with columns
                            [nodeID,dataID] (array of ints)
                    trimmedAA - dataframe with columns [lat,lon,depth,unc,etype,
                            ID,mag,S1,D1,R1,S2,D2,R2,src] representing average
                            active source and receiver functions to 
                            add to input dataset (dataframe) '''
    
    # Define lon,lat coordinate, integer ID of this node, and ID for added data
    lat = lats[j]
    lon = lons[j]
    if lat>0:
        nID = int('%i%i'%(lon*10,lat*10))
    else:
        nID = int('%i0%i'%(lon*10,lat*-10))
    maxID = maxID2+nID-1
    OGmulti = []
    
    if slab == 'sol' and lon > 149 and lon < 155:
        alen = 20.0
    if slab == 'puy' and lat < -48:
        alen = 100
        
    if slab == 'ryu' and lon > 138:
        alen = 20.0
    
    if slab == 'alu' and lat > 63:
        alen = 50.0
    
    if slab == 'man' and lat > 17.0:
        alen = 50.0
    
    if slab == 'izu' and lat < 27 and lat > 23:
        eventlist = eventlist[(eventlist.lon<lon+4.5)&(eventlist.lat<lat+4.5)&(eventlist.lon>lon-4.5)&(eventlist.lat>lat-4.5)]
    else:
        eventlist = eventlist[(eventlist.lon<lon+2.5)&(eventlist.lat<lat+2.5)&(eventlist.lon>lon-2.5)&(eventlist.lat>lat-2.5)]
    if len(eventlist)<2:
        #print ('no data here',lon,lat,nID)
        return lon, lat, np.nan, np.nan, np.nan, nID, np.nan, np.nan, np.nan, [], [], [], [], False, 0, 0

    # If in troubleshooting range, print ID and set testprint = True
    if lon>testarea[0] and lon<testarea[1]:
        if lat>testarea[2] and lat<testarea[3]:
            testprint = True
            print('___________________ %i ______________________'% nID)
        else:
            testprint = False
    else:
        testprint = False

    if testprint:
        f = open(nodeinfo, 'a')
        f.write('-%i- pre-anyfilters, len(eventlist) %i\n'%(nID, len(eventlist)))
        f.close()
    # Get depth, strike, and dip of Slab1.0 or calculate from trend at edge
    slab1, strtmp, diptmp, inside, extended, out, extlon, extlat = s2f.getslab12(slab1guide,slab1query,lon,lat,spacing,depgrid,strgrid,dipgrid,testprint,TR_data,meanBA,slab)

    if slab1 > 650 and slab == 'sam' and lat<-30:
        sdr = 60
    if testprint:
        f = open(nodeinfo, 'a')
        f.write('-%i- slab1, strtmp, diptmp, inside, extended, out %0.2f,%0.2f,%0.2f,%s,%s,%s \n'%(nID, slab1, strtmp, diptmp, inside, extended, out))
        f.close()

    if slab == 'izu' and lat < 27 and lat > 23:
        dipthresh = 45

    if diptmp > dipthresh:
        newlats = [lat,lat,lat+spacing/2.0,lat-spacing/2.0,lat+spacing/2.0,lat-spacing/2.0,lat+spacing/2.0,lat-spacing/2.0]
        newlons = [lon-spacing/2.0,lon+spacing/2.0,lon,lon,lon+spacing/2.0,lon-spacing/2.0,lon-spacing/2.0,lon+spacing/2.0]
        
        if slab == 'ker':
            newlats1 = [lat,lat,lat+spacing/4.0,lat-spacing/4.0,lat+spacing/4.0,lat-spacing/4.0,lat+spacing/4.0,lat-spacing/4.0]
            newlons1 = [lon-spacing/4.0,lon+spacing/4.0,lon,lon,lon+spacing/4.0,lon-spacing/4.0,lon-spacing/4.0,lon+spacing/4.0]

            newlats2 = [lat,lat,lat+3*spacing/4.0,lat-3*spacing/4.0,lat+3*spacing/4.0,lat-3*spacing/4.0,lat+3*spacing/4.0,lat-3*spacing/4.0]
            newlons2 = [lon-3*spacing/4.0,lon+3*spacing/4.0,lon,lon,lon+3*spacing/4.0,lon-3*spacing/4.0,lon-3*spacing/4.0,lon+3*spacing/4.0]
        
            newlats.extend(newlats1)
            newlats.extend(newlats2)
            newlons.extend(newlons1)
            newlons.extend(newlons2)
            
        newnodes = np.zeros((len(newlons),2))
        newnodes[:,0] = np.round(newlons,1)
        newnodes[:,1] = np.round(newlats,1)
        #diptmp = 90
        if slab == 'izu' and lat>23 and lat<27 and slab1 > 400:
            clen = 400
        elif slab != 'sol' and slab != 'hal' and slab != 'phi':
            clen = 200
        elif slab == 'sol' and lon>148:
            clen = 300
        else:
            clen = alen
    else:
        eventlist = eventlist[(eventlist.lon<lon+1.5)&(eventlist.lat<lat+1.5)&(eventlist.lon>lon-1.5)&(eventlist.lat>lat-1.5)]
        if testprint:
            f = open(nodeinfo, 'a')
            f.write('-%i- after reducing OGsearch range (dip<dipthresh), len(eventlist) %i\n'%(nID, len(eventlist)))
            f.close()
        clen = alen
        if len(eventlist)<2 and len(eventlist[(eventlist.etype == 'AS')|(eventlist.etype == 'RF')|(eventlist.etype == 'CP')])<1 and slab != 'ryu':
            #print ('no data here',lon,lat,nID)
            return lon, lat, slab1, strtmp, diptmp, nID, alen, blen, clen, [], [], [], [], False, 0, 0
        elif slab == 'ryu' and len(eventlist)<2 and len(eventlist[(eventlist.etype == 'AS')|(eventlist.etype == 'RF')|(eventlist.etype == 'TO')|(eventlist.etype == 'CP')])<1:
            return lon, lat, slab1, strtmp, diptmp, nID, alen, blen, clen, [], [], [], [], False, 0, 0
        else:
            newnodes = []

    if slab == 'phi' and lat > 11:
        minstk /= 2

    if extended:
        alen = s2f.refilter4(strtmp,minstk,extlon,extlat,eventlist,alen,blen,slab1query,slab,testprint)
    else:
        alen = s2f.refilter4(strtmp,minstk,lon,lat,eventlist,alen,blen,slab1query,slab,testprint)

    if lon > 170.5 and lon < 172 and lat > -22 and lat < -21:
        alen = 20
    if testprint:
        f = open(nodeinfo, 'a')
        f.write('-%i- alen, blen, clen, minstk, seismo_thick %0.2f,%0.2f,%0.2f,%0.2f,%0.2f \n'%(nID, alen, blen, clen, minstk, seismo_thick))
        f.close()
    
    if alen == 1:
    	return lon, lat,slab1, strtmp, diptmp, nID, alen, blen, clen, [], [], [], [], False, 0, 0
    else:
    	farpoint = alen
    
    if slab == 'cam' and lon < 260 and lat > 20.5 and lon > 257.5:
        alen = 20
        clen = 20
    
    if slab == 'man' and lat < 14:
        alen = 20
        if slab1 < 200:
            clen = 20
    
    if slab == 'mue':
        alen = 150
        clen = 150
        strtmp = 270
        
    # Filter data for present node
    trimmed, test, sdepth, ddepth, locstr, locdip, maxID, locdep = s2f.allFilters(
                                eventlist, lat, lon, inside, slab1, strtmp, diptmp, seismo_thick, alen, blen,
                                clen, mdist, sdr, ddr, mindip, maxID, out,
                                AA_data, TR_data, slab, maxdist, testprint, extended, datainfo, nodeinfo,nID)

    trimmed, test = s2f.doublecheckEREQ(trimmed, lon, lat)

    if not test:
        if testprint:
            print ('not test',lon,lat,nID)
        return lon, lat, locdep, locstr, locdip, nID, alen, blen, clen, [], [], [], newnodes, False, sdepth, ddepth

    trimmedAA = trimmed[trimmed.etype == 'AA']
    trimmedBA = trimmed[trimmed.etype == 'BA']
    trimmedAS = trimmed[trimmed.etype == 'AS']
    trimmedRF = trimmed[trimmed.etype == 'RF']
    trimmedTO = trimmed[trimmed.etype == 'TO']
    trimmedCP = trimmed[trimmed.etype == 'CP']
    # Save data used
    # Add this trimmed list of events to the list of used input data
    if len(trimmed) > 1 or (len(trimmed) < 2 and (len(trimmedAA)>0 or len(trimmedBA)>0 or len(trimmedAS)>0 or len(trimmedRF)>0 or len(trimmedCP)>0)):
        trimmed_index = np.array(trimmed.ID)
        node_number = np.ones(len(trimmed_index)) * nID
        used_tmp = np.array(list(zip(node_number, trimmed_index)))
    elif len(trimmed) > 1 or (len(trimmed) < 2 and len(trimmedTO)>0 and slab == 'ryu'):
        trimmed_index = np.array(trimmed.ID)
        node_number = np.ones(len(trimmed_index)) * nID
        used_tmp = np.array(list(zip(node_number, trimmed_index)))
    else:
        if testprint:
            print ('trimmed too small',lon,lat,nID)
        return lon, lat, locdep, locstr, locdip, nID, alen, blen, clen, [], [], [], newnodes, False, sdepth, ddepth

    # We only need to do the rest of this loop if there is tomography data
    if tomo_sets == 0 or slab != 'samzzz' or slab =='ryu':
        return lon, lat, locdep, locstr, locdip, nID, alen, blen, clen, [], used_tmp, trimmedAA, newnodes, True, sdepth, ddepth

    '''nothing beyond this point in this loop gets used '''
    
    # Trim dataset for PDF calculation and separate tomography data
    # We want to run PDF on non-TO, non-RF data
    TOdata = trimmed[trimmed.etype == 'TO']
    trimmed = trimmed[trimmed.etype != 'RF']
    trimmed = trimmed[trimmed.etype != 'TO']
    trimmed = trimmed[trimmed.etype != 'CP']

    # Calculate PDF
    if locdip <= mindip: #or len(trimmedAS) > 0 or len(trimmedAA) >0 or len(trimmedRF) >0:
        #print 'full',lon,lat
        peak_depth, stdv, test2, npeaks, OGmulti, centsurf = s2f.fullPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, locdep, 'pre', slab, locstr, locdip)
        peak_lon = lon
        peak_lat = lat
    else:
        peak_lon, peak_lat, peak_depth, stdv, test2, npeaks, OGmulti, centsurf, rlist = s2f.perpPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, locdep, 'pre', slab, locstr, locdip, maxthickness)
        if len(rlist) > 0:
            removeIDs = np.array(rlist.ID)
            used_tmp = used_tmp[~np.in1d(used_tmp[:, 1], removeIDs)]

    # For each node with a peak depth and tomography data, calculate differences
    # Do separately for each tomography dataset
    if test2:
        if len(TOdata) > 0:
            for idx, src in enumerate(tomo_sets):
                data = TOdata[TOdata.src == src]
                if len(data) > 0:
                    TOmean = np.mean(data.depth)
                    diff = peak_depth - TOmean
                    temp = [diff, idx]
                    temp2 = np.array(temp)
                    if used_TO is not None:
                        used_TO = np.vstack((used_TO, temp2))
                    else:
                        used_TO = temp2.copy()
                else:
                    return lon, lat, locdep, locstr, locdip, nID, farpoint, used_TO, used_tmp, trimmedAA, newnodes

    return lon, lat, locdep, locstr, locdip, nID, alen, blen, clen, used_TO, used_tmp, trimmedAA, newnodes

def loop2(testarea, lons, lats, nIDs1, locdep, locstr, locdip, used_all, eventlist, sdr, ddr, seismo_thick, slab, maxthickness, rlist, mindip, aalen, ablen, aclen, j):
    these_parameters = []

    premulti = []
    bilats, bilons, binods, bistds = np.nan, np.nan, np.nan, np.nan
    biindx, bistrs, bidips, bideps = np.nan, np.nan, np.nan, np.nan
    
    lat = lats[j]
    lon = lons[j]
    nID = nIDs1[j]
    cstr = locstr[j]
    cdip = locdip[j]
    alen = aalen[j]
    blen = ablen[j]
    clen = aclen[j]
    
    if lon>testarea[0] and lon<testarea[1] and lat>testarea[2] and lat<testarea[3]:
        testprint = True
        print('____________________________________ %i ___________________________________________'% nID)
    else:
        testprint = False

    # Collect used data associated with this node
    loc_used_data = used_all[:, 1][used_all[:, 0] == nID]

    if len(loc_used_data) > 1:
        trimmed = eventlist[eventlist['ID'].isin(loc_used_data)]
    elif len(loc_used_data) > 0:
        trimmed = eventlist[eventlist['ID'].isin(loc_used_data)]
        trimmedAA = trimmed[trimmed.etype == 'AA']
        trimmedBA = trimmed[trimmed.etype == 'BA']
        trimmedAS = trimmed[trimmed.etype == 'AS']
        trimmedRF = trimmed[trimmed.etype == 'RF']
        trimmedCP = trimmed[trimmed.etype == 'CP']
        if len(trimmedAA)<1 and len(trimmedBA)<1 and len(trimmedAS)<1 and len(trimmedRF)<1 and len(trimmedCP)<1:
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti, alen, blen, clen, 0
    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti, alen, blen, clen, 0

    # Remove RF data for PDF analysis
    if locdep[j] > seismo_thick:
        rft = trimmed[(trimmed.etype == 'RF') | (trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'CP')]
        trimmed = trimmed[(trimmed.etype != 'RF') & (trimmed.etype != 'AA') & (trimmed.etype != 'AS') & (trimmed.etype != 'CP')]
    else:
        rft = []

    if len(trimmed)<2 and len(rft)>0:
        trimmed = pd.concat([rft, trimmed])
        if abs(rft['depth'].mean()-trimmed['depth'].mean())>maxthickness:
            trimmed = trimmed[(trimmed.etype == 'RF') | (trimmed.etype == 'AA') | (trimmed.etype == 'AS') | (trimmed.etype == 'CP')]
    if len(rft)>0:
        loc_depth = np.mean(trimmed['depth'].values)
        if np.isnan(loc_depth):
            if testprint:
                print ('1 ** nan depth found for this node, returning nan output in this directory',nID)
                trimmed.to_csv('%i_faileddata.csv'%nID,header=True,index=False,na_rep=np.nan)
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti, alen, blen, clen, 0
    else:
        loc_depth = locdep[j]
        if np.isnan(loc_depth):
            if testprint:
                print ('2 ** nan depth found for this node, returning nan output in this directory',nID)
                trimmed.to_csv('%i_faileddata.csv'%nID,header=True,index=False,na_rep=np.nan)
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti, alen, blen, clen, 0

    place_holder, sdepth, ddepth, depthwritten = s2f.depthRange(loc_depth, sdr, ddr, seismo_thick, trimmed, slab, these_parameters, True)
    

    if cdip <= mindip:# or len(trimmed[trimmed.etype == 'AS' ]) > 0 or len(trimmed[trimmed.etype == 'AA' ]) > 0 or len(rft)>0:
        if cdip > mindip:
            cdip = mindip
        #try:
        peak_depth, stdv, test2, npeaks, premulti, centsurf = s2f.fullPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, 'pre2', slab, cstr, cdip)
        #except:
        #    print('full PDF didnt work! lon,lat,nID,cdip,cstr,loc_depth,trimmed', lon, lat, nID, cdip, cstr, loc_depth, trimmed)
        #    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti
        peak_lon = lon
        peak_lat = lat
    else:
        #try:
        peak_lon, peak_lat, peak_depth, stdv, test2, npeaks, premulti, centsurf, rlist = s2f.perpPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, 'pre2', slab, cstr, cdip, maxthickness)
        #except:
        #    print('perp PDF didnt work!lon,lat,nID,cdip,cstr,loc_depth,trimmed', lon, lat, nID, cdip, cstr, loc_depth, trimmed)
        #    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, rlist, premulti

    if test2:
        #tmp_res[j,0],tmp_res[j,1],tmp_res[j,2],tmp_res[j,3],tmp_res[j,4],tmp_res[j,5],tmp_res[j,6],tmp_res[j,7],tmp_res[j,8],tmp_res[j,9] = peak_lat,peak_lon,peak_depth,std,nID,lat,lon,cstr,cdip,centsurf
        jpeak_lat, jpeak_lon, jpeak_depth, jstd, jnID = peak_lat, peak_lon, peak_depth, stdv, nID
        jlat, jlon, jcstr, jcdip, jcentsurf = lat, lon, cstr, cdip, centsurf
        # If there is more than one peak
        if npeaks > 1:
            bilats, bilons, binods, bistds= lat, lon, nID, stdv
            biindx, bistrs, bidips, bideps = j, cstr, cdip, loc_depth
    else:
        jpeak_lon, jpeak_lat, jpeak_depth, jstd, jnID = np.nan, np.nan, np.nan, np.nan, np.nan
        jlat, jlon, jcstr, jcdip, jcentsurf = np.nan, np.nan, np.nan, np.nan, np.nan

    eventtypes = list(set(trimmed.etype))
    if len(eventtypes) == 1 and eventtypes[0] == 'TO':
        onlyTO = 1
    else:
        onlyTO = 0

    return jpeak_lon, jpeak_lat, jpeak_depth, jstd, jnID, jlat, jlon, jcstr, jcdip, jcentsurf, bilats, bilons, binods, bistds, biindx, bistrs, bidips, bideps, rlist, premulti, alen, blen, clen, onlyTO

def loop3(shift_out, testarea, used_all, eventlist, sdr, ddr, seismo_thick, these_parameters, slab, maxthickness, mindip, taper, nID):

    postmulti = []
    bilats, bilons, binods, bistds = np.nan, np.nan, np.nan, np.nan
    biindx, bistrs, bidips, bideps = np.nan, np.nan, np.nan, np.nan
    
    if len(shift_out)>0: # used to be if slab == 'sam' bc of variable taper
        taper = 0.0
        
    loc_shift_out = shift_out[shift_out.nID == nID]
    lat = loc_shift_out['lat'].values[0]
    lon = loc_shift_out['lon'].values[0]
    cstr = loc_shift_out['ogstr'].values[0]
    cdip = loc_shift_out['ogdip'].values[0]
    bflat = loc_shift_out['bzlat'].values[0]
    bflon = loc_shift_out['bzlon'].values[0]
    dep = loc_shift_out['psdepth'].values[0]
    
    if len(loc_shift_out) == 0:
        #print 'misnamed ID',nID
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, nID, postmulti, np.nan, np.nan
    else:
        testprint = False
    if lon>testarea[0] and lon<testarea[1] and lat>testarea[2] and lat<testarea[3]:
        testprint = True
    else:
        testprint = False
    if testprint:
        print('____________________________________ %i ___________________________________________'% nID)

    # Collect used data associated with this node
    loc_used_data = used_all[:, 1][used_all[:, 0] == nID]
    if len(loc_used_data) > 1:
        trimmed = eventlist[eventlist['ID'].isin(loc_used_data)]
    elif len(loc_used_data) > 0:
        trimmed = eventlist[eventlist['ID'].isin(loc_used_data)]
        trimmedAA = trimmed[trimmed.etype == 'AA']
        trimmedBA = trimmed[trimmed.etype == 'BA']
        trimmedAS = trimmed[trimmed.etype == 'AS']
        trimmedRF = trimmed[trimmed.etype == 'RF']
        trimmedCP = trimmed[trimmed.etype == 'CP']
        if len(trimmedAA)<1 and len(trimmedBA)<1 and len(trimmedAS)<1 and len(trimmedRF)<1 and len(trimmedCP)<1:
            #if testprint:
            print('skipping but returning something wonky', lon, lat, dep)
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, nID, postmulti, np.nan, np.nan
    else:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, nID, postmulti, np.nan, np.nan
    # Calculate average strike, rake, dip of CMTs for final output
    avg_strike, avg_dip, avg_rake = s2f.avStrDipRak(trimmed)

    if dep < seismo_thick or (len(trimmed[trimmed.etype == 'RF']) < 1 and \
                    len(trimmed[trimmed.etype == 'AA']) < 1 and \
                    len(trimmed[trimmed.etype == 'CP']) < 1 and \
                    len(trimmed[trimmed.etype == 'AS']) < 1):
        #print 'shallow',lon,lat,nID
        rstd = loc_shift_out['stdv'].values[0]
        depth = loc_shift_out['depth'].values[0]
        if dep<seismo_thick:
            return depth, rstd, avg_strike, avg_dip, avg_rake, bilats, bilons, binods, bistds, biindx, bistrs, bidips, bideps, nID, postmulti, bflon, bflat
        else:
            return depth, rstd, avg_strike, avg_dip, avg_rake, bilats, bilons, binods, bistds, biindx, bistrs, bidips, bideps, nID, postmulti, lon, lat

    if len(trimmed[(trimmed.etype != 'RF') & (trimmed.etype != 'CP')])<1 and len(trimmed[(trimmed.etype == 'RF') | (trimmed.etype == 'CP')])>0:
        rstd = loc_shift_out['stdv'].values[0]
        depth = loc_shift_out['psdepth'].values[0]
        return depth, rstd, avg_strike, avg_dip, avg_rake, bilats, bilons, binods, bistds, biindx, bistrs, bidips, bideps, nID, postmulti, bflon, bflat

    if len(trimmed[(trimmed.etype != 'RF') & (trimmed.etype != 'AS') & (trimmed.etype != 'AA') & (trimmed.etype != 'CP')]) > 0:
        # Calculate local depth difference between shifted and non-shifted results
        dep_dif = dep - loc_shift_out['depth'].values[0]
        lon_dif = bflon - lon
        lat_dif = bflat - lat
        
        # Apply shift to non-RF data
        trimmed['depth'][(trimmed.etype != 'RF') & (trimmed.etype != 'AS') & (trimmed.etype != 'AA') & (trimmed.etype != 'CP')] -= dep_dif
        trimmed['lon'][(trimmed.etype != 'RF') & (trimmed.etype != 'AS') & (trimmed.etype != 'AA') & (trimmed.etype != 'CP')] -= lon_dif
        trimmed['lat'][(trimmed.etype != 'RF') & (trimmed.etype != 'AS') & (trimmed.etype != 'AA') & (trimmed.etype != 'CP')] -= lat_dif

    # Calculate new PDF
    loc_depth = trimmed['depth'].mean()
    place_holder, sdepth, ddepth, depthwritten = s2f.depthRange(loc_depth, sdr, ddr, seismo_thick, trimmed, slab, these_parameters, True)

    if cdip <= mindip:
        if cdip > mindip:
            cdip = mindip
        try:
            new_peak_depth, new_std, test3, npeaks, postmulti, centsurf = s2f.fullPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, 'post', slab, cstr, cdip)
            peak_lon = lon
            peak_lat = lat
        except:
            print('post full pdf didnt work!!', lon, lat, nID)
            new_peak_depth, new_std, test3, npeaks, postmulti, centsurf = loc_shift_out['depth'].values[0], loc_shift_out['stdv'].values[0], False, 1, pd.DataFrame(), 2000
            peak_lon = lon
            peak_lat = lat
    else:
        try:
            peak_lon, peak_lat, new_peak_depth, new_std, test3, npeaks, postmulti, centsurf, rlist = s2f.perpPDFcalc(trimmed, sdepth, ddepth, testprint, nID, lat, lon, loc_depth, 'post', slab, cstr, cdip, maxthickness)
        except:
            print('post perp pdf didnt work!!', lon, lat, nID)
            peak_lon, peak_lat, new_peak_depth, new_std, test3, npeaks, postmulti, centsurf, rlist = lon, lat, loc_shift_out['depth'].values[0], loc_shift_out['stdv'].values[0], False, 1, pd.DataFrame(), 2000, pd.DataFrame()
    # Populate results array
    if test3:
        rdepth, rstd, rstrike, rdip, rrake = new_peak_depth, new_std, avg_strike, avg_dip, avg_rake
        if npeaks > 1:
            bilats, bilons, binods, bistds= lat, lon, nID, new_std
            biindx, bistrs, bidips, bideps = 1, cstr, cdip, loc_depth
    else:
        rdepth, rstd, rstrike, rdip, rrake = np.nan, np.nan, np.nan, np.nan, np.nan

    return rdepth, rstd, rstrike, rdip, rrake, bilats, bilons, binods, bistds, biindx, bistrs, bidips, bideps, nID, postmulti, peak_lon, peak_lat


