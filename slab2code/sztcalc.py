#!/usr/bin/env python

# import libraries
from datetime import datetime
import os.path
import argparse
import numpy as np
import pandas as pd
import warnings
import math
import mapio.gmt as gmt
from scipy import ndimage
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from copy import deepcopy
from pylab import arccos,argsort,cross,dot,double,eigh,pi,trace,zeros
from sklearn import mixture
from sklearn.metrics import mean_squared_error
import slab2functions as s2f

def main(args):

    slabsbyfile = args.slabsbyfile
    slabsbyslab = args.slabsbyslab
    inputfolder = args.inputfolder
    inputdate = args.inputdate
    origorcentl = args.origorcentl
    origorcentd = args.origorcentd
    if args.minlength is None:
        minlength = 999999999
    else:
        minlength = args.minlength
    maxdep = args.maxdep
    maxdepdiff = args.maxdepdiff
    slaborev = args.slaborev

    pd.options.mode.chained_assignment = None
    warnings.filterwarnings("ignore", message="invalid value encountered in less")
    warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
    warnings.filterwarnings("ignore", message="invalid value encountered in greater")
    warnings.filterwarnings("ignore", message="invalid value encountered in double_scalars")

    # list folder names of all slab models to calculate szt from, and move to new folder structure
    filelist = ['alu_slab2_09.20.19','cal_slab2_09.19.19','cam_slab2_09.19.19','car_slab2_09.19.19','cas_slab2_09.19.19','cot_slab2_09.20.19','hal_slab2_09.19.19','hel_slab2_09.19.19','him_slab2_09.19.19','hin_slab2_09.19.19','izu_slab2_09.20.19','ker_slab2_09.19.19','kur_slab2_09.19.19','mak_slab2_09.19.19','man_slab2_09.19.19','mue_slab2_09.20.19','pam_slab2_09.19.19','phi_slab2_09.20.19','png_slab2_09.19.19','puy_slab2_09.19.19','ryu_slab2_09.20.19','sam_slab2_09.19.19','sco_slab2_09.19.19','sol_slab2_09.19.19','sul_slab2_09.20.19','sum_slab2_09.19.19','van_slab2_09.20.19']

    # create new directory system for slab output
    os.system('rm -r %s'%slabsbyfile)
    os.system('mkdir %s'%slabsbyfile)
    os.system('mkdir %s/grids'%slabsbyfile)
    os.system('mkdir %s/supplement'%slabsbyfile)
    os.system('mkdir %s/nodes'%slabsbyfile)
    os.system('mkdir %s/surfacetext'%slabsbyfile)
    os.system('mkdir %s/clippingmasks'%slabsbyfile)
    os.system('mkdir %s/filtereddata'%slabsbyfile)
    os.system('mkdir %s/inputdata'%slabsbyfile)
    os.system('mkdir %s/parameters'%slabsbyfile)
    os.system('mkdir %s/szt'%slabsbyfile)
    os.system('mkdir %s/maps'%slabsbyfile)
    os.system('mkdir %s/crossections'%slabsbyfile)

    printtest = False

    figz = plt.figure(figsize=(15, 10))
    ax1z = figz.add_subplot(111)
    n = 0
    clist = ['yellowgreen','yellow','wheat','violet','turquoise','teal','silver','sienna','salmon','red','plum','pink','orchid','orange','olive','navy','magenta','lightgreen','lightblue','green','goldenrod','cyan','blue','darkgreen','yellowgreen','yellow','wheat','violet','turquoise']

    slabdf = []
    shaldf = []
    deepdf = []
    peakdf = []
    nevsdf = []
    savedir = '%s/szt'%slabsbyfile

    for folder in filelist:

        (slab,s2k,date) = folder.split('_')

        print ('                ---   %s   ---'%folder)
        inFile = '%s/%s_%s_input.csv'%(inputfolder,slab,inputdate)
        fname = '%s/%s/%s_slab2_dep_%s.grd'%(slabsbyslab,folder,slab,date)

        thisfolder = '%s/%s/%s_slab2'%(slabsbyslab,folder,slab)

        eventlistALL = pd.read_table('%s' % inFile, sep=',', dtype={
                'lon': np.float64, 'lat': np.float64,'depth': np.float64,
                'unc': np.float64, 'etype': str, 'ID': str, 'mag': np.float64,
                'S1': np.float64, 'D1': np.float64, 'R1': np.float64,
                'S2': np.float64, 'D2': np.float64, 'R2': np.float64,
                'src': str, 'time': str, 'mlon': np.float64, 'mlat': np.float64,
                'mdep': np.float64, 'id_no':str})

        eventlistALL['ID'] = eventlistALL['id_no'].values
        ogcolumns = ['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', \
                    'S1', 'D1', 'R1','S2', 'D2', 'R2', 'src']
        kagancols = ['lat', 'lon', 'depth', 'unc', 'etype', 'ID', 'mag', 'time', \
                    'S1', 'D1', 'R1','S2', 'D2', 'R2', 'src', 'mlon', 'mlat', 'mdep']

        eventlist = eventlistALL[kagancols]

        depgrid = gmt.GMTGrid.load(fname)
        strgrid, dipgrid = s2f.mkSDgrd(depgrid)
        slab1data = s2f.mkSlabData(depgrid, strgrid, dipgrid, printtest)

        slab1data.loc[slab1data.lon < 0, 'lon'] += 360
        eventlist.loc[eventlist.lon < 0, 'lon'] += 360
        eventlist.loc[eventlist.mlon < 0, 'mlon'] += 360

        try:
            maskdf = pd.read_csv('%s_clp_%s.csv'%(thisfolder,date), delim_whitespace=True, names=['lon','lat'])
            slab1data = s2f.getDFinMask(slab1data,maskdf)
            eventlist = s2f.getDFinMask(eventlist,maskdf)
        except:
            maskdf = pd.read_csv('%s_clp_%s.csv'%(thisfolder,date), sep=",", names=['lon','lat'])
            slab1data = s2f.getDFinMask(slab1data,maskdf)
            eventlist = s2f.getDFinMask(eventlist,maskdf)

        eventlist = s2f.getReferenceKagan(slab1data, eventlist, origorcentl, origorcentd)

        seismo_thick, taper_start, deplist, normpdfD, lendata = s2f.getSZthickness(eventlist,folder,slab,maxdep,maxdepdiff,origorcentl,origorcentd,slaborev,savedir,minlength)
        print ('slab, seismo_thick, ndata:',slab, seismo_thick, lendata)

        if lendata > 0:
            ax1z.plot(deplist, normpdfD,label='%s, s:%.1f, d:%.1f'%(slab,taper_start,seismo_thick),linewidth=2,c=clist[n])
            ax1z.plot([seismo_thick,seismo_thick],[0.09,0.1],linewidth=2,linestyle='dashed',c=clist[n])
            ax1z.plot([taper_start,taper_start],[0.09,0.1],linewidth=2,linestyle='dotted',c=clist[n]),
            n+=1

            slabdf.append(slab)
            shaldf.append(taper_start)
            deepdf.append(seismo_thick)
            peakdf.append(deplist[np.argmax(normpdfD)])
            nevsdf.append(lendata)

        try:
            interface = pd.read_csv('%s/%s_slab2_szt_%s.csv' % (savedir,slab,date))
        except:
            interface = pd.DataFrame()
        inter, upper, intra = s2f.orgEQs(interface,eventlist,maxdepdiff, seismo_thick, slab, maxdep)
        inter.to_csv('%s/%s_slab2_inter_%s.csv'%(savedir,slab,date),header=True,index=False,na_rep=np.nan)
        upper.to_csv('%s/%s_slab2_upper_%s.csv'%(savedir,slab,date),header=True,index=False,na_rep=np.nan)
        intra.to_csv('%s/%s_slab2_intra_%s.csv'%(savedir,slab,date),header=True,index=False,na_rep=np.nan)

        os.system('cp %s %s/inputdata'%(inFile,slabsbyfile))
        os.system('cp %s_dep_%s.grd %s/grids'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_str_%s.grd %s/grids'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_dip_%s.grd %s/grids'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_unc_%s.grd %s/grids'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_thk_%s.grd %s/grids'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_res_%s.csv %s/surfacetext'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_clp_%s.csv %s/clippingmasks'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_dat_%s.csv %s/filtereddata'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_nod_%s.csv %s/nodes'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_sup_%s.csv %s/supplement'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_par_%s.csv %s/parameters'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_figs_%s/* %s/maps'%(thisfolder,date,slabsbyfile))
        os.system('cp %s_xsecs_%s/* %s/crossections'%(thisfolder,date,slabsbyfile))

    ax1z.set_xlim([0,65])
    ax1z.legend(loc='best')
    ax1z.grid()
    ax1z.set_xlabel('Depths')
    ax1z.set_ylabel('P')
    ax1z.set_title('Slab depth distributions (surfacefilt = %i km, orig = %s, depth = %s, hist= %s)'%(maxdepdiff,origorcentl,origorcentd,slaborev))
    figtitle = '%s/szt/allpdf.png' % (slabsbyfile)
    figz.savefig(figtitle)
    plt.close()

    deetsdf = pd.DataFrame({'slab':slabdf,'shallow_lim':shaldf,'deep_lim':deepdf,'peak_depth':peakdf,'number_events':nevsdf})

    deetsdf = deetsdf[['slab','shallow_lim','deep_lim','peak_depth','number_events']]
    deetsdf.to_csv('%s/szt/table.csv' % (slabsbyfile),header=True,index=False,na_rep=np.nan,float_format='%0.1f')

    # fix overlapping slabs for catalogs
    ryufiltk = False
    ryufilti = False
    sulfilt = False
    phifilt = False
    cotfilt = False
    muefilt = False
    halhere = False
    japhere = False
    carhere = False

    for folder1 in filelist:
        (slab1,s2k,date1) = folder1.split('_')
        if slab1 == 'hal':
            halhere = True
            for folder2 in filelist:
                (slab2,s2k,date2) = folder2.split('_')
                if slab2 == 'sul':
                    s2f.sortoverlap(slab1, slab2, date1, date2, savedir)
                    sulfilt = True
                if slab2 == 'phi':
                    s2f.sortoverlap(slab1, slab2, date1, date2, savedir)
                    phifilt = True
                if slab2 == 'cot':
                    s2f.sortoverlap(slab1, slab2, date1, date2, savedir)
                    cotfilt = True
        if slab1 == 'car':
            carhere = True
            for folder2 in filelist:
                (slab2,s2k,date2) = folder2.split('_')
                if slab2 == 'mue':
                    s2f.sortoverlap(slab1, slab2, date1, date2, savedir)
                    muefilt = True
        if slab1 == 'kur' or slab1 == 'izu' or slab1 == 'jap':
            japhere = True
            for folder2 in filelist:
                (slab2,s2k,date2) = folder2.split('_')
                if slab2 == 'ryu':
                    s2f.sortoverlap(slab1, slab2, date1, date2, savedir)
                    if slab1 == 'kur':
                        ryufiltk = True
                    if slab1 == 'izu':
                        ryufilti = True
                    if slab1 == 'jap':
                        ryufiltk = True
                        ryufilti = True

    if not sulfilt and halhere:
        print ('halmahera and suluwesi catalogs were not filtered by eachother')
        print ('ensure that a slab model from either region is in filelist')
    if not phifilt and halhere:
        print ('halmahera and philippines catalogs were not filtered by eachother')
        print ('ensure that a slab model from either region is in filelist')
    if not cotfilt and halhere:
        print ('halmahera and cotobato catalogs were not filtered by eachother')
        print ('ensure that a slab model from either region is in filelist')
    if not muefilt and carhere:
        print ('caribbean and muertos catalogs were not filtered by eachother')
        print ('ensure that a slab model from either region is in filelist')
    if (not ryufiltk or not ryufilti) and japhere:
        print ('kur, or izu and ryukyu catalogs were not filtered by eachother')
        print ('ensure that a slab model from all three regions is in filelist')

# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        this can be used to move individual slab files (grids, parameters, 
        data, nodes, etc.) to a new file structure organized by file type 
        instead of by slab. This also calculates seismogenic zone thickness 
        for each slab using the Slab2 model.

        Required arguments include: 
            directory leading to where the original slab2 output folders are stored (-s slabsbyslab)
            a new directory to save the new file structure to (-f slabsbyfile)
            a directory listing the original input folders (-i inputfolder)
            the date of all of the input folders as MM-YY (-t inputdate)
            a flag indicating whether to use event origin or cmt origin for slab reference (-l origorcentl)
            a flag indicating whether to use event depth or cmt depth for depth histogram (-d origorcentd)
            a minimum length to use for 5th and 95th percentiles instead of 10th and 90th (-n minlength)
            depth distance around slab2 to filter events by (-m maxdepdiff)
            maximum depth to extend distribution to (-x maxdep)
            a flag indicating whether to make histogram of slab depths or event depths (-b slaborev)

        The list of slab folders/versions must be changed manually in the code.

        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', '--slabsbyslab', dest='slabsbyslab', type=str,
                        required=True, help='directory containing Slab2 output folders')

    parser.add_argument('-f', '--slabsbyfile', dest='slabsbyfile', type=str,
                        required=True, help='new directory to save file structure to')

    parser.add_argument('-i', '--inputfolder', dest='inputfolder', type=str,
                        required=True, help='directory containing Slab2 input files')
                        
    parser.add_argument('-t', '--inputdate', dest='inputdate', type=str,
                        required=True, help='date of input files (MM-YY)')

    parser.add_argument('-l', '--origorcentl', dest='origorcentl', type=str,
                        required=True, help='flag indicating origin (o) or cmt (c) lon lat for slab reference')

    parser.add_argument('-d', '--origorcentd', dest='origorcentd', type=str,
                        required=True, help='flag indicating origin (o) or cmt (c) depth for slab reference')

    parser.add_argument('-n', '--minlength', dest='minlength', type=int,
                        help='minimum length for 5th and 95th percentile calculations (optional)')

    parser.add_argument('-m', '--maxdepdiff', dest='maxdepdiff', type=int,
                        required=True, help='depth distance around slab2 to filter events by')

    parser.add_argument('-x', '--maxdep', dest='maxdep', type=int,
                        required=True, help='maximum depth to extend distribution to')

    parser.add_argument('-b','--slaborev', dest='slaborev', type=str,
                        required=True, help='make histogram of slab depths (s) or event depths (e)')

    pargs = parser.parse_args()

    main(pargs)
