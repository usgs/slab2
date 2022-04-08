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

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from copy import deepcopy
from pylab import arccos, argsort, cross, dot, double, eigh, pi, trace, zeros
from sklearn import mixture
from sklearn.metrics import mean_squared_error
import slab2functions as s2f


def main(args):

    model = args.model
    slabsbyfile = args.slabsbyfile
    slabsbyslab = args.slabsbyslab
    inputfolder = args.inputfolder
    inputdate = args.inputdate
    origorcentl = args.origorcentl
    origorcentd = args.origorcentd
    if args.box is not None:
        boxlonmin = args.box[0]
        boxlonmax = args.box[1]
        boxlatmin = args.box[2]
        boxlatmax = args.box[3]
    else:
        boxlonmin = -99999
        boxlonmax = 99999
        boxlatmin = -99999
        boxlatmax = 99999
    if args.minlength is None:
        minlength = 999999999
    else:
        minlength = args.minlength
    maxdep = args.maxdep
    maxdepdiff = args.maxdepdiff
    slaborev = args.slaborev

    pd.options.mode.chained_assignment = None
    warnings.filterwarnings("ignore", message="invalid value encountered in less")
    warnings.filterwarnings(
        "ignore", message="invalid value encountered in true_divide"
    )
    warnings.filterwarnings("ignore", message="invalid value encountered in greater")
    warnings.filterwarnings(
        "ignore", message="invalid value encountered in double_scalars"
    )

    # create new directory system for slab output
    os.system("mkdir %s" % slabsbyfile)
    printtest = False

    figz = plt.figure(figsize=(15, 10))
    ax1z = figz.add_subplot(111)
    n = 0

    slabdf = []
    shaldf = []
    deepdf = []
    peakdf = []
    nevsdf = []

    (slab, s2k, date) = model.split("_")

    inFile = "%s/%s_%s_input.csv" % (inputfolder, slab, inputdate)
    fname = slabsbyslab

    thisfolder = "%s/%s_slab2" % (slabsbyslab, slab)

    eventlistALL = pd.read_table(
        "%s" % inFile,
        sep=",",
        dtype={
            "lon": np.float64,
            "lat": np.float64,
            "depth": np.float64,
            "unc": np.float64,
            "etype": str,
            "ID": str,
            "mag": np.float64,
            "S1": np.float64,
            "D1": np.float64,
            "R1": np.float64,
            "S2": np.float64,
            "D2": np.float64,
            "R2": np.float64,
            "src": str,
            "time": str,
            "mlon": np.float64,
            "mlat": np.float64,
            "mdep": np.float64,
            "id_no": str,
        },
    )

    eventlistALL["ID"] = eventlistALL["id_no"].values

    ogcolumns = [
        "lat",
        "lon",
        "depth",
        "unc",
        "etype",
        "ID",
        "mag",
        "time",
        "S1",
        "D1",
        "R1",
        "S2",
        "D2",
        "R2",
        "src",
    ]
    kagancols = [
        "lat",
        "lon",
        "depth",
        "unc",
        "etype",
        "ID",
        "mag",
        "time",
        "S1",
        "D1",
        "R1",
        "S2",
        "D2",
        "R2",
        "src",
        "mlon",
        "mlat",
        "mdep",
    ]

    eventlist = eventlistALL[kagancols]

    depgrid = gmt.GMTGrid.load(fname)
    strgrid, dipgrid = s2f.mkSDgrd(depgrid)
    slab1data = s2f.mkSlabData(depgrid, strgrid, dipgrid, printtest)

    slab1data.loc[slab1data.lon < 0, "lon"] += 360
    eventlist.loc[eventlist.lon < 0, "lon"] += 360
    eventlist.loc[eventlist.mlon < 0, "mlon"] += 360

    if args.box is not None:
        if boxlonmin < 0:
            boxlonmin += 360
        if boxlonmax < 0:
            boxlonmax += 360
        slab1data = slab1data[slab1data.lon < boxlonmax]
        slab1data = slab1data[slab1data.lon > boxlonmin]
        slab1data = slab1data[slab1data.lat < boxlatmax]
        slab1data = slab1data[slab1data.lat > boxlatmin]

        if origorcentl == "c":
            eventlist = eventlist[eventlist.mlon < boxlonmax]
            eventlist = eventlist[eventlist.mlon > boxlonmin]
            eventlist = eventlist[eventlist.mlat < boxlatmax]
            eventlist = eventlist[eventlist.mlat > boxlatmin]
        else:
            eventlist = eventlist[eventlist.lon < boxlonmax]
            eventlist = eventlist[eventlist.lon > boxlonmin]
            eventlist = eventlist[eventlist.lat < boxlatmax]
            eventlist = eventlist[eventlist.lat > boxlatmin]

    try:
        maskdf = pd.read_csv(args.mask, delim_whitespace=True, names=["lon", "lat"])
        slab1data = s2f.getDFinMask(slab1data, maskdf)
        eventlist = s2f.getDFinMask(eventlist, maskdf)
    except:
        maskdf = pd.read_csv(args.mask, sep=",", names=["lon", "lat"])
        slab1data = s2f.getDFinMask(slab1data, maskdf)
        eventlist = s2f.getDFinMask(eventlist, maskdf)

    eventlist = s2f.getReferenceKagan(slab1data, eventlist, origorcentl, origorcentd)

    savedir = "%s" % slabsbyfile
    seismo_thick, taper_start, deplist, normpdfD, lendata = s2f.getSZthickness(
        eventlist,
        model,
        slab,
        maxdep,
        maxdepdiff,
        origorcentl,
        origorcentd,
        slaborev,
        savedir,
        minlength,
    )
    print("slab, seismo_thick:", slab, seismo_thick, lendata)

    try:
        interface = pd.read_csv("%s/%s_slab2_szt_%s.csv" % (savedir, slab, date))
    except:
        interface = pd.DataFrame()
    inter, upper, intra = s2f.orgEQs(
        interface, eventlist, maxdepdiff, seismo_thick, slab, maxdep
    )
    inter.to_csv(
        "%s/%s_slab2_inter_%s.csv" % (savedir, slab, date),
        header=True,
        index=False,
        na_rep=np.nan,
    )
    upper.to_csv(
        "%s/%s_slab2_upper_%s.csv" % (savedir, slab, date),
        header=True,
        index=False,
        na_rep=np.nan,
    )
    intra.to_csv(
        "%s/%s_slab2_intra_%s.csv" % (savedir, slab, date),
        header=True,
        index=False,
        na_rep=np.nan,
    )

    slabdf.append(slab)
    shaldf.append(taper_start)
    deepdf.append(seismo_thick)
    peakdf.append(deplist[np.argmax(normpdfD)])
    nevsdf.append(lendata)

    deetsdf = pd.DataFrame(
        {
            "slab": slabdf,
            "shallow_lim": shaldf,
            "deep_lim": deepdf,
            "peak_depth": peakdf,
            "number_events": nevsdf,
        }
    )

    deetsdf = deetsdf[
        ["slab", "shallow_lim", "deep_lim", "peak_depth", "number_events"]
    ]

    if args.box is not None or args.mask is not None:
        if args.box is not None:
            deetsdf.to_csv(
                "%s/%s_slab2_sztable_%s_%0.1f-%0.1f-%0.1f-%0.1f.csv"
                % (slabsbyfile, slab, date, boxlonmin, boxlonmax, boxlatmin, boxlatmax),
                header=True,
                index=False,
                na_rep=np.nan,
                float_format="%0.1f",
            )
        else:
            flist = args.mask.split("/")
            maskname = flist[-1]
            maskname = maskname[:-4]
            deetsdf.to_csv(
                "%s/%s_slab2_sztable_%s_%s.csv" % (slabsbyfile, slab, date, maskname),
                header=True,
                index=False,
                na_rep=np.nan,
                float_format="%0.1f",
            )
    else:
        deetsdf.to_csv(
            "%s/%s_slab2_sztable_%s.csv" % (slabsbyfile, slab, date),
            header=True,
            index=False,
            na_rep=np.nan,
            float_format="%0.1f",
        )


# Help/description and command line argument parser
if __name__ == "__main__":
    desc = """
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

        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-s",
        "--slabsbyslab",
        dest="slabsbyslab",
        type=str,
        required=True,
        help="directory containing Slab2 output folders",
    )

    parser.add_argument(
        "-f",
        "--slabsbyfile",
        dest="slabsbyfile",
        type=str,
        required=True,
        help="new directory to save output to",
    )

    parser.add_argument(
        "-i",
        "--inputfolder",
        dest="inputfolder",
        type=str,
        required=True,
        help="directory containing Slab2 input files",
    )

    parser.add_argument(
        "-r",
        "--model",
        dest="model",
        type=str,
        required=True,
        help="the slab model to use",
    )

    parser.add_argument(
        "-t",
        "--inputdate",
        dest="inputdate",
        type=str,
        required=True,
        help="date of input files (MM-YY)",
    )

    parser.add_argument(
        "-l",
        "--origorcentl",
        dest="origorcentl",
        type=str,
        required=True,
        help="flag indicating origin (o) or cmt (c) lon lat for slab reference",
    )

    parser.add_argument(
        "-d",
        "--origorcentd",
        dest="origorcentd",
        type=str,
        required=True,
        help="flag indicating origin (o) or cmt (c) depth for slab reference",
    )

    parser.add_argument(
        "-n",
        "--minlength",
        dest="minlength",
        type=int,
        help="minimum length for 5th and 95th percentile calculations (optional)",
    )

    parser.add_argument(
        "-m",
        "--maxdepdiff",
        dest="maxdepdiff",
        type=int,
        required=True,
        help="depth distance around slab2 to filter events by",
    )

    parser.add_argument(
        "-x",
        "--maxdep",
        dest="maxdep",
        type=int,
        required=True,
        help="maximum depth to extend distribution to",
    )

    parser.add_argument(
        "-b",
        "--slaborev",
        dest="slaborev",
        type=str,
        required=True,
        help="make histogram of slab depths (s) or event depths (e)",
    )

    parser.add_argument(
        "-o",
        "--box",
        metavar=("lonmin", "lonmax", "latmin", "latmax"),
        dest="box",
        type=float,
        nargs=4,
        help="bounding box [lonmin lonmax latmin latmax]",
    )

    parser.add_argument(
        "-k",
        "--mask",
        dest="mask",
        type=str,
        required=True,
        help="optional mask to calculate SZT within",
    )
    pargs = parser.parse_args()

    main(pargs)

