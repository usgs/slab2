#!/usr/bin/env python

import argparse
import cProfile
import math
import os.path
import warnings
from datetime import datetime
from functools import partial

import mapio.gmt as gmt
import matplotlib
import numpy as np
import pandas as pd
import psutil
from multiprocess import Pool
from pandas import DataFrame
from scipy import ndimage

import loops as loops
import slab2functions as s2f

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(args):

    trenches = "library/misc/trenches_usgs_2017_depths.csv"
    agesFile = "library/misc/interp_age.3.2g.nc"
    ageerrorsFile = "library/misc/interp_ageerror.3.2g.nc"
    polygonFile = "library/misc/slab_polygons.txt"
    parFile = args.parFile
    pd.options.mode.chained_assignment = None
    warnings.filterwarnings("ignore", message="invalid value encountered in less")
    warnings.filterwarnings(
        "ignore", message="invalid value encountered in true_divide"
    )
    warnings.filterwarnings("ignore", message="invalid value encountered in greater")
    warnings.filterwarnings(
        "ignore", message="invalid value encountered in double_scalars"
    )

    date = datetime.today().strftime("%m.%d.%y")
    now = datetime.now()
    time = "%s.%s" % (now.hour, now.minute)

    for line in open(parFile):
        plist = line.split()
        if len(plist) > 2:
            if plist[0] == "inFile":
                inFile = plist[2]
            if plist[0] == "use_box":
                use_box = plist[2]
            if plist[0] == "latmin":
                latmin = np.float64(plist[2])
            if plist[0] == "latmax":
                latmax = np.float64(plist[2])
            if plist[0] == "lonmin":
                lonmin = np.float64(plist[2])
            if plist[0] == "lonmax":
                lonmax = np.float64(plist[2])
            if plist[0] == "slab":
                slab = plist[2]
            if plist[0] == "grid":
                grid = np.float64(plist[2])
            if plist[0] == "radius1":
                radius1 = np.float64(plist[2])
            if plist[0] == "radius2":
                radius2 = np.float64(plist[2])
            if plist[0] == "sdr":
                sdr = np.float64(plist[2])
            if plist[0] == "ddr":
                ddr = np.float64(plist[2])
            if plist[0] == "taper":
                taper = np.float64(plist[2])
            if plist[0] == "T":
                T = np.float64(plist[2])
            if plist[0] == "node":
                node = np.float64(plist[2])
            if plist[0] == "filt":
                filt = np.float64(plist[2])
            if plist[0] == "maxdist":
                maxdist = np.float64(plist[2])
            if plist[0] == "minunc":
                minunc = np.float64(plist[2])
            if plist[0] == "mindip":
                mindip = np.float64(plist[2])
            if plist[0] == "minstk":
                minstk = np.float64(plist[2])
            if plist[0] == "maxthickness":
                maxthickness = np.float64(plist[2])
            if plist[0] == "seismo_thick":
                seismo_thick = np.float64(plist[2])
            if plist[0] == "dipthresh":
                dipthresh = np.float64(plist[2])
            if plist[0] == "fracS":
                fracS = np.float64(plist[2])
            if plist[0] == "kdeg":
                kdeg = np.float64(plist[2])
            if plist[0] == "knot_no":
                knot_no = np.float64(plist[2])
            if plist[0] == "rbfs":
                rbfs = np.float64(plist[2])

    TR_data = pd.read_csv(trenches)
    TR_data = TR_data[TR_data.slab == slab]

    fullfolder = args.folder
    dirlist = fullfolder.split("/")
    folder = dirlist[-1]
    (slab, slk, date) = folder.split("_")

    dataornodes = "nodes"
    shift_out = pd.read_csv("%s/%s_slab2_nod_%s.csv" % (fullfolder, slab, date))
    used_data = pd.read_csv("%s/%s_slab2_dat_%s.csv" % (fullfolder, slab, date))
    clip = pd.read_csv(
        "%s/%s_slab2_clp_%s.csv" % (fullfolder, slab, date), names=["lon", "lat"]
    )  # this way splits the csv file into 2 columns with names lon, lat - KLH 09/18/2019

    npass = 1

    print("    Creating surfaces...")

    surfdata = np.zeros((len(shift_out), 4))
    if dataornodes == "nodes":
        surfdata[:, 0], surfdata[:, 1], surfdata[:, 2], surfdata[:, 3] = (
            shift_out["lon"].values,
            shift_out["lat"].values,
            shift_out["depth"].values,
            shift_out["stdv"].values,
        )
    elif dataornodes == "data":
        surfdata[:, 0], surfdata[:, 1], surfdata[:, 2], surfdata[:, 3] = (
            shift_out["bzlon"].values,
            shift_out["bzlat"].values,
            shift_out["psdepth"].values,
            shift_out["stdv"].values,
        )

    errordata = np.zeros((len(shift_out), 4))
    errordata[:, 0], errordata[:, 1], errordata[:, 2], errordata[:, 3] = (
        shift_out["lon"].values,
        shift_out["lat"].values,
        shift_out["stdv"].values,
        np.ones(len(shift_out)),
    )

    errordataB = np.zeros((len(shift_out), 4))
    errordataB[:, 0], errordataB[:, 1], errordataB[:, 2], errordataB[:, 3] = (
        shift_out["lon"].values,
        shift_out["lat"].values,
        shift_out["shiftstd"].values,
        np.ones(len(shift_out)),
    )

    thickdata = np.zeros((len(shift_out), 4))
    thickdata[:, 0], thickdata[:, 1], thickdata[:, 2], thickdata[:, 3] = (
        shift_out["lon"].values,
        shift_out["lat"].values,
        shift_out["thickness"].values,
        np.ones(len(shift_out)),
    )

    filters = np.logspace(-0.2, 0.5, num=30)
    smoothers = [rbfs]
    meanBA = 5

    for b in range(len(smoothers)):
        misfits = []
        theoutput = []
        objfvs = []
        objfvmin = 0
        for i in range(len(filters)):
            filt = round(1 / filters[i], 4)
            print("filter:", filt, node, grid, meanBA, kdeg, knot_no, rbfs)

            if slab == "sum":
                Surfgrid, xi, dl = s2f.chunksurface(
                    surfdata,
                    node,
                    T,
                    slab,
                    grid,
                    "depth",
                    time,
                    "test.txt",
                    filt,
                    pd.DataFrame(),
                    npass,
                    TR_data,
                    meanBA,
                    kdeg,
                    knot_no,
                    rbfs,
                    shift_out,
                    "fin",
                    "og",
                    "lon",
                    100,
                    110,
                    105,
                )
                flipornot = "flip"
            elif slab == "jap":
                Surfgrid, xi, dl = s2f.chunksurface(
                    surfdata,
                    node,
                    T,
                    slab,
                    grid,
                    "depth",
                    time,
                    "test.txt",
                    filt,
                    pd.DataFrame(),
                    npass,
                    TR_data,
                    meanBA,
                    kdeg,
                    knot_no,
                    rbfs,
                    shift_out,
                    "fin",
                    "og",
                    "lat",
                    30,
                    40,
                    35,
                )
                flipornot = "flip"
            else:
                Surfgrid, xi, dl = s2f.pySurface3(
                    surfdata,
                    node,
                    T,
                    slab,
                    grid,
                    "depth",
                    time,
                    "test.txt",
                    filt,
                    pd.DataFrame(),
                    npass,
                    TR_data,
                    meanBA,
                    kdeg,
                    knot_no,
                    rbfs,
                    shift_out,
                    "fin",
                    "og",
                )
                flipornot = "dontflip"

            sigma = (filt / 2.0) / node

            Filtgrid = ndimage.filters.gaussian_filter(Surfgrid, sigma, mode="reflect")

            # Create output array
            output = np.zeros([len(xi), 10]) * np.nan

            output[:, 0] = xi[:, 0]  # lon Longitude at node (not shifted)
            output[:, 1] = xi[:, 1]  # lat Latitude at node
            output[
                :, 3
            ] = (
                Filtgrid.flatten()
            )  # dep_shift_smooth Post-shift surface depth after smoothing
            output[:, 0][output[:, 0] < 0] += 360

            # Want to make lon in the clip file positive (had to change how clip file was read in above):
            clip.loc[(clip["lon"] < 0), ["lon"]] += 360  # KLH 09/18/2019

            output[:, 2][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 3][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 4][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 5][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 6][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 7][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 8][output[:, 3] > shift_out["depth"].max()] = np.nan
            output[:, 9][output[:, 3] > shift_out["depth"].max()] = np.nan

            if dataornodes == "nodes":
                addextra = shift_out[shift_out.depth < 35] * 1.0
                shift_out2 = pd.concat(
                    [shift_out, addextra, addextra, addextra], sort=True
                )
                misfit, datafit = s2f.nodeMisfit(shift_out, output, clip)

            elif dataornodes == "data":
                doubleweight = used_data[
                    (used_data.etype == "AA")
                    | (used_data.etype == "AS")
                    | (used_data.etype == "BA")
                    | (used_data.etype == "RF")
                ]
                used_data2 = pd.concat(
                    [
                        used_data,
                        doubleweight,
                        doubleweight,
                        doubleweight,
                        doubleweight,
                        doubleweight,
                    ],
                    sort=True,
                )
                misfit, datafit = s2f.nodeMisfit(used_data2, output, clip)

            print(
                "----------------- kdeg, knot_no, rbfs, filt, filters[i], datafilt (%i of %i loops)"
                % (i, len(filters)),
                kdeg,
                knot_no,
                rbfs,
                filt,
                filters[i],
                datafit,
            )
            surf = str(filt)

            s2f.histofdiffs(misfit, knot_no, rbfs, filt, kdeg, slab, fullfolder, date)
            objfv = math.sqrt(datafit * datafit + filters[i] * filters[i])
            misfits.append(datafit)
            objfvs.append(objfv)

            if objfvmin != 0:
                if objfv < objfvmin:
                    finaloutput = np.copy(output)
                    minfilt = filt
                    objfvmin = objfv
            else:
                objfvmin = objfv
                minfilt = filt
                finaloutput = np.copy(output)

        filtfitdf = pd.DataFrame({"filt": filters, "dfit": misfits, "objf": objfvs})
        bestfilt, filtfitdf = s2f.plotLcurve(
            filtfitdf, ("%s/%s_slab2_lcv_%s.png" % (fullfolder, slab, date))
        )
        filtfitdf.to_csv(
            "%s/%s_slab2_lcv_%s.csv" % (fullfolder, slab, date),
            header=True,
            index=False,
            na_rep=np.nan,
            float_format="%.4f",
        )
        print("bestfilt", bestfilt)


# Help/description and command line argument parser
if __name__ == "__main__":
    desc = """
        Expected slab regions include: 
            Alaska-Aleutians         > alu
            Central America          > mex
            Cascadia                 > cas
            Izu-Bonin                > izu
            Kermadec-Tonga           > ker
            Kamchatka/Kurils/Japan   > kur
            Philippines              > phi
            Ryukyu                   > ryu
            Santa Cruz Islands       > van
            Scotia                   > sco
            Solomon Islands          > sol
            South America            > sam
            Sumatra-Java             > sum
        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-p",
        "--parFile",
        dest="parFile",
        type=str,
        required=True,
        help="file listing slab parameters",
    )
    parser.add_argument(
        "-f",
        "--folder",
        dest="folder",
        type=str,
        required=True,
        help="directory containing model to be L-Curved",
    )
    parser.add_argument(
        "-t",
        "--test",
        metavar=("lonmin", "lonmax", "latmin", "latmax"),
        dest="test",
        type=float,
        nargs=4,
        help="test box [lonmin lonmax latmin latmax]",
    )
    pargs = parser.parse_args()

    main(pargs)

