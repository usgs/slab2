#!/usr/bin/env python

import argparse
import math
import os

import numpy as np
import pandas as pd


def main(args):

    # set command line variables
    folderdir = args.folderdir
    slab = args.slab
    date = args.date
    cint = args.cint
    cnam = args.colname
    if cnam == "strike":
        cnam = "str"
    if cnam == "depth":
        cnam = "dep"
    folder = "%s_slab2_%s" % (slab, date)
    supplement = pd.read_csv("%s/%s_slab2_sup_%s.csv" % (folderdir, slab, date))

    # empty contour file if exists for rewriting
    os.system("rm %s/%s_slab2_c%s%i_%s.txt" % (folderdir, slab, cnam, cint, date))

    # make array of contours and depths in file
    contourlist = [cint]
    thisc = cint
    while thisc < supplement[args.colname].max():
        thisc += cint
        contourlist.append(thisc)

    depthlist = np.array(list((set(supplement[args.colname].values))))
    depthlist = np.sort(depthlist)

    # identify value spacing
    n = 0
    sumd = 0
    for i in range(1, len(depthlist)):
        diff = depthlist[i] - depthlist[i - 1]
        sumd += diff
        n += 1

    maxdiff = math.ceil(sumd / n)

    # open contour file for writing
    with open(
        "%s/%s_slab2_c%s%i_%s.txt" % (folderdir, slab, cnam, cint, date), "a"
    ) as f:

        # loop through contours
        for c in contourlist:

            # make array of distances between contour and depths
            distdepths = np.abs(c - depthlist)

            if cnam == "dep":
                # if the depth is within range of depth spacing, make list of xy pts
                if min(distdepths) <= maxdiff:
                    supdep = depthlist[np.argmin(distdepths)]
                    dat = supplement[supplement.depth == supdep]

                    # if there are points at the depth near this contour, sort them
                    if len(dat) > 0:
                        if slab == "izu" or slab == "man" or slab == "ker":
                            dat = dat.sort_values(by=["lat"], ascending=False)
                        if slab == "sol" or slab == "hin" or slab == "pam":
                            dat = dat.sort_values(by=["lon"], ascending=False)

                        # write contours to contour file
                        f.write("> %i \n" % c)
                        dat = dat[["lon", "lat"]]
                        dat.to_csv(f, header=False, index=False, sep=" ")
                    else:
                        continue
                else:
                    continue
            else:
                dat = supplement[
                    (supplement[args.colname].values < c + maxdiff * 0.6)
                    & (supplement[args.colname].values > c - maxdiff * 0.6)
                ]
                if len(dat) > 0:
                    # write contours to contour file
                    f.write("> %i \n" % c)
                    dat = dat[["lon", "lat"]]
                    dat.to_csv(f, header=False, index=False, sep=" ")


# Help/description and command line argument parser
if __name__ == "__main__":
    desc = """
        Expects slab (-s), model-date(-d), path to output folder(-f),
        contour interval in km (-i), and column name to make contours of (-n)

        Writes file to [-f]/[-s]_slab2_c[-n][-i]_[-d].txt

        Example 1 input:
            python makecontours.py -s man -d 03.19.18 -f Output/man_slab2_03.19.18 -i 20 -n depth

        Example 1 output:
            Output/man_slab2_03.19.18/man_slab2_cdep20_03.19.18.txt

        Example 2 input:
            python makecontours.py -s man -d 03.19.18 -f Output/man_slab2_03.19.18 -i 5 -n dip

        Example 2 output:
            Output/man_slab2_03.19.18/man_slab2_cdip5_03.19.18.txt

        Depth contours can be plotted as a continuous line, every other contour
        specified by column name is designed to be plotted as individual points.

        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-s",
        "--slab",
        dest="slab",
        type=str,
        required=True,
        help="three letter slab code",
    )
    parser.add_argument(
        "-d",
        "--date",
        dest="date",
        type=str,
        required=True,
        help="date for model (MM.DD.YY)",
    )
    parser.add_argument(
        "-f",
        "--folderdir",
        dest="folderdir",
        type=str,
        required=True,
        help="directory/to/[slab]_slab2_[date] folder",
    )
    parser.add_argument(
        "-i",
        "--cint",
        dest="cint",
        type=int,
        required=True,
        help="contour interval (km)",
    )
    parser.add_argument(
        "-n",
        "--colname",
        type=str,
        required=True,
        help="column name to make contours of",
    )

    pargs = parser.parse_args()

    main(pargs)
