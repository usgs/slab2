import argparse
import os
from functools import partial

import numpy as np
import pandas as pd
from multiprocess import Pool


def calls2d(database, filedate, slablist, i):
    slab = slablist[i]
    os.system(
        "python s2d.py -p %s -d %s -f %s_%s_input.csv"
        % (slab, database, slab, filedate)
    )


def main(args):

    filedate = args.filedate
    database = args.database

    slablist = [
        "alu",
        "cal",
        "cam",
        "car",
        "cas",
        "cot",
        "hal",
        "hel",
        "him",
        "hin",
        "izu",
        "jap",
        "ker",
        "kur",
        "mak",
        "man",
        "mue",
        "pam",
        "png",
        "phi",
        "puy",
        "ryu",
        "sam",
        "sco",
        "sol",
        "sul",
        "sum",
        "van",
    ]

    indices = range(len(slablist))
    pool1 = Pool(args.nCores)
    partial_loop1 = partial(calls2d, database, filedate, slablist)

    pts = pool1.map(partial_loop1, indices)
    pool1.close()
    pool1.join()


# Help/description and command line argument parser
if __name__ == "__main__":
    desc = """
        Expected slab regions include: 

        Aleutians               alu
        Calabria                cal
        Central America         cam
        Caribbean               car
        Cascadia                cas
        Cotabato                cot
        Halmahera               hal
        Hellenic                hel
        Himalaya                him
        Hindu Kush              hin
        Izu-Bonin               izu
        Kermadec                ker
        Kuril                   kur
        Makran                  mak
        Manila                  man
        Muertos                 mue
        Pamir                   pam
        New Guinea              png
        Philippines             phi
        Puysegur                puy
        Ryukyu                  ryu
        South America           sam
        Scotia                  sco
        Solomon Islands         sol
        Sulawesi                sul
        Sumatra/Java            sum
        Vanuatu                 van
        
        This script will make an input file for each slab region

        """
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-d",
        "--database",
        dest="database",
        type=str,
        required=True,
        help="file listing slab parameters",
    )

    parser.add_argument(
        "-f",
        "--filedate",
        dest="filedate",
        type=str,
        required=True,
        help="MM-YY of most recent database update",
    )

    parser.add_argument(
        "-c",
        "--nCores",
        dest="nCores",
        type=int,
        required=True,
        help="number of cores to run loop over",
    )

    pargs = parser.parse_args()
    main(pargs)
