import argparse
import copy
import os
import pandas as pd
import pygmt
import shutil
from plotting_funcs import *

# Defining grd types in dictionary for ease of use and access
grd_dic = {
    "dep": "Depth",
    "dip": "Dip",
    "thk": "Thickness",
    "str": "Strike",
    "unc": "Unc",
}
# Defining all slab regions
all_slabs = [
    "cal",
    "cam",
    "cot",
    "hin",
    "man",
    "sco",
    "sul",
    "sam",
    "cas",
    "him",
    "puy",
    "mak",
    "hal",
    "kur",
    "mue",
    "alu",
    "ryu",
    "phi",
    "ker",
    "van",
    "png",
    "car",
    "hel",
    "pam",
    "sol",
    "sum",
    "izu",
]
# Defining slabs with no trench
no_trench = ["hal", "hin", "pam"]
# Slabs with a complex geometry
tilted_slabs = ["izu", "ker", "man", "sol"]
# Slabs that have issues clipping
clip_fail = ["car", "mue", "sco", "sam"]

# Defining command line arguments and their help statements
slab_help = "str : 3 letter slab code"
date_help = "str : Date of model output, in the form DD.MM.YY"
color_help = "str : specify to use a color blind friendly color pallette"
parser = argparse.ArgumentParser()
parser.add_argument("slab", type=str, help=slab_help)
parser.add_argument("date", type=str, help=date_help)
parser.add_argument("--color_blind", action="store_true", help=color_help)
args = parser.parse_args()

slab = args.slab
date = args.date

# Current Directory
cwd = os.getcwd()
parent = os.path.join(cwd, os.pardir)
par_dir = os.path.abspath(parent)
# Paths to files used in plotting
ghayes_cpt = "forplotting/ghayes2.cpt"
trench = "forplotting/trenches_usgs_2017_depths.csv"
bath = "forplotting/world_lr.grd"
# Defining a perspective and projection in the global scope for ease
perspective = [180, 90]
projection = "M15c"
# Path to output directory
path = f"{par_dir}/Output/{slab}_slab2_{date}"
supp_dir = f"{path}/supplementary_Files"


def plot(file: str, region: list, z: list, out_dir: str) -> None:
    """
    Function to create 2d plots of Slab2 models.

    This function will populate a directory within a Slab2 Output directory
    with pdf files of the plot type requested

    Parameters
    ----------
    file : str
        path to a grd file to plot
    region : list of (float/int, float/int)
        The region of interest to plot
    z : list of (float, float)
        The minimum and maximum values in the grd file
    out_dir : str
        The path to the new directory for the plots to be saved to
    Returns
    ----------
    None
    """
    # Defining grid type
    grd_type = file[-16:-13]
    # Path to clipping file
    clip = f"{path}/{slab}_slab2_clp_{date}.csv"
    clp = pd.read_csv(clip)
    # Clipping grid file
    if slab not in clip_fail:
        file = clip_grd(file, clp, region, grd_type, supp_dir, projection)
    # Reading trench data
    trenchdf = pd.read_csv(trench)
    trenchdf = trenchdf.loc[trenchdf.slab == slab]
    trenchdf = trenchdf.iloc[:, [0, 1]]
    # Create another instance of the file object
    temp_file = copy.deepcopy(file)

    grd = grd_dic[grd_type]
    name = f"{grd} of {slab} {date}"
    # Name of file to store pdf in
    fname = f"{out_dir}/{slab}_{date}_{grd_type}.pdf"

    print(f"Creating plot for {name}")

    fig = pygmt.Figure()

    # Plotting world topography below
    fig.grdimage(
        grid=bath,
        nan_transparent=True,
        region=region,
        perspective=perspective,
        projection=projection,
        cmap=ghayes_cpt,
    )

    # Creating coast outlines
    fig.coast(
        region=region,
        shorelines=True,
        projection=projection,
        frame=["ag", f'+t"{name}"'],
        perspective=perspective,
    )
    # Defining color palettes and units
    cmap_dict = {
        "dep": ["seis", "Km"],
        "str": ["no_green", "Deg"],
        "dip": ["no_green", "Deg"],
        "thk": ["no_green", "Km"],
        "unc": ["gray"],
    }
    # Creating cpt below
    if args.color_blind:
        cmap = f"forplotting/hawaii.cpt"
    else:
        cmap = cmap_dict[grd_type][0]

    pygmt.grd2cpt(grid=file, cmap=cmap, reverse=True, continuous=True, output="cpt.cpt")

    # Plotting the grd below with the cpt defined above
    fig.grdimage(
        grid=file,
        nan_transparent=True,
        region=region,
        transparency=20,
        perspective=perspective,
        projection=projection,
        cmap="cpt.cpt",
    )

    # Creating cpts and plotting the areas of complex geometries
    if slab in tilted_slabs:
        tilted = f"{path}/{slab}_slab2_sup_{date}.csv"  # Supplementary file
        index_dict = {"dep": 2, "str": 3, "dip": 4, "unc": 5, "thk": 8}

        df = pd.read_csv(tilted)
        df = df.iloc[
            :, [0, 1, index_dict[grd_type]]
        ]  # Collecting relevant supplementary information

        if grd_type == "dep":
            cpt = pygmt.makecpt(
                cmap=cmap,
                series=f"{str(-1*z[1])}/{str(-1*z[0])}/10",
                background=True,
                color_model="r",
                continuous=True,
                output="cpt2.cpt",
            )

            # Plotting complex dep surface
            fig.plot(
                data=df,
                projection=projection,
                region=region,
                style="c0.08",
                cmap="cpt2.cpt",
            )

        else:
            # Plotting other grd type surfaces
            fig.plot(
                data=df,
                projection=projection,
                region=region,
                style="c0.08",
                cmap="cpt.cpt",
            )

    # Drawing contours
    fig.grdcontour(grid=file, perspective=perspective, projection=projection)

    # Plotting Slab outline
    fig.plot(data=clp, pen=".5,black", region=region, projection=projection)

    # Plotting Trench
    if slab not in no_trench:
        fig.plot(data=trenchdf, projection=projection, pen="2,yellow", label="Trench")
        fig.legend(box=True)

    # Creating colorbar
    if grd_type != "unc":
        frame = [f"x+l{grd} ({cmap_dict[grd_type][1]})"]
    else:
        frame = [f"x+l{grd}"]
    fig.colorbar(
        cmap="cpt.cpt",
        frame=frame,
    )

    # Saves figure
    print(f"Saving Plot in : {fname}")
    print("-" * 40)
    fig.savefig(fname=fname)
    os.remove("cpt.cpt")


def main():
    try:
        print("-" * 30)
        if slab not in all_slabs:

            print("Please input a valid slab code")
            exit()

        # Defining directories to store plots
        out_dir = f"{path}/{slab}_{date}_maps"

        # Checking if directory exists, if so, remove and recreate
        if os.path.exists(out_dir) == False:
            os.mkdir(out_dir)
        else:
            shutil.rmtree(out_dir)
            os.mkdir(out_dir)
        if os.path.exists(supp_dir) == False:
            os.mkdir(supp_dir)
        else:
            pass

        # Loop through all grd files and make plots
        for dir, sub, files in os.walk(path):
            for file in files:
                if ".grd" in file and "slab2" in file:
                    f = f"{path}/{file}"
                    info = get_info(f, slab)
                    region = info[0]
                    z = info[1]
                    plot(f, region, z, out_dir)

    except FileNotFoundError:
        print("File Not Found!")
        print("Please format the command line arguments correctly")
        print("e.g. python map.py alu 06.14.22")
        print("-" * 30)


if __name__ == "__main__":
    main()
"""
In order to use this file, 2 command line arguments are necessary:
1. 3 letter slab code (e.g. alu)
2. Date of model output in the form MM.DD.YY (e.g. 06.14.22)

Optionally append:
--color_blind to use a color blind friendly color palette

An example of what the command line should look like is below:

python map.py alu 06.14.22 --color_blind
"""
