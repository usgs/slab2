import argparse
import copy
import os
import shutil
import numpy as np
import pandas as pd
import pygmt
from plotting_funcs import *
from pygmt.clib import Session

# Defining Slab names
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
tilted_slabs = ["izu", "ker", "man", "sol"]
no_trench = ["hal", "hin", "pam"]

# Defining the command line arguments and their help statements
slab_help = "str : 3 letter slab code"
date_help = "str : Date of model output, in the form DD.MM.YY"
input_help = "str : Date of the input model used in the form MM-YY"
n_help = "int or str or list : Refer to Slab2Instructions for more details"
w_help = "str : Total width in km of the cross section. Defaulted to 1 if unspecified"
mul_help = "str : specify to make multiple cross sections based on what n is"
color_help = "str : specify to use a color blind friendly color pallette"
focal_help = "str : specify to plot focal mechanisms"

parser = argparse.ArgumentParser()
parser.add_argument("slab", type=str, help=slab_help)
parser.add_argument("date", type=str, help=date_help)
parser.add_argument("input_date", type=str, help=input_help)
parser.add_argument("n", help=n_help)
parser.add_argument("w", type=int, nargs="?", default=1, help=w_help)
parser.add_argument("--multiple", action="store_true", help=mul_help)
parser.add_argument("--on_map", action="store_true")
parser.add_argument("--color_blind", action="store_true", help=color_help)
parser.add_argument("--focal", action="store_true", help=focal_help)
args = parser.parse_args()

slab = args.slab
date = args.date
input_date = args.input_date
n = args.n

if n != "all" and "[" not in n:
    n = int(n)

# Divide total width by two to define width on either side of line
w = args.w / 2

# Projection to use in plots
projection = "M15c"

cwd = os.getcwd()
parent = os.path.join(cwd, os.pardir)
par_dir = os.path.abspath(parent)
# Path to output directory
path = f"{par_dir}/Output/{slab}_slab2_{date}"
# Path to trench, if slab has a trench
if slab not in no_trench:
    trench = "forplotting/trenches_usgs_2017_depths.csv"

# Defining paths to dat, clip, and input files
data = f"{path}/{slab}_slab2_dat_{date}.csv"
clip = f"{path}/{slab}_slab2_clp_{date}.csv"
infile = f"{par_dir}/Input/{input_date}/{slab}_{input_date}_input.csv"
# Reading input data
inputdf = pd.read_csv(infile, low_memory=False)
inputdf = inputdf.iloc[:, :]

bath = "forplotting/world_lr.grd"
ghayes_cpt = "forplotting/ghayes2.cpt"

# Defining the cmap to use
if args.color_blind:
    cmap = "forplotting/hawaii.cpt"
else:
    cmap = "seis"

# Creating the output directory
if not args.focal:
    xsec_dir = f"{path}/{slab}_{date}_xsec"
if args.focal:
    xsec_dir = f"{path}/{slab}_{date}_focal_xsec"
supp_dir = f"{path}/supplementary_Files"
if os.path.exists(xsec_dir) == False:
    os.mkdir(xsec_dir)
else:
    pass
if os.path.exists(supp_dir) == False:
    os.mkdir(supp_dir)
else:
    pass

# Preliminary implementation for focal mechanisms
def get_focal(data_type: str, lon: float, lat: float, az: float) -> pd.DataFrame:
    """
    Collects and projects focal mechanism data

     Parameters
    ------
    data_type : str
        type of data in input file to select
    lon : float
        longitude of the the starting point of the cross section line
    lat : float
        latitude of the the starting point of the cross section line
    az : float
        the azimuth value of the angle that is perpendicular to the trench at the cross section location
    Returns
    --------
    line2 : pandas.dataframe
        Data frame with all the relevant information to plot focal mechanisms
    """
    df = inputdf
    df = df.loc[df.etype == data_type]

    width = [-w, w]
    # Selecting relevant information from input file
    focals = df.iloc[:, [1, 0, 2, 12, 13, 14, 6]].copy()
    for i in range(len(df)):
        if focals.iloc[i, 0] < 0:  # add 360 if the longitude is negative
            focals.iloc[i, 0] += 360

    focals.depth = (focals.depth * -1).copy()
    focals.loc[len(focals.index)] = [212.5 - 360, 61.05, -25, 66, 85, 95, 9.2]
    # Trys to project the correct information
    try:
        line = pygmt.project(
            data=focals, center=[lon, lat], unit=True, width=width, azimuth=az
        )

        return line

    # If exception, return string and print error statement
    except:
        print(f"Not enough data for {data_type} to plot focal mechanisms")
        return "string"


def input_data(
    data_type: str, lon: float, lat: float, az: float, length: float
) -> pd.DataFrame:
    """
    Projects the data of a particular input type onto the cross section line of interest

    Parameters
    ------
    data_type : str
        type of data in input file to select
    lon : float
        longitude of the the starting point of the cross section line
    lat : float
        latitude of the the starting point of the cross section line
    az : float
        the azimuth value of the angle that is perpendicular to the trench at the cross section location
    length : float
        the length to constraint the projection to
    Returns
    --------
    line2 : pandas.dataframe
        a df representing the depth/distance of the data types in question for a particular line

    if exception is raised:
        line2 : str
            a string is returned for type comparison to allow for an empty df to not be plotted
    """
    # Collecting relevant input data
    df = inputdf
    df = df.loc[df.etype == data_type]
    width = [-w, w]

    df = df.iloc[:, [1, 0, 2]]
    df = df.reset_index(drop=True)

    # Correcting values in input data
    for i in range(len(df)):
        if df.iloc[i, 0] < 0:  # add 360 if the longitude is negative
            df.iloc[i, 0] += 360
        if df.iloc[i, 2] > 0:  # make the depth of the data negative if it is positive
            df.iloc[i, 2] = df.iloc[i, 2] * -1

    try:
        # Project the track on the line
        length = [0, length]

        line2 = pygmt.project(
            data=df,
            center=[lon, lat],
            unit=True,
            width=width,
            length=length,
            azimuth=az,
        )
        return line2

    except:
        print(f"Not enough data in this cross section to project {data_type}")
        line2 = "string"
        return line2


def index_info(i: int, df: pd.DataFrame) -> tuple:
    """
    Function to collect and print trench data

    Parameters
    ----------
    i : int
        Trench index to collect values at
    df : pd.DataFrame
        Trench data frame

    Returns
    -------
    tuple containing lon (float), lat (float), azimuth (float)

    """
    print("-" * 30)
    print(f"Creating XSEC plots for line {i}")
    # Finding azimuth of the angle that is perpendicular to the trench
    az = df.iloc[i, 2] + 90
    if az > 360:
        az = az - 360
    lon = df.iloc[i, 0]
    lat = df.iloc[i, 1]

    print(f"Longitude = {lon}")
    print(f"Latitude = {lat}")
    print(f"Azimuth = {az}")
    return lon, lat, az


def make_cpt(cmap: str, z: list) -> None:
    """
    Makes a CPT file for coloring
    Writes the cpt to cpt.cpt

    Parameters
    ----------
    cmap : str
        the code for the color palate to be used
    z : list of [min_z (float), max_z (float)]
        contains the min and max z values

    Returns
    -------
    None
    """
    pygmt.makecpt(
        cmap=cmap,
        series=f"{str(z[0])}/{str(z[1])}/10",
        background=True,
        color_model="r",
        continuous=True,
        reverse=True,
        output="cpt.cpt",
    )


def xsec(region: list, grd: str, z: list, az: float, lon: float, lat: float, n) -> None:
    """
    A function to plot cross sections and overview maps

    params
    ------
    region : list of [min_lon (float), max_lon (float), min_lat (float), max_lat (float)]
        the region of interest representing the min/max lon and min/max lat
    grd : str
        the path to the dep.grd file
    z : list of [min_z (float), max_z (float)]
        the min/max of the depth grid values
    az : float
        the azimuth value of the angle that is perpendicular to the trench at the cross section location

    lon : float
        longitude of the the starting point of the cross section line
    lat : float
        latitude of the the starting point of the cross section line
    n : int or list
        the number cross section currently being created
        corresponds to the row (index) of the trench file being used
        or a list containing the lat/lon/az values for a single cross section

    returns
    -------
    None

    Saves a pdf containing an overview map and plot of the cross section to a directory
    within the particular slab/date output directory
    """
    # Initial length
    length = 2000
    if slab == "izu" and 10 < lat < 15:
        length = 500

    # Axis labels
    xlabel = "Distance (Km)"
    ylabel = "Depth (Km)"
    # Plot title
    name = f"{slab} Cross Section Line {n}"
    # Reading clip file
    clip = f"{path}/{slab}_slab2_clp_{date}.csv"
    clp = pd.read_csv(clip)
    # Clipping grid
    grd = clip_grd(grd, clp, region, "dep", supp_dir, projection)
    # Read trench data, if slab has trench
    if slab not in no_trench:
        df = pd.read_csv(trench)
        df = df.loc[df.slab == slab]
        df = df.iloc[:, [0, 1]]

    # Projecting line perpendicular to the trench, with a starting point on trench
    # and a length of 2000
    line = pygmt.project(
        center=[lon, lat],
        azimuth=az,
        length=[-length, length],
        unit=True,
        generate="1",
    )

    # Changing lon values to degrees east
    for i in range(len(line.r)):
        if line.iloc[i, 0] < 0:
            line.iloc[i, 0] += 360

    # Constraining the line to be within the bounds of the slab
    line = pygmt.select(
        data=line.iloc[:, [0, 1]],
        F=f"{supp_dir}/poly.shp",
    )

    # Track the points on the line in the depth grd file
    track = pygmt.grdtrack(points=line, grid=grd, newcolname="tracks")
    track = track.dropna()

    # Redefining the starting lon/lat to the beginning of the xsec line if a lon lat az was
    # given for n
    if type(n) == type([]):
        lon = line.r[0]
        lat = line.s[0]

    # Project the track along a line perpendicular to the trench
    line2 = pygmt.project(data=track, center=[lon, lat], azimuth=az, unit=True)

    if slab in tilted_slabs:

        tilted = f"{path}/{slab}_slab2_sup_{date}.csv"  # Supplementary file
        tilted = pd.read_csv(tilted)
        tilted = tilted.iloc[:, :3]
        tilted.iloc[:, 2] = tilted.iloc[:, 2] * -1

        pw = 1  # counter for projection width increase.
        condition = None
        while condition == None and pw < 6:
            try:
                # Project tilted slab data
                line3 = pygmt.project(
                    data=tilted, center=[lon, lat], azimuth=az, unit=True, width=[0, pw]
                )

                condition = 0  # Will break the loop

                # The izu slab has length constraint issues in the first few trench indices
                # where the projection line re intersects slab. This will break the loop, as
                # interpolation is not needed in non-tilted sections, and the indices of length
                # constraint issues are outside of this region.
                if slab == "izu" and 10 < lat < 15:
                    break

                # Dropping last 10 rows of the line
                for i in range(len(line2) - 10, len(line2)):
                    line2 = line2.drop(i)
                # Filling the last rows with NaN
                empty = [np.nan for j in range(7)]
                for i in range(len(line2) - 1, len(line2) + 10):
                    line2.loc[len(line2)] = empty
                # Interpolating the last values of the normal slab to the tilted slab values
                line2 = pd.concat([line2, line3]).reset_index(drop=True)
                line2 = line2.interpolate(method="linear", limit_direction="forward")

            # Increase projection width if projection was empty
            except pd.errors.EmptyDataError:
                if pw == 1:
                    print("Issue With complex geometry")
                    print("Increasing tilted region projection width...")
                pw += 1
                if pw == 6:
                    print("Cross section not within bounds of tilted slab.")

    line2 = line2.iloc[:, [3, 2]]

    # Sometimes the lengths from the trench are negative. If so, set the zero point to the
    # first index shift all points accordingly
    if line2.iloc[0, 0] < 0:
        val = line2.iloc[0, 0] * -1
        for i in range(len(line2)):
            line2.iloc[i, 0] += val
    else:
        val = 0

    # Collecting extrema to define the axis sizes
    x1 = round((line2.iloc[:, 0]).min()) - 2
    x2 = round((line2.iloc[:, 0]).max()) + 5
    y1 = round((line2.iloc[:, 1]).min()) - 5
    x = max(abs(x2), abs(y1))
    x3 = x * -1
    if args.focal:
        # Square plot for focal mechanisms
        basemap_region = [-2, x, x3, 2]
        basemap_projection = "X8/8"
    else:
        basemap_region = [x1, x2, y1, z[1]]
        basemap_projection = "X16/6"
    if slab == "izu" and 10 < lat < 15:
        pass
    # else:
    length = (x2 - x1) + 20

    #################################
    # Code for XSEC Plotting Below  #
    #################################

    # Defines the input data types, and colors to be plotted with
    data_types = {
        "EQ": "skyblue",
        "AS": "grey",
        "ER": "hotpink",
        "BA": "red",
        "TO": "brown",
        "AA": "purple",
        "CP": "yellow",
        "RF": "green",
    }

    fig = pygmt.Figure()

    # Making basemap

    fig.basemap(
        region=basemap_region,
        projection=basemap_projection,
        frame=["WSne", "xaf+lDistance (Km)", "yaf+lDepth (Km)"],
    )
    # Plots the depth of the slab (surface)
    fig.plot(data=line2, pen="2p,blue", label="Slab")

    # Run through data types, collect the data for this particular xsec, and plot
    data_on_map = {}
    for key in data_types.keys():
        data = input_data(key, lon, lat, az, length)
        if type(data) != type(
            "string"
        ):  # input_data is set to return a str if there is no data along the xsec
            if args.on_map:
                data2 = copy.deepcopy(data)
                data_on_map[key] = data2.iloc[:, [0, 1]]
            data = data.iloc[:, [3, 2]]
            data.iloc[:, 0] += val

            if key != "BA":
                fig.plot(
                    data=data,
                    style="c.25c",
                    pen=".1,white",
                    color=data_types[key],
                    label=key,
                )
        else:
            pass

    # Plotting Focal mechanisms in the cross sectional plot below
    # DISCLAIMER - Preliminary and in development
    if args.focal:
        print(
            "DISCLAIMER -  Focal mechanism plotting within the cross sectional plot is\
 preliminary and in development"
        )
        count = 0
        focal_dict = {}
        for dat in ["EQ", "ER"]:
            # Projecting to get the correct data
            focal_data = get_focal(dat, lon, lat, az)
            focal_dict[dat] = focal_data

            if type(focal_data) != type("string"):
                for i in range(len(focal_data)):
                    if str(focal_data.iloc[i, 4]).lower() != "nan":
                        rake = float(focal_data.iloc[i, 5])
                        # Rotating rake angle to give cross sectional view
                        rake -= 90
                        focal_mechanism = dict(
                            strike=float(focal_data.iloc[i, 3]),
                            dip=float(focal_data.iloc[i, 4]),
                            rake=rake,
                            magnitude=float(focal_data.iloc[i, 6]),
                        )
                        fig.meca(
                            focal_mechanism,
                            scale=".25c",
                            longitude=float(focal_data.iloc[i, 7]),
                            latitude=float(focal_data.iloc[i, 2]),
                            depth=float(focal_data.iloc[i, 2]),
                        )

                        count += 1

                    else:
                        pass
                if count == 0:
                    print(f"No Focal Mechanism data for {dat} in this line")
                else:
                    pass
            else:
                pass

    fig.legend(box=True, position="JBL+jBL+o0.2c", region=basemap_region)

    fig.shift_origin(yshift=1)

    fig.text(text="A`", position="TR", font="28p,Helvetica-Bold,black")

    fig.text(text="A", position="Tl", font="28p,Helvetica-Bold,black")
    # Plotting Overview Map Below
    fig.grdimage(
        grid=bath,
        nan_transparent=True,
        region=region,
        projection=projection,
        cmap=ghayes_cpt,
        yshift=+11,
    )

    fig.coast(region=region, projection=projection, frame=["ag"], shorelines=True)

    if slab not in tilted_slabs:
        fig.grdimage(
            grid=grd,
            nan_transparent=True,
            region=region,
            transparency=20,
            projection=projection,
            cmap="cpt.cpt",
        )
    else:

        tilted = f"{path}/{slab}_slab2_sup_{date}.csv"  # Supplementary file
        tilted = pd.read_csv(tilted)
        tilted = tilted.iloc[:, :3]

        fig.grdimage(
            grid=grd,
            nan_transparent=True,
            region=region,
            transparency=20,
            projection=projection,
            cmap="cpt.cpt",
        )

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
            data=tilted,
            projection=projection,
            region=region,
            style="c0.08",
            cmap="cpt2.cpt",
        )
    # Plotting cross section line on map
    fig.plot(
        x=line.r,
        y=line.s,
        pen="2p,blue",
        region=region,
        projection=projection,
        frame=["ag", f"+t{name}"],
        label="Cross section",
    )
    # Adding markers on map if specified
    if args.on_map:
        for key in data_on_map.keys():
            data2 = data_on_map[key]
            fig.plot(
                data=data2,
                region=region,
                projection=projection,
                style="c.1c",
                pen=".1,black",
                color=data_types[key],
                label=key,
            )
    # Plotting Focal Mechanisms on map
    if args.focal:
        for dat in ["EQ", "ER"]:
            focal_data = focal_dict[dat]
            if type(focal_data) != type("string"):
                for i in range(len(focal_data)):
                    if str(focal_data.iloc[i, 4]) != "nan":
                        focal_mechanism = dict(
                            strike=float(focal_data.iloc[i, 3]),
                            dip=float(focal_data.iloc[i, 4]),
                            rake=float(focal_data.iloc[i, 5]),
                            magnitude=float(focal_data.iloc[i, 6]),
                        )
                        fig.meca(
                            focal_mechanism,
                            scale=".25c",
                            longitude=float(focal_data.iloc[i, 0]),
                            latitude=float(focal_data.iloc[i, 1]),
                            depth=float(focal_data.iloc[i, 2]),
                        )
                    else:
                        pass
            else:
                pass
    # Plotting clip to add outline
    fig.plot(data=clp, pen=".5,black", region=region, projection=projection)
    if slab not in no_trench:
        # Plotting the trench as a yellow line
        fig.plot(data=df, projection=projection, pen="2,yellow", label="Trench")
    # Adding contours
    fig.grdcontour(grid=grd, projection=projection)
    fig.colorbar(cmap="cpt.cpt", frame=[f"x+lDepth(Km)"])
    fig.legend(box=True)

    fig.text(
        x=[line.r[0], line.r[len(line) - 1]],
        y=[line.s[0], line.s[len(line) - 1]],
        text=["A", "A`"],
        font="10p,Helvetica-Bold,white",
        fill="black",
        transparency=10,
    )
    print(f"Saving plot in {xsec_dir}/{slab}_xsec_{n}.pdf")
    fig.savefig(f"{xsec_dir}/{slab}_xsec_{n}.pdf")


def main(n: int or str) -> None:

    file = f"{path}/{slab}_slab2_dep_{date}.grd"
    info = get_info(file, slab)
    region = info[0]
    z = info[1]
    # Creating depth CPT
    make_cpt(cmap, z)

    df = pd.read_csv(trench)
    df = df.loc[df.slab == slab]
    df = df.iloc[:, [0, 1, 2]]

    if n == "all":
        for i in range(0, len(df), 5):

            info = index_info(i, df)
            lon = info[0]
            lat = info[1]
            az = info[2]
            try:
                xsec(region, file, z, az, lon, lat, i)
            except:
                print("Could not make an XSEC, skipping this line.")
                pass
        exit()

    # Looping through trench rows
    count = 1
    adder = round(len(df) / n)
    i = 0

    while count <= n:
        try:
            info = index_info(i, df)
            lon = info[0]
            lat = info[1]
            az = info[2]
            condition = None
            while condition is None:
                try:
                    xsec(region, file, z, az, lon, lat, i)
                except:
                    print("Could not make an XSEC, skipping this line.")
                    i += 1
                    info = index_info(i, df)
                    lon = info[0]
                    lat = info[1]
                    az = info[2]
                    pass

                else:
                    i += adder
                    count += 1
                    condition = 0

        except IndexError:
            i = len(df) - 1
            info = index_info(i, df)
            lon = info[0]
            lat = info[1]
            az = info[2]
            count += n
            condition = None
            while condition is None:
                try:
                    xsec(region, file, z, az, lon, lat, i)
                except:
                    print("Could not make an XSEC, skipping this line.")
                    i -= 1
                    info = index_info(i, df)
                    lon = info[0]
                    lat = info[1]
                    az = info[2]
                    pass

                else:
                    condition = 0
    # Removing cpt files
    os.remove("cpt.cpt")
    if os.path.exists("cpt2.cpt"):
        os.remove("cpt2.cpt")


def single(n: int or list) -> None:
    """
    A function to create a single cross section
    The same function as main, but does not loop through trench rows

    params
    ------
    n : int
        the particular row of the trench file to create a cross section for

    returns
    -------
    None

    """

    file = f"{path}/{slab}_slab2_dep_{date}.grd"
    info = get_info(file, slab)
    region = info[0]
    z = info[1]
    make_cpt(cmap, z)

    if type(n) == type(1):
        try:
            df = pd.read_csv(trench)
            df = df.loc[df.slab == slab]
            df = df.iloc[:, [0, 1, 2]]
            info = index_info(n, df)
            lon = info[0]
            lat = info[1]
            az = info[2]
            xsec(region, file, z, az, lon, lat, n)

        except KeyError:  # IndexError:
            print(
                "The value you have chosen for the cross section number is out of range"
            )
            print(f"Please pick a value from 0-{len(df)}")
    else:
        n = n.strip("[]").split(",")
        lon = float(n[0])
        lat = float(n[1])
        az = float(n[2])
        if -180 <= lon <= 360 and -180 <= lat <= 180 and 0 <= az <= 360:
            n = [lon, lat, az]
            info = get_info(file, slab)
            region = info[0]
            z = info[1]
            xsec(region, file, z, az, lon, lat, n)

        else:
            print("Please input valid values for [lon,lat,az]")
            quit()

    # Removing cpt files
    os.remove("cpt.cpt")
    if os.path.exists("cpt2.cpt"):
        os.remove("cpt2.cpt")


if __name__ == "__main__":

    if args.multiple or n == "all":
        main(n)
    else:
        single(n)

"""
In order to use this file, 4 command line arguments are required:
1. 3 letter slab code (e.g. izu)
2. Date of model output in the form MM.DD.YY (e.g. 06.17.22)
3. The date corresponding to the database used for the model (e.g. 12-19)
4. N : int or str - The index of the trench or [lon,lat,az] to make the cross section at if single, 
                    number of cross sections to make, or 'all' to make multiple plots, 5 indicies apart

Optionally, Include all of the following arguments:
5. w : int - The Width of the cross section sampling. This is a positional argument.
             Default is 1, meaning the cross section will only sample along the line

6. --multiple - Argument to make multiple plots for a slab
                If specified, N will become the total number of cross sections to make, given that
                N < number of trench indices. If it is greater, cross sections for all indicies will be made.
                Specifying N as 'all' will make sections for every 5th trench index

7. --on_map - Argument to plot the data points on the overview map.

8. --color_blind - Creates the overview map with a colorblind sensitive CPT
The CPT used for this option comes from:
Crameri et al., 2020

9. --focal - Plots focal mechanisms

An example of what the command line should look like is below:

python map.py izu 06.17.22 12-19 10 --multiple
 - This will make 10 cross section plots for the izu output created on 06.17.22, using the 12-19 database
"""
