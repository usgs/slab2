import fiona
import os
import pandas as pd
import pygmt
from shapely.geometry import Polygon, mapping
from web_files import *
'''
Helper functions imported into xsec.ipynb and map.ipynb
'''
# Defining slabs
all_slabs = ['cal', 'cam', 'cot','hin', 'man', 'sco', 'sul', 'sam', 'cas', 'him', 
'puy', 'mak', 'hal', 'kur', 'mue', 'alu', 'ryu', 'phi', 'ker',
 'van', 'png', 'car', 'hel', 'pam', 'sol', 'sum', 'izu']

no_trench = ['hal', 'hin', 'pam']

tilted_slabs = ['izu', 'ker', 'man', 'sol'] # Slabs with a complex geometry which require extra work
map_issue = ['car', 'mue', 'sco', 'sam']

def get_info(grd: str, slab: str) -> tuple:
    """
    Helper function to collect relavent information from a grd file.
    Parameters
    ---------
    grd : str 
        Path to a grd file for reading
    Returns
    ---------
    region, z_scale : tuple of (list, list)
        region : list of (float, float, float, float)
            Contains in order, rounded values of
            minimum longitude, maximum longitude, minimum latitude, maximum latitude
        z_scale : list of (float, float)
            Contains in order, rounded values of
            minimum Z value, maximum Z value 
    """
    info = pygmt.grdinfo(grd)
    info = info.split()
    x_min = info[info.index('x_min:')+1]
    x_max = info[info.index('x_max:')+1]
    y_min = info[info.index('y_min:')+1]
    y_max = info[info.index('y_max:')+1]
    v_min = info[info.index('v_min:')+1]
    v_max = info[info.index('v_max:')+1]

    region = [round(float(x_min)), round(float(x_max)), round(float(y_min)), round(float(y_max))]
    z_scale = [float(v_min), float(v_max)]
    
    if slab in tilted_slabs: # if a slab has complex geometry, the region size is increased
        region[0] += -3
        region[1] += 3
        region[2] += -3
        region[3] += 3
    return region, z_scale    


def clip_grd(grid: str, clp: pd.DataFrame, region: list, gtype: str, date: str, slab: str, projection: str) -> str:
    """
    Function to clip the grids for slab files.
    When making plots of slab2 output grd files, the grd must be clipped 
    to the bounds of the slab.
    Parameters
    ----------
    grid : str
        path to the grid to clp
    clp : Pandas.DataFrame
        A dataframe of the clipping data, i.e. an outline of the slab
    region : list of [float,float]
        Contains in order, rounded values of
        minimum longitude, maximum longitude, minimum latitude, maximum latitude
    gtype : str
        grid type i.e. dip, dep, thk, unc, str
    Returns
    --------
    fname : str
        path to the clipped grid file
    """
    cwd = os.getcwd()
    parent = os.path.join(cwd, f'{os.pardir}/{os.pardir}')
    par_dir = os.path.abspath(parent)
    
    plot_par = os.path.join(cwd, os.pardir)
    plot_dir = os.path.abspath(plot_par)
    
    path = f"{par_dir}/Output/{slab}_slab2_{date}"
    supp_dir = f'{path}/supplementary_Files'
    out_dir = f'{path}/{slab}_{date}_maps'
    
    # Converting grd file to xyz format
    data = pygmt.grd2xyz(
        grid = grid,  
    ).dropna().reset_index(drop = True)

    if data.iloc[0,0] < 0:
        for i in range(len(data)):
            data.iloc[i,0] += 360
    # Looping through clp values and appending as tuple to list
    clp_vals = []
    for i in range(len(clp)):
        v1 = float(clp.iloc[i,0])
        v2 = float(clp.iloc[i,1])
        clp_vals.append((v1,v2))
    # Creating clip polygon
    clp_polygon = Polygon(clp_vals)
    
    schema = {
        'geometry': 'Polygon',
        'properties': {'id': 'float'},
    }
    # Write a new Shapefile
    with fiona.open(f'{supp_dir}/poly.shp', 'w', 'ESRI Shapefile', schema) as c:
        c.write({
            'geometry': mapping(clp_polygon),
            'properties': {'id': 123},
        })
    # Selecting grid xyz data that falls within polygon bounds
    clipped = pygmt.select(
        data = data,
        F = f'{supp_dir}/poly.shp',
    )

    # Changing xyz back to grd, and saving cliped file in fname
    grd = pygmt.xyz2grd(
        data = clipped,
        region = region,
        spacing = 0.05,
        projection = projection,
    )

    fname = f'{supp_dir}/{gtype}.grd'
    return grd


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

def map_input_data(data_type: str, input_date: str, slab : str) -> pd.DataFrame:
    '''
    Function to collect the input data of a slab for a given data type
    
    Parameters
    ----------
    data_type : str
        The type of input data to plot. e.g. TO, EQ, ER, BA, etc.
    input_date : str
        The date of the input file on github to use in the form MM-YY
    slab : str
        3 letter slab code
    Returns
    --------
    df : pd.DataFrame
        Data frame containing the input data to plot 
    '''
    df = get_input_df(input_date, slab)
    
    try:
        # Collecting relevant input data
        df = df.loc[df.etype == data_type]
        df = df.iloc[:,[0,1,2,6]]
        # Making depth negative
        df.iloc[:,2] *= -1
        # Changing to degrees east
        if df.iloc[0,1] < 0:
            df.iloc[:,1] += 360
            
    except IndexError:
        print(f'No {data_type} input data for {slab} from the {input_date} database')
        quit()
        
    return df  

def xsec_input_data(
    data_type: str, lon: float, lat: float, az: float, w: int, length: float, inputdf: pd.DataFrame, n : int or str or list, endpoint : list or None
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
    w : int
        The total width of the cross section in Km
    length : float
        the length to constraint the projection to
    inputdf : pd.DataFrame
        DataFrame containing the input data
    n : int or str or list
        The parameter which decides where to make the cross section
    endpoint: list or None
        With a list value for n, projections use an endpoint value, else they use an azimuth 
        
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
        if type(n) == type([]):
            length = [-length, length]
            line2 = pygmt.project(
            data=df,
            center=[lon, lat],
            unit=True,
            width=width,
            length=length,
            endpoint=endpoint,
        )
        else:
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


def get_focal(data_type: str, lon: float, lat: float, az: float, w: int, inputdf: pd.DataFrame) -> pd.DataFrame:
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
    w : int
        The total width of the cross section in Km
    
    Returns
    --------
    line2 : pd.DataFrame
        Data frame with all the relevant information to plot focal mechanisms
    inputdf : pd.DataFrame
        DataFrame containing the input data
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

    # If exeption, return string and print error statement
    except:
        print(f'Not enough data for {data_type} to plot focal mechanisms')
        return 'string'

def main():
    pass

if __name__ =="__main__":
    main()