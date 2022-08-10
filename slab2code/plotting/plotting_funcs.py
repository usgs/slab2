import fiona
import pandas as pd
import pygmt
from shapely.geometry import Polygon, mapping

# Defining Slab names
all_slabs = ['cal', 'cam', 'cot','hin', 'man', 'sco', 'sul', 'sam', 'cas', 'him', 
'puy', 'mak', 'hal', 'kur', 'mue', 'alu', 'ryu', 'phi', 'ker',
 'van', 'png', 'car', 'hel', 'pam', 'sol', 'sum', 'izu']
tilted_slabs = ['izu', 'ker', 'man', 'sol'] 
no_trench = ['hal', 'hin', 'pam']

def get_info(grd: str, slab : str) -> tuple:
    """
    Helper function to collect relevant information from a grd file.

    Parameters
    ---------
    grd : str 
        Path to a grd file for reading
    slab : str
        Name of slab region
    Returns
    ---------
    region, z_scale : tuple of (list, list)
        region : list of [min_lon (float), max_lon (float), min_lat (float), max_lat (float)]
            Contains in order, rounded values of
            minimum longitude, maximum longitude, minimum latitude, maximum latitude
        z_scale : list of [min_z (float), max_z (float)]
            Contains in order, rounded values of
            minimum Z value, maximum Z value 
    """
    info = pygmt.grdinfo(grd).split()

    x_min = info[info.index('x_min:')+1]
    x_max = info[info.index('x_max:')+1]

    y_min = info[info.index('y_min:')+1]
    y_max = info[info.index('y_max:')+1]

    v_min = info[info.index('v_min:')+1]
    v_max = info[info.index('v_max:')+1]

    region = [round(float(x_min)), round(float(x_max)), round(float(y_min)), round(float(y_max))]
    z_scale = [float(v_min), float(v_max)]

    # If a slab has complex geometry, the region size is increased
    if slab in tilted_slabs: 
        region[0] += -3
        region[1] += 3
        region[2] += -3
        region[3] += 3

    return region, z_scale


def clip_grd(grid : str, clp, region : list, gtype : str, supp_dir : str, projection : str) -> str:
    """
    Function to clip slab grid files.

    When making plots of slab2 output grd files, the grd must be clipped 
    to the bounds of the slab.

    Parameters
    ----------
    grid : str
        path to the grid to clp
    clp : Pandas.DataFrame
        A dataframe of the clipping data, i.e. an outline of the slab
    region : list of [min_lon (float), max_lon (float), min_lat (float), max_lat (float)]
        Contains in order, rounded values of
        minimum longitude, maximum longitude, minimum latitude, maximum latitude
    gtype : str
        grid type i.e. dip, dep, thk, unc, str
    supp_dir : str
        path to the supplementary directory which holds temporary files used for plotting
    projection : str
        type of projection used in the plots

    Returns
    --------
    fname : str
        path to the clipped grid file
    """
    # Converting grd file to xyz format
    data = pygmt.grd2xyz(
        grid = grid,  
    )
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
    # Changing xyz back to grd, and saving clipped file in fname
    pygmt.xyz2grd(
        data = clipped,
        region = region,
        spacing = 0.05,
        projection = projection,
        outgrid = f'{supp_dir}/{gtype}.grd'
    )
    fname = f'{supp_dir}/{gtype}.grd'
    return fname

def main():
    exit()

if __name__ == '__main__':
    main()