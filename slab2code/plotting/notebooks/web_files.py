import copy
import requests
from bs4 import BeautifulSoup as bs
import os
import pandas as pd
import pygmt
import io
'''
Helper functions imported into map.ipynb to webscrape ScienceBase (https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467)
to collect the proper grid files to make plots of maps.
Files from the GitHub (https://github.com/usgs/slab2) are also used.
'''
# Defining slabs
all_slabs = ['cal', 'cam', 'cot','hin', 'man', 'sco', 'sul', 'sam', 'cas', 'him', 
'puy', 'mak', 'hal', 'kur', 'mue', 'alu', 'ryu', 'phi', 'ker',
 'van', 'png', 'car', 'hel', 'pam', 'sol', 'sum', 'izu']

no_trench = ['hal', 'hin', 'pam']

tilted_slabs = ['izu', 'ker', 'man', 'sol']
map_issue = ['car', 'mue', 'sco', 'sam']


def paths(slab: str, nb_type: str) -> tuple:
    '''
    Function to create global paths for ease of use
    
    Parameters
    ----------
    slab : str
        3 letter slab code
    nb_type : str
        the type of notebook used, either map or xsec
    
    Returns 
    --------
    tuple of str, str, str, str, str
        contains the date, path to output directory, path to plotting directory within
        slab2code, path to model output directory, and path to supplementary directory
    '''
    # Defining paths
    global date
    date = get_date(slab)
    global cwd
    cwd = os.getcwd()
    parent = os.path.join(cwd, f'{os.pardir}/{os.pardir}')
    # Path to slab2/slab2code
    par_dir = os.path.abspath(parent)
    global plot_dir
    plot_par = os.path.join(cwd, os.pardir)

    # Path to slab2/slab2code/plotting
    plot_dir = os.path.abspath(plot_par)
    global path    
    path = f"{par_dir}/Output/{slab}_slab2_{date}"
    global supp_dir
    supp_dir = f'{path}/supplementary_Files'
    global out_dir
    if 'map' in nb_type:
        out_dir = f'{path}/{slab}_{date}_maps'
    else:
        out_dir = f'{path}/{slab}_{date}_xsec'

    # Creating output directories if they do not yet exist
    if os.path.exists(path) == False: 
        os.mkdir(path)
    if os.path.exists(supp_dir) == False:
        os.mkdir(supp_dir)
    if os.path.exists(out_dir) == False:
        os.mkdir(out_dir)
    
    return date, path, plot_dir, out_dir, supp_dir

def collect_slab_links() -> dict:
    '''
    Function to collect links to files on ScienceBase

    Parameters
    -----------
    None

    Returns
    --------
    slab_link_dic : dict
        Dictionary containing links to relevant files for plotting
    ''' 
    # ScienceBase url (page 1)
    url = 'https://www.sciencebase.gov/catalog/items?parentId=5aa1b00ee4b0b1c392e86467'
    r = requests.get(url)
    soup = bs(r.content,'html.parser')
    # Find all links on page 1
    p1_links = soup.find_all('a')
    # ScienceBase url (page 2)
    url = 'https://www.sciencebase.gov/catalog/items?parentId=5aa1b00ee4b0b1c392e86467&offset=20&max=20'
    r = requests.get(url)
    soup = bs(r.content,'html.parser')
    # Find all links on page 2
    p2_links = soup.find_all('a')
    
    slab_links = []
    # Finding slab2 links
    for link in p1_links:
        name = str(link.text).lower()
        if 'slab2' in name:
            slab_links.append(link)
            
    for link in p2_links:
        name = str(link.text).lower()
        if 'slab2' in name:
            slab_links.append(link)
          
    slab_link_dic = {}   
    # Separating links by slab
    for link in slab_links:
        slab = str(link.text).lower().split()[8][:3]
        # Slabs below have a 3 letter code that is not within the name of the region given on ScienceBase
        if slab == 'ala': slab = 'alu'
        if slab == 'cen': slab = 'cam'
        if slab == 'new': slab = 'png'
        if slab == 'sou': slab = 'sam'
        if slab == 'kam': slab = 'kur'
        # Adding links to dictionary, where key is slab code
        data_link = link['href']
        slab_link_dic[slab] = f'{url[:27]}{data_link}'
        
    return slab_link_dic


def get_date(slab: str) -> str:
    '''
    Function to get the output date of the files on ScienceBase

    Parameters
    ----------
    slab : str
        3 letter slab code
    
    Returns
    ---------
    date : str
        The date of the Slab2 output on ScienceBase
    '''
    # Collecting links
    links = collect_slab_links()
    url = links[slab]
    r = requests.get(url)
    # Collecting all links+text
    soup = bs(r.content,'html.parser')
    all_links = soup.find_all('span',{'class':"sb-file-get sb-download-link"})    
    f_type = str(all_links[3].text)
    # Finding date
    date = f_type[14:22]
            
    return date


def grid_links(path : str, date : str, slab : str) -> list:
    '''
    Function to write files locally and return paths to files

    Parameters
    ----------
    path : str
        Path to output directory where files will be stored
    date : str
        date of the output model
    slab : str
        3 letter slab code
    
    Returns
    --------
    files : list of [grid file list (list), clipping file (str), tilted file (str)]
        The paths to the relevant plotting files
    '''
    # Collect all links
    links = collect_slab_links()
    # Collect slab specific links
    url = links[slab]
    # Get url content
    r = requests.get(url)
    soup = bs(r.content,'html.parser')
    # Find all download links
    all_links = soup.find_all('span',{'class':"sb-file-get sb-download-link"})
    # Storing all grid and csv links in dictionary
    link_dir = {}
    for link in all_links:
        f_type = str(link.text)
        if 'grd' in f_type or 'csv' in f_type:
            link1 = str(link['data-url'])
            link_dir[f_type] = link1
                
    files = []
    for link in link_dir.keys():
        if 'csv' not in link:
            url = links[slab]
            url = f'{url[:27]}{link_dir[link]}'            
            # Opening grid file links and writing contents to local file
            r = requests.get(url)
            with open(f'{path}/{slab}_slab2_{link[10:13]}_{date}.grd','wb') as f:
                f.write(r.content)
            # Appending path to list
            f = f'{path}/{slab}_slab2_{link[10:13]}_{date}.grd'
            files.append(f)

        elif 'csv' in link:
            if 'clp' in link:
                url = links[slab]
                url = f'{url[:27]}{link_dir[link]}'            
                # Writing and appending clipping file
                r = requests.get(url)
                with open(f'{path}/{slab}_slab2_clp_{date}.csv','wb') as f:
                    f.write(r.content)
                files.append(f'{path}/{slab}_slab2_clp_{date}.csv')

            if 'sup' in link:
                url = links[slab]
                url = f'{url[:27]}{link_dir[link]}'            
                # Writing and appending supplementary file
                r = requests.get(url)
                with open(f'{path}/{slab}_slab2_sup_{date}.csv','wb') as f:
                    f.write(r.content)
                files.append(f'{path}/{slab}_slab2_sup_{date}.csv')
    for file in files:
        # Defining tilted file
        if 'sup' in file:
            tilted = file
            files.remove(file)
        # Defining clipping file
        if 'clp' in file:
            clp = file
            files.remove(file)
            
    if slab not in tilted_slabs:
        tilted = None
    
    return [files, clp, tilted]
    
def get_trench():
    '''
    Function to get trench file

    Parameters
    ----------
    None

    Returns 
    ---------
    trench : pd.DataFrame
    '''
    # URL to trench file stored on github
    url = 'https://github.com/usgs/slab2/raw/master/slab2code/plotting/forplotting/trenches_usgs_2017_depths.csv'
    # Reading link 
    x = requests.get(url=url).content 
    trench = pd.read_csv(io.StringIO(x.decode('utf8')))
    return trench


def get_ghayes(plot_dir: str) -> str:
    '''
    Function to read/write the ghayes2.cpt file stored on github

    Parameters
    ----------
    plot_dir : str
        path to slab2/slab2code/plotting

    Returns
    -------
    ghayes_cpt : str
        path to file
    '''
    ghayes_url = "https://github.com/usgs/slab2/raw/master/slab2code/plotting/forplotting/ghayes2.cpt" 
    # Write files locally if they do not exist
    ghayes_cpt =f'{plot_dir}/forplotting/ghayes2.cpt'
    if os.path.exists(ghayes_cpt) == False:
        r = requests.get(ghayes_url)
        with open(ghayes_cpt,'wb') as f:
            f.write(r.content)

    return ghayes_cpt


def get_bath(plot_dir: str) -> str:
    '''
    Function to read/write the world_lr.grd file stored on github

    Parameters
    ----------
    plot_dir : str
        path to slab2/slab2code/plotting

    Returns
    -------
    bath : str
        path to file
    '''
    # Write files locally if they do not exist
    bath_url = "https://github.com/usgs/slab2/blob/master/slab2code/plotting/forplotting/world_lr.grd?raw=true"
    bath = f'{plot_dir}/forplotting/world_lr.grd'
    if os.path.exists(bath) == False:
        r = requests.get(bath_url)
        with open(bath,'wb') as f:
            f.write(r.content)

    return bath


def get_colorblind(plot_dir: str) -> str:
    '''
    Function to read/write the colorblind cpt file stored on github

    Parameters
    ----------
    plot_dir : str
        path to slab2/slab2code/plotting

    Returns
    -------
    color : str
        path to file
    '''
    # Write files locally if they do not exist
    color_url = "https://github.com/usgs/slab2/raw/master/slab2code/plotting/forplotting/hawaii.cpt"
    color = f'{plot_dir}/forplotting/hawaii.cpt'
    if os.path.exists(color) == False:
        r = requests.get(color_url)
        with open(color,'wb') as f:
            f.write(r.content)

    return color

def get_input_df(input_date: str, slab: str) -> pd.DataFrame:
    '''
    Function to get input data from github

    Parameters
    ----------
    input_date : str
        Date of the input database
    slab : str
        3 letter slab code

    Returns
    --------
    df : pd.Dataframe
        dataframe containing the input data
    '''
    # URL of the input data
    input_url = f'https://github.com/usgs/slab2/blob/master/slab2code/Input/{input_date}/{slab}_{input_date}_input.csv?raw=true'
    x = requests.get(input_url).content 
    df = pd.read_csv(io.StringIO(x.decode('utf8')), low_memory = False)

    return df


def main():
    pass

if __name__ == '__main__':
    main()
