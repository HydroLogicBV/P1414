#########################
# File with different function definitions, used in converting HyDAMO data to a D-HyDAMO model
import numpy as np
import geopandas as gpd
import xarray as xr
import pandas as pd
import shapely 

def get_crosssection_culvert_AGV(shape: int = 1, height: float = None, 
                                width: float = None, closed:int = 1):
    '''
    Function defines a cross section based on different parameters that interact with the shape of the culvert.

    '''
    shapedict = {1:"circle",3: "rectangle", 5:"circle", 99:"rectangle"} 
    shape_str = shapedict[shape]

    # Include the diameter when the culvert is a circle
    if shape == (1 or 5): diameter = width
    else: diameter = None

    crosssection = {"shape"    : shape_str,
                    "diameter" : diameter,
                    "height"   : height,
                    "width"    : width,
                    "closed"   : closed,
                    }
    return crosssection


def getuniquecode(startstr: str, lengthlist: int) -> list:
    code = []
    for i in range(lengthlist):
        addstr = startstr+str(i)
        code.append(addstr)
    return code

def add_streefpeil_to_gemaal(streefpeilfile, inputfilename, outputfilename):

    streefpeilenxr = xr.open_dataset(streefpeilfile)
    streefpeilenxr.rio.write_crs("epsg:4326",inplace=True)

    filename = inputfilename
    weirs = gpd.read_file(filename)
    if not weirs.crs == 'epsg:4326':  weirs = weirs.to_crs('epsg:4326') 

    weirs['streefpeil'] = np.nan
    for i in range(len(weirs)):
        # Make the selection available for MultiPoints
        if type(weirs.geometry[i]) == shapely.geometry.multipoint.MultiPoint:
            lon = weirs.geometry[i][0].x
            lat = weirs.geometry[i][0].y

        # Make the selection available for single Points
        if type(weirs.geometry[i]) == shapely.geometry.point.Point:
            lon = weirs.geometry[i].x
            lat = weirs.geometry[i].y
            
        value = streefpeilenxr.sel(lon=lon, lat=lat, method="nearest").get('Band1').data
        weirs['streefpeil'].iloc[i] = value

    weirs.to_file(outputfilename)

def vormkoker_str2int(stringlistvorm: pd.Series) -> list:
    vormdict = {'Rond': 1, 'Driehoekig': 2,'Rechthoekig':3,'Eivormig':4,
                'Ellipsvormig': 5,'Paraboolvormig':6, 'Trapeziumvormig': 7,
                'Heulprofiel': 8, 'Muilprofiel': 9, 'Langwerpig': 10,
                'Scherp': 11,'Onbekend': 99}

    intlist = []

    for vorm in stringlistvorm:
        if vorm == None: vorm_to_list = None
        else: vorm_to_list = vormdict[vorm]
        intlist.append(vorm_to_list)

    return intlist