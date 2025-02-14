'''
File to adjust the bedlevel.xyz file to the right values
'''
#%%
## Imports
import numpy as np
import pandas as pd
import rasterio
from rtree import index as rtreeind
import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import box
from shapely.geometry import Point as sPoint
from rasterstats import zonal_stats
### Paths

AHN_raster = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\AHN\AHN4_WSS_filled.tif"
branches = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\HYDAMO\Combined_100m_28_januari.gpkg"
bedlevel_file_old = r"\\srv30\D\Ludo\P24050\Model_runs\Combined_V3.14_50m_export2_2025-01-29T12-27-47_TestSom4_KW_classMap\dflowfm\bedlevel.xyz"
inifile_path = r"\\srv30\D\Ludo\P24050\Model_runs\Combined_V3.14_50m_export2_2025-01-29T12-27-47_TestSom4_KW_classMap\dflowfm\InitialWaterLevel.ini"

def read_DHYDRO_file(file_path, keyword_list:list):
        # Read the crsdef and crsloc files
        with open(file_path, 'r') as _file:
            dhydro_file = _file.readlines()

        # Create a dictionary of all the csrdef entries
        item_dict = {}
        for keyword in keyword_list:
            for n, line in enumerate(dhydro_file):
                # Find the start of each definition
                if line.strip() == keyword:
                    j = n
                    item_id = dhydro_file[j+1].strip().split('#')[0].strip().split('=')[-1].strip()
                    item_dict[item_id] = {}
                    item_dict[item_id]['type_header'] = keyword

                    # Extract all information until a white line is found
                    while dhydro_file[j+1].strip() != '':
                        item_key = dhydro_file[j+1].strip().split('#')[0].strip().split('=')[0].strip()
                        item_value = dhydro_file[j+1].strip().split('#')[0].strip().split('=')[-1].strip()

                        if item_key not in item_dict[item_id].keys():
                            item_dict[item_id][item_key] = item_value
                        else:
                            i = 1
                            new_key = f"{item_key}_{i}"
                            while new_key in item_dict[item_id]:
                                i +=1 
                                new_key = f"{item_key}_{i}"

                            item_dict[item_id][new_key] = item_value

                        # Jump to the next line
                        j = j+1

        return item_dict

#%%
# Read the old bedlevel file
bedlevel_old = pd.read_csv(bedlevel_file_old, header=None)
bedlevel_old[['x','y','z']] = bedlevel_old[0].str.split(r'\s+',expand=True).loc[:,1:]
#%%

xnodes = bedlevel_old['x']
ynodes = bedlevel_old['y']

dx=50
dy=50

# Read in the initial values in a dictionary to use later 
inifields = read_DHYDRO_file(inifile_path, ['[Branch]'])
print('Initial water levels loaded to correct the AHN')

with rasterio.open(AHN_raster) as src:
    
    ahn = src.read(1)
    print('AHN loaded')
#%%
z = []

# Creeer een geodataframe van alle branches
branches_df = gpd.read_file(branches, layer='waterloop')

branch_index = rtreeind.Index()
for i, (branch_num, branch_line) in enumerate(branches_df.iterrows()):
    branch_index.insert(i,branch_line.geometry.bounds)
print('Created a spatial index for all branches')
count_cells_not_high_enough = 0
# Per node (xnodes,ynodes) check the intersection met die node. -> Van de nodes nog een box maken zodat er geen gaten vallen?
# Assign het peil + 5 cm voor die nodes
for ix, (x,y) in tqdm(enumerate(zip(xnodes, ynodes)), total=len(xnodes),desc='Adding height to mesh point'):
    
    # Check whether the gridcell touches one of the branches
    x = float(x)
    y = float(y)

    cell = box(x-dx, y-dy, x+dx, y+dy)
    candidate_indices = list(branch_index.intersection(cell.bounds))
    if len(candidate_indices) == 0:
        check_peil = -100
    else:
        check_peil_list = []
        for index in candidate_indices:
            branch_id, branch_line = branches_df.iloc[index][['code','geometry']]
            if branch_line.intersects(cell):
                try:
                    peilen_ini = inifields[branch_id]['values'].split(' ')
                    if len(peilen_ini) == 2:
                        peil_1, peil_2 = peilen_ini
                    elif len(peilen_ini) >2:
                        peil_1, peil_2 = peilen_ini[:2]
                    elif len(peilen_ini) == 1:
                        peil_1 = peil_2 = peilen_ini[0]
                    check_peil_list.append(np.mean((float(peil_1),float(peil_2))))
                except KeyError:
                    check_peil_list.append(-100)
            
        if len(check_peil_list) > 0:
            check_peil = np.max(check_peil_list) # Take the maximum of the peilen
        else:
            check_peil = -100
    
    p = sPoint(x, y)
    bp = p.buffer((dx + dy) // 2, cap_style=3)

    try:
        stats = zonal_stats(
            bp,
            ahn,
            affine=src.transform,
            stats="",
            add_stats={"nanmean": np.nanmean},
            nodata=np.nan,
        )
        heights = []
        for result in stats:
            heights.append(result["nanmean"])

        if heights[0] > check_peil:
            z.append(*[z for z in heights])
        else:
            z.append(check_peil+0.05) # Add 5 cm to the initial water level
            count_cells_not_high_enough += 1
    except:
        z.append(-10)

bedlevel_old['new_z'] = np.array(z).flatten()

print(f"Bed elevation lower than initial water level for {count_cells_not_high_enough} cells. Set to 0.05 m above initial waterlevel")
print("Added elevation to 2D mesh")
# %%
def format_row(row):
    return f"        {row[0]}        {row[1]}    {row[2]}"


# Apply formatting and write to file
with open(r"P:\HL-P24050\05_Analysis\05_Validatie\bedlevel.xyz", "w", newline="") as f:
    for row in bedlevel_old[['x','y','new_z']].itertuples(index=False):
        f.write(format_row(row) + "\n")

print("CSV file written successfully!")
# %%
# Read in the new bedlevel file, replace the old network mesh2d with the new bedlevel
new_bedlevel = pd.read_csv(r"P:\HL-P24050\05_Analysis\05_Validatie\bedlevel.xyz", header=None)
new_bedlevel[['x','y','z']] = new_bedlevel[0].str.strip().str.split(r'\s+',expand=True)

#%%
import sys
from pathlib import Path
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\HYDROLIB_adapted")
from hydrolib.core.io.net.models import Network
test_n = Network()
test_n = test_n.from_file(Path(r"P:\HL-P24050\05_Analysis\02_Model\Combined_V3.16_50m\dflowfm\network.nc"))
# %%
float_val_bedlevel = [float(val) for val in new_bedlevel['z'].values]
test_n._mesh2d.mesh2d_node_z = np.array(float_val_bedlevel)
test_n.to_file(Path(r"P:\HL-P24050\05_Analysis\02_Model\Combined_V3.16_50m\dflowfm\network_2.nc"))
# %%
