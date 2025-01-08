'''
De rijntakken hebben niet het juiste initiële waterstandsniveau. Dit script maakt een losse InitialWaterLevels.ini bestand aan voor
specifiek de Rijntakken en de RMM o.b.v. waterdieptes en profielinformatie.
'''
#%%
import pandas as pd
import numpy as np
import os
import geopandas as gpd
#### First define some helper functions

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

                    if item_id not in item_dict.keys():
                        item_id_in_dict = item_id
                            
                    else:
                        i = 1
                        item_id_in_dict = f"{item_id}_{i}"
                        while item_id_in_dict in item_dict.keys():
                            i +=1 
                            item_id_in_dict = f"{item_id}_{i}"

                    item_dict[item_id_in_dict] = {}
                    item_dict[item_id_in_dict]['type_header'] = keyword

                    # Extract all information until a white line is found
                    while dhydro_file[j+1].strip() != '':
                        item_key = dhydro_file[j+1].strip().split('#')[0].strip().split('=')[0].strip()
                        item_value = dhydro_file[j+1].strip().split('#')[0].strip().split('=')[-1].strip()

                        if item_key not in item_dict[item_id_in_dict].keys():
                            item_dict[item_id_in_dict][item_key] = item_value
                        else:
                            i = 1
                            new_key = f"{item_key}_{i}"
                            while new_key in item_dict[item_id_in_dict]:
                                i +=1 
                                new_key = f"{item_key}_{i}"

                            item_dict[item_id_in_dict][new_key] = item_value

                        # Jump to the next line
                        j = j+1

        return item_dict

def write_DHYDRO_file(new_dict: dict, 
                      folder: str, 
                      filename: str, 
                      filetype:str ='1dField', 
                      fileversion:str ='2.00',
                      globalquantity:str = 'waterlevel',
                      globalunit:str ='m',
                      globalvalue:str ='0.0'):
    new_file_name = os.path.join(folder, filename)
    
    # Write the new initial file
    with open(new_file_name, mode='w') as file:
        # Write the general part
        # file.write('[General]\n')
        # file.write(f'\tfileVersion\t= {fileversion}\n')
        # file.write(f'\tfileType\t={filetype}\n')
        # file.write('\n')
        # file.write('[Global]\n')
        # file.write(f'\tquantity\t= {globalquantity}\n')
        # file.write(f'\tunit\t= {globalunit}\n')
        # file.write(f'\tvalue\t= {globalvalue}\n')
        # file.write('\n')

        
        # Write a new entry per crosssection
        for ID, data in new_dict.items():
            file.write(f"{data['type_header']}\n")

            for keys, values in data.items():
                if keys == 'type_header': continue

                tabs = "\t" * int(np.ceil(5 - len(keys)/4))
                file.write('\t'+ keys + tabs + "= " + str(values) + '\n')

            file.write('\n')

# Set the paths to the right information
RTK_ini_waterdepth_path =r"P:\P1414\05_Analysis\03_RWS_modellen\sobek-rijn-j22_v1a1\sobek-rijn-j22_v1a\dflow1d\InitialWaterDepth.ini"
RMM_ini_waterlevel_path = r"P:\P1414\05_Analysis\03_RWS_modellen\Rijn-MaasMonding-j15_5-v2\modellen\sobek\sobek-rmm-vzm-j15_5-v4\sobek-rmm-vzm-j15_5-v4.dsproj_data\rmm_output\dflow1d\InitialWaterLevel.ini"

# Set the Rijntakken and RMM shapes path
RTK_shape_path = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\RTK_BranchesV3.shp"
RMM_shape_path = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\RMM\RMM_Branches_V3.shp"

# Set the Rijntakken profile path
RTK_prof_path = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken\ZW_cross_v3.csv"

# Load the necessary information
RTK_ini_waterdepth = read_DHYDRO_file(RTK_ini_waterdepth_path, ['[Definition]'])
RMM_ini_waterlevel = read_DHYDRO_file(RMM_ini_waterlevel_path, ['[Definition]'])
RTK_branches = gpd.read_file(RTK_shape_path)
RMM_branches = gpd.read_file(RMM_shape_path)
RTK_profs = pd.read_csv(RTK_prof_path)
# %%

#####################################################################
### MATCH THE NEW WATERWAYS WITH THE ORIGINAL INITIAL WATERDEPTHS ###
#####################################################################

# Start of with the Rijntakken
RTK_new_ini_dict = {}

# Set the original values in a dataframe
RTK_df = pd.DataFrame(RTK_ini_waterdepth).T
RTK_df['chainage'] = np.round(RTK_df['chainage'].astype(float),2)

# Lek and Nederrijn have been split up
Lek_takken = RTK_branches[RTK_branches['Name'].str.startswith('Lek')]
Nederrijn_takken = RTK_branches[RTK_branches['Name'].str.startswith('Nederrijn')]

Lek_lengte = Lek_takken.geometry.length.to_list()
Lek_lengte_cum = np.cumsum(Lek_lengte)
Lek_ids = Lek_takken['Name'].to_list()

# # Nederrijn
Nederrijn_lengte = Nederrijn_takken.geometry.length.to_list()
Nederrijn_lengte_cum = np.cumsum(Nederrijn_lengte)
Nederrijn_ids = Nederrijn_takken['Name'].to_list()

for n, row in RTK_branches.iterrows():
    RTK_new_ini_dict[row['Name']] = {
        'type_header': '[Branch]',
        'branchId': 'rijn_wl_'+row['Name'],
    }
    chainages = []
    values = []
    profiles_branch = RTK_profs[(RTK_profs['branch'] == row['Name']) & (RTK_profs['Data_type'] =='meta')]
    if row['Name'] in Lek_ids:
        Lek_idx = Lek_ids.index(row['Name'])
        # Determine the extent of the Lek reach
        if Lek_idx == 0: 
            min_chain = 0
        else:
            min_chain = Lek_lengte_cum[Lek_idx-1]
        
        for k, prof in profiles_branch.iterrows():
            min_level = np.nanmin(RTK_profs[RTK_profs['id'] == prof['id']]['level'])
            # Select the matching waterdepths with the profile
            RTK_waterdepths = RTK_df[(RTK_df['branchId']=='Lek')]
            waterdepth = RTK_waterdepths.loc[(RTK_waterdepths['chainage'] - (prof['chainage'] + min_chain)).abs().idxmin()]['value']
            
            # Add the chainages and the values 
            chainages.append(str(prof['chainage']))
            values.append(str(np.round(min_level + float(waterdepth), 2)))

    elif row['Name'] in Nederrijn_ids:
        Nederrijn_idx = Nederrijn_ids.index(row['Name'])
        # Determine the extent of the Nederrijn reach
        if Nederrijn_idx == 0: min_chain = 0
        else: min_chain = Nederrijn_lengte_cum[Nederrijn_idx-1]

        for k, prof in profiles_branch.iterrows():
            min_level = np.nanmin(RTK_profs[RTK_profs['id'] == prof['id']]['level'])
            # Select the matching waterdepths with the profile
            RTK_waterdepths = RTK_df[(RTK_df['branchId']=='Nederrijn')]
            waterdepth = RTK_waterdepths.loc[(RTK_waterdepths['chainage'] - (prof['chainage'] + min_chain)).abs().idxmin()]['value']
            
            # Add the chainages and the values 
            chainages.append(str(prof['chainage']))
            values.append(str(np.round(min_level + float(waterdepth), 2)))
    else:
        for k, prof in profiles_branch.iterrows():
            min_level = np.nanmin(RTK_profs[RTK_profs['id'] == prof['id']]['level'])

            # Select the matching waterdepths with the profile
            RTK_waterdepths = RTK_df[(RTK_df['branchId']==row['Name'])]
            waterdepth = RTK_waterdepths.loc[(RTK_waterdepths['chainage'] - prof['chainage']).abs().idxmin()]['value']
            
            # Add the chainages and the values 
            chainages.append(str(prof['chainage']))
            values.append(str(np.round(min_level + float(waterdepth), 2)))

    RTK_new_ini_dict[row['Name']]['numLocations'] = str(len(chainages))
    RTK_new_ini_dict[row['Name']]['chainage'] = ' '.join(chainages)
    RTK_new_ini_dict[row['Name']]['values'] = ' '.join(values)
   

# Add information for Lek_3
RTK_new_ini_dict['Lek_3'] = RTK_new_ini_dict['Lek_2']
RTK_new_ini_dict['Lek_3']['chainage'] = '110.3'

write_DHYDRO_file(RTK_new_ini_dict, 
                folder = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\Rijntakken",
                filename='InitialWaterLevels_Rijn.ini')

# %%

###############################################################
###### CONVERT WATERLEVELS FOR RMM ############################
###############################################################

# Start of with the Rijntakken
RMM_new_ini_dict = {}

# Set the original values in a dataframe
RMM_df = pd.DataFrame(RMM_ini_waterlevel).T
RMM_df['chainage'] = np.round(RMM_df['chainage'].astype(float),2)

# HollandseIJssel and NieuweMaas04 have been split up
HIJ_takken = RMM_branches[RMM_branches['Name'].str.startswith('HollandseIJssel1')]
NM_takken = RMM_branches[RMM_branches['Name'].str.startswith('NieuweMaas04')]

HIJ_lengte = HIJ_takken.geometry.length.to_list()
HIJ_lengte_cum = np.cumsum(HIJ_lengte)
HIJ_ids = HIJ_takken['Name'].to_list()

# NM
NM_lengte = NM_takken.geometry.length.to_list()
NM_lengte_cum = np.cumsum(NM_lengte)
NM_ids = NM_takken['Name'].to_list()

for n, row in RMM_branches.iterrows():
    RMM_new_ini_dict[row['Name']] = {
        'type_header': '[Branch]',
        'branchId': 'rijn_wl_'+row['Name'],
    }
    chainages = []
    values = []
    #profiles_branch = RMM_profs[(RMM_profs['branch'] == row['Name']) & (RTK_profs['Data_type'] =='meta')]
    if row['Name'] in HIJ_ids:
        HIJ_idx = HIJ_ids.index(row['Name'])
        # Determine the extent of the Lek reach
        if HIJ_idx == 0: 
            min_chain = 0
        else:
            min_chain = HIJ_lengte_cum[HIJ_idx-1]
        max_chain = HIJ_lengte_cum[HIJ_idx]
        
        RMM_waterlevels = RMM_df[(RMM_df['branchId']=='HollandseIJssel1') &
                                 (RMM_df['chainage'] >= min_chain) &
                                (RMM_df['chainage'] < max_chain)]
        for n_row, RMM_row in RMM_waterlevels.iterrows():
            chainages.append(str(RMM_row['chainage'] - max_chain))
            values.append(str(RMM_row['value']))
        
    elif row['Name'] in NM_ids:
        NM_idx = NM_ids.index(row['Name'])
        # Determine the extent of the Lek reach
        if NM_idx == 0: 
            min_chain = 0
        else:
            min_chain = NM_lengte_cum[NM_idx-1]
        max_chain = NM_lengte_cum[NM_idx]
        
        RMM_waterlevels = RMM_df[(RMM_df['branchId']=='NieuweMaas04') &
                                 (RMM_df['chainage'] >= min_chain) &
                                (RMM_df['chainage'] < max_chain)]
        for n_row, RMM_row in RMM_waterlevels.iterrows():
            chainages.append(str(RMM_row['chainage'] - max_chain))
            values.append(str(RMM_row['value']))
    else:
        RMM_waterlevels = RMM_df[(RMM_df['branchId']==row['Name'])]
        for n_row, RMM_row in RMM_waterlevels.iterrows():
            chainages.append(str(RMM_row['chainage']))
            values.append(str(RMM_row['value']))

    RMM_new_ini_dict[row['Name']]['numLocations'] = str(len(chainages))
    RMM_new_ini_dict[row['Name']]['chainage'] = ' '.join(chainages)
    RMM_new_ini_dict[row['Name']]['values'] = ' '.join(values)
   

# Add information for Lek_3
#RMM_new_ini_dict['Lek_3'] = RTK_new_ini_dict['Lek_2']
#RMM_new_ini_dict['Lek_3']['chainage'] = '110.3'

write_DHYDRO_file(RMM_new_ini_dict, 
                folder = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database\GIS\RMM",
                filename='InitialWaterLevels_RMM.ini')

# %%
