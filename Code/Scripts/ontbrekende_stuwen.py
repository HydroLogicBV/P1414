import sys
import os

# Specify where the scripts are located
path_code = r"C:\Work\HL-P24050\P1414\Code"
#sys.path.append("D:\Work\git\GIS_tools\Code")
sys.path.append(path_code)
sys.path.append(r'C:\Work\HL-P24050\P1414\HYDROLIB_adapted\hydrolib')
sys.path.append(r'C:\Work\HL-P24050\P1414\HYDROLIB_adapted')
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\Code")
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\HYDROLIB_adapted")
sys.path.append(r"C:\Work\Projects\P24050_ROI_voor_ROR\GitHub\P1414\HYDROLIB_adapted\hydrolib")

from data_structures.dhydro_data import DHydroData

# Specify the location where the GIS folder is located and where the models must be saved:
folder_path_GIS = r"P:\HL-P24050\05_Analysis\01_GIS\03_Complete_GIS_database"
folder_path_output = r"D:\Work\Project\P1414"
os.environ['GIS_folder_path'] = folder_path_GIS

gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Ontbrekende_stuwen.gpkg"


config = r"ontbrekende_stuwen_config"
defaults = r"defaults"

build_database = True


if build_database:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. convert raw data to hydamo data
    dhd.hydamo_from_raw_data(defaults=defaults, config=config, GIS_folder=folder_path_GIS)

    # 3. save data to gpkg
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)
