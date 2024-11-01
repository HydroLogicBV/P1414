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
folder_path_output = r"P:\HL-P24050\05_Analysis\02_Model"
#folder_path_output = r"C:\Work\Projects\P24050_ROI_voor_ROR\Testmodellen"
os.environ['GIS_folder_path'] = folder_path_GIS

#gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Buitenwater_24oktober_zonder_waterschappen.gpkg"
#gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\randvoorwaarden_test.gpkg"
gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Combined_1_november.gpkg"
#gpkg_file = folder_path_GIS + r"\GIS\HYDAMO\Buitenwater_1_november.gpkg"

gpkgs_list = [
    folder_path_GIS + r"\GIS\HYDAMO\Buitenwater_1_november.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHSK_30_oktober_no_coupuresV2.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HDSR_31_oktober.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHD_30_oktober.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\HHR_30_oktober.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\WAGV_30_oktober.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\ARKNZK.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\Rijntakken.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\RMM_HIJ_closed.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\noordzee.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\markermeer.gpkg",
    
    # folder_path_GIS + r"\GIS\HYDAMO\tunnels.gpkg",
    # folder_path_GIS + r"\GIS\HYDAMO\Ontbrekende_stuwen.gpkg",
    folder_path_GIS + r"\GIS\HYDAMO\Tunnels_Onderdoorgangen_31_oktober.gpkg",           # Should not be validated!
    # folder_path_GIS + r"\GIS\HYDAMO\Combined_21oktober_totrand.gpkg",
]
#output_folder = folder_path_output + r"\Models\Combined\V30_WBD_500"
#output_folder = folder_path_output + r"\Combined_V2.2_500m_tunnel_no_clip_V2"
output_folder = folder_path_output + r"\Combined_V2.5_100m_ALL_no_coupures"

config_dhydro = r"combined_WBD_config_100m"
config_list = [
    r"ark_nzk_config",
    r"rijntakken_config",
    r"rijnmaasmonding_open_config",
    r"noordzee_config",
    r"markermeer_config",
    # r"hhsk_config",
    # r"hdsr_config",
    # r"hhd_config",
    # r"hhr_config",
    # r"wagv_config",
    r"ontbrekende_stuwen_config",
    r"randvoorwaarden_config", 
    # r"tunnel_config",
    # r"underpass_config"                         # Should not be validated! 
]
snap_dist_list = [0, 0, 10, 10, 50, 10, 10, 100, 200, 100, 50, 0, 0 ,0]
snap_dist_list = [10, 10, 50, 100, 50, 0, 0]
defaults = r"defaults"

build_database = False
load_gpkgs = False
build_model = True

if build_database:
    dhd = DHydroData()
    for ix, config in enumerate(config_list):
        #print("\n" + config)

        dhd.hydamo_from_raw_data(
            defaults=defaults, config=config, branch_snap_dist=snap_dist_list[ix], GIS_folder=folder_path_GIS
        )
        try:
            dhd.fixed_weirs_from_raw_data(config=config, defaults=defaults)
        except AttributeError:
            pass

    dhd.clip_structures_by_branches()
    #dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    #dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults)
    #dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)
    #dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if load_gpkgs:
    dhd = DHydroData()
    for ix, gpkg in enumerate(gpkgs_list):
        print("\n" + gpkg)

        # 2. load data
        dhd.hydamo_from_gpkg(gpkg, branch_snap_dist=snap_dist_list[ix])#, GIS_folder=folder_path_GIS)

    dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=defaults)
    dhd.fixed_weirs_from_raw_data(config="noordzeekeringen_config", defaults=defaults)

    # dhd.dambreaks_from_config(config="dambreak_v0_config", defaults=defaults)
    dhd.hydamo_to_gpkg(output_gpkg=gpkg_file)

if build_model:
    # 1. initialize an instance of DHydamoData
    dhd = DHydroData()

    # 2. load data
    dhd.hydamo_from_gpkg(gpkg_file)

    # remove brug as it needs a cs
    del dhd.ddm.brug
    dhd.features.remove("brug")
    dhd.ddm.pomp["maximalecapaciteit"] = 0

    # 3. save as dhydro model
    dhd.to_dhydro(config=config_dhydro, output_folder=output_folder)
