from data_structures.dhydamo_data import DHydamoData

folder = r"D:\Work\Project\P1414"
gpkg_file = folder + r"\GIS\HDSR\HDSR_hydamo.gpkg"
output_gpkg_path = folder + r"\GIS\HDSR\HDSR_hydamo_v2.gpkg"
output_folder = folder + r"\Models\HDSR\V6"

# load data
dhd = DHydamoData()
dhd.from_gpkg(gpkg_file)

# save data to gpkg
dhd.ddm.to_gpkg(output_gpkg=output_gpkg_path)

# save as dhydro model
dhd.to_dhydro(output_folder=output_folder)
