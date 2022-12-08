# Dit script laadt de streefpeilenkaart in als NetCDF en extract de gewenste streefpeilen bij de 
# verschillende gemalen om de streefpeilen op te halen 

import data_functions as daf
import importlib
importlib.reload(daf)

streefpeilfile = "D:\work\P1414_ROI\peilgebieden_zp_25m_full2_25m.nc"

# Add streefpeil AGV
daf.add_streefpeil_to_gemaal(streefpeilfile = streefpeilfile,
                             inputfilename = "D:\work\P1414_ROI\GIS\WAGV\pomp_gemaal_v13\pomp_gemaal_v13_clipped.shp",
                             outputfilename = "D:\work\P1414_ROI\GIS\WAGV\pomp_gemaal_v13\pomp_gemaal_v13_clipped_peil.shp" )


# Add streefpeil HDSR
daf.add_streefpeil_to_gemaal(streefpeilfile = streefpeilfile,
                             inputfilename = "D:\work\P1414_ROI\GIS\HDSR\Legger\Gemalen\Gemalen.shp",
                             outputfilename = "D:\work\P1414_ROI\GIS\HDSR\Legger\Gemalen\Gemalen_peil.shp" )

# Add streefpeil HHSK
daf.add_streefpeil_to_gemaal(streefpeilfile = streefpeilfile,
                             inputfilename = "D:\work\P1414_ROI\GIS\HHSK\Legger\Gemaal.shp",
                             outputfilename = "D:\work\P1414_ROI\GIS\HHSK\Legger\Gemaal_peil.shp" )

# Add streefpeil HHDL
daf.add_streefpeil_to_gemaal(streefpeilfile = streefpeilfile,
                             inputfilename = "D:\work\P1414_ROI\GIS\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Gemaal.shp",
                             outputfilename = "D:\work\P1414_ROI\GIS\HHDelfland\Legger_Delfland_shp\Ondersteunende kunstwerken\Gemaal_peil.shp" ) 

# Add streefpeil HHRL
daf.add_streefpeil_to_gemaal(streefpeilfile = streefpeilfile,
                             inputfilename = "D:\work\P1414_ROI\GIS\HHRijnland\Legger\Gemaal\gemaal.shp",
                             outputfilename = "D:\work\P1414_ROI\GIS\HHRijnland\Legger\Gemaal\gemaal_peil.shp" )                                                                            