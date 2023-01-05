## PATHS
p_folder = r"D:\Work\Project\P1414\GIS"
branches_path = p_folder + r"\Uitgesneden watergangen\HDSR_v4_test.shp"
bridges_path = p_folder + r"\HDSR\Legger\Bruggen\Bruggen.shp"
culvert_path = p_folder + r"\HDSR\Legger\Kokers_Lijnen\Kokers_Lijnen.shp"
pump_path = p_folder + r"\HDSR\Legger\Gemalen\Gemalen_peil.shp"
weir_path = p_folder + r"\HDSR\Legger\Stuwen\BR_Stuwen.shp"

output_gpkg = p_folder + r"\HDSR\HDSR_hydamo.gpkg"

## Branches
branch_index_mapping = dict(
    [
        ("bodembreedte", "IWS_W_BODB"),
        ("bodemhoogte benedenstrooms", "IWS_W_BODH"),
        ("bodemhoogte bovenstrooms", "IWS_W_BODH"),
        ("hoogte insteek linkerzijde", "IWS_W_INST"),
        ("hoogte insteek rechterzijde", "IWS_W_IN_1"),
        ("taludhelling linkerzijde", "IWS_W_TALU"),
        ("taludhelling rechterzijde", "IWS_W_TA_1"),
        ("typeruwheid", None),
        ("ruwheidhoog", None),
        ("ruwheidlaag", None),
        ("water_width_index", "IWS_W_WATB"),
    ]
)

## Bridges
bridge_index_mapping = dict(
    [
        ("code", "CODE"),
        ("geometry", "geometry"),
        ("globalid", "globalid"),
        ("intreeverlies", None),
        ("typeruwheid", None),
        ("ruwheid", None),
        ("uittreeverlies", None),
        ("lengte", "WS_DOORVAA"),
    ]
)

## Culverts
culvert_index_mapping = dict(
    [
        ("breedteopening", "BREEDTEOPE"),
        ("code", "CODE"),
        ("geometry", "geometry"),
        ("globalid", "globalid"),
        ("hoogtebinnenonderkantbene", "HOOGTEBOKB"),
        ("hoogtebinnenonderkantbov", "HOOGTEBO_1"),
        ("hoogteopening", "HOOGTEOPEN"),
        ("intreeverlies", None),
        ("lengte", "LENGTE"),
        ("typeruwheid", None),
        ("ruwheid", None),
        ("uittreeverlies", None),
        ("vormkoker", "VORMKOKER"),
    ]
)

## Pumps
pumpingstation_index_mapping = dict(
    [
        ("code", "GEMAALID"),
        ("geometry", "geometry"),
        ("globalid", "globalid"),
    ]
)
pumps = dict()
management = dict()
