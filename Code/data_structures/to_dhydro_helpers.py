import importlib
import os
from copy import copy
from pathlib import Path
from typing import List, Union
from datetime import datetime, timedelta
import pandas as pd
from tqdm import tqdm
import geopandas as gpd
import numpy as np
import rasterio
from rtree import index as rtreeind
from geo_tools.roughness_to_mesh import create_roughness_xyz
from hydrolib.core.basemodel import DiskOnlyFileModel
from hydrolib.core.io.crosssection.models import CrossDefModel, CrossLocModel
from hydrolib.core.io.dimr.models import DIMR, FMComponent
from hydrolib.core.io.ext.models import ExtModel
from hydrolib.core.io.friction.models import FrictBranch, FrictGlobal, FrictionModel
from hydrolib.core.io.ini.models import INIBasedModel, INIGeneral, INIModel
from hydrolib.core.io.inifield.models import IniFieldModel, InitialField, ParameterField
from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.onedfield.models import OneDFieldBranch, OneDFieldGlobal, OneDFieldModel
from hydrolib.core.io.polyfile.models import Description, Metadata
from hydrolib.core.io.polyfile.models import Point as hPoint
from hydrolib.core.io.polyfile.models import PolyFile, PolyObject
from hydrolib.core.io.structure.models import Dambreak, StructureModel
from hydrolib.dhydamo.converters.df2hydrolibmodel import Df2HydrolibModel
from hydrolib.dhydamo.core.hydamo import HyDAMO
from hydrolib.dhydamo.geometry import mesh
from hydrolib.dhydamo.io.dimrwriter import DIMRWriter
from rasterstats import zonal_stats
from scipy.spatial import KDTree
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Point as sPoint
from shapely.geometry import Polygon, MultiPolygon, box

from shapely.ops import unary_union, linemerge

from data_structures.roi_data_model import ROIDataModel as DataModel

FM_FOLDER = "dflowfm"
OUTPUTDIR = "output"


def to_dhydro(
    self,
    config: str,
    extent: Union[gpd.GeoDataFrame, Polygon] = None,
    initial_1D_waterdepth: float = 1,
    load_mesh2d_path = None,
    output_folder=None,
):
    """
    Function to convert datamodel to D-HYDRO model

    Args:
        config (str): configuration file
        extent (Union[gpd.GeodataFrame, Polygon]): extent of 2D grid
        initial_1D_waterdepth (float): initial 1D waterdepth (defaults to 1)
        ouput_folder (str): output folder path

    Returns:
        None
    """

    def add_1D_initial_waterdepth(fm: FMModel, fm_path: Path, initial_1D_waterdepth: float):

        globalwd = OneDFieldGlobal(quantity="waterdepth", unit="m", value=initial_1D_waterdepth)
        onedfile_wd = OneDFieldModel(global_=globalwd)
        onedfile_wd.filepath = fm_path / "Initialwaterdepth.ini"
        onedfile_wd.save()

        initialwd = InitialField(
            quantity="waterdepth",
            datafile=DiskOnlyFileModel(filepath=Path("Initialwaterdepth.ini")),
            datafiletype="1dField",
            locationtype="1d",
        )

        # inifieldmodel = IniFieldModel(initial=[initialwd])
        # fm.geometry.inifieldfile = inifieldmodel
        if hasattr(fm.geometry, "inifieldfile") and isinstance(
            fm.geometry.inifieldfile, IniFieldModel
        ):
            fm.geometry.inifieldfile.initial.append(initialwd)
        else:
            fm.geometry.inifieldfile = IniFieldModel(initial=[initialwd])

        return fm

    def add_1D_initial_waterlevel(
        fm: FMModel, fm_path: Path, ddm: DataModel = None, initial_1D_waterlevl=0
    ):
        branches_gdf = ddm.waterloop
        branch_field_list = []
        for name, branch in branches_gdf.iterrows():
            if np.isnan(branch["peil"]):
                continue

            branch_field_list.append(
                OneDFieldBranch(
                    branchid=branch["code"],
                    numlocations=2,  # zero locations does not work
                    chainage=[
                        np.around(branch["geometry"].length * 0.1, 3),
                        np.around(branch["geometry"].length * 0.9, 3),
                    ],
                    values=[np.around(branch["peil"], 3), np.around(branch["peil"], 3)],
                )
            )

        globalwl = OneDFieldGlobal(quantity="waterlevel", unit="m", value=initial_1D_waterlevl)
        onedfile_wl = OneDFieldModel(branch=branch_field_list, global_=globalwl)
        onedfile_wl.filepath = fm_path / "Initialwaterlevels.ini"
        onedfile_wl.save()

        initialwl = InitialField(
            quantity="waterlevel",
            datafile=DiskOnlyFileModel(filepath=Path("Initialwaterlevels.ini")),
            datafiletype="1dField",
            locationtype="1d",
        )
        # inifieldmodel = IniFieldModel(initial=[initialwl])
        # fm.geometry.inifieldfile = inifieldmodel

        if hasattr(fm.geometry, "inifieldfile") and isinstance(
            fm.geometry.inifieldfile, IniFieldModel
        ):
            fm.geometry.inifieldfile.initial.append(initialwl)
        else:
            fm.geometry.inifieldfile = IniFieldModel(initial=[initialwl])

        return fm

    def add_branches(ddm: DataModel, features: List[str], hydamo: HyDAMO) -> HyDAMO:
        hydamo.branches.set_data(
            ddm.waterloop, index_col="code", check_geotype=False
        )  # using globalid leads to errors in D-HYDRO. Probably too long name

        # Check for circular features in branches
        hydamo.branches_popped = copy(hydamo.branches)  # .copy()
        for _, branch in hydamo.branches.iterrows():

            start = np.array(branch.geometry.coords)[0]
            end = np.array(branch.geometry.coords)[-1]
            if np.array_equal(start, end):
                code = branch.code
                hydamo.branches_popped = hydamo.branches_popped.drop(code)
                print(f'Dropped branch with code: {code} because it is circular')
            if branch.geometry.length < 1e-3:
                hydamo.branches_popped = hydamo.branches_popped.drop(code)

        hydamo.branches.delete_all()
        hydamo.branches.set_data(
            hydamo.branches_popped.set_geometry("geometry"), index_col="code", check_geotype=False
        )

        if (
            ("profielgroep" in features)
            and ("profielpunt" in features)
            and ("profiellijn" in features)
            and ("ruwheidsprofiel" in features)
        ):
            # hydamo.profile.read_gpkg_layer(
            #     gpkg_path,
            #     layer_name="profielpunt",
            #     groupby_column="profiellijnid",
            #     order_column="codevolgnummer",
            #     index_col="code",
            # )

            # Group profile_points by line_id and create LineStrings
            profiel_punt_gdf = ddm.profielpunt
            

            # Filter out profiles with only one point in advance
            profiel_punt_gdf_filtered = profiel_punt_gdf.groupby('profiellijnid').filter(lambda x: len(x) > 1)
            grouped = profiel_punt_gdf_filtered.groupby("profiellijnid")
            profiel_punt_lines = []

            for line, profiel_punten in tqdm(grouped, total=len(grouped), desc='Adding profiles with grouped'):
                profiel_punten = profiel_punten.sort_values(by="codevolgnummer", ascending=True)
                linestring = LineString(profiel_punten.geometry.values)
                profiel_punt_lijn = {
                    "code": line,
                    "geometry": linestring,
                    "globalid": profiel_punten.iloc[0]["globalid"],
                    "profiellijnid": line,
                    "codevolgnummer": 1,
                }
                profiel_punt_lines.append(profiel_punt_lijn)

            # unique_lines = profiel_punt_gdf["profiellijnid"].unique()
            # profiel_punt_lines = []

            # for ix, line in tqdm(enumerate(unique_lines), total=len(unique_lines),desc='Adding profiles to branches'):
            #     profiel_punten = copy(profiel_punt_gdf[profiel_punt_gdf["profiellijnid"] == line])
            #     profiel_punten.sort_values(
            #         by="codevolgnummer", axis=0, ascending=True, inplace=True
            #     )

            #     # skip profiles with only one point. MIght be result of clipping
            #     if profiel_punten.shape[0] > 1:
            #         linestring = LineString(profiel_punten.geometry.values)

            #         profiel_punt_lijn = dict(
            #             [
            #                 ("code", line),
            #                 ("geometry", linestring),
            #                 ("globalid", profiel_punten.iloc[0, :]["globalid"]),
            #                 ("profiellijnid", line),
            #                 ("codevolgnummer", 1),
            #             ]
            #         )
            #         profiel_punt_lines.append(profiel_punt_lijn)

            profiel_punt_lines_gdf = gpd.GeoDataFrame(
                data=profiel_punt_lines, geometry="geometry", crs=profiel_punt_gdf.crs
            )

            hydamo.profile.set_data(profiel_punt_lines_gdf)
            print('Added profile data')
            # Snap profiles to branch
            # hydamo.profile_roughness.read_gpkg_layer(gpkg_path, layer_name="ruwheidsprofiel")
            hydamo.profile_roughness.set_data(ddm.ruwheidsprofiel)
            hydamo.profile.snap_to_branch(hydamo.branches, snap_method='intersecting')# Was: snap_method="intersecting")
            print('Snapped profiles to branches')
            hydamo.profile.dropna(axis=0, inplace=True, subset=["branch_offset"])
            # hydamo.profile_line.read_gpkg_layer(gpkg_path, layer_name="profiellijn")
            # hydamo.profile_group.read_gpkg_layer(gpkg_path, layer_name="profielgroep")
            hydamo.profile_line.set_data(ddm.profiellijn)
            hydamo.profile_group.set_data(ddm.profielgroep)

            hydamo.profile.drop("code", axis=1, inplace=True)
            hydamo.profile["code"] = hydamo.profile["profiellijnid"]

        if ("hydroobject_normgp" in features) and ("normgeparamprofielwaarde" in features):
            hydamo.param_profile.set_data(ddm.hydroobject_normgp, index_col="globalid")
            hydamo.param_profile_values.set_data(
                ddm.normgeparamprofielwaarde, index_col="normgeparamprofielid"
            )
        return hydamo

    def add_bridges(
        ddm: DataModel, hydamo: HyDAMO, default_height=0, max_snap_dist: float = 5
    ) -> HyDAMO:

        hydamo.bridges.set_data(ddm.brug, index_col="code")
        hydamo.bridges.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
        hydamo.bridges.dropna(axis=0, inplace=True, subset=["branch_offset"])
        for i, bridge in hydamo.bridges.iterrows():
            if bridge.lengte is None:
                continue
            hydamo.structures.add_bridge(
                id=bridge.code,
                name=bridge.code,
                length=bridge.lengte,
                branchid=bridge.branch_id,
                chainage=bridge.branch_offset,
                frictiontype=bridge.typeruwheid,
                csdefid=bridge.code,  # TODO Check influence of csdefid
                shift=default_height,  # TODO Validate the use of offset
                friction=bridge.ruwheid,
                inletlosscoeff=bridge.intreeverlies,
                outletlosscoeff=bridge.uittreeverlies,
            )

        return hydamo

    def add_culverts(ddm: DataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
        hydamo.culverts.set_data(ddm.duiker)

        for name, culvert in ddm.duiker.iterrows():
            if isinstance(culvert["geometry"], LineString):
                hydamo.culverts.loc[name, "geometry"] = hydamo.culverts.loc[
                    name, "geometry"
                ].interpolate(0.5, normalized=True)

        hydamo.culverts.snap_to_branch(
            hydamo.branches, snap_method="overal", maxdist=max_snap_dist
        )
        hydamo.culverts.dropna(axis=0, inplace=True, subset=["branch_offset"])
        for i, culvert in hydamo.culverts.iterrows():
            hydamo.structures.add_culvert(
                id=culvert.code,
                name=culvert.code,
                branchid=culvert.branch_id,
                chainage=culvert.branch_offset,
                leftlevel=culvert.hoogtebinnenonderkantbov,
                rightlevel=culvert.hoogtebinnenonderkantbene,
                length=culvert.lengte,
                inletlosscoeff=culvert.intreeverlies,
                outletlosscoeff=culvert.uittreeverlies,
                crosssection=eval(
                    culvert.doorstroomopening
                ),  # dict saved as string, so eval turns it into dict
                bedfrictiontype=culvert.typeruwheid,
                bedfriction=culvert.ruwheid,
            )

        return hydamo

    def add_dambreaks(ddm: DataModel, fm=FMModel, max_dist=100, z_default=0) -> FMModel:
        db_gdf = ddm.doorbraak
        fw_gdf = ddm.keringen

        fw_point_list = []
        for ix, fw in fw_gdf.iterrows():
            coords = fw.geometry.coords[:]

            for x, y, z in coords:
                fw_point_list.append([x, y, z])

        fw_kdtree = KDTree(data=fw_point_list)

        # mesh1d = fm.geometry.netfile.network._mesh1d
        # ids = mesh1d.network1d_node_id
        # xs = mesh1d.network1d_node_x
        # ys = mesh1d.network1d_node_y

        # id_list = []
        # wl_node_list = []
        # for x, y, id in zip(xs, ys, ids):
        #     wl_node_list.append([x, y])
        #     id_list.append(id)

        # wl_kdtree = KDTree(data=wl_node_list)

        for name, db in db_gdf.iterrows():
            geom = db.geometry
            if not isinstance(geom, LineString):
                raise ValueError("I need a LineString...")

            coords = np.array(geom.coords[:])
            centroid = geom.interpolate(0.5, normalized=True)
            if len(coords[0]) == 2:
                mid_point = [centroid.x, centroid.y, z_default]
            elif len(coords[0]) == 3:
                mid_point = [centroid.x, centroid.y, centroid.z]
            else:
                raise ValueError("can't cope with {} coordinates".format(len(coords[0])))

            _, ix_point_in_kering = fw_kdtree.query(mid_point, distance_upper_bound=max_dist)
            point_in_kering = fw_point_list[ix_point_in_kering]

            # _, ix_point_in_wl = wl_kdtree.query(coords[0, :], distance_upper_bound=max_dist)
            # wl_node_id = id_list[ix_point_in_wl]

            dambreak = Dambreak(
                algorithm=db.algorithm,
                breachwidthini=db.breachwidthini,
                crestlevelini=point_in_kering[2] - db.crestlevelini,
                crestlevelmin=db.crestlevelmin,
                f1=db.f1,
                f2=db.f2,
                id=db.code,
                name=db.globalid,
                numcoordinates=coords.shape[0],
                startlocationx=point_in_kering[0],
                startlocationy=point_in_kering[1],
                t0=db.t0,
                timetobreachtomaximumdepth=db.timetobreachtomaximumdepth,
                ucrit=db.ucrit,
                waterleveldownstreamlocationx=coords[-1, 0],
                waterleveldownstreamlocationy=coords[-1, 1],
                waterlevelupstreamlocationx=coords[0, 0],
                waterlevelupstreamlocationy=coords[0, 1],
                # waterlevelupstreamlocationx=xs[ix_point_in_wl],
                # waterlevelupstreamlocationy=ys[ix_point_in_wl],
                # waterlevelupstreamnodeid=wl_node_id,
                xcoordinates=coords[:, 0].tolist(),
                ycoordinates=coords[:, 1].tolist(),
            )
            if hasattr(fm.geometry, "structurefile") and isinstance(
                fm.geometry.structurefile[0], StructureModel
            ):
                fm.geometry.structurefile[0].structure.append(dambreak)
            else:
                fm.geometry.structurefile = StructureModel(structure=[dambreak])

        return fm

    def add_fixed_weirs(ddm: DataModel, fm=FMModel, data: List = [0, 0, 5, 4, 4, 0]) -> FMModel:
        def split_line_by_point(line, max_length: float = 10000, max_seg_length: float = 10):
            length = line.length
            if length < max_length:
                return line

            n_seg = np.ceil(length / max_length).astype(int)
            l_list = []
            for n in range(0, n_seg):
                sp = n / n_seg
                ep = (n + 1) / n_seg
                points_norm = np.linspace(sp, ep, np.ceil(max_length / max_seg_length).astype(int))
                points = []
                for p in points_norm:
                    points.append(line.interpolate(p, normalized=True))
                l_list.append(LineString(points))

            mls = MultiLineString(l_list)
            return mls

        fw_list = []
        if hasattr(ddm, "keringen") and ddm.keringen is not None:
            fw_list.append(ddm.keringen)
        if hasattr(ddm, "overiglijnelement") and ddm.overiglijnelement is not None:
            fw_list.append(ddm.overiglijnelement)

        line_list = []
        for fw_gdf in fw_list:
            fw_gdf = (
                fw_gdf.assign(
                    geometry=fw_gdf.apply(
                        lambda x: split_line_by_point(
                            x.geometry,
                        ),
                        axis=1,
                    )
                )
                .explode(index_parts=False)
                .reset_index(drop=True)
            )
            #fw_gdf["geometry"] = fw_gdf["geometry"].simplify(tolerance=0.1)

            for name, row in fw_gdf.iterrows():
                coord_list = row.geometry.coords[:]
                point_list = []
                for (x, y, z) in coord_list:
                    point_list.append(hPoint(x=x, y=y, z=z, data=data))

                metadata = Metadata(name=str(row["code"]), n_rows=len(coord_list), n_columns=9)
                polyline = PolyObject(
                    description=Description(content=str(row["code"])),
                    metadata=metadata,
                    points=point_list,
                )
                line_list.append(polyline)

        fm.geometry.fixedweirfile = [PolyFile(has_z_values=True, objects=line_list)]
        return fm

    def add_pumps(ddm: DataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
        hydamo.pumpstations.set_data(ddm.gemaal)
        # hydamo.pumps.set_data(ddm.pomp)
        hydamo.pumpstations.snap_to_branch(
            hydamo.branches, snap_method="overal", maxdist=max_snap_dist
        )
        hydamo.pumpstations.dropna(axis=0, inplace=True, subset=["branch_offset"])

        hydamo.pumps.set_data(ddm.pomp[ddm.pomp["gemaalid"].isin(hydamo.pumpstations["globalid"])])
        hydamo.management.set_data(ddm.sturing)

        hydamo.structures.convert.pumps(
            hydamo.pumpstations, pumps=hydamo.pumps, management=hydamo.management
        )

        return hydamo

    def add_river_profiles(ddm: DataModel, hydamo: HyDAMO) -> HyDAMO:
        branches = ddm.waterloop
        for name, prof in ddm.rivier_profielen_data.iterrows():
            if not branches["code"].str.contains(prof["branchid"]).any():
                print(f"Branch for profile {prof['name']} not found")
                continue

            branch = branches.loc[branches["code"] == prof["branchid"], :]

            if prof["chainage"] > np.ceil(branch.geometry.length.values[:]):
                print(f"Profile {prof['name']} too far on branch {prof['branchid']}. Length: {branch.geometry.length.values[:]}, Chainage: {prof['chainage']}")
                continue

            profile = ddm.rivier_profielen[ddm.rivier_profielen["name"] == prof["name"]]
            profile = profile.sort_values(by="levels", ascending=True)

            kwargs = dict(
                [
                    ("flowWidths", profile["flowWidths"].values.tolist()),
                    ("fp1width", prof["fp1width"]),
                    ("fp2width", prof["fp2width"]),
                    ("frictiontypes", None),
                    ("frictionvalues", None),
                    ("mainwidth", prof["mainwidth"]),
                    ("name", prof["name"]),
                    ("numLevels", prof["numLevels"]),
                    ("leveebaselevel", prof["leveebaselevel"]),
                    ("leveecrestLevel", prof["leveecrestLevel"]),
                    ("leveeflowarea", prof["leveeflowarea"]),
                    ("leveetotalarea", prof["leveetotalarea"]),
                    ("levels", profile["levels"].values.tolist()),
                    ("thalweg", prof["thalweg"]),
                    ("totalWidths", profile["totalWidths"].values.tolist()),
                ]
            )
            if prof["frictionids"] is None or not isinstance(prof['frictionids'], str):
                kwargs["frictionids"] = None
            else:
                kwargs["frictionids"] = prof["frictionids"].split(",")

            hydamo.crosssections.add_zw_river_definition(**kwargs)

            # hydamo.crosssections.add_zw_river_definition(
            #     name=prof["name"],
            #     numLevels=prof["numLevels"],
            #     levels=profile["levels"].values.tolist(),
            #     flowWidths=profile["flowWidths"].values.tolist(),
            #     totalWidths=profile["totalWidths"].values.tolist(),
            #     thalweg=prof["thalweg"],
            #     leveecrestLevel=prof["leveecrestLevel"],
            #     leveebaselevel=prof["leveebaselevel"],
            #     leveeflowarea=prof["leveeflowarea"],
            #     leveetotalarea=prof["leveetotalarea"],
            #     mainwidth=prof["mainwidth"],
            #     fp1width=prof["fp1width"],
            #     fp2width=prof["fp2width"],
            #     frictionids=prof["frictionids"].split(","),
            #     frictiontypes=None,
            #     frictionvalues=None,
            # )
            hydamo.crosssections.add_crosssection_location(
                branchid=prof["branchid"], chainage=prof["chainage"], definition=prof["name"]
            )

        return hydamo

    def add_river_roughness(ddm: DataModel, fm=FMModel, default=45) -> FMModel:
        branches = ddm.waterloop
        gdf = ddm.rivier_ruwheid
        sections = gdf.section.unique()

        friction_file_list = []
        for section in sections:
            frict_list = []

            fric_def_global = FrictGlobal(
                frictionid=section, frictiontype="Chezy", frictionvalue=default
            )

            section_slice = gdf.loc[gdf["section"] == section, :]
            for ix, (name, row) in enumerate(section_slice.iterrows()):
                if not branches["code"].str.contains(row["branchid"]).any():
                    print(f"Branch for roughness {name} not found")
                    continue

                branch = branches.loc[branches["code"] == row["branchid"], :]
                if row["chainage"] is not None and not (isinstance(row["chainage"], float) and np.isnan(row["chainage"])):
                    chainage = np.amax([float(c) for c in row["chainage"].split(",")])
                    if chainage > np.ceil(branch.geometry.length.values[:]):
                        print(f"Profile {name} too far on branch {row['branchid']}. Length: {branch.geometry.length.values[:]}, Chainage: {chainage}")
                        continue

                kwargs = copy(
                    row[
                        [
                            "branchid",
                            "frictiontype",
                            "functiontype",
                            "numlevels",
                            "levels",
                            "numlocations",
                            "chainage",
                            "frictionvalues",
                        ]
                    ]
                )

                if kwargs["chainage"] is not None and not (isinstance(kwargs["chainage"], float) and np.isnan(kwargs["chainage"])):
                    chainage = [float(i) for i in kwargs["chainage"].split(",")]
                    numlocations = int(kwargs["numlocations"])
                else:
                    chainage = [1]
                    numlocations = 1
                kwargs = kwargs.drop("chainage")
                kwargs = kwargs.drop("numlocations")

                frictionvalues = [float(i) for i in kwargs["frictionvalues"].split(",")]
                if kwargs["levels"] is not None and not (isinstance(kwargs["levels"], float) and np.isnan(kwargs["levels"])):
                    levels = [float(i) for i in kwargs["levels"].split(",")]
                    numlevels = int(kwargs["numlevels"])
                    frictionvalues = (
                        np.reshape(
                            frictionvalues,
                            (numlocations, numlevels),
                        )
                        .T.astype(float)
                        .tolist()
                    )
                else:
                    levels = None
                    numlevels = None
                    frictionvalues = (
                        np.reshape(
                            frictionvalues,
                            (numlocations, 1),
                        )
                        .T.astype(float)
                        .tolist()
                    )

                kwargs = kwargs.drop("frictionvalues")
                kwargs = kwargs.drop("levels")
                kwargs = kwargs.drop("numlevels")

                frict_branch = FrictBranch(
                    **kwargs,
                    chainage=chainage,
                    frictionvalues=frictionvalues,
                    levels=levels,
                    numlevels=numlevels,
                    numlocations=numlocations,
                )

                frict_list.append(frict_branch)
            friction_model = FrictionModel(branch=frict_list, global_=fric_def_global)
            friction_model.filepath = r"roughness_{}.ini".format(section)
            fm.geometry.frictfile.append(friction_model)
        return fm

    def add_river_roughness_default(fm: FMModel, default=45, section: str = "Main") -> FMModel:
        fric_def_global = FrictGlobal(
            frictionid=section, frictiontype="Chezy", frictionvalue=default
        )
        friction_model = FrictionModel(global_=fric_def_global)
        friction_model.filepath = r"roughness_{}.ini".format(section)
        fm.geometry.frictfile.append(friction_model)
        return fm

    def add_tunnels(
        hydamo: HyDAMO,
        fm: FMModel,
        max_dist_to_struct=3,
        extent=None,
    ) -> FMModel:
        if "tunnel" in hydamo.branches.columns:
            branches = hydamo.branches.loc[hydamo.branches.tunnel, :]
            
        else:
            return fm
        # if extent is None:
        #     pass
        # elif isinstance(extent, gpd.GeoDataFrame):
        #     mask = branches['geometry'].apply(lambda line: line.within(extent['geometry'][0]))
        #     branches = branches[mask] #TODO: Profiles are not deleted!
        #     #branches = branches.overlay(gpd.GeoDataFrame(extent["geometry"]), how="within")
        # elif isinstance(extent, Polygon) or isinstance(extent, gpd.GeoSeries):
        #     mask = branches['geometry'].apply(lambda line: line.within(extent))
        #     branches = branches[mask]
        #     #branches = branches.overlay(gpd.GeoDataFrame(geometry=extent), how="within")
        network = fm.geometry.netfile.network

        
        mesh.mesh1d_add_branches_from_gdf(
            network,
            branches=branches,
            branch_name_col="code",
            node_distance=1e6,
            max_dist_to_struc=max_dist_to_struct,
            structures=None,
        )

        # if two_d:
        #     branch_list = branches["code"].values[:]
        #     # node_mask = network._mesh1d.get_node_mask(branch_list)
        #     # print(np.sum(node_mask))

        #     # # If not provided, create a box from the maximum bounds
        #     # xmin = network._mesh2d.mesh2d_node_x.min()
        #     # xmax = network._mesh2d.mesh2d_node_x.max()
        #     # ymin = network._mesh2d.mesh2d_node_y.min()
        #     # ymax = network._mesh2d.mesh2d_node_y.max()

        #     # within = box(xmin, ymin, xmax, ymax)
        #     # geometrylist = GeometryList.from_geometry(within)
        #     # network._link1d2d._link_from_1d_to_2d(node_mask, polygon=geometrylist)

        #     for ix, branchid in enumerate(branch_list):
        #         node_mask = network._mesh1d.get_node_mask(branchid)

        #         point1 = sPoint(network._mesh1d.branches[branchid].node_xy[0, :])
        #         point2 = sPoint(network._mesh1d.branches[branchid].node_xy[-1, :])
        #         _within = point1.buffer(250).union(point2.buffer(250))
        #         if ix == 0:
        #             within = _within
        #         else:
        #             within = within.union(_within)
        #         # mp = MultiPolygon(polygons=[point1.buffer(250), point2.buffer(250)])
        #         # geometrylist = GeometryList.from_geometry(mp)
        #         # network._link1d2d._link_from_1d_to_2d(node_mask, polygon=geometrylist)
        #     mesh.links1d2d_add_links_2d_to_1d_embedded(
        #         network=network, branchids=branch_list, within=within
        #     )
        #     # mesh.links1d2d_add_links_2d_to_1d_embedded(network=network, branchids=branch_list)
        #     print(network._link1d2d.link1d2d.shape)

        return fm

    def add_weirs(ddm: DataModel, hydamo: HyDAMO, max_snap_dist: float = 5) -> HyDAMO:
        hydamo.weirs.set_data(ddm.stuw)
        hydamo.opening.set_data(ddm.kunstwerkopening)
        hydamo.management_device.set_data(ddm.regelmiddel)

        hydamo.weirs.snap_to_branch(hydamo.branches, snap_method="overal", maxdist=max_snap_dist)
        hydamo.weirs.dropna(axis=0, inplace=True, subset=["branch_offset"])

        hydamo.structures.convert.weirs(
            weirs=hydamo.weirs, opening=hydamo.opening, management_device=hydamo.management_device
        )

        return hydamo

    def add_boundary_condition(fm: FMModel, hydamo: HyDAMO):
        
        hydamo.external_forcings.convert.boundaries(hydamo.boundary_conditions, mesh1d=fm.geometry.netfile.network)

        return(hydamo)

    def build_1D_model(
        features: List,
        fm: FMModel,
        hydamo: HyDAMO,
        node_distance=20,
        max_dist_to_struct=3,
    ) -> HyDAMO:

        kwargs = {}
        if "brug" in features:
            # bbrug = True
            kwargs["bridges"] = True
        else:
            kwargs["bridges"] = False

        if "duiker" in features:
            kwargs["culverts"] = True
        else:
            kwargs["culverts"] = False

        if "gemaal" in features:
            kwargs["pumps"] = True
        else:
            kwargs["pumps"] = False

        if "stuw" in features:
            kwargs["rweirs"] = True
            kwargs["uweirs"] = True
        else:
            kwargs["rweirs"] = False
            kwargs["uweirs"] = False

        # Collect all structures in a Structures dataframe
        structures = hydamo.structures.as_dataframe(**kwargs)

        # if min_dist_between_struct > 0:
        #     for branchid, branch in hydamo.branches.iterrows():
        #         struct_bool = structures["branchid"] == branchid
        #         slice = structures.loc[struct_bool, :]
        #         if slice.shape[0] > 1:
        #             chainages = slice["chainage"]

        #             branch_length = hydamo.branches.at[branchid, "geometry"].length
        #             if branch_length < (min_dist_between_struct * 2):
        #                 struct_chainage = branch_length / 2
        #                 structures.loc[struct_bool, "chainage"] = struct_chainage
        #             else:
        #                 for structid, struct in slice.iterrows():
        #                     struct_chainage = struct["chainage"]

        #                     if struct_chainage > min_dist_between_struct:
        #                         new_chainage = np.around(
        #                             (struct_chainage // min_dist_between_struct)
        #                             * min_dist_between_struct
        #                         )
        #                     else:
        #                         new_chainage = np.around(min_dist_between_struct)

        #                     structures.at[structid, "chainage"] = new_chainage

        if "tunnel" in hydamo.branches.columns:
            branches = hydamo.branches.loc[~hydamo.branches.tunnel, :]
            # structures = structures.loc[structures["branchid"].isin(branches["code"]), :]
        else:
            branches = hydamo.branches

        mesh.mesh1d_add_branches_from_gdf(
            fm.geometry.netfile.network,
            branches=branches,
            branch_name_col="code",
            node_distance=node_distance,
            max_dist_to_struc=max_dist_to_struct,
            structures=structures,
        )

        # Set the crosssections for the branches
        hydamo.crosssections.convert.profiles(
            crosssections=hydamo.profile,
            crosssection_roughness=hydamo.profile_roughness,
            # profile_groups=hydamo.profile_group, # skip as no bridges or weirs with profile
            profile_lines=hydamo.profile_line,
            param_profile=hydamo.param_profile,
            param_profile_values=hydamo.param_profile_values,
            branches=hydamo.branches,
            roughness_variant="High",
        )

        # Set a default cross section
        default = hydamo.crosssections.add_rectangle_definition(
            height=5.0,
            width=5.0,
            closed="no",
            roughnesstype="StricklerKs",
            roughnessvalue=30,
            name="default",
        )
        hydamo.crosssections.set_default_definition(definition=default, shift=-2.5)

        if ('hydrologischerandvoorwaarde') in features:
                hydamo = add_boundary_condition(fm=fm, hydamo=hydamo)

        models = Df2HydrolibModel(hydamo)

        extmodel = ExtModel()
        extmodel.boundary = models.boundaries_ext
        extmodel.lateral = models.laterals_ext

        # Add external forcings
        fm.external_forcing.extforcefilenew = extmodel

        # Assign to FM model
        fm.geometry.structurefile = [StructureModel(structure=models.structures)]
        fm.geometry.crosslocfile = CrossLocModel(crosssection=models.crosslocs)
        fm.geometry.crossdeffile = CrossDefModel(definition=models.crossdefs)

        fm.geometry.frictfile = []
        for i, fric_def in enumerate(models.friction_defs):
            fric_model = FrictionModel(global_=fric_def)
            fric_model.filepath = f"roughness_{i}.ini"
            fm.geometry.frictfile.append(fric_model)

        return fm, hydamo

    def build_2D_model(
        coupling_type: str,
        dx: float,
        dy: float,
        elevation_raster_path: str,
        extent: gpd.GeoDataFrame,
        fm: FMModel,
        fm_path,
        initial_peil_raster_path: str,
        roughness_2d_raster_path: str,
        two_d_buffer: float,
        clip_extent: gpd.GeoDataFrame = None,
        clip_buffer: int = 0,
        load_mesh2d_path: str = None, # TODO: Implementeren dat er een grid ingeladen kan worden! 
        lateral_branches: list = None,
        branch_1d2d_embed_list: list = None,
        one_d: bool = False,
        bound_ids: list = ["mark", "noor", "rijn"],
        max_lat_dist: float = 2e3,
    ) -> FMModel:

        network = fm.geometry.netfile.network

        if load_mesh2d_path is None:
            print('Start constructing the mesh2d grid, because no exisiting grid has been specified. This might take a while!')
            if two_d_buffer is None:
                polygon = extent.geometry.values[0]
            elif isinstance(two_d_buffer, (int, float)):
                polygon = extent.geometry.values[0].buffer(two_d_buffer)
            else:
                raise ValueError("not supported buffer")

            mesh.mesh2d_add_rectilinear(network=network, polygon=polygon, dx=dx, dy=dy)
            print("created 2D mesh")

            # Clip main branches from the 2D mesh and give them lateral links
            # Clip polygons from the 2D mesh
            clipped = 0
            if clip_extent is not None:
                # Explode potential MultiPolygons:
                clip_extent = clip_extent.explode(index_parts=True).reset_index(drop=True)
                
                if len(clip_extent) == 1:
                    clip_polygons = clip_extent.loc[0, 'geometry']
                    mesh.mesh2d_clip(network=network, polygon=clip_polygons)
                    print('Clipped 1 polygon from the mesh\n')
                else:
                    n_finish = 0
                    for i, poly in tqdm(enumerate(clip_extent.iterrows()),total=len(clip_extent),desc="clip progress", unit=" rows"):
                        clip_polygons = poly[1].geometry.buffer(clip_buffer)

                        # Ensure the result is a single polygon although the result should be exploded already
                        if clip_polygons.type == 'MultiPolygon':
                            merged_polygon = unary_union(clip_polygons)
                        else:
                            merged_polygon = clip_polygons

                        if merged_polygon.type == 'MultiPolygon':
                            # If still a MultiPolygon, extract the largest polygon
                            largest_polygon = max(merged_polygon, key=lambda p: p.area)
                        else:
                            largest_polygon = merged_polygon

                        clip_polygon = Polygon(largest_polygon.exterior)    # Make a polygon with the exterior = filling up potential holes
                        clipped += 1
                        mesh.mesh2d_clip(network=network, polygon=clip_polygon)
                        n_finish +=1
                    print(f'Clipped {clipped} polygons one by one from the mesh\n')

                # Add lateral 1D-2D links because the regular links are removed when clipped
                if len(lateral_branches) > 0:
                    mesh.links1d2d_add_links_2d_to_1d_lateral(
                            network=network,
                            branchids=lateral_branches,
                            dist_factor=None,   # within=extent.buffer(1e3)
                            max_length=2000,      # Set the maximum length to the grid resolution
                        )
                    
                    # Remove the links that go to end points as these deliver problems when running the simulation
                    #mesh.links1d2d_remove_1d_endpoints(network=network)
                    print('Added lateral links for the clipped branches and removed the ones to end points.')
            # ---> End add Pepijn

            if elevation_raster_path is not None:
                # mesh.mesh2d_altitude_from_raster(
                #     network=network,
                #     rasterpath=elevation_raster_path,
                #     where="node",  # Face does not work
                #     stat="nanmean",
                #     fill_option="fill_value",
                #     fill_value=-10,
                #     window_size=int(dx + dy),
                # )

                xnodes = network._mesh2d.mesh2d_node_x
                ynodes = network._mesh2d.mesh2d_node_y

                # Read in the initial values in a dictionary to use later 
                inifile_path = fm_path / fm.geometry.inifieldfile.initial[0].datafile.filepath
                inifields = read_DHYDRO_file(inifile_path, ['[Branch]'])
                print('Initial water levels loaded to correct the AHN')
                with rasterio.open(elevation_raster_path) as src:
                    #     bounds = src.bounds
                    #     _box = box(*bounds)
                    #     points = []
                    #     ixs = []
                    #     for ix, (x, y) in enumerate(zip(xnodes, ynodes)):
                    #         p = sPoint(x, y)
                    #         if p.within(_box):
                    #             bp = p.buffer((dx + dy) // 2, cap_style=3)
                    #             points.append(bp)
                    #             ixs.append(ix)

                    #     stats = zonal_stats(
                    #         bp,
                    #         src.read(1),
                    #         affine=src.transform,
                    #         stats="",
                    #         add_stats={"nanmean": np.nanmean},
                    #         nodata=np.nan,
                    #     )
                    #     heights = []
                    #     for result in stats:
                    #         heights.append(result["nanmean"])

                    #     # convert to numpy array and remove values more than two standard devations from the nanmean
                    #     z = np.array([z for z in heights]).flatten()
                    #     z_new = np.full((network._mesh2d.mesh2d_node_y.shape[0]), -10)
                    #     for ix, (x, y) in enumerate(zip(xnodes, ynodes)):
                    #         if ix in ixs:
                    #             z_new[ix] = z
                    ahn = src.read(1)
                    print('AHN loaded')
                    z = []

                    # Creeer een geodataframe van alle branches
                    branches_list = [(branch_id, LineString(coords.geometry)) for branch_id, coords in network._mesh1d.branches.items()]
                    
                    branch_index = rtreeind.Index()
                    for i, (branch_id, branch_line) in enumerate(branches_list):
                        branch_index.insert(i,branch_line.bounds)
                    print('Created a spatial index for all branches')
                    count_cells_not_high_enough = 0
                    # Per node (xnodes,ynodes) check the intersection met die node. -> Van de nodes nog een box maken zodat er geen gaten vallen?
                    # Assign het peil + 5 cm voor die nodes
                    for ix, (x, y) in tqdm(enumerate(zip(xnodes, ynodes)), total=len(xnodes),desc='Adding height to mesh point'):
                        
                        # Check whether the gridcell touches one of the branches
                    
                        cell = box(x-dx, y-dy, x+dx, y+dy)
                        candidate_indices = list(branch_index.intersection(cell.bounds))
                        if len(candidate_indices) == 0:
                            check_peil = -100
                        else:
                            check_peil_list = []
                            for index in candidate_indices:
                                branch_id, branch_line = branches_list[index]
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

                    network._mesh2d.mesh2d_node_z = np.array(z).flatten()

                print(f"Bed elevation lower than initial water level for {count_cells_not_high_enough} cells. Set to 0.05 m above initial waterlevel")
                print("Added elevation to 2D mesh")
        else:
            # Use a pre-specified grid
            network._mesh2d.read_file(Path(load_mesh2d_path))
            print('Pre-existing mesh has been loaded!')

        if initial_peil_raster_path is not None:
            initial_wl = InitialField(
                averagingtype="nearestNb",
                datafile=DiskOnlyFileModel(filepath=Path(initial_peil_raster_path)),
                datafiletype="GeoTIFF",
                interpolationmethod="averaging",
                locationtype="2d",
                quantity="waterlevel",
            )
            if hasattr(fm.geometry, "inifieldfile") and isinstance(
                fm.geometry.inifieldfile, IniFieldModel
            ):
                fm.geometry.inifieldfile.initial.append(initial_wl)
            else:
                fm.geometry.inifieldfile = IniFieldModel(initial=[initial_wl])
            print("added initial peilen to 2D mesh")
        # else:
        # initial_wd = InitialField(
        #     quantity="waterdepth",
        #     datafile=DiskOnlyFileModel(filepath=Path(initial_peil_raster_path)),
        #     datafiletype="polygon",
        #     interpolationmethod="constant",
        #     locationtype="2d",
        #     value=0,
        # )
        # if hasattr(fm.geometry, "inifieldfile") and isinstance(
        #     fm.geometry.inifieldfile, IniFieldModel
        # ):
        #     fm.geometry.inifieldfile.initial.append(initial_wd)
        # else:
        #     fm.geometry.inifieldfile = IniFieldModel(initial=[initial_wd])
        # print("added initial peilen to 2D mesh")
        # fm.geometry.waterlevini = -10

        if roughness_2d_raster_path is not None:
            path_f = os.path.split(roughness_2d_raster_path)[0]
            xyz_path = path_f + "\\roughness_sample.xyz"
            create_roughness_xyz(
                xnodes=network._mesh2d.mesh2d_face_x,
                ynodes=network._mesh2d.mesh2d_face_y,
                dx=dx,
                dy=dy,
                roughness_tif_path=roughness_2d_raster_path,
                roughness_mesh_name=path_f + "\\roughness_mesh_v2.tif",
                roughness_xyz_name=xyz_path,
                convert_to_manning=True,
            )

            initial_rn = ParameterField(
                quantity="frictioncoefficient",
                datafile=DiskOnlyFileModel(filepath=Path(xyz_path)),
                datafiletype="sample",
                interpolationmethod="averaging",
                locationtype="2d",
                ifrctyp=1,
            )

            if hasattr(fm.geometry, "inifieldfile") and isinstance(
                fm.geometry.inifieldfile, IniFieldModel
            ):
                fm.geometry.inifieldfile.parameter.append(initial_rn)
            else:
                fm.geometry.inifieldfile = IniFieldModel(parameter=[initial_rn])

            print("added roughness to 2D mesh")

        # branchids = network._mesh1d.network1d_branch_id
        # branchids = [x for x in branchids if not x.startswith(r"rijn_")]
        branchids = None
        polygon = None
        if one_d:
            if coupling_type == "1Dto2D":
                mesh.links1d2d_add_links_1d_to_2d(
                    branchids=branch_1d2d_embed_list, network=network, within=polygon
                )
                print('1Dto2D embedded links added')
            elif coupling_type == "2Dto1D":
                mesh.links1d2d_add_links_2d_to_1d_embedded(
                    branchids=branch_1d2d_embed_list, network=network, within=polygon
                )
                print('2Dto1D embedded links added')
            else:
                print("No 1D to 2D links have been set")

            branchids = network._mesh1d.network1d_branch_id
            branch_list = []
            for id in branchids:
                for prefix in bound_ids:
                    if prefix in id:
                        branch_list.append(id)
            if len(branch_list) > 0:
                mesh.links1d2d_add_links_2d_to_1d_lateral(
                    network=network,
                    branchids=branch_list,
                    dist_factor=None,  # within=extent.buffer(1e3)
                    max_length=max_lat_dist,
                )
                print('Lateral 1d2d links added for the boundaries (Rijntakken, RMM, Noordzee, Markermeer')

            # Add lateral 1D-2D links because the regular links are removed when clipped
            if (len(lateral_branches) > 0) and (clip_extent is not None):
                mesh.links1d2d_add_links_2d_to_1d_lateral(
                        network=network,
                        branchids=lateral_branches,
                        dist_factor=None,   # within=extent.buffer(1e3)
                        max_length=2000,      # Set the maximum length to the grid resolution
                    )
                
                # Remove the links that go to end points as these deliver problems when running the simulation
                #mesh.links1d2d_remove_1d_endpoints(network=network)
                print('Added lateral links for the clipped branches and removed the ones to end points.')

            network._link1d2d.link1d2d, unique_idx = np.unique(network._link1d2d.link1d2d, axis=0, return_index=True)
            network._link1d2d.link1d2d_id = network._link1d2d.link1d2d_id[unique_idx]
            network._link1d2d.link1d2d_long_name = network._link1d2d.link1d2d_long_name[unique_idx]
            network._link1d2d.link1d2d_contact_type = network._link1d2d.link1d2d_contact_type[unique_idx]
        return fm

    def merge_ext_files(ext_force_model, extforcefile: ExtModel) -> None:
        '''
        Merge existing ExtModel files if already exisiting
        '''
        #assert isinstance(fmmodel, FMModel)

        if (not hasattr(ext_force_model, "extforcefilenew")) or (
            ext_force_model.extforcefilenew is None
        ):
            ext_force_model.extforcefilenew = extforcefile
        else:
            if hasattr(extforcefile, "lateral"):
                if (not hasattr(ext_force_model.extforcefilenew, "lateral")) or (
                    ext_force_model.extforcefilenew.lateral is None
                ):
                    ext_force_model.extforcefilenew.lateral = extforcefile.lateral
                else:
                    for lateral in extforcefile.lateral:
                        ext_force_model.extforcefilenew.lateral.append(lateral)
            if hasattr(extforcefile, "boundary"):
                if (not hasattr(ext_force_model.extforcefilenew, "boundary")) or \
                    (ext_force_model.extforcefilenew.boundary is None) or \
                    (isinstance(ext_force_model.extforcefilenew.boundary, list) and not \
                    ext_force_model.extforcefilenew.boundary): # Check for an empty list

                    ext_force_model.extforcefilenew.boundary = extforcefile.boundary
                else:
                    new_forcingmodel = ext_force_model.extforcefilenew.boundary[0].forcingfile.forcing

                    for boundary in extforcefile.boundary:
                        new_forcingmodel.append(boundary.forcingfile.forcing[0])
                        ext_force_model.extforcefilenew.boundary.append(boundary)

                    # Add the updated forcing model to all boundaries    
                    for n, bound in enumerate(ext_force_model.extforcefilenew.boundary):
                        ext_force_model.extforcefilenew.boundary[n].forcingfile.forcing = new_forcingmodel
                        
        return ext_force_model

    def set_hydrolib_core_options(fmmodel: FMModel, options):
        """ """

        def set_options(obj_loc, options):
            """
            sets options recursively
            """
            # loop over all attributes of object
            for attribute, value in options.__dict__.items():
                # Skip python special attributes
                # if attribute.startswith("__"):
                if "__" in attribute:
                    continue
                # if an option is found, try to set it
                # options should be int, float, list, or str type
                if isinstance(value, (int, float, list, str, INIBasedModel, INIGeneral, INIModel)):
                    if isinstance(value, (ExtModel)):
                        obj_loc = merge_ext_files(obj_loc, value)
                    else:
                        try:
                            setattr(obj_loc, attribute, value)
                        except Exception as e:
                            print("failed to set {}".format(attribute))
                            print(e)
                # if not, a subclass was encounterd and we need to traverse deeper
                else:
                    set_options(obj_loc=getattr(obj_loc, attribute), options=value)

            return obj_loc

        return set_options(obj_loc=fmmodel, options=options)

    def simple_fm_model(
        dtmax: int = 60,
        start_time: int = 20160601,
        stop_time: int = 2 * 86400,
        statsinterval: int = 3600,
    ) -> FMModel:
        fm = FMModel()
        fm.time.refdate = start_time
        fm.time.tstop = stop_time
        fm.time.dtmax = dtmax
        fm.output.statsinterval = [statsinterval]
        return fm


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


    def write_DHYDRO_file(new_dict: dict, 
                        folder: str, 
                        filename: str):
        new_file_name = os.path.join(folder, filename)
        
        # Write the new initial file
        with open(new_file_name, mode='w') as file:           
            # Write a new entry per crosssection
            for ID, data in new_dict.items():
                file.write(f"{data['type_header']}\n")

                for keys, values in data.items():
                    if keys == 'type_header': continue

                    tabs = "\t" * int(np.ceil(5 - len(keys)/4))
                    file.write('\t'+ keys + tabs + "= " + str(values) + '\n')

                file.write('\n')

    # load configuration file
    model_config = getattr(importlib.import_module("dataset_configs." + config), "Models")

    # build FM model if set in config file
    if hasattr(model_config, "FM"):
        # initialize a simple FM model
        self.fm = simple_fm_model(
            start_time=model_config.FM.start_time, stop_time=model_config.FM.stop_time
        )
        self.fm.output.outputdir = OUTPUTDIR

        self.hydamo = HyDAMO()

        if model_config.FM.two_d_bool:  # compute extent of 1D network
            # Add to gpgk?

            if hasattr(model_config.FM.two_d, "extent_path") and (
                model_config.FM.two_d.extent_path is not None
            ):
                extent = gpd.read_file(model_config.FM.two_d.extent_path).dissolve(by=None)
            else:
                extent = self.ddm.waterloop.dissolve(by=None).convex_hull.buffer(
                    model_config.FM.two_d.two_d_buffer
                )
            # if "doorbraak" in self.features:
            #     buf_doorbraak = copy(self.ddm.doorbraak)
            #     buf_doorbraak["geometry"] = buf_doorbraak["geometry"].buffer(
            #         (model_config.FM.two_d.dx + model_config.FM.two_d.dy) / 2
            #     )
            #     extent = gpd.overlay(extent, buf_doorbraak, how="union").dissolve(by=None)

        # build 1D model
        if model_config.FM.one_d_bool:
            # check if there are branches in the data
            if "waterloop" not in self.features:
                raise ValueError("Missing branches")

            # add branches, profiles and norm profiles
            print("\nworking on 1D branches\n")
            self.hydamo = add_branches(ddm=self.ddm, features=self.features, hydamo=self.hydamo)

            if ("rivier_profielen_data" in self.features) and (
                "rivier_profielen" in self.features
            ):
                print("\nworking on river profiles\n")
                self.hydamo = add_river_profiles(ddm=self.ddm, hydamo=self.hydamo)

            # Loop over structures and add when present in the data
            struct_functions = {
                "brug": add_bridges,
                "duiker": add_culverts,
                "gemaal": add_pumps,
                "stuw": add_weirs,
            }

            for structure, function in struct_functions.items():
                if structure in self.features:
                    print("\nworking on {}\n".format(structure))
                    self.hydamo = function(
                        ddm=self.ddm,
                        hydamo=self.hydamo,
                        max_snap_dist=model_config.FM.one_d.max_snap_dist,
                    )

            if ("tunnel" in self.ddm.waterloop.columns) and (model_config.FM.two_d_bool):
                self.fm = add_tunnels(
                    hydamo=self.hydamo,
                    fm=self.fm,
                    max_dist_to_struct=model_config.FM.one_d.max_dist_to_struct,
                    extent=extent,
                )

            if ('hydrologischerandvoorwaarde') in self.features:
                self.hydamo.boundary_conditions.set_data(self.ddm.hydrologischerandvoorwaarde)
                self.hydamo.boundary_conditions.snap_to_branch(self.hydamo.branches, snap_method='overal', maxdist=10)


            # add 1D model
            print("\nCompleting 1D model\n")
            self.fm, self.hydamo = build_1D_model(
                fm=self.fm,
                features=self.features,
                hydamo=self.hydamo,
                node_distance=model_config.FM.one_d.node_distance,
                max_dist_to_struct=model_config.FM.one_d.max_dist_to_struct,
            )
            if ("rivier_profielen_data" in self.features) and (
                "rivier_profielen" in self.features
            ):
                if "rivier_ruwheid" in self.features:
                    self.fm = add_river_roughness(ddm=self.ddm, fm=self.fm)
                else:
                    self.fm = add_river_roughness_default(fm=self.fm)

            if ("peil" in self.ddm.waterloop.columns) and (output_folder is not None):
                fm_path = Path(output_folder) / FM_FOLDER
                self.fm = add_1D_initial_waterlevel(
                    ddm=self.ddm,
                    fm=self.fm,
                    fm_path=fm_path,
                )
                
                inifile_path = fm_path / self.fm.geometry.inifieldfile.initial[0].datafile.filepath
                inifields = read_DHYDRO_file(inifile_path, ['[General]', '[Global]','[Branch]'])

                # Corrigeer de initiële peilen op de Rijntakken en RMM
                if hasattr(model_config.FM.one_d, "rijntakken_ini_path"):
                    rijntakken_ini = read_DHYDRO_file(model_config.FM.one_d.rijntakken_ini_path, ['[Branch]'])
                    for key in rijntakken_ini.keys():
                        if key in inifields.keys():
                            inifields[key] = rijntakken_ini[key]
                        else:
                            print(f'{key} not present in original branches')
                    print('Updated initial waterlevels for the Rijntakken')
                
                if hasattr(model_config.FM.one_d, 'rmm_ini_path'):
                    rmm_ini = read_DHYDRO_file(model_config.FM.one_d.rmm_ini_path, ['[Branch]'])
                    for key in rmm_ini.keys():
                        if key in inifields.keys():
                            inifields[key] = rmm_ini[key]
                        else:
                            print(f'{key} not present in original branches')
                    print('Updated initial waterlevels for the RMM')

                if hasattr(model_config.FM.one_d, 'correctie_overig_ini_path'):
                    corr_ini = read_DHYDRO_file(model_config.FM.one_d.correctie_overig_ini_path, ['[Branch]'])
                    for key in corr_ini.keys():
                        if key in inifields.keys():
                            inifields[key] = corr_ini[key]
                        else:
                            print(f'{key} not present in original branches')
                    print('Updated initial waterlevels for the corrections')
                
                write_DHYDRO_file(inifields, fm_path, self.fm.geometry.inifieldfile.initial[0].datafile.filepath)

            else:
                self.fm = add_1D_initial_waterdepth(
                    fm=self.fm,
                    fm_path=fm_path,
                    initial_1D_waterdepth=initial_1D_waterdepth,
                )

        # build 2D model
        if model_config.FM.two_d_bool:  # compute extent of 1D network
            # Add to gpgk?

            if hasattr(model_config.FM.two_d, "extent_path") and (
                model_config.FM.two_d.extent_path is not None
            ):
                extent = gpd.read_file(model_config.FM.two_d.extent_path).dissolve(by=None)
            else:
                extent = self.ddm.waterloop.dissolve(by=None).convex_hull.buffer(
                    model_config.FM.two_d.two_d_buffer
                )
            # if "doorbraak" in self.features:
            #     buf_doorbraak = copy(self.ddm.doorbraak)
            #     buf_doorbraak["geometry"] = buf_doorbraak["geometry"].buffer(
            #         (model_config.FM.two_d.dx + model_config.FM.two_d.dy) / 2
            #     )
            #     extent = gpd.overlay(extent, buf_doorbraak, how="union").dissolve(by=None)

            option_list = [
                "coupling_type",
                "elevation_raster_path",
                "initial_peil_raster_path",
                "roughness_2d_raster_path",
                "two_d_buffer",
            ]
            kwargs = {}
            for option in option_list:
                if hasattr(model_config.FM.two_d, option) and (
                    getattr(model_config.FM.two_d, option) is not None
                ):
                    kwargs[option] = getattr(model_config.FM.two_d, option)
                else:
                    kwargs[option] = None


            # Get the geometries that must be clipped
            CRS = 'EPSG:28992'   

            # Stuk uit code Pepijn ter test:
            if hasattr(model_config.FM.two_d, "clip_extent_path") and (model_config.FM.two_d.clip_extent_path is not None):
                # Concat possible multiple clip_extent_paths
                if isinstance(model_config.FM.two_d.clip_extent_path, dict):
                    for key, value in model_config.FM.two_d.clip_extent_path.items():
                        if "base" in key:
                            gdf = gpd.read_file(value, crs=CRS)
                        elif "concat" in key:
                            gdf = gpd.GeoDataFrame(
                                data=pd.concat([gdf, gpd.read_file(value, crs=CRS)]),
                                geometry="geometry",
                                crs=gdf.crs,
                            )
                    clip_extent_gdf = gdf
                else:
                    clip_extent_gdf = gpd.read_file(model_config.FM.two_d.clip_extent_path, crs=CRS)
                
                #clip_extent = clip_extent_gdf[clip_extent_gdf['geometry'].apply(lambda geom: isinstance(geom, Polygon))]    # Only select the ones with a polygon, because that can only be clipped
                clip_extent = clip_extent_gdf.copy()

                # Apply selection criteria:
                #   1. all shapes that have 'meer, plas' in their 'typewater' column must be clipped (must remain in clip_extent)
                #   2. all shapes that meet the clip_selection criteria for the 'breedtekla' column must be clipped
                if hasattr(model_config.FM.two_d, "clip_selection") and (model_config.FM.two_d.clip_selection is not None):
                    condition1 = (clip_extent['typewater'] == 'meer, plas')
                    condition2 = (clip_extent['typewater'] == 'waterloop') & (clip_extent['breedtekla'].isin(model_config.FM.two_d.clip_selection))
                    clip_extent = clip_extent[condition1 | condition2]

                # Ensure that only unique column names exist
                if 'index_left' in clip_extent.columns:
                    clip_extent = clip_extent.rename(columns={'index_left': 'clip_index_left'})
                if 'index_right' in clip_extent.columns:
                    clip_extent = clip_extent.rename(columns={'index_right': 'clip_index_right'})
            else:
                clip_extent = None
                clip_buffer = 0
                lateral_branches = []
            
            # Get the clip buffer, because it is not defined everywhere, to prevent errors set the buffer to 100m if not defined.
            if hasattr(model_config.FM.two_d, "clip_buffer"):
                clip_buffer = model_config.FM.two_d.clip_buffer
            else:
                clip_buffer = 100
                print('** Assumed default clip buffer of 100m because no buffer was given.')
            
            # Get a list of the branchid's that are clipped and now need 1D-2D lateral links
            if hasattr(model_config.FM.two_d, "clip_extent_path") and (model_config.FM.two_d.clip_extent_path is not None):
                all_branches_gdf = gpd.GeoDataFrame(self.hydamo.branches, geometry='geometry', crs=CRS)
                clip_extent_buffer = clip_extent.copy()
                               
                clip_extent_buffer = clip_extent_buffer.set_geometry('geometry')
                all_branches_gdf = all_branches_gdf.set_geometry('geometry')
                clip_extent_buffer = clip_extent_buffer.to_crs(CRS)
                all_branches_gdf = all_branches_gdf.to_crs(CRS)

                clipped_branches_gdf = gpd.sjoin(all_branches_gdf, clip_extent_buffer, how='inner', op='intersects')
                lst_temp = clipped_branches_gdf.code.tolist()
                lateral_branches = lst_temp
                if 'wagv_wl_6440_2668' in lateral_branches:
                    lateral_branches.remove('wagv_wl_6440_2668')
                    print('\n Removed wagv_wl_6440_2668 from list of branches that get lateral links')

            ## EINDE STUK PEPIJN 
            #      
            # if hasattr(model_config.FM.two_d, "clip_extent_path") and (model_config.FM.two_d.clip_extent_path is not None):
            #     clip_extent_gdf = gpd.read_file(model_config.FM.two_d.clip_extent_path, crs=CRS)
            #     clip_extent = clip_extent_gdf[clip_extent_gdf['geometry'].apply(lambda geom: isinstance(geom, Polygon))]    # Only select the ones with a polygon, because that can only be clipped
                
            #     # Apply selection criteria:
            #     #   1. all shapes that have 'meer, plas' in their 'typewater' column must be clipped (must remain in clip_extent)
            #     #   2. all shapes that meet the clip_selection criteria for the 'breedtekla' column must be clipped
            #     if hasattr(model_config.FM.two_d, "clip_selection") and (model_config.FM.two_d.clip_selection is not None):
            #         condition1 = (clip_extent['typewater'] == 'meer, plas')
            #         condition2 = (clip_extent['typewater'] == 'waterloop') & (clip_extent['breedtekla'].isin(model_config.FM.two_d.clip_selection))
            #         clip_extent = clip_extent[condition1 | condition2]

            #     # Ensure that only unique column names exist
            #     if 'index_left' in clip_extent.columns:
            #         clip_extent = clip_extent.rename(columns={'index_left': 'clip_index_left'})
            #     if 'index_right' in clip_extent.columns:
            #         clip_extent = clip_extent.rename(columns={'index_right': 'clip_index_right'})

            #     # Get the clip buffer, because it is not defined everywhere, to prevent errors set the buffer to 100m if not defined.
            #     if hasattr(model_config.FM.two_d, "clip_buffer"):
            #         clip_buffer = model_config.FM.two_d.clip_buffer
            #     else:
            #         clip_buffer = 100
            #         print('** Assumed default clip buffer of 100m because no buffer was given.')
                
            #     # Get a list of the branchid's that are clipped and now need 1D-2D lateral links
            #     all_branches_gdf = gpd.GeoDataFrame(self.hydamo.branches, geometry='geometry', crs=CRS)
            #     clip_extent_buffer = clip_extent.copy()
                            
            #     clip_extent_buffer = clip_extent_buffer.set_geometry('geometry')
            #     all_branches_gdf = all_branches_gdf.set_geometry('geometry')
            #     clip_extent_buffer = clip_extent_buffer.to_crs(CRS)
            #     all_branches_gdf = all_branches_gdf.to_crs(CRS)

            #     clipped_branches_gdf = gpd.sjoin(all_branches_gdf, clip_extent_buffer, how='inner', op='intersects')
            #     lst_temp = clipped_branches_gdf.code_left.tolist()
            #     lateral_branches = lst_temp 

            # # If no extent is giving, don't use it
            # else:
            #     clip_extent = None
            #     clip_buffer = 0
            #     lateral_branches = []
            
            
            # add 2D model
            print("\nBuilding 2D model grid\n")

            # Get the list of branches that should not get a 1D2D link because they are duikers
            all_branches_code_list = list(self.ddm.waterloop['code'])
            duiker_list = list(self.ddm.waterloop[self.ddm.waterloop['is_duiker'] =="JA"]['code'])
            branch_1d2d_embed_list = [item for item in all_branches_code_list if item not in duiker_list]
            fm_path = Path(output_folder) / FM_FOLDER
            self.fm = build_2D_model(
                dx=model_config.FM.two_d.dx,
                dy=model_config.FM.two_d.dy,
                extent=extent,
                fm=self.fm,
                fm_path = fm_path,
                one_d=model_config.FM.one_d_bool,
                clip_extent=clip_extent,
                clip_buffer=clip_buffer,
                load_mesh2d_path = load_mesh2d_path,
                lateral_branches=lateral_branches,
                branch_1d2d_embed_list = branch_1d2d_embed_list,
                **kwargs,
            )

            if ("keringen" in self.features) or ("overiglijnelement" in self.features):
                self.fm = add_fixed_weirs(ddm=self.ddm, fm=self.fm)

            if "doorbraak" in self.features:
                self.fm = add_dambreaks(ddm=self.ddm, fm=self.fm)

        if hasattr(model_config.FM, "hydrolib_core_options"):
            self.fm = set_hydrolib_core_options(
                fmmodel=self.fm, options=model_config.FM.hydrolib_core_options
            )


def write_dimr(fm: FMModel, output_folder: str):
    """ """

    output_path = Path(output_folder)
    output_path.mkdir(exist_ok=True, parents=True)
    fm_path = Path(output_folder) / FM_FOLDER
    fm.filepath = fm_path / "test.mdu"
    # for bound_idx in range(0, len(fm.external_forcing.extforcefilenew.boundary)):
    #     if fm.external_forcing.extforcefilenew.boundary[bound_idx].forcingfile.filepath is None:
    #         fm.external_forcing.extforcefilenew.boundary[bound_idx].forcingfile.filepath = fm_path / 'boundaryconditions.bc'
    dimr = DIMR()
    dimr.component.append(
        FMComponent(
            name="DFM", workingDir=output_path / "dflowfm", model=fm, inputfile=fm.filepath
        )
    )
    dimr.save(recurse=True)
    dimr_path = r"C:\Program Files\Deltares\D-HYDRO Suite 2023.01 1D2D\plugins\DeltaShell.Dimr\kernels\x64\dimr\scripts\run_dimr.bat"
    dimr = DIMRWriter(dimr_path=dimr_path, output_path=output_path)
    dimr.write_dimrconfig(fm=fm)  # , rr_model=drrmodel, rtc_model=drtcmodel)
    dimr.write_runbat()
