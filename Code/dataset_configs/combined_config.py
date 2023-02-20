class Models:
    class FM:
        one_d_bool = True
        two_d_bool = True
        start_time = 20000101
        stop_time = 86400

        class one_d:
            max_dist_to_struct = 3
            max_snap_dist = 1
            node_distance = 500

        class two_d:
            coupling_type = "2Dto1D"
            dx = 500
            dy = 500
            elevation_raster_path = "D:\Work\Project\P1414\GIS\AHN\AHN_merged.TIF"
            extent_path = "D:\Work\Project\P1414\GIS\Randstad_shape\dijkringen_randstad_merged.shp"
            two_d_buffer = 100

        class hydrolib_core_options:
            class external_forcing:
                pass

            class geometry:
                dxmin1d = 500
                usecaching = 1

            class numerics:
                cflmax = 0.7

            class output:
                hisinterval = [0]

            class time:
                dtmax = 60
