import os
import geopandas as gpd
import rasterio
from rasterio.features import shapes
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape
import numpy as np
import pyproj
from pysheds.grid import Grid
from pathlib import Path

def calculate_upstream_catchments(discharges_path, id_field, dem_path, channels_path, output_path):
    """
    Calculate upstream catchments for a set of given discharge points, considering watercourse paths.

    Parameters:
        discharges_path (str): Path to the shapefile with discharge points.
        id_field (str): Name of the ID field in the shapefile.
        dem_path (str): Path to the DEM GeoTIFF file.
        channels_path (str): Path to the shapefile with watercourse paths.
        output_path (str): Path to save the output shapefile.

    Returns:
        None
    """
    # Read discharge points
    discharge_points = gpd.read_file(discharges_path)

    # Read watercourse paths
    channels = gpd.read_file(channels_path)

    # Load DEM
    grid = Grid.from_raster(dem_path, data_name='dem')

    # Compute the flow direction grid
    grid.flowdir(data='dem', out_name='dir')

    catchments = []

    for _, point in discharge_points.iterrows():
        point_geom = point.geometry
        id_value = point[id_field]

        # Convert point geometry to grid coordinates
        x, y = point_geom.x, point_geom.y
        col, row = grid.nearest_cell(x, y)

        # Delineate the catchment upstream of the point
        grid.catchment(data='dir', x=col, y=row, dirmap=grid.dirmap, out_name='catch', xytype='index')

        # Mask and convert catchment grid to shapes
        shapes_generator = shapes(grid.view('catch').astype(np.uint8), mask=grid.view('catch') == 1, transform=grid.affine)
        for geom, value in shapes_generator:
            if value == 1:
                catchment_geom = shape(geom)

                # Check intersection with channels
                if channels.intersects(catchment_geom).any():
                    catchments.append({
                        'geometry': catchment_geom,
                        'id': id_value
                    })

    # Create GeoDataFrame for catchments
    catchments_gdf = gpd.GeoDataFrame(catchments, crs=discharge_points.crs)

    # Save the catchments as a shapefile
    catchments_gdf.to_file(output_path, driver='ESRI Shapefile')

if __name__ == "__main__":
    
    # Paths to input files and output
    discharges_path = r"c:\GITHUB\blokkendoosLimburg\data\knelpunt.shp"
    id_field = "id"
    dem_path = r"c:\GITHUB\blokkendoosLimburg\data\AHN_DTM_5M.tif"
    channels_path = r"c:\GITHUB\blokkendoosLimburg\data\waterlopen_WL.shp"
    output_path = r"c:\GITHUB\blokkendoosLimburg\data\knelpunt_strgeb.shp"

    calculate_upstream_catchments(discharges_path, id_field, dem_path, channels_path, output_path)
