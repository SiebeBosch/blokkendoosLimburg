import os
import sys
from pathlib import Path

# Set up GDAL environment variables
if os.name == 'nt':  # Windows
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        # Add Conda Library paths to system PATH
        lib_path = str(Path(conda_prefix) / 'Library/bin')
        if lib_path not in os.environ['PATH']:
            os.environ['PATH'] = lib_path + os.pathsep + os.environ['PATH']
        
        # Set GDAL specific environment variables
        os.environ['GDAL_DATA'] = str(Path(conda_prefix) / 'Library/share/gdal')
        os.environ['PROJ_LIB'] = str(Path(conda_prefix) / 'Library/share/proj')
        
        # Additional environment variables that might help with DLL loading
        os.environ['GDAL_DRIVER_PATH'] = str(Path(conda_prefix) / 'Library/lib/gdalplugins')
        os.environ['GDAL_VERSION'] = '3.6.2'  # Match version in environment.yml

# Now import the rest of the packages
try:
    from osgeo import gdal
except ImportError as e:
    print("Error loading GDAL:", str(e))
    sys.exit(1)

try:
    import rasterio
except ImportError as e:
    print("Error loading rasterio:", str(e))
    sys.exit(1)

try:
    import geopandas as gpd
except ImportError as e:
    print("Error loading geopandas:", str(e))
    sys.exit(1)

from rasterio.features import shapes
from rasterio.warp import reproject, Resampling
from shapely.geometry import shape
import numpy as np
import pyproj
from pysheds.grid import Grid

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
    grid = Grid.from_raster(dem_path)
    
    # Read DEM data
    dem = grid.read_raster(dem_path)
    
    # Condition DEM
    print("Conditioning DEM...")
    # Fill pits in DEM
    pit_filled_dem = grid.fill_pits(dem)
    # Fill depressions in DEM
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    # Resolve flats in DEM
    inflated_dem = grid.resolve_flats(flooded_dem)
    
    grid.dem = inflated_dem

    # Define flow direction mapping (D8 flow directions)
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

    print("Calculating flow direction...")
    # Compute the flow direction grid
    fdir = grid.flowdir(inflated_dem, dirmap=dirmap)

    # Get the DEM's transform
    with rasterio.open(dem_path) as dem_src:
        dem_transform = dem_src.transform
        dem_crs = dem_src.crs

    catchments = []

    # Ensure discharge points are in the same CRS as the DEM
    if discharge_points.crs != dem_crs:
        discharge_points = discharge_points.to_crs(dem_crs)

    print("Processing discharge points...")
    for _, point in discharge_points.iterrows():
        point_geom = point.geometry
        id_value = point[id_field]

        # Convert point geometry to grid coordinates using the DEM's transform
        col, row = ~dem_transform * (point_geom.x, point_geom.y)
        col, row = int(col), int(row)

        # Check if the point is within bounds
        if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1]:
            try:
                print(f"Delineating catchment for point ID {id_value}...")
                # Delineate the catchment upstream of the point
                catch = grid.catchment(fdir=fdir, x=col, y=row, dirmap=dirmap, xytype='index')
                
                # Print the number of cells in the catchment
                num_cells = np.sum(catch == 1)
                print(f"Catchment for point ID {id_value} has {num_cells} cells")

                # Mask and convert catchment grid to shapes
                shapes_generator = shapes(catch.astype(np.uint8), 
                                       mask=catch == 1, 
                                       transform=grid.affine)
                for geom, value in shapes_generator:
                    if value == 1:
                        catchment_geom = shape(geom)

                        # Check intersection with channels
                        if channels.intersects(catchment_geom).any():
                            catchments.append({
                                'geometry': catchment_geom,
                                'id': id_value
                            })
            except Exception as e:
                print(f"Warning: Could not process point with ID {id_value}: {str(e)}")
        else:
            print(f"Warning: Point with ID {id_value} is outside DEM bounds")

    if not catchments:
        print("Warning: No valid catchments were generated")
        return

    print(f"Generated {len(catchments)} catchments")
    
    # Create GeoDataFrame for catchments
    catchments_gdf = gpd.GeoDataFrame(catchments, crs=dem_crs)

    # Save the catchments as a shapefile
    catchments_gdf.to_file(output_path, driver='ESRI Shapefile')
    print(f"Saved catchments to {output_path}")

if __name__ == "__main__":
    
    # Paths to input files and output
    discharges_path = r"c:\GITHUB\blokkendoosLimburg\data\knelpunt.shp"
    id_field = "id"
    dem_path = r"c:\GITHUB\blokkendoosLimburg\data\AHN_DTM_25M.tif"
    channels_path = r"c:\GITHUB\blokkendoosLimburg\data\waterlopen_WL.shp"
    output_path = r"c:\GITHUB\blokkendoosLimburg\data\knelpunt_strgeb_25m.shp"

    calculate_upstream_catchments(discharges_path, id_field, dem_path, channels_path, output_path)
