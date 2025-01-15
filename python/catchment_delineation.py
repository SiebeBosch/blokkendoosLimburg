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
    Calculate upstream catchments for a set of given discharge polygons, considering watercourse paths.

    Parameters:
        discharges_path (str): Path to the shapefile with discharge polygons.
        id_field (str): Name of the ID field in the shapefile.
        dem_path (str): Path to the DEM GeoTIFF file.
        channels_path (str): Path to the shapefile with watercourse paths.
        output_path (str): Path to save the output shapefile.

    Returns:
        None
    """
    # Read discharge polygons
    discharge_areas = gpd.read_file(discharges_path)

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

    # Get the DEM's transform and metadata
    with rasterio.open(dem_path) as dem_src:
        dem_transform = dem_src.transform
        dem_crs = dem_src.crs
        dem_shape = dem_src.shape

    catchments = []

    # Ensure discharge areas are in the same CRS as the DEM
    if discharge_areas.crs != dem_crs:
        discharge_areas = discharge_areas.to_crs(dem_crs)

    print("Processing discharge areas...")
    for _, area in discharge_areas.iterrows():
        area_geom = area.geometry
        id_value = area[id_field]
        print(f"Processing area with ID {id_value}...")

        # Create a mask for the current polygon
        mask = np.zeros(dem_shape, dtype=bool)
        
        # Get bounds of the polygon in pixel coordinates
        minx, miny, maxx, maxy = area_geom.bounds
        col_min, row_max = ~dem_transform * (minx, miny)
        col_max, row_min = ~dem_transform * (maxx, maxy)
        
        # Convert to integers and ensure within bounds
        col_min = max(0, int(col_min))
        col_max = min(dem_shape[1], int(col_max) + 1)
        row_min = max(0, int(row_min))
        row_max = min(dem_shape[0], int(row_max) + 1)

        # Create a grid of coordinates for the bounding box
        rows, cols = np.mgrid[row_min:row_max, col_min:col_max]
        x_coords, y_coords = dem_transform * (cols.flatten(), rows.flatten())
        
        # Create points and check which ones are within the polygon
        from shapely.geometry import Point
        points = [Point(x, y) for x, y in zip(x_coords, y_coords)]
        for point, row, col in zip(points, rows.flatten(), cols.flatten()):
            if area_geom.contains(point):
                mask[row, col] = True

        try:
            # Initialize combined catchment
            combined_catch = np.zeros_like(dem, dtype=bool)
            
            # Get all cells within the polygon
            start_cells = np.column_stack(np.where(mask))
            
            # Delineate catchment from each cell
            for row, col in start_cells:
                print(f"Delineating catchment for cell ({row}, {col}) in area {id_value}...")
                catch = grid.catchment(fdir=fdir, x=col, y=row, dirmap=dirmap, xytype='index')
                combined_catch = np.logical_or(combined_catch, catch == 1)

            # Print the number of cells in the combined catchment
            num_cells = np.sum(combined_catch)
            print(f"Combined catchment for area ID {id_value} has {num_cells} cells")

            if num_cells > 0:
                # Mask and convert catchment grid to shapes
                shapes_generator = shapes(combined_catch.astype(np.uint8), 
                                       mask=combined_catch, 
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
            print(f"Warning: Could not process area with ID {id_value}: {str(e)}")

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
    # Notice that the discharge locations are polygons, not points
    discharges_path = r"c:\GITHUB\blokkendoosLimburg\data\knelpunt_polygon.shp"
    id_field = "NAAM"
    dem_path = r"c:\GITHUB\blokkendoosLimburg\data\eudem_dem_4258_europe_clip.tif"
    channels_path = r"c:\GITHUB\blokkendoosLimburg\data\waterlopen_WL.shp"
    output_path = r"c:\GITHUB\blokkendoosLimburg\data\strgeb_knelpunt_EU_DEM.shp"

    calculate_upstream_catchments(discharges_path, id_field, dem_path, channels_path, output_path)
