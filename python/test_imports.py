import sys
print("Python version:", sys.version)
print("Python path:", sys.executable)

try:
    import rasterio
    print("Rasterio version:", rasterio.__version__)
except Exception as e:
    print("Rasterio import error:", str(e))

try:
    from osgeo import gdal
    print("GDAL version:", gdal.__version__)
except Exception as e:
    print("GDAL import error:", str(e))