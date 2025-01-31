import imod
import os
import xarray as xr
import numpy as np
from osgeo import gdal

# Definieer paden
input_file = r"c:\SYNC\PROJECTEN\H3147.WRL\02.GIS\Landgebruik\SOPP.idf"
output_file = r"c:\SYNC\PROJECTEN\H3147.WRL\02.GIS\Landgebruik\SOPP.tif"

# Check of input bestaat
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input bestand niet gevonden: {input_file}")

# Lees het IDF bestand
data, *metadata = imod.idf.read(input_file)
meta = metadata[0]  # Pak de metadata dictionary

# Creëer een GeoTIFF
driver = gdal.GetDriverByName('GTiff')
out_ds = driver.Create(output_file, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)

# Stel de geotransformatie in
# (linksboven x, pixel breedte, 0, linksboven y, 0, pixel hoogte)
geotransform = (meta['xmin'], meta['dx'], 0, meta['ymax'], 0, meta['dy'])
out_ds.SetGeoTransform(geotransform)

# Schrijf de data
out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(data)

# Schoon op
out_ds = None

print(f"Conversie voltooid: {output_file}")