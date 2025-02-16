import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

# Define the input GeoTIFF file and output file
input_tif = "../data/dem/cdem-canada-hillshade.tif"  # Replace with actual file path
output_tif = "../data/dem/axel_heiberg.tif"

# Define the center GPS position (latitude, longitude)
center_lat, center_lon = 79.417183, -90.767350


# Define the crop size in meters
half_size = 6000

# Open the original raster file
with rasterio.open(input_tif) as src:
    # Get the CRS of the raster (assume it's in a projected coordinate system like UTM)
    raster_crs = src.crs

    # Transform GPS lat/lon to the raster's coordinate system
    transformer = Transformer.from_crs("EPSG:4326", raster_crs, always_xy=True)
    center_x, center_y = transformer.transform(center_lon, center_lat)

    # Define the bounding box (min_x, min_y, max_x, max_y)
    min_x, max_x = center_x - half_size, center_x + half_size
    min_y, max_y = center_y - half_size, center_y + half_size

    # Get the window corresponding to the bounding box
    window = from_bounds(min_x, min_y, max_x, max_y, src.transform)

    # Read the data from the window
    data = src.read(window=window)

    # Update the transform for the cropped area
    new_transform = src.window_transform(window)

    # Write the cropped section to a new GeoTIFF file
    profile = src.profile
    profile.update(
        {"height": data.shape[1], "width": data.shape[2], "transform": new_transform}
    )

    with rasterio.open(output_tif, "w", **profile) as dst:
        dst.write(data)

print(f"Cropped GeoTIFF saved as {output_tif}")
