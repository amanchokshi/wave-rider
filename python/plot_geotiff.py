import cmasher as cmr
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer
from scipy.ndimage import gaussian_filter

# Define the cropped GeoTIFF file
cropped_tif = "../data/dem/axel_heiberg.tif"

# Open the cropped raster
with rasterio.open(cropped_tif) as src:
    # Read image data (assuming single band for simplicity)
    # print(dir(src))
    img = src.read(1)  # Read the first band

    # Get raster's CRS
    raster_crs = src.crs

    # Get the bounds in the raster's CRS
    left, bottom, right, top = src.bounds

    # Compute the center of the image in the raster CRS
    center_x, center_y = (left + right) / 2, (bottom + top) / 2

    # Convert bounds to distance in kilometers from the center
    def meters_to_km(x, y):
        return (x - center_x) / 1000, (y - center_y) / 1000

    x_left, y_bottom = meters_to_km(left, bottom)
    x_right, y_top = meters_to_km(right, top)

    # Define the extent in km
    extent = [x_left, x_right, y_bottom, y_top]

# Plot the image with lat/lon axis labels
# contour_levels = [0, 50, 100, 150, 200]
contour_levels = np.arange(0, 275, 25)

# **Apply Gaussian filter for smoothing**
sigma = 7  # Adjust for more or less smoothing
smoothed_img = gaussian_filter(img, sigma=sigma)

fig, ax = plt.subplots()
im = ax.imshow(smoothed_img, extent=extent, cmap=cmr.pride, origin="upper", alpha=0.4)
contour = ax.contour(
    smoothed_img,
    extent=extent,
    levels=contour_levels,
    origin="upper",
    cmap=cmr.pride,
    linewidths=0.3,
)
ax.set_xlabel("East Distance from Center (km)")
ax.set_ylabel("North Distance from Center (km)")
ax.set_title("Cropped GeoTIFF with Latitude/Longitude")
fig.colorbar(im, label="Elevation (m)")
# plt.grid(True)
plt.show()
