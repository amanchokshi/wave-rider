import csv

import cmasher as cmr
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage
from skyfield.api import wgs84

# Set font styles for plots
plt.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsfonts}",
        "font.family": "serif",
        "axes.labelsize": 9,
        "axes.titlesize": 11,
        "font.size": 8,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    }
)

# List of coordinates (lat[deg], lon[deg], altitude[m])
coords = [
    ("ANT01", 79.417183, -90.767350, 189),
    ("ANT02", 79.417217, -90.758850, 176),
    ("ANT03", 79.415417, -90.773083, 175),
    ("ANT04", 79.388467, -91.019267, 22),
    ("ANT05", 79.418317, -90.667467, 53),
    ("ANT06", 79.398000, -90.799933, 43),
]

# Convert coordinates to Skyfield wgs84 objects
gps_objects = {name: wgs84.latlon(lat, lon, alt) for name, lat, lon, alt in coords}

# Compute ENH offsets relative to ANT01
ref_position = gps_objects["ANT01"].itrs_xyz.m
enh_offsets = {ant: gps.itrs_xyz.m - ref_position for ant, gps in gps_objects.items()}

antenna_data = []
for ant, enh in enh_offsets.items():
    antenna_data.append([ant, enh[0], enh[1], enh[2]])


# Writing to the CSV file
with open("../data/layouts/albatros.csv", mode="w", newline="") as file:
    writer = csv.writer(file, delimiter="\t")  # Use tab as the delimiter
    # Write the header
    writer.writerow(["Name", "Number", "BeamID", "E", "N", "U"])

    # Write the data with f-string formatting for 4 significant figures
    for i, row in enumerate(antenna_data):
        formatted_row = [
            row[0],  # Name
            i,  # Beam number
            0,  # Beam type?
            f"{row[1]:.4f}",  # E with 4 significant figures
            f"{row[2]:.4f}",  # N with 4 significant figures
            f"{row[3]:.4f}",  # U with 4 significant figures
        ]
        writer.writerow(formatted_row)


# Load the cropped GeoTIFF file
with rasterio.open("../data/dem/axel_heiberg.tif") as src:
    img = ndimage.zoom(src.read(1), 3)  # Read and scale the first band
    smoothed_img = ndimage.gaussian_filter(img, sigma=14)
    left, bottom, right, top = src.bounds
    center_x, center_y = (left + right) / 2, (bottom + top) / 2
    extent = [left - center_x, right - center_x, bottom - center_y, top - center_y]

# Plot the image with contours
fig, ax = plt.subplots(figsize=(7, 4))
im = ax.imshow(smoothed_img, extent=extent, cmap=cmr.pride, origin="upper", alpha=0.7)

contour_levels = np.arange(15, 270, 15)
ax.contour(
    smoothed_img,
    extent=extent,
    levels=contour_levels,
    origin="upper",
    cmap=cmr.pride,
    linewidths=0.3,
    alpha=1.0,
)

# Use make_axes_locatable to create a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Elevation (m)")

# Plot antenna positions
plt.rcParams["text.usetex"] = False  # Disable LaTeX for markers
for ant, (e, n, _) in enh_offsets.items():
    ax.scatter(e, n, marker=r"$✠$", color="black", ec="white", lw=0.03, s=63, zorder=42)
    ax.text(
        e - 540,
        n + 90,
        ant,
        fontsize=6,
        ha="left",
        va="bottom",
        color="black",
        rotation=-45,
        zorder=63,
    )

# Formatting and labels
ax.set_aspect("equal")
ax.set_xlabel("East Distance from Center (m)")
ax.set_ylabel("North Distance from Center (m)")
ax.set_title("Albatros Antennas on Axel Heiberg Island")
ax.set_xlim([-5900, 2500])
ax.set_ylim([-3600, 900])

plt.tight_layout()
plt.savefig("../data/plots/albatros_antennas.png", dpi=300)
plt.show()
