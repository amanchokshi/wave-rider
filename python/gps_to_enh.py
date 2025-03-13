import csv
from typing import Dict, Union

import cmasher as cmr
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from astropy.coordinates import EarthLocation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage

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


def gps_to_enu(coords: Dict[Union[int, str], list]) -> np.ndarray:
    """
    Converts set of GPS coords (lat, lon, height) into ENU relative to first entry

    The first antenna coordinate in the `coords` dictionary is considered the
    reference point. All other coordinates are converted to ENU relative to
    this reference

    Parameters
    ----------
    coords : dict
        A dictionary where keys are indices (int or str) and values are lists
        of [latitude (deg), longitude (deg), height (m)].

    Returns
    -------
    np.ndarray
        An (n, 3) ndarray containing ENU coordinates relative to first location
    """

    # Convert Geodetic (lat, lon, height) to ITRS Cartesian (ECEF)
    antenna_locations = {
        key: EarthLocation.from_geodetic(
            lat=lat, lon=lon, height=height
        ).itrs.cartesian.xyz.value
        for key, (lat, lon, height) in coords.items()
    }

    # Identify the reference location (first entry in the dictionary)
    ref_key = next(iter(coords))
    ref_location = antenna_locations[ref_key]

    # Compute baselines in ITRS (relative to reference antenna)
    baselines_itrs = {key: loc - ref_location for key, loc in antenna_locations.items()}

    # Define ENU Coordinate System at Reference Antenna, and normalize
    up_vector = ref_location / np.linalg.norm(ref_location)

    # Compute and normalize East vector (unit vector in the local eastward direction)
    lon_0_rad = np.deg2rad(coords[ref_key][1])
    east_vector = np.array([-np.sin(lon_0_rad), np.cos(lon_0_rad), 0])
    east_vector /= np.linalg.norm(east_vector)

    # Compute and normalize North vector (cross product of Up and East)
    north_vector = np.cross(up_vector, east_vector)
    north_vector /= np.linalg.norm(north_vector)

    # Convert ITRS Baselines to ENU
    baselines_enu = np.array(
        [
            np.array(
                [
                    np.dot(baseline, east_vector),  # East component
                    np.dot(baseline, north_vector),  # North component
                    np.dot(baseline, up_vector),  # Up component
                ]
            )
            for baseline in baselines_itrs.values()
        ]
    )

    return baselines_enu


if __name__ == "__main__":
    # List of coordinates (lat[deg], lon[deg], altitude[m])
    gps_coords = {
        "ANT01": [79 + 25.031 / 60, -90 - 46.041 / 60, 189],
        "ANT02": [79 + 25.033 / 60, -90 - 45.531 / 60, 176],
        "ANT03": [79 + 24.925 / 60, -90 - 46.385 / 60, 175],
        "ANT04": [79 + 23.308 / 60, -91 - 1.156 / 60, 22],
        "ANT05": [79 + 25.099 / 60, -90 - 40.048 / 60, 53],
        "ANT06": [79 + 23.880 / 60, -90 - 47.996 / 60, 43],
    }

    enu_coords = gps_to_enu(gps_coords)

    # Writing to the CSV file
    with open("../data/layouts/albatros.csv", mode="w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")  # Use tab as the delimiter
        # Write the header
        writer.writerow(["Name", "Number", "BeamID", "E", "N", "U"])

        # Write the data with f-string formatting for 4 significant figures
        for i, row in enumerate(enu_coords):
            formatted_row = [
                list(gps_coords.keys())[i],  # Name
                i,  # Beam number
                0,  # Beam type?
                f"{row[0]}",  # E
                f"{row[1]}",  # N
                f"{row[2]}",  # U
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
    im = ax.imshow(
        smoothed_img, extent=extent, cmap=cmr.pride, origin="upper", alpha=0.7
    )

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
    # for ant, (e, n, _) in enh_offsets.items():
    for i, row in enumerate(enu_coords):
        e, n, _ = row
        ant = (list(gps_coords.keys())[i],)  # Name
        ax.scatter(
            e, n, marker=r"$✠$", color="black", ec="white", lw=0.03, s=63, zorder=42
        )
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
