import cmasher as cmr
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.constants import c
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyuvdata import UVData

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

# Load UVFITS file
uvd = UVData.from_file(
    "../data/uvfits/albatros_sim_zenith_uniform_beam_4_16MHz_24h.uvfits"
)
uvw_meters = uvd.uvw_array[..., np.newaxis]
wavelengths = c.value / uvd.freq_array
uvw_wavelengths = uvw_meters / wavelengths
uu, vv, ww = uvw_wavelengths.transpose(1, 0, 2).reshape(3, -1)

# Create figure and main u-v plot
fig, ax_uv = plt.subplots(figsize=(6, 5))
gridsize = 256
vmax = 5e2
vmin = 1e-1

# Create side plots using make_axes_locatable
divider = make_axes_locatable(ax_uv)
ax_uw = divider.append_axes("bottom", size="25%", pad=0.05, sharex=ax_uv)
ax_wv = divider.append_axes("left", size="25%", pad=0.05, sharey=ax_uv)


# Function to compute density
def compute_density(uu, vv, gridsize, ax):
    hb = ax.hexbin(
        uu, vv, gridsize=gridsize, cmap=cmr.pride, norm=LogNorm(vmax=vmax, vmin=vmin)
    )

    # Compute hexagon radius and area
    xlim = ax.get_xlim()
    width = xlim[1] - xlim[0]
    hex_radius = width / (gridsize * np.sqrt(3))
    hex_area = (3 * np.sqrt(3) / 2) * hex_radius**2

    # Convert hex counts to density
    hex_counts = hb.get_array()
    hex_density = hex_counts / hex_area
    hb.set_array(hex_density)

    return hb


# Main u-v coverage plot
hb = compute_density(uu, vv, gridsize, ax_uv)
ax_uv.set_aspect("equal")

# u-w plot (below main)
hb_uw = compute_density(uu, ww, gridsize, ax_uw)
ax_uw.set_xlabel("$u$ [$\lambda$]")
ax_uw.set_ylabel("")
ax_uw.set_yticks([-100, 0])

# w-v plot (left of main)
hb_wv = compute_density(ww, vv, gridsize, ax_wv)
ax_wv.set_xlabel("")
ax_wv.set_ylabel("$v$ [$\lambda$]")
ax_wv.set_xticks([-100, 0])

# Add colorbar covering all plots
cbar = plt.colorbar(
    hb,
    ax=[ax_uv, ax_uw, ax_wv],
    orientation="vertical",
    pad=0.02,
    fraction=0.1,
    aspect=28,
    extend="both",
)
cbar.set_label("Density [count/$\lambda^2$]")

# Axis labels
fig.text(0.18, 0.18, "$w$ [$\lambda$]", va="center", rotation=45)

# Title
fig.suptitle("Albatros Sim $uvw$ coverage: 24h@60s, 4-16MHz@100kHz", y=0.93)
plt.savefig("../data/plots/albatros_uvw.png", dpi=300)
# plt.show()
