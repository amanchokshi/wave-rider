import cmasher as cmr
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
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
uvd = UVData.from_file("../data/uvfits/albatros_sim_diffuse_ateam.uvfits")

enu = uvd.get_enu_data_ants()


def baseline_length_meters(enu, ant1, ant2):
    return np.linalg.norm(enu[0][ant1] - enu[0][ant2])


for ant1 in enu[1]:
    for ant2 in enu[1]:
        if ant2 > ant1:
            print(
                f"V_{ant1}{ant2}: [ {baseline_length_meters(enu, ant1, ant2):8.3f}m ]"
            )

            waterfall_data = uvd.get_data((ant1, ant2, uvd.polarization_array[0]))
            waterfall_times = Time(
                uvd.get_times((1, 2, uvd.polarization_array[0])), format="jd"
            ).to_datetime()
            waterfall_freqs = uvd.freq_array / 1e6

            extent = [
                waterfall_freqs[0],
                waterfall_freqs[-1],
                mdates.date2num(waterfall_times[-1]),
                mdates.date2num(waterfall_times[0]),
            ]

            fig, ax = plt.subplots(figsize=(7, 4))
            im = ax.imshow(
                np.abs(waterfall_data), cmap=cmr.pride, aspect="auto", extent=extent
            )

            ax.yaxis_date()
            ax.yaxis.set_major_locator(mdates.HourLocator(interval=1))
            ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

            ax.set_xlabel("Frequency [ MHz ]")
            ax.set_ylabel("Time [ HH:MM ]")
            ax.set_title(
                rf"$\mathcal{{V}}_{{{ant1}{ant2}}}$ : {{[ {baseline_length_meters(enu, ant1, ant2):6.1f}m, 24h@60s, 4-16MHz@100kHz ]}}"
            )
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("Visibility Amplitude")

            plt.tight_layout()
            plt.savefig(
                f"../data/plots/diffuse_ateam/diffuse_ateam_v_{ant1}{ant2}_amp_{baseline_length_meters(enu, ant1, ant2):06.1f}m.png",
                dpi=300,
            )
            # plt.show()
            plt.close()

            fig, ax = plt.subplots(figsize=(7, 4))
            im = ax.imshow(
                np.angle(waterfall_data), cmap=cmr.pride, aspect="auto", extent=extent
            )

            ax.yaxis_date()
            ax.yaxis.set_major_locator(mdates.HourLocator(interval=1))
            ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

            ax.set_xlabel("Frequency [ MHz ]")
            ax.set_ylabel("Time [ HH:MM ]")
            ax.set_title(
                rf"$\mathcal{{V}}_{{{ant1}{ant2}}}$ : {{[ {baseline_length_meters(enu, ant1, ant2):6.1f}m, 24h@60s, 4-16MHz@100kHz ]}}"
            )
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label("Visibility Phase")

            plt.tight_layout()
            plt.savefig(
                f"../data/plots/diffuse_ateam/diffuse_ateam_v_{ant1}{ant2}_phase_{baseline_length_meters(enu, ant1, ant2):06.1f}m.png",
                dpi=300,
            )
            # plt.show()
            plt.close()
