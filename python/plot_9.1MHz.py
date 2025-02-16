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


def extract_visi(uvfits_path, ant1, ant2, freq):
    """Extract visibility data for a given frequency and antenna pair."""
    uvd = UVData()
    uvd.read(uvfits_path)

    # Ensure frequency is properly indexed
    freqs = uvd.freq_array.flatten()  # Flatten in case freq_array is 2D
    freq_idx = np.argmin(np.abs(freqs - freq))  # Find nearest frequency index

    # Extract visibility data
    pol = uvd.polarization_array[0]  # Assuming single polarization
    times = Time(uvd.get_times((ant1, ant2, pol)), format="jd").to_datetime()
    visi = uvd.get_data((ant1, ant2, pol))[:, freq_idx]

    return visi, times


# Define parameters
ant1, ant2 = 0, 1
freq = 9.1e6  # Hz

# Extract visibility data
ateam, times = extract_visi(
    "../data/uvfits/albatros_sim_ateam.uvfits", ant1, ant2, freq
)
diff, _ = extract_visi("../data/uvfits/albatros_sim_diffuse.uvfits", ant1, ant2, freq)
diff_ateam, _ = extract_visi(
    "../data/uvfits/albatros_sim_diffuse_ateam.uvfits", ant1, ant2, freq
)


c = cmr.pride([0.7, 0.3, 0.4])

# Plot results
fig, ax = plt.subplots(figsize=(7, 4))

ax.plot(times, np.unwrap(np.angle(ateam)), label="A-Team", color=c[0])
ax.plot(times, np.unwrap(np.angle(diff)), label="Diffuse", color=c[1])
ax.plot(times, np.unwrap(np.angle(diff_ateam)), label="Diffuse+A-Team", color=c[2])

# Format time axis
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_xlabel("Time (UTC)")
ax.set_ylabel("Visibility Phase")
ax.legend()

plt.tight_layout()
plt.savefig("../data/plots/174m/unwrapped_phase_9.1MHz_174m.png", dpi=300)
# plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(7, 4))

ax.plot(times, np.abs(ateam), label="A-Team", color=c[0])
ax.plot(times, np.abs(diff), label="Diffuse", color=c[1])
ax.plot(times, np.abs(diff_ateam), label="Diffuse+A-Team", color=c[2])

# Format time axis
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax.set_xlabel("Time (UTC)")
ax.set_ylabel("Visibility Amplitude")
ax.legend()

plt.tight_layout()
plt.savefig("../data/plots/174m/amp_9.1MHz_174m.png", dpi=300)
# plt.show()
