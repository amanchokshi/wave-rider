import sys

import astropy.units as units
import healpy as hp
import numpy as np
import pyradiosky
from astropy.units import Quantity
from pygdsm import GlobalSkyModel


def generate_pygdsm_map(
    freq_mhz,
    nside,
    output_frame="equatorial",  # Convert to celestial if needed
    spec_index=-2.5,
):
    """
    Generate a HEALPix sky map using pyGDSM at a given frequency and nside.
    Converts from Galactic to Equatorial (ICRS) if required.
    """
    gsm = GlobalSkyModel(freq_unit="MHz")
    hpx_gsm = gsm.generate(freq_mhz)  # HEALPix map in Galactic coordinates
    hpx_gsm_nside = hp.pixelfunc.ud_grade(hpx_gsm, nside, pess=True)

    history_str = f"Generated using pyGDSM at {freq_mhz} MHz, nside={nside}"

    # Default coordinate system from pyGDSM
    coordsys = "G"

    # Convert to Equatorial (Celestial) if requested
    if output_frame == "equatorial":
        print("Converting map from Galactic to Equatorial (ICRS) frame...")

        # Get pixel theta, phi in Galactic frame
        npix = hp.nside2npix(nside)
        theta_gal, phi_gal = hp.pix2ang(nside, np.arange(npix), nest=False)

        # Define Galactic → Equatorial rotation
        rot = hp.rotator.Rotator(coord=["G", "C"], inv=True)  # Inverse rotation (G → C)

        # Apply rotation
        theta_eq, phi_eq = rot(theta_gal, phi_gal)

        # Interpolate HEALPix map at new positions
        hpx_gsm_nside = hp.get_interp_val(hpx_gsm_nside, theta_eq, phi_eq)

        history_str += ", Rotated from Galactic to Equatorial frame"
        coordsys = "C"

    # Define the frame string for pyradiosky
    if coordsys == "C":
        frame_str = "icrs"
    elif coordsys == "G":
        frame_str = "galactic"
    else:
        sys.exit("ERROR: Unsupported coordsys.")

    # Create pyradiosky SkyModel
    skymodel = pyradiosky.SkyModel()
    skymodel.component_type = "healpix"
    skymodel.nside = nside
    skymodel.hpx_order = "ring"  # Default order for pyGDSM
    skymodel.frame = frame_str
    skymodel.hpx_frame = frame_str
    skymodel.Nfreqs = 1
    skymodel.Ncomponents = hp.nside2npix(nside)
    skymodel.stokes = Quantity(
        np.zeros((4, skymodel.Nfreqs, skymodel.Ncomponents)), "Kelvin"
    )
    skymodel.stokes[0, 0, :] = hpx_gsm_nside * units.Kelvin
    skymodel.hpx_inds = np.arange(skymodel.Ncomponents)
    skymodel.spectral_type = "spectral_index"
    skymodel.reference_frequency = Quantity(
        np.full(skymodel.Ncomponents, freq_mhz * 1e6), "hertz"
    )
    skymodel.spectral_index = np.full(skymodel.Ncomponents, spec_index)
    skymodel.history = history_str

    return skymodel


if __name__ == "__main__":
    freq_mhz = 14
    nside = 1024
    output_skyh5 = f"./gsm_{freq_mhz}MHz_nside{nside}_equatorial.skyh5"

    skymodel_new = generate_pygdsm_map(freq_mhz, nside, output_frame="equatorial")

    skymodel_new.write_skyh5(
        output_skyh5,
        run_check=True,
        clobber=True,
    )

    print(f"Saved pyGDSM sky model to {output_skyh5}")
