import os

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
from pyradiosky import SkyModel

albatros = EarthLocation(lat=79.417183, lon=-90.767350, height=189)
time = Time("2025-01-01 00:00:00", scale="utc", location=albatros)

############################
# Zenith Source
############################
source_coord = SkyCoord(
    alt=Angle(90, unit=u.deg),
    az=Angle(0, unit=u.deg),
    obstime=time,
    frame="altaz",
    location=albatros,
)
stokes = [1.0, 0, 0, 0] * u.Jy
icrs_coord = source_coord.transform_to("icrs")
sm = SkyModel(
    name=["ZEN"],
    skycoord=icrs_coord,
    stokes=stokes,
    spectral_type="flat",
    history="1Jy unpolarised source at zenith from Albatros site at 20250101 12:00:00",
)
write_file = os.path.join("../data/srclists/", "zenith.txt")
sm.write_text_catalog(write_file)

############################
# Off Zenith Source
############################
source_coord = SkyCoord(
    alt=Angle(80, unit=u.deg),
    az=Angle(0, unit=u.deg),
    obstime=time,
    frame="altaz",
    location=albatros,
)
stokes = [1.0, 0, 0, 0] * u.Jy
icrs_coord = source_coord.transform_to("icrs")
sm = SkyModel(
    name=["10N"],
    skycoord=icrs_coord,
    stokes=stokes,
    spectral_type="flat",
    history="1Jy unpolarised source 10deg North of zenith from Albatros site at 20250101 12:00:00",
)
write_file = os.path.join("../data/srclists/", "off_zenith.txt")
sm.write_text_catalog(write_file)


############################
# North Sources
############################
source_coord = SkyCoord(
    alt=Angle([90, 85, 80, 75, 70, 65, 60], unit=u.deg),
    az=Angle([0, 0, 0, 0, 0, 0, 0], unit=u.deg),
    obstime=time,
    frame="altaz",
    location=albatros,
)
stokes = np.zeros((4, 1, 7)) * u.Jy
stokes[0, ...] = 1.0 * u.Jy
icrs_coord = source_coord.transform_to("icrs")
sm = SkyModel(
    name=["ZEN", "5N", "10N", "15N", "20N", "25N", "30N"],
    skycoord=icrs_coord,
    stokes=stokes,
    spectral_type="flat",
    history="Line of 1Jy unpolarised source at zenith from Albatros site at 20250101 12:00:00",
)
write_file = os.path.join("../data/srclists/", "north_line.txt")
sm.write_text_catalog(write_file)
