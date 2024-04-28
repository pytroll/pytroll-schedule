#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2019 PyTroll community

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""The SwathBoundary class."""


import logging
import logging.handlers

import numpy as np
from pyorbital import geoloc, geoloc_instrument_definitions
from pyresample.boundary import Boundary

logger = logging.getLogger(__name__)

INSTRUMENT = {"avhrr/3": "avhrr",
              "avhrr/2": "avhrr",
              "avhrr-3": "avhrr",
              "mwhs-2": "mwhs2"}


class SwathBoundary(Boundary):
    """Boundaries for satellite overpasses."""

    def get_instrument_points(self, overpass, utctime,
                              scans_nb, scanpoints, scan_step=1):
        """Get the boundary points for a given overpass."""
        instrument, scan_angle = self.get_instrument_and_angle(overpass)

        sgeom = self.create_instrument_geometry(instrument, scans_nb, scanpoints, scan_step, scan_angle)

        times = sgeom.times(utctime)

        pixel_pos = geoloc.compute_pixels((self.orb.tle._line1,
                                           self.orb.tle._line2),
                                          sgeom, times)
        lons, lats, alts = geoloc.get_lonlatalt(pixel_pos, times)

        del alts
        return (lons.reshape(-1, len(scanpoints)),
                lats.reshape(-1, len(scanpoints)))

    def get_instrument_and_angle(self, overpass):
        """Get the instrument and angle for an overpass."""
        instrument = overpass.instrument
        if instrument == "modis":
            scan_angle = 55.0
            instrument = "avhrr"
        elif instrument == "viirs":
            scan_angle = 55.84
            instrument = "viirs"
        elif instrument == "iasi":
            scan_angle = 48.3
            instrument = "avhrr"
        elif overpass.satellite == "noaa 16":
            scan_angle = 55.25
            instrument = "avhrr"
        elif instrument.startswith("mersi"):
            scan_angle = 55.4
            instrument = "avhrr"
        elif overpass.satellite.name.startswith("aws"):
            scan_angle = 55.25
            instrument = "avhrr"
        else:
            scan_angle = 55.25
        return instrument, scan_angle

    def create_instrument_geometry(self, instrument, scans_nb, scanpoints, scan_step, scan_angle):
        """Create an instrument geometry object."""
        instrument_fun = getattr(geoloc_instrument_definitions,
                                 INSTRUMENT.get(instrument, instrument))

        if instrument.startswith("avhrr"):
            sgeom = instrument_fun(scans_nb, scanpoints, scan_angle=scan_angle, frequency=100)
        elif instrument in ["ascat", ]:
            sgeom = instrument_fun(scans_nb, scanpoints)
        elif instrument in ["amsua", "mhs"]:
            sgeom = instrument_fun(scans_nb, scanpoints)
        elif instrument in ["mwhs2", ]:
            sgeom = instrument_fun(scans_nb, scanpoints)
        elif instrument in ["olci", ]:
            sgeom = instrument_fun(scans_nb, scanpoints)
        elif instrument == "viirs":
            sgeom = instrument_fun(scans_nb, scanpoints, scan_step=scan_step)
        elif instrument in ["mhs", "atms", "mwhs-2"]:
            sgeom = instrument_fun(scans_nb, scanpoints)
        else:
            logger.warning("Instrument not tested: %s", instrument)
            sgeom = instrument_fun(scans_nb)
        return sgeom

    def __init__(self, overpass, scan_step=50, frequency=200):
        """Initialize the boundary.

        Arguments:
            overpass: the overpass to use
            scan_step: how many scans we should skip for a smaller boundary
            frequency: how much to decimate the top and bottom rows of the boundary.
        """
        # compute area covered by pass
        super().__init__()

        self.overpass = overpass
        self.orb = overpass.orb

        # compute sides

        scanlength_seconds = ((overpass.falltime - overpass.risetime).seconds +
                              (overpass.falltime - overpass.risetime).microseconds / 1000000.0)

        logger.debug("Instrument = %s", self.overpass.instrument)
        scan_step, sec_scan_duration, along_scan_reduce_factor = self.get_steps_and_duration(scan_step)

        # From pass length in seconds and the seconds for one scan derive the number of scans in the swath:
        scans_nb = scanlength_seconds / sec_scan_duration * along_scan_reduce_factor
        # Devide by the scan step to a reduced number of scans:
        scans_nb = np.floor(scans_nb / scan_step)
        scans_nb = int(max(scans_nb, 1))

        sides_lons, sides_lats = self.get_instrument_points(self.overpass,
                                                            overpass.risetime,
                                                            scans_nb,
                                                            np.array([0, self.overpass.number_of_fovs - 1]),
                                                            scan_step=scan_step)

        side_shape = sides_lons[::-1, 0].shape[0]
        nmod = 1

        if side_shape != scans_nb:
            nmod = side_shape // scans_nb
            logger.debug("Number of scan lines (%d) does not match number of scans (%d)",
                         side_shape, scans_nb)
            logger.info("Take every %d th element on the sides...", nmod)

        self.left_lons = sides_lons[::-1, 0][::nmod]
        self.left_lats = sides_lats[::-1, 0][::nmod]
        self.right_lons = sides_lons[:, 1][::nmod]
        self.right_lats = sides_lats[:, 1][::nmod]

        # compute bottom
        maxval = self.overpass.number_of_fovs
        rest = maxval % frequency
        mid_range = np.arange(rest / 2, maxval, frequency)
        if mid_range[0] == 0:
            start_idx = 1
        else:
            start_idx = 0

        reduced = np.hstack([0, mid_range[start_idx::], maxval - 1]).astype("int")

        lons, lats = self.get_instrument_points(self.overpass,
                                                overpass.falltime,
                                                1,
                                                reduced)

        self.bottom_lons = lons[0][::-1]
        self.bottom_lats = lats[0][::-1]

        # compute top
        lons, lats = self.get_instrument_points(self.overpass,
                                                overpass.risetime,
                                                1,
                                                reduced)

        self.top_lons = lons[0]
        self.top_lats = lats[0]

        return

    def get_steps_and_duration(self, scan_step):
        """Get the steps and duration for the instrument."""
        if self.overpass.instrument == "viirs":
            sec_scan_duration = 1.779166667
            along_scan_reduce_factor = 1
        elif self.overpass.instrument.startswith("avhrr"):
            sec_scan_duration = 1. / 6.
            along_scan_reduce_factor = 0.1
        elif self.overpass.instrument == "ascat":
            sec_scan_duration = 3.74747474747
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 1
        elif self.overpass.instrument == "amsua":
            sec_scan_duration = 8.
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 1
        elif self.overpass.instrument == "mhs":
            sec_scan_duration = 8./3.
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 1
        elif self.overpass.instrument == "mwhs2":
            sec_scan_duration = 8./3.
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 1
        elif self.overpass.instrument == "olci":
            # 3 minutes of data is 4091 300meter lines:
            sec_scan_duration = 0.04399902224395014
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 100
        elif self.overpass.instrument == "atms":
            sec_scan_duration = 8/3.
            along_scan_reduce_factor = 1
            # Overwrite the scan step
            scan_step = 1

        else:
            # Assume AVHRR!
            logmsg = ("Instrument scan duration not known. Setting it to AVHRR. Instrument: ")
            logger.info(logmsg + "%s", str(self.overpass.instrument))
            sec_scan_duration = 1. / 6.
            along_scan_reduce_factor = 0.1
        return scan_step, sec_scan_duration, along_scan_reduce_factor

    def decimate(self, ratio):
        """Remove points from the boundary."""
        length = len(self.top_lons)
        start = (length % ratio) / 2
        points = np.concatenate(([0], np.arange(start, length, ratio), [length - 1]))

        self.top_lons = self.top_lons[points]
        self.top_lats = self.top_lats[points]
        self.bottom_lons = self.bottom_lons[points]
        self.bottom_lats = self.bottom_lats[points]

        length = len(self.right_lons)
        start = (length % ratio) / 2
        points = np.concatenate(([0], np.arange(start, length, ratio), [length - 1]))

        self.right_lons = self.right_lons[points]
        self.right_lats = self.right_lats[points]
        self.left_lons = self.left_lons[points]
        self.left_lats = self.left_lats[points]

        return

    def contour(self):
        """Get the contour lon/lats."""
        lons = np.concatenate((self.top_lons,
                               self.right_lons[1:-1],
                               self.bottom_lons,
                               self.left_lons[1:-1]))
        lats = np.concatenate((self.top_lats,
                               self.right_lats[1:-1],
                               self.bottom_lats,
                               self.left_lats[1:-1]))
        return lons, lats
