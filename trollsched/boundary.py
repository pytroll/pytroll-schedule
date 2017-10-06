#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015, 2017 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""The Boundary classes.
"""


import logging
import logging.handlers

import numpy as np

from pyorbital import geoloc, geoloc_instrument_definitions
from trollsched.spherical import SphPolygon

logger = logging.getLogger(__name__)


class Boundary(object):

    """Boundary objects.
    """

    def __init__(self, lons=None, lats=None, frequency=1):
        self._contour_poly = None
        if lons is not None:
            self.lons = lons[::frequency]
        if lats is not None:
            self.lats = lats[::frequency]

    def contour(self):
        return self.lons, self.lats

    @property
    def contour_poly(self):
        """Get the Spherical polygon corresponding to the Boundary
        """
        if self._contour_poly is None:
            self._contour_poly = SphPolygon(
                np.deg2rad(np.vstack(self.contour()).T))
        return self._contour_poly

    def draw(self, mapper, options):
        """Draw the current boundary on the *mapper*
        """
        self.contour_poly.draw(mapper, options)


class AreaBoundary(Boundary):

    """Area boundary objects.
    """

    def __init__(self, *sides):
        Boundary.__init__(self)
        self.sides_lons, self.sides_lats = zip(*sides)
        self.sides_lons = list(self.sides_lons)
        self.sides_lats = list(self.sides_lats)

    def decimate(self, ratio):
        """Remove some points in the boundaries, but never the corners.
        """
        for i in range(len(self.sides_lons)):
            l = len(self.sides_lons[i])
            start = int((l % ratio) / 2)
            points = np.concatenate(([0], np.arange(start, l, ratio), [l - 1]))
            if points[1] == 0:
                points = points[1:]
            if points[-2] == (l - 1):
                points = points[:-1]
            self.sides_lons[i] = self.sides_lons[i][points]
            self.sides_lats[i] = self.sides_lats[i][points]

    def contour(self):
        """Get the (lons, lats) tuple of the boundary object.
        """
        lons = np.concatenate([lns[:-1] for lns in self.sides_lons])
        lats = np.concatenate([lts[:-1] for lts in self.sides_lats])

        return lons, lats


class AreaDefBoundary(AreaBoundary):

    """Boundaries for area definitions (pyresample)
    """

    def __init__(self, area, frequency=1):
        lons, lats = area.get_boundary_lonlats()
        AreaBoundary.__init__(self,
                              (lons.side1, lats.side1),
                              (lons.side2, lats.side2),
                              (lons.side3, lats.side3),
                              (lons.side4, lats.side4))

        if frequency != 1:
            self.decimate(frequency)


class SwathBoundary(Boundary):

    """Boundaries for satellite overpasses.
    """

    def get_instrument_points(self, overpass, utctime,
                              scans_nb, scanpoints, frequency=1):
        """Get the boundary points for a given overpass.
        """
        instrument = overpass.instrument
        # cheating at the moment.
        scan_angle = 55.37
        if instrument == "modis":
            scan_angle = 55.0
        elif instrument == "viirs":
            scan_angle = 55.84
        elif instrument == "iasi":
            scan_angle = 48.3
        elif overpass.satellite == "noaa 16":
            scan_angle = 55.25
        instrument = "avhrr"
        instrument_fun = getattr(geoloc_instrument_definitions, instrument)
        sgeom = instrument_fun(scans_nb, scanpoints,
                               scan_angle=scan_angle, frequency=frequency)
        times = sgeom.times(utctime)

        pixel_pos = geoloc.compute_pixels((self.orb.tle._line1,
                                           self.orb.tle._line2),
                                          sgeom, times)
        lons, lats, alts = geoloc.get_lonlatalt(pixel_pos, times)

        del alts
        return (lons.reshape(-1, len(scanpoints)),
                lats.reshape(-1, len(scanpoints)))

    def __init__(self, overpass, frequency=100.0):
        # compute area covered by pass

        Boundary.__init__(self)

        self.overpass = overpass
        self.orb = overpass.orb

        # compute sides

        scans_nb = np.ceil(((overpass.falltime - overpass.risetime).seconds +
                            (overpass.falltime -
                             overpass.risetime).microseconds
                            / 1000000.0) / frequency)

        scans_nb = int(max(scans_nb, 1))

        sides_lons, sides_lats = self.get_instrument_points(self.overpass,
                                                            overpass.risetime,
                                                            scans_nb,
                                                            np.array(
                                                                [0, 2047]),
                                                            frequency=frequency)

        self.left_lons = sides_lons[::-1, 0]
        self.left_lats = sides_lats[::-1, 0]
        self.right_lons = sides_lons[:, 1]
        self.right_lats = sides_lats[:, 1]

        # compute bottom

        # avhrr
        maxval = 2048
        rest = maxval % frequency
        reduced = np.hstack(
            [0, np.arange(rest / 2, maxval, frequency), maxval - 1])

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

    def decimate(self, ratio):
        l = len(self.top_lons)
        start = (l % ratio) / 2
        points = np.concatenate(([0], np.arange(start, l, ratio), [l - 1]))

        self.top_lons = self.top_lons[points]
        self.top_lats = self.top_lats[points]
        self.bottom_lons = self.bottom_lons[points]
        self.bottom_lats = self.bottom_lats[points]

        l = len(self.right_lons)
        start = (l % ratio) / 2
        points = np.concatenate(([0], np.arange(start, l, ratio), [l - 1]))

        self.right_lons = self.right_lons[points]
        self.right_lats = self.right_lats[points]
        self.left_lons = self.left_lons[points]
        self.left_lats = self.left_lats[points]

    def contour(self):
        lons = np.concatenate((self.top_lons,
                               self.right_lons[1:-1],
                               self.bottom_lons,
                               self.left_lons[1:-1]))
        lats = np.concatenate((self.top_lats,
                               self.right_lats[1:-1],
                               self.bottom_lats,
                               self.left_lats[1:-1]))
        return lons, lats
