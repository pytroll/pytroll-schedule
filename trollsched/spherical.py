#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014, 2015, 2018 Martin Raspaud

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

"""Some generalized spherical functions.

base type is a numpy array of size (n, 2) (2 for lon and lats)

"""

import logging

import numpy as np
import pyresample.spherical

logger = logging.getLogger(__name__)

EPSILON = 0.0000001


def modpi(val, mod=np.pi):
    """Put *val* between -*mod* and *mod*."""
    return (val + mod) % (2 * mod) - mod


class SphPolygon(pyresample.spherical.SphPolygon):
    """A spherical polygon with drawing capabilities."""

    def draw(self, mapper, options, **more_options):
        """Draw the polygon."""
        lons = np.rad2deg(self.lon.take(np.arange(len(self.lon) + 1),
                                        mode="wrap"))
        lats = np.rad2deg(self.lat.take(np.arange(len(self.lat) + 1),
                                        mode="wrap"))
        rx, ry = mapper(lons, lats)
        mapper.plot(rx, ry, options, **more_options)


def get_twilight_poly(utctime):
    """Return a polygon enclosing the sunlit part of the globe at *utctime*."""
    from pyorbital import astronomy
    ra, dec = astronomy.sun_ra_dec(utctime)
    lon = modpi(ra - astronomy.gmst(utctime))
    lat = dec

    vertices = np.zeros((4, 2))

    vertices[0, :] = modpi(lon - np.pi / 2), 0
    if lat <= 0:
        vertices[1, :] = lon, np.pi / 2 + lat
        vertices[3, :] = modpi(lon + np.pi), -(np.pi / 2 + lat)
    else:
        vertices[1, :] = modpi(lon + np.pi), np.pi / 2 - lat
        vertices[3, :] = lon, -(np.pi / 2 - lat)

    vertices[2, :] = modpi(lon + np.pi / 2), 0

    return SphPolygon(vertices)
