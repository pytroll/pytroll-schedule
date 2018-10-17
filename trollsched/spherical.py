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

import numpy as np
import pyresample.spherical
import logging

logger = logging.getLogger(__name__)


class SCoordinate(object):

    """Spherical coordinates
    """

    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat

    def cross2cart(self, point):
        """Compute the cross product, and convert to cartesian coordinates
        """

        lat1 = self.lat
        lon1 = self.lon
        lat2 = point.lat
        lon2 = point.lon

        ad = np.sin(lat1 - lat2) * np.cos((lon1 - lon2) / 2.0)
        be = np.sin(lat1 + lat2) * np.sin((lon1 - lon2) / 2.0)
        c = np.sin((lon1 + lon2) / 2.0)
        f = np.cos((lon1 + lon2) / 2.0)
        g = np.cos(lat1)
        h = np.cos(lat2)
        i = np.sin(lon2 - lon1)
        res = CCoordinate(np.array([-ad * c + be * f,
                                    ad * f + be * c,
                                    g * h * i]))

        return res

    def to_cart(self):
        """Convert to cartesian.
        """
        return CCoordinate(np.array([np.cos(self.lat) * np.cos(self.lon),
                                     np.cos(self.lat) * np.sin(self.lon),
                                     np.sin(self.lat)]))

    def distance(self, point):
        """Vincenty formula.
        """

        dlambda = self.lon - point.lon
        num = ((np.cos(point.lat) * np.sin(dlambda)) ** 2 +
               (np.cos(self.lat) * np.sin(point.lat) -
                np.sin(self.lat) * np.cos(point.lat) *
                np.cos(dlambda)) ** 2)
        den = (np.sin(self.lat) * np.sin(point.lat) +
               np.cos(self.lat) * np.cos(point.lat) * np.cos(dlambda))

        return np.arctan2(num ** .5, den)

    def hdistance(self, point):
        """Haversine formula
        """

        return 2 * np.arcsin((np.sin((point.lat - self.lat) / 2.0) ** 2.0 +
                              np.cos(point.lat) * np.cos(self.lat) *
                              np.sin((point.lon - self.lon) / 2.0) ** 2.0) ** .5)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return np.allclose((self.lon, self.lat), (other.lon, other.lat))

    def __str__(self):
        return str((np.rad2deg(self.lon), np.rad2deg(self.lat)))

    def __repr__(self):
        return str((np.rad2deg(self.lon), np.rad2deg(self.lat)))

    def __iter__(self):
        return [self.lon, self.lat].__iter__()


class CCoordinate(object):

    """Cartesian coordinates
    """

    def __init__(self, cart):
        self.cart = np.array(cart)

    def norm(self):
        """Euclidean norm of the vector.
        """
        return np.sqrt(np.einsum('...i, ...i', self.cart, self.cart))

    def normalize(self):
        """normalize the vector.
        """

        self.cart /= np.sqrt(np.einsum('...i, ...i', self.cart, self.cart))

        return self

    def cross(self, point):
        """cross product with another vector.
        """
        return CCoordinate(np.cross(self.cart, point.cart))

    def dot(self, point):
        """dot product with another vector.
        """
        return np.inner(self.cart, point.cart)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return np.allclose(self.cart, other.cart)

    def __str__(self):
        return str(self.cart)

    def __repr__(self):
        return str(self.cart)

    def __add__(self, other):
        try:
            return CCoordinate(self.cart + other.cart)
        except AttributeError:
            return CCoordinate(self.cart + np.array(other))

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        try:
            return CCoordinate(self.cart * other.cart)
        except AttributeError:
            return CCoordinate(self.cart * np.array(other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def to_spherical(self):
        return SCoordinate(np.arctan2(self.cart[1], self.cart[0]),
                           np.arcsin(self.cart[2]))


EPSILON = 0.0000001


def modpi(val, mod=np.pi):
    """Puts *val* between -*mod* and *mod*.
    """
    return (val + mod) % (2 * mod) - mod


class Arc(object):

    """An arc of the great circle between two points.
    """
    start = None
    end = None

    def __init__(self, start, end):
        self.start, self.end = start, end

    def __eq__(self, other):
        if(self.start == other.start and self.end == other.end):
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return (str(self.start) + " -> " + str(self.end))

    def __repr__(self):
        return (str(self.start) + " -> " + str(self.end))

    def angle(self, other_arc):
        """Oriented angle between two arcs.
        """
        if self.start == other_arc.start:
            a__ = self.start
            b__ = self.end
            c__ = other_arc.end
        elif self.start == other_arc.end:
            a__ = self.start
            b__ = self.end
            c__ = other_arc.start
        elif self.end == other_arc.end:
            a__ = self.end
            b__ = self.start
            c__ = other_arc.start
        elif self.end == other_arc.start:
            a__ = self.end
            b__ = self.start
            c__ = other_arc.end
        else:
            raise ValueError("No common point in angle computation.")

        ua_ = a__.cross2cart(b__)
        ub_ = a__.cross2cart(c__)

        val = ua_.dot(ub_) / (ua_.norm() * ub_.norm())
        if abs(val - 1) < EPSILON:
            angle = 0
        elif abs(val + 1) < EPSILON:
            angle = np.pi
        else:
            angle = np.arccos(val)

        n__ = ua_.normalize()
        if n__.dot(c__.to_cart()) > 0:
            return -angle
        else:
            return angle

    def intersections(self, other_arc):
        """Gives the two intersections of the greats circles defined by the
       current arc and *other_arc*.
       From http://williams.best.vwh.net/intersect.htm
        """

        if self.end.lon - self.start.lon > np.pi:
            self.end.lon -= 2 * np.pi
        if other_arc.end.lon - other_arc.start.lon > np.pi:
            other_arc.end.lon -= 2 * np.pi
        if self.end.lon - self.start.lon < -np.pi:
            self.end.lon += 2 * np.pi
        if other_arc.end.lon - other_arc.start.lon < -np.pi:
            other_arc.end.lon += 2 * np.pi

        ea_ = self.start.cross2cart(self.end).normalize()
        eb_ = other_arc.start.cross2cart(other_arc.end).normalize()

        cross = ea_.cross(eb_)
        lat = np.arctan2(cross.cart[2],
                         np.sqrt(cross.cart[0] ** 2 + cross.cart[1] ** 2))
        lon = np.arctan2(cross.cart[1], cross.cart[0])

        return (SCoordinate(lon, lat),
                SCoordinate(modpi(lon + np.pi), -lat))

    def intersects(self, other_arc):
        """Says if two arcs defined by the current arc and the *other_arc*
        intersect. An arc is defined as the shortest tracks between two points.
        """

        return bool(self.intersection(other_arc))

    def intersection(self, other_arc):
        """Says where, if two arcs defined by the current arc and the
        *other_arc* intersect. An arc is defined as the shortest tracks between
        two points.
        """
        if self == other_arc:
            return None
        # if (self.end == other_arc.start or
        #     self.end == other_arc.end or
        #     self.start == other_arc.start or
        #     self.start == other_arc.end):
        #     return None

        for i in self.intersections(other_arc):
            a__ = self.start
            b__ = self.end
            c__ = other_arc.start
            d__ = other_arc.end

            ab_ = a__.hdistance(b__)
            cd_ = c__.hdistance(d__)

            if(((i in (a__, b__)) or
                (abs(a__.hdistance(i) + b__.hdistance(i) - ab_) < EPSILON)) and
               ((i in (c__, d__)) or
                    (abs(c__.hdistance(i) + d__.hdistance(i) - cd_) < EPSILON))):
                return i
        return None

    def get_next_intersection(self, arcs, known_inter=None):
        """Get the next intersection between the current arc and *arcs*
        """
        res = []
        for arc in arcs:
            inter = self.intersection(arc)
            if (inter is not None and
                    inter != arc.end and
                    inter != self.end):
                res.append((inter, arc))

        def dist(args):
            """distance key.
            """
            return self.start.distance(args[0])

        take_next = False
        for inter, arc in sorted(res, key=dist):
            if known_inter is not None:
                if known_inter == inter:
                    take_next = True
                elif take_next:
                    return inter, arc
            else:
                return inter, arc

        return None, None


class SphPolygon(pyresample.spherical.SphPolygon):

    def draw(self, mapper, options, **more_options):
        lons = np.rad2deg(self.lon.take(np.arange(len(self.lon) + 1),
                                        mode="wrap"))
        lats = np.rad2deg(self.lat.take(np.arange(len(self.lat) + 1),
                                        mode="wrap"))
        rx, ry = mapper(lons, lats)
        mapper.plot(rx, ry, options, **more_options)


def get_twilight_poly(utctime):
    """Return a polygon enclosing the sunlit part of the globe at *utctime*.
    """
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
