#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll developers

# Author(s):

#   Adam Dybbroe <Firstname.Lastname at smhi.se>

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

"""Common utilities for satellite pass and schedule calculations."""


from pathlib import Path
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
import pyproj
from shapely.ops import transform

class SatScore:
    """docstring for SatScore."""

    def __init__(self, day, night):
        """Initialize the score."""
        self.day = day
        self.night = night


class Satellite:
    """docstring for Satellite."""

    def __init__(self, name, day, night,
                 schedule_name=None, international_designator=None):
        """Initialize the satellite."""
        self.name = name
        self.international_designator = international_designator
        self.score = SatScore(day, night)
        self.schedule_name = schedule_name or name


def get_polygon_from_contour(contour_poly):
    """From a pytroll-schedule contour-poly return a shapely Polygon."""
    geodata = np.vstack((contour_poly.lon, contour_poly.lat)).T
    return Polygon(np.rad2deg(geodata))


def create_shapefile_from_pass(sat_pass, outpath):
    """From a satellite overpass (instrument scanning outline) create a shapefile and save."""
    sat_poly = get_polygon_from_contour(sat_pass.boundary.contour_poly)

    wgs84 = pyproj.CRS('EPSG:4326')  # WGS 84
    mycrs = {'init': 'epsg:4326'}

    project = pyproj.Transformer.from_crs(wgs84, mycrs, always_xy=True).transform
    new_shapes = transform(project, sat_poly)

    satname = sat_pass.satellite.name.replace(' ', '-')
    outname = f'{sat_pass.instrument}_{satname}_{sat_pass.risetime:%Y%m%d%H%M}_{sat_pass.falltime:%Y%m%d%H%M}_outline.shp'
    output_filepath = Path(outpath) / outname

    gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                     crs=mycrs,
                     geometry=[new_shapes]).to_file(output_filepath)
