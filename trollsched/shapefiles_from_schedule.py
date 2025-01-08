#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll community

# Author(s):

#   Adam Dybbroe <Firstname.Lastname @ smhi.se>

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

"""Create shapefiles of satellite instrument swath outlines from an XML pass reception plan."""


import argparse
import logging
from datetime import datetime
from pathlib import Path

import defusedxml.ElementTree as ET
import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
from shapely.geometry import Polygon
from shapely.ops import transform

from trollsched import INSTRUMENT, SATELLITE_NAMES
from trollsched.logger import setup_logging
from trollsched.satpass import create_pass

logger = logging.getLogger(__name__)


def get_shapely_polygon_from_lonlat(lons, lats):
    """From a arrays of lons and lats return a shapely Polygon."""
    geodata = np.vstack((lons, lats)).T
    return Polygon(np.rad2deg(geodata))


def create_shapefile_filename(satellite_pass_obj):
    """From a trollsched.satpass instance create the shapefile filename."""
    satname = satellite_pass_obj.satellite.name.replace(" ", "-")
    prefix = f"{satellite_pass_obj.instrument}_{satname}"
    return f"{prefix}_{satellite_pass_obj.risetime:%Y%m%d%H%M}_{satellite_pass_obj.falltime:%Y%m%d%H%M}_outline.shp"


def create_shapefile_from_pass(sat_pass, output_filepath):
    """From a satellite overpass (instrument scanning outline) create a shapefile and save."""
    sat_poly = get_shapely_polygon_from_lonlat(sat_pass.boundary.contour_poly.lon,
                                               sat_pass.boundary.contour_poly.lat)

    wgs84 = pyproj.CRS("EPSG:4326")  # WGS 84
    a_crs = {"init": "epsg:4326"}

    project = pyproj.Transformer.from_crs(wgs84, a_crs, always_xy=True).transform
    new_shapes = transform(project, sat_poly)

    gpd.GeoDataFrame(pd.DataFrame(["p1"], columns=["geom"]), crs=a_crs,
                     geometry=[new_shapes]).to_file(output_filepath)


def shapefiles_from_schedule_xml_requests(filename, satellites, tle_file, output_dir):
    """From XML file with planned satellite passes create shapefiles for each Pass outline.

    Open the XML file with the scheduled satellite pass plan, and create
    shapefiles for each Pass outline assuming the instrument as provided in the
    INSTRUMENTS dictionary.
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    for child in root:
        if child.tag == "pass":
            platform_name = SATELLITE_NAMES.get(child.attrib["satellite"], child.attrib["satellite"])
            instrument = INSTRUMENT.get(platform_name)
            if not instrument:
                logger.debug("Instrument not defined for platform %s", platform_name)
                continue

            if platform_name not in satellites:
                logger.debug("Satellite platform name not requested: %s", platform_name)
                continue

            overpass = create_pass(platform_name, instrument,
                                   datetime.strptime(child.attrib["start-time"],
                                                     "%Y-%m-%d-%H:%M:%S"),
                                   datetime.strptime(child.attrib["end-time"],
                                                     "%Y-%m-%d-%H:%M:%S"),
                                   tle_filename=tle_file)

            output_filepath = Path(output_dir) / create_shapefile_filename(overpass)
            create_shapefile_from_pass(overpass, output_filepath)


def parse_args(args):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log-config",
                        dest="log_config",
                        default=None,
                        required=False,
                        help="Log config file to use instead of the standard logging.")
    parser.add_argument("-s", "--satellites",
                        dest="list_of_satellites",
                        nargs="*",
                        default=[],
                        help="Complete file path to the TLE file to use.",
                        required=True)
    parser.add_argument("-t", "--tle_filename",
                        dest="tle_filepath",
                        help="Complete file path to the TLE file to use.",
                        required=True)
    parser.add_argument("-x", "--xml_filename",
                        dest="xml_filepath",
                        help="Complete path to the XML satellite schedule file.",
                        required=True)
    parser.add_argument("-o", "--output_dir",
                        dest="output_dir",
                        help="Complete path to the XML satellite schedule file.",
                        default="./")
    parser.add_argument("-v", "--verbose",
                        dest="verbosity",
                        action="count",
                        default=0,
                        help="Verbosity (between 1 and 2 occurrences with more leading to more "
                        "verbose logging). WARN=0, INFO=1, "
                        "DEBUG=2. This is overridden by the log config file if specified.")

    return parser.parse_args(args)


def run(args=None):
    """Create shapefiles from an XML schedule of satellite passes."""
    opts = parse_args(args)

    setup_logging(opts)
    filename = opts.xml_filepath
    satellites = opts.list_of_satellites
    tle_filename = opts.tle_filepath
    outdir = opts.output_dir

    shapefiles_from_schedule_xml_requests(filename,
                                          satellites,
                                          tle_filename,
                                          outdir)
