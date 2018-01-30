#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Alexander Maul
#
# Author(s):
#
#   Alexander Maul <alexander.maul@dwd.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Combine several graphs.
"""
import os
import yaml
import logging
from collections import Mapping
from ConfigParser import ConfigParser

import schedule
import satpass
from pyresample import utils as resample_utils

logger = logging.getLogger("trollsched")

def read_yaml_file(file_name):
    """Read one or more files in to a single dict object."""
    if isinstance(file_name, str):
        file_name = [file_name]
    conf_dict = {}
    for file_obj in file_name:
        with open(file_obj) as fp:
            tmp_dict = yaml.load(fp)
        conf_dict = recursive_dict_update(conf_dict, tmp_dict)
    return conf_dict

def recursive_dict_update(d, u):
    """Recursive dictionary update.
    Copied from:
        http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            r = recursive_dict_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

def read_config(filename):
    try:
        return read_config_yaml(filename)
    except yaml.parser.ParserError as e:
        return read_config_cfg(filename)

def read_config_cfg(filename):
    """Read the config file *filename* and replace the values in global
    variables.
    """
    station_list = []
    cfg = ConfigParser()
    cfg.read(filename)

    stations = cfg.get("default", "station").split(",")
    forward = cfg.getint("default", "forward")
    start = cfg.getfloat("default", "start")

    pattern = {}
    for k, v in cfg.items("pattern"):
        pattern[k] = v

    for station in stations:
        station_name = cfg.get(station, "name")
        station_lon = cfg.getfloat(station, "longitude")
        station_lat = cfg.getfloat(station, "latitude")
        station_alt = cfg.getfloat(station, "altitude")

        area = resample_utils.parse_area_file(cfg.get(station, "area_file"),
                                                        cfg.get(station, "area"))[0]

        satellites = cfg.get(station, "satellites").split(",")

        sat_scores = {}
        for sat in satellites:
            sat_scores[sat] = (cfg.getfloat(sat, "night"),
                               cfg.getfloat(sat, "day"))

        station_list.append(((station_lon, station_lat, station_alt),
                station_name, area, sat_scores))

    return (station_list, forward, start, pattern)

def read_config_yaml(filename):
    """Read the yaml file *filename* and replace the values in global
    variables.
    """
    station_list = []
    cfg = read_yaml_file(filename)

    if cfg["default"].get("center_id"):
        schedule.CENTER_ID = cfg["default"].get("center_id")
    if cfg["default"].get("dump_url"):
        satpass.HOST = cfg["default"].get("dump_url")
    if cfg["default"].get("min_pass"):
        satpass.MIN_PASS = cfg["default"].get("min_pass")

    stations = cfg["default"]["station"]
    forward = cfg["default"]["forward"]
    start = cfg["default"]["start"]

    pattern = {}
    for k, v in cfg["pattern"].items():
        pattern[k] = v

    for station in stations:
        station_name = cfg["stations"][station]["name"]
        station_lon = cfg["stations"][station]["longitude"]
        station_lat = cfg["stations"][station]["latitude"]
        station_alt = cfg["stations"][station]["altitude"]

        area = resample_utils.parse_area_file(cfg["stations"][station]["area_file"],
                                              cfg["stations"][station]["area"])[0]

        satellites = cfg["stations"][station]["satellites"]

        sat_scores = {}
        for sat in satellites:
            sat_scores[sat] = (cfg["satellites"][sat]["night"],
                                             cfg["satellites"][sat]["day"])

        station_list.append(((station_lon, station_lat, station_alt),
                station_name, area, sat_scores))

    return (station_list, forward, start, pattern)

