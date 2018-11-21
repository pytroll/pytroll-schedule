#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, 2018 Alexander Maul
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

"""Utility functions and config reading for the pytroll-scheduler
"""

import yaml
import logging
from collections import Mapping
from six.moves.configparser import ConfigParser

try:
    from trollsched import schedule
except ImportError:
    import schedule


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
    cfg = ConfigParser()
    cfg.read(filename)

    def read_cfg_opts(section):
        """Read the option:value pairs in one section,
        converting value to int/float if applicable.
        """
        kv_dict = {}
        for k, v in cfg.items(section):
            try:
                kv_dict[k] = int(v)
            except:
                try:
                    kv_dict[k] = float(v)
                except:
                    kv_dict[k] = v
        return kv_dict

    default_params = read_cfg_opts("default")
    pattern = {}
    for k, v in cfg.items("pattern"):
        pattern[k] = v
    station_list = []
    for station_id in default_params["station"].split(","):
        station_params = read_cfg_opts(station_id)
        satellites = cfg.get(station_id, "satellites").split(",")
        sat_list = []
        for sat_name in satellites:
            sat_list.append(schedule.Satellite(sat_name,
                                               **read_cfg_opts(sat_name)
                                               ))
        new_station = schedule.Station(station_id, **station_params)
        new_station.satellites = sat_list
        station_list.append(new_station)
    scheduler = schedule.Scheduler(stations=station_list,
                                   min_pass=default_params.get("min_pass", 4),
                                   forward=default_params.get("forward"),
                                   start=default_params.get("start"),
                                   dump_url=default_params.get("dump_url", None),
                                   patterns=pattern,
                                   center_id=default_params.get("center_id", "unknown"))
    return scheduler


def read_config_yaml(filename):
    """Read the yaml file *filename* and create a scheduler."""

    cfg = read_yaml_file(filename)
    satellites = {sat_name: schedule.Satellite(sat_name, **sat_params)
                  for (sat_name, sat_params) in cfg["satellites"].items()}

    stations = {}
    for station_id, station in cfg["stations"].items():
        if isinstance(station['satellites'], dict):
            sat_list = []
            for (sat_name, sat_params) in station["satellites"].items():
                if sat_params is None:
                    sat_list.append(satellites[sat_name])
                else:
                    sat_list.append(schedule.Satellite(sat_name, **sat_params))
        else:
            sat_list = [satellites[sat_name] for sat_name in station['satellites']]
        new_station = schedule.Station(station_id, **station)
        new_station.satellites = sat_list
        stations[station_id] = new_station

    pattern = {}
    for k, v in cfg["pattern"].items():
        pattern[k] = v

    sched_params = cfg['default']
    scheduler = schedule.Scheduler(stations=[stations[st_id]
                                             for st_id in sched_params['station']],
                                   min_pass=sched_params.get('min_pass', 4),
                                   forward=sched_params['forward'],
                                   start=sched_params['start'],
                                   dump_url=sched_params.get('dump_url'),
                                   patterns=pattern,
                                   center_id=sched_params.get('center_id', 'unknown'))

    return scheduler
