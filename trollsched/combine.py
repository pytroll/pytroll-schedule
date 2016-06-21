#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

# Author(s):

#   Alexander Maul <alexander.maul@dwd.de>

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

"""Combine several graphs.
"""
import logging
import logging.handlers
import urlparse
import os
from datetime import datetime, timedelta
from pprint import pformat
import numpy as np
from pyresample import utils
from pyorbital import astronomy
from trollsched.spherical import get_twilight_poly
from trollsched.graph import Graph
from trollsched.satpass import get_next_passes, SimplePass
from trollsched.boundary import AreaDefBoundary
from trollsched.schedule import get_passes_from_xml_file

def get_best_sched(overpasses, area_of_interest, scores, delay, avoid_list=None):
    """Get the best schedule based on *area_of_interest*.
    """
    avoid_list = avoid_list or []
    print avoid_list
    raw_input()
    passes = sorted(overpasses, key=lambda x: x.risetime)
    grs = conflicting_passes(passes, delay)
    logger.debug("conflicting %s", str(grs))
    ncgrs = [get_non_conflicting_groups(gr, delay) for gr in grs]
    logger.debug("non conflicting %s", str(ncgrs))
    n_vertices = len(passes)

    graph = Graph(n_vertices=n_vertices + 2)

    def add_arc(graph, p1, p2, hook=None):
        logger.debug("Adding arc between " + str(p1) +
                     " and " + str(p2) + "...")
        if p1 in avoid_list or p2 in avoid_list:
            w = 0
            logger.debug("...0 because in the avoid_list!")
        else:
            w = combine(p1, p2, area_of_interest, scores)
        logger.debug("...with weight " + str(w))

#         with open("/tmp/schedule.gv", "a") as fp_:
#             fp_.write('        "' + str(p1) + '" -> "' + str(p2) +
#                       '" [ label = "' + str(w) + '" ];\n')

        graph.add_arc(passes.index(p1) + 1,
                      passes.index(p2) + 1, w)
        if hook is not None:
            hook()

    prev = set()
    for ncgr in ncgrs:
        for pr in prev:
            foll = set(gr[0] for gr in ncgr)
            for f in foll:
                add_arc(graph, pr, f)

        prev = set(sorted(gr, key=lambda x: x.falltime)[-1] for gr in ncgr)
        for gr in ncgr:
            if len(gr) > 1:
                for p1, p2 in zip(gr[:-1], gr[1:]):
                    add_arc(graph, p1, p2)

    for pr in prev:
        graph.add_arc(passes.index(pr) + 1, n_vertices + 1)
    for first in ncgrs[0][0]:
        graph.add_arc(0, passes.index(first) + 1)

    dist, path = graph.dag_longest_path(0, n_vertices + 1)

    del dist
    return [passes[idx - 1] for idx in path[1:-1]], (graph, passes)

def add_graphs(graph_set, passes_set):
    """add all graphs to one combined graph.
    """
    pass


if __name__ == '__main__':
    try:
        from trollsched.schedule import read_config
        import argparse
        #global logger
        logger = logging.getLogger("trollsched")
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--config", help="configuration file to use", default=None)
        parser.add_argument("-g", "--graph", help="save graph info to this directory", default=None)
        parser.add_argument("-r", "--report", help="xml report file", default=None)
        opts = parser.parse_args()

        if opts.config is None:
            parser.error("Configuration file required.")
        if opts.graph is None:
            parser.error("Graph directory required.")
        if opts.report is None:
            parser.error("Report file required.")
        
        # [coords, station, area, scores], forward, start
        station_list, forward, start = read_config(opts.config)
        graph = {}
        allpasses = {}
        for coords, station, area, scores in station_list:
            graph[station] = Graph()
            graph[station].load(os.path.join(opts.graph, "graph." + station + ".npz"))
            print graph[station].adj_matrix
            print graph[station].weight_matrix
            allpasses[station] = get_passes_from_xml_file(os.path.join(opts.report, "acquisition-schedule-report." + station + ".xml"))
            print allpasses[station]
            
            for v in graph[station].neighbours(1):
                print v, " : ", allpasses[station][v].risetime, "->", graph[station].weight(1, v)

    except:
        logger.exception("Something wrong happened!")
        raise


