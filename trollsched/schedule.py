#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014, 2015, 2016 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
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

"""Scheduling
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
from trollsched.combine import get_combined_sched

from ConfigParser import ConfigParser

logger = logging.getLogger(__name__)

# shortest allowed pass in minutes
MIN_PASS = 4

# name/id for centre/org creating schedules
# CENTER_ID = "SMHI"
CENTER_ID = "DWD-OF"

def conflicting_passes(allpasses, delay=timedelta(seconds=0)):
    """Get the passes in groups of conflicting passes.
    """

    passes = sorted(allpasses, key=lambda x: x.risetime)

    overpass = passes[0]
    last_time = overpass.falltime
    group = [overpass]
    groups = []
    for overpass in passes[1:]:
        if overpass.risetime - delay < last_time:
            group.append(overpass)
            if last_time < overpass.falltime:
                last_time = overpass.falltime

        else:
            groups.append(group)
            group = [overpass]
            last_time = overpass.falltime

    groups.append(group)
    return groups


def get_non_conflicting_groups(passes, delay=timedelta(seconds=0)):
    """Get the different non-conflicting solutions in a group of conflicting
    passes.
    """
    # Uses graphs and maximal clique finding with the Bron-Kerbosch algorithm.

    order = len(passes)

    if order == 1:
        return [passes]

    graph = Graph(order)

    for i, overpass in enumerate(sorted(passes, key=lambda x: x.risetime)):
        for j in range(i + 1, order):
            if not overpass.overlaps(passes[j], delay):
                graph.add_edge(i, j)

    groups = []
    for res in graph.bron_kerbosch(set(), set(graph.vertices), set()):
        grp = []
        for vertex in res:
            grp.append(passes[vertex])
        groups.append(sorted(grp))

    return groups


def fermia(t):
    a = 0.25
    b = a / 4
    k = b * np.log(1 / 3.0) + a
    sh = k - 0.25
    return 0.5 / (np.exp(((t + sh) - a) / b) + 1) + 0.5


def fermib(t):
    a = 0.25
    b = a / 4
    return 1 / (np.exp((t - a) / b) + 1)

combination = {}


def combine(p1, p2, area_of_interest, scores):
    """Combine passes together.
    """

    try:
        return combination[p1, p2]
    except KeyError:
        pass

    area = area_of_interest.poly.area()

    def pscore(poly, coeff=1):
        if poly is None:
            return 0
        else:
            return poly.area() * coeff

    twi1 = get_twilight_poly(p1.uptime)
    twi2 = get_twilight_poly(p2.uptime)

    ip1, sip1 = p1.score.get(area_of_interest, (None, None))
    if sip1 is None:
        ip1 = p1.boundary.contour_poly.intersection(area_of_interest.poly)
        # FIXME: ip1 or ip2 could be None if the pass is entirely inside the
        # area (or vice versa)
        if ip1 is None:
            return 0

        ip1d = ip1.intersection(twi1)
        if ip1d is None:
            lon, lat = np.rad2deg(ip1.vertices[0, :])
            theta = astronomy.cos_zen(p1.uptime,
                                      lon, lat)
            if np.sign(theta) > 0:
                ip1d = ip1
                ip1n = None
            else:
                ip1n = ip1
        else:
            twi1.invert()
            ip1n = ip1.intersection(twi1)
            twi1.invert()

        ns1 = pscore(ip1n, scores[p1.satellite][0] / area)
        ds1 = pscore(ip1d, scores[p1.satellite][1] / area)
        sip1 = ns1 + ds1
        p1.score[area_of_interest] = (ip1, sip1)

    ip2, sip2 = p2.score.get(area_of_interest, (None, None))
    if sip2 is None:
        ip2 = p2.boundary.contour_poly.intersection(area_of_interest.poly)
        if ip2 is None:
            return 0

        ip2d = ip2.intersection(twi2)
        if ip2d is None:
            lon, lat = np.rad2deg(ip2.vertices[0, :])
            theta = astronomy.cos_zen(p2.uptime,
                                      lon, lat)
            if np.sign(theta) > 0:
                ip2d = ip2
                ip2n = None
            else:
                ip2n = ip2
        else:
            twi2.invert()
            ip2n = ip2.intersection(twi2)
            twi2.invert()

        ns2 = pscore(ip2n, scores[p2.satellite][0] / area)
        ds2 = pscore(ip2d, scores[p2.satellite][1] / area)
        sip2 = ns2 + ds2
        p2.score[area_of_interest] = (ip2, sip2)

    ip1p2 = ip1.intersection(ip2)

    if ip1p2 is None:
        sip1p2 = 0
    else:
        ip1p2da = ip1p2.intersection(twi1)
        twi1.invert()
        ip1p2na = ip1p2.intersection(twi1)
        twi1.invert()

        ip1p2db = ip1p2.intersection(twi2)
        twi2.invert()
        ip1p2nb = ip1p2.intersection(twi2)
        twi2.invert()

        ns12a = pscore(ip1p2na, scores[p1.satellite][0] / area)
        ds12a = pscore(ip1p2da, scores[p1.satellite][1] / area)
        ns12b = pscore(ip1p2nb, scores[p2.satellite][0] / area)
        ds12b = pscore(ip1p2db, scores[p2.satellite][1] / area)

        sip1p2a = ns12a + ds12a
        sip1p2b = ns12b + ds12b
        sip1p2 = (sip1p2a + sip1p2b) / 2.0

    if p2 > p1:
        tdiff = (p2.uptime - p1.uptime).seconds / 3600.
    else:
        tdiff = (p1.uptime - p2.uptime).seconds / 3600.

    res = fermia(tdiff) * (sip1 + sip2) - fermib(tdiff) * sip1p2
    combination[p1, p2] = res

    return res


def get_best_sched(overpasses, area_of_interest, scores, delay, avoid_list=None):
    """Get the best schedule based on *area_of_interest*.
    """
    avoid_list = avoid_list or []
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


def argmax(iterable):
    return max((x, i) for i, x in enumerate(iterable))[1]


def get_max(groups, fun):
    """Get the best group of *groups* using the score function *fun*
    """
    scores = []
    for grp in groups:
        scores.append(sum([fun(p) for p in grp]))
    return groups[argmax(scores)]


def generate_sch_file(output_file, overpasses, coords):

    with open(output_file, "w") as out:
        # create epochs
        out.write("#Orbital elements\n#\n#SCName           Epochtime\n#\n")
        satellites = set()

        for overpass in overpasses:
            epoch = "!{0:<16} {1}".format(overpass.satellite.upper(),
                                          overpass.orb.tle.epoch.strftime("%Y%m%d %H%M%S"))
            satellites |= set([epoch])
        sats = "\n".join(satellites) + "\n"
        out.write(sats)
        out.write("#\n#\n#Pass List\n#\n")

        out.write(
            "#SCName          RevNum Risetime        Falltime        Elev Dura ANL   Rec Dir Man Ovl OvlSCName        OvlRev OvlRisetime     OrigRisetime    OrigFalltime    OrigDuration\n#\n")

        for overpass in sorted(overpasses):
            out.write(overpass.print_vcs(coords) + "\n")


def generate_xml_requests(sched, start, end, station_name, report_mode=False):
    """Create xml requests.
    """
    import xml.etree.ElementTree as ET

    sats = {"noaa 15": "noaa15",
            "noaa 16": "noaa16",
            "noaa 18": "noaa18",
            "noaa 19": "noaa19",
            "metop-a": "metop-a",
            "metop-b": "metop-b",
            "terra": "terra",
            "aqua": "aqua",
            "suomi npp": "npp",
            }

    reqtime = datetime.utcnow()
    eum_format = "%Y-%m-%d-%H:%M:%S"

    root = ET.Element("acquisition-schedule")
    props = ET.SubElement(root, "properties")
    proj = ET.SubElement(props, "project")
    proj.text = "Pytroll"
    typep = ET.SubElement(props, "type")
    if report_mode:
        typep.text = "report"
    else:
        typep.text = "request"
    station = ET.SubElement(props, "station")
    station.text = station_name
    file_start = ET.SubElement(props, "file-start")
    file_start.text = start.strftime(eum_format)
    file_end = ET.SubElement(props, "file-end")
    file_end.text = end.strftime(eum_format)
    reqby = ET.SubElement(props, "requested-by")
    reqby.text = CENTER_ID
    reqon = ET.SubElement(props, "requested-on")
    reqon.text = reqtime.strftime(eum_format)
    for overpass in sorted(sched):
        if (overpass.rec or report_mode) and overpass.risetime > start:
            ovpass = ET.SubElement(root, "pass")
            ovpass.set("satellite", sats.get(overpass.satellite,
                                             overpass.satellite))
            ovpass.set("start-time", overpass.risetime.strftime(eum_format))
            ovpass.set("end-time", overpass.falltime.strftime(eum_format))
            if report_mode:
                if overpass.fig is not None:
                    ovpass.set("img", overpass.fig)
                ovpass.set("rec", str(overpass.rec))

    return root, reqtime


def generate_xml_file(sched, start, end, xml_file, station, report_mode=False):
    """Create an xml request file.
    """
    import xml.etree.ElementTree as ET
    tree, reqtime = generate_xml_requests(sched,
                                          start, end,
                                          station, report_mode)
    filename = xml_file
    tmp_filename = xml_file + reqtime.strftime("%Y-%m-%d-%H-%M-%S") + ".tmp"
    with open(tmp_filename, "w") as fp_:
        if report_mode:
            fp_.write("<?xml version='1.0' encoding='utf-8'?>"
                      "<?xml-stylesheet type='text/xsl' href='reqreader.xsl'?>")
        fp_.write(ET.tostring(tree))
    os.rename(tmp_filename, filename)
    return filename


def parse_datetime(strtime):
    """Parse the time string *strtime*
    """
    return datetime.strptime(strtime, "%Y%m%d%H%M%S")


def read_config(filename):
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

        area = utils.parse_area_file(cfg.get(station, "area_file"), 
                                                        cfg.get(station, "area"))[0]

        satellites = cfg.get(station, "satellites").split(",")

        sat_scores = {}
        for sat in satellites:
            sat_scores[sat] = (cfg.getfloat(sat, "night"),
                               cfg.getfloat(sat, "day"))

        station_list.append(((station_lon, station_lat, station_alt),
                station_name, area, sat_scores))

    return (station_list, forward, start, pattern)


def save_passes(allpasses, poly, output_dir):
    for passage in allpasses:
        passage.save_fig(poly, directory=output_dir)

def get_passes_from_xml_file(filename):
    """Read passes from aquisition xml file."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(filename)
    root = tree.getroot()
    pass_list = []
    for overpass in root.iter('pass'):
        start_time = datetime.strptime(overpass.attrib['start-time'], '%Y-%m-%d-%H:%M:%S')
        end_time = datetime.strptime(overpass.attrib['end-time'], '%Y-%m-%d-%H:%M:%S')
        pass_list.append(SimplePass(overpass.attrib['satellite'], start_time, end_time))
    return pass_list


def build_filename(pattern_name, pattern_dict, kwargs):
    """Build absolute path from pattern dictionary."""
    for k in pattern_dict.keys():
        for v in pattern_dict.values():
            if "{" + k + "}" in v:
                kwargs[k] = pattern_dict[k].format(**kwargs)
    return pattern_dict[pattern_name].format(**kwargs)


def send_file(url, file):
    pathname, filename = os.path.split(file)
    del pathname
    if url.scheme in ["file", ""]:
        pass
    elif url.scheme == "ftp":
        import ftplib
        session = ftplib.FTP(url.hostname, url.username, url.password)
        with open(file, "rb") as xfile:
            session.storbinary('STOR ' + str(filename), xfile)
        session.quit()
    else:
        logger.error("Cannot save to %s, but file is there:", str(url.scheme), str(file))


def single_station(opts, pattern, station, coords, area, scores, start_time, start, forward, tle_file):
    """Calculate passes, graph, and schedule for one station."""

    logger.debug("station: %s coords: %s area: %s scores: %s", station, coords, area.area_id, scores)

    pattern_args = {
            "station":station,
            "output_dir":opts.output_dir,
            "date":start_time.strftime("%Y%m%d"),
            "time":start_time.strftime("%H%M%S")
            }
    if opts.xml:
        pattern_args['mode'] = "request"
    elif opts.report:
        pattern_args['mode'] = "report"

    satellites = scores.keys()

    if opts.lon and opts.lat and opts.alt:
        coords = (opts.lon, opts.lat, opts.alt)


    logger.info("Computing next satellite passes")
    allpasses = get_next_passes(satellites, start_time,
                                forward, coords, tle_file, aqua_terra_dumps=opts.no_aqua_terra_dump)
    logger.info("Computation of next overpasses done")

    logger.debug(str(sorted(allpasses, key=lambda x: x.risetime)))

    area_boundary = AreaDefBoundary(area, frequency=500)
    area.poly = area_boundary.contour_poly

    if opts.plot:
        logger.info("Saving plots to %s", build_filename("dir_plots", pattern, pattern_args))
        from threading import Thread
        image_saver = Thread(
                             target=save_passes,
                             args=(allpasses,
                                   area.poly,
                                   build_filename("dir_plots", pattern, pattern_args)
                                   )
                             )
        image_saver.start()


    if opts.avoid is not None:
        avoid_list = get_passes_from_xml_file(opts.avoid)
    else:
        avoid_list = None

    logger.info("computing best schedule for area %s" % area.area_id)
    schedule, (graph, labels) = get_best_sched(allpasses, area, scores,
                                            timedelta(seconds=opts.delay),
                                            avoid_list)

    logger.debug(pformat(schedule))
    for opass in schedule:
        opass.rec = True
    logger.info("generating file")

    if opts.scisys:
        generate_sch_file(build_filename("file_sci", pattern, pattern_args), allpasses, coords)

    if opts.xml or opts.report:
        url = urlparse.urlparse(opts.output_url or opts.output_dir)
        if url.scheme not in ["file", ""]:
            directory = "/tmp"
        else:
            directory = url.path
        if opts.report:
            logger.info("Waiting for images to be saved...")
            image_saver.join()
            logger.info("Done!")
        if opts.xml or opts.report:
            """Allways create xml-file in request-mode"""
            pattern_args['mode'] = "request"
            xmlfile = generate_xml_file(allpasses,
                                    start_time + timedelta(hours=start),
                                    start_time + timedelta(hours=forward),
                                    build_filename("file_xml", pattern, pattern_args),
                                    station,
                                    False
                                    )
            logger.info("Generated " + str(xmlfile))
            send_file(url, xmlfile)
        if opts.report:
            """'If report-mode was set"""
            pattern_args['mode'] = "report"
            xmlfile = generate_xml_file(allpasses,
                                    start_time + timedelta(hours=start),
                                    start_time + timedelta(hours=forward),
                                    build_filename("file_xml", pattern, pattern_args),
                                    station,
                                    True
                                    )
            logger.info("Generated " + str(xmlfile))

    if opts.graph or opts.comb:
        graph.save(build_filename("file_graph", pattern, pattern_args))
        graph.export(
                    labels=[str(label) for label in labels],
                    filename=build_filename("file_graph", pattern, pattern_args) + ".gv"
                    )
    if opts.comb:
        import pickle
        ph = open(os.path.join(build_filename("dir_output", pattern, pattern_args), "allpasses.%s.pkl" % station), "wb")
        pickle.dump(allpasses, ph)
        ph.close()

    return graph, allpasses


def combined_stations(opts, pattern, station_list, graph, allpasses, start_time, start, forward):
    """The works around the combination of schedules for two or more stations."""

    logger.info("Generating coordinated schedules ...")

    pattern_args = {
            "output_dir":opts.output_dir,
            "date":start_time.strftime("%Y%m%d"),
            "time":start_time.strftime("%H%M%S")
            }
    if opts.xml:
        pattern_args['mode'] = "request"
    elif opts.report:
        pattern_args['mode'] = "report"

    passes = {}
    # reset flag "rec" for all passes.
    try:
        for s, ap in allpasses.items():
            passes[s] = list(ap)
            for p in passes[s]:
                p.rec = False
    except:
        print "s", s
        print "ap", ap
        print "passes[s]", passes[s]
        print "p", p
        raise

    station_meta = {}
    for coords, station, area, scores in station_list:
        station_meta[station] = {'coords':coords, 'area':area, 'scores':scores}

    stats, schedule, (newgraph, newpasses) = get_combined_sched(graph, passes)
#     logger.debug(pformat(schedule))

    for opass in schedule:
        for i, ipass in zip(range(len(opass)), opass):
            if ipass[0] is None:
                continue
            ipass[0].rec = True

    logger.info("generating files")

    if opts.graph:
        # save graph as npz file.
        pattern_args["station"] = "comb"
        newgraph.save(build_filename("file_graph", pattern, pattern_args))

        # collect labels, each with one pass per station.
        # TODO: is there a simpler way?
        clabels = []
        if sys.version_info < (2, 7):
            npasses = dict((s, set()) for s in stats)
        else:
            npasses = {s:set() for s in stats}
        for npass in newpasses:
            cl = []
            for i, s in zip(range(len(stats)), stats):
                if npass[i][0] is None:
                    cl.append("---")
                else:
                    npasses[s].add(npass[i][0])
                    if npass[i][0].rec:
                        cl.append("+ " + str(npass[i][0]))
                    else:
                        cl.append("  " + str(npass[i][0]))
            clabels.append("\\n".join(cl))
        # save graph as gv file for "dot"-plot
        newgraph.export(labels=[str(label) for label in clabels],
                     filename=build_filename("file_graph", pattern, pattern_args) + ".gv")

    for station in passes.keys():
        pattern_args["station"] = station + "-comb"
        logger.info("Create schedule file(s) for %s", station)
        if opts.scisys:
            generate_sch_file(build_filename("file_sci", pattern, pattern_args), passes[station],
                              station_meta[station]['coords'])
        if opts.xml or opts.report:
            pattern_args['mode'] = "request"
            xmlfile = generate_xml_file(passes[station], start_time + timedelta(hours=start),
                                    start_time + timedelta(hours=forward),
                                    build_filename("file_xml", pattern, pattern_args),
                                    station, False)
            logger.info("Generated " + str(xmlfile))
            url = urlparse.urlparse(opts.output_url or opts.output_dir)
            send_file(url, xmlfile)
        if opts.report:
            pattern_args['mode'] = "report"
            xmlfile = generate_xml_file(passes[station], start_time + timedelta(hours=start),
                                    start_time + timedelta(hours=forward),
                                    build_filename("file_xml", pattern, pattern_args),
                                    station, True)
            logger.info("Generated " + str(xmlfile))

    logger.info("Finished coordinated schedules.")


def run():
    """The schedule command
    """
    import argparse
    global logger

    parser = argparse.ArgumentParser()
    # general arguments
    parser.add_argument("-c", "--config", default=None,
                        help="configuration file to use")
    parser.add_argument("-t", "--tle", default=None,
                        help="tle file to use")
    parser.add_argument("-l", "--log", default=None,
                        help="File to log to (defaults to stdout)")
    parser.add_argument("-m", "--mail", nargs="*", default=None,
                        help="mail address(es) to send error messages to.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="print debug messages too")
    # argument group: coordinates and times
    group_postim = parser.add_argument_group(title="start-parameter",
                        description="(or set values in the configuration file)")
    group_postim.add_argument("--lat", type=float,
                        help="Latitude, degrees north")
    group_postim.add_argument("--lon", type=float,
                        help="Longitude, degrees east")
    group_postim.add_argument("--alt", type=float,
                        help="Altitude, km")
    group_postim.add_argument("-f", "--forward", type=float,
                        help="time ahead to compute the schedule")
    group_postim.add_argument("-s", "--start-time", type=parse_datetime,
                        help="start time of the schedule to compute")
    group_postim.add_argument("-d", "--delay", default=60, type=float,
                        help="delay (in seconds) needed between two "
                        + "consecutive passes (60 seconds by default)")
    # argument group: special behaviour
    group_spec = parser.add_argument_group(title="special",
                        description="(additional parameter changing behaviour)")
    group_spec.add_argument("-a", "--avoid",
                        help="xml request file with passes to avoid")
    group_spec.add_argument("--no-aqua-terra-dump", action="store_false",
                       help="do not consider Aqua/Terra-dumps")
    group_spec.add_argument("--multiproc", action="store_true",
                        help="use multiple parallel processes")
    # argument group: output-related
    group_outp = parser.add_argument_group(title="output",
                        description="(file pattern are taken from configuration file)")
    group_outp.add_argument("-o", "--output-dir", default=None,
                        help="where to put generated files")
    group_outp.add_argument("-u", "--output-url", default=None,
                        help="URL where to put generated schedule file(s)"
                        + ", otherwise use output-dir")
    group_outp.add_argument("-x", "--xml", action="store_true",
                       help="generate an xml request file (schedule)"
                       )
    group_outp.add_argument("-r", "--report", action="store_true",
                       help="generate an xml report file (schedule)")
    group_outp.add_argument("--scisys", action="store_true",
                       help="generate a SCISYS schedule file")
    group_outp.add_argument("-p", "--plot", action="store_true",
                        help="generate plot images")
    group_outp.add_argument("-g", "--graph", action="store_true",
                        help="save graph info")
    opts = parser.parse_args()

    if opts.config:
        # read_config() returns:
        #     [(coords, station, area, scores)], forward, start, {pattern}
        station_list, forward, start, pattern = read_config(opts.config)

    if (not opts.config) and (not (opts.lon or opts.lat or opts.alt)):
        parser.error("Coordinates must be provided in the absence of "
                     "configuration file.")

    if not (opts.xml or opts.scisys or opts.report):
        parser.error("No output specified, use '--scisys' or '-x/--xml'")

    if opts.output_dir is None:
        opts.output_dir = os.path.curdir
    if "dir_output" not in pattern:
        pattern["dir_output"] = opts.output_dir

    if opts.log:
        previous = os.path.exists(opts.log)
        handler = logging.handlers.RotatingFileHandler(opts.log, backupCount=7)
        if previous:
            handler.doRollover()
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[%(levelname)s: %(asctime)s :"
                                           " %(name)s] %(message)s",
                                           '%Y-%m-%d %H:%M:%S'))
    if opts.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    handler.setLevel(loglevel)
    logging.getLogger('').setLevel(loglevel)
    logging.getLogger('').addHandler(handler)

    if opts.mail:
        mhandler = logging.handlers.SMTPHandler("localhost",
                                                "martin.raspaud@smhi.se",
                                                opts.mail,
                                                "Scheduler")
        mhandler.setLevel(logging.WARNING)
        logging.getLogger('').addHandler(mhandler)

    logger = logging.getLogger("trollsched")

    tle_file = opts.tle
    if opts.forward:
        forward = opts.forward
    if opts.start_time:
        start_time = opts.start_time
    else:
        start_time = datetime.utcnow()

    allpasses = {}
    graph = {}

    logger.debug("start: %s forward: %s" % (start, forward))

    pattern_args = {
                    "output_dir":opts.output_dir,
                    "date":start_time.strftime("%Y%m%d"),
                    "time":start_time.strftime("%H%M%S")
                    }
    dir_output = build_filename("dir_output", pattern, pattern_args)
    if not os.path.exists(dir_output):
        logger.debug("Create output dir " + dir_output)
        os.makedirs(dir_output)

    if len(station_list) > 1:
        opts.comb = True
        import pickle
        ph = open(os.path.join(dir_output, "opts.pkl"), "wb")
        pickle.dump(opts, ph)
        ph.close()
    else:
        opts.comb = False

    # single- or multi-processing?
    if not opts.multiproc:
        # sequential processing all stations' single schedule.
        for coords, station, area, scores in station_list:
            graph[station], allpasses[station] = single_station(opts, pattern, station, coords,
                                                                area, scores, start_time, start,
                                                                forward, tle_file)

    else:
        # processing the stations' single schedules with multiprocessing.

        process_single = {}
        statlst_ordered = []

        # first round through the stations, forking sub-processes to do the
        # "single station calculations" in parallel.
        # the pickling of passes and graphs is done inside single_station().
        for coords, station, area, scores in station_list:
            statlst_ordered.append(station)
            from multiprocessing import Process
            process_single[station] = Process(
                    target=single_station,
                    args=(
                          opts, pattern, station, coords,
                          area, scores, start_time, start, forward, tle_file
                          )
                    )
            process_single[station].start()

        # second round through the stations, collecting the sub-processes and their results.
        for station in statlst_ordered:
            process_single[station].join()
            pattern_args["station"] = station
            # load graph for station
            graph[station] = Graph()
            graph[station].load(build_filename("file_graph", pattern, pattern_args) + ".npz")
            # load pickled passes for station
            ph = open(os.path.join(dir_output, "allpasses.%s.pkl" % station), "rb")
            allpasses[station] = pickle.load(ph)
            ph.close()

    if opts.comb:
        combined_stations(opts, pattern, station_list, graph, allpasses, start_time, start, forward)


if __name__ == '__main__':
    try:
        run()
    except:
        logger.exception("Something wrong happened!")
        raise
