#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

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

"""Scheduling 
"""
from tempfile import mkstemp
import logging
import logging.handlers
import operator
import os
from datetime import datetime, timedelta
from pprint import pformat
from scipy.optimize import brentq
import numpy as np
from pyresample import utils
from pyorbital import (orbital, geoloc, geoloc_instrument_definitions,
                       astronomy, tlefile)
from trollsched.spherical import SphPolygon, get_twilight_poly
from trollsched.graph import Graph

from ConfigParser import ConfigParser
import glob

logger = logging.getLogger(__name__)

# shortest allowed pass in minutes
MIN_PASS = 4

class Mapper(object):
    """A class to generate nice plots with basemap.
    """
    def __init__(self):
        from mpl_toolkits.basemap import Basemap

        self.map = Basemap(projection='nsper', lat_0 = 58, lon_0 = 16,
                           resolution = 'l', area_thresh = 1000.)

        self.map.drawcoastlines()
        self.map.drawcountries()

        self.map.drawmapboundary(fill_color='white')

        self.map.drawmeridians(np.arange(0, 360, 30))
        self.map.drawparallels(np.arange(-90, 90, 30))

    def __enter__(self):
        return self.map

    def __exit__(self, etype, value, tb):
        pass

class Boundary(object):
    """Area boundary objects.
    """
    def __init__(self, *sides):
        self.sides_lons, self.sides_lats = zip(*sides)
        self.sides_lons = list(self.sides_lons)
        self.sides_lats = list(self.sides_lats)
        
        self._contour_poly = None

    def decimate(self, ratio):
        """Remove some points in the boundaries, but never the corners.
        """
        for i in range(len(self.sides_lons)):
            l = len(self.sides_lons[i])
            start = (l % ratio) / 2
            points = np.concatenate(([0], np.arange(start, l, ratio), [l-1]))
            
            self.sides_lons[i] = self.sides_lons[i][points]
            self.sides_lats[i] = self.sides_lats[i][points]



    def contour(self):
        """Get the (lons, lats) tuple of the boundary object.
        """
        lons = np.concatenate([lns[:-1] for lns in self.sides_lons])
        lats = np.concatenate([lts[:-1] for lts in self.sides_lats])

        return lons, lats

    def contour_poly(self):
        """Get the Spherical polygon corresponding to the Boundary
        """
        if self._contour_poly is None:
            self._contour_poly = SphPolygon(np.deg2rad(np.vstack(self.contour()).T))
        return self._contour_poly
    
    def draw(self, mapper, options):
        """Draw the current boundary on the *mapper*
        """
        self.contour_poly().draw(mapper, options)

    def show(self):
        """Show the current boundary.
        """
        
        with Mapper() as mapper:
            self.draw(mapper, "-r")

class SwathBoundary(Boundary):
    """Boundaries for satellite overpasses.
    """
    def get_instrument_points(self, overpass, utctime,
                              scans_nb, scanpoints, decimate=1):
        """Get the boundary points for a given overpass.
        """
        instrument = overpass.instrument
        # cheating at the moment.
        scan_angle = 55.37
        if instrument == "modis":
            scan_angle = 55.0
        elif instrument == "viirs":
            scan_angle = 55.84
        elif overpass.satellite == "noaa 16":
            scan_angle = 55.25
        instrument = "avhrr"
        instrument_fun = getattr(geoloc_instrument_definitions, instrument)
        sgeom = instrument_fun(scans_nb, scanpoints,
                               scan_angle=scan_angle, decimate=decimate)
        times = sgeom.times(utctime)
        pixel_pos = geoloc.compute_pixels((self.orb.tle._line1,
                                           self.orb.tle._line2),
                                          sgeom, times)
        lons, lats, alts = geoloc.get_lonlatalt(pixel_pos, times)
        
        del alts
        return (lons.reshape(-1, len(scanpoints)),
                lats.reshape(-1, len(scanpoints)))

    def __init__(self, overpass):
        # compute area covered by pass

        self.overpass = overpass
        self.orb = overpass.orb

        decimate = 500.0

        ## compute sides

        scans_nb = np.ceil(((overpass.falltime - overpass.risetime).seconds +
                            (overpass.falltime - overpass.risetime).microseconds
                            / 1000000.0) * 6 / decimate)

        sides_lons, sides_lats = self.get_instrument_points(self.overpass,
                                                            overpass.risetime,
                                                            scans_nb,
                                                            np.array([0, 2047]),
                                                            decimate=decimate)

        self.left_lons = sides_lons[::-1, 0]
        self.left_lats = sides_lats[::-1, 0]
        self.right_lons = sides_lons[:, 1]
        self.right_lats = sides_lats[:, 1]

        ## compute bottom

        # avhrr
        maxval = 2048
        rest = maxval % decimate
        reduced = np.hstack([0, np.arange(rest/2, maxval, decimate), maxval -1])

        lons, lats = self.get_instrument_points(self.overpass,
                                                overpass.falltime,
                                                1,
                                                reduced)

        self.bottom_lons = lons[0][::-1]
        self.bottom_lats = lats[0][::-1]

        ## compute top

        lons, lats = self.get_instrument_points(self.overpass,
                                                overpass.risetime,
                                                1,
                                                reduced)

        self.top_lons = lons[0]
        self.top_lats = lats[0]

        self._contour_poly = None

    def decimate(self, ratio):
        l = len(self.top_lons)
        start = (l % ratio) / 2
        points = np.concatenate(([0], np.arange(start, l, ratio), [l-1]))

        self.top_lons = self.top_lons[points]
        self.top_lats = self.top_lats[points]
        self.bottom_lons = self.bottom_lons[points]
        self.bottom_lats = self.bottom_lats[points]

        l = len(self.right_lons)
        start = (l % ratio) / 2
        points = np.concatenate(([0], np.arange(start, l, ratio), [l-1]))

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



class Pass(object):
    """A pass: satellite, risetime, falltime, (orbital)
    """

    buffer = timedelta(minutes=2)

    def __init__(self, satellite, risetime, falltime, orb=None, uptime=None, instrument=None):
        """
        
        Arguments:
        - `satellite`:
        - `risetime`:
        - `falltime`:
        """
        self.satellite = satellite
        self.risetime = risetime
        self.falltime = falltime
        self.uptime = uptime
        self.instrument = instrument
        self.orb = orb
        self.score = {}
        self.old_bound = None
        self.boundary = SwathBoundary(self)
        # make boundary lighter.
        #self.boundary.decimate(100)
        self.subsattrack = {"start": None,
                            "end": None}
        self.rec = False
        
    def overlaps(self, other, delay=timedelta(seconds=0)):
        """Check if two passes overlap.
        """
        return ((self.risetime < other.falltime + delay) and
                (self.falltime + delay > other.risetime))

    def __cmp__(self, other):
        if self.uptime < other.uptime:
            return -1
        if self.uptime > other.uptime:
            return 1
        else:
            return 0

    def __eq__(self, other):
        return (self.risetime == other.risetime and
                self.falltime == other.falltime and
                self.satellite == other.satellite)

    def __str__(self):
        return (self.satellite + " "
                + self.risetime.isoformat() + " " + self.falltime.isoformat())

    def __repr__(self):
        return str(self)

    def duration(self):
        """Get the duration of an overpass.
        """
        return self.falltime - self.risetime

    def seconds(self):
        """Get the duration of an overpass.
        """
        duration = self.duration()
        return (duration.days * 24 * 60 * 60
                + duration.seconds
                + duration.microseconds * 1e-6)

    def pass_direction(self):
        """Get the direction of the pass in (ascending, descending).
        """
        start_lat = self.orb.get_lonlatalt(self.risetime)[1]
        end_lat = self.orb.get_lonlatalt(self.falltime)[1]

        if start_lat > end_lat:
            return "descending"
        else:
            return "ascending"

    def slsearch(self, sublat):
        """Find sublatitude.
        """
    
        def nadirlat(minutes):
            return self.orb.get_lonlatalt(self.risetime +
                                          timedelta(minutes=minutes))[1] - sublat

        arr = np.array([nadirlat(m) for m in range(15)])
        a = np.where(np.diff(np.sign(arr)))[0]
        for guess in a:
            sublat_mins = brentq(nadirlat, guess, guess + 1)
            return self.risetime + timedelta(minutes=sublat_mins)
        
    def area_coverage(self, area_of_interest):
        """Get the score depending on the coverage of the area of interest.
        """
        inter = self.boundary.contour_poly().intersection(area_of_interest.poly)
        return inter.area() / area_of_interest.poly.area()


    def save_fig(self, poly=None, directory="/tmp/plots",
                 overwrite=False, labels=None):
        """Save the pass as a figure. Filename is automatically generated.
        """
        logger.debug("Save fig " + str(self))
        filename = os.path.join(directory,
                                (self.risetime.isoformat()
                                 + self.satellite
                                 + self.falltime.isoformat()) + ".pdf")
        if not overwrite and os.path.exists(filename):
            return filename
        
        import matplotlib.pyplot as plt
        plt.clf()
        with Mapper() as mapper:
            mapper.nightshade(self.uptime, alpha=0.2)
            self.draw(mapper, "-r")
            if poly is not None:
                poly.draw(mapper, "-b")
        plt.title(str(self))
        for label in (labels or []):
            plt.figtext(*label[0], **label[1])
        plt.savefig(filename)
        return filename

    def show(self, poly=None, labels=None):
        """Show the current pass on screen (matplotlib, basemap).
        """
        import matplotlib.pyplot as plt
        plt.clf()
        with Mapper() as mapper:
            mapper.nightshade(self.uptime, alpha=0.2)
            self.draw(mapper, "-r")
            if poly is not None:
                poly.draw(mapper, "-b")
        plt.title(str(self))
        for label in (labels or []):
            plt.figtext(*label[0], **label[1])
        plt.show()

    def draw(self, mapper, options):
        """Draw the pass to the *mapper* object (basemap).
        """
        self.boundary.contour_poly().draw(mapper, options)

    def print_vcs(self, coords):
        """Should look like this::


#SCName          RevNum Risetime        Falltime        Elev Dura ANL   Rec Dir Man Ovl OvlSCName        OvlRev OvlRisetime     OrigRisetime    OrigFalltime    OrigDuration
#
NOAA 19           24845 20131204 001450 20131204 003003 32.0 15.2 225.6 Y   Des N   N   none                  0 19580101 000000 20131204 001450 20131204 003003 15.2

        
        """

        max_elevation = self.orb.get_observer_look(self.uptime, *coords)[1]
        anl = self.orb.get_observer_look(self.risetime, *coords)[0]
        if self.rec:
            rec = "Y"
        else:
            rec = "N"
        line_list = ["{satellite:<16}",
                     "{orbit:>6}",
                     "{risetime}",
                     "{falltime}",
                     "{elevation:>4.1f}",
                     "{duration:>4.1f}",
                     "{anl:>5.1f}",
                     "{rec:<3}",
                     "{direction}",
                     "N   N   none                  0 19580101 000000",
                     "{risetime}",
                     "{falltime}",
                     "{duration:>4.1f}",
                     ]                     
        line = " ".join(line_list).format(
            satellite=self.satellite.upper(),
            orbit=self.orb.get_orbit_number(self.risetime),
            risetime=self.risetime.strftime("%Y%m%d %H%M%S"),
            falltime=self.falltime.strftime("%Y%m%d %H%M%S"),
            elevation=max_elevation,
            duration=(self.falltime - self.risetime).seconds / 60.0,
            anl=anl,
            rec=rec,
            direction=self.pass_direction().capitalize()[:3])
        return line
        
def conflicting_passes(allpasses, delay=timedelta(seconds=0)):
    """Get the passes in groups of conflicting passes.
    """

    passes = sorted(allpasses)

    last_time = None
    group = []
    groups = []
    for overpass in passes:
        if last_time is None:
            last_time = overpass.falltime
            group.append(overpass)
            continue
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

    for i, overpass in enumerate(passes):
        for j in range(i+1, order):
            if not overpass.overlaps(passes[j], delay):
                graph.add_edge(i, j)
    
    groups = []
    for res in graph.bron_kerbosch(set(), set(graph.vertices), set()):
        grp = []
        for v in res:
            grp.append(passes[v])
        groups.append(grp)

    return groups


def fermia(t):
    a = 0.25
    b = a / 4
    k = b * np.log(1/3.0) + a
    sh = k - 0.25
    return 0.5/(np.exp(((t + sh) - a) / b) + 1) + 0.5

def fermib(t):
    a = 0.25
    b = a / 4
    return 1/(np.exp((t - a) / b) + 1)

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
        ip1 = p1.boundary.contour_poly().intersection(area_of_interest.poly)
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
        ip2 = p2.boundary.contour_poly().intersection(area_of_interest.poly)
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


def get_best_sched(overpasses, area_of_interest, scores, delay):
    """Get the best schedule based on *area_of_interest*.
    """
    passes = sorted(overpasses)
    grs = conflicting_passes(passes, delay)
    ncgrs = [get_non_conflicting_groups(gr, delay) for gr in grs]

    n_vertices = len(passes)
    
    graph = Graph(n_vertices=n_vertices + 2)


    def add_arc(graph, p1, p2, hook=None):
        logger.debug("Adding arc between " + str(p1) +
                     " and " + str(p2) + "...")
        w = combine(p1, p2, area_of_interest, scores)
        logger.debug("...with weight " + str(w))
        
        with open("/tmp/schedule.gv", "a") as fp_:
            fp_.write('        "' + str(p1) + '" -> "' + str(p2) +
                      '" [ label = "' + str(w) + '" ];\n')

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

        prev = set(gr[-1] for gr in ncgr)
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
    return [passes[idx - 1] for idx in path[1:-1]], graph
        

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

        out.write("#SCName          RevNum Risetime        Falltime        Elev Dura ANL   Rec Dir Man Ovl OvlSCName        OvlRev OvlRisetime     OrigRisetime    OrigFalltime    OrigDuration\n#\n")
    
        for overpass in sorted(overpasses):
            out.write(overpass.print_vcs(coords) + "\n")
            
HOST = "ftp://is.sci.gsfc.nasa.gov/ancillary/ephemeris/schedule/aqua/downlink/"
import urlparse
import ftplib
import socket

def get_aqua_dumps_from_ftp(start_time, end_time, satorb):
    url = urlparse.urlparse(HOST)
    logger.debug("Connect to ftp server")
    try:
        f = ftplib.FTP(url.netloc)
    except (socket.error, socket.gaierror), e:
        logger.error('cannot reach to %s ' % HOST + str(e))
        f = None


    if f is not None:
        try:
            f.login('anonymous','guest')
            logger.debug("Logged in")
        except ftplib.error_perm:
            logger.error('cannot login anonymously')
            f.quit()
            f = None


    if f is not None:
        data = []        
        try:
            f.dir(url.path, data.append)
        except socket.error, e:
            logger.error("Can't get any data: " + str(e))
            f.quit()
            f = None
        else:
            filenames = [line.split()[-1] for line in data]

    if f is None:
        logger.info("Can't access ftp server, using cached data")
        filenames = glob.glob("/tmp/*.rpt")
        
    dates = [datetime.strptime("".join(filename.split(".")[2:4]), "%Y%j%H%M%S")
             for filename in filenames]
    filedates = dict(zip(dates, filenames))

    dumps = []

    for date in sorted(dates):
        lines = []
        if not filedates[date].endswith(".rpt"):
            continue
        if not os.path.exists(os.path.join("/tmp", filedates[date])):
            f.retrlines('RETR ' + os.path.join(url.path, filedates[date]), lines.append)
            with open(os.path.join("/tmp", filedates[date]), "w") as fd_:
                for line in lines:
                    fd_.write(line + "\n")
        else:
            with open(os.path.join("/tmp", filedates[date]), "r") as fd_:
                for line in fd_:
                    lines.append(line)

        for line in lines[7::2]:
            if line.strip() == '':
                break
            station, aos, elev, los = line.split()[:4]
            aos = datetime.strptime(aos, "%Y:%j:%H:%M:%S")
            los = datetime.strptime(los, "%Y:%j:%H:%M:%S")
            if los >= start_time and aos <= end_time:
                uptime = aos + (los - aos) / 2
                overpass = Pass("aqua", aos, los, satorb, uptime, "modis")
                overpass.station = station
                overpass.max_elev = elev
                dumps.append(overpass)
    if f is not None:
        f.quit()
    return dumps
            

def get_next_passes(satellites, utctime, forward, coords, tle_file=None):
    passes = {}
    orbitals = {}

    if tle_file is None:
        fp_, tle_file = mkstemp(prefix="tle", dir="/tmp")
        os.close(fp_)
        logger.info("Fetch tle info from internet")
        tlefile.fetch(tle_file)
        
    if not os.path.exists(tle_file):
        logger.info("Fetch tle info from internet")
        tlefile.fetch(tle_file)

    for sat in satellites:
        satorb = orbital.Orbital(sat, tle_file=tle_file)
        orbitals[sat] = satorb
        passlist = satorb.get_next_passes(utctime,
                                          forward,
                                          *coords,
                                          tol=0.001)
        if sat.startswith("metop") or sat.startswith("noaa"):
            instrument = "avhrr"
        elif sat in ["aqua", "terra"]:
            instrument = "modis"
        elif sat.endswith("npp"):
            instrument = "viirs"
        else:
            instrument = "unknown"
        # take care of metop-a
        if sat == "metop-a":
            metop_passes = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                            for rtime, ftime, uptime in passlist]

            passes["metop-a"] = []
            for overpass in metop_passes:
                if overpass.pass_direction() == "descending":
                    new_rise = overpass.slsearch(60)
                    if new_rise is not None and new_rise < overpass.falltime:
                        overpass.risetime = new_rise
                        overpass.boundary = SwathBoundary(overpass)
                        if overpass.seconds() > MIN_PASS * 60:
                            passes["metop-a"].append(overpass)
        # take care of aqua (dumps in svalbard and poker flat)
        elif sat == "aqua":
            
            wpcoords = (-75.457222, 37.938611, 0)
            passlist_wp = satorb.get_next_passes(utctime - timedelta(minutes=30),
                                                 forward + 1,
                                                 *wpcoords,
                                                 tol=0.001)
            wp_passes = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                         for rtime, ftime, uptime in passlist_wp]
            
            svcoords = (15.399, 78.228, 0)
            passlist_sv = satorb.get_next_passes(utctime - timedelta(minutes=30),
                                                 forward + 1,
                                                 *svcoords,
                                                 tol=0.001)
            sv_passes = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                         for rtime, ftime, uptime in passlist_sv]
            pfcoords = (-147.43, 65.12, 0.51)
            passlist_pf = satorb.get_next_passes(utctime - timedelta(minutes=30),
                                                 forward + 1,
                                                 *pfcoords,
                                                 tol=0.001)
            pf_passes = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                         for rtime, ftime, uptime in passlist_pf]
            
            aqua_passes = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                           for rtime, ftime, uptime in passlist]

            dumps = get_aqua_dumps_from_ftp(utctime - timedelta(minutes=30),
                                            utctime + timedelta(hours=forward+0.5),
                                            satorb)

            # remove the known dumps
            for dump in dumps:
                #print "*", dump.station, dump, dump.max_elev
                logger.debug("dump from ftp: " + str((dump.station, dump,
                                                      dump.max_elev)))
                for i, sv_pass in enumerate(sv_passes):
                    if sv_pass.overlaps(dump, timedelta(minutes=40)):
                        sv_elevation = sv_pass.orb.get_observer_look(sv_pass.uptime,
                                                                     *svcoords)[1]
                        logger.debug("Computed " +str(("SG", sv_pass,
                                                       sv_elevation)))
                        del sv_passes[i]
                for i, pf_pass in enumerate(pf_passes):
                    if pf_pass.overlaps(dump, timedelta(minutes=40)):
                        pf_elevation = pf_pass.orb.get_observer_look(pf_pass.uptime,
                                                                     *pfcoords)[1]
                        logger.debug("Computed " +str(("PF", pf_pass,
                                                       pf_elevation)))
                        del pf_passes[i]
                for i, wp_pass in enumerate(wp_passes):
                    if wp_pass.overlaps(dump, timedelta(minutes=40)):
                        wp_elevation = wp_pass.orb.get_observer_look(wp_pass.uptime,
                                                                     *wpcoords)[1]
                        logger.debug("Computed " +str(("WP", wp_pass,
                                                       wp_elevation)))
                        del wp_passes[i]
                        
            # sort out dump passes first
            # between sv an pf, we take the one with the highest elevation if
            # pf < 20°, pf otherwise
            # I think wp is also used if sv is the only other alternative
            used_pf = []
            for sv_pass in sv_passes:
                found_pass = False
                for pf_pass in pf_passes:
                    if sv_pass.overlaps(pf_pass):
                        found_pass = True
                        used_pf.append(pf_pass)
                        sv_elevation = sv_pass.orb.get_observer_look(sv_pass.uptime,
                                                                     *svcoords)[1]
                        pf_elevation = pf_pass.orb.get_observer_look(pf_pass.uptime,
                                                                     *pfcoords)[1]
                        if pf_elevation > 20:
                            dumps.append(pf_pass)
                        elif sv_elevation > pf_elevation:
                            dumps.append(sv_pass)
                        else:
                            dumps.append(pf_pass)
                        break
                if not found_pass:
                    dumps.append(sv_pass)

            for pf_pass in pf_passes:
                if pf_pass not in used_pf:
                    dumps.append(pf_pass)


            passes["aqua"] = []
            for overpass in aqua_passes:
                add = True
                for dump_pass in dumps:
                    if dump_pass.overlaps(overpass):
                        if (dump_pass.uptime < overpass.uptime and
                            dump_pass.falltime > overpass.risetime):
                            logger.debug("adjusting " + str(overpass)
                                         + " to new risetime " +
                                         str(dump_pass.falltime))
                            overpass.risetime = dump_pass.falltime
                            overpass.boundary = SwathBoundary(overpass)
                        elif (dump_pass.uptime >= overpass.uptime and
                              dump_pass.risetime < overpass.falltime):
                            logger.debug("adjusting " + str(overpass)
                                         + " to new falltime " +
                                         str(dump_pass.risetime))
                            overpass.falltime = dump_pass.risetime
                            overpass.boundary = SwathBoundary(overpass)
                        if overpass.falltime <= overpass.risetime:
                            add = False
                            logger.debug("skipping " + str(overpass))
                if add and overpass.seconds() > MIN_PASS * 60:                    
                    passes["aqua"].append(overpass)

        else:
            passes[sat] = [Pass(sat, rtime, ftime, satorb, uptime, instrument)
                           for rtime, ftime, uptime in passlist
                           if ftime - rtime > timedelta(minutes=MIN_PASS)]

            
    return set(reduce(operator.concat, passes.values()))


def generate_xml_requests(sched, start, end, station_name):
    """Create xml requests.
    """
    import xml.etree.ElementTree as ET

    sats = {"noaa 15": "noaa15",
            "noaa 16": "noaa16",
            "noaa 18": "noaa18",
            "noaa 19": "noaa19",
            "metop-a": "metopa",
            "metop-b": "metopb",
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
    typep.text = "request"
    station = ET.SubElement(props, "station")
    station.text = station_name
    file_start = ET.SubElement(props, "file-start")
    file_start.text = start.strftime(eum_format)
    file_end = ET.SubElement(props, "file-end")
    file_end.text = end.strftime(eum_format)
    reqby = ET.SubElement(props, "requested-by")
    reqby.text = "SMHI"
    reqon = ET.SubElement(props, "requested-on")
    reqon.text = reqtime.strftime(eum_format)
    for overpass in sorted(sched):
        if overpass.rec:
            ovpass = ET.SubElement(root, "pass")
            ovpass.set("satellite", sats[overpass.satellite])
            ovpass.set("start-time", overpass.risetime.strftime(eum_format))
            ovpass.set("end-time", overpass.falltime.strftime(eum_format))

    return root, reqtime

def generate_xml_file(sched, start, end, directory, station):
    """Create an xml request file.
    """
    import xml.etree.ElementTree as ET
    tree, reqtime = generate_xml_requests(sched, start, end, station)
    filename = (reqtime.strftime("%Y-%m-%d-%H-%M-%S")
                + "-acquisition-schedule-request-"
                + station + ".xml")
    filename = os.path.join(directory, filename)
    with open(filename, "w") as fp_:
        fp_.write(ET.tostring(tree))
    return filename

def parse_datetime(strtime):
    return datetime.strptime(strtime, "%Y%m%d%H%M%S")

def read_config(filename):
    """Read the config file *filename* and replace the values in global
    variables.
    """
    cfg = ConfigParser()
    cfg.read(filename)

    station = cfg.get("default", "station")
    satellites = cfg.get("default", "satellites").split(",")
    forward = cfg.getfloat("default", "forward")
    start = cfg.getfloat("default", "start")

    station_name = cfg.get(station, "name")
    station_lon = cfg.getfloat(station, "longitude")
    station_lat = cfg.getfloat(station, "latitude")
    station_alt = cfg.getfloat(station, "altitude")

    sat_scores = {}
    for sat in satellites:
        sat_scores[sat] = (cfg.getfloat(sat, "night"),
                           cfg.getfloat(sat, "day"))
    
    area = utils.parse_area_file(cfg.get(station, "area_file"),
                                 cfg.get(station, "area"))[0]

    
    return ((station_lon, station_lat, station_alt),
            sat_scores, station_name, area, forward, start)
        
def run():
    """The schedule command
    """
    import argparse
    global logger
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--lon", help="Longitude, degrees east", type=float)
    parser.add_argument("--lat", help="Latitude, degrees north", type=float)
    parser.add_argument("--alt", help="Altitude, km", type=float)
    parser.add_argument("-l", "--log",
                        help="File to log to (defaults to stdout)",
                        default=None)
    parser.add_argument("-m", "--mail", nargs="*",
                        help="mail address(es) to send error messages to.",
                        default=None)
    parser.add_argument("-v", "--verbose", help="print debug messages too",
                        action="store_true")
    parser.add_argument("-t", "--tle", help="tle file to use", default=None)
    parser.add_argument("-f", "--forward", type=float,
                        help="time ahead to compute the schedule")
    parser.add_argument("-s", "--start-time",
                        type=parse_datetime,
                        help="start time of the schedule to compute")
    parser.add_argument("-d", "--delay", default=60, type=float,
                        help="delay (in seconds) needed between two "
                        + "consecutive passes (60 seconds by default)")
    parser.add_argument("-c", "--config", help="configuration file to use",
                        default=None)
    group = parser.add_argument_group(title="output")
    group.add_argument("-x", "--xml", default=".",
                        help="generate an xml request file and"
                       " put it in this directory. Could be a url")
    group.add_argument("--scisys", default=None,
                        help="path to the schedule file")
                        
    opts = parser.parse_args()

    if opts.config:
        coords, scores, station, area, forward, start = read_config(opts.config)

    if (not opts.config) and (not (opts.lon or opts.lat or opts.alt)):
        parser.error("Coordinates must be provided in the absence of "
                     "configuration file.")
        
    if not (opts.xml or opts.scisys):
        parser.error("No output specified, use '--scisys' or '-x/--xml'")

    if opts.log:
        previous = os.path.exists(opts.log)
        handler = logging.handlers.RotatingFileHandler(opts.log,
                                                       backupCount=7)
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

    if opts.lon and opts.lat and opts.alt:
        coords = (opts.lon, opts.lat, opts.alt)


    # test line
    # python schedule.py -v 16.148649 58.581844 0.052765 -f 216 -s 20140118140000 -t tle_20140120.txt -x . --scisys myched.txt

    satellites = ["noaa 19", "noaa 18", "noaa 16", "noaa 15",
                  "metop-a", "metop-b",
                  "terra", "aqua",
                  "suomi npp"]


    logger.info("Computing next satellite passes")
    
    tle_file = opts.tle
    if opts.forward:
        forward = opts.forward
    if opts.start_time:
        start_time = opts.start_time
    else:
        start_time = datetime.utcnow() + timedelta(hours=start)
    allpasses = get_next_passes(satellites, start_time,
                                forward, coords, tle_file)


    logger.info("Computation of next overpasses done")

    logger.debug(str(sorted(allpasses)))

    lons, lats = area.get_boundary_lonlats()
    area_boundary = Boundary((lons.side1, lats.side1),
                             (lons.side2, lats.side2),
                             (lons.side3, lats.side3),
                             (lons.side4, lats.side4))
    area_boundary.decimate(500)
    area.poly = area_boundary.contour_poly()

    logger.info("computing best schedule for area euron1")
    schedule, graph = get_best_sched(allpasses, area, scores,
                                     timedelta(seconds=opts.delay))
    

    logger.debug(pformat(schedule))
    for opass in schedule:
        opass.rec = True
    logger.info("generating file")
    
    if opts.scisys:
        generate_sch_file(opts.scisys, allpasses, coords)

    if opts.xml:
        url = urlparse.urlparse(opts.xml)
        if url.scheme not in ["file", ""]:
            directory = "/tmp"
        else:
            directory = url.path
        xmlfile = generate_xml_file(allpasses, start_time,
                                    start_time + timedelta(hours=forward),
                                    directory, station)
        logger.info("Generated " + str(xmlfile))
        pathname, filename = os.path.split(xmlfile)
        del pathname
        if url.scheme in ["file", ""]:
            pass
        elif url.scheme == "ftp":
            session = ftplib.FTP(url.hostname, url.username, url.password)
            with open(xmlfile, "rb") as xfile:
                session.storbinary('STOR ' + str(filename), xfile)
            session.quit()
        else:
            logger.error("Cannot save to " + str(url.scheme)
                         + ", but file is there" + str(xmlfile))
        
    #graph.save("my_graph")
if __name__ == '__main__':
    try:
        run()
    except:
        logger.exception("Something wrong happened!")
        raise

