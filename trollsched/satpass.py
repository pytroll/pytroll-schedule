#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 - 2019 PyTroll Community
# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Alexander Maul <alexander.maul@dwd.de>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

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
"""Satellite passes."""

import ftplib
import glob
import logging
import logging.handlers
import operator
import os
import socket
from datetime import datetime, timedelta
from functools import reduce as fctools_reduce
from tempfile import gettempdir, mkstemp
from urllib.parse import urlparse

import numpy as np
from pyorbital import orbital, tlefile
from pyresample.boundary import AreaDefBoundary

from trollsched import MIN_PASS, NOAA20_NAME, NUMBER_OF_FOVS
from trollsched.boundary import SwathBoundary

logger = logging.getLogger(__name__)

VIIRS_PLATFORM_NAMES = ["SUOMI NPP", "SNPP",
                        "NOAA-20", "NOAA 20"]
MERSI_PLATFORM_NAMES = ["FENGYUN 3C", "FENGYUN-3C", "FY-3C"]
MERSI2_PLATFORM_NAMES = ["FENGYUN 3D", "FENGYUN-3D", "FY-3D",
                         "FENGYUN 3E", "FENGYUN-3E", "FY-3E"]


class SimplePass:
    """A pass: satellite, risetime, falltime, (orbital)."""

    buffer = timedelta(minutes=2)

    def __init__(self, satellite, risetime, falltime):
        """Initialize the simple pass."""
        if not hasattr(satellite, "name"):
            from trollsched.schedule import Satellite
            self.satellite = Satellite(satellite, 0, 0)
        else:
            self.satellite = satellite
        self.risetime = risetime
        self.falltime = falltime
        self.score = {}
        self.subsattrack = {"start": None, "end": None}
        self.rec = False
        self.fig = None

    def __hash__(self):
        """Hash the pass."""
        return super.__hash__(self)

    def overlaps(self, other, delay=None):
        """Check if two passes overlap in time."""
        if delay is None:
            delay = timedelta(seconds=0)
        return ((self.risetime < other.falltime + delay) and (self.falltime + delay > other.risetime))

    def __lt__(self, other):
        """Check if this pass starts earlier than the other one."""
        return self.uptime < other.uptime

    def __gt__(self, other):
        """Check if this pass startss later than the other one."""
        return self.uptime > other.uptime

    def __cmp__(self, other):
        """Compare two passes."""
        if self.uptime < other.uptime:
            return -1
        if self.uptime > other.uptime:
            return 1
        else:
            return 0

    def __eq__(self, other):
        """Determine if two satellite passes are the same."""
        # Two passes, maybe observed from two distinct stations, are compared by
        # a) satellite name and orbit number,
        # or if the later is not available
        # b) the time difference between rise- and fall-times.
        if other is not None and isinstance(self, Pass) and isinstance(
                other, Pass):
            return (self.satellite.name == other.satellite.name and
                    self.orb.get_orbit_number(self.risetime) == other.orb.get_orbit_number(other.risetime))
        timedelta(seconds=1)
        return (other is not None and
                self.satellite.name == other.satellite.name and
                self.overlaps(other))

    def __str__(self):
        """Give a string version of the pass."""
        return (self.satellite.name + " " + self.risetime.isoformat() + " " +
                self.falltime.isoformat())

    def __repr__(self):
        """Represent the pass."""
        return str(self)

    def duration(self):
        """Get the duration of an overpass."""
        return self.falltime - self.risetime

    def seconds(self):
        """Get the duration of an overpass."""
        duration = self.duration()
        return (duration.days * 24 * 60 * 60 + duration.seconds +
                duration.microseconds * 1e-6)


class Pass(SimplePass):
    """A pass: satellite, risetime, falltime, (orbital)."""

    def __init__(self, satellite, risetime, falltime, **kwargs):
        """Initialize the pass."""
        SimplePass.__init__(self, satellite, risetime, falltime)

        logger.debug("kwargs: %s", str(kwargs))
        orb = kwargs.get("orb", None)
        uptime = kwargs.get("uptime", None)
        instrument = kwargs.get("instrument", None)
        tle1 = kwargs.get("tle1", None)
        tle2 = kwargs.get("tle2", None)
        logger.debug("instrument: %s", str(instrument))

        if isinstance(instrument, (list, set)):
            if "avhrr" in instrument:
                logger.warning("Instrument is a sequence! Assume avhrr...")
                instrument = "avhrr"
            elif "viirs" in instrument:
                logger.warning("Instrument is a sequence! Assume viirs...")
                instrument = "viirs"
            elif "modis" in instrument:
                logger.warning("Instrument is a sequence! Assume modis...")
                instrument = "modis"
            elif "mersi" in instrument:
                logger.warning("Instrument is a sequence! Assume mersi...")
                instrument = "mersi"
            elif "mersi-2" in instrument:
                logger.warning("Instrument is a sequence! Assume mersi-2...")
                instrument = "mersi-2"
            else:
                raise TypeError("Instrument is a sequence! Don't know which one to choose!")

        default = NUMBER_OF_FOVS.get(instrument, 2048)
        self.number_of_fovs = kwargs.get("number_of_fovs", default)
        # The frequency shouldn't actualy depend on the number of FOVS along a scanline should it!?
        # frequency = kwargs.get('frequency', int(self.number_of_fovs / 4))
        frequency = kwargs.get("frequency", 300)

        self.station = None
        self.max_elev = None
        self.uptime = uptime or (risetime + (falltime - risetime) / 2)
        self.instrument = instrument
        self.frequency = frequency
        if orb:
            self.orb = orb
        else:
            try:
                self.orb = orbital.Orbital(satellite, line1=tle1, line2=tle2)
            except KeyError as err:
                logger.debug("Failed in PyOrbital: %s", str(err))
                self.orb = orbital.Orbital(
                    NOAA20_NAME.get(satellite, satellite),
                    line1=tle1,
                    line2=tle2)
                logger.info("Using satellite name %s instead",
                            str(NOAA20_NAME.get(satellite, satellite)))

        self._boundary = None

    @property
    def boundary(self):
        """Get the boundary of the swath."""
        if not self._boundary:
            self._boundary = SwathBoundary(self, frequency=self.frequency)
        return self._boundary

    @boundary.setter
    def boundary(self, value):
        self._boundary = SwathBoundary(self, frequency=self.frequency)

    def pass_direction(self):
        """Get the direction of the pass in (ascending, descending)."""
        start_lat = self.orb.get_lonlatalt(self.risetime)[1]
        end_lat = self.orb.get_lonlatalt(self.falltime)[1]

        if start_lat > end_lat:
            return "descending"
        else:
            return "ascending"

    def slsearch(self, sublat):
        """Find sublatitude."""

        def nadirlat(minutes):
            return self.orb.get_lonlatalt(self.risetime + timedelta(
                minutes=np.float64(minutes)))[1] - sublat

        def get_root(fun, start, end):
            p = np.polyfit(
                [start, (start + end) / 2.0, end],
                [fun(start), fun((start + end) / 2),
                 fun(end)], 2)
            for root in np.roots(p):
                if root <= end and root >= start:
                    return root

        arr = np.array([nadirlat(m) for m in range(15)])
        a = np.where(np.diff(np.sign(arr)))[0]
        for guess in a:
            sublat_mins = get_root(nadirlat, guess, guess + 1)
            return self.risetime + timedelta(minutes=sublat_mins)

    def area_coverage(self, area_of_interest):
        """Get the ratio of coverage (between 0 and 1) of the pass with the area of interest."""
        try:
            area_boundary = area_of_interest.poly
        except AttributeError:
            area_boundary = AreaDefBoundary(area_of_interest, frequency=100)
            area_boundary = area_boundary.contour_poly

        inter = self.boundary.contour_poly.intersection(area_boundary)

        if inter is None:
            return 0
        return inter.area() / area_boundary.area()

    def generate_metno_xml(self, coords, root):
        """Generate a metno xml schedule."""
        import xml.etree.ElementTree as ET  # noqa because defusedxml has no SubElement

        asimuth_at_max_elevation, max_elevation = self.orb.get_observer_look(self.uptime, *coords)
        pass_direction = self.pass_direction().capitalize()[:1]
        # anl = self.orb.get_lonlatalt(self.orb.get_last_an_time(self.risetime))[0] % 360
        asimuth_at_aos, aos_elevation = self.orb.get_observer_look(self.risetime, *coords)
        orbit = self.orb.get_orbit_number(self.risetime)
        # aos_epoch=int((self.risetime-datetime(1970,1,1)).total_seconds())
        sat_lon, sat_lat, alt = self.orb.get_lonlatalt(self.risetime)

        ovpass = ET.SubElement(root, "pass")
        ovpass.set("satellite", self.satellite.name)
        ovpass.set("aos", self.risetime.strftime("%Y%m%d%H%M%S"))
        ovpass.set("los", self.falltime.strftime("%Y%m%d%H%M%S"))
        ovpass.set("orbit", "{:d}".format(orbit))
        ovpass.set("max-elevation", "{:.3f}".format(max_elevation))
        ovpass.set("asimuth-at-max-elevation", "{:.3f}".format(asimuth_at_max_elevation))
        ovpass.set("asimuth-at-aos", "{:.3f}".format(asimuth_at_aos))
        ovpass.set("pass-direction", pass_direction)
        ovpass.set("satellite-lon-at-aos", "{:.3f}".format(sat_lon))
        ovpass.set("satellite-lat-at-aos", "{:.3f}".format(sat_lat))
        ovpass.set("tle-epoch", self.orb.orbit_elements.epoch.astype(datetime).strftime("%Y%m%d%H%M%S.%f"))
        if self.fig:
            ovpass.set("figure", self.fig)

        return True

    def print_meos(self, coords, line_no):
        """No. Date    Satellite  Orbit Max EL  AOS      Ovlp  LOS      Durtn  Az(AOS/MAX)."""
        asimuth_at_max_elevation, max_elevation = self.orb.get_observer_look(self.uptime, *coords)
        pass_direction = self.pass_direction().capitalize()[:1]
        # anl = self.orb.get_lonlatalt(self.orb.get_last_an_time(self.risetime))[0] % 360
        asimuth_at_aos, aos_elevation = self.orb.get_observer_look(self.risetime, *coords)
        orbit = self.orb.get_orbit_number(self.risetime)
        aos_epoch = int((self.risetime - datetime(1970, 1, 1)).total_seconds())
        sat_lon, sat_lat, alt = self.orb.get_lonlatalt(self.risetime)

        dur_secs = (self.falltime - self.risetime).seconds
        dur_hours, dur_reminder = divmod(dur_secs, 3600)
        dur_minutes, dur_seconds = divmod(dur_reminder, 60)
        duration = "{:0>2}:{:0>2}".format(dur_minutes, dur_seconds)

        satellite_meos_translation = {"NOAA 19": "NOAA_19",
                                      "NOAA 18": "NOAA_18",
                                      "NOAA 15": "NOAA_15",
                                      "METOP-A": "M02",
                                      "METOP-B": "M01",
                                      "FENGYUN 3A": "FENGYUN-3A",
                                      "FENGYUN 3B": "FENGYUN-3B",
                                      "FENGYUN 3C": "FENGYUN-3C",
                                      "SUOMI NPP": "NPP"}

        import hashlib

        pass_key = hashlib.md5(("{:s}|{:d}|{:d}|{:.3f}|{:.3f}".  # noqa : md5 is insecure, but not sensitive here.
                                format(satellite_meos_translation.get(self.satellite.name.upper(),
                                                                      self.satellite.name.upper()),
                                       int(orbit),
                                       aos_epoch,
                                       sat_lon,
                                       sat_lat)).encode("utf-8")).hexdigest()

        line_list = [" {line_no:>2}",
                     "{date}",
                     "{satellite:<10}",
                     "{orbit:>5}",
                     "{elevation:>6.3f} ",
                     "{risetime}",
                     "{overlap:<5s}",
                     "{falltime}",
                     "{duration}",
                     "{asimuth_at_aos:>5.1f}",
                     "{asimuth_at_max:>5.1f}",
                     "-- Undefined(Scheduling not done {aos_epoch} )",
                     "{passkey}",
                     "{pass_direction}"
                     ]

        line = " ".join(line_list).format(
            # line_no=line_no,
            line_no=1,
            date=self.risetime.strftime("%Y%m%d"),
            satellite=satellite_meos_translation.get(self.satellite.name.upper(),
                                                     self.satellite.name.upper()),
            orbit=orbit,
            elevation=max_elevation,
            risetime=self.risetime.strftime("%H:%M:%S"),
            overlap="n/a",
            falltime=self.falltime.strftime("%H:%M:%S"),
            duration=duration,
            asimuth_at_aos=asimuth_at_aos,
            asimuth_at_max=asimuth_at_max_elevation,
            aos_epoch=aos_epoch,
            passkey=pass_key,
            pass_direction=pass_direction)
        return line

    def print_vcs(self, coords):
        """Print a vcs/scisys/cgi schedule.

        Should look like this::

        # SCName          RevNum Risetime        Falltime        Elev Dura ANL   Rec Dir Man Ovl OvlSCName
        #      OvlRev OvlRisetime     OrigRisetime    OrigFalltime    OrigDuration
        # NOAA 19           24845 20131204 001450 20131204 003003 32.0 15.2 225.6 Y   Des N   N   none
        #      0 19580101 000000 20131204 001450 20131204 003003 15.2


        """
        max_elevation = self.orb.get_observer_look(self.uptime, *coords)[1]
        anl = self.orb.get_lonlatalt(self.orb.get_last_an_time(
            self.risetime))[0] % 360
        # anl = self.orb.get_observer_look(self.risetime, *coords)[0]
        if self.rec:
            rec = "Y"
        else:
            rec = "N"
        line_list = [
            "{satellite:<16}",
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
            satellite=self.satellite.name.upper(),
            orbit=self.orb.get_orbit_number(self.risetime),
            risetime=self.risetime.strftime("%Y%m%d %H%M%S"),
            falltime=self.falltime.strftime("%Y%m%d %H%M%S"),
            elevation=max_elevation,
            duration=(self.falltime - self.risetime).seconds / 60.0,
            anl=anl,
            rec=rec,
            direction=self.pass_direction().capitalize()[:3])
        return line


HOST = "ftp://is.sci.gsfc.nasa.gov/ancillary/ephemeris/schedule/%s/downlink/"


def get_aqua_terra_dumps(start_time,
                         end_time,
                         satorb,
                         sat,
                         dump_url=None):
    """Get the Terra and Aqua overpasses.

     We take into account the fact that when
    there are global dumps there is no direct broadcast.
    """
    # Get the list of aqua/terra dump info:
    dump_info_list = get_aqua_terra_dumpdata_from_ftp(sat, dump_url)

    dumps = []
    for elem in dump_info_list:
        if elem["los"] >= start_time and elem["aos"] <= end_time:
            uptime = elem["aos"] + (elem["los"] - elem["aos"]) / 2
            overpass = Pass(sat, elem["aos"], elem["los"],
                            orb=satorb, uptime=uptime, instrument="modis")
            overpass.station = elem["station"]
            overpass.max_elev = elem["elev"]
            dumps.append(overpass)

    return dumps


def get_aqua_terra_dumpdata_from_ftp(sat, dump_url):
    """Get the information on the internet on the actual global dumps of Terra and Aqua."""
    logger.info("Fetch %s dump info from internet", str(sat.name))
    if isinstance(dump_url, str):
        url = urlparse(dump_url % sat.name)
    else:
        url = urlparse(HOST % sat.name)
    logger.debug("Connect to ftp server")
    try:
        f = ftplib.FTP_TLS(url.netloc)
    except (socket.error, socket.gaierror) as e:
        logger.error("cannot reach to %s " % HOST + str(e))
        f = None

    if f is not None:
        try:
            f.login("anonymous", "guest")
            logger.debug("Logged in")
        except ftplib.error_perm:
            logger.error("cannot login anonymously")
            f.quit()
            f = None

    if f is not None:
        data = []
        try:
            f.prot_p()  # explicitly call for protected transfer
            f.dir(url.path, data.append)
        except socket.error as e:
            logger.error("Can't get any data: " + str(e))
            f.quit()
            f = None
        else:
            filenames = [line.split()[-1] for line in data]

    if f is None:
        logger.info("Can't access ftp server, using cached data")
        filenames = glob.glob(os.path.join(gettempdir(), "*.rpt"))

    filenames = [
        x for x in filenames if x.startswith("wotis.") and x.endswith(".rpt")
    ]
    dates = [
        datetime.strptime("".join(filename.split(".")[2:4]), "%Y%j%H%M%S")
        for filename in filenames
    ]
    filedates = dict(zip(dates, filenames))

    dumps = []

    for date in sorted(dates):
        lines = []
        if not os.path.exists(os.path.join(gettempdir(), filedates[date])):
            try:
                f.prot_p()  # explicitly call for protected transfer
                f.retrlines("RETR " + os.path.join(url.path, filedates[date]),
                            lines.append)
            except ftplib.error_perm:
                logger.info("Permission error (???) on ftp server, skipping.")
                continue
            with open(os.path.join(gettempdir(), filedates[date]), "w") as fd_:
                for line in lines:
                    fd_.write(line + "\n")

        else:
            with open(os.path.join(gettempdir(), filedates[date]), "r") as fd_:
                for line in fd_:
                    lines.append(line)

        # for line in lines[7::2]:
        #     if line.strip() == '':
        #         break
        #     station, aos, elev, los = line.split()[:4]
        #     aos = datetime.strptime(aos, "%Y:%j:%H:%M:%S")
        #     los = datetime.strptime(los, "%Y:%j:%H:%M:%S")
        #     if los >= start_time and aos <= end_time:
        #         uptime = aos + (los - aos) / 2
        #         overpass = Pass(sat, aos, los, orb=satorb, uptime=uptime, instrument="modis")
        #         overpass.station = station
        #         overpass.max_elev = elev
        #         dumps.append(overpass)

        for line in lines[7::2]:
            if line.strip() == "":
                break
            station, aos, elev, los = line.split()[:4]
            aos = datetime.strptime(aos, "%Y:%j:%H:%M:%S")
            los = datetime.strptime(los, "%Y:%j:%H:%M:%S")
            dumps.append({"station": station, "aos": aos, "los": los, "elev": elev})

    if f is not None:
        f.quit()
    return dumps


def get_next_passes(satellites,
                    utctime,
                    forward,
                    coords,
                    tle_file=None,
                    aqua_terra_dumps=None,
                    min_pass=MIN_PASS,
                    local_horizon=0):
    """Get the next passes for *satellites*.

    Get the next passes for *satellites* , starting at *utctime*, for a
    duration of *forward* hours, with observer at *coords* ie lon (°E), lat
    (°N), altitude (km). Uses *tle_file* if provided, downloads from celestrack
    otherwise.

    Metop-A, Terra and Aqua need special treatment due to downlink restrictions.
    """
    passes = {}

    if tle_file is None and "TLES" not in os.environ:
        fp_, tle_file = mkstemp(prefix="tle", dir=gettempdir())
        os.close(fp_)
        logger.info("Fetch tle info from internet")
        tlefile.fetch(tle_file)

    if not os.path.exists(tle_file) and "TLES" not in os.environ:
        logger.info("Fetch tle info from internet")
        tlefile.fetch(tle_file)

    for sat in satellites:
        if not hasattr(sat, "name"):
            from trollsched.schedule import Satellite
            sat = Satellite(sat, 0, 0)

        satorb = orbital.Orbital(sat.name, tle_file=tle_file)
        passlist = satorb.get_next_passes(utctime,
                                          forward,
                                          *coords,
                                          horizon=local_horizon,
                                          )

        if sat.name.lower() == "metop-a":
            # Take care of metop-a special case
            passes["metop-a"] = get_metopa_passes(sat, passlist, satorb)
        elif sat.name.lower() in ["aqua", "terra"] and aqua_terra_dumps:
            # Take care of aqua (dumps in svalbard and poker flat)
            # Get the Terra/Aqua passes and fill the passes dict:
            get_terra_aqua_passes(passes, utctime, forward, sat, passlist, satorb, aqua_terra_dumps)
        else:
            if sat.name.upper() in VIIRS_PLATFORM_NAMES:
                instrument = "viirs"
            elif sat.name.lower().startswith("metop") or sat.name.lower().startswith("noaa"):
                instrument = "avhrr"
            elif sat.name.lower() in ["aqua", "terra"]:  # when aqua_terra_dumps=False
                instrument = "modis"
            elif sat.name.upper() in MERSI_PLATFORM_NAMES:
                instrument = "mersi"
            elif sat.name.upper() in MERSI2_PLATFORM_NAMES:
                instrument = "mersi-2"
            else:
                instrument = "unknown"

            passes[sat.name] = [
                Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument=instrument)
                for rtime, ftime, uptime in passlist
                if ftime - rtime > timedelta(minutes=min_pass)
            ]

    return set(fctools_reduce(operator.concat, list(passes.values())))


def get_metopa_passes(sat, passlist, satorb):
    """Get the Metop-A passes, taking care that Metop-A doesn't transmit to ground everywhere."""
    metop_passes = [
        Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument="avhrr")
        for rtime, ftime, uptime in passlist if rtime < ftime
    ]

    passes = []
    for overpass in metop_passes:
        if overpass.pass_direction() == "descending":
            new_rise = overpass.slsearch(60)
            if new_rise is not None and new_rise < overpass.falltime:
                overpass.risetime = new_rise
                # overpass has a boundary property, and it is not really needed here anyways!
                # overpass.boundary = SwathBoundary(overpass)
                if overpass.seconds() > MIN_PASS * 60:
                    passes.append(overpass)

    return passes


def get_terra_aqua_passes(passes, utctime, forward, sat, passlist, satorb, aqua_terra_dumps):
    """Get the Terra/Aqua passes.

    We take care that Terra and Aqua do not have direct broadcast when there are global dumps.

    Args:
        passes: The dictionary of satellite passes which is being built

        utctime: The start time (datetime object)

        forward: The number of hours ahead for which we will get the coming passes

        sat: The Satellite platform considered

        passlist: List of Pass objects

        satorb: Orbital instance for the actual satellite and tles considered

        aqua_terra_dumps: True or False or the actual URL to get info on Terra/Aqua
           dumps.  If True, the default URL will be used. If False or None, no dump
           info will be considered.

    """
    instrument = "modis"

    wpcoords = (-75.457222, 37.938611, 0)
    passlist_wp = satorb.get_next_passes(
        utctime - timedelta(minutes=30), forward + 1, *wpcoords)
    wp_passes = [
        Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument=instrument)
        for rtime, ftime, uptime in passlist_wp if rtime < ftime
    ]

    svcoords = (15.399, 78.228, 0)
    passlist_sv = satorb.get_next_passes(
        utctime - timedelta(minutes=30), forward + 1, *svcoords)
    sv_passes = [
        Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument=instrument)
        for rtime, ftime, uptime in passlist_sv if rtime < ftime
    ]
    pfcoords = (-147.43, 65.12, 0.51)
    passlist_pf = satorb.get_next_passes(
        utctime - timedelta(minutes=30), forward + 1, *pfcoords)
    pf_passes = [
        Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument=instrument)
        for rtime, ftime, uptime in passlist_pf if rtime < ftime
    ]

    aqua_passes = [
        Pass(sat, rtime, ftime, orb=satorb, uptime=uptime, instrument=instrument)
        for rtime, ftime, uptime in passlist if rtime < ftime
    ]

    dumps = get_aqua_terra_dumps(utctime - timedelta(minutes=30),
                                 utctime + timedelta(hours=forward + 0.5),
                                 satorb, sat, aqua_terra_dumps)

    # remove the known dumps
    for dump in dumps:
        # print "*", dump.station, dump, dump.max_elev
        logger.debug("dump from ftp: " + str((dump.station, dump,
                                              dump.max_elev)))
        for i, sv_pass in enumerate(sv_passes):
            if sv_pass.overlaps(dump, timedelta(minutes=40)):
                sv_elevation = sv_pass.orb.get_observer_look(
                    sv_pass.uptime, *svcoords)[1]
                logger.debug("Computed " + str(("SG", sv_pass,
                                                sv_elevation)))
                del sv_passes[i]
        for i, pf_pass in enumerate(pf_passes):
            if pf_pass.overlaps(dump, timedelta(minutes=40)):
                pf_elevation = pf_pass.orb.get_observer_look(
                    pf_pass.uptime, *pfcoords)[1]
                logger.debug("Computed " + str(("PF", pf_pass,
                                                pf_elevation)))
                del pf_passes[i]
        for i, wp_pass in enumerate(wp_passes):
            if wp_pass.overlaps(dump, timedelta(minutes=40)):
                wp_elevation = wp_pass.orb.get_observer_look(
                    wp_pass.uptime, *wpcoords)[1]
                logger.debug("Computed " + str(("WP", wp_pass,
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
                sv_elevation = sv_pass.orb.get_observer_look(
                    sv_pass.uptime, *svcoords)[1]
                pf_elevation = pf_pass.orb.get_observer_look(
                    pf_pass.uptime, *pfcoords)[1]
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

    passes[sat.name] = []
    for overpass in aqua_passes:
        add = True
        for dump_pass in dumps:
            if dump_pass.overlaps(overpass):
                if (dump_pass.uptime < overpass.uptime and
                        dump_pass.falltime > overpass.risetime):
                    logger.debug("adjusting " + str(overpass) +
                                 " to new risetime " +
                                 str(dump_pass.falltime))
                    overpass.risetime = dump_pass.falltime
                    overpass.boundary = SwathBoundary(overpass)
                elif (dump_pass.uptime >= overpass.uptime and
                      dump_pass.risetime < overpass.falltime):
                    logger.debug("adjusting " + str(overpass) +
                                 " to new falltime " +
                                 str(dump_pass.risetime))
                    overpass.falltime = dump_pass.risetime
                    overpass.boundary = SwathBoundary(overpass)
                if overpass.falltime <= overpass.risetime:
                    add = False
                    logger.debug("skipping " + str(overpass))
        if add and overpass.seconds() > MIN_PASS * 60:
            passes[sat.name].append(overpass)

    return
