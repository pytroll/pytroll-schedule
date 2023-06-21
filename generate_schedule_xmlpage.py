#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016, 2018, 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c20671.ad.smhi.se>

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

"""From a schedule request xml file generate png's with swath coverage outline and an xml page for visualisation.

It uses posttroll to listen for incoming schedule request xml files and then triggers the png and xml output generation.
"""

import logging
import os.path
import sys
from datetime import datetime

import defusedxml.ElementTree as ET
import posttroll.subscriber
from posttroll.publisher import Publish
from six.moves.configparser import RawConfigParser
from six.moves.urllib.parse import urlparse

from trollsched import INSTRUMENT, SATELLITE_NAMES
from trollsched.drawing import save_fig
from trollsched.satpass import Pass

LOG = logging.getLogger(__name__)

CFG_DIR = os.environ.get("PYTROLL_SCHEDULE_CONFIG_DIR", "./")
CONF = RawConfigParser()
CFG_FILE = os.path.join(CFG_DIR, "pytroll_schedule_config.cfg")
LOG.debug("Config file = " + str(CFG_FILE))
if not os.path.exists(CFG_FILE):
    raise IOError("Config file %s does not exist!" % CFG_FILE)

CONF.read(CFG_FILE)
OPTIONS = {}
for option, value in CONF.items("DEFAULT"):
    OPTIONS[option] = value


#: Default time format
_DEFAULT_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

#: Default log format
_DEFAULT_LOG_FORMAT = "[%(levelname)s: %(asctime)s : %(name)s] %(message)s"


def process_xmlrequest(filename, plotdir, output_file, excluded_satellites):
    """Process the xml request."""
    tree = ET.parse(filename)
    root = tree.getroot()

    for child in root:
        if child.tag == "pass":
            LOG.debug("Pass: %s", str(child.attrib))
            platform_name = SATELLITE_NAMES.get(child.attrib["satellite"], child.attrib["satellite"])
            instrument = INSTRUMENT.get(platform_name)
            if not instrument:
                LOG.error("Instrument unknown! Platform = %s", platform_name)
                continue

            if platform_name in excluded_satellites:
                LOG.debug("Platform name excluded: %s", platform_name)
                continue
            try:
                overpass = Pass(platform_name,
                                datetime.strptime(child.attrib["start-time"],
                                                  "%Y-%m-%d-%H:%M:%S"),
                                datetime.strptime(child.attrib["end-time"],
                                                  "%Y-%m-%d-%H:%M:%S"),
                                instrument=instrument)
            except KeyError as err:
                LOG.warning("Failed on satellite %s: %s", platform_name, str(err))
                continue

            save_fig(overpass, directory=plotdir)
            child.set("img", overpass.fig)
            child.set("rec", "True")
            LOG.debug("Plot saved - plotdir = %s, platform_name = %s", plotdir, platform_name)

    tree.write(output_file, encoding="utf-8", xml_declaration=True)

    with open(output_file) as fpt:
        lines = fpt.readlines()
        lines.insert(
            1, "<?xml-stylesheet type='text/xsl' href='reqreader.xsl'?>")

    with open(output_file, "w") as fpt:
        fpt.writelines(lines)


def start_plotting(jobreg, message, **kwargs):
    """Make a web-publishable version of the xml schedule.

    Read the xmlschedule request file and make the png images of swath outlines
    and generate the output xml file for web publication.
    """
    excluded_satellites = kwargs.get("excluded_satellites", [])

    LOG.info("")
    LOG.info("job-registry dict: " + str(jobreg))
    LOG.info("\tMessage:")
    LOG.info(message)
    urlobj = urlparse(message.data["uri"])
    # path, fname = os.path.split(urlobj.path)

    process_xmlrequest(urlobj.path,
                       OPTIONS["path_plots"], OPTIONS["xmlfilepath"],
                       excluded_satellites)

    return jobreg


def schedule_page_generator(excluded_satellite_list=None):
    """Listens and triggers processing."""
    LOG.info(
        "*** Start the generation of the schedule xml page with swath outline plots")
    with posttroll.subscriber.Subscribe("", [OPTIONS["posttroll_topic"], ],
                                        True) as subscr:
        with Publish("schedule_page_generator", 0) as publisher:
            job_registry = {}
            for msg in subscr.recv():
                job_registry = start_plotting(
                    job_registry, msg, publisher=publisher, excluded_satellites=excluded_satellite_list)
                # Cleanup in registry (keep only the last 5):
                keys = job_registry.keys()
                if len(keys) > 5:
                    keys.sort()
                    job_registry.pop(keys[0])


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--excluded_satellites", nargs="*",
                        help="List of platform names to exclude",
                        default=[])
    opts = parser.parse_args()

    no_sats = opts.excluded_satellites

    handler = logging.StreamHandler(sys.stderr)

    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(fmt=_DEFAULT_LOG_FORMAT,
                                  datefmt=_DEFAULT_TIME_FORMAT)
    handler.setFormatter(formatter)
    logging.getLogger("").addHandler(handler)
    logging.getLogger("").setLevel(logging.DEBUG)
    logging.getLogger("posttroll").setLevel(logging.INFO)

    LOG = logging.getLogger("schedule_page_generator")
    LOG.info("Exclude the following satellite platforms: %s", str(no_sats))

    schedule_page_generator(no_sats)

    # uri = "/data/temp/AdamD/xxx/2018-10-22-00-42-28-acquisition-schedule-confirmation-nrk.xml"
    # urlobj = urlparse(uri)
    # process_xmlrequest(urlobj.path,
    #                    OPTIONS['path_plots'], OPTIONS['xmlfilepath'], no_sats)
