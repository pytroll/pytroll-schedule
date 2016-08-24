#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Adam.Dybbroe

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

"""From a schedule request xml file generate png's with swath coverage outline
and an xml page for visualisation. It uses posttroll to listen for incoming
schedule request xml files and then triggers the png and xml output generation.

"""

import os
from ConfigParser import RawConfigParser
import logging
LOG = logging.getLogger(__name__)

CFG_DIR = os.environ.get('PYTROLL_SCHEDULE_CONFIG_DIR', './')
CONF = RawConfigParser()
CFG_FILE = os.path.join(CFG_DIR, "pytroll_schedule_config.cfg")
LOG.debug("Config file = " + str(CFG_FILE))
if not os.path.exists(CFG_FILE):
    raise IOError('Config file %s does not exist!' % CFG_FILE)

CONF.read(CFG_FILE)
OPTIONS = {}
for option, value in CONF.items("DEFAULT"):
    OPTIONS[option] = value


#: Default time format
_DEFAULT_TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

#: Default log format
_DEFAULT_LOG_FORMAT = '[%(levelname)s: %(asctime)s : %(name)s] %(message)s'

import sys
from urlparse import urlparse
import posttroll.subscriber
from posttroll.publisher import Publish
import xml.etree.ElementTree as ET
from datetime import datetime
import os.path

from trollsched.satpass import Pass

sat_dict = {'npp': 'Suomi NPP',
            'noaa19': 'NOAA 19',
            'noaa18': 'NOAA 18',
            'noaa15': 'NOAA 15',
            'aqua': 'Aqua',
            'terra': 'Terra',
            'metop-b': 'Metop-B',
            'metop-a': 'Metop-A',
            }


def process_xmlrequest(filename, plotdir, output_file):

    tree = ET.parse(filename)
    root = tree.getroot()

    for child in root:
        if child.tag == 'pass':
            print child.attrib
            overpass = Pass(sat_dict.get(child.attrib['satellite'],
                                         child.attrib['satellite']),
                            datetime.strptime(child.attrib['start-time'],
                                              '%Y-%m-%d-%H:%M:%S'),
                            datetime.strptime(child.attrib['end-time'],
                                              '%Y-%m-%d-%H:%M:%S'))

            overpass.save_fig(directory=plotdir)
            child.set('img', overpass.fig)
            child.set('rec', 'True')

    tree.write(output_file, encoding='utf-8', xml_declaration=True)

    with open(output_file) as fpt:
        lines = fpt.readlines()
        lines.insert(
            1, "<?xml-stylesheet type='text/xsl' href='reqreader.xsl'?>")

    with open(output_file, 'w') as fpt:
        fpt.writelines(lines)


def start_plotting(jobreg, message, **kwargs):
    """Read the xmlschedule request file and make the png images of swath outlines
    and generate the output xml file for web publication

    """
    LOG.info("")
    LOG.info("job-registry dict: " + str(jobreg))
    LOG.info("\tMessage:")
    LOG.info(message)
    urlobj = urlparse(message.data['uri'])
    # path, fname = os.path.split(urlobj.path)

    process_xmlrequest(urlobj.path,
                       OPTIONS['path_plots'], OPTIONS['xmlfilepath'])

    return jobreg


def schedule_page_generator():
    """Listens and triggers processing"""

    LOG.info(
        "*** Start the generation of the schedule xml page with swath outline plots")
    with posttroll.subscriber.Subscribe('', [OPTIONS['posttroll_topic'], ],
                                        True) as subscr:
        with Publish('schedule_page_generator', 0) as publisher:
            job_registry = {}
            for msg in subscr.recv():
                job_registry = start_plotting(
                    job_registry, msg, publisher=publisher)
                # Cleanup in registry (keep only the last 5):
                keys = job_registry.keys()
                if len(keys) > 5:
                    keys.sort()
                    job_registry.pop(keys[0])

if __name__ == "__main__":

    handler = logging.StreamHandler(sys.stderr)

    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(fmt=_DEFAULT_LOG_FORMAT,
                                  datefmt=_DEFAULT_TIME_FORMAT)
    handler.setFormatter(formatter)
    logging.getLogger('').addHandler(handler)
    logging.getLogger('').setLevel(logging.DEBUG)
    logging.getLogger('posttroll').setLevel(logging.INFO)

    LOG = logging.getLogger('schedule_page_generator')

    schedule_page_generator()
