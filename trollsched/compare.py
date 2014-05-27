#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

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

"""Compare the request file and the confirmation file.
"""

import logging
import logging.handlers
import sys
import os
import glob

logger = logging.getLogger(__name__)

def xml_compare(x1_, x2_, reporter=None, skiptags=None):
    """Compare xml objects.
    """
    if x1_.tag != x2_.tag:
        if reporter:
            reporter('Tags do not match: %s and %s' % (x1_.tag, x2_.tag))
        return False
    for name, value in x1_.attrib.items():
        if x2_.attrib.get(name) != value:
            if reporter:
                reporter('Attributes do not match: %s=%r, %s=%r'
                         % (name, value, name, x2_.attrib.get(name)))
            return False
    for name in x2_.attrib.keys():
        if name not in x1_.attrib:
            if reporter:
                reporter('x2_ has an attribute x1_ is missing: %s'
                         % name)
            return False
    if not text_compare(x1_.text, x2_.text):
        if reporter:
            reporter('text: %r != %r' % (x1_.text, x2_.text))
        return False
    if not text_compare(x1_.tail, x2_.tail):
        if reporter:
            reporter('tail: %r != %r' % (x1_.tail, x2_.tail))
        return False
    cl1 = x1_.getchildren()
    cl2 = x2_.getchildren()
    if len(cl1) != len(cl2):
        if reporter:
            reporter('not the same number of passes, %i != %i'
                     % (len(cl1), len(cl2)))
        return False
    i = 0
    for c1, c2 in zip(cl1, cl2):
        i += 1
        if skiptags and (c1.tag in skiptags):
            continue
        if not xml_compare(c1, c2, reporter=reporter):
            if reporter:
                reporter('element %i do not match: %s'
                         % (i, c1.tag))
            return False
    return True


def text_compare(t1_, t2_):
    """Compare text fields.
    """
    if not t1_ and not t2_:
        return True
    if t1_ == '*' or t2_ == '*':
        return True
    return (t1_ or '').strip() == (t2_ or '').strip()

def compare(file1, file2):
    """Compare two xml files, request and confirmation.
    """
    import xml.etree.ElementTree as ET
    xml1 = ET.parse(file1).getroot()
    xml2 = ET.parse(file2).getroot()
    if xml_compare(xml1, xml2, logger.error,
                   ["confirmed-by", "confirmed-on", "properties"]):
        logger.info("All passes confirmed.")
        return True
    else:
        return False


# import fnmatch
# import pyinotify

# class EventHandler(pyinotify.ProcessEvent):
#     """Manage events.
#     """
#     def process_IN_CLOSE_WRITE(self, event):
#         if not fnmatch.fnmatch(event.pathname, "*-acquisition-schedule-confirmation-???.xml"):
#             return
#         logger.info("Processing: %s", event.pathname)
#         reqname = event.pathname[:-20] + "request" +  event.pathname[-8:]
#         logger.info("Validating against: %s", reqname)
#         compare(reqname, event.pathname)
#         sys.exit(0)

#     def process_IN_MOVED_TO(self, event):
#         self.process_IN_CLOSE_WRITE(event)

def run():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="req/conf files to compare",
                        nargs=2)

    parser.add_argument("-m", "--mail", nargs="*",
                        help="mail address(es) to send error messages to.",
                        default=None)
    parser.add_argument("-v", "--verbose", help="activate debug messages",
                        action="store_true")
    parser.add_argument("-l", "--log", help="file to log to")
#    parser.add_argument("-w", "--watch",
#                        help="directory to watch for new confirmation files")
    parser.add_argument("-r", "--most-recent",
                        help="check the most recent request against the" +
                        " corresponding confirmation, from the given directory")
    parser.add_argument("-c", "--confirmation",
                        help="directory for the confirmation files")
    
    opts = parser.parse_args()

    if opts.log:
        handler = logging.handlers.TimedRotatingFileHandler(opts.log,
                                                            "midnight",
                                                            backupCount=7)
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
                                                "satsateknik@smhi.se",
                                                opts.mail,
                                                "Scheduler confirmation")
        mhandler.setLevel(logging.WARNING)
        logging.getLogger('').addHandler(mhandler)
    logger = logging.getLogger("compare")

    logger.debug("DEBUG on")

    if opts.file:
        compare(opts.file[0], opts.file[1])

    # if opts.watch:
    #     wm = pyinotify.WatchManager() # Watch Manager
    #     mask = pyinotify.IN_CLOSE_WRITE | pyinotify.IN_MOVED_TO
    #     handler = EventHandler()
    #     notifier = pyinotify.Notifier(wm, handler)
    #     wdd = wm.add_watch(opts.watch, mask, rec=False)

    #     notifier.loop()



    if opts.most_recent:
        logger.debug("looking for most recent file in " +
                     os.path.join(opts.most_recent, "*request*.xml"))
        filelist = glob.glob(os.path.join(opts.most_recent, "*request*.xml"))
        newest = max(filelist, key=lambda x: os.stat(x).st_mtime)
        logger.debug("checking " + newest)
        reqdir, newfile = os.path.split(newest)
        confdir = opts.confirmation or reqdir
        confname = os.path.join(confdir,
                                newfile[:-15] + "confirmation" +  newfile[-8:])
        logger.debug("against " + confname)
        try:
            compare(newest, confname)
        except IOError:
            logger.exception("Something went wrong!") 

if __name__ == '__main__':
    run()
