.. pytroll-schedule documentation master file, created by
   sphinx-quickstart on Mon Feb 24 23:43:03 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pytroll-schedule's documentation!
============================================

Usage of the schedule script::

  usage: schedule [-h] [--lon LON] [--lat LAT] [--alt ALT] [-l LOG] [-v]
                  [-t TLE] [-f FORWARD] [-s START_TIME] [-d DELAY] [-c CONFIG]
                  [-x XML] [--scisys SCISYS]

  optional arguments:
    -h, --help            show this help message and exit
    --lon LON             Longitude, degrees east
    --lat LAT             Latitude, degrees north
    --alt ALT             Altitude, km
    -l LOG, --log LOG     File to log to (defaults to stdout)
    -v, --verbose         print debug messages too
    -t TLE, --tle TLE     tle file to use
    -f FORWARD, --forward FORWARD
                          time ahead to compute the schedule
    -s START_TIME, --start-time START_TIME
                          start time of the schedule to compute
    -d DELAY, --delay DELAY
                          delay (in seconds) needed between two consecutive
                          passes (60 seconds by default)
    -c CONFIG, --config CONFIG
                          configuration file to use

  output:
    -x XML, --xml XML     generate an xml request file and put it in this
                          directory
    --scisys SCISYS       path to the schedule file


Contents:

.. toctree::
   :maxdepth: 2

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

