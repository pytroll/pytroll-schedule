#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2018 Martin Raspaud

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

"""Test the schedule module.
"""

import unittest
import numpy as np
from datetime import datetime, timedelta

from trollsched.schedule import fermia, fermib, conflicting_passes
from pyresample.boundary import AreaBoundary

# class TestPass(unittest.TestCase):

# def test_day(self):
#     satellite = "noaa 16"
#     tle1 = "1 26536U 00055A   13076.42963155  .00000201  00000-0  13237-3 0  1369"
#     tle2 = "2 26536  99.0540 128.2392 0010826  39.9070  85.2960 14.12848373643614"
#     orb = Orbital(satellite, line1=tle1, line2=tle2)
#     tstart = datetime(2013, 3, 18, 8, 15, 22, 352000)
#     tup = datetime(2013, 3, 18, 8, 22, 52, 352000)
#     tend = datetime(2013, 3, 18, 8, 30, 22, 352000)
#     overp = Pass(satellite, tstart, tend, orb, tup)

# a little night

#     day = overp.day()

#     self.assertEquals(0.99735685408290298, day)

# on the area of interest there is no night

#     area_of_interest = get_area_def("euron1")
#     day = overp.day(area_of_interest)
#     self.assertEquals(1.0, day)

#     tstart = datetime(2013, 3, 18, 8, 16, 22, 352000)
#     overp = Pass(satellite, tstart, tend, orb, tup)

# an entire pass without night

#     day = overp.day()
#     self.assertEquals(1.0, day)


class TestTools(unittest.TestCase):

    def test_conflicting_passes(self):

        class MyPass(object):

            def __init__(self, rise, fall):
                self.risetime = rise
                self.falltime = fall

        ref_time = datetime.utcnow()
        passes = [MyPass(ref_time, ref_time + timedelta(minutes=10)),
                  MyPass(ref_time + timedelta(minutes=10.01),
                         ref_time + timedelta(minutes=20))]
        self.assertEquals(
            len(conflicting_passes(passes, timedelta(seconds=0))), 2)
        self.assertEquals(
            len(conflicting_passes(passes, timedelta(seconds=60))), 1)


class TestAreaBoundary(unittest.TestCase):

    def test_contour(self):

        side1_lons = np.arange(4)
        side1_lats = np.arange(4) + 20
        side2_lons = np.arange(4) + 3
        side2_lats = np.arange(4) + 20 + 3
        side3_lons = np.arange(4) + 6
        side3_lats = np.arange(4) + 20 + 6
        side4_lons = np.arange(4) + 9
        side4_lats = np.arange(4) + 20 + 9

        bond = AreaBoundary((side1_lons, side1_lats),
                            (side2_lons, side2_lats),
                            (side3_lons, side3_lats),
                            (side4_lons, side4_lats))

        lons, lats = bond.contour()
        self.assertTrue(np.allclose(lons, np.arange(12)))
        self.assertTrue(np.allclose(lats, np.arange(12) + 20))

    def test_decimate(self):

        side1_lons = np.arange(8)
        side1_lats = np.arange(8) + 30
        side2_lons = np.arange(8) + 7
        side2_lats = np.arange(8) + 30 + 7
        side3_lons = np.arange(8) + 14
        side3_lats = np.arange(8) + 30 + 14
        side4_lons = np.arange(8) + 21
        side4_lats = np.arange(8) + 30 + 21

        bond = AreaBoundary((side1_lons, side1_lats),
                            (side2_lons, side2_lats),
                            (side3_lons, side3_lats),
                            (side4_lons, side4_lats))

        bond.decimate(5)
        lons, lats = bond.contour()

        self.assertTrue(np.allclose(lons,
                                    np.array([0, 1, 6, 7, 8,
                                              13, 14, 15, 20, 21, 22, 27])))
        self.assertTrue(np.allclose(lats,
                                    np.array([30, 31, 36, 37, 38, 43, 44, 45,
                                              50, 51, 52, 57])))


class TestUtils(unittest.TestCase):

    def test_fermi(self):
        self.assertEquals(fermia(0.25), 0.875)
        self.assertEquals(fermib(0.25), 0.5)


class TestAll(unittest.TestCase):

    def test_gen(self):
        utctime = datetime(2014, 1, 18, 14, 20)
        satellites = ["noaa 19", "noaa 18", "noaa 16", "noaa 15",
                      "metop-a", "metop-b",
                      "terra", "aqua",
                      "suomi npp"]
        tle_file = "./tle_20140120.txt"
        #allpasses = get_next_passes(satellites, utctime, 2, tle_file)


def suite():
    """The suite for test_schedule
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestUtils))
    mysuite.addTest(loader.loadTestsFromTestCase(TestAreaBoundary))
    mysuite.addTest(loader.loadTestsFromTestCase(TestTools))

    return mysuite
