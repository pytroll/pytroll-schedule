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

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

import numpy as np
from datetime import datetime, timedelta

from trollsched.schedule import fermia, fermib, conflicting_passes
from pyresample.boundary import AreaBoundary
from pyorbital import orbital
from trollsched.satpass import get_next_passes
from trollsched.satpass import Pass


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

    def setUp(self):
        """Set up"""
        from pyorbital import orbital

        self.utctime = datetime(2018, 11, 28, 10, 0)
        self.satellites = ["noaa-20", ]
        self.tles = {'noaa-20': {}}
        self.tles['noaa-20']['line1'] = "1 43013U 17073A   18331.00000000  .00000048  00000-0  22749-4 0  3056"
        self.tles['noaa-20']['line2'] = "2 43013 098.7413 267.0121 0001419 108.5818 058.1314 14.19552981053016"

        self.aquas = ["aqua", ]
        self.tles['aqua'] = {}
        self.tles['aqua']['line1'] = "1 27424U 02022A   18332.21220389  .00000093  00000-0  30754-4 0  9994"
        self.tles['aqua']['line2'] = "2 27424  98.2121 270.9368 0001045 343.9225 155.8703 14.57111538881313"

        self.orb = orbital.Orbital('NOAA 20',
                                   line1=self.tles['noaa-20']['line1'],
                                   line2=self.tles['noaa-20']['line2'])
        self.aqua_orb = orbital.Orbital('AQUA',
                                        line1=self.tles['aqua']['line1'],
                                        line2=self.tles['aqua']['line2'])

        self.dumpdata = [
            {'los': datetime(2018, 11, 28, 10, 0, 30), 'station': 'USAK05',
             'aos': datetime(2018, 11, 28, 9, 50, 24), 'elev': '11.188'},
            {'los': datetime(2018, 11, 28, 11, 39, 47), 'station': 'AS2',
             'aos': datetime(2018, 11, 28, 11, 28, 51), 'elev': '39.235'},
            {'los': datetime(2018, 11, 28, 13, 19, 8), 'station': 'USAK05',
             'aos': datetime(2018, 11, 28, 13, 6, 36), 'elev': '58.249'},
            {'los': datetime(2018, 11, 28, 14, 54, 25), 'station': 'AS2',
             'aos': datetime(2018, 11, 28, 14, 44, 37), 'elev': '22.403'},
            {'los': datetime(2018, 11, 28, 16, 27, 22), 'station': 'SG1',
             'aos': datetime(2018, 11, 28, 16, 16, 58), 'elev': '9.521'}
        ]

    @patch('os.path.exists')
    def test_get_next_passes_viirs(self, exists):

        exists.return_code = True

        # mymock:
        with patch('pyorbital.orbital.Orbital') as mymock:
            instance = mymock.return_value
            instance.get_next_passes = self.orb.get_next_passes

            allpasses = get_next_passes(self.satellites, self.utctime,
                                        4, (16, 58, 0), tle_file='nonexisting')

            self.assertEqual(len(allpasses), 2)

            n20pass1 = allpasses.pop()

            rt1 = datetime(2018, 11, 28, 10, 53, 42, 79483)
            ft1 = datetime(2018, 11, 28, 11, 9, 6, 916787)
            rt2 = datetime(2018, 11, 28, 12, 34, 44, 667963)
            ft2 = datetime(2018, 11, 28, 12, 49, 25, 134067)

            dt_ = n20pass1.risetime - rt1
            self.assertAlmostEqual(dt_.seconds, 0)

            dt_ = n20pass1.falltime - ft1
            self.assertAlmostEqual(dt_.seconds, 0)

            n20pass2 = allpasses.pop()

            dt_ = n20pass2.risetime - rt2
            self.assertAlmostEqual(dt_.seconds, 0)

            dt_ = n20pass2.falltime - ft2
            self.assertAlmostEqual(dt_.seconds, 0)

    @patch('os.path.exists')
    @patch('trollsched.satpass.get_aqua_terra_dumpdata_from_ftp')
    def test_get_next_passes_with_aquadumps(self, dumps_from_ftp, exists):
        dumps_from_ftp.return_value = self.dumpdata
        exists.return_code = True

        # mymock:
        with patch('pyorbital.orbital.Orbital') as mymock:
            instance = mymock.return_value
            instance.get_next_passes = self.aqua_orb.get_next_passes

            allpasses = get_next_passes(self.aquas, self.utctime,
                                        6, (16, 58, 0), tle_file='nonexisting',
                                        aqua_terra_dumps=True)

            self.assertEqual(len(allpasses), 3)

            rt1 = datetime(2018, 11, 28, 11, 12, 8, 728455)
            ft1 = datetime(2018, 11, 28, 11, 26, 8, 250021)
            rt2 = datetime(2018, 11, 28, 12, 50, 46, 574975)
            ft2 = datetime(2018, 11, 28, 13, 3, 53, 262440)
            rt3 = datetime(2018, 11, 28, 14, 33, 33, 973194)
            ft3 = datetime(2018, 11, 28, 14, 40, 10, 761405)

            for mypass in allpasses:
                dtmin = timedelta(seconds=10000000)
                for risetime in [rt1, rt2, rt3]:
                    dt_ = abs(mypass.risetime - risetime)
                    if dt_ < dtmin:
                        dtmin = dt_

                self.assertAlmostEqual(dtmin.seconds, 0)

                dtmin = timedelta(seconds=10000000)
                for falltime in [ft1, ft2, ft3]:
                    dt_ = abs(mypass.falltime - falltime)
                    if dt_ < dtmin:
                        dtmin = dt_

                self.assertAlmostEqual(dtmin.seconds, 0)

    def tearDown(self):
        """Clean up"""
        pass


def suite():
    """The suite for test_schedule
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestUtils))
    mysuite.addTest(loader.loadTestsFromTestCase(TestAreaBoundary))
    mysuite.addTest(loader.loadTestsFromTestCase(TestTools))

    return mysuite
