#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 - 2019 PyTroll

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
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

"""Test the schedule module.
"""

import numpy as np
from datetime import datetime, timedelta

from trollsched.schedule import fermia, fermib, conflicting_passes
from trollsched.schedule import parse_datetime, build_filename
from pyresample.boundary import AreaBoundary
from trollsched.satpass import get_next_passes
from trollsched.satpass import get_aqua_terra_dumps
from trollsched.satpass import get_metopa_passes

import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch


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

    def test_parse_datetime(self):

        dtobj = parse_datetime('20190104110059')
        self.assertEqual(dtobj, datetime(2019, 1, 4, 11, 0, 59))

    def test_build_filename(self):

        pattern_name = "dir_output"
        pattern_dict = {'file_xml': '{dir_output}/{date}-{time}-aquisition-schedule-{mode}-{station}.xml', 'file_sci': '{dir_output}/scisys-schedule-{station}.txt',
                        'dir_plots': '{dir_output}/plots.{station}', 'dir_output': '/tmp', 'file_graph': '{dir_output}/graph.{station}'}
        kwargs = {'date': '20190104', 'output_dir': '.', 'dir_output': '/tmp', 'time': '122023'}

        res = build_filename(pattern_name, pattern_dict, kwargs)
        self.assertEqual(res, '/tmp')

        pattern_name = "file_xml"
        kwargs = {'station': 'nrk', 'mode': 'request', 'time': '125334',
                  'date': '20190104', 'dir_output': '/tmp', 'output_dir': '.'}
        res = build_filename(pattern_name, pattern_dict, kwargs)
        self.assertEqual(res, '/tmp/20190104-125334-aquisition-schedule-request-nrk.xml')


class TestAll(unittest.TestCase):

    def setUp(self):
        """Set up"""
        from pyorbital import orbital
        from trollsched.schedule import Satellite

        self.utctime = datetime(2018, 11, 28, 10, 0)
        self.satellites = ["noaa-20", ]
        self.tles = {'noaa-20': {}}
        self.tles['noaa-20']['line1'] = "1 43013U 17073A   18331.00000000  .00000048  00000-0  22749-4 0  3056"
        self.tles['noaa-20']['line2'] = "2 43013 098.7413 267.0121 0001419 108.5818 058.1314 14.19552981053016"

        self.aquas = ["aqua", ]
        self.terras = ["terra", ]
        self.terra = Satellite('terra', 0, 0)
        self.metopa = Satellite('metop-a', 0, 0)

        self.tles['aqua'] = {}
        self.tles['aqua']['line1'] = "1 27424U 02022A   18332.21220389  .00000093  00000-0  30754-4 0  9994"
        self.tles['aqua']['line2'] = "2 27424  98.2121 270.9368 0001045 343.9225 155.8703 14.57111538881313"
        self.tles['terra'] = {}
        self.tles['terra']['line1'] = "1 25994U 99068A   18338.20920286  .00000076  00000-0  26867-4 0  9999"
        self.tles['terra']['line2'] = "2 25994  98.2142  50.5750 0000577 102.5211 257.6060 14.57132862  8586"
        self.tles['metop-a'] = {}
        self.tles['metop-a']['line1'] = "1 29499U 06044A   18338.30873671  .00000000  00000+0  31223-4 0 00013"
        self.tles['metop-a']['line2'] = "2 29499  98.6045  31.7725 0001942  91.8780 346.4884 14.21536046629175"

        self.orb = orbital.Orbital('NOAA 20',
                                   line1=self.tles['noaa-20']['line1'],
                                   line2=self.tles['noaa-20']['line2'])
        self.aqua_orb = orbital.Orbital('AQUA',
                                        line1=self.tles['aqua']['line1'],
                                        line2=self.tles['aqua']['line2'])
        self.terra_orb = orbital.Orbital('TERRA',
                                         line1=self.tles['terra']['line1'],
                                         line2=self.tles['terra']['line2'])
        self.metopa_orb = orbital.Orbital('Metop-A',
                                          line1=self.tles['metop-a']['line1'],
                                          line2=self.tles['metop-a']['line2'])

        # These values were used to generate the get_next_passes list mock:
        # utctime = datetime(2018, 12, 4, 9, 0)
        # forward = 6
        # coords = (16, 58, 0)
        self.metopa_passlist = [(datetime(2018, 12, 4, 9, 10, 4, 574801),
                                 datetime(2018, 12, 4, 9, 25, 29, 157194),
                                 datetime(2018, 12, 4, 9, 17, 48, 530484)),
                                (datetime(2018, 12, 4, 10, 50, 23, 899232),
                                 datetime(2018, 12, 4, 11, 4, 2, 335184),
                                 datetime(2018, 12, 4, 10, 57, 13, 691637)),
                                (datetime(2018, 12, 4, 12, 30, 24, 97160),
                                 datetime(2018, 12, 4, 12, 40, 42, 403698),
                                 datetime(2018, 12, 4, 12, 35, 33, 317647)),
                                (datetime(2018, 12, 4, 14, 9, 1, 937869),
                                 datetime(2018, 12, 4, 14, 17, 20, 556654),
                                 datetime(2018, 12, 4, 14, 13, 11, 247497))]

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
        self.dumpdata_terra = [{'los': datetime(2018, 11, 20, 23, 24, 41), 'station': 'SG2',
                                'aos': datetime(2018, 11, 20, 23, 12, 32), 'elev': '17.4526'},
                               {'los': datetime(2018, 11, 22, 23, 19, 21), 'station': 'AS3',
                                'aos': datetime(2018, 11, 22, 23, 8, 55), 'elev': '28.9558'},
                               {'los': datetime(2018, 11, 22, 23, 19, 21), 'station': 'AS3',
                                'aos': datetime(2018, 11, 22, 23, 8, 55), 'elev': '28.9558'},
                               {'los': datetime(2018, 11, 26, 22, 47, 34), 'station': 'SG1',
                                'aos': datetime(2018, 11, 26, 22, 34, 58), 'elev': '21.5694'},
                               {'los': datetime(2018, 11, 26, 22, 47, 34), 'station': 'SG1',
                                'aos': datetime(2018, 11, 26, 22, 34, 58), 'elev': '21.5694'},
                               {'los': datetime(2018, 11, 26, 22, 47, 34), 'station': 'SG1',
                                'aos': datetime(2018, 11, 26, 22, 34, 58), 'elev': '21.5694'},
                               {'los': datetime(2018, 11, 27, 23, 30, 44), 'station': 'SG2',
                                'aos': datetime(2018, 11, 27, 23, 18, 39), 'elev': '16.8795'},
                               {'los': datetime(2018, 11, 27, 23, 30, 44), 'station': 'SG2',
                                'aos': datetime(2018, 11, 27, 23, 18, 39), 'elev': '16.8795'},
                               {'los': datetime(2018, 11, 28, 22, 43, 53), 'station': 'USAK05',
                                'aos': datetime(2018, 11, 28, 22, 31, 57), 'elev': '40.9264'},
                               {'los': datetime(2018, 11, 28, 22, 43, 53), 'station': 'USAK05',
                                'aos': datetime(2018, 11, 28, 22, 31, 57), 'elev': '40.9264'},
                               {'los': datetime(2018, 11, 29, 23, 25, 11), 'station': 'USAK05',
                                'aos': datetime(2018, 11, 29, 23, 14, 47), 'elev': '26.9937'},
                               {'los': datetime(2018, 11, 29, 23, 25, 11), 'station': 'USAK05',
                                'aos': datetime(2018, 11, 29, 23, 14, 47), 'elev': '26.9937'},
                               {'los': datetime(2018, 11, 30, 22, 31, 3), 'station': 'AS2',
                                'aos': datetime(2018, 11, 30, 22, 19, 48), 'elev': '47.8599'},
                               {'los': datetime(2018, 12, 1, 1, 29, 2), 'station': 'WG1',
                                'aos': datetime(2018, 12, 1, 1, 21, 11), 'elev': '8.0543'},
                               {'los': datetime(2018, 11, 30, 22, 31, 3), 'station': 'AS2',
                                'aos': datetime(2018, 11, 30, 22, 19, 48), 'elev': '47.8599'},
                               {'los': datetime(2018, 12, 1, 1, 29, 2), 'station': 'WG1',
                                'aos': datetime(2018, 12, 1, 1, 21, 11), 'elev': '8.0543'},
                               {'los': datetime(2018, 12, 3, 1, 28, 14), 'station': 'SG2',
                                'aos': datetime(2018, 12, 3, 1, 17, 53), 'elev': '9.2428'},
                               {'los': datetime(2018, 12, 3, 22, 53, 35), 'station': 'SG1',
                                'aos': datetime(2018, 12, 3, 22, 41, 5), 'elev': '20.8371'},
                               {'los': datetime(2018, 12, 3, 22, 53, 35), 'station': 'SG1',
                                'aos': datetime(2018, 12, 3, 22, 41, 5), 'elev': '20.8371'},
                               {'los': datetime(2018, 12, 4, 23, 43, 5), 'station': 'AS2',
                                'aos': datetime(2018, 12, 4, 23, 33, 8), 'elev': '23.546'}]

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

    @patch('trollsched.satpass.get_aqua_terra_dumpdata_from_ftp')
    def test_get_aqua_terra_dumps(self, dumps_from_ftp):
        dumps_from_ftp.return_value = self.dumpdata_terra

        # mymock:
        with patch('pyorbital.orbital.Orbital') as mymock:
            instance = mymock.return_value
            instance.get_next_passes = self.terra_orb.get_next_passes

            dumps = get_aqua_terra_dumps(datetime(2018, 12, 3, 0, 0),
                                         datetime(2018, 12, 10, 0, 0),
                                         self.terra_orb,
                                         self.terra)

            self.assertEqual(len(dumps), 4)
            self.assertEqual(dumps[0].station, 'SG2')
            self.assertEqual(dumps[0].max_elev, '9.2428')
            self.assertEqual(dumps[0].pass_direction(), 'ascending')
            self.assertEqual((dumps[0].risetime - datetime(2018, 12, 3, 1, 17, 53)).seconds, 0)
            self.assertEqual((dumps[0].falltime - datetime(2018, 12, 3, 1, 28, 14)).seconds, 0)

            self.assertEqual(dumps[3].station, 'AS2')
            self.assertEqual(dumps[3].max_elev, '23.546')
            self.assertEqual(dumps[3].pass_direction(), 'descending')
            self.assertEqual((dumps[3].risetime - datetime(2018, 12, 4, 23, 33, 8)).seconds, 0)
            self.assertEqual((dumps[3].falltime - datetime(2018, 12, 4, 23, 43, 5)).seconds, 0)

    @patch('os.path.exists')
    def test_get_metopa_passes(self, exists):

        exists.return_code = True

        # mymock:
        with patch('pyorbital.orbital.Orbital') as mymock:
            instance = mymock.return_value
            instance.get_next_passes = self.metopa_orb.get_next_passes

            metopa_passes = get_metopa_passes(self.metopa, self.metopa_passlist, self.metopa_orb)

            self.assertEqual(len(metopa_passes), 2)
            self.assertEqual(metopa_passes[0].pass_direction(), 'descending')
            self.assertEqual(metopa_passes[0].seconds(), 462.466119)
            self.assertEqual((metopa_passes[0].uptime - datetime(2018, 12, 4, 9, 17, 48, 530484)).seconds, 0)
            self.assertEqual((metopa_passes[0].risetime - datetime(2018, 12, 4, 9, 17, 46, 691075)).seconds, 0)

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
