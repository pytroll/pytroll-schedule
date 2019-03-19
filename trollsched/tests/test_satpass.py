#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2018 - 2019 PyTroll

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>

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

"""Test the satellite pass and swath boundary classes
"""

import unittest
import numpy as np
from datetime import datetime, timedelta
from trollsched.satpass import Pass
from trollsched.boundary import SwathBoundary
from pyorbital.orbital import Orbital


LONS1 = np.array([-122.29913729160562, -131.54385362589042, -155.788034272281,
                  143.1730880418349, 105.69172088208997, 93.03135571771092,
                  87.26010432019743, 83.98598584966442, 81.86683434871546,
                  80.37175346216411, 79.2509798668123, 78.37198926578984,
                  77.65800714027662, 77.06147400915819, 76.55132566889495,
                  76.10637628220547, 75.71164306799828, 75.35619180525052,
                  75.03181505238287, 74.73218847041143, 74.45231256197947,
                  74.18813012461848, 73.9362540393912, 73.69376447765231,
                  73.45804804675883, 73.22665809422263, 72.99717692793544,
                  72.76705638792168, 72.53339841609603, 72.2925978414254,
                  72.03965795937306, 71.76661774368146, 71.45957469190316,
                  57.97687872167697, 45.49802548616658, 34.788857347919546,
                  25.993525469714424, 18.88846123000295, 13.14317179269443,
                  8.450362684274728, -0.27010733525252295, -3.0648326302431794,
                  -5.116189000358824, -6.73429807721795, -8.072680053386163,
                  -9.21696007364773, -10.220171884036919, -11.11762132513045,
                  -11.934120125548072, -12.687881125682765, -13.392781001351315,
                  -14.059756511026736, -14.69771138916782, -15.314133712696703,
                  -15.915536409615735, -16.507788289068856, -17.09637839792269,
                  -17.686643087306685, -18.283978247944123, -18.894056410060063,
                  -19.523069195727878, -20.17801994245519, -20.867100607022966,
                  -21.600204055760642, -22.389653641849733, -23.251288693929943,
                  -24.206153922914886, -25.283264445138713, -26.524411381004743,
                  -27.993172418988525, -29.79361072725673, -32.11515837055801,
                  -35.36860848223405, -35.38196057933595, -35.96564490844792,
                  -37.14469461070555, -39.34032289002443, -43.49756191648018,
                  -52.140150361811244, -73.32968630186114], dtype='float64')

LATS1 = np.array([84.60636067724808, 86.98555849233523, 88.49911967556697,
                  88.90233393880413, 88.23555365613707, 87.41630911481282,
                  86.64939216187459, 85.94959841469182, 85.30839167814023,
                  84.71507625588431, 84.16010931725756, 83.63544438659248,
                  83.13431099825148, 82.65092034888734, 82.18020003036649,
                  81.71757084925224, 81.25875743723827, 80.79962022032255,
                  80.33599602524967, 79.86353436733512, 79.37751495806062,
                  78.87262831355378, 78.3426942980262, 77.78028071690198,
                  77.17616119674511, 76.51850934329316, 75.79164459606967,
                  74.97397797613992, 74.03443588562436, 72.92573674313518,
                  71.57038280824118, 69.82683886377178, 67.40109717220513,
                  67.03242839212335, 65.54326755696877, 63.11784822611803,
                  59.98023069591168, 56.32647323215378, 52.30373268534935,
                  48.01531077177335, 36.33799056582854, 37.200362356448125,
                  37.78169598891329, 38.210308430109684, 38.54535234179983,
                  38.8181101172057, 39.0470359762339, 39.24386487280032,
                  39.41648482997921, 39.57043267820405, 39.70973443234515,
                  39.83740623634436, 39.955767569171485, 40.06664498984812,
                  40.17150923539549, 40.271570238680745, 40.36784473887322,
                  40.46120553672548, 40.55241811035527, 40.64216822927882,
                  40.7310828091462, 40.819745180454284, 40.90870492549053,
                  40.99848114410508, 41.08955592221846, 41.18235086149538,
                  41.27717142920562, 41.37408580927609, 41.472661254399455,
                  41.57136466452366, 41.66608254408796, 41.745942562974314,
                  41.77850750277849, 54.62516158367828, 59.69624962433962,
                  64.7365168572082, 69.72588498397877, 74.61859631181376,
                  79.2863412851444, 83.25136141880888], dtype='float64')

LONS2 = np.array([-174.41109502,  167.84584132,  148.24213696,  130.10334782,
                  115.7074828,  105.07369809,   97.28481583,   91.4618503,
                  86.98024241,   83.4283141,   80.53652225,   78.1253594,
                  76.07228855,   74.29143113,   72.72103408,   71.31559576,
                  70.04080412,   68.87020177,   67.78293355,   66.76218577,
                  65.79407472,   64.86682945,   63.97016605,   63.09478077,
                  62.23190558,   61.37287373,   60.50863405,   59.62912286,
                  58.72232744,   57.77268809,   56.75796498,   55.6419694,
                  54.36007027,   41.41762911,   41.15660793,   40.9331126,
                  40.73252665,   40.54677784,   40.37092304,   40.20150965,
                  40.0358693,   39.87175642,   39.70713409,   39.54002703,
                  39.36840323,   39.1900621,   39.00251256,   38.80282499,
                  38.58743647,   38.35188019,   38.09039231,   37.79531831,
                  37.45618154,   37.05815986,   36.57947382,   35.98665163,
                  35.22533847,   34.20085643,   32.73220377,   30.42514135,
                  26.23397747,   16.29417395,  -23.91719576, -102.71481425,
                  -122.5294795, -129.09284487], dtype='float64')

LATS2 = np.array([83.23214786, 84.90973645, 85.62529048, 85.74243351, 85.52147568,
                  85.13874302, 84.69067959, 84.22338069, 83.75720094, 83.30023412,
                  82.85480916, 82.42053485, 81.9957309, 81.57810129, 81.16504231,
                  80.75376801, 80.34133891, 79.92463458, 79.50028749, 79.0645828,
                  78.61332046, 78.14162813, 77.64370408, 77.11245516, 76.5389713,
                  75.91173559, 75.21538754, 74.42869094, 73.52099029, 72.44554294,
                  71.12561977, 69.42093758, 67.03973793, 67.40770791, 69.8341456,
                  71.57844446, 72.93459921, 74.04414258, 74.98457279, 75.80317362,
                  76.53102217, 77.1897121, 77.79492994, 78.3585095, 78.88968633,
                  79.39590402, 79.88335693, 80.35737249, 80.8226939, 81.28370137,
                  81.74459732, 82.20957417, 82.68298027, 83.16949849, 83.67435372,
                  84.20356848, 84.76429067, 85.36521771, 86.01711637, 86.73327122,
                  87.5286869, 88.40887156, 89.21959299, 88.71884272, 87.09172665,
                  84.6670132], dtype='float64')

LONS3 = np.array([-8.66259458, -6.20984986, 15.99813586, 25.41134052, 33.80598414,
                  48.28641356, 49.55596283, 45.21769275, 43.95449327, 30.04053601,
                  22.33028017, 13.90584249, -5.59290326, -7.75625031], dtype='float64')

LATS3 = np.array([66.94713585, 67.07854554, 66.53108388, 65.27837805, 63.50223596,
                  58.33858588, 57.71210872, 55.14964148, 55.72506407, 60.40889798,
                  61.99561474, 63.11425455, 63.67173255, 63.56939058], dtype='float64')


def assertNumpyArraysEqual(self, other):
    if self.shape != other.shape:
        raise AssertionError("Shapes don't match")
    if not np.allclose(self, other):
        raise AssertionError("Elements don't match!")


def get_n20_orbital():
    """Return the orbital instance for a given set of TLEs for NOAA-20.
    From 16 October 2018.
    """
    tle1 = "1 43013U 17073A   18288.00000000  .00000042  00000-0  20142-4 0  2763"
    tle2 = "2 43013 098.7338 224.5862 0000752 108.7915 035.0971 14.19549169046919"
    return Orbital('NOAA-20', line1=tle1, line2=tle2)


def get_n19_orbital():
    """Return the orbital instance for a given set of TLEs for NOAA-19.
    From 16 October 2018.
    """
    tle1 = "1 33591U 09005A   18288.64852564  .00000055  00000-0  55330-4 0  9992"
    tle2 = "2 33591  99.1559 269.1434 0013899 353.0306   7.0669 14.12312703499172"
    return Orbital('NOAA-19', line1=tle1, line2=tle2)


def get_region(areaid):
    try:
        from satpy.resample import get_area_def
    except ImportError:
        from mpop.projector import get_area_def

    return get_area_def(areaid)


class TestPass(unittest.TestCase):

    def setUp(self):
        """Set up"""
        self.n20orb = get_n20_orbital()
        self.n19orb = get_n19_orbital()

    def test_pass_instrument_interface(self):

        tstart = datetime(2018, 10, 16, 2, 48, 29)
        tend = datetime(2018, 10, 16, 3, 2, 38)

        instruments = set(('viirs', 'avhrr', 'modis'))
        overp = Pass('NOAA-20', tstart, tend, orb=self.n20orb, instrument=instruments)
        self.assertEqual(overp.instrument, 'avhrr')

        instruments = set(('viirs', 'modis'))
        overp = Pass('NOAA-20', tstart, tend, orb=self.n20orb, instrument=instruments)
        self.assertEqual(overp.instrument, 'viirs')

        instruments = set(('amsu-a', 'mhs'))
        self.assertRaises(TypeError, Pass, self,
                          'NOAA-20', tstart, tend, orb=self.n20orb, instrument=instruments)

    def tearDown(self):
        """Clean up"""
        pass


class TestSwathBoundary(unittest.TestCase):

    def setUp(self):
        """Set up"""
        self.n20orb = get_n20_orbital()
        self.n19orb = get_n19_orbital()
        self.euron1 = get_region('euron1')

    def test_swath_boundary(self):

        tstart = datetime(2018, 10, 16, 2, 48, 29)
        tend = datetime(2018, 10, 16, 3, 2, 38)

        overp = Pass('NOAA-20', tstart, tend, orb=self.n20orb, instrument='viirs')
        overp_boundary = SwathBoundary(overp)

        cont = overp_boundary.contour()

        assertNumpyArraysEqual(cont[0], LONS1)
        assertNumpyArraysEqual(cont[1], LATS1)

        tstart = datetime(2018, 10, 16, 4, 29, 4)
        tend = datetime(2018, 10, 16, 4, 30, 29, 400000)

        overp = Pass('NOAA-20', tstart, tend, orb=self.n20orb, instrument='viirs')
        overp_boundary = SwathBoundary(overp, frequency=200)

        cont = overp_boundary.contour()

        assertNumpyArraysEqual(cont[0], LONS2)
        assertNumpyArraysEqual(cont[1], LATS2)

        # NOAA-19 AVHRR:
        tstart = datetime.strptime('20181016 04:00:00', '%Y%m%d %H:%M:%S')
        tend = datetime.strptime('20181016 04:01:00', '%Y%m%d %H:%M:%S')

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr')
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        assertNumpyArraysEqual(cont[0], LONS3)
        assertNumpyArraysEqual(cont[1], LATS3)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr/3')
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        assertNumpyArraysEqual(cont[0], LONS3)
        assertNumpyArraysEqual(cont[1], LATS3)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr-3')
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        assertNumpyArraysEqual(cont[0], LONS3)
        assertNumpyArraysEqual(cont[1], LATS3)

    def test_swath_coverage(self):

        # NOAA-19 AVHRR:
        tstart = datetime.strptime('20181016 03:54:13', '%Y%m%d %H:%M:%S')
        tend = datetime.strptime('20181016 03:55:13', '%Y%m%d %H:%M:%S')

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr')

        cov = overp.area_coverage(self.euron1)
        self.assertEqual(cov, 0)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr', frequency=80)

        cov = overp.area_coverage(self.euron1)
        self.assertEqual(cov, 0)

        tstart = datetime.strptime('20181016 04:00:00', '%Y%m%d %H:%M:%S')
        tend = datetime.strptime('20181016 04:01:00', '%Y%m%d %H:%M:%S')

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr')

        cov = overp.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.103526, 5)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr', frequency=100)

        cov = overp.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.103526, 5)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr/3', frequency=133)

        cov = overp.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.103526, 5)

        overp = Pass('NOAA-19', tstart, tend, orb=self.n19orb, instrument='avhrr', frequency=300)

        cov = overp.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.103526, 5)

        # ASCAT and AVHRR on Metop-B:
        tstart = datetime.strptime("2019-01-02T10:19:39", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=180)
        tle1 = '1 38771U 12049A   19002.35527803  .00000000  00000+0  21253-4 0 00017'
        tle2 = '2 38771  98.7284  63.8171 0002025  96.0390 346.4075 14.21477776326431'

        mypass = Pass('Metop-B', tstart, tend, instrument='ascat', tle1=tle1, tle2=tle2)
        cov = mypass.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.322812, 5)

        mypass = Pass('Metop-B', tstart, tend, instrument='avhrr', tle1=tle1, tle2=tle2)
        cov = mypass.area_coverage(self.euron1)
        self.assertAlmostEqual(cov, 0.357324, 5)

        tstart = datetime.strptime("2019-01-05T01:01:45", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=60*15.5)

        tle1 = '1 43010U 17072A   18363.54078832 -.00000045  00000-0 -79715-6 0  9999'
        tle2 = '2 43010  98.6971 300.6571 0001567 143.5989 216.5282 14.19710974 58158'

        mypass = Pass('FENGYUN 3D', tstart, tend, instrument='mersi2', tle1=tle1, tle2=tle2)
        cov = mypass.area_coverage(self.euron1)

        self.assertAlmostEqual(cov, 0.786836, 5)

    def tearDown(self):
        """Clean up"""
        pass


def suite():
    """The suite for test_satpass
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestSwathBoundary))
    mysuite.addTest(loader.loadTestsFromTestCase(TestPass))

    return mysuite
