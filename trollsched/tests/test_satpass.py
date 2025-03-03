#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2018 - 2024 Pytroll-schedule developers

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

"""Test the satellite pass and swath boundary classes."""

from datetime import datetime, timedelta

import numpy as np
import numpy.testing
import pytest
from pyorbital.orbital import Orbital
from pyresample.geometry import AreaDefinition, create_area_def

from trollsched.boundary import InstrumentNotSupported, SwathBoundary
from trollsched.satpass import Pass, create_pass

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
                  -52.140150361811244, -73.32968630186114], dtype="float64")

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
                  79.2863412851444, 83.25136141880888], dtype="float64")

LONS2 = np.array([-174.41109502, 167.84584132, 148.24213696, 130.10334782,
                  115.7074828, 105.07369809, 97.28481583, 91.4618503,
                  86.98024241, 83.4283141, 80.53652225, 78.1253594,
                  76.07228855, 74.29143113, 72.72103408, 71.31559576,
                  70.04080412, 68.87020177, 67.78293355, 66.76218577,
                  65.79407472, 64.86682945, 63.97016605, 63.09478077,
                  62.23190558, 61.37287373, 60.50863405, 59.62912286,
                  58.72232744, 57.77268809, 56.75796498, 55.6419694,
                  54.36007027, 41.41762911, 41.15660793, 40.9331126,
                  40.73252665, 40.54677784, 40.37092304, 40.20150965,
                  40.0358693, 39.87175642, 39.70713409, 39.54002703,
                  39.36840323, 39.1900621, 39.00251256, 38.80282499,
                  38.58743647, 38.35188019, 38.09039231, 37.79531831,
                  37.45618154, 37.05815986, 36.57947382, 35.98665163,
                  35.22533847, 34.20085643, 32.73220377, 30.42514135,
                  26.23397747, 16.29417395, -23.91719576, -102.71481425,
                  -122.5294795, -129.09284487], dtype="float64")

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
                  84.6670132], dtype="float64")

LONS3 = np.array([-8.66259458, -6.20984986, 15.99813586, 25.41134052, 33.80598414,
                  48.28641356, 49.55596283, 45.21769275, 43.95449327, 30.04053601,
                  22.33028017, 13.90584249, -5.59290326, -7.75625031], dtype="float64")

LATS3 = np.array([66.94713585, 67.07854554, 66.53108388, 65.27837805, 63.50223596,
                  58.33858588, 57.71210872, 55.14964148, 55.72506407, 60.40889798,
                  61.99561474, 63.11425455, 63.67173255, 63.56939058], dtype="float64")

AREA_DEF_EURON1 = AreaDefinition("euron1", "Northern Europe - 1km",
                                 "", {"proj": "stere", "ellps": "WGS84",
                                      "lat_0": 90.0, "lon_0": 0.0, "lat_ts": 60.0},
                                 3072, 3072, (-1000000.0, -4500000.0, 2072000.0, -1428000.0))


def get_n20_orbital():
    """Return the orbital instance for a given set of TLEs for NOAA-20.

    From 16 October 2018.
    """
    tle1 = "1 43013U 17073A   18288.00000000  .00000042  00000-0  20142-4 0  2763"
    tle2 = "2 43013 098.7338 224.5862 0000752 108.7915 035.0971 14.19549169046919"
    return Orbital("NOAA-20", line1=tle1, line2=tle2)


def get_n19_orbital():
    """Return the orbital instance for a given set of TLEs for NOAA-19.

    From 16 October 2018.
    """
    tle1 = "1 33591U 09005A   18288.64852564  .00000055  00000-0  55330-4 0  9992"
    tle2 = "2 33591  99.1559 269.1434 0013899 353.0306   7.0669 14.12312703499172"
    return Orbital("NOAA-19", line1=tle1, line2=tle2)


def get_mb_orbital():
    """Return orbital for a given set of TLEs for MetOp-B.

    From 2021-02-04
    """
    tle1 = "1 38771U 12049A   21034.58230818 -.00000012  00000-0  14602-4 0 9998"
    tle2 = "2 38771  98.6992  96.5537 0002329  71.3979  35.1836 14.21496632434867"
    return Orbital("Metop-B", line1=tle1, line2=tle2)


def get_s3a_orbital():
    """From 2022-06-06."""
    tle1 = "1 41335U 16011A   22157.82164820  .00000041  00000-0  34834-4 0  9994"
    tle2 = "2 41335  98.6228 225.2825 0001265  95.7364 264.3961 14.26738817328255"
    return Orbital("Sentinel-3A", line1=tle1, line2=tle2)


class TestPass:
    """Tests for the Pass object."""

    def setup_method(self):
        """Set up."""
        self.n20orb = get_n20_orbital()
        self.n19orb = get_n19_orbital()

    def test_pass_instrument_interface(self):
        """Test the intrument interface."""
        tstart = datetime(2018, 10, 16, 2, 48, 29)
        tend = datetime(2018, 10, 16, 3, 2, 38)

        instruments = set(("viirs", "avhrr", "modis", "mersi", "mersi-2"))
        for instrument in instruments:
            overp = Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument=instrument)
            assert overp.instrument == instrument

        instruments = set(("viirs", "avhrr", "modis"))
        overp = Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument=instruments)
        assert overp.instrument == "avhrr"

        instruments = set(("viirs", "modis"))
        overp = Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument=instruments)
        assert overp.instrument == "viirs"

        instruments = set(("amsu-a", "mhs"))
        with pytest.raises(TypeError):
            Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument=instruments)


class TestSwathBoundary:
    """Test the swath boundary object."""

    def setup_method(self):
        """Set up."""
        self.n20orb = get_n20_orbital()
        self.n19orb = get_n19_orbital()
        self.mborb = get_mb_orbital()
        self.s3aorb = get_s3a_orbital()
        self.euron1 = AREA_DEF_EURON1
        self.antarctica = create_area_def(
            "antarctic",
            {"ellps": "WGS84", "lat_0": "-90", "lat_ts": "-60",
             "lon_0": "0", "no_defs": "None", "proj": "stere",
             "type": "crs", "units": "m", "x_0": "0", "y_0": "0"},
            width=1000, height=1000,
            area_extent=(-4008875.4031, -4000855.294,
                         4000855.9937, 4008874.7048))
        self.arctica = create_area_def(
            "arctic",
            {"ellps": "WGS84", "lat_0": "90", "lat_ts": "60",
             "lon_0": "0", "no_defs": "None", "proj": "stere",
             "type": "crs", "units": "m", "x_0": "0", "y_0": "0"},
            width=1000, height=1000,
            area_extent=(-4008875.4031, -4000855.294,
                         4000855.9937, 4008874.7048))

    def test_swath_boundary(self):
        """Test generating a swath boundary."""
        tstart = datetime(2018, 10, 16, 2, 48, 29)
        tend = datetime(2018, 10, 16, 3, 2, 38)

        overp = Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument="viirs")
        overp_boundary = SwathBoundary(overp)

        cont = overp_boundary.contour()

        numpy.testing.assert_array_almost_equal(cont[0], LONS1)
        numpy.testing.assert_array_almost_equal(cont[1], LATS1)

        tstart = datetime(2018, 10, 16, 4, 29, 4)
        tend = datetime(2018, 10, 16, 4, 30, 29, 400000)

        overp = Pass("NOAA-20", tstart, tend, orb=self.n20orb, instrument="viirs")
        overp_boundary = SwathBoundary(overp, frequency=200)

        cont = overp_boundary.contour()

        numpy.testing.assert_array_almost_equal(cont[0], LONS2)
        numpy.testing.assert_array_almost_equal(cont[1], LATS2)

        # NOAA-19 AVHRR:
        tstart = datetime.strptime("20181016 04:00:00", "%Y%m%d %H:%M:%S")
        tend = datetime.strptime("20181016 04:01:00", "%Y%m%d %H:%M:%S")

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr")
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        numpy.testing.assert_array_almost_equal(cont[0], LONS3)
        numpy.testing.assert_array_almost_equal(cont[1], LATS3)

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr/3")
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        numpy.testing.assert_array_almost_equal(cont[0], LONS3)
        numpy.testing.assert_array_almost_equal(cont[1], LATS3)

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr-3")
        overp_boundary = SwathBoundary(overp, frequency=500)

        cont = overp_boundary.contour()

        numpy.testing.assert_array_almost_equal(cont[0], LONS3)
        numpy.testing.assert_array_almost_equal(cont[1], LATS3)

    def test_swath_coverage_does_not_cover_data_outside_area(self):
        """Test that swath covergate is 0 when the data is outside the area of interest."""
        # NOAA-19 AVHRR:
        tstart = datetime.strptime("20181016 03:54:13", "%Y%m%d %H:%M:%S")
        tend = datetime.strptime("20181016 03:55:13", "%Y%m%d %H:%M:%S")

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr")

        cov = overp.area_coverage(self.euron1)
        assert cov == 0

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr", frequency=80)

        cov = overp.area_coverage(self.euron1)
        assert cov == 0

    def test_swath_coverage_over_area(self):
        """Test that swath coverage matches when covering a part of the area of interest."""
        tstart = datetime.strptime("20181016 04:00:00", "%Y%m%d %H:%M:%S")
        tend = datetime.strptime("20181016 04:01:00", "%Y%m%d %H:%M:%S")

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr")

        cov = overp.area_coverage(self.euron1)
        assert cov == pytest.approx(0.103526, 1e-5)

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr", frequency=100)

        cov = overp.area_coverage(self.euron1)
        assert cov == pytest.approx(0.103526, 1e-5)

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr/3", frequency=133)

        cov = overp.area_coverage(self.euron1)
        assert cov == pytest.approx(0.103526, 1e-5)

        overp = Pass("NOAA-19", tstart, tend, orb=self.n19orb, instrument="avhrr", frequency=300)

        cov = overp.area_coverage(self.euron1)
        assert cov == pytest.approx(0.103526, 1e-5)

    def test_swath_coverage_metop(self):
        """Test ascat and avhrr coverages."""
        # ASCAT and AVHRR on Metop-B:
        tstart = datetime.strptime("2019-01-02T10:19:39", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=180)
        tle1 = "1 38771U 12049A   19002.35527803  .00000000  00000+0  21253-4 0 00017"
        tle2 = "2 38771  98.7284  63.8171 0002025  96.0390 346.4075 14.21477776326431"

        mypass = Pass("Metop-B", tstart, tend, instrument="ascat", tle1=tle1, tle2=tle2)
        cov = mypass.area_coverage(self.euron1)
        assert cov == pytest.approx(0.322812, 1e-5)

        mypass = Pass("Metop-B", tstart, tend, instrument="avhrr", tle1=tle1, tle2=tle2)
        cov = mypass.area_coverage(self.euron1)
        assert cov == pytest.approx(0.357324, 1e-5)

    def test_swath_coverage_slstr_not_supported(self):
        """Test Sentinel-3 SLSTR swath coverage - SLSTR is currently not supported!"""
        # Sentinel 3A slstr
        tstart = datetime(2022, 6, 6, 19, 58, 0)
        tend = tstart + timedelta(seconds=60)

        tle1 = "1 41335U 16011A   22156.83983125  .00000043  00000-0  35700-4 0  9996"
        tle2 = "2 41335  98.6228 224.3150 0001264  95.7697 264.3627 14.26738650328113"
        mypass = Pass("SENTINEL 3A", tstart, tend, instrument="slstr", tle1=tle1, tle2=tle2)

        with pytest.raises(InstrumentNotSupported) as exec_info:
            mypass.area_coverage(self.euron1)

        assert str(exec_info.value) == "SLSTR is a conical scanner, and currently not supported!"


    def test_swath_coverage_fy3(self):
        """Test FY3 coverages."""
        tstart = datetime.strptime("2019-01-05T01:01:45", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=60*15.5)

        tle1 = "1 43010U 17072A   18363.54078832 -.00000045  00000-0 -79715-6 0  9999"
        tle2 = "2 43010  98.6971 300.6571 0001567 143.5989 216.5282 14.19710974 58158"

        mypass = Pass("FENGYUN 3D", tstart, tend, instrument="mersi2", tle1=tle1, tle2=tle2, frequency=100)
        cov = mypass.area_coverage(self.euron1)
        assert cov == pytest.approx(0.786836, 1e-5)

        mypass = Pass("FENGYUN 3D", tstart, tend, instrument="mersi-2", tle1=tle1, tle2=tle2, frequency=100)
        cov = mypass.area_coverage(self.euron1)
        assert cov == pytest.approx(0.786836, 1e-5)

        tstart = datetime.strptime("2025-03-03T11:53:01", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=60*12.2)
        tle1 = "1 57490U 23111A   25061.77275788  .00000196  00000+0  11297-3 0  9996"
        tle2 = "2 57490  98.7372 134.1450 0001865  84.6719 275.4671 14.20001545 81964"
        mypass = Pass("FENGYUN 3F", tstart, tend, instrument="mersi-3", tle1=tle1, tle2=tle2, frequency=100)
        cov = mypass.area_coverage(self.euron1)
        assert cov == pytest.approx(0.70125, 1e-5)

    def test_arctic_is_not_antarctic(self):
        """Test that artic and antarctic are not mixed up."""
        tstart = datetime(2021, 2, 3, 16, 28, 3)
        tend = datetime(2021, 2, 3, 16, 31, 3)

        overp = Pass("Metop-B", tstart, tend, orb=self.mborb, instrument="avhrr")

        cov_south = overp.area_coverage(self.antarctica)
        cov_north = overp.area_coverage(self.arctica)

        assert cov_north == 0
        assert cov_south != 0


class TestPassList:
    """Tests for the pass list."""

    def test_meos_pass_list(self):
        """Test generating a meos pass list."""
        orig = ("  1 20190105 FENGYUN 3D  5907 52.943  01:01:45 n/a   01:17:15 15:30  18.6 107.4 -- "
                "Undefined(Scheduling not done 1546650105 ) a3d0df0cd289244e2f39f613f229a5cc D")

        tstart = datetime.strptime("2019-01-05T01:01:45", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=60 * 15.5)

        tle1 = "1 43010U 17072A   18363.54078832 -.00000045  00000-0 -79715-6 0  9999"
        tle2 = "2 43010  98.6971 300.6571 0001567 143.5989 216.5282 14.19710974 58158"

        mypass = Pass("FENGYUN 3D", tstart, tend, instrument="mersi2", tle1=tle1, tle2=tle2)
        coords = (10.72, 59.942, 0.1)
        meos_format_str = mypass.print_meos(coords, line_no=1)
        assert meos_format_str == orig

        mypass = Pass("FENGYUN 3D", tstart, tend, instrument="mersi-2", tle1=tle1, tle2=tle2)
        coords = (10.72, 59.942, 0.1)
        meos_format_str = mypass.print_meos(coords, line_no=1)
        assert meos_format_str == orig

    def test_generate_metno_xml(self):
        """Test generating a metno xml."""
        import xml.etree.ElementTree as ET  # noqa because defusedxml has no Element, see defusedxml#48
        root = ET.Element("acquisition-schedule")

        orig = ('<acquisition-schedule><pass satellite="FENGYUN 3D" aos="20190105010145" los="20190105011715" '
                'orbit="5907" max-elevation="52.943" asimuth-at-max-elevation="107.385" asimuth-at-aos="18.555" '
                'pass-direction="D" satellite-lon-at-aos="76.204" satellite-lat-at-aos="80.739" '
                'tle-epoch="20181229125844.110848" /></acquisition-schedule>')

        tstart = datetime.strptime("2019-01-05T01:01:45", "%Y-%m-%dT%H:%M:%S")
        tend = tstart + timedelta(seconds=60 * 15.5)

        tle1 = "1 43010U 17072A   18363.54078832 -.00000045  00000-0 -79715-6 0  9999"
        tle2 = "2 43010  98.6971 300.6571 0001567 143.5989 216.5282 14.19710974 58158"

        mypass = Pass("FENGYUN 3D", tstart, tend, instrument="mersi2", tle1=tle1, tle2=tle2)

        coords = (10.72, 59.942, 0.1)
        mypass.generate_metno_xml(coords, root)

        # Dictionaries don't have guaranteed ordering in Python 3.7, so convert the strings to sets and compare them
        res = set(ET.tostring(root).decode("utf-8").split())
        assert res == set(orig.split())

    def tearDown(self):
        """Clean up."""
        pass


@pytest.mark.usefixtures("fake_tle_file")
def test_create_pass(fake_tle_file):
    """Test creating a pass given a start and an end-time, platform, instrument and TLE-filepath."""
    starttime = datetime(2024, 9, 17, 1, 25, 52)
    endtime = starttime + timedelta(minutes=15)
    apass = create_pass("NOAA-20", "viirs", starttime, endtime, str(fake_tle_file))

    assert isinstance(apass, Pass)
    assert apass.risetime == datetime(2024, 9, 17, 1, 25, 52)
    assert apass.falltime == datetime(2024, 9, 17, 1, 40, 52)
    contours = apass.boundary.contour()

    np.testing.assert_array_almost_equal(contours[0], np.array([
        -70.36110203, -67.46084331, -37.31525344, 15.33742429,
        45.70673859, 58.20758754, 64.53232683, 68.33939709,
        70.91122142, 72.79385654, 74.25646605, 75.44689972,
        76.45347206, 77.33271412, 78.12312657, 78.85259555,
        79.5427448, 80.21178505, 80.87674043, 81.55576564,
        82.27163678, 83.05937724, 83.99111296, 84.1678927 ,
        71.14095118, 59.65602344, 50.0921629 , 42.33689333,
        36.07906912, 30.99289682, 26.80566985, 23.30761782,
        17.5195915, 16.76660935, 13.3562137 , 11.05669417,
         9.31517975, 7.90416741,  6.70609822, 5.65154688,
         4.69542071, 3.80603134,  2.95939777, 2.13592964,
         1.31825077, 0.48953714, -0.3680201, -1.27496102,
        -2.25693369, -3.34841058, -4.59918301, -6.08698296,
        -7.94529693, -10.43647307, -14.21005445, -15.07433768,
       -14.49280142, -14.53229624, -14.87821968, -15.67970647,
        -17.21558546, -20.06333116, -25.60407459, -37.83432182], dtype="float64"))

    np.testing.assert_array_almost_equal(contours[1], np.array([
        84.33463271, 84.996855, 87.51595932, 87.91311185, 87.0788011 ,
        86.08609579, 85.16108707, 84.31800847, 83.54109034, 82.81224557,
        82.11521918, 81.4355823 , 80.75994129, 80.07498982, 79.36644307,
        78.61771056, 77.80801807, 76.90942731, 75.88160126, 74.66160838,
        73.14141261, 71.10852475, 68.03499509, 67.34603232, 66.29582898,
        64.22360654, 61.34855651, 57.88332007, 53.99536906, 49.80456636,
        45.39366513, 40.81954787, 30.76199763, 30.98222964, 31.91529108,
        32.48607791, 32.88793494, 33.1947272 , 33.44223661, 33.65034541,
        33.8311831 , 33.99269414, 34.14039982, 34.27835064, 34.40968082,
        34.53695383, 34.66239272, 34.78804435, 34.91590262, 35.047996  ,
        35.18641366, 35.33315856, 35.48940778, 35.65238903, 35.79943487,
        35.81635723, 46.51491314, 51.61438345, 56.69722991, 61.75596299,
        66.7771717 , 71.73310959, 76.55555826, 81.03181048], dtype="float64"))
