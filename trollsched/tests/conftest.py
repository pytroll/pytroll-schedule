#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c22526.ad.smhi.se>

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

"""Fixtures for unit testing."""


from datetime import datetime, timedelta

import pytest

from trollsched.satpass import create_pass

# Create a fake TLE file for testing.
TEST_TLEFILE1 = """NOAA-20
1 43013U 17073A   24260.23673431  .00000000  00000+0  20543-3 0 00012
2 43013  98.7102 196.6451 0000981  27.6583 181.3633 14.19588904353826
"""

FAKE_XML_SCHEDULE = """<acquisition-schedule><properties><project>Pytroll</project><type>request</type><station>nrk</station><file-start>2024-12-02-21:54:34</file-start><file-end>2024-12-02-23:54:34</file-end><requested-by>SMHI</requested-by><requested-on>2024-12-02-19:55:05</requested-on></properties><pass satellite="noaa18" start-time="2024-12-02-22:01:31" end-time="2024-12-02-22:16:17" /><pass satellite="npp" start-time="2024-12-02-22:38:18" end-time="2024-12-02-22:49:39" /><pass satellite="noaa20" start-time="2024-12-02-23:01:26" end-time="2024-12-02-23:14:17" /></acquisition-schedule>"""  # noqa

TEST_TLEFILE2 = """NOAA-20
1 43013U 17073A   24336.73680280  .00000000  00000+0  21070-3 0 00017
2 43013  98.7240 271.9468 0002361  91.7021 243.0836 14.19579317364671
NOAA 18
1 28654U 05018A   24336.77887833  .00000603  00000+0  34401-3 0  9991
2 28654  98.8621  52.2610 0013351 310.3593  49.6413 14.13439071  6932
SUOMI NPP
1 37849U 11061A   24336.86623462  .00000353  00000+0  18825-3 0  9992
2 37849  98.7391 272.4265 0002092 106.9131 253.2275 14.19509233678657
"""

@pytest.fixture
def fake_tle_file(tmp_path):
    """Write fake TLE file."""
    file_path = tmp_path / "sometles.txt"
    with open(file_path, "w") as fpt:
        fpt.write(TEST_TLEFILE1)

    return file_path

@pytest.fixture
def fake_long_tle_file(tmp_path):
    """Write fake TLE file."""
    file_path = tmp_path / "some_more_tles.txt"
    with open(file_path, "w") as fpt:
        fpt.write(TEST_TLEFILE2)

    return file_path

@pytest.fixture
def fake_noaa20_viirs_pass_instance(fake_tle_file):
    """Create a fake trollsched.,satpass instance for NOAA-20 VIIRS."""
    starttime = datetime(2024, 9, 17, 1, 25, 52)
    endtime = starttime + timedelta(minutes=15)
    return create_pass("NOAA-20", "viirs", starttime, endtime, str(fake_tle_file))

@pytest.fixture
def fake_xml_schedule_file(tmp_path):
    """Create a fake XML schedule file."""
    file_path = tmp_path / "20241202-195434-aquisition-schedule-request-nrk.xml"
    with open(file_path, "w") as fpt:
        fpt.write(FAKE_XML_SCHEDULE)

    return file_path
