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

@pytest.fixture
def fake_tle_file(tmp_path):
    """Write fake TLE file."""
    file_path = tmp_path / "sometles.txt"
    with open(file_path, "w") as fpt:
        fpt.write(TEST_TLEFILE1)

    return file_path

@pytest.fixture
def fake_noaa20_viirs_pass_instance(fake_tle_file):
    """Create a fake trollsched.,satpass instance for NOAA-20 VIIRS."""
    starttime = datetime(2024, 9, 17, 1, 25, 52)
    endtime = starttime + timedelta(minutes=15)
    return create_pass("NOAA-20", "viirs", starttime, endtime, str(fake_tle_file))
