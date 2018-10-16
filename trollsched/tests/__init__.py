#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 - 2018 PyTroll Community

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

"""Tests for scheduler.
"""

from trollsched.tests import (test_schedule, test_spherical, test_satpass)

import unittest


def suite():
    """The global test suite.
    """
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_schedule.suite())
    mysuite.addTests(test_spherical.suite())
    mysuite.addTests(test_satpass.suite())

    return mysuite
