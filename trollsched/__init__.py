#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 - 2019 PyTroll Community

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

"""Package file."""

from . import version

__version__ = version.get_versions()['version']


# shortest allowed pass in minutes
MIN_PASS = 4

# DRL still use the name JPSS-1 in the TLEs:
NOAA20_NAME = {'NOAA-20': 'JPSS-1'}

NUMBER_OF_FOVS = {
    'avhrr': 2048,
    'mhs': 90,
    'amsua': 30,
    'mwhs2': 98,
    'atms': 96,
    'ascat': 42,
    'viirs': 6400,
    'atms': 96,
    'mwhs-2': 98
}

SATELLITE_NAMES = {'npp': 'Suomi NPP',
                   'noaa19': 'NOAA 19',
                   'noaa18': 'NOAA 18',
                   'noaa15': 'NOAA 15',
                   'aqua': 'Aqua',
                   'terra': 'Terra',
                   'metopc': 'Metop-C',
                   'metopb': 'Metop-B',
                   'metopa': 'Metop-A',
                   'noaa20': 'NOAA-20',
                   'fengyun3d': 'FY-3D',
                   'fengyun3c': 'FY-3C'
                   }

INSTRUMENT = {'Suomi NPP': 'viirs',
              'NOAA-20': 'viirs',
              'Aqua': 'modis',
              'Terra': 'modis',
              'NOAA 19': 'avhrr',
              'NOAA 18': 'avhrr',
              'NOAA 15': 'avhrr',
              'Metop-A': 'avhrr',
              'Metop-B': 'avhrr',
              'Metop-C': 'avhrr',
              'FY-3D': 'avhrr',
              'FY-3C': 'avhrr'}

from . import version
__version__ = version.get_versions()['version']
