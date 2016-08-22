#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2016 Martin Raspaud

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

"""
"""
# workaround python bug: http://bugs.python.org/issue15881#msg170215
import multiprocessing

from setuptools import setup
import imp
import sys

version = imp.load_source('trollsched.version', 'trollsched/version.py')

requires = ['numpy', 'pyresample', 'pyorbital']

if sys.version_info < (2, 7):
    # multiprocessing is not in the standard library
    requires.append('argparse')

setup(name='pytroll-schedule',
      version=version.__version__,
      description='Scheduling satellite passes in Python',
      author='Martin Raspaud',
      author_email='martin.raspaud@smhi.se',
      classifiers=["Development Status :: 4 - Beta",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: Astronomy"],
      test_suite='trollsched.tests.suite',
      entry_points={
          'console_scripts': ['schedule = trollsched.schedule:run',
                              'compare_scheds = trollsched.compare:run']},
      scripts=['generate_schedule_xmlpage.py'],
      packages=['trollsched'],
      install_requires=requires,
      zip_safe=False,
      )
