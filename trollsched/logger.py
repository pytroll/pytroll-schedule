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

"""The log handling."""

import logging
import logging.config
import logging.handlers

import yaml

LOG_FORMAT = "[%(asctime)s %(levelname)-8s] %(message)s"

log_levels = {
    0: logging.WARN,
    1: logging.INFO,
    2: logging.DEBUG,
}


def setup_logging(cmd_args):
    """Set up logging."""
    if cmd_args.log_config is not None:
        with open(cmd_args.log_config) as fd:
            log_dict = yaml.safe_load(fd.read())
            logging.config.dictConfig(log_dict)
            return

    root = logging.getLogger("")
    root.setLevel(log_levels[cmd_args.verbosity])

    fh_ = logging.StreamHandler()

    formatter = logging.Formatter(LOG_FORMAT)
    fh_.setFormatter(formatter)

    root.addHandler(fh_)
