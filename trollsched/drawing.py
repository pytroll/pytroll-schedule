#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2018 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c20671.ad.smhi.se>

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

"""Drawing satellite overpass outlines on maps
"""

import os
import logging
import logging.handlers
import numpy as np
import matplotlib as mpl
MPL_BACKEND = mpl.get_backend()

logger = logging.getLogger(__name__)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    BASEMAP_NOT_CARTOPY = False
except ImportError:
    logger.warning("Failed loading Cartopy, will try Basemap instead")
    BASEMAP_NOT_CARTOPY = True


class MapperBasemap(object):
    """A class to generate nice plots with basemap.
    """

    def __init__(self, **proj_info):

        from mpl_toolkits.basemap import Basemap

        if not proj_info:
            proj_info = {
                'projection': 'nsper',
                'lat_0': 58,
                'lon_0': 16,
                'resolution': 'l',
                'area_thresh': 1000.
            }

        self.map = Basemap(**proj_info)

        self.map.drawcoastlines()
        self.map.drawcountries()

        self.map.drawmapboundary(fill_color='white')

        self.map.drawmeridians(np.arange(0, 360, 30))
        self.map.drawparallels(np.arange(-90, 90, 30))

    def __enter__(self):
        return self.map

    def __exit__(self, etype, value, tb):
        pass


class MapperCartopy(object):
    """A class to generate nice plots with Cartopy.
    """

    def __init__(self, **proj_info):

        mpl.use(MPL_BACKEND)
        import matplotlib.pyplot as plt

        if not proj_info:
            proj_info = {
                'central_latitude': 58,
                'central_longitude': 16,
                'satellite_height': 35785831,
                'false_easting': 0,
                'false_northing': 0,
                'globe': None
            }

        fig = plt.figure(figsize=(8, 6))

        self._ax = fig.add_subplot(
            1, 1, 1, projection=ccrs.NearsidePerspective(**proj_info))

        self._ax.add_feature(cfeature.OCEAN, zorder=0)
        self._ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
        self._ax.add_feature(cfeature.BORDERS, zorder=0)

        self._ax.set_global()
        self._ax.gridlines()

    def plot(self, *args, **kwargs):
        mpl.use(MPL_BACKEND)
        import matplotlib.pyplot as plt

        kwargs['transform'] = ccrs.Geodetic()
        return plt.plot(*args, **kwargs)

    def nightshade(self, utctime, **kwargs):

        from trollsched.helper_functions import fill_dark_side

        color = kwargs.get('color', 'black')
        alpha = kwargs.get('alpha', 0.4)
        fill_dark_side(self._ax, time=utctime, color=color, alpha=alpha)

    def __call__(self, *args):
        return args

    def __enter__(self):
        return self

    def __exit__(self, etype, value, tb):
        pass


if BASEMAP_NOT_CARTOPY:
    Mapper = MapperBasemap
else:
    Mapper = MapperCartopy


def save_fig(pass_obj,
             poly=None,
             directory="/tmp/plots",
             overwrite=False,
             labels=None,
             extension=".png",
             outline='-r'):
    """Save the pass as a figure. Filename is automatically generated.
    """
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.clf()

    logger.debug("Save fig " + str(pass_obj))
    rise = pass_obj.risetime.strftime("%Y%m%d%H%M%S")
    fall = pass_obj.falltime.strftime("%Y%m%d%H%M%S")
    if not os.path.exists(directory):
        logger.debug("Create plot dir " + directory)
        os.makedirs(directory)
    filename = os.path.join(
        directory,
        (rise + pass_obj.satellite.name.replace(" ", "_") + fall + extension))

    pass_obj.fig = filename
    if not overwrite and os.path.exists(filename):
        return filename

    logger.debug("Filename = <%s>", filename)
    with Mapper() as mapper:
        mapper.nightshade(pass_obj.uptime, alpha=0.2)
        logger.debug("Draw: outline = <%s>", outline)
        draw(pass_obj.boundary.contour_poly, mapper, outline)
        if poly is not None:
            draw(poly, mapper, "-b")

    logger.debug("Title = %s", str(pass_obj))
    plt.title(str(pass_obj))
    for label in labels or []:
        plt.figtext(*label[0], **label[1])
    logger.debug("Save plot...")
    plt.savefig(filename)
    logger.debug("Return...")
    return filename


def show(pass_obj,
         poly=None,
         labels=None,
         other_poly=None,
         proj=None,
         outline='-r'):
    """Show the current pass on screen (matplotlib, basemap).
    """
    mpl.use(MPL_BACKEND)
    import matplotlib.pyplot as plt

    proj = proj or {}
    with Mapper(**proj) as mapper:
        mapper.nightshade(pass_obj.uptime, alpha=0.2)
        draw(pass_obj.boundary.contour_poly, mapper, outline)
        if poly is not None:
            draw(poly, mapper, "-b")
        if other_poly is not None:
            draw(other_poly, mapper, "-g")
    plt.title(str(pass_obj))
    for label in (labels or []):
        plt.figtext(*label[0], **label[1])
    plt.show()


def draw(poly, mapper, options, **more_options):
    lons = np.rad2deg(poly.lon.take(np.arange(len(poly.lon) + 1), mode="wrap"))
    lats = np.rad2deg(poly.lat.take(np.arange(len(poly.lat) + 1), mode="wrap"))
    rx, ry = mapper(lons, lats)
    mapper.plot(rx, ry, options, **more_options)


def main():
    from trollsched.satpass import get_next_passes
    from datetime import datetime

    passes = get_next_passes(["noaa 19", "suomi npp"], datetime.now(), 24, (16, 58, 0))
    for p in passes:
        save_fig(p, directory="/tmp/plots/")


if __name__ == '__main__':
    main()
