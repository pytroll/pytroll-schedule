#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll Community


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

"""Test the creation and storing of satellite pass instrument swath outlines."""


import numpy as np
import pytest
import shapely

from trollsched.shapefiles_from_schedule import (
    create_shapefile_filename,
    create_shapefile_from_pass,
    get_shapely_polygon_from_lonlat,
    run,
)

CONT_POLY_LON = np.array([-2.13452262, -2.29587336, -2.71901413, 2.49884179, 1.84466852,
                          1.62370346, 1.52297613, 1.46583198, 1.42884581, 1.40275172,
                          1.38319053, 1.36784925, 1.35538792, 1.34497645, 1.33607268,
                          1.32830685, 1.32141745, 1.31521366, 1.30955222, 1.30432275,
                          1.29943799, 1.29482714, 1.29043107, 1.28619883, 1.2820848 ,
                          1.27804628, 1.27404108, 1.27002472, 1.26594662, 1.26174386,
                          1.25732922, 1.25256377, 1.24720486, 1.01188742, 0.79409035,
                          0.6071801, 0.4536726, 0.32966584, 0.22939162, 0.14748665,
                          -0.00471426, -0.05349142, -0.08929434, -0.11753567, -0.14089485,
                          -0.1608663 , -0.17837565, -0.1940391, -0.20828969, -0.2214453,
                          -0.23374812, -0.24538904, -0.25652346, -0.26728206, -0.27777851,
                          -0.28811526, -0.29838809, -0.30869016, -0.31911562, -0.32976349,
                          -0.34074184, -0.35217288, -0.36419961, -0.37699468, -0.39077317,
                          -0.40581154, -0.42247709, -0.44127621, -0.46293831, -0.48857303,
                          -0.5199966 , -0.56051525, -0.61729867, -0.61753171, -0.62771892,
                          -0.64829722, -0.68661816, -0.75917567, -0.9100173, -1.27984447],
                         "float64")
CONT_POLY_LAT = np.array([1.47665956, 1.5181844 , 1.54460102, 1.55163844, 1.54000093,
                          1.52570241, 1.51231719, 1.50010348, 1.48891231, 1.47855701,
                          1.46887101, 1.45971388, 1.45096745, 1.44253069, 1.43431507,
                          1.42624067, 1.41823286, 1.41021941, 1.40212764, 1.39388163,
                          1.38539899, 1.37658705, 1.36733796, 1.35752199, 1.34697812,
                          1.33549993, 1.32281374, 1.30854277, 1.29214467, 1.27279422,
                          1.24913883, 1.21870824, 1.17637107, 1.16993658, 1.14394582,
                          1.10161427, 1.04685251, 0.98308241, 0.91287235, 0.83802526,
                          0.63421758, 0.64926881, 0.65941499, 0.66689569, 0.67274331,
                          0.67750383, 0.68149934, 0.68493465, 0.68794744, 0.69063434,
                          0.69306561, 0.6952939 , 0.6973597 , 0.69929488, 0.7011251 ,
                          0.7028715 , 0.7045518 , 0.70618126, 0.70777322, 0.70933965,
                          0.7108915 , 0.71243895, 0.71399159, 0.71555848, 0.71714804,
                          0.71876762, 0.72042255, 0.72211402, 0.72383449, 0.72555719,
                          0.72721033, 0.72860415, 0.72917251, 0.95338892, 1.04189611,
                          1.12986537, 1.21694627, 1.30234019, 1.38380771, 1.45301036],
                         "float64")


def test_get_shapely_polygon_from_lonlat():
    """Test getting a polygon from a swath contour polygon."""
    sat_poly = get_shapely_polygon_from_lonlat(CONT_POLY_LON, CONT_POLY_LAT)

    assert isinstance(sat_poly, shapely.geometry.polygon.Polygon)
    assert sat_poly.is_valid
    assert sat_poly.is_simple
    assert sat_poly.is_closed is False
    np.testing.assert_almost_equal(sat_poly.length, 649.84664, decimal=5)


@pytest.mark.usefixtures("fake_noaa20_viirs_pass_instance")
def test_create_shapefile_filename(fake_noaa20_viirs_pass_instance):
    """Test creating the shapefile filename from a NOAA-20/VIIRS trollsched.satpass instance."""
    result = create_shapefile_filename(fake_noaa20_viirs_pass_instance)
    assert result == "viirs_NOAA-20_202409170125_202409170140_outline.shp"


@pytest.mark.usefixtures("fake_noaa20_viirs_pass_instance")
def test_create_shapefile_from_pass(fake_noaa20_viirs_pass_instance, tmp_path):
    """Test create a shapefile from a fake NOAA-20/VIIRS pass."""
    tmp_filename = tmp_path / "fake_noaa20_viirs_shapefile.shp"
    create_shapefile_from_pass(fake_noaa20_viirs_pass_instance, tmp_filename)
    assert tmp_filename.exists()
    nfiles = len([f for f in tmp_filename.parent.glob("fake_noaa20_viirs_shapefile.*")])
    assert nfiles == 5


@pytest.mark.usefixtures("fake_long_tle_file")
@pytest.mark.usefixtures("fake_xml_schedule_file")
def test_create_shapefiles_from_schedule(fake_xml_schedule_file, fake_long_tle_file, tmp_path):
    """Test running the script to create shapefiles from an xml schedule request file."""
    output_dir = tmp_path / "results"
    output_dir.mkdir()

    args = ["-s", "Suomi NPP", "NOAA 18", "NOAA-20", "-t", str(fake_long_tle_file),
            "-x", str(fake_xml_schedule_file),
            "-o", str(output_dir)]
    run(args)

    n18files = [str(name) for name in output_dir.glob("avhrr_NOAA-18_202412022201_202412022216_outline.*")]
    assert len(n18files) == 5
    nppfiles = [str(name) for name in output_dir.glob("viirs_Suomi-NPP_202412022238_202412022249_outline.*")]
    assert len(nppfiles) == 5
    n20files = [str(name) for name in output_dir.glob("viirs_NOAA-20_202412022301_202412022314_outline.*")]
    assert len(n20files) == 5
