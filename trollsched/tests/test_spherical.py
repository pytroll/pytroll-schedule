#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014, 2015, 2018 Martin Raspaud

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

"""Test cases for spherical geometry.
"""

from trollsched.spherical import SphPolygon, Arc, SCoordinate, CCoordinate
import unittest
import numpy as np


class TestSCoordinate(unittest.TestCase):

    """Test SCoordinates.
    """

    def test_distance(self):
        """Test Vincenty formula
        """
        d = SCoordinate(0, 0).distance(SCoordinate(1, 1))
        self.assertEquals(d, 1.2745557823062943)

    def test_hdistance(self):
        """Test Haversine formula
        """
        d = SCoordinate(0, 0).hdistance(SCoordinate(1, 1))
        self.assertTrue(np.allclose(d, 1.2745557823062943))

    def test_str(self):
        """Check the string representation
        """
        d = SCoordinate(0, 0)
        self.assertEqual(str(d), "(0.0, 0.0)")

    def test_repr(self):
        """Check the representation
        """
        d = SCoordinate(0, 0)
        self.assertEqual(repr(d), "(0.0, 0.0)")


class TestCCoordinate(unittest.TestCase):

    """Test SCoordinates.
    """

    def test_str(self):
        """Check the string representation
        """
        d = CCoordinate((0, 0, 0))
        self.assertEqual(str(d), "[0 0 0]")

    def test_repr(self):
        """Check the representation
        """
        d = CCoordinate((0, 0, 0))
        self.assertEqual(repr(d), "[0 0 0]")

    def test_norm(self):
        """Euclidean norm of a cartesian vector
        """
        d = CCoordinate((1, 0, 0))
        self.assertEqual(d.norm(), 1.0)

    def test_normalize(self):
        """Normalize a cartesian vector
        """
        d = CCoordinate((2., 0., 0.))
        self.assertTrue(np.allclose(d.normalize().cart, [1, 0, 0]))

    def test_cross(self):
        """Test cross product in cartesian coordinates
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        self.assertTrue(np.allclose(d.cross(c).cart, [0., 0., 1.]))

    def test_dot(self):
        """Test the dot product of two cartesian vectors.
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        self.assertEqual(d.dot(c), 0)

    def test_ne(self):
        """Test inequality of two cartesian vectors.
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        self.assertTrue(c != d)

    def test_eq(self):
        """Test equality of two cartesian vectors.
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        self.assertFalse(c == d)

    def test_add(self):
        """Test adding cartesian vectors.
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        b = CCoordinate((1., 1., 0.))
        self.assertTrue(np.allclose((d + c).cart, b.cart))

        self.assertTrue(np.allclose((d + (0, 1, 0)).cart, b.cart))

        self.assertTrue(np.allclose(((0, 1, 0) + d).cart, b.cart))

    def test_mul(self):
        """Test multiplying (element-wise) cartesian vectors.
        """
        d = CCoordinate((1., 0., 0.))
        c = CCoordinate((0., 1., 0.))
        b = CCoordinate((0., 0., 0.))
        self.assertTrue(np.allclose((d * c).cart, b.cart))
        self.assertTrue(np.allclose((d * (0, 1, 0)).cart, b.cart))

        self.assertTrue(np.allclose(((0, 1, 0) * d).cart, b.cart))

    def test_to_spherical(self):
        """Test converting to spherical coordinates.
        """
        d = CCoordinate((1., 0., 0.))
        c = SCoordinate(0, 0)
        self.assertEqual(d.to_spherical(), c)


class TestArc(unittest.TestCase):

    """Test arcs
    """

    def test_eq(self):
        arc1 = Arc(SCoordinate(0, 0),
                   SCoordinate(np.deg2rad(10), np.deg2rad(10)))
        arc2 = Arc(SCoordinate(0, np.deg2rad(10)),
                   SCoordinate(np.deg2rad(10), 0))

        self.assertFalse(arc1 == arc2)

        self.assertTrue(arc1 == arc1)

    def test_ne(self):
        arc1 = Arc(SCoordinate(0, 0),
                   SCoordinate(np.deg2rad(10), np.deg2rad(10)))
        arc2 = Arc(SCoordinate(0, np.deg2rad(10)),
                   SCoordinate(np.deg2rad(10), 0))

        self.assertTrue(arc1 != arc2)

        self.assertFalse(arc1 != arc1)

    def test_str(self):
        arc1 = Arc(SCoordinate(0, 0),
                   SCoordinate(np.deg2rad(10), np.deg2rad(10)))
        self.assertEqual(str(arc1), str(arc1.start) + " -> " + str(arc1.end))
        self.assertEqual(repr(arc1), str(arc1.start) + " -> " + str(arc1.end))

    def test_intersection(self):
        arc1 = Arc(SCoordinate(0, 0),
                   SCoordinate(np.deg2rad(10), np.deg2rad(10)))
        arc2 = Arc(SCoordinate(0, np.deg2rad(10)),
                   SCoordinate(np.deg2rad(10), 0))
        lon, lat = arc1.intersection(arc2)

        self.assertTrue(np.allclose(np.rad2deg(lon), 5))
        self.assertEquals(np.rad2deg(lat), 5.0575148968282093)

        arc1 = Arc(SCoordinate(0, 0),
                   SCoordinate(np.deg2rad(10), np.deg2rad(10)))

        self.assertTrue(arc1.intersection(arc1) is None)

        arc1 = Arc(SCoordinate(np.deg2rad(24.341215776575297),
                               np.deg2rad(44.987819588259327)),
                   SCoordinate(np.deg2rad(18.842727517611817),
                               np.deg2rad(46.512483610284178)))
        arc2 = Arc(SCoordinate(np.deg2rad(20.165961750361905),
                               np.deg2rad(46.177305385810541)),
                   SCoordinate(np.deg2rad(20.253297585831707),
                               np.deg2rad(50.935830837274324)))
        inter = SCoordinate(np.deg2rad(20.165957021925202),
                            np.deg2rad(46.177022633103398))
        self.assertEquals(arc1.intersection(arc2), inter)

        arc1 = Arc(SCoordinate(np.deg2rad(-2.4982818108326734),
                               np.deg2rad(48.596644847869655)),
                   SCoordinate(np.deg2rad(-2.9571441235622835),
                               np.deg2rad(49.165688435261394)))
        arc2 = Arc(SCoordinate(np.deg2rad(-3.4976667413531688),
                               np.deg2rad(48.562704872921373)),
                   SCoordinate(np.deg2rad(-5.893976312685715),
                               np.deg2rad(48.445795283217116)))

        self.assertTrue(arc1.intersection(arc2) is None)

    def test_angle(self):
        arc1 = Arc(SCoordinate(np.deg2rad(157.5),
                               np.deg2rad(89.234600944314138)),
                   SCoordinate(np.deg2rad(90),
                               np.deg2rad(89)))
        arc2 = Arc(SCoordinate(np.deg2rad(157.5),
                               np.deg2rad(89.234600944314138)),
                   SCoordinate(np.deg2rad(135),
                               np.deg2rad(89)))

        self.assertAlmostEqual(np.rad2deg(arc1.angle(arc2)), -44.996385007218926, 13)

        arc1 = Arc(SCoordinate(np.deg2rad(112.5),
                               np.deg2rad(89.234600944314138)),
                   SCoordinate(np.deg2rad(90), np.deg2rad(89)))
        arc2 = Arc(SCoordinate(np.deg2rad(112.5),
                               np.deg2rad(89.234600944314138)),
                   SCoordinate(np.deg2rad(45), np.deg2rad(89)))

        self.assertAlmostEqual(np.rad2deg(arc1.angle(arc2)), 44.996385007218883, 13)

        arc1 = Arc(SCoordinate(0, 0), SCoordinate(1, 0))
        self.assertEqual(arc1.angle(arc1), 0)

        arc2 = Arc(SCoordinate(1, 0), SCoordinate(0, 0))
        self.assertEqual(arc1.angle(arc2), 0)

        arc2 = Arc(SCoordinate(0, 0), SCoordinate(-1, 0))
        self.assertEqual(arc1.angle(arc2), np.pi)

        arc2 = Arc(SCoordinate(2, 0), SCoordinate(1, 0))
        self.assertEqual(arc1.angle(arc2), np.pi)

        arc2 = Arc(SCoordinate(2, 0), SCoordinate(3, 0))
        self.assertRaises(ValueError, arc1.angle, arc2)


class TestSphericalPolygon(unittest.TestCase):

    """Test the spherical polygon.
    """

    def test_area(self):
        """Test the area function
        """
        vertices = np.array([[1, 2, 3, 4, 3, 2],
                             [3, 4, 3, 2, 1, 2]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.00121732523118, polygon.area())

        vertices = np.array([[1, 2, 3, 2],
                             [3, 4, 3, 2]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.000608430665842, polygon.area())

        vertices = np.array([[0, 0, 1, 1],
                             [0, 1, 1, 0]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.000304609684862, polygon.area())

        # Across the dateline

        vertices = np.array([[179.5, -179.5, -179.5, 179.5],
                             [1, 1, 0, 0]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.000304609684862, polygon.area())

        vertices = np.array([[0, 90, 90, 0],
                             [1, 1, 0, 0]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.0349012696772, polygon.area())

        vertices = np.array([[90, 0, 0],
                             [0, 0, 90]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(np.pi / 2, polygon.area())

        # Around the north pole

        vertices = np.array([[0, -90, 180, 90],
                             [89, 89, 89, 89]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.000609265770322, polygon.area())

        # Around the south pole

        vertices = np.array([[0, 90, 180, -90],
                             [-89, -89, -89, -89]]).T
        polygon = SphPolygon(np.deg2rad(vertices))

        self.assertAlmostEqual(0.000609265770322, polygon.area())

    def test_is_inside(self):
        """Test checking if a polygon is inside of another.
        """

        vertices = np.array([[1, 1, 20, 20],
                             [1, 20, 20, 1]]).T

        polygon1 = SphPolygon(np.deg2rad(vertices))

        vertices = np.array([[0, 0, 30, 30],
                             [0, 30, 30, 0]]).T

        polygon2 = SphPolygon(np.deg2rad(vertices))

        self.assertTrue(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))
        self.assertTrue(polygon2.area() > polygon1.area())

        polygon2.invert()
        self.assertFalse(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))

        vertices = np.array([[0, 0, 30, 30],
                             [21, 30, 30, 21]]).T

        polygon2 = SphPolygon(np.deg2rad(vertices))
        self.assertFalse(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))

        polygon2.invert()

        self.assertTrue(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))

        vertices = np.array([[100, 100, 130, 130],
                             [41, 50, 50, 41]]).T

        polygon2 = SphPolygon(np.deg2rad(vertices))

        self.assertFalse(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))

        polygon2.invert()

        self.assertTrue(polygon1._is_inside(polygon2))
        self.assertFalse(polygon2._is_inside(polygon1))

        vertices = np.array([[-1.54009253, 82.62402855],
                             [3.4804808, 82.8105746],
                             [20.7214892, 83.00875812],
                             [32.8857629, 82.7607758],
                             [41.53844302, 82.36024339],
                             [47.92062759, 81.91317164],
                             [52.82785062, 81.45769791],
                             [56.75107895, 81.00613046],
                             [59.99843787, 80.56042986],
                             [62.76998034, 80.11814453],
                             [65.20076209, 79.67471372],
                             [67.38577498, 79.22428],
                             [69.39480149, 78.75981318],
                             [71.28163984, 78.27283234],
                             [73.09016378, 77.75277976],
                             [74.85864685, 77.18594725],
                             [76.62327682, 76.55367303],
                             [78.42162204, 75.82918893],
                             [80.29698409, 74.97171721],
                             [82.30538638, 73.9143231],
                             [84.52973107, 72.53535661],
                             [87.11696138, 70.57600156],
                             [87.79163209, 69.98712409],
                             [72.98142447, 67.1760143],
                             [61.79517279, 63.2846272],
                             [53.50600609, 58.7098766],
                             [47.26725347, 53.70533139],
                             [42.44083259, 48.42199571],
                             [38.59682041, 42.95008531],
                             [35.45189206, 37.3452509],
                             [32.43435578, 30.72373327],
                             [31.73750748, 30.89485287],
                             [29.37284023, 31.44344415],
                             [27.66001308, 31.81016309],
                             [26.31358296, 32.08057499],
                             [25.1963477, 32.29313986],
                             [24.23118049, 32.46821821],
                             [23.36993508, 32.61780082],
                             [22.57998837, 32.74952569],
                             [21.8375532, 32.86857867],
                             [21.12396693, 32.97868717],
                             [20.42339605, 33.08268331],
                             [19.72121983, 33.18284728],
                             [19.00268283, 33.28113306],
                             [18.2515215, 33.3793305],
                             [17.4482606, 33.47919405],
                             [16.56773514, 33.58255576],
                             [15.57501961, 33.6914282],
                             [14.4180087, 33.8080799],
                             [13.01234319, 33.93498577],
                             [11.20625437, 34.0742239],
                             [8.67990371, 34.22415978],
                             [7.89344478, 34.26018768],
                             [8.69446485, 41.19823568],
                             [9.25707165, 47.17351118],
                             [9.66283477, 53.14128114],
                             [9.84134875, 59.09937166],
                             [9.65054241, 65.04458004],
                             [8.7667375, 70.97023122],
                             [6.28280904, 76.85731403]])
        polygon1 = SphPolygon(np.deg2rad(vertices))

        vertices = np.array([[49.94506701, 46.52610743],
                             [51.04293649, 46.52610743],
                             [62.02163129, 46.52610743],
                             [73.0003261, 46.52610743],
                             [83.9790209, 46.52610743],
                             [85.05493299, 46.52610743],
                             [85.05493299, 45.76549301],
                             [85.05493299, 37.58315571],
                             [85.05493299, 28.39260587],
                             [85.05493299, 18.33178739],
                             [85.05493299, 17.30750918],
                             [83.95706351, 17.30750918],
                             [72.97836871, 17.30750918],
                             [61.9996739, 17.30750918],
                             [51.0209791, 17.30750918],
                             [49.94506701, 17.30750918],
                             [49.94506701, 18.35262921],
                             [49.94506701, 28.41192025],
                             [49.94506701, 37.60055422],
                             [49.94506701, 45.78080831]])
        polygon2 = SphPolygon(np.deg2rad(vertices))

        self.assertFalse(polygon2._is_inside(polygon1))
        self.assertFalse(polygon1._is_inside(polygon2))

    def test_bool(self):
        """Test the intersection and union functions.
        """
        vertices = np.array([[180, 90, 0, -90],
                             [89, 89, 89, 89]]).T
        poly1 = SphPolygon(np.deg2rad(vertices))
        vertices = np.array([[-45, -135, 135, 45],
                             [89, 89, 89, 89]]).T
        poly2 = SphPolygon(np.deg2rad(vertices))

        uni = np.array([[157.5,   89.23460094],
                        [-225.,   89.],
                        [112.5,   89.23460094],
                        [90.,   89.],
                        [67.5,   89.23460094],
                        [45.,   89.],
                        [22.5,   89.23460094],
                        [0.,   89.],
                        [-22.5,   89.23460094],
                        [-45.,   89.],
                        [-67.5,   89.23460094],
                        [-90.,   89.],
                        [-112.5,   89.23460094],
                        [-135.,   89.],
                        [-157.5,   89.23460094],
                        [-180.,   89.]])
        inter = np.array([[157.5,   89.23460094],
                          [112.5,   89.23460094],
                          [67.5,   89.23460094],
                          [22.5,   89.23460094],
                          [-22.5,   89.23460094],
                          [-67.5,   89.23460094],
                          [-112.5,   89.23460094],
                          [-157.5,   89.23460094]])
        poly_inter = poly1.intersection(poly2)
        poly_union = poly1.union(poly2)

        self.assertTrue(poly_inter.area() <= poly_union.area())

        self.assertTrue(np.allclose(poly_inter.vertices,
                                    np.deg2rad(inter)))
        self.assertTrue(np.allclose(poly_union.vertices,
                                    np.deg2rad(uni)))

        # Test 2 polygons sharing 2 contiguous edges.

        vertices1 = np.array([[-10,  10],
                              [-5,  10],
                              [0,  10],
                              [5,  10],
                              [10,  10],
                              [10, -10],
                              [-10, -10]])

        vertices2 = np.array([[-5,  10],
                              [0,  10],
                              [5,  10],
                              [5,  -5],
                              [-5,  -5]])

        vertices3 = np.array([[5,  10],
                              [5,  -5],
                              [-5,  -5],
                              [-5,  10],
                              [0,  10]])

        poly1 = SphPolygon(np.deg2rad(vertices1))
        poly2 = SphPolygon(np.deg2rad(vertices2))
        poly_inter = poly1.intersection(poly2)

        self.assertTrue(np.allclose(poly_inter.vertices,
                                    np.deg2rad(vertices3)))

        # Test when last node of the intersection is the last vertice of the
        # second polygon.

        swath_vertices = np.array([[-115.32268301,   66.32946139],
                                   [-61.48397172,  58.56799254],
                                   [-60.25004314, 58.00754686],
                                   [-71.35057076,   49.60229517],
                                   [-113.746486,  56.03008985]])
        area_vertices = np.array([[-68.32812107,  52.3480829],
                                  [-67.84993896,  53.07015692],
                                  [-55.54651296,  64.9254637],
                                  [-24.63341856,  74.24628796],
                                  [-31.8996363,  27.99907764],
                                  [-39.581043,  37.0639821],
                                  [-50.90185988,  45.56296169],
                                  [-67.43022017,  52.12399581]])

        res = np.array([[-62.77837918,   59.12607053],
                        [-61.48397172,   58.56799254],
                        [-60.25004314,   58.00754686],
                        [-71.35057076,   49.60229517],
                        [-113.746486,     56.03008985],
                        [-115.32268301,   66.32946139]])

        poly1 = SphPolygon(np.deg2rad(swath_vertices))
        poly2 = SphPolygon(np.deg2rad(area_vertices))

        poly_inter = poly1.intersection(poly2)
        self.assertTrue(np.allclose(poly_inter.vertices,
                                    np.deg2rad(res)))

        poly_inter = poly2.intersection(poly1)
        self.assertTrue(np.allclose(poly_inter.vertices,
                                    np.deg2rad(res)))

        # vertices = np.array([[ -84.54058691,   71.80094043],
        #                      [ -74.68557932,   72.16812631],
        #                      [ -68.06987203,   72.1333064 ],
        #                      [ -63.17961469,   71.96265   ],
        #                      [ -59.33392061,   71.73824792],
        #                      [ -56.16798418,   71.49047832],
        #                      [ -53.46489053,   71.231076  ],
        #                      [ -51.08551155,   70.96395329],
        #                      [ -48.93484325,   70.68929276],
        #                      [ -46.94415494,   70.40519826],
        #                      [ -45.06071892,   70.10832093],
        #                      [ -43.24140861,   69.7939738 ],
        #                      [ -41.44830671,   69.45591086],
        #                      [ -39.64527217,   69.08578252],
        #                      [ -37.79474271,   68.6721527 ],
        #                      [ -35.85408829,   68.1987858 ],
        #                      [ -33.7705704 ,   67.64156121],
        #                      [ -31.47314483,   66.9625364 ],
        #                      [ -28.85703847,   66.09736791],
        #                      [ -25.74961912,   64.92465312],
        #                      [ -21.81516555,   63.17261421],
        #                      [ -18.62398733,   62.28633798],
        #                      [ -16.93359509,   62.89011263],
        #                      [ -15.17161807,   63.47161418],
        #                      [ -13.33621801,   64.02936211],
        #                      [ -11.42593772,   64.56180886],
        #                      [  -9.43979715,   65.0673476 ],
        #                      [  -7.37739816,   65.54432277],
        #                      [  -5.23903263,   65.99104411],
        #                      [  -3.02579085,   66.40580433],
        #                      [  -0.73966571,   66.78690012],
        #                      [   1.61635637,   67.13265703],
        #                      [   4.03822468,   67.44145758],
        #                      [   6.52078043,   67.71177166],
        #                      [   9.05775043,   67.94218891],
        #                      [  11.64178394,   68.13145134],
        #                      [  14.26453542,   68.27848476],
        #                      [  16.9167971 ,   68.38242749],
        #                      [  19.58867724,   68.44265471],
        #                      [  22.26981526,   68.45879658],
        #                      [  24.94962586,   68.43074943],
        #                      [  27.61755654,   68.35867876],
        #                      [  30.26334172,   68.24301426],
        #                      [  32.87724117,   68.08443684],
        #                      [  35.45024798,   67.88385879],
        #                      [  37.97425437,   67.64239838],
        #                      [  40.44217258,   67.36135027],
        #                      [  42.84800609,   67.04215364],
        #                      [  45.18687531,   66.68635947],
        #                      [  47.45500013,   66.2955988 ],
        #                      [  49.64965026,   65.87155246],
        #                      [  52.34514841,   66.28428851],
        #                      [  56.04377347,   68.57914951],
        #                      [  59.05474396,   70.10401937],
        #                      [  61.66799965,   71.23110288],
        #                      [  64.02929638,   72.12002156],
        #                      [  66.22835251,   72.85391032],
        #                      [  68.32829893,   73.48143318],
        #                      [  70.37866226,   74.03347161],
        #                      [  72.42237212,   74.53085444],
        #                      [  74.50035309,   74.98833047],
        #                      [  76.65524775,   75.41675945],
        #                      [  78.93517067,   75.824363  ],
        #                      [  81.39826053,   76.21741056],
        #                      [  84.11897279,   76.600482  ],
        #                      [  87.19757467,   76.97627542],
        #                      [  90.77537201,   77.3447072 ],
        #                      [  95.06035831,   77.70058684],
        #                      [ 100.37229526,   78.02797258],
        #                      [ 107.22498444,   78.28582497],
        #                      [ 116.481466  ,   78.36746171],
        #                      [ 129.66805239,   77.96163057],
        #                      [ 134.67038545,   78.4115401 ],
        #                      [ 136.40302873,   79.30544125],
        #                      [ 138.4763311 ,   80.18558961],
        #                      [ 140.98282558,   81.04796485],
        #                      [ 144.04700981,   81.88693584],
        #                      [ 147.83664747,   82.6944745 ],
        #                      [ 152.57512293,   83.45896996],
        #                      [ 158.54810167,   84.16352558],
        #                      [ 166.0844409 ,   84.78383882],
        #                      [ 175.46720475,   85.28657382],
        #                      [-173.27937931,   85.6309921 ],
        #                      [-160.67741256,   85.77820349],
        #                      [-147.84352095,   85.70789809],
        #                      [-136.01435526,   85.4301266 ],
        #                      [-125.94447471,   84.97922118],
        #                      [-117.77450148,   84.39683471],
        #                      [-111.28213275,   83.71944226],
        #                      [-106.1391311 ,   82.97447237],
        #                      [-102.03983076,   82.18121521],
        #                      [ -98.73868716,   81.3529452 ],
        #                      [ -96.04944891,   80.49880811],
        #                      [ -93.83359781,   79.62518236],
        #                      [ -91.98834044,   78.73659234],
        #                      [ -90.43691725,   77.83630659],
        #                      [ -89.12142407,   76.92672961],
        #                      [ -87.99766337,   76.0096614 ],
        #                      [ -87.03148527,   75.08647127],
        #                      [ -86.19618441,   74.15821627],
        #                      [ -85.47063566,   73.22572391],
        #                      [ -84.83794555,   72.28964996]])
        # polygon = SphPolygon(np.deg2rad(vertices))
        # polygon.invert()
        # from datetime import datetime
        # utctime = datetime(2013, 12, 12, 9, 31, 54, 485719)
        # utctime = datetime(2013, 11, 11, 11, 11)
        # twi = get_twilight_poly(utctime)
        # poly_inter_day = twi.intersection(polygon)
        # twi.invert()
        # poly_inter_night = twi.intersection(polygon)
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.basemap import Basemap
        # map = Basemap(projection='nsper', lat_0 = 58, lon_0 = 16,
        #               resolution = 'l', area_thresh = 1000.)
        # map = Basemap(resolution = "l")
        # map.drawcoastlines()
        # map.drawcountries()
        # map.drawmapboundary(fill_color='white')
        # map.drawmeridians(np.arange(0, 360, 30))
        # map.drawparallels(np.arange(-90, 90, 30))

        # poly_inter_day.draw(map, "-r")
        # poly_inter_night.draw(map, "-b")
        # plt.show()

    # def test_twilight(self):
    #     """Test the twilight polygon.
    #     """
    #     from datetime import datetime
    #     utctime = datetime(2013, 3, 20, 12, 0)

    # print np.rad2deg(get_twilight_poly(utctime).vertices)

    #     vertices = np.array([[0, -90, 180, 90],
    #                          [89, 89, 89, 89]]).T


def suite():
    """The suite for test_spherical
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestSCoordinate))
    mysuite.addTest(loader.loadTestsFromTestCase(TestCCoordinate))
    mysuite.addTest(loader.loadTestsFromTestCase(TestArc))
    mysuite.addTest(loader.loadTestsFromTestCase(TestSphericalPolygon))

    return mysuite


if __name__ == '__main__':
    unittest.main()
