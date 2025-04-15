import unittest
import numpy as np
from gecatsim.pyfiles.rimg import rimg

class TestRimg(unittest.TestCase):

    def test_rimg_default(self):
        result = rimg(5)
        expected = np.array([
            [2.82842712, 2.23606798, 2.0, 2.23606798, 2.82842712],
            [2.23606798, 1.41421356, 1.0, 1.41421356, 2.23606798],
            [2.0, 1.0, 0.0, 1.0, 2.0],
            [2.23606798, 1.41421356, 1.0, 1.41421356, 2.23606798],
            [2.82842712, 2.23606798, 2.0, 2.23606798, 2.82842712]
        ])
        np.testing.assert_almost_equal(result, expected, decimal=6)

    def test_rimg_custom_center(self):
        result = rimg(5, 5, 1, 3, 3)
        expected = np.array([
            [2.82842712, 2.23606798, 2.0, 2.23606798, 2.82842712],
            [2.23606798, 1.41421356, 1.0, 1.41421356, 2.23606798],
            [2.0, 1.0, 0.0, 1.0, 2.0],
            [2.23606798, 1.41421356, 1.0, 1.41421356, 2.23606798],
            [2.82842712, 2.23606798, 2.0, 2.23606798, 2.82842712]
        ])
        np.testing.assert_almost_equal(result, expected, decimal=6)

    def test_rimg_non_square(self):
        result = rimg(4, 3)
        expected = np.array([
            [1.80277564, 1.11803399, 1.11803399, 1.80277564],
            [1.5, 0.5, 0.5, 1.5],
            [1.80277564, 1.11803399, 1.11803399, 1.80277564]
        ])
        np.testing.assert_almost_equal(result, expected, decimal=6)

if __name__ == '__main__':
    unittest.main()