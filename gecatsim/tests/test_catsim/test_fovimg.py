import unittest
import numpy as np
from numpy.testing import assert_array_equal
from gecatsim.pyfiles.fovimg import fovimg

class TestFovimg(unittest.TestCase):

    def test_fovimg_default_radius_center(self):
        # Test case with default radius and center
        result = fovimg(5, 5)
        expected = np.array([
            [False,  True,  True,  True, False],
            [ True,  True,  True,  True,  True],
            [ True,  True,  True,  True,  True],
            [ True,  True,  True,  True,  True],
            [False,  True,  True,  True, False]
        ])
        assert_array_equal(result, expected)

    def test_fovimg_custom_radius_center(self):
        # Test case with custom radius and center
        result = fovimg(5, 5, radius=1, centerrow=3, centercol=3)
        expected = np.array([
            [False, False, False, False, False],
            [False, False,  True, False, False],
            [False,  True,  True,  True, False],
            [False, False,  True, False, False],
            [False, False, False, False, False]
        ])
        assert_array_equal(result, expected)

    def test_fovimg_3d_volume(self):
        # Test case with 3D volume (cylinder)
        result = fovimg(3, 3, nrplanes=2)
        expected = np.array([
            [[True, True], [ True,  True], [ True,  True]],
            [[True, True], [ True,  True], [ True,  True]],
            [[True, True], [ True,  True], [ True,  True]]
        ])
        print("Actual result for 3D volume:")
        print(result)
        assert_array_equal(result, expected)

    def test_fovimg_non_square(self):
        # Test case with non-square dimensions
        result = fovimg(3, 5)
        expected = np.array([
            [False,  True,  True,  True, False],
            [False,  True,  True,  True, False],
            [False,  True,  True,  True, False]
        ])
        assert_array_equal(result, expected)

if __name__ == '__main__':
    unittest.main()