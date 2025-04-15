import unittest
import numpy as np
from gecatsim.pyfiles.find_R import find_R

class TestFindR(unittest.TestCase):

    def test_find_R(self):
        x1 = np.array([1, 0, 0])
        x2 = np.array([0, 1, 0])
        y1 = np.array([0, 1, 0])
        y2 = np.array([-1, 0, 0])
        expected_R = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        R = find_R(x1, x2, y1, y2)
        np.testing.assert_almost_equal(R, expected_R, decimal=5)

if __name__ == '__main__':
    unittest.main()