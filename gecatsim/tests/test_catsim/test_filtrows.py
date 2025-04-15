import unittest
import numpy as np
from gecatsim.pyfiles.filtrows import filtrows

class TestFiltrows(unittest.TestCase):

    def test_filtrows(self):
        x = np.array([[1, 2, 3], [4, 5, 6]])
        h = np.array([0.2, 0.5, 0.3])
        expected_y = np.array([[2.3, 1.8, 1.9], [5.3, 4.8, 4.9]])
        y = filtrows(x, h)
        np.testing.assert_almost_equal(y, expected_y, decimal=5)

if __name__ == '__main__':
    unittest.main()