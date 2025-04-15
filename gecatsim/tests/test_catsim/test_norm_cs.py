import unittest
import numpy as np
from gecatsim.pyfiles.norm_cs import norm_cs

class TestNormCS(unittest.TestCase):

    def test_norm_cs_vector(self):
        vec = np.array([1, 2, 3])
        result = norm_cs(vec)
        expected = np.linalg.norm(vec)
        self.assertAlmostEqual(result, expected, places=6)

    def test_norm_cs_single_element(self):
        vec = np.array([5])
        result = norm_cs(vec)
        expected = np.linalg.norm(vec)
        self.assertAlmostEqual(result, expected, places=6)

    def test_norm_cs_2d_vector(self):
        vec = np.array([[1], [2], [3]])
        result = norm_cs(vec)
        expected = np.linalg.norm(vec)
        self.assertAlmostEqual(result, expected, places=6)

    def test_norm_cs_invalid_input(self):
        mat = np.array([[1, 2], [3, 4]])
        with self.assertRaises(ValueError):
            norm_cs(mat)

if __name__ == '__main__':
    unittest.main()