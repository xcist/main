import unittest
import numpy as np
from gecatsim.pyfiles.vectornorm import vectornorm

class TestVectorNorm(unittest.TestCase):
    def test_vectornorm_valid_input(self):
        xyz = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        expected_norms = np.sqrt(np.sum(xyz * xyz, axis=0))
        actual_norms = vectornorm(xyz)
        np.testing.assert_array_equal(actual_norms, expected_norms)

    def test_vectornorm_invalid_input(self):
        xyz = np.array([[1, 2], [3, 4]])  # Invalid shape
        result = vectornorm(xyz)
        self.assertIsNone(result)

    def test_vectornorm_empty_input(self):
        xyz = np.array([[], [], []])  # Empty input
        expected_norms = np.array([])
        actual_norms = vectornorm(xyz)
        np.testing.assert_array_equal(actual_norms, expected_norms)

if __name__ == '__main__':
    unittest.main()