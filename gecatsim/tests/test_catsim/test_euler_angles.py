import unittest
import numpy as np
from gecatsim.pyfiles.euler_angles import euler_angles

class TestEulerAngles(unittest.TestCase):

    def test_euler_angles(self):
        R = np.array([
            [np.cos(0.5) * np.cos(0.3) - np.sin(0.5) * np.cos(0.2) * np.sin(0.3), np.cos(0.5) * np.sin(0.3) + np.sin(0.5) * np.cos(0.2) * np.cos(0.3), np.sin(0.5) * np.sin(0.2)],
            [-np.sin(0.5) * np.cos(0.3) - np.cos(0.5) * np.cos(0.2) * np.sin(0.3), -np.sin(0.5) * np.sin(0.3) + np.cos(0.5) * np.cos(0.2) * np.cos(0.3), np.cos(0.5) * np.sin(0.2)],
            [np.sin(0.2) * np.sin(0.3), -np.sin(0.2) * np.cos(0.3), np.cos(0.2)]
        ])
        expected_ea = np.array([0.5, 0.2, 0.3])
        ea = euler_angles(R)
        np.testing.assert_almost_equal(ea, expected_ea, decimal=5)

if __name__ == '__main__':
    unittest.main()