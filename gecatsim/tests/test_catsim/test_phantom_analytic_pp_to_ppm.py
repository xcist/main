import unittest
import os
from unittest.mock import MagicMock
import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.Phantom_Analytic_pp_to_ppm import *

class TestPhantomAnalyticConversion(unittest.TestCase):
    def setUp(self):
        # Create a temporary .pp file for testing
        self.test_filename = 'test_phantom'
        with open(f'{self.test_filename}.pp', 'w') as f:
            f.write("# Material Table\n")
            f.write("Material1\n")
            f.write("# Object Data\n")
            f.write("Object1\n")
            f.write("0 0 0 1\n")
            f.write("1 0 0 0\n")
            f.write("0 1 0 0\n")
            f.write("0 0 1 0\n")
            f.write("0 0 0 1\n")
            f.write("C 0 0 0 0\n")
            f.write("Material 1\n")
            f.write("Density 1.0\n")

        # Create a temporary file for testing read_text_lines2
        self.test_text_filename = 'test_file.tmp'
        with open(self.test_text_filename, 'w') as f:
            f.write("Line 1\nLine 2\nLine 3\n")

        # Create a temporary file for testing read_object
        self.test_object_filename = 'test_object.tmp'
        with open(self.test_object_filename, 'w') as f:
            f.write("type 1\n")
            f.write("0 0 0 1\n")
            f.write("1 0 0 0\n")
            f.write("0 1 0 0\n")
            f.write("0 0 1 0\n")
            f.write("0 0 0 1\n")
            f.write("Material 1\n")
            f.write("Density 1.0\n")

    def tearDown(self):
        # Remove the temporary files after testing
        os.remove(f'{self.test_filename}.pp')
        if os.path.exists(f'{self.test_filename}.ppm'):
            os.remove(f'{self.test_filename}.ppm')
        if os.path.exists(f'{self.test_filename}.tmp'):
            os.remove(f'{self.test_filename}.tmp')
        os.remove(self.test_text_filename)
        os.remove(self.test_object_filename)

    def test_phantom_analytic_pp_to_ppm(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")
        cfg.clib.TranslatePhantom_FORBILD_to_tmp = MagicMock()

        # Simulate the creation of the .tmp file
        with open(f'{self.test_filename}.tmp', 'w') as f:
            f.write("# Material Table\n")
            f.write("Material1\n")
            f.write("# Object Data\n")
            f.write("Object1\n")

        phantom_analytic_pp_to_ppm(cfg, self.test_filename, scale=1)
        self.assertTrue(os.path.exists(f'{self.test_filename}.ppm'))

    def test_read_text_lines2(self):
        expected_lines = ['Line 1', '\nLine 2', '\nLine 3']
        actual_lines = read_text_lines2(self.test_text_filename)
        self.assertEqual(actual_lines, expected_lines)

    def test_euler_angs(self):
        R = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ]
        expected_angles = [0, 0, 0]
        actual_angles = euler_angs(R)
        self.assertEqual(actual_angles, expected_angles)

if __name__ == '__main__':
    unittest.main()