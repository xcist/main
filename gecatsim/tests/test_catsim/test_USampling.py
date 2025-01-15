import unittest
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.USampling import USampling

class TestUSampling(unittest.TestCase):

    def setUp(self):
        self.cfg = CFG("../examples/cfg/Phantom_Sample",
                       "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")
        # Add missing attributes
        self.cfg.recon_slice_thickness = 1.0
        self.cfg.recon_mu = 0.02
        self.cfg.col_size = 1.0
        self.cfg.col_intensive_oversample_length = 0.1
        self.cfg.plate_thickness = 0.0
        self.cfg.plate_height = 0.0
        self.cfg.plate_airgap = 0.0
        self.cfg.sdd = 1000.0
        self.cfg.col_cast = 0.0
        self.cfg.detection_on_kerf = 0
        self.cfg.callback_detector = 'Detector_SVCT'
        self.cfg.plate_location = 0

    def test_usampling_no_intensive(self):
        total_n = 10
        intensive_n = 0
        us, steplen = USampling(self.cfg, total_n, intensive_n)

        # Check the length of the output arrays
        self.assertEqual(len(us), total_n)
        self.assertEqual(len(steplen), total_n)

        # Check the values of the output arrays
        self.assertTrue(np.allclose(steplen, steplen[0]))  # All step lengths should be equal

    def test_usampling_with_intensive(self):
        total_n = 20
        intensive_n = 5
        us, steplen = USampling(self.cfg, total_n, intensive_n)

        # Check the length of the output arrays
        self.assertEqual(len(us), total_n)
        self.assertEqual(len(steplen), total_n)

        # Check the values of the output arrays
        self.assertTrue(np.allclose(steplen[:intensive_n], steplen[0]))  # Intensive region step lengths should be equal
        self.assertTrue(
            np.allclose(steplen[-intensive_n:], steplen[0]))  # Intensive region step lengths should be equal


if __name__ == '__main__':
    unittest.main()