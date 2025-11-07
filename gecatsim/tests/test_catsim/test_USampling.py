import unittest
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.USampling import USampling

cfg = CFG("../examples/cfg/Phantom_Sample",
          "../examples/cfg/Scanner_Sample_generic",
          "../examples/cfg/Protocol_Sample_axial")

cfg.scanner.detectorColSize = 1.0
cfg.scanner.colIntensiveOversampleLength = 0.1
cfg.scanner.plateThickness = 0.0
cfg.scanner.plateHeight = 0.0
cfg.scanner.plateAirgap = 0.0
cfg.scanner.sdd = 1000.0
cfg.scanner.colCast = 0.0
cfg.scanner.detectionOnKerf = 0
cfg.scanner.callbackDetector = 'Detector_SVCT'
cfg.scanner.plateLocation = 0

class TestUSampling(unittest.TestCase):

    def test_usampling_no_intensive(self):
        total_n = 10
        intensive_n = 0
        us, steplen = USampling(cfg, total_n, intensive_n)

        # Check the length of the output arrays
        self.assertEqual(len(us), total_n)
        self.assertEqual(len(steplen), total_n)

        # Check the values of the output arrays
        self.assertTrue(np.allclose(steplen, steplen[0]))  # All step lengths should be equal

    def test_usampling_with_intensive(self):
        total_n = 20
        intensive_n = 5
        us, steplen = USampling(cfg, total_n, intensive_n)

        # Check the length of the output arrays
        self.assertEqual(len(us), total_n)
        self.assertEqual(len(steplen), total_n)

        # Check the values of the output arrays
        self.assertTrue(np.allclose(steplen[:intensive_n], steplen[0]))  # Intensive region step lengths should be equal
        self.assertTrue(np.allclose(steplen[-intensive_n:], steplen[0]))  # Intensive region step lengths should be equal

if __name__ == '__main__':
    unittest.main()