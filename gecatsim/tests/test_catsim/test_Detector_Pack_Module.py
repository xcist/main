import unittest
import numpy as np
from unittest.mock import patch
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.Detector_Pack_Module import detector_pack_module

class TestDetectorPackModule(unittest.TestCase):

    @patch('gecatsim.pyfiles.USampling.USampling')
    def test_detector_pack_module(self, mock_USampling):
        cfg = CFG("../examples/cfg/Phantom_Sample",
                  "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.scanner.detectorColSize = 1.0
        cfg.scanner.detectorRowSize = 1.0
        cfg.scanner.detectorColsPerMod = 4
        cfg.scanner.detectorRowsPerMod = 4
        cfg.scanner.colPacksPerMod = 2
        cfg.scanner.rowPacksPerMod = 2
        cfg.scanner.colPackGap = 0.1
        cfg.scanner.rowPackGap = 0.1
        cfg.scanner.colOversample = 4
        cfg.scanner.colIntensiveOversample = 1
        cfg.scanner.colIntensiveOversampleLength = 0.2
        cfg.scanner.colCrosstalk = 0.1
        cfg.scanner.rowOversample = 4
        cfg.scanner.rowCast = 0.1
        cfg.scanner.colCast = 0.0
        cfg.scanner.rowCrosstalk = 0.1
        cfg.scanner.callbackDetector = 'Detector_SVCT'

        det = type('', (), {})()  # Create a simple object to hold attributes
        mock_USampling.return_value = (np.array([0.25, 0.75]), np.array([0.5, 0.5]))

        result = detector_pack_module(cfg, det)

        self.assertEqual(result.n_cells, 16)
        self.assertEqual(result.n_samples, 36)
        self.assertIsNotNone(result.cellcoords)
        self.assertIsNotNone(result.samplecoords)
        self.assertIsNotNone(result.weights)
        self.assertIsNotNone(result.activearea)
        self.assertIsNotNone(result.width)
        self.assertIsNotNone(result.height)
        self.assertIsNotNone(result.cols)
        self.assertIsNotNone(result.rows)
        self.assertIsNotNone(result.sample_du)
        self.assertIsNotNone(result.sample_dv)

if __name__ == '__main__':
    unittest.main()