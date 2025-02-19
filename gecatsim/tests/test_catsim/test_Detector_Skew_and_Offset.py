import unittest
import numpy as np
from gecatsim.pyfiles.Detector_Skew_and_Offset import detector_skew_and_offset
from gecatsim.pyfiles.CommonTools import CFG

class TestDetectorSkewAndOffset(unittest.TestCase):

    def setUp(self):
        cfg = CFG("../examples/cfg/Phantom_Sample",
                  "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")
        cfg.individual_cell_offset = 1
        cfg.col_count = 2
        cfg.row_count = 2
        cfg.individual_pack_skew = 1
        cfg.col_packs_per_mod = 1
        cfg.cols_per_pack = 2
        cfg.rows_per_pack = 2
        cfg.col_size = 1
        cfg.col_pack_gap = 0
        cfg.individual_module_skew = 1
        cfg.module_skew_w = [0]
        cfg.pack_skew_w = [0]  # 1 module with 1 pack

        self.cfg = cfg

        det = type('', (), {})()  # Create a simple object to hold attributes
        det.n_cells = 4
        det.n_modules = 1
        det.modtypes = [0]  # Ensure modtypes match the dimension of cellcoords
        det.cellcoords = np.zeros((2, 4, 1))
        det.n_moddefs = 1
        det.n_samples = np.zeros((1, 1))
        det.samplecoords = np.zeros((2, 1, 1))
        det.weights = np.zeros((1, 1, 1))
        det.activearea = np.zeros((1, 1))
        det.width = np.zeros((1, 1))
        det.height = np.zeros((1, 1))
        det.uvecs = np.zeros((3, 1))
        det.vvecs = np.zeros((3, 1))

        self.det = det

    def test_detector_skew_and_offset(self):
        updated_det = detector_skew_and_offset(self.cfg, self.det)
        self.assertIsNotNone(updated_det)
        self.assertTrue(hasattr(updated_det, 'cellcoords'))
        self.assertTrue(hasattr(updated_det, 'uvecs'))
        self.assertTrue(hasattr(updated_det, 'vvecs'))

if __name__ == '__main__':
    unittest.main()