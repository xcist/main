import unittest
from unittest.mock import patch
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.ModularFlatPanelDetector import modular_flat_panel_detector

class TestModularFlatPanelDetector(unittest.TestCase):

    def test_modular_flat_panel_detector(self):

        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.scanner.sid = 541
        cfg.scanner.sdd = 949

        cfg.scanner.detectorRowsPerMod = 16
        cfg.scanner.detectorColsPerMod = 16
        cfg.scanner.detectorColCount = 32
        cfg.scanner.colOffset = 0.0
        cfg.scanner.rowOffset = 0.0
        cfg.scanner.detectorColSize = 1.0239
        cfg.scanner.detectorRowSize = 1.096349
        cfg.scanner.colOversample = 4
        cfg.scanner.rowOversample = 4
        cfg.scanner.colFillFraction = 0.8
        cfg.scanner.rowFillFraction = 0.8
        cfg.scanner.colCrosstalk = 0.1
        cfg.scanner.rowCrosstalk = 0.1

        det = modular_flat_panel_detector(cfg)

        self.assertIsNotNone(det.modcoords)
        self.assertIsNotNone(det.uvecs)
        self.assertIsNotNone(det.vvecs)
        self.assertEqual(det.n_modules, 2)
        self.assertEqual(det.total_n_cells, 512)  # 32 * 16
        self.assertEqual(len(det.startindices), 2)
        self.assertEqual(det.n_moddefs, 1)
        self.assertTrue((det.modtypes == 1).all())

if __name__ == '__main__':
    unittest.main()