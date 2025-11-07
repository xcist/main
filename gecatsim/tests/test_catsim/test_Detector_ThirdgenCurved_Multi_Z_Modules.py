import unittest
from unittest.mock import patch
from gecatsim.pyfiles.Detector_ThirdgenCurved_Multi_Z_Modules import detector_thirdgen_curved_multi_z_modules
from gecatsim.pyfiles.CommonTools import CFG

class TestDetectorThirdgenCurvedMultiZModules(unittest.TestCase):

    @patch('gecatsim.pyfiles.Detector_Pack_Module.detector_pack_module')
    def test_detector_thirdgen_curved_multi_z_modules(self, mock_detector_pack_module):

        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.scanner.sid = 541
        cfg.scanner.sdd = 949

        cfg.scanner.detectorRowsPerMod = 16
        cfg.scanner.detectorColsPerMod = 16
        cfg.scanner.colOffset = 0.0
        cfg.scanner.modRowOffset = [0, 0]
        cfg.scanner.detectorColSize = 1.0239
        cfg.scanner.detectorRowSize = 1.096349
        cfg.scanner.modsPerDetX = 2
        cfg.scanner.modsPerDetZ = 2
        cfg.scanner.modGapU = 0.5
        cfg.scanner.modGapV = 0.5
        cfg.scanner.modOffsetW = 0.1
        cfg.scanner.modRotationU = 0.0
        cfg.scanner.colPackGap = 0.1
        cfg.scanner.colPacksPerMod = 1
        cfg.scanner.rowPacksPerMod = 1
        cfg.scanner.rowPackGap = 0.1
        cfg.scanner.modDeltaAlpha = 0.0
        cfg.scanner.callbackDetector = 'Detector_SVCT'

        cfg.physics.colSampleCount = 4
        cfg.physics.rowSampleCount = 3

        cfg.scanner.colOversample = 4
        cfg.scanner.rowOversample = 4
        cfg.scanner.colCrosstalk = 0.1
        cfg.scanner.rowCrosstalk = 0.1
        cfg.scanner.rowCast = 0.1
        cfg.scanner.colCast = 0.0

        det = detector_thirdgen_curved_multi_z_modules(cfg)

        self.assertIsNotNone(det.modcoords)
        self.assertIsNotNone(det.uvecs)
        self.assertIsNotNone(det.vvecs)
        self.assertEqual(det.n_modules, 4)
        self.assertEqual(det.total_n_cells, 1024)  # 16 * 16 * 4
        self.assertEqual(len(det.startindices), 4)
        self.assertEqual(det.n_moddefs, 1)
        self.assertTrue((det.modtypes == 1).all())

if __name__ == '__main__':
    unittest.main()