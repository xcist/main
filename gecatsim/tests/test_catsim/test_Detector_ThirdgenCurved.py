from gecatsim.pyfiles.Detector_ThirdgenCurved import Detector_ThirdgenCurved
import unittest.mock
from unittest.mock import patch
from gecatsim.pyfiles.CommonTools import *
import numpy as np

class test_Detector_ThirdgenCurved(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_Detector_ThirdgenCurved(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.scanner.sid = 541
        cfg.scanner.sdd = 949
        cfg.scanner.detectorRowsPerMod = 16
        cfg.scanner.detectorColsPerMod = 16
        cfg.scanner.detectorColOffset = -1.25
        cfg.scanner.detectorRowOffset = 0
        cfg.scanner.detectorColSize = 1.0239
        cfg.scanner.detectorRowSize = 1.096349
        cfg.scanner.detectorColFillFraction = 0.8
        cfg.scanner.detectorRowFillFraction = 0.8
        cfg.scanner.detectorColCount = 880
        cfg.scanner.detectorRowCount = cfg.scanner.detectorRowsPerMod

        cfg.physics.colSampleCount = 4
        cfg.physics.rowSampleCount = 3

        cfg = Detector_ThirdgenCurved(cfg)

        assert cfg.det.cellCoords is not None
        assert cfg.det.sampleCoords is not None
        assert cfg.det.modCoords is not None
        assert cfg.det.uvecs is not None
        assert cfg.det.vvecs is not None