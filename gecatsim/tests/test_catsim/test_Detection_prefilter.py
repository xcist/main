from gecatsim.pyfiles.Detection_prefilter import Detection_prefilter
import unittest.mock
from unittest.mock import patch
from gecatsim.pyfiles.CommonTools import *
import numpy as np

class test_Detection_prefilter(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_Detection_prefilter(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.detector = CFG()
        cfg.scanner.detectorPrefilter = ['al', 0.1, 'water', 2]

        cfg.det.totalNumCells = 3

        cfg = feval(cfg.protocol.spectrumCallback, cfg)
        Wvec = Detection_prefilter(cfg)
        assert Wvec is not None