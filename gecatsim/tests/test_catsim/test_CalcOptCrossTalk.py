from gecatsim.pyfiles.CalcOptCrossTalk import CalcOptCrossTalk
import gecatsim as xc
import unittest.mock
from unittest.mock import patch
import gecatsim.pyfiles.CommonTools as c
import numpy as np

class test_CalcOptCrossTalk(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_CalcOptCrossTalk(self, signal_mock):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")
        thisView = np.zeros([900, 16], dtype=np.double)
        cfg.physics.row_crosstalk_opt = 0.045
        cfg.physics.col_crosstalk_opt = 0.04
        signal_mock.return_value = np.zeros([cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount], dtype=np.double)
        CalcOptCrossTalk(thisView, cfg)
        signal_mock.call_count == cfg.physics.energyCount
