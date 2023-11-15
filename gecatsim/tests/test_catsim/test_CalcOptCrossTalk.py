import unittest.mock
from unittest.mock import patch

import numpy as np

import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.CalcOptCrossTalk import CalcOptCrossTalk


class TestCalcOptCrossTalk(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_calcoptcrosstalk(self, signal_mock):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        thisView = np.zeros([900, 16], dtype=np.double)

        cfg.physics.row_crosstalk_opt = 0.045
        cfg.physics.col_crosstalk_opt = 0.04

        signal_mock.return_value = np.zeros([cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount],
                                            dtype=np.double)

        CalcOptCrossTalk(thisView, cfg)

        signal_mock.call_count == cfg.physics.energyCount
