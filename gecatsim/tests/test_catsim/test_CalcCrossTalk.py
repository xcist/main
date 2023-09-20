import unittest.mock
from unittest.mock import patch

import numpy as np

import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.CalcCrossTalk import CalcCrossTalk


class TestCalcCrossTalk(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_calccrosstalk(self, signal_mock):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        thisView = np.zeros([14400, 20], dtype=np.double)

        cfg.physics.row_crosstalk = 0.02
        cfg.physics.col_crosstalk = 0.025

        signal_mock.return_value = np.zeros([cfg.scanner.detectorColCount, cfg.scanner.detectorRowCount],
                                            dtype=np.double)
        CalcCrossTalk(thisView, cfg)

        signal_mock.call_count == cfg.physics.energyCount
