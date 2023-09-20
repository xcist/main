import unittest.mock
from unittest.mock import patch

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detection_Lag import Detection_Lag


class TestDetection_Lag(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_detection_lag(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.physics.lag_taus = [0.9, 6]
        cfg.physics.lag_alphas = [0.930, 0.070]

        cfg.protocol.rotationTime = 1.0

        cfg.protocol.viewsPerRotation = 1000
        cfg.physics.viewSampleCount = 2

        cfg.sim.startViewId = 0
        cfg.sim.isAirScan = 1

        Detection_Lag(0, 0, 0, cfg)

        assert cfg.memview1 is not None
        assert cfg.memview2 is not None
