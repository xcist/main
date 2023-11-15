import unittest.mock
from unittest.mock import patch

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detection_Flux import Detection_Flux


class TestDetection_Flux(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_detection_flux(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.sim.isOffsetScan = 0
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)

        cfg = feval(cfg.physics.rayAngleCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)
        cfg = feval(cfg.protocol.filterCallback, cfg)

        cfg = Detection_Flux(cfg)

        assert cfg.detFlux is not None
