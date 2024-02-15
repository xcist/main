from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Polygonal import phantom_polygonal
from unittest.mock import patch, MagicMock

from unittest.mock import patch, call

class Test_Phantom_Polygonal(unittest.TestCase):

    @patch('gecatsim.pyfiles.Phantom_Polygonal.feval', create=True)
    def test_Phantom_Polygonal(self, feval_mock):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample_Analytic")

        cfg.phantom.filename = '../phantom/W20.ppm'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_src_info_vox = MagicMock()
        cfg = phantom_polygonal(cfg)

        expected = [call("C_Phantom_Polygonal_Clear", cfg), call("C_Phantom_Polygonal_SetPolygon", cfg)]
        assert feval_mock.mock_calls == expected

        assert cfg.clib.set_src_info_vox.call_count == 1