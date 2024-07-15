from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Analytic import Phantom_Analytic
from unittest.mock import patch, MagicMock

class Test_Phantom_Analytic(unittest.TestCase):
    def test_Phantom_Analytic(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample_Analytic", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/W20.ppm'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_material_info = MagicMock()
        cfg.clib.set_phantom_info = MagicMock()

        cfg = Phantom_Analytic(cfg)

        assert cfg.clib.set_material_info.call_count == 1
        assert cfg.clib.set_phantom_info.call_count == 1