from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Analytic import Phantom_Analytic
from unittest.mock import patch, MagicMock

#FIXME: need to fix this
## Getting 'UnboundLocalError: local variable 'ppmPhantomFilename' referenced before assignment


class Test_Phantom_Analytic(unittest.TestCase):
    def test_Phantom_Analytic(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/CatSimLogo_1024/CatSim_logo_1024.json'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_phantom_info = MagicMock()
        cfg.clib.set_bounding_info = MagicMock()

        cfg = Phantom_Analytic(cfg)

        assert cfg.clib.set_phantom_info.call_count == 1
        assert cfg.clib.set_bounding_info.call_count == 1