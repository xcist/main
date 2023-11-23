from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_NCAT import Phantom_NCAT
from unittest.mock import patch, MagicMock

#FIXME: need to fix this

#Getting 'allocation error in cp_matrix'


class Test_Phantom_NCAT(unittest.TestCase):
    def test_Phantom_NCAT(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/CatSimLogo_1024/CatSim_logo_1024.json'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_src_info_NCAT = MagicMock()
        cfg.clib.set_material_info_NCAT = MagicMock()
        cfg.clib.set_module_info_NCAT = MagicMock()

        cfg = Phantom_NCAT(cfg)

        assert cfg.clib.set_src_info_NCAT.call_count == 1
        assert cfg.clib.set_material_info_NCAT.call_count == 1
        assert cfg.clib.set_module_info_NCAT == 1

