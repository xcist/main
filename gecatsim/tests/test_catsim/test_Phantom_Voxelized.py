from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Voxelized import Phantom_Voxelized
from unittest.mock import patch, MagicMock

class Test_Phantom_Voxelized(unittest.TestCase):
    def test_Phantom_Voxelized(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/CatSimLogo_1024/CatSim_logo_1024.json'

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_src_info_vox = MagicMock()
        cfg = Phantom_Voxelized(cfg)

        assert cfg.clib.set_src_info_vox.call_count == 1