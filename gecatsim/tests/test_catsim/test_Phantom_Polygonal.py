from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Polygonal import Phantom_Polygonal
from unittest.mock import patch, MagicMock

class Test_Phantom_NCAT(unittest.TestCase):
    def test_Phantom_NCAT(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample_Polygonal", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/female_adult_average_lung_lesions_reduced.nrb'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_material_info_polygon = MagicMock()
        cfg.clib.clear_polygonalized_phantom = MagicMock()
        cfg.clib.pass_polygon_to_c = MagicMock()

        cfg = Phantom_Polygonal(cfg)

        assert cfg.clib.set_material_info_polygon.call_count == 1
        assert cfg.clib.clear_polygonalized_phantom.call_count == 1
        assert cfg.clib.pass_polygon_to_c.call_count > 0

