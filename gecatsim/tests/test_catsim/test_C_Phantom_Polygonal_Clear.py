from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.C_Phantom_Polygonal_Clear import C_Phantom_Polygonal_Clear
import unittest

from unittest.mock import patch, call, MagicMock


class Test_C_Phantom_Polygonal_Clear(unittest.TestCase):
    def test_C_Phantom_Polygonal_Clear(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        num_polygons = 5
        cfg.clib.clear_polygonalized_phantom = MagicMock()
        C_Phantom_Polygonal_Clear(cfg, num_polygons)

        assert cfg.clib.clear_polygonalized_phantom.call_count == 1