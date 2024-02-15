import unittest.mock
from unittest.mock import ANY
from unittest.mock import MagicMock

import ctypes
import numpy as np

import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.C_Phantom_Polygonal_SetPolygon import C_Phantom_Polygonal_SetPolygon
class test_C_Phantom_Polygonal_SetPolygon(unittest.TestCase):

    def test_c_phantom_polygonal_setpolygon_calls_pass_polygon_to_c(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")
        vertices = np.array([[0.25, 0.25, 0.25, 0.25]], dtype=np.single)

        cfg.clib.pass_polygon_to_c = MagicMock()

        C_Phantom_Polygonal_SetPolygon(cfg, vertices, 4, 100.0, 10)

        assert cfg.clib.pass_polygon_to_c.call_count == 1
        cfg.clib.pass_polygon_to_c.assert_called_once_with(ANY, 4, 100.0, 10)
        assert (cfg.clib.pass_polygon_to_c.call_args[0][0] == np.array([[0.25, 0.25, 0.25, 0.25]])).all()

        assert cfg.clib.pass_polygon_to_c.argtypes == [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_int]
        assert cfg.clib.pass_polygon_to_c.restype == None
