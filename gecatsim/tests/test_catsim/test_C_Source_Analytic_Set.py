from gecatsim.pyfiles.C_Source_Analytic_Set import C_Source_Analytic_Set
import numpy as np
import unittest.mock
import gecatsim.pyfiles.CommonTools as c
from unittest.mock import patch, MagicMock
import os

class test_C_Source_Analytic_Set(unittest.TestCase):

    def test_C_Projector_Analytic_calls_clib_Projector(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic", "../examples/cfg/Protocol_Sample_axial")
        cfg.src.weights = np.array([[0.25, 0.25, 0.25, 0.25]], dtype=np.single)

        cfg.clib.set_src_info = MagicMock()
        C_Source_Analytic_Set(cfg,0,0)
        assert cfg.clib.set_src_info.call_count == 1
