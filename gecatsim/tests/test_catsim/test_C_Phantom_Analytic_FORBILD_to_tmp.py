import unittest.mock
from unittest.mock import MagicMock
import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.C_Phantom_Analytic_FORBILD_to_tmp import C_Phantom_Analytic_FORBILD_to_tmp


class TestC_Phantom_Analytic_FORBILD_to_tmp(unittest.TestCase):

    def test_C_Phantom_Analytic_FORBILD_to_tmp(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.clib.TranslatePhantom_FORBILD_to_tmp = MagicMock()

        scale = 1.0
        pp_phantom_filename = "input.pp"
        tmp_phantom_filename = "output.tmp"

        C_Phantom_Analytic_FORBILD_to_tmp(cfg, scale, pp_phantom_filename, tmp_phantom_filename)

        assert cfg.clib.TranslatePhantom_FORBILD_to_tmp.call_count == 1


if __name__ == '__main__':
    unittest.main()