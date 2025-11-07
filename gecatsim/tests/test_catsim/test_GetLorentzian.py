import unittest
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.GetLorentzian import GetLorentzian

class TestGetLorentzian(unittest.TestCase):

    def test_GetLorentzian(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.lorentzian_a = 1.0
        cfg.lorentzian_b = 0.5
        cfg.lorentzian_c = 0.3

        ptch = 1.0
        kernelsize = 8
        expected_H = GetLorentzian(ptch, kernelsize, cfg)
        H = GetLorentzian(ptch, kernelsize, cfg)
        np.testing.assert_almost_equal(H, expected_H, decimal=5)

if __name__ == '__main__':
    unittest.main()