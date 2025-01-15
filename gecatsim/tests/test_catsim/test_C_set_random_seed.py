import unittest
from unittest.mock import MagicMock, patch
from gecatsim.pyfiles.C_set_random_seed import C_set_random_seed
from gecatsim.pyfiles.CommonTools import *
class TestCSetRandomSeed(unittest.TestCase):
    def setUp(self):

        self.cfg = CFG("../examples/cfg/Phantom_Sample_Analytic",
                       "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")
        self.cfg.clib.setall = MagicMock()

    def test_C_set_random_seed_with_seeds(self):
        seed1 = 123
        seed2 = 456

        C_set_random_seed(self.cfg, seed1, seed2)

        self.cfg.clib.setall.assert_called_once_with(seed1, seed2)

    def test_C_set_random_seed_without_seeds(self):
        C_set_random_seed(self.cfg)

        self.cfg.clib.setall.assert_called_once_with(None, None)

if __name__ == '__main__':
    unittest.main()