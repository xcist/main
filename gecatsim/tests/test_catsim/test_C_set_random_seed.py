import unittest
from unittest.mock import MagicMock, patch
from gecatsim.pyfiles.C_set_random_seed import C_set_random_seed
from gecatsim.pyfiles.CommonTools import CFG

cfg = CFG("../examples/cfg/Phantom_Sample_Analytic",
          "../examples/cfg/Scanner_Sample_generic",
          "../examples/cfg/Protocol_Sample_axial")
cfg.clib = MagicMock()
cfg.clib.setall = MagicMock()

def test_C_set_random_seed_with_seeds():
    # Reset the mock
    cfg.clib.setall.reset_mock()

    seed1 = 123
    seed2 = 456

    C_set_random_seed(cfg, seed1, seed2)

    cfg.clib.setall.assert_called_once_with(seed1, seed2)

def test_C_set_random_seed_without_seeds():
    # Reset the mock
    cfg.clib.setall.reset_mock()

    C_set_random_seed(cfg)

    cfg.clib.setall.assert_called_once_with(None, None)

if __name__ == '__main__':
    unittest.main()