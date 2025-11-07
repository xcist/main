import unittest
from unittest.mock import MagicMock, patch
import ctypes
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.set_rndseed import set_rndseed

cfg = CFG("../examples/cfg/Phantom_Sample_Analytic",
          "../examples/cfg/Scanner_Sample_generic",
          "../examples/cfg/Protocol_Sample_axial")
cfg.clib = MagicMock()
cfg.clib.setall = MagicMock()

@patch('gecatsim.pyfiles.CommonTools.load_C_lib')
def test_set_rndseed_positive(mock_load_C_lib):
    # Call the function with a positive seed
    set_rndseed(cfg, 42)
    args, kwargs = cfg.clib.setall.call_args
    assert args[0].value == 430
    assert args[1].value == 43

@patch('gecatsim.pyfiles.CommonTools.load_C_lib')
def test_set_rndseed_negative(mock_load_C_lib):
    # Call the function with a negative seed
    with patch('builtins.print') as mocked_print:
        set_rndseed(cfg, -42)
        mocked_print.assert_called_with('Warning: Seed is set to negative, will be converted to positive')
    args, kwargs = cfg.clib.setall.call_args
    assert args[0].value == 430
    assert args[1].value == 43

@patch('gecatsim.pyfiles.CommonTools.load_C_lib')
def test_set_rndseed_zero(mock_load_C_lib):
    # Call the function with a seed of zero
    set_rndseed(cfg, 0)
    args, kwargs = cfg.clib.setall.call_args
    assert args[0].value == 10
    assert args[1].value == 1

if __name__ == '__main__':
    unittest.main()