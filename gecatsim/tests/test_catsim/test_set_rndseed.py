import unittest
from unittest.mock import MagicMock, patch
import ctypes
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.set_rndseed import set_rndseed

class TestSetRndSeed(unittest.TestCase):
    def setUp(self):

        self.cfg = CFG("../examples/cfg/Phantom_Sample_Analytic",
                       "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")
        self.cfg.clib = MagicMock()
        self.cfg.clib.setall = MagicMock()

    @patch('gecatsim.pyfiles.CommonTools.load_C_lib')
    def test_set_rndseed_positive(self, mock_load_C_lib):
        # Call the function with a positive seed
        set_rndseed(self.cfg, 42)
        args, kwargs = self.cfg.clib.setall.call_args
        self.assertEqual(args[0].value, 430)
        self.assertEqual(args[1].value, 43)

    @patch('gecatsim.pyfiles.CommonTools.load_C_lib')
    def test_set_rndseed_negative(self, mock_load_C_lib):
        # Call the function with a negative seed
        with patch('builtins.print') as mocked_print:
            set_rndseed(self.cfg, -42)
            mocked_print.assert_called_with('Warning: Seed is set to negative, will be converted to positive')
        args, kwargs = self.cfg.clib.setall.call_args
        self.assertEqual(args[0].value, 430)
        self.assertEqual(args[1].value, 43)

    @patch('gecatsim.pyfiles.CommonTools.load_C_lib')
    def test_set_rndseed_zero(self, mock_load_C_lib):
        # Call the function with a seed of zero
        set_rndseed(self.cfg, 0)
        args, kwargs = self.cfg.clib.setall.call_args
        self.assertEqual(args[0].value, 10)
        self.assertEqual(args[1].value, 1)

if __name__ == '__main__':
    unittest.main()