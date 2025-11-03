import unittest
from unittest.mock import MagicMock, patch
from gecatsim.pyfiles.C_Phantom_Analytic_FORBILD_to_tmp import C_Phantom_Analytic_FORBILD_to_tmp

class TestC_Phantom_Analytic_FORBILD_to_tmp(unittest.TestCase):

    @patch('gecatsim.pyfiles.C_Phantom_Analytic_FORBILD_to_tmp.load_C_lib')
    def test_C_Phantom_Analytic_FORBILD_to_tmp(self, mock_load_C_lib):

        mock_clib = MagicMock()
        mock_func = MagicMock()
        mock_clib.TranslatePhantom_FORBILD_to_tmp = mock_func
        mock_load_C_lib.return_value = mock_clib

        scale = 1.0
        pp_phantom_filename = "input.pp"
        tmp_phantom_filename = "output.tmp"

        C_Phantom_Analytic_FORBILD_to_tmp(scale, pp_phantom_filename, tmp_phantom_filename)

        mock_func.assert_called_once_with(scale, b"input.pp", b"output.tmp")

if __name__ == '__main__':
    unittest.main()