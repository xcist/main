import unittest
from unittest.mock import patch, call
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.catdkvp import catdkvp

class TestCatdkvp(unittest.TestCase):

    @patch('gecatsim.CatSim')
    @patch('builtins.print')
    def test_catdkvp(self, mock_print, mock_CatSim):

        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.rotate_rotate = 1
        cfg.results_basename = "test_basename"
        configfilename = 'configfile'

        catdkvp(configfilename, cfg, 'adjust', 'preadjust', 'silent')

        expected_print_calls = [
            call('\n Simulating low kVp scan...\n'),
            call('\n... done simulating low kVp scan.\n'),
            call('\n Simulating high kVp scan...\n'),
            call('\n... done simulating high kVp scan.\n')
        ]
        mock_print.assert_has_calls(expected_print_calls, any_order=False)

        # Check CatSim calls
        expected_CatSim_calls = [
            call(configfilename, cfg, 'cfg["rr_kvp"]="low_kvp"; cfg["results_basename"]="test_basename_lowkvp";', 'preadjust', 'silent'),
            call(configfilename, cfg, 'cfg["rr_kvp"]="high_kvp"; cfg["results_basename"]="test_basename_highkvp";', 'preadjust', 'silent')
        ]
        mock_CatSim.assert_has_calls(expected_CatSim_calls, any_order=False)

if __name__ == '__main__':
    unittest.main()