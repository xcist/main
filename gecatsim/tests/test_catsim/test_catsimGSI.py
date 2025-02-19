import os
import unittest
from gecatsim.pyfiles.catsimGSI import catsimGSI
from unittest.mock import MagicMock, patch
from gecatsim.pyfiles.CommonTools import *

class TestCatSimGSI(unittest.TestCase):
    @patch('gecatsim.pyfiles.CommonTools.PathHelper.find')
    @patch('gecatsim.pyfiles.CommonTools.source_cfg')
    def test_catsimGSI(self, mock_source_cfg, mock_find):

        cfg = CFG("../examples/cfg/Phantom_Sample",
                          "../examples/cfg/Scanner_Sample_generic",
                          "../examples/cfg/Protocol_Sample_axial")
        cfg.kVp = 120
        cfg.rotation_period = 1.0
        cfg.views_per_rotation = 1000
        cfg.start_time = 0.0
        cfg.start_view = 1
        cfg.duty_cycle_per_view = [0.5, 0.5]
        cfg.subphase_grouping = [[1], [2]]
        cfg.kVp_cycle_per_view = [120, 100]
        cfg.number_Ebins = 100
        cfg.spectrum_filename = "spectrum_120"
        cfg.mA_cycle_per_view = [200, 150]
        cfg.this_is_an_airscan = False
        cfg.this_is_an_offsetscan = False
        cfg.total_n_cells = 100
        cfg.airscan_total_n_views = 10
        cfg.offsetscan_total_n_views = 10
        cfg.total_n_views = 10
        cfg.convert_to_prep_per_phase = 0
        cfg.results_basename = "test"
        cfg.start_time_GSI_mode = 1
        cfg.start_view = 1
        cfg.stop_view = 10

        # Mock other arguments
        adjust = MagicMock()
        preadjust = MagicMock()
        silent = MagicMock()

        catsimGSI(cfg, adjust, preadjust, silent)

        self.assertTrue(os.path.exists('test_subphase_1_subphase_2_phase_1.air'))
        self.assertTrue(os.path.exists('test_subphase_1_subphase_2_phase_1.offset'))
        self.assertTrue(os.path.exists('test_subphase_1_subphase_2_phase_1.scan'))

if __name__ == '__main__':
    unittest.main()