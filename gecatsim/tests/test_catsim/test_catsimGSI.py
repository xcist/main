import unittest
from unittest.mock import patch, MagicMock
import numpy as np
from gecatsim.pyfiles.catsimGSI import catsimGSI
from gecatsim.pyfiles.CommonTools import CFG

class TestCatsimGSI(unittest.TestCase):

    @patch('gecatsim.pyfiles.CommonTools.rawread')
    @patch('gecatsim.CatSim')
    @patch('gecatsim.pyfiles.PrepView.prep_view')
    def test_catsimGSI(self, mock_prep_view, mock_CatSim, mock_rawread):

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
        cfg.kVp_cycle_per_view = [80, 100]
        cfg.number_Ebins = 100
        cfg.spectrum_filename = "spectrum_120"
        cfg.mA_cycle_per_view = [200, 300]
        cfg.this_is_an_airscan = False
        cfg.this_is_an_offsetscan = False
        cfg.total_n_cells = 512
        cfg.airscan_total_n_views = 1000
        cfg.offsetscan_total_n_views = 1000
        cfg.total_n_views = 1000
        cfg.convert_to_prep_per_phase = 1
        cfg.results_basename = "test"
        cfg.start_time_GSI_mode = 1
        cfg.start_view = 0
        cfg.stop_view = 999

        # Mock return values for rawread
        mock_rawread.return_value = np.zeros((cfg.total_n_cells, 1, cfg.total_n_views))

        catsimGSI(cfg, adjust=None, preadjust=None, silent=True)

        self.assertEqual(mock_CatSim.call_count, 2)
        self.assertEqual(mock_rawread.call_count, 6)
        self.assertEqual(mock_prep_view.call_count, 2000)

if __name__ == '__main__':
    unittest.main()