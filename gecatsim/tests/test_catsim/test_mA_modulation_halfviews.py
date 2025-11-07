import unittest
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.mA_modulation_halfviews import mA_modulation_halfviews

class TestMAModulationHalfviews(unittest.TestCase):

    def test_mA_modulation_halfviews(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.start_view = 0
        cfg.views_per_rotation = 360

        viewnr_in_range = cfg.start_view + cfg.views_per_rotation // 4
        viewnr_out_of_range = cfg.start_view + cfg.views_per_rotation // 4 + cfg.views_per_rotation // 2 + 1

        expected_scale_factor_in_range = 2
        expected_scale_factor_out_of_range = 1

        scale_factor_in_range = mA_modulation_halfviews(cfg, viewnr_in_range)
        scale_factor_out_of_range = mA_modulation_halfviews(cfg, viewnr_out_of_range)

        self.assertEqual(scale_factor_in_range, expected_scale_factor_in_range)
        self.assertEqual(scale_factor_out_of_range, expected_scale_factor_out_of_range)

if __name__ == '__main__':
    unittest.main()