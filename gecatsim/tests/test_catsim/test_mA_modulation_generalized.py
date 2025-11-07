import unittest
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.mA_modulation_generalized import *


class TestMAModulationGeneralized(unittest.TestCase):

    def test_mA_modulation_with_modulation(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.max_mA = 300
        cfg.min_mA = 100
        cfg.start_angle = 0
        cfg.views_per_rotation = 360
        cfg.rotation_direction = 1
        cfg.mA_modulation = 0.5
        cfg.mA = 200

        viewnr = 1
        expected_scale_factor = 1.25
        scale_factor = mA_modulation_generalized(cfg, viewnr)
        self.assertAlmostEqual(scale_factor, expected_scale_factor, places=5)

    def test_mA_modulation_with_table(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.modulation_table_file_path = 'mod_table.npz'
        cfg.views_per_rotation = 360
        cfg.mA = 200

        mod_table = np.ones((360, 1)) * 200
        np.savez(cfg.modulation_table_file_path, mod_table=mod_table)

        viewnr = 1
        expected_scale_factor = 1.0
        scale_factor = mA_modulation_generalized(cfg, viewnr)
        self.assertAlmostEqual(scale_factor, expected_scale_factor, places=5)

if __name__ == '__main__':
    unittest.main()