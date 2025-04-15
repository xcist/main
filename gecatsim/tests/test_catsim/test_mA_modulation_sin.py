import unittest
import numpy as np
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.mA_modulation_sin import mA_modulation_sin

class TestMAModulationSin(unittest.TestCase):

    def test_mA_modulation_sin(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.start_angle = 0
        cfg.views_per_rotation = 360
        cfg.rotation_direction = 1
        cfg.mA_modulation = 0.5

        viewnr = 1
        expected_scale_factor = 1.0 + cfg.mA_modulation * np.sin(cfg.start_angle * np.pi / 180.0)
        scale_factor = mA_modulation_sin(cfg, viewnr)
        self.assertAlmostEqual(scale_factor, expected_scale_factor, places=5)

if __name__ == '__main__':
    unittest.main()