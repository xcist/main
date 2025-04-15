import unittest
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.mA_modulation_Rev2Rot import mA_modulation_Rev2Rot

class TestMAModulationRev2Rot(unittest.TestCase):

    def test_mA_modulation_Rev2Rot_within_rotation(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.views_per_rotation = 360
        cfg.mA = 200

        viewnr = 180
        expected_scale_factor = 1
        scale_factor = mA_modulation_Rev2Rot(cfg, viewnr)
        self.assertEqual(scale_factor, expected_scale_factor)

    def test_mA_modulation_Rev2Rot_beyond_rotation(self):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.views_per_rotation = 360
        cfg.mA = 200

        viewnr = 361
        expected_scale_factor = 100 / cfg.mA
        scale_factor = mA_modulation_Rev2Rot(cfg, viewnr)
        self.assertEqual(scale_factor, expected_scale_factor)

if __name__ == '__main__':
    unittest.main()