from gecatsim.pyfiles.CommonTools import *
import unittest.mock
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.Phantom_Analytic import Phantom_Analytic,Rmat, euler_angs, Phantom_Analytic_ConvertTori,get_clip_dD
from unittest.mock import patch, MagicMock

class Test_Phantom_Analytic(unittest.TestCase):
    def test_Phantom_Analytic(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample_Analytic", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = '../phantom/W20.ppm'
        cfg.sim.subViewCount = 1

        cfg = feval(cfg.scanner.detectorCallback, cfg)
        cfg = feval(cfg.scanner.focalspotCallback, cfg)
        cfg = feval(cfg.protocol.spectrumCallback, cfg)

        cfg.clib.set_material_info = MagicMock()
        cfg.clib.set_phantom_info = MagicMock()

        cfg = Phantom_Analytic(cfg)

        assert cfg.clib.set_material_info.call_count == 1
        assert cfg.clib.set_phantom_info.call_count == 1

    def test_get_clip_dD(self):
        c0 = [0, 0, 0]
        c1 = [0, 0, 0]
        result = get_clip_dD(c0, c1)
        expected = 0.0
        self.assertEqual(result, expected)
    def test_Rmat(self):
        p = [30, 45, 60]
        expected_output = np.array([
            [0.12682648, 0.9267767, 0.35355339],
            [-0.78033009, -0.12682648, 0.61237244],
            [0.61237244, -0.35355339, 0.70710678]
        ])
        np.testing.assert_almost_equal(Rmat(p), expected_output, decimal=6)

    def test_euler_angs(self):
        R = np.array([
            [0.36, -0.48, 0.8],
            [0.8, 0.6, 0],
            [-0.48, 0.64, 0.6]
        ])
        expected_ea = [np.arctan2(0.8, 0.6), np.arccos(0.6), np.arctan2(-0.48, 0.64)]
        result = euler_angs(R)
        np.testing.assert_almost_equal(result, expected_ea, decimal=5)