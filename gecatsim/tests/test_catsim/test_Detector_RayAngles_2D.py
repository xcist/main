import unittest.mock
from unittest.mock import patch

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Detector_RayAngles_2D import Detector_RayAngles_2D


class TestDetector_RayAngles_2D(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_detector_ray_angles_2d(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg = feval("Detector_ThirdgenCurved", cfg)
        #cfg = feval("Source_Uniform", cfg)
        cfg.src.samples = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                       [2.500000e-01, 5.420361e+02, -2.500000e-01],
                                       [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                       [2.500000e-01, 5.379639e+02, 2.500000e-01]], dtype=np.single)

        cfg.src.nSamples = 4

        cfg = Detector_RayAngles_2D(cfg)

        assert cfg.det.rayDistance is not None
        assert cfg.det.alphas is not None
        assert cfg.det.cosBetas is not None
        assert cfg.det.gammas is not None
