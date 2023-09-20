import unittest.mock
from unittest.mock import patch

import numpy as np

from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Gantry_Helical import Gantry_Helical


class test_Gantry_Helical(unittest.TestCase):

    @patch('gecatsim.pyfiles.CalcCrossTalk.signal.convolve2d', create=True)
    def test_Gantry_Helical(self, signal_mock):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.src.corners = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                    [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                    [2.500000e-01, 5.379639e+02, 2.500000e-01],
                                    [2.500000e-01, 5.420361e+02, -2.500000e-01]], dtype=np.single)
        cfg.src.samples = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                    [2.500000e-01, 5.420361e+02, -2.500000e-01],
                                    [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                    [2.500000e-01, 5.379639e+02, 2.500000e-01]], dtype=np.single)
        cfg.src.front = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                  [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                  [2.500000e-01, 5.379639e+02, 2.500000e-01],
                                  [2.500000e-01, 5.420361e+02, -2.500000e-01]], dtype=np.single)
        cfg.src.lateral = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                    [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                    [2.500000e-01, 5.379639e+02, 2.500000e-01],
                                    [2.500000e-01, 5.420361e+02, -2.500000e-01]], dtype=np.single)
        cfg.src.long = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                 [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                 [2.500000e-01, 5.379639e+02, 2.500000e-01],
                                 [2.500000e-01, 5.420361e+02, -2.500000e-01]], dtype=np.single)

        cfg.det.modCoords = np.zeros((900, 3), dtype=np.single)

        cfg.det.uvecs = np.zeros((900, 3), dtype=np.single)
        cfg.det.vvecs = np.zeros((900, 3), dtype=np.single)
        cfg.det.nMod = 900

        cfg.time = 100

        cfg.src.nCorners = 4
        cfg.src.nSamples = 4

        cfg = Gantry_Helical(cfg, 0)

        assert cfg.srcNew is not None
        assert cfg.srcNew.corners is not None
        assert cfg.srcNew.samples is not None
        assert cfg.srcNew.front is not None
        assert cfg.srcNew.lateral is not None
        assert cfg.srcNew.long is not None

        assert cfg.detNew is not None
        assert cfg.detNew.modCoords is not None
        assert cfg.detNew.uvecs is not None
        assert cfg.detNew.vvecs is not None
