import unittest.mock
from unittest.mock import MagicMock

import numpy as np

import gecatsim.pyfiles.CommonTools as c
from gecatsim.pyfiles.C_Projector_Voxelized import C_Projector_Voxelized


class TestC_Projector_Voxelized(unittest.TestCase):

    def test_c_projector_analytic_calls_clib_voxelized_projector(self):
        cfg = c.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                    "../examples/cfg/Protocol_Sample_axial")

        cfg.detNew.totalNumCells = 14400
        cfg.spec.nEbin = 20

        cfg.srcNew.weights = np.array([[0.25, 0.25, 0.25, 0.25]], dtype=np.single)
        cfg.srcNew.samples = np.array([[-2.500000e-01, 5.420361e+02, -2.500000e-01],
                                       [2.500000e-01, 5.420361e+02, -2.500000e-01],
                                       [-2.500000e-01, 5.379639e+02, 2.500000e-01],
                                       [2.500000e-01, 5.379639e+02, 2.500000e-01]], dtype=np.single)
        cfg.srcNew.nSamples = 4

        cfg.detNew.modTypes = np.zeros((900, 1), dtype=np.single)
        cfg.detNew.vvecs = np.zeros((900, 3), dtype=np.single)
        cfg.detNew.uvecs = np.zeros((900, 3), dtype=np.single)
        cfg.detNew.modCoords = np.zeros((900, 3), dtype=np.single)

        cfg.phantom.numberOfMaterials = 2
        cfg.sim.stopViewId = 49

        for i in range(0, 900):
            cfg.detNew.modTypes[i] = np.array([0], dtype=np.single)
            cfg.detNew.vvecs[i] = np.array([0, 0, 1], dtype=np.single)
            cfg.detNew.uvecs[i] = np.array([0.8902536, -0.4554652, 0.0], dtype=np.single)
            cfg.detNew.modCoords[i] = np.array([-432.69196, -305.7409, 0.0], dtype=np.single)

        cfg.detNew.nMod = 900
        cfg.thisSubView = np.zeros((14400, 20), dtype=np.single)

        for i in range(0, 14400):
            cfg.thisSubView[i] = [0.0, 0.0, 0.0, 5.8932144e-11, 0.0022313918, 2.8641064, 85.15998, 503.93176, 1335.6624,
                                  3907.3804, 3119.3367, 4286.113, 3069.4067, 3123.5225, 2996.3645, 2717.0347, 2307.0815,
                                  1784.581, 1145.0819, 386.48895]

        cfg.clib.voxelized_projector = MagicMock()

        C_Projector_Voxelized(cfg, 0, 0)

        assert cfg.clib.voxelized_projector.call_count == 8
