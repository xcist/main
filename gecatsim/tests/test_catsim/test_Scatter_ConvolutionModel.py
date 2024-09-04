import unittest.mock
from unittest.mock import patch

import numpy
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Scatter_ConvolutionModel import Scatter_ConvolutionModel
import gecatsim as xc


class TestScatter_ConvolutionModel(unittest.TestCase):

    #@patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_scatter_convolutionmodel(self):
        ct = xc.CatSim()  # initialization

        ##--------- Make changes to parameters (optional)
        ct.resultsName = "test_sc"
        ct.physics.scatterCallback = "Scatter_ConvolutionModel"  # scatter model
        ct.physics.scatterKernelCallback = ""  # using default
        ct.physics.scatterScaleFactor = 1

        ct.physics.callback_pre_log = "Scatter_Correction"  # scatter correction
        ct.physics.scatterCorrectionScaleFactor = 1
        ct.physics.scatterCorrectionSaveView = 1  # save the corrected phantom scan

        ct.protocol.viewsPerRotation = 1
        ct.protocol.viewCount = ct.protocol.viewsPerRotation
        ct.protocol.stopViewId = ct.protocol.viewCount - 1
        ct.sim = emptyCFG
        ct.sim.startViewId = 0

        ct.thisSubView = np.random.random([57600, 20])
        ct.detFlux  = np.random.random([57600, 20])
        ct.scanner.detectorRowsPerMod = 64
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

        cfg = Scatter_ConvolutionModel(ct, 0, 0)

        for i in range(57600):
            for j in range(20):
                assert cfg.thisSubView[i,j] != 0
