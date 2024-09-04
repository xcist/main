import unittest.mock
from unittest.mock import patch

import numpy
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Scatter_Correction import Scatter_Correction
import gecatsim as xc


class TestScatter_ConvolutionModel(unittest.TestCase):

    @patch('gecatsim.pyfiles.Scatter_Correction.rawwrite', create=True)
    def test_scatter_convolutionmodel(self, scatter_raw_write):
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

        airscan = np.random.random([1, 57600])
        offsetScan  = np.random.random([1, 57600])
        phantomScan  = np.random.random([1, 57600])
        ct.scanner.detectorRowsPerMod = 64
        ct.scanner.detectorRowCount = ct.scanner.detectorRowsPerMod

        airscan, offsetScan, phantomScan = Scatter_Correction(ct, airscan, offsetScan, phantomScan)

        for i in range(57600):
            assert airscan[0, i] != 0
            assert offsetScan[0, i] != 0
            assert phantomScan[0, i] != 0

        assert scatter_raw_write.call_count == 1
