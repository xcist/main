import unittest.mock
from unittest.mock import patch

import numpy
from gecatsim.pyfiles.CommonTools import *
from gecatsim.pyfiles.Prep_BHC_Accurate import Prep_BHC_Accurate
import gecatsim as xc


class TestPrep_BHC_Accurate(unittest.TestCase):

    #@patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_prep_bhc_accurate(self):
        cfg = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial")

        cfg.phantom.filename = 'CTDI_16cm_WaterAirPEBoneChambers.ppm'

        cfg.resultsName = "test_Analytic"
        cfg.protocol.viewsPerRotation = 500
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1

        cfg.physics.enableQuantumNoise = 0
        cfg.physics.enableElectronicNoise = 0

        cfg.physics.callback_post_log = 'Prep_BHC_Accurate'
        cfg.physics.EffectiveMu = 0.2
        cfg.physics.BHC_poly_order = 5
        cfg.physics.BHC_max_length_mm = 220
        cfg.physics.BHC_length_step_mm = 10

        cfg.physics.colSampleCount = 2
        cfg.physics.rowSampleCount = 2
        cfg.physics.srcXSampleCount = 2
        cfg.physics.srcYSampleCount = 2
        cfg.physics.viewSampleCount = 1
        
        prep = numpy.zeros((500, 14400))

        cfg = cfg.air_scan(doPrint=False)

        prep = Prep_BHC_Accurate(cfg, prep)

        for i in range(500):
            for j in range(14400):
                assert prep[i,j] != 0

        ndet = cfg.scanner.detectorColCount * cfg.scanner.detectorRowCount
        poly_order = cfg.physics.BHC_poly_order
        num_poly_coefs = poly_order + 1

        if os.path.exists(cfg.resultsName+".air"):
            os.remove(cfg.resultsName+".air")