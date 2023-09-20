import unittest.mock
from unittest.mock import patch, call

import gecatsim.pyfiles.OneScan as o
from gecatsim.pyfiles import CommonTools


class TestOneScan(unittest.TestCase):

    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_air_scan(self, faval_mock):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.sim.thisScanType = [1, 0, 0]

        # dummy values for the test to continue without failures
        cfg.detFlux = [[[1, 0, 0]], [[1, 0, 0]], [[1, 0, 0]]]
        faval_mock.return_value = cfg

        o.one_scan(cfg)

        expected = [call(cfg.scanner.detectorCallback, cfg), call(cfg.scanner.focalspotCallback, cfg),
                    call(cfg.physics.rayAngleCallback, cfg), call(cfg.protocol.spectrumCallback, cfg),
                    call(cfg.protocol.filterCallback, cfg), call(cfg.physics.fluxCallback, cfg),
                    call(cfg.scanner.detectionCallback, cfg, 0, 0), call(cfg.physics.outputCallback, cfg, 0)]
        assert faval_mock.mock_calls == expected

    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_offset_scan(self, faval_mock):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")
        cfg.sim.thisScanType = [0, 1, 0]
        # dummy values for the test to continue without failures
        cfg.detFlux = [[[1, 0, 0]], [[1, 0, 0]], [[1, 0, 0]]]
        faval_mock.return_value = cfg

        o.one_scan(cfg)

        expected = [call(cfg.scanner.detectorCallback, cfg), call(cfg.scanner.focalspotCallback, cfg),
                    call(cfg.physics.rayAngleCallback, cfg), call(cfg.protocol.spectrumCallback, cfg),
                    call(cfg.protocol.filterCallback, cfg), call(cfg.physics.fluxCallback, cfg),
                    call(cfg.scanner.detectionCallback, cfg, 0, 0), call(cfg.physics.outputCallback, cfg, 0)]
        assert faval_mock.mock_calls == expected

    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_phantom_scan(self, faval_mock):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")
        cfg.sim.thisScanType = [0, 0, 1]
        # set this a lesser value
        cfg.protocol.stopViewId = 5
        # dummy values for the test to continue without failures
        cfg.detFlux = [[[1, 0, 0]], [[1, 0, 0]], [[1, 0, 0]]]
        faval_mock.return_value = cfg

        o.one_scan(cfg)

        expected = [call(cfg.phantom.callback, cfg),
                    call(cfg.scanner.detectorCallback, cfg),
                    call(cfg.scanner.focalspotCallback, cfg),
                    call(cfg.physics.rayAngleCallback, cfg),
                    call(cfg.protocol.spectrumCallback, cfg),
                    call(cfg.protocol.filterCallback, cfg),
                    call(cfg.physics.fluxCallback, cfg)]

        for viewId in range(cfg.protocol.stopViewId + 1):
            for subViewId in range(2):
                expected.append(call(cfg.protocol.scanTrajectory, cfg, viewId))
                expected.append(call(cfg.phantom.projectorCallback, cfg, viewId, subViewId))
                expected.append(call(cfg.scanner.detectionCallback, cfg, viewId, subViewId))
            expected.append(call(cfg.physics.outputCallback, cfg, viewId))

        assert faval_mock.mock_calls == expected

    def test_update_scan_time(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                              "../examples/cfg/Protocol_Sample_axial")

        cfg.time = 100
        cfg.subViewTime = 50
        cfg.protocol.dutyRatio = 10
        cfg.sim.subViewCount = 20

        o.update_scan_time(cfg, 10)

        assert cfg.time == 150

        cfg.time = 100
        cfg.subViewTime = 50
        cfg.protocol.dutyRatio = 0.5
        cfg.sim.subViewCount = 11
        cfg.viewTime = 25

        o.update_scan_time(cfg, 10)

        assert cfg.time == 162.5
