import unittest.mock
from unittest.mock import patch

import gecatsim as xc


class TestCatSim(unittest.TestCase):

    @patch('gecatsim.CatSim.prep_view', create=True)
    @patch('gecatsim.CatSim.phantom_scan', create=True)
    @patch('gecatsim.CatSim.offset_scan', create=True)
    @patch('gecatsim.CatSim.air_scan', create=True)
    def test_run_all(self, air_scan_mock, offset_scan_mock, phantom_scan_mock, prep_view_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.cfg.protocol.scanTypes = [1, 1, 1]
        air_scan_mock.return_value = ct.cfg

        ct.run_all()

        assert air_scan_mock.call_count == 1
        assert offset_scan_mock.call_count == 1
        assert phantom_scan_mock.call_count == 1
        assert prep_view_mock.call_count == 1

    @patch('gecatsim.CatSim.prep_view', create=True)
    @patch('gecatsim.CatSim.phantom_scan', create=True)
    @patch('gecatsim.CatSim.offset_scan', create=True)
    @patch('gecatsim.CatSim.air_scan', create=True)
    def test_run_all_air_scan_only(self, air_scan_mock, offset_scan_mock, phantom_scan_mock, prep_view_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.cfg.protocol.scanTypes = [1, 0, 0]
        air_scan_mock.return_value = ct.cfg

        ct.run_all()

        assert air_scan_mock.call_count == 1
        assert offset_scan_mock.call_count == 0
        assert phantom_scan_mock.call_count == 0
        assert prep_view_mock.call_count == 0

    @patch('gecatsim.CatSim.prep_view', create=True)
    @patch('gecatsim.CatSim.phantom_scan', create=True)
    @patch('gecatsim.CatSim.offset_scan', create=True)
    @patch('gecatsim.CatSim.air_scan', create=True)
    def test_run_all_offset_scan_only(self, air_scan_mock, offset_scan_mock, phantom_scan_mock, prep_view_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.cfg.protocol.scanTypes = [0, 1, 0]
        air_scan_mock.return_value = ct.cfg

        ct.run_all()

        assert air_scan_mock.call_count == 0
        assert offset_scan_mock.call_count == 1
        assert phantom_scan_mock.call_count == 0
        assert prep_view_mock.call_count == 0

    @patch('gecatsim.CatSim.prep_view', create=True)
    @patch('gecatsim.CatSim.phantom_scan', create=True)
    @patch('gecatsim.CatSim.offset_scan', create=True)
    @patch('gecatsim.CatSim.air_scan', create=True)
    def test_run_all_phantom_scan_only(self, air_scan_mock, offset_scan_mock, phantom_scan_mock, prep_view_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.cfg.protocol.scanTypes = [0, 0, 1]
        air_scan_mock.return_value = ct.cfg

        ct.run_all()

        assert air_scan_mock.call_count == 0
        assert offset_scan_mock.call_count == 0
        assert phantom_scan_mock.call_count == 1
        assert prep_view_mock.call_count == 0

    @patch('gecatsim.pyfiles.CatSim.one_scan', create=True)
    def test_air_scan(self, one_scan_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.air_scan(ct.get_current_cfg())

        one_scan_mock.assert_called_once_with(ct.get_current_cfg())

        assert ct.get_current_cfg().sim.thisScanType == [1, 0, 0]

    @patch('gecatsim.pyfiles.CatSim.one_scan', create=True)
    def test_offset_scan(self, one_scan_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.offset_scan(ct.get_current_cfg())

        one_scan_mock.assert_called_once_with(ct.get_current_cfg())
        assert ct.get_current_cfg().sim.thisScanType == [0, 1, 0]

    @patch('gecatsim.pyfiles.PrepView.prep_view', create=True)
    def test_phantom_scan(self, prep_view_mock):
        ct = xc.CatSim("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                       "../examples/cfg/Protocol_Sample_axial")  # initialization

        ct.prep_view(ct.get_current_cfg())

        prep_view_mock.assert_called_once_with(ct.get_current_cfg())
