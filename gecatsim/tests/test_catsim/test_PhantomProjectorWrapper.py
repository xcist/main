import unittest.mock
from unittest.mock import patch, call

from gecatsim.pyfiles.PhantomProjectorWrapper import CopyCfgPhantom, SplitCfgPhantom, PhantomWrapper, ProjectorWrapper
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.CommonTools import *


class TestOneScan(unittest.TestCase):

    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_copy_cfg_phantom_creates_copy(self, feval_mock):
        cfg = self.setup()

        copiedCfg = CopyCfgPhantom(cfg)

        assert cfg != copiedCfg
        assert cfg.phantom != copiedCfg.phantom

        cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
        for i in range(len(cfg.phantom.filename)):
            for thisattr in cfgphantom_attrs:
                self.assertListEqual(getattr(cfg.phantom, thisattr), getattr(copiedCfg.phantom, thisattr))



    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_copy_cfg_phantom_updates_copy(self, feval_mock):
        cfg = self.setup()
        copiedCfg = CFG()

        returned_config = CopyCfgPhantom(cfg, copiedCfg)

        assert cfg != copiedCfg
        assert returned_config == copiedCfg
        assert cfg.phantom != copiedCfg.phantom

        cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
        for i in range(len(cfg.phantom.filename)):
            for thisattr in cfgphantom_attrs:
                self.assertListEqual(getattr(cfg.phantom, thisattr), getattr(copiedCfg.phantom, thisattr))

    @patch('gecatsim.pyfiles.OneScan.feval', create=True)
    def test_split_cfg_phantom(self, feval_mock):
        cfg = self.setup()

        cfg_list = SplitCfgPhantom(cfg)

        cfgphantom_attrs = [x for x in dir(cfg.phantom) if not x.startswith('__')]
        for i in range(len(cfg.phantom.filename)):
            for thisattr in cfgphantom_attrs:
                for split_config in cfg_list:
                    assert getattr(split_config.phantom, thisattr) in getattr(cfg.phantom, thisattr)

    @patch('gecatsim.pyfiles.PhantomProjectorWrapper.feval', create=True)
    def test_phantom_wrapper(self, faval_mock):
        cfg = self.setup()

        cfg = PhantomWrapper(cfg)

        cfg_list = SplitCfgPhantom(cfg)

        expected = []
        actual = []

        for split_cfg in cfg_list:
            expected.append(split_cfg.phantom.callback)


        for mock_call in faval_mock.mock_calls:
            actual.append(mock_call.args[0])

        self.assertListEqual(actual, expected)

    @patch('gecatsim.pyfiles.PhantomProjectorWrapper.feval', create=True)
    def test_phantom_wrapper_single_projector(self, faval_mock):
        cfg = self.setup()

        cfg_list = SplitCfgPhantom(cfg)

        cfg = PhantomWrapper(cfg_list[0])

        assert len(faval_mock.mock_calls) == 1

    @patch('gecatsim.pyfiles.PhantomProjectorWrapper.sys.exit', create=True)
    def test_phantom_wrapper_duplicate_projector(self, sys_exit_mock):
        cfg = self.setup()

        if not isinstance(cfg.phantom.callback, list):
            cfg.phantom.callback = self.get_as_list(cfg.phantom.callback)
        cfg.phantom.callback.append(cfg.phantom.callback[0])

        try:
            cfg = PhantomWrapper(cfg)
            assert False
        except BaseException:
            assert True

        assert len(sys_exit_mock.mock_calls) == 1


    @patch('gecatsim.pyfiles.PhantomProjectorWrapper.feval', create=True)
    def test_projector_wrapper(self, faval_mock):
        cfg = self.setup()

        cfg = ProjectorWrapper(cfg, 0, 0)

        cfg_list = SplitCfgPhantom(cfg)

        expected = []
        actual = []

        for split_cfg in cfg_list:
            expected.append(split_cfg.phantom.projectorCallback)


        for mock_call in faval_mock.mock_calls:
            actual.append(mock_call.args[0])

        self.assertListEqual(actual, expected)

    @patch('gecatsim.pyfiles.PhantomProjectorWrapper.feval', create=True)
    def test_projector_wrapper_single_projector(self, faval_mock):
        cfg = self.setup()

        cfg_list = SplitCfgPhantom(cfg)

        cfg = ProjectorWrapper(cfg_list[0], 0, 0)

        assert len(faval_mock.mock_calls) == 1




    def get_as_list(self, object):
        if isinstance(object, list):
            return object
        return [object]
    def setup(self):
        cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample_Hybrid")

        cfg.resultsName = "test_Hybrid"
        cfg.protocol.viewsPerRotation = 500
        cfg.protocol.viewCount = cfg.protocol.viewsPerRotation
        cfg.protocol.stopViewId = cfg.protocol.viewCount - 1

        cfg.physics.enableQuantumNoise = 0
        cfg.physics.enableElectronicNoise = 0

        cfg.physics.callback_post_log = 'Prep_BHC_Accurate'
        cfg.physics.EffectiveMu = 0.2
        cfg.physics.BHC_poly_order = 5
        cfg.physics.BHC_max_length_mm = 300
        cfg.physics.BHC_length_step_mm = 10

        cfg.physics.colSampleCount = 2
        cfg.physics.rowSampleCount = 2
        cfg.physics.srcXSampleCount = 2
        cfg.physics.srcYSampleCount = 2
        cfg.physics.viewSampleCount = 1

        return cfg

