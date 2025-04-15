import unittest
from unittest.mock import patch, mock_open
import numpy as np
from gecatsim.pyfiles import CommonTools
from gecatsim.pyfiles.CreateHeelEffect import CreateHeelEffect, HeelEffectIntensity

cfg = CommonTools.CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                      "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

class TestCreateHeelEffect(unittest.TestCase):

    @patch('gecatsim.pyfiles.Spectrum.spectrum_read')
    @patch('gecatsim.pyfiles.GetMu.GetMu')
    @patch('gecatsim.pyfiles.GetMu.ReadMaterialFile')
    @patch('builtins.open', new_callable=mock_open, read_data='3\n10,0.9048374180359595\n20,0.8187307530779818\n30,0.7408182206817179\n')
    def test_create_heel_effect(self, mock_open, mock_ReadMaterialFile, mock_GetMu, mock_spectrum_read):

        mock_spectrum_read.return_value = (np.array([1, 2, 3]), np.array([10, 20, 30]), 0)
        mock_GetMu.return_value = np.array([0.1, 0.2, 0.3])
        mock_ReadMaterialFile.return_value = (1, 19.3, [74], [1.0])


        cfg.heel_effect_limit_angle = 10
        cfg.heel_effect_angle_decimation = 5
        cfg.target_angle = 15
        cfg.spectrum_filename = 'test_spectrum.dat'
        cfg.spectrum_dir = 'test_dir'
        cfg.anode_electron_penetration_in_mm = 1
        cfg.reference_spectrum_angle = 5
        cfg.tube_tilt = 2

        CreateHeelEffect(cfg)

        # Filter out irrelevant paths
        actual_filenames = [call[0][0] for call in mock_open.call_args_list if 'test_dir' in call[0][0]]

        expected_filenames = [
            'test_dir/test_spectrum_000_-10.dat',
            'test_dir/test_spectrum_001_-5.dat',
            'test_dir/test_spectrum_002_0.dat',
            'test_dir/test_spectrum_003_5.dat',
            'test_dir/test_spectrum_004_10.dat'
        ]
        self.assertEqual(expected_filenames, actual_filenames)

        handle = mock_open()

        handle.write.assert_any_call('3\n')
        handle.write.assert_any_call('0.904837429523468,10.0\n')
        handle.write.assert_any_call('0.8187307715415955,20.0\n')
        handle.write.assert_any_call('0.740818202495575,30.0\n')
        handle.write.assert_any_call('7.00')


class TestHeelEffectIntensity(unittest.TestCase):
    @patch('gecatsim.pyfiles.GetMu.GetMu')
    def test_heel_effect_intensity(self, mock_GetMu):
        mock_GetMu.return_value = np.array([0.1, 0.2, 0.3])

        theta = np.array([-10, 0, 10])
        Energy = np.array([10, 20, 30])

        result = HeelEffectIntensity(theta, Energy, cfg)

        expected_result = np.array([
            [1.00000e+000, 1.00000e+000, 2.81366e+106],
            [8.90066e-049, 2.44266e+019, 5.74422e+119],
            [2.54395e-216, 4.97641e+006, 2.32803e+041]
        ])

        self.assertTrue(np.allclose(result, expected_result, rtol=1e-5, atol=1e-8))

if __name__ == '__main__':
    unittest.main()