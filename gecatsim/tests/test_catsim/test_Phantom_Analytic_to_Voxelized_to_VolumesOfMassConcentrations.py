import unittest
import numpy as np
from unittest.mock import MagicMock, patch
from gecatsim.pyfiles.CommonTools import CFG
from gecatsim.pyfiles.Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations import Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations

class TestPhantomAnalyticToVoxelizedToVolumesOfMassConcentrations(unittest.TestCase):

    @patch('os.path.exists')
    @patch('os.path.getmtime')
    @patch('Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.catvoxel')
    @patch('Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations.Phantom_Voxelized_to_VolumesOfMassConcentrations')
    def test_phantom_analytic_to_voxelized_to_volumes_of_mass_concentrations(self, mock_phantom_voxelized, mock_catvoxel, mock_getmtime, mock_exists):
        cfg = CFG("../examples/cfg/Phantom_Sample", "../examples/cfg/Scanner_Sample_generic",
                  "../examples/cfg/Protocol_Sample_axial", "../examples/cfg/Scanner_Sample_generic.cfg")

        cfg.phantom_filename = 'phantom.pp'
        cfg.phantom_samples_xy = 128
        cfg.phantom_samples_voxelsize = 1.0
        cfg.phantom_samples_z = 64
        cfg.force_phantom_conversion = True

        mock_exists.return_value = False
        mock_getmtime.return_value = 0
        mock_catvoxel.return_value = (np.zeros((cfg.phantom_samples_xy, cfg.phantom_samples_xy, cfg.phantom_samples_z)), [])

        Phantom_Analytic_to_Voxelized_to_VolumesOfMassConcentrations(cfg)

        self.assertEqual(cfg.phantom_filename, 'phantom.vp')
        mock_catvoxel.assert_called_once()
        mock_phantom_voxelized.assert_called_once()

if __name__ == '__main__':
    unittest.main()